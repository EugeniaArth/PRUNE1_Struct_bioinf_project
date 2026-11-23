#!/usr/bin/env python3
"""Extract AlphaFold Server complex metrics for one selected chain pair.

Designed for the PRUNE1-GSK3B analysis, but works with other two-chain
selections in an AlphaFold Server job directory (including a pair from a
trimer). It reads all five summary/full-data JSONs and model CIFs; it does not
submit jobs or modify the input directory.

Example:
  python3 analyze_prune1_gsk3b.py fold_2026_10_22_16_06 \
      --chain-a A --chain-b B --label-a PRUNE1 --label-b GSK3B \
      --cut 0.5 --out results_PRUNE1_GSK3B

For PRUNE1 homodimer, use --chain-a A --chain-b B and labels PRUNE1-A/PRUNE1-B.
For the trimer, run separately for A:C and B:C. Output directories must not
already exist; this avoids silently overwriting prior results.
"""

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


# ---------------------------- Input and parsing

def find_samples(job_dir):
    """Find the five AlphaFold Server output sets and verify matching files."""
    samples = []
    for summary_path in sorted(job_dir.glob("*_summary_confidences_*.json")):
        match = re.search(r"_summary_confidences_(\d+)\.json$", summary_path.name)
        if not match:
            continue
        index = int(match.group(1))
        prefix = summary_path.name[:match.start()]
        samples.append({
            "index": index,
            "summary_path": summary_path,
            "full_path": job_dir / f"{prefix}_full_data_{index}.json",
            "cif_path": job_dir / f"{prefix}_model_{index}.cif",
        })
    samples.sort(key=lambda s: s["index"])
    if len(samples) != 5:
        raise ValueError(f"Expected 5 AlphaFold Server samples; found {len(samples)}")
    for sample in samples:
        for key in ("full_path", "cif_path"):
            if not sample[key].is_file():
                raise FileNotFoundError(sample[key])
        with sample["summary_path"].open() as handle:
            sample["summary"] = json.load(handle)
        with sample["full_path"].open() as handle:
            sample["full"] = json.load(handle)
    return samples


_TOKEN = re.compile(r"'(?:[^']|'(?=\S))*'|\"(?:[^\"]|\"(?=\S))*\"|\S+")


def parse_cif_atoms(path):
    """Read heavy atoms from an AlphaFold Server mmCIF atom_site loop."""
    columns = []
    atoms = []
    in_loop = False
    reading = False
    with path.open() as handle:
        for line in handle:
            text = line.strip()
            if text == "loop_":
                in_loop, reading, columns = True, False, []
                continue
            if in_loop and text.startswith("_atom_site."):
                columns.append(text.split()[0][len("_atom_site."):])
                continue
            if columns and in_loop and not reading:
                reading = True
            if not reading:
                continue
            if not text or text.startswith("#") or text.startswith("_") or text == "loop_":
                break
            fields = [token.strip("'\"") for token in _TOKEN.findall(text)]
            row = dict(zip(columns, fields))
            if row.get("group_PDB") not in ("ATOM", "HETATM") or row.get("type_symbol") == "H":
                continue
            seq = row.get("label_seq_id", ".")
            if seq in (".", "?"):
                seq = row.get("auth_seq_id", "1")
            atoms.append({
                "chain": row["label_asym_id"],
                "res": int(seq),
                "comp": row.get("label_comp_id", "UNK"),
                "atom": row.get("label_atom_id", ""),
                "xyz": (float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])),
            })
    if not atoms:
        raise ValueError(f"No atom_site coordinates found in {path}")
    return atoms


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def token_coordinates(atoms, token_chains, token_residues):
    """Return C-alpha and C-beta/pseudo-beta coordinates in token order."""
    by_residue = defaultdict(dict)
    for atom in atoms:
        by_residue[(atom["chain"], atom["res"])].setdefault(atom["atom"], atom["xyz"])
    ca = np.full((len(token_chains), 3), np.nan, dtype=np.float32)
    cb = np.full_like(ca, np.nan)
    for i, (chain, residue) in enumerate(zip(token_chains, token_residues)):
        names = by_residue.get((chain, int(residue)), {})
        if not names:
            continue
        ca[i] = names.get("CA", next(iter(names.values())))
        cb[i] = names.get("CB", ca[i])  # glycine: pseudo-beta is CA
    return ca, cb


def residue_plddt(atoms, atom_plddts):
    if len(atoms) != len(atom_plddts):
        raise ValueError(f"CIF atom count ({len(atoms)}) != atom_plddts length ({len(atom_plddts)})")
    values = defaultdict(list)
    for atom, score in zip(atoms, atom_plddts):
        values[(atom["chain"], atom["res"])].append(float(score))
    return {key: float(np.mean(scores)) for key, scores in values.items()}


def distances_for_sample(sample, token_chains, token_residues, mask_a, mask_b):
    atoms = parse_cif_atoms(sample["cif_path"])
    ca, cb = token_coordinates(atoms, token_chains, token_residues)
    ia, ib = np.where(mask_a)[0], np.where(mask_b)[0]
    d_ca = np.linalg.norm(ca[ia, None, :] - ca[None, ib, :], axis=2)
    d_cb = np.linalg.norm(cb[ia, None, :] - cb[None, ib, :], axis=2)
    return d_ca, d_cb, residue_plddt(atoms, sample["full"]["atom_plddts"])


def block_summary(block, cutoffs):
    profile_a = block.max(axis=1)
    profile_b = block.max(axis=0)
    return {
        "shape": [int(x) for x in block.shape],
        "max_contact_probability": float(block.max()),
        "cutoffs": {
            str(cut): {
                "pairs": int((block > cut).sum()),
                "chain_a_residues": int((profile_a > cut).sum()),
                "chain_b_residues": int((profile_b > cut).sum()),
            } for cut in cutoffs
        },
    }


def tsv_write(path, columns, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------- Analysis

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("job_dir", type=Path, help="AlphaFold Server result folder")
    parser.add_argument("--chain-a", default="A", help="first chain in contact/PAE blocks (default: A)")
    parser.add_argument("--chain-b", default="B", help="second chain in contact/PAE blocks (default: B)")
    parser.add_argument("--label-a", default=None, help="display label for chain A; defaults to its chain ID")
    parser.add_argument("--label-b", default=None, help="display label for chain B; defaults to its chain ID")
    parser.add_argument("--cut", type=float, default=0.5, help="contact_probs threshold (default: 0.5)")
    parser.add_argument("--distance-cut", type=float, default=8.0, help="distance threshold in Angstrom (default: 8)")
    parser.add_argument("--out", type=Path, default=None, help="new output directory; must not already exist")
    parser.add_argument("--dry-run", action="store_true", help="calculate and print the summary without writing files")
    args = parser.parse_args()

    job_dir = args.job_dir.resolve()
    if not job_dir.is_dir():
        parser.error(f"Not a directory: {job_dir}")
    if args.chain_a == args.chain_b:
        parser.error("Select two different chains")
    label_a, label_b = args.label_a or args.chain_a, args.label_b or args.chain_b
    samples = find_samples(job_dir)
    # Sort for top-sample analysis; retain sample indices for coordinate comparisons.
    top = max(samples, key=lambda s: float(s["summary"]["ranking_score"]))
    full = top["full"]
    token_chains = np.asarray(full["token_chain_ids"])
    token_residues = np.asarray(full["token_res_ids"], dtype=int)
    cp = np.asarray(full["contact_probs"], dtype=float)
    pae = np.asarray(full["pae"], dtype=float)
    if cp.shape != (len(token_chains), len(token_chains)) or pae.shape != cp.shape:
        raise ValueError("contact_probs/pae dimensions do not match token_chain_ids")
    mask_a, mask_b = token_chains == args.chain_a, token_chains == args.chain_b
    if not mask_a.any() or not mask_b.any():
        present = list(dict.fromkeys(token_chains.tolist()))
        raise ValueError(f"Requested chains {args.chain_a}/{args.chain_b}; found {present}")
    if mask_a.sum() != len(set(token_residues[mask_a])) or mask_b.sum() != len(set(token_residues[mask_b])):
        raise ValueError("This script expects one token per protein residue for the selected chains")

    block_cp = cp[np.ix_(mask_a, mask_b)]
    block_pae_ab = pae[np.ix_(mask_a, mask_b)]
    block_pae_ba = pae[np.ix_(mask_b, mask_a)].T
    residues_a, residues_b = token_residues[mask_a], token_residues[mask_b]
    cutoffs = [0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.7, 0.9]
    counts = block_summary(block_cp, cutoffs)

    # Measure all selected-chain residue pairs in all five coordinate models.
    distance_ca, distance_cb, plddt_by_sample = [], [], []
    for sample in samples:
        ca, cb, plddt = distances_for_sample(sample, token_chains.tolist(), token_residues.tolist(), mask_a, mask_b)
        distance_ca.append(ca)
        distance_cb.append(cb)
        plddt_by_sample.append(plddt)
    distance_ca = np.stack(distance_ca)
    distance_cb = np.stack(distance_cb)
    freq_ca = np.mean(distance_ca < args.distance_cut, axis=0)
    freq_cb = np.mean(distance_cb < args.distance_cut, axis=0)
    cp_diff_by_sample = []
    for sample in samples:
        f = sample["full"]
        tc = np.asarray(f["token_chain_ids"])
        sample_cp = np.asarray(f["contact_probs"], dtype=float)
        sample_block = sample_cp[np.ix_(tc == args.chain_a, tc == args.chain_b)]
        cp_diff_by_sample.append(float(np.max(np.abs(sample_block - block_cp))))

    chain_ids = list(dict.fromkeys(top["summary"].get("chain_ids", token_chains.tolist())))
    summary = top["summary"]
    pair_i, pair_j = chain_ids.index(args.chain_a), chain_ids.index(args.chain_b)
    models_table = []
    for sample in sorted(samples, key=lambda s: s["index"]):
        q = sample["summary"]
        models_table.append({
            "sample": sample["index"], "ptm": float(q["ptm"]), "iptm": float(q["iptm"]),
            "ranking_score": float(q["ranking_score"]), "fraction_disordered": float(q["fraction_disordered"]),
            "has_clash": q["has_clash"],
        })

    # Per-residue summary above selected cutoff; best partner is on chain B.
    profile_a = block_cp.max(axis=1)
    profile_b = block_cp.max(axis=0)
    top_position = next(i for i, sample in enumerate(samples) if sample is top)
    top_plddt = plddt_by_sample[top_position]
    residue_rows = []
    for i in np.where(profile_a > args.cut)[0]:
        j = int(np.argmax(block_cp[i]))
        ra, rb = int(residues_a[i]), int(residues_b[j])
        residue_rows.append({
            "chain_a": args.chain_a, "residue_a": ra,
            "max_contact_probs": round(float(profile_a[i]), 4),
            "best_chain_b": args.chain_b, "best_residue_b": rb,
            "best_pair_pae_a_aligned": round(float(block_pae_ab[i, j]), 2),
            "best_pair_pae_b_aligned": round(float(block_pae_ba[i, j]), 2),
            "chain_a_mean_pLDDT_top": round(top_plddt.get((args.chain_a, ra), float("nan")), 1),
            "chain_b_mean_pLDDT_top": round(top_plddt.get((args.chain_b, rb), float("nan")), 1),
        })

    # All pairs above cutoff, with coordinate ranges across all five samples.
    pair_rows = []
    for i, j in np.argwhere(block_cp > args.cut):
        ra, rb = int(residues_a[i]), int(residues_b[j])
        pair_rows.append({
            "chain_a": args.chain_a, "residue_a": ra, "chain_b": args.chain_b, "residue_b": rb,
            "contact_probs": round(float(block_cp[i, j]), 4),
            "pae_a_aligned": round(float(block_pae_ab[i, j]), 2),
            "pae_b_aligned": round(float(block_pae_ba[i, j]), 2),
            "ca_min_A": round(float(np.nanmin(distance_ca[:, i, j])), 2),
            "ca_max_A": round(float(np.nanmax(distance_ca[:, i, j])), 2),
            "cb_min_A": round(float(np.nanmin(distance_cb[:, i, j])), 2),
            "cb_max_A": round(float(np.nanmax(distance_cb[:, i, j])), 2),
            "fraction_samples_CA_under_cut": round(float(freq_ca[i, j]), 2),
            "ca_distances_by_sample_A": ",".join(f"{x:.2f}" for x in distance_ca[:, i, j]),
            "cb_distances_by_sample_A": ",".join(f"{x:.2f}" for x in distance_cb[:, i, j]),
        })
    pair_rows.sort(key=lambda row: -row["contact_probs"])

    # Largest mismatches between the predicted probability and five-coordinate C-alpha ruler.
    gap = np.abs(block_cp - freq_ca)
    flat_order = np.argsort(gap.ravel())[::-1][:20]
    disagreements = []
    for flat in flat_order:
        i, j = np.unravel_index(flat, gap.shape)
        disagreements.append({
            "chain_a": args.chain_a, "residue_a": int(residues_a[i]),
            "chain_b": args.chain_b, "residue_b": int(residues_b[j]),
            "contact_probs": round(float(block_cp[i, j]), 4),
            "fraction_samples_CA_under_cut": round(float(freq_ca[i, j]), 2),
            "fraction_samples_CB_under_cut": round(float(freq_cb[i, j]), 2),
            "ca_distances_by_sample_A": [round(float(x), 2) for x in distance_ca[:, i, j]],
            "cb_distances_by_sample_A": [round(float(x), 2) for x in distance_cb[:, i, j]],
            "pae_a_aligned": round(float(block_pae_ab[i, j]), 2),
            "pae_b_aligned": round(float(block_pae_ba[i, j]), 2),
            "absolute_cp_vs_ca_frequency_gap": round(float(gap[i, j]), 3),
        })

    report = {
        "job_directory": job_dir.name,
        "top_sample": top["index"],
        "selected_pair": {"chain_a": args.chain_a, "label_a": label_a, "chain_b": args.chain_b, "label_b": label_b},
        "chain_order_in_summary": chain_ids,
        "top_scores": {k: float(summary[k]) for k in ("ptm", "iptm", "ranking_score", "fraction_disordered")},
        "all_samples": models_table,
        "chain_ptm": {c: float(v) for c, v in zip(chain_ids, summary["chain_ptm"])},
        "chain_iptm": {c: float(v) for c, v in zip(chain_ids, summary["chain_iptm"])},
        "chain_pair_iptm": summary["chain_pair_iptm"],
        "chain_pair_pae_min_A": summary["chain_pair_pae_min"],
        "selected_pair_chain_pair_iptm": float(summary["chain_pair_iptm"][pair_i][pair_j]),
        "selected_pair_chain_pair_pae_min_a_aligned_A": float(summary["chain_pair_pae_min"][pair_i][pair_j]),
        "selected_pair_chain_pair_pae_min_b_aligned_A": float(summary["chain_pair_pae_min"][pair_j][pair_i]),
        "selected_pair_contact_block": counts,
        "selected_pair_pae_block_medians_A": {
            "a_aligned_b_scored": float(np.median(block_pae_ab)),
            "b_aligned_a_scored": float(np.median(block_pae_ba)),
        },
        "contact_probs_max_abs_difference_vs_top_by_sample": cp_diff_by_sample,
        "selected_cutoff": args.cut,
        "distance_cutoff_A": args.distance_cut,
        "called_residue_count_chain_a": len(residue_rows),
        "called_pair_count": len(pair_rows),
        "called_residues_chain_a": residue_rows,
        "top_20_ruler_disagreements": disagreements,
    }

    # Concise human-readable console summary; dry-run produces no files.
    print(f"Job: {job_dir.name}; selected pair: {label_a}({args.chain_a}) x {label_b}({args.chain_b})")
    print(f"Top sample: model_{top['index']} | pTM {summary['ptm']:.2f} | ipTM {summary['iptm']:.2f} | ranking_score {summary['ranking_score']:.2f}")
    print(f"Selected-pair ipTM {report['selected_pair_chain_pair_iptm']:.2f}; chain_pair_pae_min: aligned on {label_a} {report['selected_pair_chain_pair_pae_min_a_aligned_A']:.2f} A, aligned on {label_b} {report['selected_pair_chain_pair_pae_min_b_aligned_A']:.2f} A")
    print("chain_ptm:", report["chain_ptm"])
    print("chain_iptm:", report["chain_iptm"])
    print("chain_pair_iptm, row/column order", chain_ids, ":", summary["chain_pair_iptm"])
    print("chain_pair_pae_min (A), row = aligned on; row/column order", chain_ids, ":", summary["chain_pair_pae_min"])
    print(f"Contact block {counts['shape']}; max contact_probs {counts['max_contact_probability']:.2f}")
    print("cutoff | pairs | residues A | residues B")
    for cutoff in cutoffs:
        row = counts["cutoffs"][str(cutoff)]
        print(f"{cutoff:>6} | {row['pairs']:>5} | {row['chain_a_residues']:>10} | {row['chain_b_residues']:>10}")
    print(f"Contact_probs max differences across samples: {cp_diff_by_sample}")
    print("Top residue calls:")
    for row in residue_rows[:20]:
        print(f"  {row['residue_a']} -> {row['best_residue_b']}: cp={row['max_contact_probs']:.2f}, PAE={row['best_pair_pae_a_aligned']:.1f}/{row['best_pair_pae_b_aligned']:.1f} A")
    print("Largest ruler/probability differences (C-alpha ruler vs five-sample fraction):")
    for row in disagreements[:5]:
        print(f"  {row['residue_a']} x {row['residue_b']}: cp={row['contact_probs']:.2f}, CA distances={row['ca_distances_by_sample_A']}, CA<cut fraction={row['fraction_samples_CA_under_cut']:.2f}")

    if args.dry_run:
        return

    out_dir = args.out or Path(f"analysis_{args.chain_a}_{args.chain_b}")
    if out_dir.exists():
        raise FileExistsError(f"Output path already exists; choose a new --out directory: {out_dir}")
    out_dir.mkdir(parents=True)
    (out_dir / "results.json").write_text(json.dumps(report, indent=2))

    residue_columns = list(residue_rows[0].keys()) if residue_rows else [
        "chain_a", "residue_a", "max_contact_probs", "best_chain_b", "best_residue_b",
        "best_pair_pae_a_aligned", "best_pair_pae_b_aligned", "chain_a_mean_pLDDT_top", "chain_b_mean_pLDDT_top",
    ]
    pair_columns = list(pair_rows[0].keys()) if pair_rows else [
        "chain_a", "residue_a", "chain_b", "residue_b", "contact_probs", "pae_a_aligned", "pae_b_aligned",
        "ca_min_A", "ca_max_A", "cb_min_A", "cb_max_A", "fraction_samples_CA_under_cut",
        "ca_distances_by_sample_A", "cb_distances_by_sample_A",
    ]
    tsv_write(out_dir / "residues_above_cutoff.tsv", residue_columns, residue_rows)
    tsv_write(out_dir / "pairs_above_cutoff.tsv", pair_columns, pair_rows)
    tsv_write(out_dir / "ruler_disagreements.tsv", list(disagreements[0].keys()), disagreements)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.8), constrained_layout=True)
    matrices = [
        (block_cp, "contact_probs", 0, 1, "magma", "probability"),
        (block_pae_ab, f"PAE: aligned on {label_a}", 0, 32, "viridis", "PAE (Å)"),
        (block_pae_ba, f"PAE: aligned on {label_b}", 0, 32, "viridis", "PAE (Å)"),
    ]
    extent = (residues_b[0] - 0.5, residues_b[-1] + 0.5, residues_a[0] - 0.5, residues_a[-1] + 0.5)
    for ax, (matrix, title, vmin, vmax, cmap, color_label) in zip(axes, matrices):
        image = ax.imshow(matrix, origin="lower", aspect="auto", interpolation="nearest",
                          extent=extent, vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_title(title)
        ax.set_xlabel(f"{label_b} residue ({args.chain_b})")
        ax.set_ylabel(f"{label_a} residue ({args.chain_a})")
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03, label=color_label)
    fig.suptitle(f"{label_a}–{label_b}, model_{top['index']} (top-ranked); contact cutoff={args.cut}")
    fig.savefig(out_dir / "contact_and_directional_pae.png", dpi=180)
    plt.close(fig)
    print(f"\nWrote results to {out_dir}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, KeyError, IndexError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
