#!/usr/bin/env python3
"""
Usage:

  python3 structure_analysis.py --uniprot Q86TP1
  python3 structure_analysis.py --pdb 1DPX --compare-pdb 1LYZ

WHAT THIS SCRIPT DOES
    Reads one structure, measures the things Lecture 2 was about, and writes
        report_<ID>.md          a short text report with every number in it
        panels_<ID>.png         five figures: Ramachandran, secondary structure,
                                B-factor/pLDDT, burial, and contacts
    Nothing is hard-coded to the example protein.

WHAT YOU NEED
    python -m pip install biopython matplotlib numpy
    The file rama_allowed.npz and the module sbtools.py must sit next to this file.
"""

import argparse
import os

# ─────────────── Parse CLI ───────────────
parser = argparse.ArgumentParser(description="Analyse a protein structure")
parser.add_argument("--pdb", default=None, help="PDB entry ID (e.g. 1DPX)")
parser.add_argument("--uniprot", default=None, help="UniProt accession for AlphaFold model")
parser.add_argument("--chain", default=None, help="Chain ID (default: first chain >=20 residues)")
parser.add_argument("--compare-pdb", default=None, help="Optional second PDB to compare")
parser.add_argument("--compare-chain", default=None, help="Chain for comparison structure")
parser.add_argument("--compare-uniprot", default=None, help="UniProt for comparison (AlphaFold)")
parser.add_argument("--buried-cutoff", type=float, default=0.05, help="Buried RSA threshold")
parser.add_argument("--exposed-cutoff", type=float, default=0.25, help="Exposed RSA threshold")
args = parser.parse_args()

PDB_ID = args.pdb
UNIPROT = args.uniprot
CHAIN = args.chain
COMPARE_PDB_ID = args.compare_pdb
COMPARE_CHAIN = args.compare_chain
COMPARE_ALPHAFOLD = args.compare_uniprot
BURIED_CUTOFF = args.buried_cutoff
EXPOSED_CUTOFF = args.exposed_cutoff

if not PDB_ID and not UNIPROT:
    parser.error("Provide either --pdb or --uniprot (or both).")

# ─────────────── Everything below is unchanged from the original ───────────────

import textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sbtools as sb

ACC, ORA, DK, GREY = "#156082", "#E97132", "#0E2841", "#9aa3ad"
GREEN, BLUE2, PURPLE, RED = "#196B24", "#0F9ED5", "#A02B93", "#C0392B"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "savefig.dpi": 200, "figure.facecolor": "white"})

report = []      # every line we print also goes into the report file


def say(line=""):
    print(line)
    report.append(line)


# -------------------------------------------------------------- 1. get the file
if PDB_ID:
    path = sb.fetch_pdb(PDB_ID)
    label = PDB_ID.upper()
else:
    if not UNIPROT:
        raise SystemExit("Set either PDB_ID or UNIPROT in the CONFIG block.")
    path = sb.fetch_alphafold(UNIPROT)
    label = f"AF-{UNIPROT}"

model = sb.load(path)
facts = sb.header_facts(path)
is_model = label.startswith("AF-")
chains = sb.chain_ids(model)
chain = CHAIN or (chains[0] if chains else None)
if chain is None:
    raise SystemExit("No chain with at least 20 amino acids was found in this file.")

res = sb.residues(model, chain)
say(f"# Homework 2 — {label}, chain {chain}")
say()
say(f"- file: `{os.path.basename(path)}`")
if facts.get("title"):
    say(f"- title: {facts['title'].title()}")
say(f"- method: {facts.get('method', 'predicted model (AlphaFold DB)' if is_model else 'not stated')}")
if facts.get("resolution"):
    say(f"- resolution: {facts['resolution']}")
say(f"- models in the file: {sb.n_models(path)}  (more than one means an NMR ensemble)")
say(f"- chains with >= 20 residues: {', '.join(chains)}")
say(f"- residues analysed in chain {chain}: {len(res)}")
if facts.get("ssbond"):
    say(f"- SSBOND records (disulfide bridges): {facts['ssbond']}")
say()

# --------------------------------------------- 2. backbone angles, Ramachandran
import os, sys
_script_dir = os.path.dirname(os.path.abspath(__file__))
_rama_candidates = [
    "rama_allowed.npz",
    os.path.join(_script_dir, "rama_allowed.npz"),
    os.path.join(_script_dir, "..", "data", "rama_allowed.npz"),
]
_rama_path = next((p for p in _rama_candidates if os.path.exists(p)), "rama_allowed.npz")
rama = sb.Rama(_rama_path)
pp = sb.phi_psi(model, chain)
gly = pp["resname"] == "GLY"
pro = pp["resname"] == "PRO"
gen = ~gly & ~pro
allowed = rama.is_allowed(pp["phi"], pp["psi"])
bad = rama.is_outlier(pp["phi"], pp["psi"]) & gen
outliers = [(int(r), n, round(float(a)), round(float(b)), round(float(v), 2))
            for r, n, a, b, v in zip(pp["resid"][bad], pp["resname"][bad],
                                     pp["phi"][bad], pp["psi"][bad],
                                     rama.overlap(pp["phi"], pp["psi"])[bad])]

say("## Backbone conformation")
say()
say(f"- residues with both φ and ψ defined: {len(pp['phi'])}")
say(f"- inside the strict hard-sphere band from the lecture (0.30 Å tolerance): "
    f"{100 * allowed.mean():.1f}%")
say(f"- non-Gly, non-Pro residues outside the loose band (0.50 Å) — the ones actually "
    f"worth inspecting: {len(outliers)}")
if 100 * allowed.mean() < 90:
    say("    - note: a low number on the strict band usually means your protein is "
        "helix-rich, not badly built. Hard spheres clip the lower edge of the α-helix "
        "basin, which is why the outlier call uses the looser band.")
for r, n, a, b, v in outliers[:10]:
    say(f"    - {n}{r}: φ = {a:+.0f}°, ψ = {b:+.0f}°, overlap {v} Å")
gly_left = [(int(r), round(float(a)), round(float(b)))
            for r, a, b in zip(pp["resid"][gly], pp["phi"][gly], pp["psi"][gly]) if a > 0]
say(f"- glycines in a left-handed conformation (φ > 0): {len(gly_left)} of {int(gly.sum())}"
    + (f" — {', '.join('Gly' + str(r) for r, _, _ in gly_left[:8])}" if gly_left else ""))
cis = [(r, n, round(float(w))) for r, n, w in sb.omega(model, chain) if abs(w) < 90]
say(f"- cis peptide bonds (|ω| < 90°): {len(cis)}"
    + (f" — {', '.join(f'{n}{r} (ω = {w}°)' for r, n, w in cis)}" if cis else ""))
say()

# ------------------------------------------------------ 3. secondary structure
ss, ss_resids, hbonds = sb.secondary_structure(model, chain)
frac = {k: ss.count(k) / max(len(ss), 1) for k in "HEL"}
say("## Secondary structure")
say()
say(f"- α-helix: {100 * frac['H']:.0f}%   β-strand: {100 * frac['E']:.0f}%   "
    f"neither: {100 * frac['L']:.0f}%")
runs, i = [], 0
while i < len(ss):
    j = i
    while j < len(ss) and ss[j] == ss[i]:
        j += 1
    if ss[i] in "HE":
        runs.append((ss[i], int(ss_resids[i]), int(ss_resids[j - 1]), j - i))
    i = j
say(f"- elements: {sum(1 for r in runs if r[0] == 'H')} helices, "
    f"{sum(1 for r in runs if r[0] == 'E')} strands")
longest = max(runs, key=lambda r: r[3]) if runs else None
if longest:
    say(f"- longest element: {'helix' if longest[0] == 'H' else 'strand'} "
        f"{longest[1]}–{longest[2]} ({longest[3]} residues)")
say(f"- backbone hydrogen bonds found: {len(hbonds)}")
say()
say("Secondary structure string (H = helix, E = strand, L = loop):")
say()
say("```")
for k in range(0, len(ss), 60):
    say(f"{int(ss_resids[k]):>5}  {ss[k:k + 60]}")
say("```")
say()

# -------------------------------------------------------- 4. burial and the core
rsa = sb.relative_sasa(model, chain)
keys = list(rsa)
vals = np.array([rsa[k] for k in keys])
buried = [k for k, v in zip(keys, vals) if v < BURIED_CUTOFF]
exposed = [k for k, v in zip(keys, vals) if v > EXPOSED_CUTOFF]
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY", "CYS"}
say("## Burial")
say()
say(f"- buried (relative SASA < {BURIED_CUTOFF}): {len(buried)} residues "
    f"({100 * len(buried) / len(keys):.0f}%)")
say(f"- exposed (relative SASA > {EXPOSED_CUTOFF}): {len(exposed)} residues "
    f"({100 * len(exposed) / len(keys):.0f}%)")
if buried:
    hb_frac = sum(1 for k in buried if k[2] in HYDROPHOBIC) / len(buried)
    ex_frac = (sum(1 for k in exposed if k[2] in HYDROPHOBIC) / len(exposed)
               if exposed else float("nan"))
    say(f"- hydrophobic fraction: {100 * hb_frac:.0f}% of buried residues, "
        f"{100 * ex_frac:.0f}% of exposed ones")
order = np.argsort(vals)
say("- most buried residues: "
    + ", ".join(f"{keys[i][2]}{keys[i][1]} ({vals[i]:.2f})" for i in order[:5]))
say("- most exposed residues: "
    + ", ".join(f"{keys[i][2]}{keys[i][1]} ({vals[i]:.2f})" for i in order[::-1][:5]))
say()

# ------------------------------------------------------- 5. contacts / interface
iface = []
if len(chains) > 1:
    partner = [c for c in chains if c != chain][0]
    iface = sb.interface_residues(model, chain, partner)
    say("## Interface")
    say()
    say(f"- chain {chain} against chain {partner}: {len(iface)} residues within 4.5 Å")
    if iface:
        rsa_alone = sb.relative_sasa(model, chain, isolate_chain=True)
        in_cplx = [rsa.get((chain, r, n), float("nan")) for r, n, _ in iface]
        alone = [rsa_alone.get((chain, r, n), float("nan")) for r, n, _ in iface]
        say(f"- mean relative SASA of those residues in the complex: "
            f"{np.nanmean(in_cplx):.2f}")
        say(f"- mean relative SASA of the same residues on the isolated chain: "
            f"{np.nanmean(alone):.2f}")
        say(f"- so the interface buries {100 * (1 - np.nanmean(in_cplx) / max(np.nanmean(alone), 1e-9)):.0f}% "
            f"of the surface those residues expose when the chain is alone. On a monomer "
            f"they look like ordinary surface — which is why a predictor that only sees "
            f"one chain calls substitutions there tolerated.")
        say("- interface residues: "
            + ", ".join(f"{n}{r}" for r, n, _ in iface[:25])
            + (" …" if len(iface) > 25 else ""))
    say()
else:
    say("## Interface")
    say()
    say(f"- only one chain in this file, so there is no interface to measure. If your "
        f"protein is known to act as a dimer or a larger complex, say so, and check "
        f"whether the PDB entry lists a biological assembly larger than what the file "
        f"contains (REMARK 350).")
    say()

# ---------------------------------------------------------- 6. comparison, if any
cmp_rows = []
if COMPARE_PDB_ID or COMPARE_ALPHAFOLD:
    if COMPARE_PDB_ID:
        cpath = sb.fetch_pdb(COMPARE_PDB_ID)
        clabel = COMPARE_PDB_ID.upper()
    else:
        cpath = sb.fetch_alphafold(COMPARE_ALPHAFOLD)
        clabel = f"AF-{COMPARE_ALPHAFOLD}"
    cmodel = sb.load(cpath)
    cchain = COMPARE_CHAIN or sb.chain_ids(cmodel)[0]
    core = [k[1] for k, v in rsa.items() if v < 0.25]     # the well-ordered core
    try:
        all_rms, n_all = sb.rmsd(model, chain, cmodel, cchain)
        core_rms, n_core = sb.rmsd(model, chain, cmodel, cchain, resid_subset=core)
        cmp_rows = [("all matched residues", all_rms, n_all),
                    ("buried core only (rel. SASA < 0.25)", core_rms, n_core)]
        say(f"## Comparison with {clabel} chain {cchain}")
        say()
        for name, v, n in cmp_rows:
            say(f"- RMSD over {name}: {v:.2f} Å  (n = {n} Cα)")
        say(f"- difference: {all_rms - core_rms:+.2f} Å. If this is large, the "
            f"whole-chain number is being driven by a few mobile residues, not by the fold.")
        say()
    except ValueError as e:
        say(f"## Comparison with {clabel}")
        say()
        say(f"- could not match residues: {e}")
        say()

# ------------------------------------------------------------------- figures
bvals = np.array([r["CA"].get_bfactor() for r in res])
bresid = np.array([r.id[1] for r in res])

fig = plt.figure(figsize=(16, 9.5))
gs = fig.add_gridspec(3, 3, height_ratios=[1.25, 1, 1], hspace=0.45, wspace=0.28)

# panel A — Ramachandran
ax = fig.add_subplot(gs[0, 0])
rama.background(ax)
ax.scatter(pp["phi"][gen], pp["psi"][gen], s=9, color=DK, alpha=0.7, linewidths=0,
           label="general", zorder=2)
ax.scatter(pp["phi"][gly], pp["psi"][gly], s=16, color=BLUE2, alpha=0.9, linewidths=0,
           label="Gly", zorder=3)
ax.scatter(pp["phi"][pro], pp["psi"][pro], s=16, color=ORA, alpha=0.9, linewidths=0,
           label="Pro", zorder=3)
if outliers:
    ax.scatter([o[2] for o in outliers], [o[3] for o in outliers], s=60,
               facecolor="none", edgecolor=RED, lw=1.4, zorder=4, label="outlier")
ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.text(-175, -170, "shaded = sterically allowed\ndarker: hard spheres · lighter: +0.2 Å",
        fontsize=8.5, color=DK, va="bottom")
ax.set_title(f"Ramachandran — {label} chain {chain}\n"
             f"{100 * allowed.mean():.1f}% inside the allowed region", fontsize=11)

# panel B — secondary structure track
ax = fig.add_subplot(gs[0, 1:])
colmap = {"H": "#d8563c", "E": "#e8b23a", "L": "#c8ced6"}
for i, s in enumerate(ss):
    ax.axvspan(ss_resids[i] - 0.5, ss_resids[i] + 0.5, color=colmap[s], lw=0)
ax.set_yticks([])
ax.set_xlim(ss_resids.min() - 1, ss_resids.max() + 1)
ax.set_xlabel("residue number")
for k, (lab, col) in enumerate([("helix", colmap["H"]), ("strand", colmap["E"]),
                                ("loop", colmap["L"])]):
    ax.text(0.01 + 0.085 * k, 1.04, lab, transform=ax.transAxes, color=col,
            fontweight="bold", fontsize=10.5)
ax.set_title(f"Secondary structure along the chain — "
             f"{100 * frac['H']:.0f}% H, {100 * frac['E']:.0f}% E, "
             f"{100 * frac['L']:.0f}% neither", fontsize=11, pad=24)

# panel C — B-factor or pLDDT
ax = fig.add_subplot(gs[1, :2])
if is_model:
    ax.plot(bresid, bvals, color=DK, lw=1.3)
    for lo, hi, col in [(90, 100, "#0053D6"), (70, 90, "#65CBF3"),
                        (50, 70, "#FFDB13"), (0, 50, "#FF7D45")]:
        ax.axhspan(lo, hi, color=col, alpha=0.18, zorder=0)
    ax.set_ylim(0, 100); ax.set_ylabel("pLDDT")
    ax.set_title(f"Model confidence — {100 * float((bvals < 70).mean()):.0f}% of residues "
                 f"below 70", fontsize=11)
elif np.ptp(bvals) < 1e-6:
    ax.text(0.5, 0.5, f"this file carries no B-factor values\n(all {bvals[0]:.2f})",
            transform=ax.transAxes, ha="center", va="center", fontsize=12, color=GREY)
    ax.set_yticks([]); ax.set_ylabel("Cα B-factor  (Å²)")
    ax.set_title("No per-atom uncertainty in this entry — older depositions sometimes "
                 "omit it", fontsize=11)
else:
    ax.plot(bresid, bvals, color=ACC, lw=1.3)
    ax.fill_between(bresid, bvals, color=ACC, alpha=0.15)
    k = int(np.argmax(bvals))
    ax.annotate(f"{res[k].get_resname()}{bresid[k]}  {bvals[k]:.0f} Å²",
                xy=(bresid[k], bvals[k]), xytext=(0.7, 0.85), textcoords="axes fraction",
                fontsize=10, color=ORA,
                arrowprops=dict(arrowstyle="-|>", color=ORA, lw=1.4))
    ax.set_ylabel("Cα B-factor  (Å²)")
    ax.set_title(f"Where the model is least certain — B-factor "
                 f"{bvals.min():.0f}–{bvals.max():.0f} Å²", fontsize=11)
ax.set_xlabel("residue number")
ax.spines[["top", "right"]].set_visible(False)

# panel D — burial histogram
ax = fig.add_subplot(gs[1, 2])
ax.hist(vals, bins=24, color=GREY, edgecolor="white")
ax.axvline(BURIED_CUTOFF, color=ACC, lw=2)
ax.axvline(EXPOSED_CUTOFF, color=ORA, lw=2)
ax.text(BURIED_CUTOFF + 0.02, ax.get_ylim()[1] * 0.92, "buried", color=ACC, fontsize=10)
ax.text(EXPOSED_CUTOFF + 0.02, ax.get_ylim()[1] * 0.92, "exposed", color=ORA, fontsize=10)
ax.set_xlabel("relative SASA"); ax.set_ylabel("residues")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Burial distribution", fontsize=11)

# panel E — burial along the chain, with the core marked
ax = fig.add_subplot(gs[2, :2])
sasa_resid = np.array([k[1] for k in keys])
o = np.argsort(sasa_resid)
ax.plot(sasa_resid[o], vals[o], color=DK, lw=1.0)
ax.fill_between(sasa_resid[o], vals[o], BURIED_CUTOFF,
                where=vals[o] < BURIED_CUTOFF, color=ACC, alpha=0.5, step="mid")
if iface:
    for r, n, _ in iface:
        ax.axvline(r, color=ORA, lw=0.8, alpha=0.35, zorder=0)
    ax.text(0.01, 0.92, "orange lines: interface residues", transform=ax.transAxes,
            color=ORA, fontsize=10)
ax.axhline(BURIED_CUTOFF, color=ACC, lw=1, ls="--")
ax.set_xlabel("residue number"); ax.set_ylabel("relative SASA")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Burial along the chain — the shaded troughs are your hydrophobic core",
             fontsize=11)

# panel F — comparison, or amino-acid composition of the core
ax = fig.add_subplot(gs[2, 2])
if cmp_rows:
    names = ["all\nresidues", "core\nonly"]
    ax.bar([0, 1], [c[1] for c in cmp_rows], color=[ORA, ACC], width=0.55)
    for x, c in enumerate(cmp_rows):
        ax.text(x, c[1], f"{c[1]:.2f} Å", ha="center", va="bottom", fontsize=11,
                fontweight="bold", color=DK)
    ax.set_xticks([0, 1]); ax.set_xticklabels(names, fontsize=10)
    ax.set_ylabel("RMSD  (Å)")
    ax.set_ylim(0, max(c[1] for c in cmp_rows) * 1.3)
    ax.set_title("Same fold, two numbers", fontsize=11)
else:
    comp = {}
    for k, v in zip(keys, vals):
        if v < BURIED_CUTOFF:
            comp[k[2]] = comp.get(k[2], 0) + 1
    items = sorted(comp.items(), key=lambda t: -t[1])[:10]
    cols = [ACC if n in HYDROPHOBIC else ORA for n, _ in items]
    ax.barh([n for n, _ in items][::-1], [c for _, c in items][::-1],
            color=cols[::-1])
    ax.set_xlabel("residues in the core")
    ax.set_title("What is in the core (blue = nonpolar)", fontsize=10.5)
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle(f"{label} chain {chain} — Homework 2 panels", fontsize=14, y=0.985,
             color=DK)
figname = f"panels_{label}.png"
fig.savefig(figname, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ write report
say("## Figures")
say()
say(f"![panels]({figname})")
say()
say("## Questions to answer in your own words")
say()
for q in [
    "Which secondary-structure class does your protein belong to (all-α, all-β, α/β, "
    "α+β), and which numbers above support that?",
    "Pick one glycine and one proline in your structure. What would happen if each were "
    "substituted, and why — use their φ/ψ values, not general statements.",
    "Name one buried and one exposed residue from the lists above. For each, say what you "
    "would expect a hydrophobic-to-polar substitution to do, and why they differ.",
    "If your file has more than one chain: are the interface residues of the isolated "
    "monomer buried or exposed? If it has one chain: is that the biological unit?",
    "Look at the B-factor (or pLDDT) profile. Where is the structure least reliable, and "
    "would you trust a structural argument about a variant in that region?",
    "Any Ramachandran outliers? For each, decide whether it is a real strained "
    "conformation or a modelling error, and say how you decided.",
]:
    say(textwrap.fill(f"{q}", 96, initial_indent="- ", subsequent_indent="  "))
say()

out = f"report_{label}.md"
open(out, "w", encoding="utf-8").write("\n".join(report) + "\n")
print()
print(f"written: {out} and {figname}")
