#!/usr/bin/env python3
"""hw3_variant.py — variant report from a mechanistic argument.

Structural bioinformatics in the aspect of medical genetics, Lecture 3.

Now generic via CLI args. Usage:
  python3 hw3_variant.py --gene PRUNE1 --uniprot Q86TP1 \
      --variant c.316G>A p.Asp106Asn \
      --variant c.1252G>T p.Glu418Ter \
      --hpo HP:0010864 HP:0002500 HP:0002059 HP:0002187 HP:0000252

WHAT IT DOES
    Takes a gene and coding-level variant(s) and works out, from live databases:
      1. the variant class;
      2. whether the transcript is NMD-competent (the 50-nt rule);
      3. where the affected residue sits in the folded protein;
      4. whether the HPO terms you chose actually exist.
    It writes a short report and a two-panel figure.

WHAT YOU NEED
    python 3.9+, requests, numpy, matplotlib.
    Network access to rest.ensembl.org, alphafold.ebi.ac.uk, www.ebi.ac.uk.
"""

import argparse
import datetime
import json
import os
import re
import sys
import urllib.parse

import numpy as np
import requests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ─────────────── Parse CLI ───────────────
parser = argparse.ArgumentParser(description="Variant report from coding HGVS")
parser.add_argument("--gene", default="PRUNE1", help="HGNC gene symbol")
parser.add_argument("--transcript", default=None, help="Ensembl transcript ID (default: canonical)")
parser.add_argument("--uniprot", default="Q86TP1", help="UniProt accession")
parser.add_argument("--variant", action="append", nargs=2, metavar=("C_HGVS", "P_HGVS"),
                    default=[], help="Variant pair: coding and protein HGVS (repeatable)")
parser.add_argument("--hpo", action="append", default=[], help="HPO term ID (repeatable)")
parser.add_argument("--outdir", default=".", help="Output directory")
args = parser.parse_args()

GENE = args.gene
TRANSCRIPT = args.transcript
UNIPROT = args.uniprot
VARIANTS = args.variant or []
HPO_TERMS = args.hpo or []
OUTDIR = args.outdir

# ─────────────── Everything below is unchanged from the original ───────────────

ACCESS = datetime.date.today().isoformat()
S = requests.Session()
S.headers.update({"Accept": "application/json"})

DK, ACC, ORA, GREY = "#0E2841", "#156082", "#E97132", "#595F69"
RED, GREEN, PL = "#B3261E", "#15803D", "#AFC7DA"
plt.rcParams.update({"font.family": ["DejaVu Sans"], "savefig.dpi": 200,
                     "savefig.bbox": "tight"})

NOTES = []          # anything that could not be determined goes here


def warn(msg):
    NOTES.append(msg)
    print("  ! " + msg)


def get(url, **params):
    for _ in range(4):
        try:
            r = S.get(url, params=params, timeout=90)
            if r.status_code == 200:
                return r.json()
        except requests.RequestException:
            pass
    return None


# ------------------------------------------------------------ 1. transcript
def transcript_geometry(gene, tx_id=None):
    """Exon structure and the NMD boundary, from Ensembl."""
    if tx_id is None:
        g = get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}", expand=1)
        if not g:
            sys.exit(f"Could not look up {gene} in Ensembl. Check the symbol and the network.")
        can = [t for t in g["Transcript"] if t.get("is_canonical") == 1]
        tx_id = (can or g["Transcript"])[0]["id"]
    j = get(f"https://rest.ensembl.org/lookup/id/{tx_id}", expand=1)
    if not j or "Exon" not in j:
        sys.exit(f"Could not fetch the exon structure of {tx_id}.")

    strand = j["strand"]
    ex = sorted(j["Exon"], key=lambda e: e["start"], reverse=(strand == -1))
    lens = [e["end"] - e["start"] + 1 for e in ex]
    cdna_len = sum(lens)

    def cdna(gpos):
        acc = 0
        for e, L in zip(ex, lens):
            if e["start"] <= gpos <= e["end"]:
                return acc + ((gpos - e["start"] + 1) if strand == 1
                              else (e["end"] - gpos + 1))
            acc += L
        return None

    tr = j["Translation"]
    a, b = cdna(tr["start"]), cdna(tr["end"])
    cds_start, cds_end = min(a, b), max(a, b)
    cds_len = cds_end - cds_start + 1
    last_junction_cdna = cdna_len - lens[-1]

    bounds, acc = [], 0
    for k, L in enumerate(lens, start=1):
        bounds.append((k, acc + 1 - cds_start + 1, acc + L - cds_start + 1))
        acc += L

    lj_c = last_junction_cdna - cds_start + 1
    return dict(
        transcript=tx_id, n_exons=len(ex), cdna_len=cdna_len, cds_len=cds_len,
        protein_len=cds_len // 3 - 1, utr5=cds_start - 1, utr3=cdna_len - cds_end,
        exon_lens=lens, exon_bounds_c=bounds,
        last_exon_coding=cds_end - last_junction_cdna,
        last_junction_c=lj_c, nmd_boundary_c=lj_c - 50,
        pct_nmd_competent=round(100 * (lj_c - 50) / cds_len, 1))


# --------------------------------------------------------- 2. variant class
CLASS_RULES = [
    (r"c\.-?\d+[+-]\d+[ACGT]>[ACGT]$", "intronic / splice-region"),
    (r"c\.\*\d+", "3-prime UTR"),
    (r"c\.-\d+", "5-prime UTR or promoter-proximal"),
    (r"c\.\d+_\d+del$", "deletion"),
    (r"c\.\d+del", "deletion"),
    (r"c\.\d+_\d+ins", "insertion"),
    (r"c\.\d+dup", "duplication"),
    (r"c\.\d+[ACGT]>[ACGT]$", "substitution"),
]


def classify(cvar, pvar, geo):
    """Variant class, affected codon and residue. Coarse, but explicit about it."""
    out = {"hgvs_c": cvar, "hgvs_p": pvar or "(not supplied)"}
    kind = next((lab for pat, lab in CLASS_RULES if re.match(pat, cvar)), "unrecognised")
    out["record_form"] = kind

    m = re.match(r"c\.(\d+)", cvar)
    pos = int(m.group(1)) if m else None
    out["cds_position"] = pos

    sp = re.match(r"c\.(\d+)([+-])(\d+)", cvar)
    if sp:
        off = int(sp.group(3))
        out["variant_class"] = ("canonical splice site (+/-1, +/-2)" if off <= 2
                                else "extended splice region / deep intronic")
        out["note"] = ("Canonical-site variants are treated as null variants; deeper "
                       "positions need a splice predictor (Lecture 3, slide 13).")
        return out

    if pos is None or not (1 <= pos <= geo["cds_len"]):
        out["variant_class"] = kind + " (outside the coding sequence)"
        return out

    out["codon"] = (pos - 1) // 3 + 1
    out["residue"] = out["codon"]

    pv = (pvar or "").lower()
    if kind == "substitution":
        if "ter" in pv or "*" in pv:
            out["variant_class"] = "nonsense"
        elif "=" in pv:
            out["variant_class"] = "synonymous"
        elif re.search(r"p\.\(?[a-z]{3}\d+[a-z]{3}", pv):
            out["variant_class"] = "missense"
        else:
            out["variant_class"] = "substitution (protein consequence not supplied)"
    elif kind in ("deletion", "insertion", "duplication"):
        span = re.findall(r"\d+", cvar)
        n = (int(span[1]) - int(span[0]) + 1) if len(span) > 1 else 1
        if "ins" in cvar:
            n = len(re.sub(r".*ins", "", cvar))
        out["variant_class"] = ("in-frame " + kind if n % 3 == 0
                                else "frameshift (" + kind + ")")
        out["bases_changed"] = n
    else:
        out["variant_class"] = kind

    for k, lo, hi in geo["exon_bounds_c"]:
        if lo <= pos <= hi:
            out["exon"] = k
            break
    return out


# ------------------------------------------------------------------- 3. NMD
def nmd_call(v, geo):
    """Apply the 50-nt rule. Only meaningful for variants that create a PTC."""
    makes_ptc = any(k in v.get("variant_class", "")
                    for k in ("nonsense", "frameshift", "canonical splice"))
    if not makes_ptc:
        return dict(applies=False,
                    verdict="not applicable -- no premature stop codon",
                    reasoning="NMD is a question about premature termination. This "
                              "variant class does not create one, so the protein-level "
                              "question is about the residue, not the transcript.")
    pos = v.get("cds_position")
    if pos is None:
        return dict(applies=True, verdict="cannot be determined",
                    reasoning="No coding position could be parsed from the record.")
    b, lj = geo["nmd_boundary_c"], geo["last_junction_c"]
    if pos > lj:
        verdict = "ESCAPES NMD -- the stop falls in the last exon"
        why = (f"The last exon-exon junction is at c.{lj}; the variant is at c.{pos}, "
               f"downstream of it. No junction complexes remain after the stop, so a "
               f"truncated protein IS made. Consider a dominant-negative mechanism, and "
               f"downgrade PVS1.")
    elif pos > b:
        verdict = "BORDERLINE -- inside the 50-nt window"
        why = (f"The variant is at c.{pos}, between the boundary c.{b} and the last "
               f"junction c.{lj}. The terminating ribosome may cover the last junction "
               f"complex. Treat NMD as uncertain and downgrade PVS1 accordingly.")
    else:
        verdict = "NMD-COMPETENT -- the transcript is degraded"
        why = (f"The variant is at c.{pos}, which is {b - pos} nt upstream of the "
               f"boundary c.{b} (itself 50 nt before the last junction at c.{lj}). "
               f"Junction complexes remain downstream of the stop, so the transcript is "
               f"degraded and this allele makes essentially no protein.")
    return dict(applies=True, verdict=verdict, reasoning=why,
                position_pct_of_cds=round(100.0 * pos / geo["cds_len"], 1),
                protein_removed_pct=(round(100 - 100.0 * v["residue"] / geo["protein_len"], 1)
                                     if v.get("residue") else None))


# ---------------------------------------------------- 4. structural context
def structure_context(uniprot, residues):
    """Fetch the AlphaFold model and describe each residue's local environment."""
    if not uniprot:
        warn("No UniProt accession given, so the structural question was skipped.")
        return None
    meta = get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}")
    if not meta:
        warn(f"No AlphaFold model found for {uniprot}; structural question skipped.")
        return None
    path = os.path.join(OUTDIR, f"AF-{uniprot}.pdb")
    if not os.path.exists(path):
        r = S.get(meta[0]["pdbUrl"], timeout=180)
        if r.status_code != 200:
            warn("Could not download the AlphaFold model.")
            return None
        open(path, "wb").write(r.content)

    ca, plddt, resid, three = [], [], [], {}
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            ca.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            plddt.append(float(line[60:66]))
            rid = int(line[22:26])
            resid.append(rid)
            three[rid] = line[17:20].strip()
    ca, plddt, resid = np.asarray(ca), np.asarray(plddt), np.asarray(resid)
    nb = np.array([(np.linalg.norm(ca - p, axis=1) < 10).sum() - 1 for p in ca])

    out = dict(uniprot=uniprot, model=meta[0]["entryId"], length=len(ca),
               mean_plddt=round(float(plddt.mean()), 1),
               nb_median=int(np.median(nb)), residues={},
               _plddt=plddt, _resid=resid)
    known = set(resid.tolist())
    for r in residues:
        if r is None or r not in known:
            continue
        i = int(np.where(resid == r)[0][0])
        pct = float((nb < nb[i]).mean() * 100)
        out["residues"][r] = dict(
            aa=three.get(r, "?"), plddt=round(float(plddt[i]), 1),
            neighbours=int(nb[i]), burial_percentile=round(pct),
            environment=("buried core" if pct >= 75 else
                         "partly buried" if pct >= 40 else "surface-exposed"),
            model_confidence=("confident -- local geometry can be argued from"
                              if plddt[i] >= 70 else
                              "LOW -- do not argue local geometry from this model"))
    return out


# ------------------------------------------------------------------- 5. HPO
def check_hpo(terms):
    """Confirm each identifier exists, and report its label and parents."""
    out = []
    for t in terms:
        iri = urllib.parse.quote_plus(urllib.parse.quote_plus(
            f"http://purl.obolibrary.org/obo/{t.strip().replace(':', '_')}"))
        j = get(f"https://www.ebi.ac.uk/ols4/api/ontologies/hp/terms/{iri}")
        if not j or not j.get("label"):
            warn(f"{t} did not resolve in the current HPO release -- check the identifier.")
            out.append(dict(id=t, label=None, parents=[]))
            continue
        pj = get(f"https://www.ebi.ac.uk/ols4/api/ontologies/hp/terms/{iri}/parents")
        parents = [(p["short_form"].replace("_", ":"), p["label"])
                   for p in (pj or {}).get("_embedded", {}).get("terms", [])]
        out.append(dict(id=t, label=j["label"], parents=parents,
                        definition=(j.get("description") or [""])[0]))
    return out


# ----------------------------------------------------------------- 6. figure
def figure(geo, analyses, struct, gene, path):
    n = 2 if struct else 1
    fig, axes = plt.subplots(n, 1, figsize=(11.5, 3.1 * n + 0.5), squeeze=False)
    ax = axes[0][0]

    lens, total = geo["exon_lens"], geo["cdna_len"]
    sc = 100.0 / total
    cds_s, cds_e = geo["utr5"] + 1, geo["utr5"] + geo["cds_len"]
    lj = total - lens[-1]
    nbd = lj - 50

    ax.add_patch(plt.Rectangle((nbd * sc, 0.355), max((cds_e - nbd) * sc, 0.2), 0.31,
                               facecolor=GREEN, alpha=0.30, lw=0, zorder=1))
    acc = 0
    for L in lens:
        # colour the coding and UTR parts of each exon separately, so a final
        # exon that is half 3' UTR is drawn as such
        a0, a1 = acc, acc + L
        for lo, hi, col in ((a0, min(a1, cds_s - 1), PL),
                            (max(a0, cds_s - 1), min(a1, cds_e), ORA),
                            (max(a0, cds_e), a1, PL)):
            if hi > lo:
                ax.add_patch(plt.Rectangle((lo * sc, 0.38), max((hi - lo) * sc, 0.04),
                                           0.24, facecolor=col, edgecolor="none",
                                           lw=0, zorder=2))
        ax.plot([a1 * sc] * 2, [0.38, 0.62], color="white", lw=1.2, zorder=3)
        acc = a1
    ax.plot([nbd * sc] * 2, [0.33, 0.72], color=RED, lw=1.8, zorder=4)
    ax.text(nbd * sc, 0.74, f"NMD boundary c.{geo['nmd_boundary_c']}", color=RED,
            fontsize=9, ha="center", va="bottom", fontweight="bold")
    ax.plot([lj * sc] * 2, [0.34, 0.66], color=DK, lw=1.1, ls=":", zorder=4)
    ax.text(101, 0.50, "orange coding\nblue UTR\ngreen NMD escape", fontsize=7.5,
            color=GREY, ha="left", va="center")

    for k, (v, nm) in enumerate(analyses):
        if v.get("cds_position") is None:
            continue
        x = (cds_s - 1 + v["cds_position"]) * sc
        y = 0.20 - 0.10 * k
        ax.annotate(f"{v['hgvs_c']}  ({v.get('variant_class', '?')})",
                    xy=(x, 0.36), xytext=(x, y), ha="center", va="top", fontsize=9,
                    color=DK, fontweight="bold",
                    arrowprops=dict(arrowstyle="-|>", color=DK, lw=1.3))
        ax.text(x, y - 0.052, nm["verdict"].split(" -- ")[0], ha="center", va="top",
                fontsize=8.5,
                color=GREEN if "ESCAPE" in nm["verdict"] else
                      RED if "COMPETENT" in nm["verdict"] else ORA)
    ax.set_xlim(-2, 113)
    ax.set_ylim(-0.16, 0.92)
    ax.axis("off")
    ax.set_title(f"{gene}  {geo['transcript']} -- {geo['n_exons']} exons, "
                 f"{geo['pct_nmd_competent']} % of the CDS is NMD-competent",
                 fontsize=11.5, color=DK, loc="left")

    if struct:
        a2 = axes[1][0]
        a2.fill_between(struct["_resid"], 0, struct["_plddt"], color=PL, lw=0)
        a2.plot(struct["_resid"], struct["_plddt"], color=ACC, lw=0.6)
        a2.axhline(70, color=GREY, ls="--", lw=1.0)
        a2.text(struct["length"], 72, "pLDDT 70", fontsize=8, color=GREY,
                ha="right", va="bottom")
        for r, d in struct["residues"].items():
            a2.axvline(r, color=RED, lw=1.6)
            a2.annotate(f"{d['aa']}{r}\n{d['environment']}", xy=(r, d["plddt"]),
                        xytext=(r, 110), ha="center", va="bottom", fontsize=8.5,
                        color=RED, fontweight="bold",
                        arrowprops=dict(arrowstyle="-", color=RED, lw=1.0))
        a2.set_ylim(0, 132)
        a2.set_xlim(1, struct["length"])
        a2.set_yticks([0, 25, 50, 75, 100])
        a2.set_xlabel("residue")
        a2.set_ylabel("pLDDT")
        a2.spines[["top", "right"]].set_visible(False)
        a2.set_title(f"{struct['model']} -- model confidence along the sequence "
                     f"(mean {struct['mean_plddt']})", fontsize=11.5, color=DK,
                     loc="left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


# ------------------------------------------------------------------ 7. report
def report(gene, geo, analyses, struct, hpo, path):
    L = []
    A = L.append
    A(f"# Homework 3 — {gene}\n")
    A(f"Generated {ACCESS} by `hw3_variant.py`. Every number below was fetched or "
      f"computed at run time.\n")

    A("## 1. The transcript\n")
    A("| | |")
    A("|---|---|")
    A(f"| gene | *{gene}* |")
    A(f"| transcript | `{geo['transcript']}` (Ensembl canonical) |")
    A(f"| exons | {geo['n_exons']} |")
    A(f"| coding sequence | {geo['cds_len']} nt = {geo['protein_len']} amino acids |")
    A(f"| 5' UTR / 3' UTR | {geo['utr5']} nt / {geo['utr3']} nt |")
    A(f"| last exon-exon junction | `c.{geo['last_junction_c']}` |")
    A(f"| **NMD boundary** | **`c.{geo['nmd_boundary_c']}`** (50 nt before it) |")
    A(f"| NMD-competent share of the CDS | **{geo['pct_nmd_competent']} %** |")
    A(f"\n*Source: Ensembl REST, GRCh38 — accessed {ACCESS}.*\n")
    A("> Sanity check to do yourself: the coding length divided by three, minus the "
      "stop codon, must equal the known protein length. If it does not, something is "
      "wrong with the transcript you chose.\n")

    A("## 2. The variants\n")
    for v, nm in analyses:
        A(f"### {v['hgvs_c']}  {v['hgvs_p']}\n")
        A(f"- record form: {v['record_form']}")
        A(f"- variant class: **{v.get('variant_class', '?')}**")
        if v.get("exon"):
            A(f"- falls in exon {v['exon']} of {geo['n_exons']}")
        if v.get("cds_position"):
            A(f"- coding position `c.{v['cds_position']}`"
              + (f", codon {v['codon']}, residue {v['residue']}" if v.get("codon") else ""))
        A(f"- **NMD: {nm['verdict']}**")
        A(f"- reasoning: {nm['reasoning']}")
        if nm.get("protein_removed_pct") is not None:
            A(f"- a stop here removes about {nm['protein_removed_pct']} % of the protein "
              f"(the PVS1 tree asks whether this exceeds ~10 %)")
        A("")

    A("## 3. Where the residue sits in the folded protein\n")
    if struct:
        A(f"Model `{struct['model']}`, {struct['length']} residues, mean pLDDT "
          f"{struct['mean_plddt']}. Median Ca neighbour count within 10 A across the "
          f"chain: {struct['nb_median']}.\n")
        if struct["residues"]:
            A("| residue | pLDDT | neighbours <10 A | burial percentile | environment | "
              "argue from the model? |")
            A("|---|---|---|---|---|---|")
            for r, d in struct["residues"].items():
                A(f"| {d['aa']}{r} | {d['plddt']} | {d['neighbours']} | "
                  f"{d['burial_percentile']}th | {d['environment']} | "
                  f"{d['model_confidence']} |")
        else:
            A("None of the affected residues could be located in the model.")
        A(f"\n*Source: AlphaFold DB — accessed {ACCESS}.*\n")
        A("> For a nonsense or frameshift variant this section is of secondary interest: "
          "if the transcript is degraded there is no protein to localise in. It matters "
          "most for missense and in-frame variants.\n")
    else:
        A("Not determined — see the notes at the end.\n")

    A("## 4. The phenotype, as identifiers\n")
    if hpo:
        A("| term | label | parents |")
        A("|---|---|---|")
        for h in hpo:
            par = "; ".join(f"{i} {l}" for i, l in h["parents"][:2]) or "—"
            A(f"| `{h['id']}` | {h['label'] or '**did not resolve**'} | {par} |")
        A(f"\n*Source: HPO via EBI OLS4 — accessed {ACCESS}. State the release you used "
          f"in your write-up.*\n")
        A("> Specific and true beats many. A profile padded with vague ancestors carries "
          "less information than three precise terms.\n")
    else:
        A("No HPO terms were supplied.\n")

    A("## 5. What you have to write\n")
    A("The script answered determinants 1, 2, 3 and part of 6. You write the argument:\n")
    A("1. Which **mechanism** do these variants imply — haploinsufficiency, "
      "dominant-negative, gain of function, or recessive loss of function? Justify it "
      "from the NMD call, not from the variant class alone.")
    A("2. What **inheritance pattern** follows from that mechanism, and does the gene's "
      "constraint in gnomAD agree?")
    A("3. Does the phenotype you annotated in section 4 fit the tissue in which this "
      "gene is **expressed**?")
    A("4. State the **PVS1 strength** you would assign, and which branch of the tree "
      "decided it.\n")

    if NOTES:
        A("## Notes and failures\n")
        for n in NOTES:
            A(f"- {n}")
        A("")
    open(path, "w", encoding="utf-8").write("\n".join(L))
    return path


# -------------------------------------------------------------------- driver
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print(f"Homework 3 — {GENE}")
    print("  fetching the transcript ...")
    geo = transcript_geometry(GENE, TRANSCRIPT)
    print(f"    {geo['transcript']}: {geo['n_exons']} exons, {geo['cds_len']} nt CDS "
          f"= {geo['protein_len']} aa")
    print(f"    NMD boundary c.{geo['nmd_boundary_c']}  "
          f"({geo['pct_nmd_competent']} % of the CDS is NMD-competent)")

    analyses = []
    for cv, pv in VARIANTS:
        v = classify(cv, pv, geo)
        nm = nmd_call(v, geo)
        analyses.append((v, nm))
        print(f"  {cv:14s} {v.get('variant_class', '?'):22s} -> {nm['verdict']}")

    print("  fetching the structure ...")
    struct = structure_context(UNIPROT, [v.get("residue") for v, _ in analyses])
    print("  checking the HPO terms ...")
    hpo = check_hpo(HPO_TERMS)

    png = figure(geo, analyses, struct, GENE,
                 os.path.join(OUTDIR, f"panels_{GENE}.png"))
    md = report(GENE, geo, analyses, struct, hpo,
                os.path.join(OUTDIR, f"report_{GENE}.md"))
    json.dump({"gene": GENE,
               "geometry": {k: v for k, v in geo.items() if not k.startswith("_")},
               "variants": [{"variant": v, "nmd": n} for v, n in analyses],
               "hpo": hpo, "accessed": ACCESS},
              open(os.path.join(OUTDIR, f"result_{GENE}.json"), "w"), indent=1)
    print(f"\n  wrote {md}\n        {png}\n        result_{GENE}.json")
    if NOTES:
        print("\n  things that did not work (they are in the report too):")
        for n in NOTES:
            print("   -", n)


if __name__ == "__main__":
    main()
