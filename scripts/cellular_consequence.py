#!/usr/bin/env python3
"""cellular_consequence.py — determine molecular consequence from variant data.

Input:  JSON from classify_variant.py + gene symbol.
Output: tab-separated tables with:
        - NMD-based consequence per transcript
        - gnomAD constraint metrics (pLI, LOF z-score)
        - ClinGen dosage sensitivity
        - Conclusion: haploinsufficiency / dominant-negative / gain-of-function

Usage:
    python3 cellular_consequence.py --gene PRUNE1 \
        --classify-json ../data/variant_c.521_2A_gt_G.json \
        --outdir .
"""

import argparse
import datetime
import json
import os
import sys

import requests

ACCESS = datetime.date.today().isoformat()
S = requests.Session()
S.headers.update({"Accept": "application/json"})


def get(url, **params):
    for _ in range(4):
        try:
            r = S.get(url, params=params, timeout=90)
            if r.status_code == 200:
                return r.json()
        except requests.RequestException:
            pass
    return None


# ────────────────────────────────────────────── gnomAD constraint
def fetch_gnomad_constraint(gene):
    """Fetch gnomAD v2.1 constraint metrics for a gene via GraphQL."""
    query = """
    {
      gene(gene_symbol: "%s", reference_genome: GRCh37) {
        gene_id
        symbol
        exac_constraint {
          exp_lof
          obs_lof
          mu_lof
          pli
          lof_z
          exp_syn
          obs_syn
          exp_mis
          obs_mis
        }
      }
    }
    """ % gene
    url = "https://gnomad.broadinstitute.org/api"
    try:
        r = S.post(url, json={"query": query, "variables": {}}, timeout=90)
        if r.status_code == 200:
            data = r.json()
            gene_data = data.get("data", {}).get("gene", {})
            if gene_data and gene_data.get("exac_constraint"):
                c = gene_data["exac_constraint"]
                return {
                    "gene_id": gene_data.get("gene_id"),
                    "symbol": gene_data.get("symbol"),
                    "exp_lof": c.get("exp_lof"),
                    "obs_lof": c.get("obs_lof"),
                    "pLI": c.get("pli"),
                    "lof_z": c.get("lof_z"),
                    "source": "gnomAD v2.1 exac_constraint (GRCh37)",
                }
    except requests.RequestException:
        pass
    return None


# ────────────────────────────────────────────── ClinGen / ClinVar
def fetch_clingen(gene):
    """Fetch ClinGen dosage sensitivity for a gene."""
    url = f"https://clinicalgenome.org/api/curation/gene/{gene}"
    try:
        r = S.get(url, timeout=90)
        if r.status_code == 200:
            data = r.json()
            out = {"source": "ClinGen"}
            hap = data.get("haploinsufficiency")
            if hap:
                out["haploinsufficiency_score"] = hap.get("score")
                out["haploinsufficiency_description"] = hap.get("description")
            trip = data.get("triplosensitivity")
            if trip:
                out["triplosensitivity_score"] = trip.get("score")
                out["triplosensitivity_description"] = trip.get("description")
            return out
    except requests.RequestException:
        pass
    return None


def fetch_clinvar_gene(gene):
    """Fetch ClinVar variant count for a gene."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "clinvar", "term": f"{gene}[gene]", "retmax": 1, "retmode": "json"}
    try:
        r = S.get(url, params=params, timeout=90)
        if r.status_code == 200:
            data = r.json()
            count = data.get("esearchresult", {}).get("count", "?")
            return {"clinvar_variant_count": int(count) if count != "?" else count}
    except requests.RequestException:
        pass
    return None


# ────────────────────────────────────────── Determine consequence
def determine_consequence(classify_data):
    """From classify_variant.py JSON, determine molecular consequence per transcript."""
    splice = classify_data.get("splice_analysis") or {}
    coding = classify_data.get("coding_analysis") or {}
    consequences = []

    if splice:
        for c in splice.get("splicing_consequences", []):
            nmd = c.get("nmd", {})
            verdict = nmd.get("verdict", "unknown")
            exon_label = c.get("skipped_exons", str(c.get("skipped_exon", "")))

            if verdict == "NMD_not_applicable":
                cons = {
                    "transcript": f"\u0394exon{exon_label}",
                    "nmd_verdict": verdict,
                    "molecular_consequence": "truncated_protein",
                    "detail": c.get("description", ""),
                    "protein_made": True,
                    "mechanism": "hypomorph_or_dominant_negative",
                }
                consequences.append(cons)
            elif verdict == "NMD_competent":
                cons = {
                    "transcript": f"\u0394exon{exon_label}",
                    "nmd_verdict": verdict,
                    "molecular_consequence": "null",
                    "detail": c.get("description", ""),
                    "protein_made": False,
                    "mechanism": "haploinsufficiency_contributor",
                }
                consequences.append(cons)

    if coding:
        vc = coding.get("variant_class")
        nmd = coding.get("nmd") or {}
        nmd_v = nmd.get("verdict", "unknown")
        if vc in ("nonsense", "frameshift"):
            if "ESCAPES" in nmd_v:
                consequences.append({
                    "transcript": coding.get("hgvs_c", "?"),
                    "nmd_verdict": nmd_v,
                    "molecular_consequence": "truncated_protein",
                    "protein_made": True,
                    "mechanism": "dominant_negative_or_hypomorph",
                })
            else:
                consequences.append({
                    "transcript": coding.get("hgvs_c", "?"),
                    "nmd_verdict": nmd_v,
                    "molecular_consequence": "null",
                    "protein_made": False,
                    "mechanism": "haploinsufficiency",
                })
        elif vc == "missense":
            consequences.append({
                "transcript": coding.get("hgvs_c", "?"),
                "nmd_verdict": "not_applicable",
                "molecular_consequence": "substituted_protein",
                "protein_made": True,
                "mechanism": "depends_on_residue_function",
            })
    return consequences


# ────────────────────────────────────────── Interpret constraint
def interpret_constraint(constraint):
    if not constraint:
        return {"verdict": "no_data"}
    pLI = constraint.get("pLI")
    result = {"pLI": pLI}
    if pLI is not None and pLI >= 0.9:
        result["haploinsufficiency_evidence"] = "strong"
        result["message"] = f"pLI = {pLI:.3f} — highly LOF-intolerant; haploinsufficiency plausible"
    elif pLI is not None and pLI >= 0.5:
        result["haploinsufficiency_evidence"] = "moderate"
        result["message"] = f"pLI = {pLI:.3f} — moderately LOF-intolerant"
    elif pLI is not None:
        result["haploinsufficiency_evidence"] = "weak"
        result["message"] = f"pLI = {pLI:.3f} — tolerant to LOF; haploinsufficiency unlikely for this gene"
    return result


# ───────────────────────────────────────────────────── Main
def main():
    parser = argparse.ArgumentParser(
        description="Determine cellular consequence of a variant.")
    parser.add_argument("--gene", default="PRUNE1")
    parser.add_argument("--classify-json", default=None)
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    print(f"# cellular_consequence.py — {args.gene}")
    print(f"# Accessed: {ACCESS}\n")

    classify_data = {}
    if args.classify_json and os.path.exists(args.classify_json):
        with open(args.classify_json) as f:
            classify_data = json.load(f)

    # 1. Molecular consequence
    print("== Molecular consequence (from NMD analysis)")
    cons = determine_consequence(classify_data)
    if cons:
        print("\t".join(["Transcript", "NMD verdict", "Consequence",
                         "Protein?", "Mechanism"]))
        for c in cons:
            print("\t".join([c.get("transcript", "?"), c.get("nmd_verdict", "?"),
                             c.get("molecular_consequence", "?"),
                             str(c.get("protein_made", "?")), c.get("mechanism", "?")]))
        print()
        made = [c for c in cons if c.get("protein_made")]
        null = [c for c in cons if not c.get("protein_made")]
        if null and made:
            print("  Mixed: in-frame product made, out-of-frame product degraded.\n"
                  "  Net effect depends on in-frame product function.\n")
        elif null and not made:
            print("  All transcripts degraded. Allele is functionally NULL.\n")
        elif made and not null:
            print("  All transcripts make protein.\n")

    # 2. gnomAD constraint
    print("== gnomAD constraint (exac)")
    constraint = fetch_gnomad_constraint(args.gene)
    if constraint:
        print(f"Symbol\t{constraint['symbol']}")
        print(f"Gene ID\t{constraint['gene_id']}")
        print(f"pLI\t{constraint['pLI']}")
        print(f"LOF z-score\t{constraint.get('lof_z', '?')}")
        print(f"Observed LOF\t{constraint['obs_lof']}")
        print(f"Expected LOF\t{constraint['exp_lof']:.2f}")
        print()
    else:
        print("  (not available)\n")

    # 3. Interpret constraint
    print("== Constraint interpretation")
    interp = interpret_constraint(constraint)
    if interp.get("message"):
        print(f"  {interp['message']}\n")

    # 4. ClinGen
    print("== ClinGen dosage sensitivity")
    clingen = fetch_clingen(args.gene)
    if clingen:
        for k, v in clingen.items():
            print(f"  {k}\t{v}")
    else:
        print("  (not available)")
    print()

    # 5. ClinVar
    print("== ClinVar")
    clinvar = fetch_clinvar_gene(args.gene)
    if clinvar:
        for k, v in clinvar.items():
            print(f"  {k}\t{v}")
    print()

    # 6. Summary
    print("== Summary")
    print(f"  Gene: {args.gene}")
    print(f"  Inheritance: autosomal recessive")
    print(f"  Allele state: homozygous")
    if constraint:
        pLI = constraint.get("pLI")
        if pLI is not None and pLI < 0.5:
            print(f"  pLI = {pLI:.3f}: gene is NOT LOF-constrained")
            print(f"  → Carriers (heterozygotes) tolerate haploinsufficiency — consistent with recessive inheritance")
    if cons:
        has_inframe = any(c.get("protein_made") for c in cons)
        has_null = any(not c.get("protein_made") for c in cons)
        if has_inframe and has_null:
            print(f"  Allele produces TWO transcripts:")
            print(f"    - In-frame \u0394exon5 (53 aa deletion) → truncated PRUNE1 protein")
            print(f"    - Out-of-frame \u0394exons4+5 → NMD, degraded")
            print(f"  Cellular consequence: partial loss-of-function (hypomorph) vs null")
            print(f"  The paper hypothesizes: in-frame product retains residual activity,")
            print(f"  producing a degenerative phenotype rather than purely developmental.")
    print()

    # 7. Save JSON
    result = {
        "accessed": ACCESS, "gene": args.gene,
        "molecular_consequences": cons,
        "gnomad_constraint": constraint,
        "constraint_interpretation": interp,
        "clingen": clingen, "clinvar": clinvar,
    }
    out_path = os.path.join(outdir, f"consequence_{args.gene}.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"# Data saved to {out_path}")


if __name__ == "__main__":
    main()