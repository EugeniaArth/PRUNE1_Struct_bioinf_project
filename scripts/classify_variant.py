#!/usr/bin/env python3
"""classify_variant.py — generic variant classifier.

Input: gene symbol + coding HGVS notation.
Output: tab-separated tables with Exon geometry, variant classification,
        splicing consequences (if applicable), residue context (if applicable).

No interpretation, no paper-specific text — just data fetched from public APIs
or computed from transcript/exon geometry.

Usage:
    python3 classify_variant.py [--gene GENE] [--transcript TX_ID] [--variant C_HGVS]

Example:
    python3 classify_variant.py --gene PRUNE1 --variant c.521-2A>G
    python3 classify_variant.py --gene PRUNE1 --variant c.316G>A --protein p.Asp106Asn
"""

import argparse
import datetime
import json
import os
import re
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


# ──────────────────────────────────────────────────── Transcript geometry
def fetch_transcript_geometry(gene, tx_id=None):
    """Fetch exon geometry from Ensembl REST."""
    if tx_id is None:
        g = get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}", expand=1)
        if not g:
            sys.exit(f"Could not look up {gene} in Ensembl.")
        can = [t for t in g["Transcript"] if t.get("is_canonical") == 1]
        tx_id = (can or g["Transcript"])[0]["id"]
    j = get(f"https://rest.ensembl.org/lookup/id/{tx_id}", expand=1)
    if not j or "Exon" not in j:
        sys.exit(f"Could not fetch exon structure for {tx_id}.")

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

    bounds = []
    acc = 0
    for k, L in enumerate(lens, start=1):
        start_c = acc + 1 - cds_start + 1
        end_c = acc + L - cds_start + 1
        bounds.append({"exon": k, "start_c": start_c, "end_c": end_c,
                        "len_nt": end_c - start_c + 1})
        acc += L

    lj_c = last_junction_cdna - cds_start + 1
    nmd_boundary = lj_c - 50

    return {
        "gene": gene,
        "transcript_id": tx_id,
        "strand": strand,
        "n_exons": len(ex),
        "cdna_len": cdna_len,
        "cds_len": cds_len,
        "cds_start": cds_start,
        "cds_end": cds_end,
        "protein_len": cds_len // 3 - 1,
        "utr5_len": cds_start - 1,
        "utr3_len": cdna_len - cds_end,
        "last_junction_c": lj_c,
        "nmd_boundary_c": nmd_boundary,
        "exons": bounds,
    }


# ────────────────────────────────────────────── Splice acceptor analysis
SPLICE_PAT = re.compile(r"c\.(\d+)([+-])(\d+)([ACGT])>([ACGT])$")

def analyse_splice_variant(hgvs_c, geo):
    """Analyse a splice-site variant. Return structured dict."""
    m = SPLICE_PAT.match(hgvs_c)
    if not m:
        return None

    ref_pos, sign, offset, ref_base, alt_base = m.groups()
    ref_pos, offset = int(ref_pos), int(offset)

    # Determine site type
    if sign == "-":
        site_type = "acceptor"
    else:
        site_type = "donor"

    if offset <= 2:
        site_class = "canonical_splice_site"
    else:
        site_class = "extended_splice_region"

    # Find affected exon and intron
    affected_exon = None
    affected_intron = None
    if sign == "-":
        for ex in geo["exons"]:
            if ex["start_c"] == ref_pos:
                affected_exon = ex["exon"]
                affected_intron = ex["exon"] - 1
                break
    else:
        for ex in geo["exons"]:
            if ex["start_c"] <= ref_pos <= ex["end_c"]:
                affected_exon = ex["exon"]
                affected_intron = ex["exon"]
                break

    # Compute exon-skipping consequences
    consequences = []
    if affected_exon is not None:
        # Single exon skipping
        ex = geo["exons"][affected_exon - 1]
        skip_len = ex["end_c"] - ex["start_c"] + 1
        in_frame = (skip_len % 3 == 0)
        conseq = {
            "type": "single_exon_skip",
            "skipped_exon": affected_exon,
            "skipped_region_c": f"c.{ex['start_c']}–c.{ex['end_c']}",
            "skipped_nt": skip_len,
            "reading_frame_maintained": in_frame,
            "codons_removed": skip_len // 3 if in_frame else None,
        }
        # NMD analysis for this consequence
        nmd = {
            "applies": not in_frame,
            "nmd_boundary_c": geo["nmd_boundary_c"],
            "last_junction_c": geo["last_junction_c"],
        }
        if in_frame:
            nmd["verdict"] = "NMD_not_applicable"
            nmd["reason"] = "in-frame deletion; no premature stop created"
            conseq["description"] = f"in-frame deletion of {skip_len // 3} codons"
        else:
            nmd["verdict"] = "NMD_competent"
            nmd["reason"] = "frameshift at c.{} creates premature stop well before NMD boundary (c.{})".format(
                ex["start_c"], geo["nmd_boundary_c"])
            conseq["description"] = (f"frameshift of {skip_len % 3} nt, "
                                     f"premature stop likely")
        conseq["nmd"] = nmd
        consequences.append(conseq)

    # Multi-exon skipping: current + upstream exon (if acceptor)
    if sign == "-" and affected_exon is not None and affected_exon > 1:
        ex_a = geo["exons"][affected_exon - 1]
        ex_b = geo["exons"][affected_exon - 2]
        combined_skip_len = (ex_a["end_c"] - ex_a["start_c"] + 1
                             + ex_b["end_c"] - ex_b["start_c"] + 1)
        in_frame = (combined_skip_len % 3 == 0)
        conseq2 = {
            "type": "multi_exon_skip",
            "skipped_exons": f"{ex_b['exon']}+{ex_a['exon']}",
            "skipped_region_c": f"c.{ex_b['start_c']}–c.{ex_a['end_c']}",
            "skipped_nt": combined_skip_len,
            "reading_frame_maintained": in_frame,
        }
        if in_frame:
            conseq2["description"] = (f"in-frame deletion of {combined_skip_len // 3} codons")
        else:
            conseq2["description"] = (f"frameshift of {combined_skip_len % 3} nt, "
                                      f"premature stop likely")
        # NMD
        nmd = {
            "applies": not in_frame,
            "nmd_boundary_c": geo["nmd_boundary_c"],
            "last_junction_c": geo["last_junction_c"],
        }
        if in_frame:
            nmd["verdict"] = "NMD_not_applicable"
            nmd["reason"] = "in-frame deletion; no premature stop created"
        else:
            nmd["verdict"] = "NMD_competent"
            nmd["reason"] = "frameshift at c.{} creates premature stop well before NMD boundary (c.{})".format(
                ex_b["start_c"], geo["nmd_boundary_c"])
        conseq2["nmd"] = nmd
        consequences.append(conseq2)

    return {
        "hgvs_c": hgvs_c,
        "site_type": site_type,
        "site_class": site_class,
        "ref_pos_c": ref_pos,
        "intron_offset": offset,
        "ref_base": ref_base,
        "alt_base": alt_base,
        "affected_exon": affected_exon,
        "affected_intron": affected_intron,
        "splicing_consequences": consequences,
    }


# ────────────────────────────────────────────── Coding variant analysis
CODON_PAT = re.compile(r"c\.(\d+)([ACGT]>[ACGT])?")
NMD_BOUNDARY_NT = 50

def analyse_coding_variant(hgvs_c, hgvs_p, geo):
    """Analyse a coding (missense/nonsense/synonymous) variant."""
    m = re.match(r"c\.(\d+)", hgvs_c)
    if not m:
        return None
    pos = int(m.group(1))
    if not (1 <= pos <= geo["cds_len"]):
        return None

    codon = (pos - 1) // 3 + 1
    residue = codon

    # Determine variant class from hgvs_p if provided
    vclass = None
    if hgvs_p:
        pl = hgvs_p.lower()
        if "ter" in pl or "*" in pl:
            vclass = "nonsense"
        elif "=" in pl or "syn" in pl:
            vclass = "synonymous"
        elif re.search(r"[a-z]{3}\d+[a-z]{3}", pl):
            vclass = "missense"
        elif "fs" in pl or "frameshift" in pl:
            vclass = "frameshift"
        elif "del" in pl:
            vclass = "deletion"
        elif "ins" in pl:
            vclass = "insertion"
        elif "ext" in pl:
            vclass = "extension"
        else:
            vclass = "unclassified_coding"

    # Find which exon
    affected_exon = None
    for ex in geo["exons"]:
        if ex["start_c"] <= pos <= ex["end_c"]:
            affected_exon = ex["exon"]
            break

    # NMD analysis (for nonsense/frameshift only)
    nmd = None
    if vclass in ("nonsense", "frameshift"):
        stop_in_last_exon = (pos > geo["last_junction_c"])
        nmd_competent = not stop_in_last_exon
        nmd = {
            "premature_stop_c": pos,
            "last_junction_c": geo["last_junction_c"],
            "nmd_boundary_c": geo["nmd_boundary_c"],
            "stop_downstream_of_last_junction": stop_in_last_exon,
            "nmd_competent": nmd_competent,
        }
        if stop_in_last_exon:
            nmd["verdict"] = "ESCAPES_NMD"
            nmd["protein_removed_pct"] = round(
                (1 - pos / geo["cds_len"]) * 100, 1)
        else:
            nmd["verdict"] = "NMD_TARGET"

    return {
        "hgvs_c": hgvs_c,
        "hgvs_p": hgvs_p,
        "variant_class": vclass,
        "cds_position": pos,
        "codon": codon,
        "residue": residue,
        "affected_exon": affected_exon,
        "nmd": nmd,
    }


# ─────────────────────────────────────────────────────────────── Main
def main():
    parser = argparse.ArgumentParser(
        description="Classify a human gene variant (splice or coding).")
    parser.add_argument("--gene", default="PRUNE1",
                        help="HGNC gene symbol (default: PRUNE1)")
    parser.add_argument("--transcript", default=None,
                        help="Ensembl transcript ID (default: canonical)")
    parser.add_argument("--variant", default="c.521-2A>G",
                        help="Coding HGVS (default: c.521-2A>G)")
    parser.add_argument("--protein", default=None,
                        help="Protein HGVS (optional, helps variant class)")
    parser.add_argument("--outdir", default=".",
                        help="Output directory (default: .)")
    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    print(f"# classify_variant.py — {args.gene} {args.variant}")
    print(f"# Accessed: {ACCESS}")
    print(f"# Source: Ensembl REST API (rest.ensembl.org)")
    print()

    # 1. Transcript geometry
    geo = fetch_transcript_geometry(args.gene, args.transcript)

    print("== Transcript geometry")
    print(f"Gene\t{geo['gene']}")
    print(f"Transcript\t{geo['transcript_id']}")
    print(f"Exons\t{geo['n_exons']}")
    print(f"CDS length (nt)\t{geo['cds_len']}")
    print(f"Protein length (aa)\t{geo['protein_len']}")
    print(f"Last exon-exon junction (c.)\t{geo['last_junction_c']}")
    print(f"NMD boundary (c.)\t{geo['nmd_boundary_c']}")
    print(f"5' UTR (nt)\t{geo['utr5_len']}")
    print(f"3' UTR (nt)\t{geo['utr3_len']}")
    print()

    print("== Exon layout")
    print("\t".join(["Exon", "Start (c.)", "End (c.)", "Length (nt)"]))
    for ex in geo["exons"]:
        print("\t".join([
            str(ex["exon"]),
            str(ex["start_c"]),
            str(ex["end_c"]),
            str(ex["len_nt"]),
        ]))
    print()

    # 2. Detect variant type
    splice_result = analyse_splice_variant(args.variant, geo)

    if splice_result:
        # Splice variant
        print("== Variant classification (splice)")
        r = splice_result
        print(f"HGVS (c.)\t{r['hgvs_c']}")
        print(f"Site type\t{r['site_type']}")
        print(f"Site class\t{r['site_class']}")
        print(f"Intron offset (nt)\t{r['intron_offset']}")
        print(f"Reference base\t{r['ref_base']}")
        print(f"Alternate base\t{r['alt_base']}")
        print(f"Affected intron\t{r['affected_intron']}")
        print(f"Affected exon\t{r['affected_exon']}")
        print()

        if r["splicing_consequences"]:
            print("== Splicing consequences (if acceptor disrupted)")
            print("\t".join(["Type", "Skipped exon(s)", "Region (c.)",
                             "Skipped (nt)", "Reading frame", "NMD", "Description"]))
            for c in r["splicing_consequences"]:
                nmd_v = c.get("nmd", {}).get("verdict", "-")
                print("\t".join([
                    c["type"],
                    str(c.get("skipped_exons", str(c.get("skipped_exon", "")))),
                    c["skipped_region_c"],
                    str(c["skipped_nt"]),
                    "maintained" if c["reading_frame_maintained"] else "disrupted",
                    nmd_v,
                    c["description"],
                ]))
            print()

    # 3. Also try as coding variant (for completeness)
    coding_result = analyse_coding_variant(args.variant, args.protein, geo)

    if coding_result:
        print("== Variant classification (coding)")
        r = coding_result
        print(f"HGVS (c.)\t{r['hgvs_c']}")
        print(f"HGVS (p.)\t{r['hgvs_p'] or '(not provided)'}")
        print(f"Variant class\t{r['variant_class']}")
        print(f"CDS position\t{r['cds_position']}")
        print(f"Codon\t{r['codon']}")
        print(f"Residue\t{r['residue']}")
        print(f"Affected exon\t{r['affected_exon']}")
        print()

        if r["nmd"]:
            print("== NMD analysis")
            n = r["nmd"]
            print(f"Premature stop at (c.)\t{n['premature_stop_c']}")
            print(f"Last junction at (c.)\t{n['last_junction_c']}")
            print(f"NMD boundary (c.)\t{n['nmd_boundary_c']}")
            print(f"Stop in last exon\t{n['stop_downstream_of_last_junction']}")
            print(f"NMD competent\t{n['nmd_competent']}")
            print(f"Verdict\t{n['verdict']}")
            if "protein_removed_pct" in n:
                print(f"Protein removed (%)\t{n['protein_removed_pct']}")
            print()

    # 4. Save raw JSON
    result = {
        "accessed": ACCESS,
        "gene": args.gene,
        "transcript": args.transcript or geo["transcript_id"],
        "variant": args.variant,
        "transcript_geometry": geo,
        "splice_analysis": splice_result,
        "coding_analysis": coding_result,
    }
    safe_name = args.variant.replace(">", "_gt_").replace("-", "_").replace("*", "X")
    out_path = os.path.join(outdir, f"variant_{safe_name}.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"# Data saved to {out_path}")


if __name__ == "__main__":
    main()