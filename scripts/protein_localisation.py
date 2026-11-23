#!/usr/bin/env python3
"""protein_localisation.py — map variant-affected residues onto protein domains.

Input:  gene symbol + UniProt accession + variant HGVS + JSON from classify_variant.py
Output: tab-separated tables with:
        - domain annotations (fetched from InterPro)
        - residue-level pLDDT (from AlphaFold confidence JSON)
        - overlap of affected region with domains / linkers / active sites

Usage:
    python3 protein_localisation.py \\
        --gene PRUNE1 --uniprot Q86TP1 \\
        --variant c.521-2A>G \\
        --classify-json ../Working_dir/variant_c.521_2A_gt_G.json \\
        --plddt-json ../Working_dir/data/AF-Q86TP1-F1-confidence_v6.json \\
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


# ──────────────────────────────────────────────── Fetch domain annotations
def fetch_domains(uniprot):
    """Fetch InterPro domain annotations for a UniProt accession.

    Returns a list of dicts: {accession, name, type, start, end}.
    """
    url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{uniprot}/"
    data = get(url, page_size=100)
    if not data:
        return []

    domains = []
    for result in data.get("results", []):
        meta = result.get("metadata", {})
        acc = meta.get("accession", "?")
        name = meta.get("name", "?")
        domain_type = meta.get("type", "?")
        integrated = meta.get("integrated")

        for prot in result.get("proteins", []):
            for loc in prot.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    start = frag.get("start")
                    end = frag.get("end")
                    if start and end:
                        domains.append({
                            "accession": acc,
                            "name": name,
                            "type": domain_type,
                            "integrated": integrated,
                            "start": start,
                            "end": end,
                            "source": "InterPro",
                        })
    return domains


# ────────────────────────────────────────────── Load AlphaFold confidence
def fetch_alphafold_plddt(uniprot):
    """Fetch AlphaFold confidence JSON from the EBI API.

    Returns {residue: score} dict, or None if unavailable.
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    try:
        meta = S.get(url, timeout=90).json()
        if meta and len(meta) > 0:
            conf_url = meta[0].get("confidenceUrl")
            if conf_url:
                r = S.get(conf_url, timeout=90)
                if r.status_code == 200:
                    return parse_plddt(r.json())
    except Exception:
        pass
    # Fallback: try direct file URL
    direct = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-confidence_v6.json"
    try:
        r = S.get(direct, timeout=90)
        if r.status_code == 200:
            return parse_plddt(r.json())
    except Exception:
        pass
    return None


def parse_plddt(d):
    """Parse pLDDT from AlphaFold JSON into {residue: score}."""
    if isinstance(d, list) and len(d) > 0 and "residueNumber" in d[0]:
        return {e["residueNumber"]: e["confidenceScore"] for e in d}
    if "residueNumber" in d and "confidenceScore" in d:
        return dict(zip(d["residueNumber"], d["confidenceScore"]))
    if "plddt" in d:
        return {i + 1: v for i, v in enumerate(d["plddt"])}
    return None


def load_plddt(path, uniprot=None):
    """Load AlphaFold confidence JSON from local path or fetch from API."""
    if path and os.path.exists(path):
        with open(path) as f:
            d = json.load(f)
        result = parse_plddt(d)
        if result:
            return result
    if uniprot:
        result = fetch_alphafold_plddt(uniprot)
        if result:
            return result
    return {}


# ────────────────────────────────────────────────── Classify environment
def classify_residue_environment(res_start, res_end, domains, protein_len):
    """Classify a residue range: which domains, linkers, termini it falls in.

    Uses only InterPro domain-type entries (IPR*) for domain boundary detection;
    whole-protein family annotations are excluded to avoid masking linkers.

    Returns a list of overlapping regions.
    """
    regions = []

    # Collect domain boundaries: prefer IPR entries with type "domain"
    domain_spans = []
    for d in domains:
        if d["accession"].startswith("IPR") and d["type"] == "domain":
            domain_spans.append((d["start"], d["end"], d["name"], d["accession"]))

    # If no IPR domains found, fall back to any domain/family entry
    if not domain_spans:
        for d in domains:
            if d["type"] in ("domain", "family"):
                # Skip whole-length annotations (>80% of protein)
                span_len = d["end"] - d["start"] + 1
                if span_len < protein_len * 0.8:
                    domain_spans.append((d["start"], d["end"], d["name"], d["accession"]))

    # Sort by start
    domain_spans.sort()

    # Compute linker regions (gaps between domains)
    linker_spans = []
    prev_end = 1
    for ds, de, dn, da in domain_spans:
        if ds > prev_end:
            linker_spans.append((prev_end, ds - 1, f"linker_{prev_end}-{ds-1}"))
        prev_end = max(prev_end, de)
    if prev_end < protein_len:
        linker_spans.append((prev_end, protein_len, f"c_terminal_{prev_end}-{protein_len}"))
    if domain_spans and domain_spans[0][0] > 1:
        linker_spans.insert(0, (1, domain_spans[0][0] - 1, f"n_terminal_1-{domain_spans[0][0]-1}"))

    # Check overlap
    overlap_start = max(1, res_start)
    overlap_end = min(protein_len, res_end)

    # Check domain overlap
    for ds, de, dn, da in domain_spans:
        if overlap_start <= de and overlap_end >= ds:
            o_start = max(overlap_start, ds)
            o_end = min(overlap_end, de)
            regions.append({
                "region_type": "domain",
                "name": dn,
                "accession": da,
                "overlap_start": o_start,
                "overlap_end": o_end,
                "overlap_len": o_end - o_start + 1,
                "description": f"residues {o_start}-{o_end} fall in {dn} ({da})",
            })

    # Check linker overlap
    for ls, le, ln in linker_spans:
        if overlap_start <= le and overlap_end >= ls:
            o_start = max(overlap_start, ls)
            o_end = min(overlap_end, le)
            regions.append({
                "region_type": "linker",
                "name": ln,
                "accession": None,
                "overlap_start": o_start,
                "overlap_end": o_end,
                "overlap_len": o_end - o_start + 1,
                "description": f"residues {o_start}-{o_end} in inter-domain linker",
            })

    return regions


# ──────────────────────────────────────── Compute affected residues from exons
def compute_affected_residues(classify_data):
    """From classify_variant.py JSON, compute which residues are affected.

    For the in-frame product (Δexon), returns the residue range removed.
    For coding variants, returns the single residue affected.
    """
    result = []

    splice = classify_data.get("splice_analysis")
    geo = classify_data.get("transcript_geometry", {})
    protein_len = geo.get("protein_len", 0)

    if splice:
        for c in splice.get("splicing_consequences", []):
            if c.get("reading_frame_maintained") and c.get("type") == "single_exon_skip":
                # In-frame exon skipping: compute residue range removed
                exon_num = c["skipped_exon"]
                if exon_num and geo:
                    ex = geo["exons"][exon_num - 1]
                    start_c = ex["start_c"]
                    end_c = ex["end_c"]
                    # Residue n spans c.(3n-2) to c.(3n)
                    # First residue affected: ceil(start_c / 3)
                    # Last residue affected: floor(end_c / 3)
                    res_start = (start_c + 2) // 3
                    res_end = end_c // 3
                    result.append({
                        "consequence": f"Δexon{exon_num}",
                        "affected_type": "deleted_region",
                        "residue_start": res_start,
                        "residue_end": res_end,
                        "length_aa": res_end - res_start + 1,
                        "description": f"deleted residues {res_start}-{res_end} ({res_end - res_start + 1} aa)",
                    })

    coding = classify_data.get("coding_analysis")
    if coding and coding.get("variant_class") in ("missense", "nonsense"):
        residue = coding.get("residue")
        if residue:
            result.append({
                "consequence": coding.get("hgvs_p", "?"),
                "affected_type": "point_mutation",
                "residue_start": residue,
                "residue_end": residue,
                "length_aa": 1,
                "description": f"residue {residue} ({coding.get('variant_class', '?')})",
            })

    return result


# ───────────────────────────────────────────────────── Main
def main():
    parser = argparse.ArgumentParser(
        description="Map variant-affected residues onto protein domains.")
    parser.add_argument("--gene", default="PRUNE1")
    parser.add_argument("--uniprot", default="Q86TP1")
    parser.add_argument("--variant", default="c.521-2A>G")
    parser.add_argument("--classify-json", default="variant_c.521_2A_gt_G.json",
                        help="JSON output from classify_variant.py")
    parser.add_argument("--plddt-json", default=None,
                        help="AlphaFold confidence JSON (optional; fetched from API if not given)")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    print(f"# protein_localisation.py — {args.gene} {args.variant}")
    print(f"# UniProt: {args.uniprot}")
    print(f"# Accessed: {ACCESS}")
    print()

    # 1. Load classify data
    with open(args.classify_json) as f:
        classify_data = json.load(f)
    geo = classify_data.get("transcript_geometry", {})

    # 2. Fetch domains
    print("== Domain annotations (InterPro)")
    domains = fetch_domains(args.uniprot)
    if not domains:
        print("  ! Could not fetch domains from InterPro; check UniProt ID")
    else:
        # Deduplicate: prefer InterPro canonical (IPR*) over member DB (PF*, SM*, PTHR*)
        seen_ranges = {}
        uniq = []
        for d in domains:
            key = (d["start"], d["end"])
            if d["type"] not in ("domain", "family"):
                continue
            is_ipr = d["accession"].startswith("IPR")
            existing = seen_ranges.get(key)
            if existing is None:
                seen_ranges[key] = d
                uniq.append(d)
            elif is_ipr and not existing["accession"].startswith("IPR"):
                # Replace with IPR version
                seen_ranges[key] = d
                for i, u in enumerate(uniq):
                    if (u["start"], u["end"]) == key:
                        uniq[i] = d
                        break

        print("\t".join(["Accession", "Name", "Type", "Start", "End"]))
        for d in uniq:
            print("\t".join([d["accession"], d["name"], d["type"],
                             str(d["start"]), str(d["end"])]))
        domains = uniq
    print()

    # 3. Load pLDDT
    plddt = {}
    if args.plddt_json and os.path.exists(args.plddt_json):
        plddt = load_plddt(args.plddt_json, uniprot=args.uniprot)
    elif args.uniprot:
        # Try fetching from AlphaFold API
        plddt = load_plddt(None, uniprot=args.uniprot)
        print(f"== pLDDT confidence: {len(plddt)} residues loaded")
        if plddt:
            print("\t".join(["Residue", "pLDDT"]))
            for r in sorted(plddt)[:5]:
                print(f"\t{r}\t{plddt[r]}")
            print("\t...")
        print()

    # 4. Compute affected residues
    affected = compute_affected_residues(classify_data)
    if not affected:
        print("  ! Could not determine affected residue range from variant")
        # Fallback: try to use coding analysis
        print()

    print("== Affected protein region(s)")
    print("\t".join(["Consequence", "Type", "Start (res)", "End (res)", "Length (aa)", "Description"]))
    for a in affected:
        print("\t".join([a["consequence"], a["affected_type"],
                         str(a["residue_start"]), str(a["residue_end"]),
                         str(a["length_aa"]), a["description"]]))
    print()

    # 5. Classify residue environment
    print("== Domain / linker overlap")
    protein_len = geo.get("protein_len", max(plddt.keys()) if plddt else 0)
    for a in affected:
        env = classify_residue_environment(a["residue_start"], a["residue_end"],
                                            domains, protein_len)
        print(f"--- {a['consequence']} ({a['description']}) ---")
        if env:
            print("\t".join(["Region type", "Name", "Overlap start", "Overlap end",
                             "Overlap (aa)", "Description"]))
            for e in env:
                print("\t".join([
                    e["region_type"],
                    e["name"],
                    str(e["overlap_start"]),
                    str(e["overlap_end"]),
                    str(e["overlap_len"]),
                    e["description"],
                ]))
        else:
            print("  No domain/linker overlap found")
        print()

    # 6. pLDDT at key residues
    if plddt:
        print("== pLDDT at affected residues")
        for a in affected:
            print(f"--- {a['consequence']} ---")
            print("\t".join(["Residue", "pLDDT"]))
            for r in range(a["residue_start"], a["residue_end"] + 1):
                if r in plddt:
                    print(f"\t{r}\t{plddt[r]:.1f}")
            # Also show flanking residues
            flank_start = max(1, a["residue_start"] - 3)
            flank_end = min(protein_len, a["residue_end"] + 3)
            print(f"  Flanking: residues {flank_start}-{a['residue_start']-1}, "
                  f"{a['residue_end']+1}-{flank_end}")
        print()

    # 7. Summary
    print("== Summary")
    for a in affected:
        env = classify_residue_environment(a["residue_start"], a["residue_end"],
                                            domains, protein_len)
        dom_parts = [e for e in env if e["region_type"] == "domain"]
        link_parts = [e for e in env if e["region_type"] == "linker"]
        print(f"  Affected region: residues {a['residue_start']}-{a['residue_end']} "
              f"({a['length_aa']} aa)")
        if dom_parts:
            for d in dom_parts:
                print(f"    Domain: {d['name']} ({d['accession']}), "
                      f"residues {d['overlap_start']}-{d['overlap_end']}")
        if link_parts:
            for lk in link_parts:
                print(f"    Linker: {lk['name']}, "
                      f"residues {lk['overlap_start']}-{lk['overlap_end']}")
        if not dom_parts and not link_parts:
            print(f"    No domain/linker overlap (outside annotated regions)")
        # pLDDT stats
        if plddt:
            scores = [plddt[r] for r in range(a["residue_start"], a["residue_end"]+1)
                      if r in plddt]
            if scores:
                print(f"    pLDDT across region: mean {sum(scores)/len(scores):.1f}, "
                      f"min {min(scores):.1f}, max {max(scores):.1f}")
        print()

    # 8. Save JSON
    result = {
        "accessed": ACCESS,
        "gene": args.gene,
        "uniprot": args.uniprot,
        "variant": args.variant,
        "domains": domains,
        "affected_regions": affected,
        "protein_length": protein_len,
    }
    safe_name = args.variant.replace(">", "_gt_").replace("-", "_").replace("*", "X")
    out_path = os.path.join(outdir, f"localisation_{args.gene}_{safe_name}.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"# Data saved to {out_path}")


if __name__ == "__main__":
    main()