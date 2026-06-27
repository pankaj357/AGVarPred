#!/usr/bin/env python3
"""
Extract VEP annotations from gnomAD VCF for the PanelApp benign benchmark.
Uses tabix for fast region queries.
"""

import os
import re
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
EXT_DIR = ROOT / "external_validation"
GNOMAD_VCF = os.environ.get(
    "GNOMAD_VCF",
    "external_data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
)
BENCH_PATH = EXT_DIR / "benchmarks" / "benchmark_panelapp_benign_conf3.csv"
OUT_PATH = EXT_DIR / "panelapp_benign_vep.parquet"

VEP_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "ALLELE_NUM", "DISTANCE", "STRAND",
    "FLAGS", "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID",
    "CANONICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL",
    "UNIPARC", "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET",
    "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF",
    "EA_MAF", "ExAC_MAF", "ExAC_Adj_MAF", "ExAC_AFR_MAF", "ExAC_AMR_MAF",
    "ExAC_EAS_MAF", "ExAC_FIN_MAF", "ExAC_NFE_MAF", "ExAC_OTH_MAF",
    "ExAC_SAS_MAF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME",
    "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter",
    "LoF_flags", "LoF_info"
]
VEP_IDX = {f: i for i, f in enumerate(VEP_FIELDS)}


def parse_vep(vep_string):
    if not vep_string:
        return {}
    first_ann = vep_string.split(",")[0]
    parts = first_ann.split("|")
    if len(parts) < 65:
        return {}

    def get(idx):
        val = parts[idx] if idx < len(parts) else ""
        return val if val != "" else None

    sift_raw = get(VEP_IDX["SIFT"])
    sift_score = None
    sift_pred = None
    if sift_raw:
        m = re.search(r"([^(]+)\(([\d.]+)\)", sift_raw)
        if m:
            try:
                sift_score = float(m.group(2))
            except ValueError:
                pass
            sift_pred = m.group(1)

    polyphen_raw = get(VEP_IDX["PolyPhen"])
    polyphen_score = None
    polyphen_pred = None
    if polyphen_raw:
        m = re.search(r"([^(]+)\(([\d.]+)\)", polyphen_raw)
        if m:
            try:
                polyphen_score = float(m.group(2))
            except ValueError:
                pass
            polyphen_pred = m.group(1)

    consequence = get(VEP_IDX["Consequence"])
    impact = get(VEP_IDX["IMPACT"])
    lof = get(VEP_IDX["LoF"])
    protein_pos = get(VEP_IDX["Protein_position"])
    impact_map = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
    impact_score = impact_map.get(impact, 0) if impact else 0
    protein_pos_num = None
    if protein_pos:
        m = re.search(r"^(\d+)", protein_pos)
        if m:
            try:
                protein_pos_num = int(m.group(1))
            except ValueError:
                pass
    cons_list = consequence.split("&") if consequence else []
    return {
        "vep_SIFT_score": sift_score,
        "vep_SIFT_pred": sift_pred,
        "vep_PolyPhen_score": polyphen_score,
        "vep_PolyPhen_pred": polyphen_pred,
        "vep_IMPACT": impact,
        "vep_IMPACT_score": impact_score,
        "vep_Consequence": consequence,
        "vep_is_missense": 1 if "missense_variant" in cons_list else 0,
        "vep_is_synonymous": 1 if "synonymous_variant" in cons_list else 0,
        "vep_is_stop_gained": 1 if "stop_gained" in cons_list else 0,
        "vep_is_frameshift": 1 if "frameshift_variant" in cons_list else 0,
        "vep_is_splice": 1 if any(c in cons_list for c in ["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"]) else 0,
        "vep_LoF": lof,
        "vep_is_LoF_HC": 1 if lof == "HC" else 0,
        "vep_Protein_position": protein_pos_num,
    }


def main():
    print("Loading PanelApp benign benchmark...")
    df = pd.read_csv(BENCH_PATH)
    df["variant_id"] = (
        "chr" + df["chrom"].astype(str).str.replace("chr", "", case=False)
        + "_" + df["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
        + "_" + df["ref"].astype(str).str.upper()
        + "_" + df["alt"].astype(str).str.upper()
    )
    target_variants = df["variant_id"].tolist()
    print(f"PanelApp benign variants to query: {len(target_variants):,}")

    # Build BED
    bed_entries = []
    vid_lookup = {}
    for vid in target_variants:
        parts = vid.split("_")
        if len(parts) != 4:
            continue
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
        bed_entries.append((chrom, pos - 1, pos, vid))
        vid_lookup[(chrom, pos, ref, alt)] = vid

    bed_entries.sort(key=lambda x: (x[0], x[1]))
    bed_path = EXT_DIR / "raw_data" / "panelapp_benign_variants.bed"
    bed_path.parent.mkdir(exist_ok=True, parents=True)
    with open(bed_path, "w") as fh:
        for chrom, start, end, vid in bed_entries:
            fh.write(f"{chrom}\t{start}\t{end}\n")
    print(f"BED file: {bed_path} ({len(bed_entries):,} entries)")

    # Run tabix
    print("Running tabix...")
    cmd = ["tabix", "-R", str(bed_path), GNOMAD_VCF]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"tabix failed: {proc.stderr}")
    lines = proc.stdout.strip().split("\n")
    print(f"Tabix returned {len(lines):,} lines")

    # Parse VEP
    vep_map = {}
    for line in lines:
        cols = line.split("\t")
        if len(cols) < 8:
            continue
        chrom, pos, _, ref, alt = cols[0], int(cols[1]), cols[2], cols[3], cols[4]
        info = cols[7]
        alts = alt.split(",")
        for i, a in enumerate(alts):
            key = (chrom, pos, ref, a)
            if key in vid_lookup:
                vid = vid_lookup[key]
                m = re.search(r'vep=([^;]+)', info)
                if m:
                    vep_map[vid] = parse_vep(m.group(1))

    print(f"VEP annotations found for {len(vep_map):,} / {len(target_variants):,} variants")

    rows = []
    for vid in target_variants:
        row = {"variant_id": vid}
        row.update(vep_map.get(vid, {}))
        rows.append(row)
    out_df = pd.DataFrame(rows).set_index("variant_id")
    out_df.to_parquet(OUT_PATH)
    print(f"Saved {OUT_PATH} ({len(out_df):,} variants, {len(out_df.columns)} VEP features)")


if __name__ == "__main__":
    main()
