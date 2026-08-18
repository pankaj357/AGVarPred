#!/usr/bin/env python3
"""
Extract VEP annotations from gnomAD VCF for all external validation benchmarks.
Saves VEP features as parquet files that scoring scripts can merge.
"""

import os
import re
import pandas as pd
import numpy as np
import pysam
from pathlib import Path
from tqdm import tqdm

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
GNOMAD_VCF = os.environ.get(
    "GNOMAD_VCF",
    "external_data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
)
EXT_DIR = ROOT / "external_validation"

# VEP parsing (same as feature_labbling.py)
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
        m = re.search(r"([^(]+)\(([\\d.]+)\)", sift_raw)
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
        m = re.search(r"([^(]+)\(([\\d.]+)\)", polyphen_raw)
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

BENCHMARKS = [
    ("benchmarks/benchmark_independent_humsavar.csv", "humsavar_vep.parquet"),
    ("benchmarks/benchmark_mave_independent.csv", "mave_independent_vep.parquet"),
    ("benchmarks/benchmark_gnomad_benign.csv", "gnomad_benign_vep.parquet"),
]

def extract_variant_ids(bench_path):
    df = pd.read_csv(bench_path)
    df["variant_id"] = (
        "chr" + df["chrom"].astype(str).str.replace("chr", "", case=False)
        + "_" + df["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
        + "_" + df["ref"].astype(str).str.upper()
        + "_" + df["alt"].astype(str).str.upper()
    )
    return set(df["variant_id"])

def main():
    print("Loading gnomAD VCF...")
    vcf = pysam.VariantFile(GNOMAD_VCF)
    
    # Collect all variant IDs
    all_variants = set()
    for bench_file, _ in BENCHMARKS:
        bench_path = EXT_DIR / bench_file
        if bench_path.exists():
            vids = extract_variant_ids(bench_path)
            all_variants.update(vids)
            print(f"  {bench_file}: {len(vids):,} variants")
    print(f"Total unique variants to query: {len(all_variants):,}")
    
    # Query VCF
    vep_map = {}
    count = 0
    for rec in tqdm(vcf.fetch(), desc="Querying gnomAD VCF"):
        try:
            chrom = rec.chrom
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            pos = rec.pos
            ref = rec.ref
            if rec.alts is None:
                continue
            alt = rec.alts[0]
            variant_id = f"{chrom}_{pos}_{ref}_{alt}"
            if variant_id not in all_variants:
                continue
            vep = rec.info.get("vep")
            if vep:
                vep_map[variant_id] = parse_vep(vep[0])
            count += 1
        except Exception:
            continue
    print(f"VEP annotations found for {len(vep_map):,} / {len(all_variants):,} variants")
    
    # Save per-benchmark VEP files
    for bench_file, out_file in BENCHMARKS:
        bench_path = EXT_DIR / bench_file
        if not bench_path.exists():
            continue
        vids = extract_variant_ids(bench_path)
        rows = []
        for vid in vids:
            row = {"variant_id": vid}
            row.update(vep_map.get(vid, {}))
            rows.append(row)
        df = pd.DataFrame(rows).set_index("variant_id")
        out_path = EXT_DIR / out_file
        df.to_parquet(out_path)
        print(f"Saved {out_path} ({len(df):,} variants, {len(df.columns)} VEP features)")
    
    print("Done!")

if __name__ == "__main__":
    main()
