#!/usr/bin/env python3
"""
Build MAVE Benchmark from MaveDB mapping archives — NON-TRAINING GENES ONLY.

Uses pre-mapped MAVE data from MaveDB mapping archive.
For each target gene NOT present in train/cal/test, extracts variants with
genomic coordinates (GRCh38), converts VRS Alleles to VCF-style, assigns
binary labels based on score percentiles, and deduplicates against train/cal/test.

BRCA1, PTEN, TP53 and any other training genes are EXCLUDED to maintain
strict gene-holdout independence.
"""

import pandas as pd
import numpy as np
import json
import os
from pathlib import Path
from collections import defaultdict

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
EXT = ROOT / "external_validation"
EXTRACT_DIR = EXT / "source_datasets" / "mave_multi_gene" / "mapped_data_extracted" / "mappings"

NC_TO_CHROM = {
    "NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3",
    "NC_000004.12": "4", "NC_000005.10": "5", "NC_000006.12": "6",
    "NC_000007.14": "7", "NC_000008.11": "8", "NC_000009.12": "9",
    "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
    "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15",
    "NC_000016.10": "16", "NC_000017.11": "17", "NC_000018.10": "18",
    "NC_000019.10": "19", "NC_000020.11": "20", "NC_000021.9": "21",
    "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y",
    "NC_012920.1": "MT",
}

MIN_VARIANTS = 50


def load_training_genes():
    """Load all gene symbols from train/cal/test for strict holdout."""
    genes = set()
    for split in ["train.csv", "cal.csv", "test.csv"]:
        df = pd.read_csv(ROOT / split, usecols=["GeneSymbol"], dtype=str)
        genes.update(df["GeneSymbol"].dropna().unique())
    print(f"   Training genes to exclude: {len(genes):,}")
    return genes


def load_dedup_set():
    train_ids = set()
    for split in ["train.csv", "cal.csv", "test.csv"]:
        df = pd.read_csv(ROOT / split, usecols=["Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"])
        df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False)
        vids = (df["Chromosome"] + "_" + df["PositionVCF"].astype(str).str.replace(r"\.0$", "", regex=True) +
                "_" + df["ReferenceAlleleVCF"].str.upper() + "_" + df["AlternateAlleleVCF"].str.upper())
        train_ids.update(vids)
    return train_ids


def build_refget_to_chrom(data):
    mapping = {}
    mrs = data.get("mapped_reference_sequence", {})
    seq_id = mrs.get("sequence_id", "")
    if seq_id.startswith("ga4gh:"):
        refget_acc = seq_id.split(":", 1)[1]
    else:
        refget_acc = seq_id
    for nc_acc in mrs.get("sequence_accessions", []):
        chrom = NC_TO_CHROM.get(nc_acc)
        if chrom:
            mapping[refget_acc] = chrom
            break
    return mapping


def vrs_to_vcf(post_mapped, refget_to_chrom):
    loc = post_mapped.get("location", {})
    seq_ref = loc.get("sequenceReference", {})
    refget_acc = seq_ref.get("refgetAccession", "")
    chrom = refget_to_chrom.get(refget_acc)
    if chrom is None:
        return None
    start = loc.get("start")
    end = loc.get("end")
    ref_seq = loc.get("sequence", "")
    alt_seq = post_mapped.get("state", {}).get("sequence", "")
    if start is None or end is None:
        return None
    pos = start + 1
    if not ref_seq or alt_seq is None:
        return None
    return {"chrom": chrom, "pos": int(pos), "ref": ref_seq, "alt": alt_seq}


def match_gene(target_genes, training_genes):
    """Find matching target gene, but reject if it's in training."""
    for g in target_genes:
        name = g.get("name", "").upper()
        # Extract first token as the gene symbol (e.g., "BRCA1 RING domain" -> "BRCA1")
        symbol = name.split()[0] if name else ""
        if symbol and symbol not in training_genes:
            return symbol
    return None


def _alleles_from_post_mapped(post):
    post_type = post.get("type")
    if post_type == "Allele":
        return [post]
    if post_type == "CisPhasedBlock":
        return [m for m in post.get("members", []) if m.get("type") == "Allele"]
    return []


def process_json_file(json_path, train_ids, training_genes):
    if os.path.basename(json_path).startswith("._"):
        return []
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError, UnicodeDecodeError):
        return []

    target_genes = data.get("metadata", {}).get("targetGenes", [])
    gene = match_gene(target_genes, training_genes)
    if gene is None:
        return []

    refget_to_chrom = build_refget_to_chrom(data)
    if not refget_to_chrom:
        return []

    scores = []
    for ms in data.get("mapped_scores", []):
        post = ms.get("post_mapped", {})
        score = ms.get("score")
        if score is None or (isinstance(score, float) and np.isnan(score)):
            continue
        for allele in _alleles_from_post_mapped(post):
            vcf = vrs_to_vcf(allele, refget_to_chrom)
            if vcf is None:
                continue
            # Keep SNVs only
            if len(vcf["ref"]) != 1 or len(vcf["alt"]) != 1:
                continue
            vid = f"{vcf['chrom']}_{vcf['pos']}_{vcf['ref']}_{vcf['alt']}"
            if vid in train_ids:
                continue
            scores.append({
                "chrom": vcf["chrom"],
                "pos": vcf["pos"],
                "ref": vcf["ref"],
                "alt": vcf["alt"],
                "gene": gene,
                "score": float(score),
                "mavedb_id": ms.get("mavedb_id", ""),
            })
    return scores


def extract_and_process():
    print("=" * 70)
    print("MAVE MULTI-GENE BENCHMARK BUILDER (NON-TRAINING GENES ONLY)")
    print("=" * 70)

    training_genes = load_training_genes()
    train_ids = load_dedup_set()
    print(f"   Train/cal/test dedup set: {len(train_ids):,} variants")

    if not EXTRACT_DIR.exists():
        print(f"❌ Extraction directory not found: {EXTRACT_DIR}")
        print("   Please extract the tar.gz first.")
        return

    json_files = [f for f in EXTRACT_DIR.glob("*.json") if not f.name.startswith("._")]
    print(f"   Found {len(json_files):,} JSON files")

    all_scores = []
    for i, json_path in enumerate(json_files, 1):
        scores = process_json_file(str(json_path), train_ids, training_genes)
        if scores:
            all_scores.extend(scores)
            print(f"   [{i}/{len(json_files)}] {json_path.name}: {len(scores)} variants for {scores[0]['gene']}")

    if not all_scores:
        print("   ⚠️ No variants found for non-training genes")
        return

    df = pd.DataFrame(all_scores)
    print(f"\n   Total variants before dedup: {len(df):,}")
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    print(f"   After exact dedup: {len(df):,}")

    # Assign binary labels per gene using score percentiles
    labeled = []
    for gene, gdf in df.groupby("gene"):
        if len(gdf) < MIN_VARIANTS:
            print(f"   ⚠️ Skipping {gene}: only {len(gdf)} variants")
            continue

        low_thr = gdf["score"].quantile(0.20)
        high_thr = gdf["score"].quantile(0.80)

        pathogenic = gdf[gdf["score"] <= low_thr].copy()
        pathogenic["label"] = 1
        pathogenic["dataset"] = f"mave_{gene.lower()}_pathogenic"

        benign = gdf[gdf["score"] >= high_thr].copy()
        benign["label"] = 0
        benign["dataset"] = f"mave_{gene.lower()}_benign"

        labeled.append(pathogenic)
        labeled.append(benign)

        print(f"   {gene}: {len(pathogenic)} pathogenic (score≤{low_thr:.3f}), {len(benign)} benign (score≥{high_thr:.3f})")

    if not labeled:
        print("   ⚠️ No labeled variants")
        return

    final = pd.concat(labeled, ignore_index=True)
    final = final[["chrom", "pos", "ref", "alt", "gene", "label", "score", "dataset"]]

    out = EXT / "benchmarks/benchmark_mave_independent.csv"
    final.to_csv(out, index=False)

    print(f"\n   ✅ Saved {len(final):,} variants to {out}")
    print(f"   Pathogenic: {(final['label']==1).sum():,}")
    print(f"   Benign: {(final['label']==0).sum():,}")
    for gene in sorted(final["gene"].unique()):
        g = final[final["gene"] == gene]
        print(f"      {gene}: {(g['label']==1).sum()} pathogenic, {(g['label']==0).sum()} benign")


if __name__ == "__main__":
    extract_and_process()
