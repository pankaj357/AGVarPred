#!/usr/bin/env python3
"""
Build an independent benchmark from ClinGen Variant Pathogenicity curations.

Source:
    ClinGen Evidence Repository Variant Pathogenicity summary
    https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download?type=csv

Processing:
  1. Download / load ClinGen curated variant classifications
  2. Parse GRCh38 genomic HGVS expressions for simple SNVs
  3. Assign binary labels from expert-panel assertions:
     - Pathogenic / Likely Pathogenic -> 1
     - Benign / Likely Benign -> 0
     - Uncertain Significance -> excluded
  4. Apply strict gene holdout (exclude genes in train/cal/test)
  5. Deduplicate against train/cal/test by chrom/pos/ref/alt
  6. Remove variants already present in other active benchmarks
  7. Save benchmark CSV

Independence note:
  This benchmark now applies the same strict gene + variant holdout used by the
  other independent benchmarks. Because ClinGen expert panels curate established
  disease genes, this typically leaves a much smaller benchmark than the
  variant-only holdout version.

Output:
  external_validation/benchmarks/benchmark_clingen.csv
"""

import os
import re
import json
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
EXT_DIR = PROJECT_ROOT / "external_validation"
SOURCE_DIR = EXT_DIR / "raw_data"
SOURCE_FILE = SOURCE_DIR / "clingen_variant_pathogenicity.csv"
SOURCE_URL = "https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download?type=csv"
OUT_DIR = EXT_DIR / "benchmarks"
OUT_BENCHMARK = OUT_DIR / "benchmark_clingen.csv"

TRAIN_VARIANTS_CACHE = EXT_DIR / "processing" / "clingen_train_variants.pkl"
TRAIN_GENES_CACHE = EXT_DIR / "processing" / "clingen_train_genes.pkl"

# GRCh38 RefSeq chromosome accessions
NC_TO_CHR = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.10": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "M",
}


def collect_training_variants():
    """Collect all variant IDs from train/cal/test parquet parts."""
    if TRAIN_VARIANTS_CACHE.exists():
        print(f"Loading cached training variants from {TRAIN_VARIANTS_CACHE}")
        return pd.read_pickle(TRAIN_VARIANTS_CACHE)

    import pyarrow.dataset as ds

    variants = set()
    for split in ["train", "cal", "test"]:
        parts_dir = PROJECT_ROOT / f"final_dataset_parts_{split}"
        if not parts_dir.exists():
            print(f"Warning: {parts_dir} not found, skipping")
            continue
        print(f"Scanning {split} parquet dataset...")
        dataset = ds.dataset(parts_dir, format="parquet")
        table = dataset.to_table(columns=["variant_id"])
        split_variants = set(table.column("variant_id").to_pylist())
        variants.update(split_variants)
        print(f"  {split}: {len(split_variants):,} unique variants")

    TRAIN_VARIANTS_CACHE.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(variants, TRAIN_VARIANTS_CACHE)
    print(f"Cached {len(variants):,} total training variants")
    return variants


def collect_training_genes():
    """Collect all gene symbols from train/cal/test CSVs."""
    if TRAIN_GENES_CACHE.exists():
        print(f"Loading cached training genes from {TRAIN_GENES_CACHE}")
        return pd.read_pickle(TRAIN_GENES_CACHE)

    genes = set()
    for split in ["train", "cal", "test"]:
        csv_path = PROJECT_ROOT / f"{split}.csv"
        if not csv_path.exists():
            print(f"Warning: {csv_path} not found, skipping")
            continue
        df = pd.read_csv(csv_path, usecols=["GeneSymbol"], dtype=str)
        split_genes = set(df["GeneSymbol"].dropna().unique())
        genes.update(split_genes)
        print(f"  {split}: {len(split_genes):,} unique genes")

    TRAIN_GENES_CACHE.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(genes, TRAIN_GENES_CACHE)
    print(f"Cached {len(genes):,} total training genes")
    return genes


def parse_simple_snv(hgvs_expr):
    """
    Parse a GRCh38 genomic HGVS expression for a simple substitution SNV.
    Returns (chrom, pos, ref, alt) or None.
    """
    if not isinstance(hgvs_expr, str):
        return None

    # Find a genomic expression on a GRCh38 chromosome
    for expr in hgvs_expr.split(","):
        expr = expr.strip()
        m = re.match(r"(NC_\d+\.\d+):g\.(\d+)([ACGT])>([ACGT])$", expr)
        if not m:
            continue
        nc, pos, ref, alt = m.groups()
        chrom = NC_TO_CHR.get(nc)
        if chrom is None:
            continue
        return f"chr{chrom}", int(pos), ref.upper(), alt.upper()

    return None


def load_existing_benchmark_variant_ids():
    """Load variant IDs from other active benchmarks to avoid overlap."""
    existing = set()
    bench_names = ["humsavar", "mave_independent", "gnomad_benign", "vip", "grimm2015"]
    for name in bench_names:
        path = OUT_DIR / f"benchmark_{name}.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path, dtype=str)
        if "variant_id" in df.columns:
            existing.update(df["variant_id"].dropna().unique())
        elif all(c in df.columns for c in ["chrom", "pos", "ref", "alt"]):
            for _, row in df.iterrows():
                vid = f"chr{row['chrom']}_{row['pos']}_{row['ref']}_{row['alt']}"
                existing.add(vid)
    print(f"Existing benchmark variant IDs: {len(existing):,}")
    return existing


def build_benchmark():
    print("=" * 70)
    print("Building ClinGen independent benchmark")
    print("=" * 70)

    # Ensure source file exists
    if not SOURCE_FILE.exists():
        print(f"Downloading ClinGen variant pathogenicity data from {SOURCE_URL}")
        SOURCE_DIR.mkdir(parents=True, exist_ok=True)
        import urllib.request
        urllib.request.urlretrieve(SOURCE_URL, SOURCE_FILE)
        print(f"Saved to {SOURCE_FILE}")

    # Load ClinGen data
    print(f"\nLoading {SOURCE_FILE}")
    df = pd.read_csv(SOURCE_FILE, sep="\t", dtype=str)
    print(f"  Total rows: {len(df):,}")

    # Load training genes and variants
    print("\nLoading training holdout sets...")
    train_genes = collect_training_genes()
    train_variants = collect_training_variants()

    # Parse coordinates
    print("\nParsing GRCh38 genomic coordinates...")
    parsed = df["HGVS Expressions"].apply(parse_simple_snv)
    df["chrom"] = parsed.apply(lambda x: x[0] if x else None)
    df["pos"] = parsed.apply(lambda x: x[1] if x else None)
    df["ref"] = parsed.apply(lambda x: x[2] if x else None)
    df["alt"] = parsed.apply(lambda x: x[3] if x else None)

    before = len(df)
    df = df[df["chrom"].notna()].copy()
    print(f"  Kept simple GRCh38 SNVs: {len(df):,} / {before:,}")

    # Build variant_id
    df["variant_id"] = (
        df["chrom"].astype(str) + "_" +
        df["pos"].astype(int).astype(str) + "_" +
        df["ref"].astype(str) + "_" +
        df["alt"].astype(str)
    )

    # Map assertions to labels
    assertion_map = {
        "Pathogenic": 1,
        "Likely Pathogenic": 1,
        "Likely Benign": 0,
        "Benign": 0,
    }
    df["label"] = df["Assertion"].map(assertion_map)

    before = len(df)
    df = df[df["label"].notna()].copy()
    print(f"  After excluding VUS/retracted: {len(df):,} / {before:,}")
    print(f"    Pathogenic: {(df['label'] == '1').sum():,}")
    print(f"    Benign: {(df['label'] == '0').sum():,}")

    # Convert label to int
    df["label"] = df["label"].astype(int)

    # Retracted filter
    if "Retracted" in df.columns:
        before = len(df)
        df = df[df["Retracted"].astype(str).str.lower() != "true"].copy()
        print(f"  After removing retracted: {len(df):,} / {before:,}")

    # Gene holdout (strict: exclude genes in train/cal/test)
    before = len(df)
    df = df[~df["HGNC Gene Symbol"].isin(train_genes)].copy()
    print(f"  After training-gene exclusion: {len(df):,} / {before:,}")

    # Variant holdout (strict: exact chrom/pos/ref/alt)
    before = len(df)
    df = df[~df["variant_id"].isin(train_variants)].copy()
    print(f"  After training-variant exclusion: {len(df):,} / {before:,}")

    # Remove overlaps with existing benchmarks
    existing_ids = load_existing_benchmark_variant_ids()
    before = len(df)
    df = df[~df["variant_id"].isin(existing_ids)].copy()
    print(f"  After existing-benchmark exclusion: {len(df):,} / {before:,}")

    # Final deduplication within ClinGen
    before = len(df)
    df = df.drop_duplicates(subset=["variant_id"]).copy()
    print(f"  After internal deduplication: {len(df):,} / {before:,}")

    # Build output
    out = pd.DataFrame({
        "chrom": df["chrom"].str.replace("chr", ""),
        "pos": df["pos"].astype(int),
        "ref": df["ref"],
        "alt": df["alt"],
        "variant_id": df["variant_id"],
        "gene": df["HGNC Gene Symbol"],
        "disease": df["Disease"],
        "mondo_id": df["Mondo Id"],
        "assertion": df["Assertion"],
        "expert_panel": df["Expert Panel"],
        "label": df["label"],
        "category": df["label"].map({1: "P/LP", 0: "B/LB"}),
        "approval_date": df["Approval Date"],
        "uuid": df["Uuid"],
    })

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_BENCHMARK, index=False)

    print("\n" + "=" * 70)
    print(f"Saved benchmark to {OUT_BENCHMARK}")
    print(f"Final benchmark size: {len(out):,}")
    print(f"  Pathogenic (P/LP): {(out['label'] == 1).sum():,}")
    print(f"  Benign (B/LB): {(out['label'] == 0).sum():,}")
    print(f"  Genes: {out['gene'].nunique():,}")
    print(f"  Expert panels: {out['expert_panel'].nunique():,}")
    print("=" * 70)

    # Save summary
    summary_path = OUT_BENCHMARK.with_suffix(".summary.txt")
    with open(summary_path, "w") as f:
        f.write("ClinGen Independent Benchmark Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total variants: {len(out)}\n")
        f.write(f"Pathogenic (P/LP): {(out['label'] == 1).sum()}\n")
        f.write(f"Benign (B/LB): {(out['label'] == 0).sum()}\n")
        f.write(f"Genes: {out['gene'].nunique()}\n")
        f.write(f"Expert panels: {out['expert_panel'].nunique()}\n")
        f.write("\nExpert panel counts:\n")
        f.write(out["expert_panel"].value_counts().to_string())
        f.write("\n\nAssertion counts:\n")
        f.write(out["assertion"].value_counts().to_string())
    print(f"Saved summary to {summary_path}")


if __name__ == "__main__":
    build_benchmark()
