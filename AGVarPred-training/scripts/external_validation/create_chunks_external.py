#!/usr/bin/env python3
"""
Create AlphaGenome input chunks for external validation benchmarks.

Each benchmark gets its own chunk directory:
  external_validation/processing/chunks/external_humsavar/
  external_validation/processing/chunks/mave_independent/
  external_validation/processing/chunks/gnomad_benign/

Inside each: run_1/ ... run_20/ folders with per-gene .txt files.
Format: tab-separated, no header, columns = chrom, pos, ref, alt, gene
Compatible with code_external.py.
"""

import os
import math
import pandas as pd
from pathlib import Path

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
OUT_BASE = ROOT / "external_validation"
N_RUNS = 20

# Dataset configs: (csv_file, output_dir_name)
DATASETS = [
    ("benchmarks/benchmark_independent_humsavar.csv", "processing/chunks/external_humsavar"),
    ("benchmarks/benchmark_mave_independent.csv", "processing/chunks/mave_independent"),
    ("benchmarks/benchmark_gnomad_benign.csv", "processing/chunks/gnomad_benign"),
]


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    if chrom.lower() in ("chrmt", "chrmitochondria"):
        chrom = "chrM"
    return chrom


def create_chunks_for_dataset(csv_name: str, out_dir_name: str):
    csv_path = OUT_BASE / csv_name
    out_dir = OUT_BASE / out_dir_name

    if not csv_path.exists():
        print(f"\n⚠️  Skipping {csv_name} (not found yet)")
        return

    print(f"\n{'='*60}")
    print(f"Processing: {csv_name}")
    print(f"Output dir: {out_dir}")
    print(f"{'='*60}")

    df = pd.read_csv(csv_path, dtype=str)
    print(f"  Raw rows: {len(df):,}")

    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt", "gene"])
    print(f"  After dedup: {len(df):,}")

    df["chrom"] = df["chrom"].apply(normalize_chrom)
    df["pos"] = df["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()

    gene_groups = df.groupby("gene")
    genes = sorted(gene_groups.groups.keys())
    print(f"  Genes: {len(genes):,}")

    for r in range(1, N_RUNS + 1):
        (out_dir / f"run_{r}").mkdir(parents=True, exist_ok=True)

    files_created = 0
    variants_written = 0
    for i, gene in enumerate(genes):
        run_num = (i % N_RUNS) + 1
        gene_df = gene_groups.get_group(gene)
        gene_file = out_dir / f"run_{run_num}" / f"{gene}.txt"

        with open(gene_file, "w") as fh:
            for _, row in gene_df.iterrows():
                fh.write(f"{row['chrom']}\t{row['pos']}\t{row['ref']}\t{row['alt']}\t{gene}\n")

        files_created += 1
        variants_written += len(gene_df)

    print(f"  Files created: {files_created:,}")
    print(f"  Variants written: {variants_written:,}")

    print("  Per-run distribution:")
    for r in range(1, N_RUNS + 1):
        run_path = out_dir / f"run_{r}"
        n_files = len(list(run_path.glob("*.txt")))
        if n_files > 0:
            print(f"    run_{r:2d}: {n_files:4d} files")


def main():
    print("🚀 Creating AlphaGenome input chunks for external validation")

    for csv_name, out_dir_name in DATASETS:
        create_chunks_for_dataset(csv_name, out_dir_name)

    print("\n" + "="*60)
    print("✅ All chunks created successfully!")
    print("="*60)
    print("\nTo run feature extraction, set INPUT_BASE env var to one of:")
    for _, out_dir_name in DATASETS:
        print(f"  - external_validation/{out_dir_name}")


if __name__ == "__main__":
    main()
