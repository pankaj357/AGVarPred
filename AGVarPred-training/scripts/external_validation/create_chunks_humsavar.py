#!/usr/bin/env python3
"""
Create AlphaGenome input chunks for the Humsavar independent benchmark.

Output: external_validation/processing/chunks/external_humsavar/
        run_1/ ... run_20/ folders with per-gene .txt files
Format: tab-separated, no header, columns = chrom, pos, ref, alt, gene
"""

import pandas as pd
from pathlib import Path

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
OUT_BASE = ROOT / "external_validation"
CSV_PATH = OUT_BASE / "benchmarks/benchmark_independent_humsavar.csv"
OUT_DIR = OUT_BASE / "processing/chunks/external_humsavar"
N_RUNS = 20


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    if chrom.lower() in ("chrmt", "chrmitochondria"):
        chrom = "chrM"
    return chrom


def main():
    print("🚀 Creating AlphaGenome chunks for Humsavar benchmark")
    print(f"   Source: {CSV_PATH}")
    print(f"   Output: {OUT_DIR}")

    df = pd.read_csv(CSV_PATH, dtype=str)
    print(f"   Raw rows: {len(df):,}")

    # Deduplicate exact same variant-gene combinations
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt", "gene"])
    print(f"   After dedup: {len(df):,}")

    df["chrom"] = df["chrom"].apply(normalize_chrom)
    df["pos"] = df["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
    df["ref"] = df["ref"].astype(str).str.upper()
    df["alt"] = df["alt"].astype(str).str.upper()

    gene_groups = df.groupby("gene")
    genes = sorted(gene_groups.groups.keys())
    print(f"   Genes: {len(genes):,}")

    for r in range(1, N_RUNS + 1):
        (OUT_DIR / f"run_{r}").mkdir(parents=True, exist_ok=True)

    files_created = 0
    variants_written = 0
    for i, gene in enumerate(genes):
        run_num = (i % N_RUNS) + 1
        gene_df = gene_groups.get_group(gene)
        gene_file = OUT_DIR / f"run_{run_num}" / f"{gene}.txt"

        with open(gene_file, "w") as fh:
            for _, row in gene_df.iterrows():
                fh.write(f"{row['chrom']}\t{row['pos']}\t{row['ref']}\t{row['alt']}\t{gene}\n")

        files_created += 1
        variants_written += len(gene_df)

    print(f"   Files created: {files_created:,}")
    print(f"   Variants written: {variants_written:,}")
    print("✅ Done!")


if __name__ == "__main__":
    main()
