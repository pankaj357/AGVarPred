#!/usr/bin/env python3
"""
Split ClinGen benchmark into per-gene chunks for AlphaGenome feature extraction.

Creates:
  external_validation/processing/chunks/external_clingen/run_{1..N}/GENE.txt
"""

import os
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
EXT_DIR = PROJECT_ROOT / "external_validation"
BENCH_PATH = EXT_DIR / "benchmarks" / "benchmark_clingen.csv"
CHUNK_DIR = EXT_DIR / "processing" / "chunks" / "external_clingen"
N_RUNS = 20


def normalize_chrom(chrom: str) -> str:
    """Ensure chromosomes have the 'chr' prefix required by AlphaGenome."""
    chrom = str(chrom).strip()
    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    if chrom.lower() in ("chrmt", "chrmitochondria"):
        chrom = "chrM"
    return chrom


def main():
    print("Preparing ClinGen chunks for AlphaGenome feature extraction...")
    bench = pd.read_csv(BENCH_PATH)
    print(f"Total variants: {len(bench):,}")
    print(f"Genes: {bench['gene'].nunique():,}")

    # Create run directories
    for i in range(1, N_RUNS + 1):
        (CHUNK_DIR / f"run_{i}").mkdir(parents=True, exist_ok=True)

    # Assign genes to runs round-robin
    genes = sorted(bench["gene"].unique())
    for idx, gene in enumerate(genes):
        run_id = (idx % N_RUNS) + 1
        gene_df = bench[bench["gene"] == gene][["chrom", "pos", "ref", "alt", "gene"]]
        # Ensure chr prefix (AlphaGenome API requires it)
        gene_df["chrom"] = gene_df["chrom"].apply(normalize_chrom)
        out_path = CHUNK_DIR / f"run_{run_id}" / f"{gene}.txt"
        gene_df.to_csv(out_path, sep="\t", header=False, index=False)

    print(f"Created chunks for {len(genes)} genes across {N_RUNS} runs in {CHUNK_DIR}")


if __name__ == "__main__":
    main()
