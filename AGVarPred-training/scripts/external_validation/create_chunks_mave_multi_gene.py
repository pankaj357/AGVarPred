#!/usr/bin/env python3
"""
Create chunk files for MAVE multi-gene benchmark feature extraction.
Splits variants across 20 parallel runs.
"""

import os
import pandas as pd
from pathlib import Path

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
BENCH = ROOT / "external_validation/benchmarks/benchmark_mave_multi_gene.csv"
OUT_DIR = ROOT / "external_validation/processing/chunks/mave_multi_gene"
N_RUNS = 20

def normalize_chrom(chrom):
    chrom = str(chrom).strip()
    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    return chrom

print("🚀 Creating MAVE multi-gene chunks...")
df = pd.read_csv(BENCH)
print(f"   Variants: {len(df):,} | Genes: {df.gene.nunique():,}")

# Normalize
df["chrom"] = df["chrom"].apply(normalize_chrom)
df["pos"] = df["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
df["ref"] = df["ref"].astype(str).str.upper()
df["alt"] = df["alt"].astype(str).str.upper()

# Remove duplicates
df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt", "gene"])
print(f"   After dedup: {len(df):,}")

# Split into N_RUNS chunks
chunk_size = (len(df) + N_RUNS - 1) // N_RUNS

for run_id in range(1, N_RUNS + 1):
    run_dir = OUT_DIR / f"run_{run_id}"
    run_dir.mkdir(parents=True, exist_ok=True)
    
    start = (run_id - 1) * chunk_size
    end = min(run_id * chunk_size, len(df))
    chunk = df.iloc[start:end]
    
    chunk_file = run_dir / f"chunk_{run_id}.txt"
    with open(chunk_file, "w") as fh:
        for _, row in chunk.iterrows():
            fh.write(f"{row['chrom']}\t{row['pos']}\t{row['ref']}\t{row['alt']}\t{row['gene']}\n")
    
    print(f"   Run {run_id}: {len(chunk):,} variants -> {chunk_file}")

print(f"\n✅ Chunks saved to {OUT_DIR}")
