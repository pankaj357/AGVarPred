#!/usr/bin/env python3
"""Combine all per-gene DVD AlphaGenome feature parquet files into one combined.parquet."""
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FEATURE_DIR = ROOT / "external_validation" / "processing" / "features" / "dvd"
OUT_PATH = FEATURE_DIR / "combined.parquet"

parquet_files = sorted(FEATURE_DIR.rglob("*_ALL_VEP_RAW_SCORE_MATRIX.parquet"))
print(f"Found {len(parquet_files):,} parquet files")

dfs = []
for f in parquet_files:
    try:
        df = pd.read_parquet(f)
        dfs.append(df)
    except Exception as e:
        print(f"⚠ Error reading {f}: {e}")

if not dfs:
    raise FileNotFoundError(f"No parquet files found in {FEATURE_DIR}")

# Align columns across all files
all_cols = sorted(set().union(*(df.columns for df in dfs)))
combined = pd.concat([df.reindex(columns=all_cols) for df in dfs], ignore_index=False)

# Deduplicate rows by variant_id (should not happen, but be safe)
n_before = len(combined)
combined = combined[~combined.index.duplicated(keep="first")]
n_deduped = n_before - len(combined)
if n_deduped:
    print(f"⚠ Deduplicated {n_deduped:,} rows")

# DVD benchmark uses chrMT for mitochondrial variants; AlphaGenome normalizes to chrM.
combined.index = combined.index.str.replace(r"^chrM_", "chrMT_", regex=True)
combined.index.name = "variant_id"
combined.to_parquet(OUT_PATH)
print(f"Saved combined features: {OUT_PATH}")
print(f"  Variants: {len(combined):,}")
print(f"  Features: {len(combined.columns):,}")
