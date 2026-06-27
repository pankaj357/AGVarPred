#!/usr/bin/env python3
"""
Combine per-gene AlphaGenome feature parquet files for Grimm2015 benchmark.
Writes external_validation/processing/features/grimm2015/combined.parquet
"""

import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FEATURE_DIR = ROOT / "external_validation" / "processing" / "features" / "grimm2015"
OUT_PATH = FEATURE_DIR / "combined.parquet"


def clean_name(s):
    import re
    s = str(s)
    s = s.replace(":", "_").replace("-", "_").replace(" ", "_")
    s = re.sub(r'[^A-Za-z0-9_]', '', s)
    s = re.sub(r'_+', '_', s)
    return s


def main():
    print("=" * 60)
    print("Combining Grimm2015 AlphaGenome features")
    print("=" * 60)

    parquet_files = sorted(FEATURE_DIR.rglob("*.parquet"))
    print(f"Found {len(parquet_files):,} parquet files")

    if not parquet_files:
        raise FileNotFoundError(f"No parquet files in {FEATURE_DIR}")

    # Read in batches to avoid memory spikes
    batch_size = 200
    combined_parts = []
    total_rows = 0

    n_batches = (len(parquet_files) + batch_size - 1) // batch_size
    for i in range(0, len(parquet_files), batch_size):
        batch_num = i // batch_size + 1
        print(f"Reading parquet batch {batch_num}/{n_batches}...")
        batch_files = parquet_files[i:i + batch_size]
        batch_dfs = []
        for f in batch_files:
            try:
                df = pd.read_parquet(f)
                batch_dfs.append(df)
            except Exception as e:
                print(f"\n⚠ Error reading {f}: {e}")
        if batch_dfs:
            batch_combined = pd.concat(batch_dfs, ignore_index=False)
            combined_parts.append(batch_combined)
            total_rows += len(batch_combined)

    print(f"\nConcatenating {len(combined_parts)} batches ({total_rows:,} total rows)...")
    df = pd.concat(combined_parts, ignore_index=False)
    df.columns = [clean_name(c) for c in df.columns]
    if df.index.name is None:
        df.index.name = "variant_id"

    n_before = len(df)
    df = df[~df.index.duplicated(keep="first")]
    n_deduped = n_before - len(df)
    if n_deduped > 0:
        print(f"⚠ Deduplicated: removed {n_deduped:,} duplicate variant rows")

    print(f"Final combined: {len(df):,} variants, {len(df.columns):,} features")
    df.to_parquet(OUT_PATH)
    print(f"Saved: {OUT_PATH}")
    print("\n✅ Done!")


if __name__ == "__main__":
    main()
