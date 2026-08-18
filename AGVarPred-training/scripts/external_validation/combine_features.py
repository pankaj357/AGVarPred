#!/usr/bin/env python3
"""Combine per-gene parquet files into a single combined.parquet for faster scoring."""

import sys
from pathlib import Path
import pyarrow.dataset as ds
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[3]

BENCHMARKS = [
    ("humsavar", ROOT / "external_validation/processing/features/humsavar"),
    ("gnomad_benign", ROOT / "external_validation/processing/features/gnomad_benign"),
    ("mave_independent", ROOT / "external_validation/processing/features/mave_independent"),
    ("grimm2015", ROOT / "external_validation/processing/features/grimm2015"),
    ("clingen", ROOT / "external_validation/processing/features/clingen"),
]


def combine_benchmark(name, feature_dir):
    print(f"\n{'='*70}")
    print(f"Combining features for: {name}")
    print(f"{'='*70}")
    feature_dir = Path(feature_dir)
    out_path = feature_dir / "combined.parquet"

    if out_path.exists():
        print(f"  combined.parquet already exists, skipping.")
        return

    # Find run directories (output_test_run_*) — these contain only parquet files
    run_dirs = sorted(feature_dir.glob("output_test_run_*"))
    if not run_dirs:
        # Fallback: search all parquet files directly
        parquet_files = list(feature_dir.rglob("*.parquet"))
        print(f"  Found {len(parquet_files):,} parquet files (no run dirs)")
        if not parquet_files:
            print(f"  No parquet files found!")
            return
        dataset = ds.dataset([str(f) for f in parquet_files], format="parquet")
        table = dataset.to_table()
        pq.write_table(table, str(out_path))
        print(f"  Written: {out_path} ({table.num_rows:,} rows, {table.num_columns:,} cols)")
        return

    print(f"  Found {len(run_dirs)} run directories")
    tables = []
    total_rows = 0
    for run_dir in run_dirs:
        parquet_files = list(run_dir.rglob("*.parquet"))
        if not parquet_files:
            continue
        dataset = ds.dataset([str(f) for f in parquet_files], format="parquet")
        table = dataset.to_table()
        tables.append(table)
        total_rows += table.num_rows
        print(f"    {run_dir.name}: {table.num_rows:,} rows")

    print(f"  Concatenating {len(tables)} tables ({total_rows:,} total rows)...")
    combined = pa.concat_tables(tables)
    print(f"  Writing combined.parquet...")
    pq.write_table(combined, str(out_path))
    print(f"  Done: {out_path} ({combined.num_rows:,} rows, {combined.num_columns:,} cols)")


if __name__ == "__main__":
    import pyarrow as pa

    for name, fdir in BENCHMARKS:
        combine_benchmark(name, fdir)

    print("\n" + "="*70)
    print("All benchmarks combined!")
    print("="*70)
