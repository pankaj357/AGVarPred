#!/usr/bin/env python3
"""Combine per-gene parquet files into a single combined.parquet for PanelApp benign."""

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FEATURE_DIR = ROOT / "external_validation" / "processing" / "features" / "panelapp_benign"
OUT_PATH = FEATURE_DIR / "combined.parquet"


def main():
    print(f"🔄 Combining PanelApp benign features from {FEATURE_DIR}")

    if OUT_PATH.exists():
        print(f"  {OUT_PATH} already exists, skipping.")
        return

    run_dirs = sorted(FEATURE_DIR.glob("output_test_run_*"))
    if not run_dirs:
        print("  No run directories found. Looking for parquet files directly...")
        parquet_files = list(FEATURE_DIR.rglob("*.parquet"))
        if not parquet_files:
            print("  No parquet files found!")
            return
        dataset = ds.dataset([str(f) for f in parquet_files], format="parquet")
        table = dataset.to_table()
        pq.write_table(table, str(OUT_PATH))
        print(f"  Written: {OUT_PATH} ({table.num_rows:,} rows, {table.num_columns:,} cols)")
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

    if not tables:
        print("  No tables to combine!")
        return

    print(f"  Concatenating {len(tables)} tables ({total_rows:,} total rows)...")
    combined = pa.concat_tables(tables)
    print(f"  Writing combined.parquet...")
    pq.write_table(combined, str(OUT_PATH))
    print(f"  Done: {OUT_PATH} ({combined.num_rows:,} rows, {combined.num_columns:,} cols)")


if __name__ == "__main__":
    main()
