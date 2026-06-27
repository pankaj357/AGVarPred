#!/usr/bin/env python3
"""Combine per-gene parquet files into per-run combined.parquet for faster scoring."""

from pathlib import Path
import pyarrow.parquet as pq
import pyarrow as pa
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT = Path(__file__).resolve().parents[3]

BENCHMARKS = [
    ("humsavar", ROOT / "external_validation/processing/features/humsavar"),
    ("gnomad_benign", ROOT / "external_validation/processing/features/gnomad_benign"),
    ("mave_independent", ROOT / "external_validation/processing/features/mave_independent"),
]


def concat_batch(tables):
    """Concatenate tables in small batches to avoid O(n^2) schema promotion."""
    if not tables:
        return None
    if len(tables) == 1:
        return tables[0]
    # Concatenate in chunks of 10
    batch_size = 10
    batches = []
    for i in range(0, len(tables), batch_size):
        chunk = tables[i:i+batch_size]
        if len(chunk) == 1:
            batches.append(chunk[0])
        else:
            batches.append(pa.concat_tables(chunk, promote_options="default"))
    # Recursively concat batches
    while len(batches) > 1:
        new_batches = []
        for i in range(0, len(batches), batch_size):
            chunk = batches[i:i+batch_size]
            if len(chunk) == 1:
                new_batches.append(chunk[0])
            else:
                new_batches.append(pa.concat_tables(chunk, promote_options="default"))
        batches = new_batches
    return batches[0]


def combine_run(run_dir):
    run_dir = Path(run_dir)
    out_path = run_dir.with_name(run_dir.name + "_combined.parquet")
    if out_path.exists():
        return (run_dir.name, "already_exists", 0)

    parquet_files = list(run_dir.rglob("*.parquet"))
    if not parquet_files:
        return (run_dir.name, "no_files", 0)

    tables = []
    for f in parquet_files:
        try:
            tables.append(pq.read_table(str(f)))
        except Exception:
            continue
    if not tables:
        return (run_dir.name, "empty", 0)

    combined = concat_batch(tables)
    pq.write_table(combined, str(out_path))
    return (run_dir.name, "done", combined.num_rows)


def combine_benchmark(name, feature_dir):
    print(f"\n{'='*70}")
    print(f"Combining features for: {name}")
    print(f"{'='*70}")
    feature_dir = Path(feature_dir)
    run_dirs = sorted([d for d in feature_dir.glob("output_test_run_*") if d.is_dir()])
    if not run_dirs:
        print("  No run directories found.")
        return

    print(f"  Found {len(run_dirs)} run directories")
    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(combine_run, rd): rd for rd in run_dirs}
        for future in as_completed(futures):
            run_name, status, n_rows = future.result()
            print(f"    {run_name}: {status} ({n_rows:,} rows)")

    print(f"  Done with {name}")


if __name__ == "__main__":
    for name, fdir in BENCHMARKS:
        combine_benchmark(name, fdir)
    print("\n" + "="*70)
    print("All benchmarks combined!")
    print("="*70)
