#!/usr/bin/env python3
"""
Fix external validation benchmarks:
1. Remove cross-benchmark variant overlaps using precedence rules
2. Rebuild variant_benchmark_map.csv
"""

import json
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BENCH_DIR = ROOT / "external_validation" / "benchmarks"
BACKUP_DIR = BENCH_DIR / "original_backup"

# Precedence: higher number = higher priority (keep in this benchmark)
PRECEDENCE = {
    "gnomad": 1,
    "mave": 2,
    "humsavar": 3,
}

BENCHMARKS = {
    "gnomad": ("benchmark_gnomad_benign.csv", "label"),
    "humsavar": ("benchmark_independent_humsavar.csv", "label"),
    "mave": ("benchmark_mave_independent.csv", "label"),
}


def build_variant_id(df):
    return (
        df["chrom"].astype(str) + "_"
        + df["pos"].astype(str) + "_"
        + df["ref"].astype(str).str.upper() + "_"
        + df["alt"].astype(str).str.upper()
    )


def fix_cross_benchmark_overlaps():
    """Remove variants from lower-priority benchmarks when they overlap."""
    # Load all benchmarks
    bench_data = {}
    for name, (fname, label_col) in BENCHMARKS.items():
        df = pd.read_csv(BENCH_DIR / fname, dtype=str)
        df["variant_id"] = build_variant_id(df)
        bench_data[name] = df

    # Build a unified registry: for each variant, which benchmarks does it belong to?
    registry = {}
    for name, df in bench_data.items():
        for vid in df["variant_id"]:
            if vid not in registry:
                registry[vid] = []
            registry[vid].append(name)

    # Find overlaps and decide which benchmark keeps each variant
    removed_counts = {name: 0 for name in BENCHMARKS}
    overlap_report = []

    for vid, bnames in registry.items():
        if len(bnames) > 1:
            # Keep in highest-priority benchmark
            keep = max(bnames, key=lambda x: PRECEDENCE[x])
            remove_from = [b for b in bnames if b != keep]
            overlap_report.append({
                "variant_id": vid,
                "in_benchmarks": ",".join(bnames),
                "kept_in": keep,
                "removed_from": ",".join(remove_from),
            })
            for b in remove_from:
                removed_counts[b] += 1

    # Apply removals
    for name, df in bench_data.items():
        if removed_counts[name] > 0:
            # Identify variants to keep (those not removed)
            keep_mask = ~df["variant_id"].isin([
                r["variant_id"] for r in overlap_report
                if name in r["removed_from"].split(",")
            ])
            df_filtered = df[keep_mask].copy()
            n_removed = len(df) - len(df_filtered)
            fname, _ = BENCHMARKS[name]
            # Drop variant_id before saving
            df_filtered = df_filtered.drop(columns=["variant_id"])
            df_filtered.to_csv(BENCH_DIR / fname, index=False)
            print(f"[{name}] Removed {n_removed:,} overlapping variants, {len(df_filtered):,} remaining")
        else:
            # Still clean up variant_id column
            df = df.drop(columns=["variant_id"])
            fname, _ = BENCHMARKS[name]
            df.to_csv(BENCH_DIR / fname, index=False)
            print(f"[{name}] No overlaps removed, {len(df):,} variants")

    # Save overlap report
    overlap_df = pd.DataFrame(overlap_report)
    overlap_df.to_csv(BENCH_DIR / "cross_benchmark_overlap_report.csv", index=False)
    print(f"\n[cross-benchmark] Total overlapping variants: {len(overlap_df):,}")
    print(f"[cross-benchmark] Overlap report saved to cross_benchmark_overlap_report.csv")

    return overlap_df


def rebuild_variant_map():
    """Rebuild variant_benchmark_map.csv with all benchmarks."""
    rows = []
    for name, (fname, label_col) in BENCHMARKS.items():
        df = pd.read_csv(BENCH_DIR / fname, dtype=str)
        df["variant_id"] = build_variant_id(df)
        for _, row in df.iterrows():
            rows.append({
                "chrom": row.get("chrom"),
                "pos": row.get("pos"),
                "ref": row.get("ref"),
                "alt": row.get("alt"),
                "gene": row.get("gene"),
                "rsid": row.get("rsid", ""),
                "label": row.get(label_col),
                "variant_id": row["variant_id"],
                "benchmark_source": name,
            })

    map_df = pd.DataFrame(rows)
    map_df.to_csv(BENCH_DIR / "variant_benchmark_map.csv", index=False)
    print(f"\n[map] Rebuilt variant_benchmark_map.csv with {len(map_df):,} variants")
    print(f"[map] Benchmarks covered: {map_df['benchmark_source'].value_counts().to_dict()}")
    return map_df


def verify():
    """Run verification checks."""
    print("\n" + "="*70)
    print("VERIFICATION")
    print("="*70)

    # 1. Check duplicates in each benchmark
    for name, (fname, label_col) in BENCHMARKS.items():
        df = pd.read_csv(BENCH_DIR / fname, dtype=str)
        df["variant_id"] = build_variant_id(df)
        dups = df.duplicated(subset=["variant_id"], keep=False).sum()
        print(f"[{name}] Duplicates: {dups} | Total: {len(df)}")

    # 2. Check cross-benchmark overlaps
    sets = {}
    for name, (fname, label_col) in BENCHMARKS.items():
        df = pd.read_csv(BENCH_DIR / fname, dtype=str)
        df["variant_id"] = build_variant_id(df)
        sets[name] = set(df["variant_id"])

    names = list(sets.keys())
    total_overlaps = 0
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            overlap = sets[names[i]] & sets[names[j]]
            if overlap:
                total_overlaps += len(overlap)
                print(f"[OVERLAP] {names[i]} <-> {names[j]}: {len(overlap)} variants")
    if total_overlaps == 0:
        print("[OK] Zero cross-benchmark overlaps")

    # 3. Check map
    map_df = pd.read_csv(BENCH_DIR / "variant_benchmark_map.csv")
    dups = map_df.duplicated(subset=["variant_id", "benchmark_source"], keep=False).sum()
    print(f"[map] Duplicates (variant+source): {dups}")
    print(f"[map] Benchmarks: {map_df['benchmark_source'].value_counts().to_dict()}")


if __name__ == "__main__":
    print("="*70)
    print("FIXING EXTERNAL VALIDATION BENCHMARKS")
    print("="*70)

    print("\n--- Step 1: Fix cross-benchmark overlaps ---")
    fix_cross_benchmark_overlaps()

    print("\n--- Step 2: Rebuild variant_benchmark_map.csv ---")
    rebuild_variant_map()

    print("\n--- Step 3: Verify ---")
    verify()

    print("\n" + "="*70)
    print("DONE")
    print("="*70)
