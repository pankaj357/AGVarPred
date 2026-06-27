#!/usr/bin/env python3
"""
Build DVD (Deafness Variation Database) benchmark.

Source: https://deafnessvariationdatabase.org/
Processing:
  1. Download per-gene CSVs (run download_dvd_variants_parallel.py first)
  2. Parse GRCh38 genomic description chrom:pos:ref>alt
  3. Filter for Pathogenic/Likely pathogenic (label=1) and Benign/Likely benign (label=0)
  4. Filter for genes NOT in train/cal/test (strict gene holdout)
  5. Deduplicate against train/cal/test by chrom/pos/ref/alt
  6. Save benchmark CSV

Output:
  external_validation/benchmarks/benchmark_dvd.csv
"""
import os
import re
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
EXT = ROOT / "external_validation"
DVD_CSV_DIR = EXT / "dvd_data/csv"
BENCH_DIR = EXT / "benchmarks"
BENCH_DIR.mkdir(exist_ok=True)

TRAIN_CSV = ROOT / "train.csv"
CAL_CSV = ROOT / "cal.csv"
TEST_CSV = ROOT / "test.csv"

PATHOLOGIC = {"Pathogenic", "Likely pathogenic"}
BENIGN = {"Benign", "Likely benign"}

# ---------------------------------------------------------------------------
# 1. LOAD TRAINING DATA FOR HOLDOUT
# ---------------------------------------------------------------------------
def load_training_holdout():
    """Load training genes and variant coordinates for strict holdout."""
    print("=" * 70)
    print("🔄 Loading train/cal/test for gene + variant holdout")
    print("=" * 70)

    genes = set()
    variants = set()
    for csv_path, name in [(TRAIN_CSV, "train"), (CAL_CSV, "cal"), (TEST_CSV, "test")]:
        df = pd.read_csv(csv_path, dtype=str, low_memory=False)
        df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False)
        df["PositionVCF"] = df["PositionVCF"].astype(str).str.replace(r"\.0$", "", regex=True)
        genes.update(df["GeneSymbol"].dropna().unique())
        for _, row in df.iterrows():
            key = (
                str(row["Chromosome"]),
                str(row["PositionVCF"]),
                str(row["ReferenceAlleleVCF"]).upper(),
                str(row["AlternateAlleleVCF"]).upper(),
            )
            variants.add(key)
        print(f"   {name}: {len(df):,} variants | {df['GeneSymbol'].nunique():,} genes")
    print(f"   Total training genes: {len(genes):,}")
    print(f"   Total training variants: {len(variants):,}")
    return genes, variants


# ---------------------------------------------------------------------------
# 2. PARSE DVD CSVs
# ---------------------------------------------------------------------------
def parse_genomic_description(desc):
    """Parse 'chr:pos:ref>alt' into (chrom, pos, ref, alt)."""
    if pd.isna(desc) or not desc:
        return None
    desc = str(desc).strip()
    # Format: 13:20187466:T>A or 13:20187468:AATT>A
    m = re.match(r"^([^:]+):(\d+):([^>]+)>(.*)$", desc)
    if not m:
        return None
    chrom, pos, ref, alt = m.groups()
    # Normalize chromosome
    chrom = chrom.replace("chr", "").replace("Chr", "")
    return chrom, int(pos), ref.upper(), alt.upper()


def load_dvd_csvs():
    """Load and concatenate all per-gene DVD CSVs."""
    print("\n" + "=" * 70)
    print("📖 Loading DVD per-gene CSVs")
    print("=" * 70)

    if not DVD_CSV_DIR.exists():
        raise FileNotFoundError(f"DVD CSV directory not found: {DVD_CSV_DIR}")

    csv_files = sorted(DVD_CSV_DIR.glob("*.csv"))
    print(f"   Found {len(csv_files)} CSV files")

    frames = []
    for p in csv_files:
        try:
            df = pd.read_csv(p, dtype=str, low_memory=False)
            df["_source_file"] = p.stem
            frames.append(df)
        except Exception as e:
            print(f"   ⚠️ Could not read {p}: {e}")

    combined = pd.concat(frames, ignore_index=True)
    print(f"   Total rows loaded: {len(combined):,}")
    return combined


def filter_and_label(df):
    """Assign binary labels and parse coordinates."""
    print("\n" + "=" * 70)
    print("🏷️  Filtering classifications and parsing coordinates")
    print("=" * 70)

    # Classification column can have values like "Pathogenic", "Likely pathogenic", etc.
    df["classification"] = df["Variant_Classification"].astype(str).str.strip()

    rows = []
    n_total = len(df)
    n_parsed = 0
    n_labeled = 0
    n_path = 0
    n_benign = 0

    for _, row in df.iterrows():
        cls = row["classification"]
        if cls in PATHOLOGIC:
            label = 1
            n_path += 1
        elif cls in BENIGN:
            label = 0
            n_benign += 1
        else:
            continue

        parsed = parse_genomic_description(row.get("Genomic_Description_GRCh38"))
        if parsed is None:
            continue
        n_parsed += 1

        chrom, pos, ref, alt = parsed
        gene = row.get("_source_file", "")  # filename is gene symbol

        rows.append({
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "gene": gene,
            "label": label,
            "classification": cls,
            "phenotype": row.get("Phenotype", ""),
            "consequence": row.get("Consequence", ""),
            "cadd_phred": row.get("CADD_Phred", ""),
            "revel": row.get("REVEL", ""),
        })
        n_labeled += 1

    print(f"   Total rows: {n_total:,}")
    print(f"   With P/LP or B/LB classification: {n_path + n_benign:,}")
    print(f"   Parsed GRCh38 coordinates: {n_parsed:,}")
    print(f"   Pathogenic: {n_path:,}")
    print(f"   Benign: {n_benign:,}")

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 3. GENE HOLDOUT
# ---------------------------------------------------------------------------
def apply_gene_holdout(df, training_genes):
    """Remove variants in genes seen during training/calibration/testing."""
    print("\n" + "=" * 70)
    print("🧬 Applying strict gene holdout")
    print("=" * 70)
    before = len(df)
    df = df[~df["gene"].isin(training_genes)].copy()
    after = len(df)
    removed_genes = before - after
    print(f"   Before: {before:,}")
    print(f"   After:  {after:,} (removed {removed_genes:,} variants)")
    print(f"   Remaining genes: {df['gene'].nunique()}")
    return df


# ---------------------------------------------------------------------------
# 4. VARIANT DEDUPLICATION
# ---------------------------------------------------------------------------
def deduplicate(df, training_variants):
    """Remove exact coordinate overlaps with train/cal/test."""
    print("\n" + "=" * 70)
    print("🧹 Deduplicating against train/cal/test variants")
    print("=" * 70)
    before = len(df)
    df["_key"] = list(zip(
        df["chrom"].astype(str),
        df["pos"].astype(str),
        df["ref"].astype(str).str.upper(),
        df["alt"].astype(str).str.upper(),
    ))
    df = df[~df["_key"].isin(training_variants)].drop(columns=["_key"])
    after = len(df)
    print(f"   Before: {before:,}")
    print(f"   After:  {after:,} (removed {before - after:,})")
    return df


# ---------------------------------------------------------------------------
# 5. FINALIZE AND SAVE
# ---------------------------------------------------------------------------
def finalize_and_save(df):
    """Finalize benchmark and save CSV."""
    print("\n" + "=" * 70)
    print("🏁 Finalizing benchmark")
    print("=" * 70)

    # Drop duplicates on chrom/pos/ref/alt (keep first label)
    before = len(df)
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    print(f"   After exact dedup: {len(df):,} (removed {before - len(df):,})")

    df = df.sort_values(["chrom", "pos", "ref", "alt"]).reset_index(drop=True)

    n_path = (df["label"] == 1).sum()
    n_benign = (df["label"] == 0).sum()
    n_genes = df["gene"].nunique()

    print(f"\n   Final benchmark statistics:")
    print(f"      Total variants: {len(df):,}")
    print(f"      Pathogenic:     {n_path:,}")
    print(f"      Benign:         {n_benign:,}")
    print(f"      Genes:          {n_genes:,}")

    out_cols = [
        "chrom", "pos", "ref", "alt", "gene", "label",
        "classification", "phenotype", "consequence", "cadd_phred", "revel"
    ]
    out = df[out_cols].copy()
    out["dataset"] = "dvd"

    out_path = BENCH_DIR / "benchmark_dvd.csv"
    out.to_csv(out_path, index=False)
    print(f"\n   ✅ Saved: {out_path}")

    # Update variant_benchmark_map.csv
    map_path = BENCH_DIR / "variant_benchmark_map.csv"
    df_map = out[["chrom", "pos", "ref", "alt", "gene", "label"]].copy()
    df_map["variant_id"] = (
        df_map["chrom"].astype(str) + "_" +
        df_map["pos"].astype(str) + "_" +
        df_map["ref"].astype(str) + "_" +
        df_map["alt"].astype(str)
    )
    df_map["benchmark_source"] = "dvd"

    if map_path.exists():
        existing = pd.read_csv(map_path)
        combined = pd.concat([existing, df_map], ignore_index=True)
        combined = combined.drop_duplicates(subset=["variant_id"])
        combined.to_csv(map_path, index=False)
        print(f"   ✅ Appended {len(df_map):,} variants to {map_path}")
    else:
        df_map.to_csv(map_path, index=False)
        print(f"   ✅ Created {map_path}")

    return out


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("🧬 DVD BENCHMARK BUILDER")
    print("=" * 70)
    print(f"   CSV dir: {DVD_CSV_DIR}")
    print(f"   Output:  {BENCH_DIR}/benchmark_dvd.csv")

    training_genes, training_variants = load_training_holdout()
    raw = load_dvd_csvs()
    if raw.empty:
        print("\n⚠️ No DVD CSVs found. Please run download step first.")
        return

    labeled = filter_and_label(raw)
    if labeled.empty:
        print("\n⚠️ No labeled variants found.")
        return

    heldout_genes = apply_gene_holdout(labeled, training_genes)
    if heldout_genes.empty:
        print("\n⚠️ No variants survived gene holdout.")
        return

    clean = deduplicate(heldout_genes, training_variants)
    if clean.empty:
        print("\n⚠️ No variants survived variant deduplication.")
        return

    finalize_and_save(clean)

    print("\n" + "=" * 70)
    print("✅ DVD benchmark build complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
