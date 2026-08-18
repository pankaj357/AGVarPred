#!/usr/bin/env python3
"""
Build VIP (Variant Interpretation for clinical Practice) benchmark.

Source: BCM-HGSC Neptune VIP database (GitLab)
URL: https://gitlab.com/bcm-hgsc/neptune/-/raw/master/data/vip/2020-10-13.vip.txt

Processing:
  1. Parse GRCh37 VCF, keep SNVs only
  2. Filter for genes NOT in train/cal/test (strict gene holdout)
  3. Assign binary labels from VIP curation fields:
     - Pathogenic (1): eMerge_category in {Path, Lpath, Path_NR}
       OR (eMerge_category == 'review' AND ClinVar Path/LP)
     - Benign (0): eMerge_category in {B, LB}
  4. Lift over GRCh37 -> GRCh38 using pyliftover
  5. Deduplicate against train/cal/test by chrom/pos/ref/alt
  6. Save benchmark CSV

Output:
  external_validation/benchmarks/benchmark_vip.csv
"""

import os
import re
import gzip
from pathlib import Path
from collections import Counter

import pandas as pd
from pyliftover import LiftOver

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
EXT = ROOT / "external_validation"
VIP_VCF = EXT / "source_datasets/vip/vip_2020-10-13.vcf"
BENCH_DIR = EXT / "benchmarks"
BENCH_DIR.mkdir(exist_ok=True)

TRAIN_CSV = ROOT / "train.csv"
CAL_CSV = ROOT / "cal.csv"
TEST_CSV = ROOT / "test.csv"

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
        df = pd.read_csv(csv_path, dtype=str)
        df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False)
        df["PositionVCF"] = df["PositionVCF"].astype(str).str.replace(r"\.0$", "", regex=True)
        genes.update(df["GeneSymbol"].dropna().unique())
        for _, row in df.iterrows():
            key = (row["Chromosome"],
                   row["PositionVCF"],
                   str(row["ReferenceAlleleVCF"]).upper(),
                   str(row["AlternateAlleleVCF"]).upper())
            variants.add(key)
        print(f"   {name}: {len(df):,} variants | {df['GeneSymbol'].nunique():,} genes")
    print(f"   Total training genes: {len(genes):,}")
    print(f"   Total training variants: {len(variants):,}")
    return genes, variants


# ---------------------------------------------------------------------------
# 2. PARSE VIP VCF AND EXTRACT CANDIDATES
# ---------------------------------------------------------------------------
def parse_vip_vcf(training_genes):
    """Parse VIP VCF and extract labeled SNVs for non-training genes."""
    print("\n" + "=" * 70)
    print("📖 Parsing VIP VCF and extracting candidates")
    print("=" * 70)

    rows = []
    n_total = 0
    n_snv = 0
    n_not_train = 0
    n_labeled = 0

    pathogenic_emerge = {"Path", "Lpath", "Path_NR"}
    benign_emerge = {"B", "LB"}

    with open(VIP_VCF, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_total += 1
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, id_, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]

            # Keep SNVs only
            if len(ref) != 1 or len(alt) != 1:
                continue
            n_snv += 1

            info = parts[7]

            # Extract Gene
            gene_match = re.search(r"Gene=([^;]+)", info)
            if not gene_match:
                continue
            gene = gene_match.group(1)

            # Gene holdout
            if gene in training_genes:
                continue
            n_not_train += 1

            # Extract eMerge_category
            cat_match = re.search(r"eMerge_category=([^;]+)", info)
            if not cat_match:
                continue
            emerge_cat = cat_match.group(1)

            # Extract ClinVar category (if present)
            clinvar_match = re.search(r"Clinvar_cat=([^;]+)", info)
            clinvar_cat = clinvar_match.group(1) if clinvar_match else ""

            # Assign label
            label = None
            label_source = ""

            if emerge_cat in pathogenic_emerge:
                label = 1
                label_source = f"VIP_{emerge_cat}"
            elif emerge_cat == "review" and clinvar_cat and (
                "Pathogenic" in clinvar_cat or "Likely pathogenic" in clinvar_cat
            ):
                label = 1
                label_source = f"VIP_review+ClinVar_{clinvar_cat}"
            elif emerge_cat in benign_emerge:
                label = 0
                label_source = f"VIP_{emerge_cat}"

            if label is None:
                continue
            n_labeled += 1

            rows.append({
                "chrom_grch37": chrom,
                "pos_grch37": int(pos),
                "ref": ref.upper(),
                "alt": alt.upper(),
                "gene": gene,
                "label": label,
                "label_source": label_source,
                "emerge_category": emerge_cat,
                "clinvar_category": clinvar_cat,
                "info": info,
            })

    print(f"   Total VIP variants: {n_total:,}")
    print(f"   SNVs: {n_snv:,}")
    print(f"   In non-training genes: {n_not_train:,}")
    print(f"   With clear labels: {n_labeled:,}")
    print(f"      Pathogenic: {sum(1 for r in rows if r['label'] == 1):,}")
    print(f"      Benign:     {sum(1 for r in rows if r['label'] == 0):,}")

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 3. LIFT OVER GRCh37 -> GRCh38
# ---------------------------------------------------------------------------
def liftover_grch37_to_grch38(df):
    """Lift over GRCh37 coordinates to GRCh38 using pyliftover."""
    print("\n" + "=" * 70)
    print("🗺️  Lifting over GRCh37 -> GRCh38")
    print("=" * 70)

    lo = LiftOver("hg19", "hg38")

    lifted_rows = []
    n_failed = 0
    n_success = 0
    n_multiple = 0

    for _, row in df.iterrows():
        chrom = row["chrom_grch37"]
        pos = row["pos_grch37"]

        # pyliftover uses UCSC-style chromosome names with "chr" prefix
        ucsc_chrom = f"chr{chrom}"
        # VCF is 1-based; pyliftover expects 0-based coordinates
        result = lo.convert_coordinate(ucsc_chrom, pos - 1)

        if not result:
            n_failed += 1
            continue
        if len(result) > 1:
            n_multiple += 1
            # Take the first mapping
        n_success += 1

        new_chrom, new_pos_0based, strand, _ = result[0]
        # Convert back from UCSC style to plain chromosome number
        new_chrom = new_chrom.replace("chr", "")
        # Convert to 1-based VCF coordinate
        new_pos = new_pos_0based + 1

        lifted_rows.append({
            "chrom": new_chrom,
            "pos": int(new_pos),
            "ref": row["ref"],
            "alt": row["alt"],
            "gene": row["gene"],
            "label": row["label"],
            "label_source": row["label_source"],
            "emerge_category": row["emerge_category"],
            "clinvar_category": row["clinvar_category"],
            "chrom_grch37": row["chrom_grch37"],
            "pos_grch37": row["pos_grch37"],
        })

    print(f"   Lifted successfully: {n_success:,}")
    print(f"   Failed: {n_failed:,}")
    print(f"   Multiple mappings (used first): {n_multiple:,}")

    return pd.DataFrame(lifted_rows)


# ---------------------------------------------------------------------------
# 4. DEDUPLICATE AGAINST TRAIN/CAL/TEST
# ---------------------------------------------------------------------------
def deduplicate(df, train_variants):
    """Remove variants that overlap train/cal/test by chrom/pos/ref/alt."""
    print("\n" + "=" * 70)
    print("🧹 Deduplicating against train/cal/test")
    print("=" * 70)

    before = len(df)
    df["_key"] = list(zip(
        df["chrom"].astype(str),
        df["pos"].astype(str),
        df["ref"].astype(str).str.upper(),
        df["alt"].astype(str).str.upper()
    ))
    df = df[~df["_key"].isin(train_variants)].drop(columns=["_key"])
    after = len(df)
    print(f"   Before dedup: {before:,}")
    print(f"   After dedup:  {after:,} (removed {before - after:,})")
    return df


# ---------------------------------------------------------------------------
# 5. FINALIZE AND SAVE
# ---------------------------------------------------------------------------
def finalize_and_save(df):
    """Finalize benchmark and save CSV."""
    print("\n" + "=" * 70)
    print("🏁 Finalizing benchmark")
    print("=" * 70)

    # Drop duplicates on chrom/pos/ref/alt (can happen if liftover produces same coord)
    before = len(df)
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    print(f"   After exact dedup: {len(df):,} (removed {before - len(df):,})")

    # Sort
    df = df.sort_values(["chrom", "pos", "ref", "alt"]).reset_index(drop=True)

    # Summary
    n_path = (df["label"] == 1).sum()
    n_benign = (df["label"] == 0).sum()
    n_genes = df["gene"].nunique()

    print(f"\n   Final benchmark statistics:")
    print(f"      Total variants: {len(df):,}")
    print(f"      Pathogenic:     {n_path:,}")
    print(f"      Benign:         {n_benign:,}")
    print(f"      Genes:          {n_genes:,}")
    print(f"\n   Label source breakdown:")
    for src, count in df["label_source"].value_counts().items():
        print(f"      {src}: {count:,}")

    # Select output columns
    out = df[["chrom", "pos", "ref", "alt", "gene", "label", "label_source",
              "emerge_category", "clinvar_category", "chrom_grch37", "pos_grch37"]].copy()
    out["dataset"] = "vip"

    out_path = BENCH_DIR / "benchmark_vip.csv"
    out.to_csv(out_path, index=False)
    print(f"\n   ✅ Saved: {out_path}")

    # Also update variant_benchmark_map.csv for chunking
    map_path = BENCH_DIR / "variant_benchmark_map.csv"
    df_map = out[["chrom", "pos", "ref", "alt", "gene", "label"]].copy()
    df_map["variant_id"] = df_map["chrom"] + "_" + df_map["pos"].astype(str) + "_" + df_map["ref"] + "_" + df_map["alt"]
    df_map["benchmark_source"] = "vip"

    # Append to existing map if present
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
    print("🏥 VIP DATABASE BENCHMARK BUILDER")
    print("=" * 70)
    print(f"   VIP VCF: {VIP_VCF}")
    print(f"   Output:  {BENCH_DIR}/benchmark_vip.csv")

    if not VIP_VCF.exists():
        print(f"\n❌ VIP VCF not found: {VIP_VCF}")
        print("   Please run the download step first.")
        return

    training_genes, training_variants = load_training_holdout()
    candidates = parse_vip_vcf(training_genes)
    if candidates.empty:
        print("\n⚠️ No candidate variants found.")
        return

    lifted = liftover_grch37_to_grch38(candidates)
    if lifted.empty:
        print("\n⚠️ No variants survived liftover.")
        return

    clean = deduplicate(lifted, training_variants)
    if clean.empty:
        print("\n⚠️ No variants survived deduplication.")
        return

    finalize_and_save(clean)

    print("\n" + "=" * 70)
    print("✅ VIP benchmark build complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
