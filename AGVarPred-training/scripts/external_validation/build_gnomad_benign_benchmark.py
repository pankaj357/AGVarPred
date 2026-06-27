#!/usr/bin/env python3
"""
Build a gnomAD Common Benign benchmark restricted to genes NOT in training.

Strategy:
- Source: gnomAD exomes liftover VCF (GRCh38)
- For each gene not present in train/cal/test, query the VCF for common SNVs
  (AF > 1%) within the transcript region.
- All labels = 0 (genuinely benign — purifying selection removes pathogenic
  common variants from the population).
- Strict variant deduplication against train/cal/test.

Note: The liftover gnomAD VCF does not contain VEP consequence annotations.
We therefore include all common exonic SNVs; downstream VEP annotation will
classify them. For scoring, one may optionally filter to missense only.
"""

import os
import subprocess
import gzip
import pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
EXT = ROOT / "external_validation"
GNOMAD_VCF = os.environ.get(
    "GNOMAD_VCF",
    "/data/kanaka/pankaj/Raw_Dataset/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
)
REFGENE = EXT / "raw_data" / "refGene.txt.gz"

# Sampling limits to keep runtime reasonable
MAX_GENES = 800          # Number of non-training genes to sample
MAX_PER_GENE = 15        # Max benign variants to keep per gene
MIN_AF = 0.001           # AF > 0.1% = common benign


def load_training_genes_and_variants():
    print("🔄 Loading train/cal/test for holdout...")
    seen_vars = set()
    seen_genes = set()
    for split in ["train.csv", "cal.csv", "test.csv"]:
        df = pd.read_csv(ROOT / split, dtype=str)
        df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False)
        df["PositionVCF"] = df["PositionVCF"].astype(str).str.replace(".0", "", regex=False)
        for _, row in df.iterrows():
            key = (row["Chromosome"], row["PositionVCF"],
                   str(row["ReferenceAlleleVCF"]).upper(),
                   str(row["AlternateAlleleVCF"]).upper())
            seen_vars.add(key)
        seen_genes.update(df["GeneSymbol"].dropna().unique())
    print(f"  Training genes: {len(seen_genes):,}")
    print(f"  Training variants: {len(seen_vars):,}")
    return seen_vars, seen_genes


def load_refgene():
    print("📚 Loading refGene coordinates...")
    # refGene columns: bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames
    cols = [0, 1, 2, 4, 5, 12]
    names = ["name", "chrom", "txStart", "txEnd", "gene"]
    df = pd.read_csv(REFGENE, sep="\t", header=None, usecols=cols,
                     names=["bin", "name", "chrom", "txStart", "txEnd", "gene"],
                     dtype={"chrom": str, "txStart": int, "txEnd": int, "gene": str})
    # Keep "chr" prefix for tabix compatibility
    df["chrom"] = df["chrom"].astype(str)
    # Keep canonical / longest transcript per gene
    df = df.sort_values("txEnd", ascending=False).drop_duplicates(subset=["gene"], keep="first")
    print(f"  refGene genes: {len(df):,}")
    return df


def parse_vcf_line(line):
    cols = line.strip().split("\t")
    if len(cols) < 8:
        return None
    chrom, pos, _, ref, alt, _, _, info = cols[:8]
    if len(ref) != 1 or len(alt) != 1:
        return None
    af = 0.0
    for entry in info.split(";"):
        if entry.startswith("AF="):
            try:
                af = float(entry.split("=")[1])
            except ValueError:
                pass
            break
    return {
        "chrom": chrom.replace("chr", ""),
        "pos": int(pos),
        "ref": ref.upper(),
        "alt": alt.upper(),
        "af": af,
    }


def query_gnomad_by_region(chrom, start, end, max_variants=1000):
    region = f"chr{chrom.replace('chr', '')}:{start}-{end}"
    cmd = f"tabix {GNOMAD_VCF} {region}"
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
    except subprocess.TimeoutExpired:
        return []
    variants = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        v = parse_vcf_line(line)
        if v and v["af"] > MIN_AF:
            variants.append(v)
        if len(variants) >= max_variants:
            break
    return variants


def build_benchmark():
    print("=" * 70)
    print("GNOMAD COMMON BENIGN BENCHMARK BUILDER")
    print("=" * 70)

    if not Path(GNOMAD_VCF).exists():
        raise FileNotFoundError(
            f"gnomAD VCF not found: {GNOMAD_VCF}\n"
            "Set GNOMAD_VCF env var or place the VCF at the default path."
        )
    if not REFGENE.exists():
        raise FileNotFoundError(
            f"refGene not found: {REFGENE}\n"
            "Download from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
        )

    seen_vars, seen_genes = load_training_genes_and_variants()
    refgene = load_refgene()

    # Exclude training genes
    refgene = refgene[~refgene["gene"].isin(seen_genes)].copy()
    print(f"  Non-training genes in refGene: {len(refgene):,}")

    # Sample genes for feasibility
    if len(refgene) > MAX_GENES:
        refgene = refgene.sample(MAX_GENES, random_state=42)
        print(f"  Sampled {MAX_GENES} genes for querying")

    benign_rows = []
    skipped_no_data = 0
    skipped_no_common = 0

    for idx, row in refgene.iterrows():
        chrom = str(row["chrom"])
        start = max(1, row["txStart"] - 1000)
        end = row["txEnd"] + 1000
        variants = query_gnomad_by_region(chrom, start, end, max_variants=MAX_PER_GENE * 2)

        if not variants:
            skipped_no_data += 1
            continue

        # Keep up to MAX_PER_GENE common variants per gene
        kept = 0
        for v in variants:
            key = (v["chrom"], str(v["pos"]), v["ref"], v["alt"])
            if key in seen_vars:
                continue
            benign_rows.append({
                "chrom": v["chrom"],
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"],
                "gene": row["gene"],
                "label": 0,
                "af": v["af"],
                "dataset": "gnomad_common_benign",
            })
            kept += 1
            if kept >= MAX_PER_GENE:
                break

        if kept == 0:
            skipped_no_common += 1

        if (idx + 1) % 100 == 0 or idx == len(refgene) - 1:
            print(f"  Processed {min(idx + 1, len(refgene))}/{len(refgene)} genes | "
                  f"Variants: {len(benign_rows):,} | No data: {skipped_no_data} | No common: {skipped_no_common}")

    df = pd.DataFrame(benign_rows)
    if len(df) == 0:
        print("⚠️  No benign variants found.")
        return

    # Deduplicate exact variants
    before = len(df)
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    print(f"\n  Before exact dedup: {before:,}")
    print(f"  After exact dedup:  {len(df):,}")

    out = EXT / "benchmarks/benchmark_gnomad_benign.csv"
    df.to_csv(out, index=False)
    print(f"\n  ✅ Saved {len(df):,} variants to {out}")
    print(f"  Genes: {df['gene'].nunique():,}")
    print(f"  Mean AF: {df['af'].mean():.4f}")


if __name__ == "__main__":
    build_benchmark()
