#!/usr/bin/env python3
"""
Build a PanelApp-based common benign benchmark.

Strategy:
- Use PanelApp expert-curated genes (green/amber confidence) that are NOT in training.
- For each such gene, query gnomAD exomes liftover VCF (GRCh38) for common SNVs
  (AF > 0.1%) within the transcript region.
- All labels = 0 (genuinely benign — common variants in disease genes are under
  strong purifying selection and would be removed if pathogenic).
- Strict variant deduplication against train/cal/test.

This benchmark is clinically relevant because it tests benign variant
classification specifically in expert-curated disease genes.
"""

import os
import subprocess
import gzip
import json
import pandas as pd
from pathlib import Path
from time import time

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
EXT = ROOT / "external_validation"
PANELAPP_GENES_PATH = EXT / "panelapp_data" / "panelapp_novel_genes.csv"
TRAIN_GENES_PATH = EXT / "train_genes_upper.json"
GNOMAD_VCF = os.environ.get(
    "GNOMAD_VCF",
    "/data/kanaka/pankaj/Raw_Dataset/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
)
REFGENE = EXT / "raw_data" / "refGene.txt.gz"
OUTDIR = EXT / "benchmarks"
OUTDIR.mkdir(exist_ok=True, parents=True)

# Parameters
MIN_CONFIDENCE = 3           # PanelApp confidence level: 3=green, 2=amber, 1=red, 0=unknown
MAX_PER_GENE = 15            # Max benign variants per gene
MIN_AF = 0.001               # AF > 0.1%
INCLUDE_MISSENSE_ONLY = False  # If True, rely on downstream VEP to filter missense


def load_training_variants():
    """Load train/cal/test variants for deduplication."""
    print("🔄 Loading train/cal/test variants for deduplication...")
    seen_vars = set()
    for split in ["train.csv", "cal.csv", "test.csv"]:
        path = ROOT / split
        if not path.exists():
            print(f"   ⚠ {path} not found, skipping")
            continue
        df = pd.read_csv(path, dtype=str)
        df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False)
        df["PositionVCF"] = df["PositionVCF"].astype(str).str.replace(".0", "", regex=False)
        for _, row in df.iterrows():
            key = (
                row["Chromosome"],
                row["PositionVCF"],
                str(row["ReferenceAlleleVCF"]).upper(),
                str(row["AlternateAlleleVCF"]).upper(),
            )
            seen_vars.add(key)
    print(f"   Training variants: {len(seen_vars):,}")
    return seen_vars


def load_panelapp_genes():
    """Load PanelApp genes not in training, filtered by confidence."""
    print(f"\n🔄 Loading PanelApp novel genes (confidence >= {MIN_CONFIDENCE})...")
    df = pd.read_csv(PANELAPP_GENES_PATH)
    df = df[df["confidence_level"] >= MIN_CONFIDENCE].copy()
    genes = sorted(df["gene_symbol"].str.upper().unique())
    print(f"   PanelApp genes passing confidence filter: {len(genes):,}")
    return genes


def load_refgene_for_genes(genes):
    """Load refGene coordinates for the specified genes."""
    print("\n🔄 Loading refGene coordinates...")
    cols = ["bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
            "exonCount", "exonStarts", "exonEnds", "score", "gene", "cdsStartStat", "cdsEndStat", "exonFrames"]
    refgene = pd.read_csv(
        REFGENE, sep="\t", header=None, names=cols,
        dtype={"chrom": str, "txStart": int, "txEnd": int, "gene": str}
    )
    refgene["gene_upper"] = refgene["gene"].str.upper()
    refgene = refgene[refgene["gene_upper"].isin({g.upper() for g in genes})]
    # Keep longest transcript per gene
    refgene = refgene.sort_values("txEnd", ascending=False).drop_duplicates(subset=["gene_upper"], keep="first")
    print(f"   RefGene coverage: {refgene['gene_upper'].nunique()}/{len(genes)} genes")
    return refgene


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
    print("PANELAPP-BASED COMMON BENIGN BENCHMARK BUILDER")
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
    if not PANELAPP_GENES_PATH.exists():
        raise FileNotFoundError(
            f"PanelApp genes not found: {PANELAPP_GENES_PATH}\n"
            "Run external_validation/scripts/download_panelapp_genes.py first."
        )

    seen_vars = load_training_variants()
    panelapp_genes = load_panelapp_genes()
    refgene = load_refgene_for_genes(panelapp_genes)

    benign_rows = []
    skipped_no_data = 0
    skipped_no_common = 0
    gene_records = refgene.to_dict("records")

    t0 = time()
    for i, row in enumerate(gene_records):
        gene = row["gene_upper"]
        chrom = str(row["chrom"])
        start = max(1, row["txStart"] - 1000)
        end = row["txEnd"] + 1000
        variants = query_gnomad_by_region(chrom, start, end, max_variants=MAX_PER_GENE * 2)

        if not variants:
            skipped_no_data += 1
            continue

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
                "gene": gene,
                "label": 0,
                "af": v["af"],
                "dataset": "panelapp_benign",
            })
            kept += 1
            if kept >= MAX_PER_GENE:
                break

        if kept == 0:
            skipped_no_common += 1

        if (i + 1) % 50 == 0 or i == len(gene_records) - 1:
            elapsed = time() - t0
            print(f"   Processed {i+1}/{len(gene_records)} genes | "
                  f"Variants: {len(benign_rows):,} | No data: {skipped_no_data} | "
                  f"No common: {skipped_no_common} | Elapsed: {elapsed:.1f}s")

    df = pd.DataFrame(benign_rows)
    if len(df) == 0:
        print("⚠️  No benign variants found.")
        return

    before = len(df)
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    print(f"\n   Before exact dedup: {before:,}")
    print(f"   After exact dedup:  {len(df):,}")

    out = OUTDIR / f"benchmark_panelapp_benign_conf{MIN_CONFIDENCE}.csv"
    df.to_csv(out, index=False)
    print(f"\n✅ Saved {len(df):,} variants to {out}")
    print(f"   Genes: {df['gene'].nunique():,}")
    print(f"   Mean variants per gene: {len(df) / df['gene'].nunique():.1f}")


if __name__ == "__main__":
    build_benchmark()
