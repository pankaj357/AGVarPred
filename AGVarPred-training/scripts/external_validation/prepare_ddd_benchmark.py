#!/usr/bin/env python3
"""
Prepare DDD Clinical Benchmark (Kaplanis/Gunning et al., J Med Genet 2021)

Paper: "Assessing performance of pathogenicity predictors using clinically
        relevant variant datasets"
PMC ID: PMC8327323
PubMed: 32843488

Supplementary Table S1 contains 1,757 clinically curated variants:
  - DDD study: 687 pathogenic + 533 benign
  - Diagnostic lab: 452 pathogenic + 28 benign  
  - Amish population: 57 benign

All coordinates are GRCh37/hg19 (in HGVS strings); lifted to GRCh38.
Labels: BENIGN = 0, PATHOGENIC = 1
"""

import os
import sys
from pathlib import Path
from urllib.request import urlretrieve
from urllib.error import HTTPError

import pandas as pd
try:
    from pyliftover import LiftOver
except ImportError:
    LiftOver = None

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]          # external_validation/
PROJECT_ROOT = ROOT.parent                           # project root (train.csv)
RAW_DIR = ROOT / "raw_data"
BENCHMARK_DIR = ROOT / "benchmarks"
RAW_DIR.mkdir(parents=True, exist_ok=True)
BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)

XLSX_PATH = RAW_DIR / "ddd_clinical_benchmark.xlsx"
OUTPUT_CSV = BENCHMARK_DIR / "benchmark_ddd_clinical.csv"

# ── Download URLs ────────────────────────────────────────────────────────────
# NCBI PMC often blocks curl; manual download may be required.
URLS = [
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8327323/bin/jmedgenet-2020-107003supp001.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC8327323/bin/jmedgenet-2020-107003supp001.xlsx",
]

MANUAL_INSTRUCTIONS = """
╔══════════════════════════════════════════════════════════════════════════════╗
║  AUTOMATIC DOWNLOAD FAILED                                                   ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Please download the DDD Clinical Benchmark supplement manually:             ║
║                                                                              ║
║  1. Go to: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8327323/             ║
║  2. Scroll to "Supplementary Material"                                       ║
║  3. Download:                                                                ║
║     "jmedgenet-2020-107003supp001.xlsx"                                      ║
║     (Supplementary data - clinical dataset)                                  ║
║                                                                              ║
║  4. Place the file at:                                                       ║
║     {xlsx_path}                                                              ║
║                                                                              ║
║  5. Re-run this script.                                                      ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""


def download() -> bool:
    """Try to download automatically; return True on success."""
    if XLSX_PATH.exists():
        print(f"✅ Found existing xlsx: {XLSX_PATH}")
        return True

    for url in URLS:
        print(f"⬇️  Trying: {url}")
        try:
            urlretrieve(url, XLSX_PATH)
            # Sanity check: xlsx files start with PK
            with open(XLSX_PATH, "rb") as f:
                header = f.read(2)
            if header == b"PK":
                print(f"✅ Downloaded successfully")
                return True
            else:
                print(f"   ⚠️  Returned HTML instead of xlsx (likely blocked)")
                XLSX_PATH.unlink()
        except HTTPError as e:
            print(f"   ❌ HTTP {e.code}")
        except Exception as e:
            print(f"   ❌ {e}")
    return False


def load_training_genes_and_variants():
    """Load genes and variants from train/cal/test for strict holdout."""
    print("🔄 Loading train/cal/test for gene + variant holdout...")
    seen_vars = set()
    seen_genes = set()

    for split in ["train.csv", "cal.csv", "test.csv"]:
        path = PROJECT_ROOT / split
        if not path.exists():
            print(f"   ⚠️  {split} not found at {path}; skipping")
            continue
        df = pd.read_csv(path, dtype=str, low_memory=False)
        df["Chromosome"] = (
            df["Chromosome"].astype(str).str.replace("chr", "", case=False)
        )
        df["PositionVCF"] = (
            df["PositionVCF"].astype(str).str.replace(r"\.0$", "", regex=True)
        )
        for _, row in df.iterrows():
            key = (
                str(row["Chromosome"]),
                str(row["PositionVCF"]),
                str(row["ReferenceAlleleVCF"]).upper(),
                str(row["AlternateAlleleVCF"]).upper(),
            )
            seen_vars.add(key)
        seen_genes.update(df["GeneSymbol"].dropna().unique())

    print(f"   Training genes:    {len(seen_genes):,}")
    print(f"   Training variants: {len(seen_vars):,}")
    return seen_vars, seen_genes


def annotate_genes_from_refgene(df: pd.DataFrame) -> pd.DataFrame:
    """Map variant positions to gene symbols using refGene exons via bedtools."""
    print("📚 Annotating genes from refGene exons...")
    refgene_path = RAW_DIR / "refGene_hg19.txt.gz"
    if not refgene_path.exists():
        print(f"   ⚠️  refGene not found at {refgene_path}; skipping gene annotation")
        df["gene"] = "UNKNOWN"
        return df

    import subprocess
    import tempfile

    # 1. Build variant BED (0-based)
    var_bed = []
    for idx, row in df.iterrows():
        chrom = str(row["chrom"]).replace("chr", "")
        pos = int(row["pos"])
        var_bed.append(f"{chrom}\t{pos-1}\t{pos}\t{idx}")

    # 2. Build exon BED from refGene
    exon_bed_lines = []
    with subprocess.Popen(
        ["zcat", str(refgene_path)], stdout=subprocess.PIPE, text=True
    ) as proc:
        for line in proc.stdout:
            cols = line.strip().split("\t")
            if len(cols) < 16:
                continue
            chrom = cols[2].replace("chr", "")
            gene = cols[12]
            starts = cols[9].rstrip(",").split(",")
            ends = cols[10].rstrip(",").split(",")
            for s, e in zip(starts, ends):
                try:
                    s_int = int(s)
                    e_int = int(e)
                    if e_int > s_int:
                        exon_bed_lines.append(f"{chrom}\t{s_int}\t{e_int}\t{gene}")
                except ValueError:
                    continue

    # 3. Write temp files and run bedtools intersect
    with tempfile.TemporaryDirectory() as tmpdir:
        var_path = Path(tmpdir) / "variants.bed"
        exon_path = Path(tmpdir) / "exons.bed"
        out_path = Path(tmpdir) / "intersect.bed"

        var_path.write_text("\n".join(var_bed) + "\n")
        exon_path.write_text("\n".join(exon_bed_lines) + "\n")

        subprocess.run(
            ["bedtools", "intersect", "-a", str(var_path), "-b", str(exon_path), "-wa", "-wb"],
            stdout=open(out_path, "w"), check=True
        )

        idx_to_gene = {}
        with open(out_path) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 8:
                    idx = int(cols[3])
                    gene = cols[7]
                    if idx not in idx_to_gene:
                        idx_to_gene[idx] = gene

    genes = []
    annotated = 0
    for idx in range(len(df)):
        gene = idx_to_gene.get(idx, "UNKNOWN")
        genes.append(gene)
        if gene != "UNKNOWN":
            annotated += 1

    df["gene"] = genes
    print(f"   Annotated: {annotated:,} / {len(df):,}")
    return df



def liftover_hg19_to_hg38(df: pd.DataFrame) -> pd.DataFrame:
    """Liftover chrom/pos from hg19 to GRCh38 using pyliftover."""
    print("🔄 Lifting over hg19 → GRCh38...")
    lo = LiftOver("hg19", "hg38")

    new_chroms = []
    new_poss = []
    lifted = 0
    failed = 0

    for _, row in df.iterrows():
        chrom = str(row["chrom"]).replace("chr", "")
        pos = int(row["pos"])
        # pyliftover uses 0-based coordinates
        result = lo.convert_coordinate(f"chr{chrom}", pos - 1)
        if result and len(result) > 0:
            new_chrom = result[0][0].replace("chr", "")
            new_pos = int(result[0][1]) + 1  # back to 1-based
            new_chroms.append(new_chrom)
            new_poss.append(new_pos)
            lifted += 1
        else:
            new_chroms.append(None)
            new_poss.append(None)
            failed += 1

    df["chrom"] = new_chroms
    df["pos"] = new_poss
    df = df.dropna(subset=["chrom", "pos"]).copy()
    df["pos"] = df["pos"].astype(int)
    print(f"   Lifted: {lifted:,} | Failed: {failed:,}")
    return df


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    # 1. Download (or use existing)
    if not download():
        print(MANUAL_INSTRUCTIONS.format(xlsx_path=XLSX_PATH))
        sys.exit(1)

    # 2. Load xlsx
    print(f"📖 Reading DDD Clinical Benchmark...")
    df = pd.read_excel(XLSX_PATH, dtype=str)
    print(f"   Total rows: {len(df):,}")
    print(f"   Columns: {list(df.columns)}")

    # 3. Parse #Variant column (format: Chr1(GRCh37):g.1167764C>A)
    import re
    print("   Parsing #Variant column...")
    parsed = df["#Variant"].str.extract(r"[Cc]hr([0-9XYM]+)\(GRCh37\):g\.(\d+)([ACGTN])>([ACGTN])")
    parsed.columns = ["chrom", "pos", "ref", "alt"]
    parsed["chrom"] = parsed["chrom"].astype(str)
    parsed["pos"] = pd.to_numeric(parsed["pos"], errors="coerce")
    parsed["ref"] = parsed["ref"].str.upper()
    parsed["alt"] = parsed["alt"].str.upper()

    # Drop rows that didn't parse
    before = len(df)
    df = df.join(parsed)
    df = df.dropna(subset=["chrom", "pos", "ref", "alt"]).copy()
    print(f"   Parsed variants: {len(df):,} / {before:,}")

    # 4. Standardize pathogenicity
    df["pathogenicity"] = df["PATHOGENICITY"].str.strip().str.upper()
    df["source"] = df["GROUP"].str.strip()

    # 5. Annotate genes from refGene (GRCh37/hg19 coordinates)
    df = annotate_genes_from_refgene(df)

    # 5. Filter to BENIGN / PATHOGENIC only
    if "pathogenicity" in df.columns:
        df = df[df["pathogenicity"].isin(["BENIGN", "PATHOGENIC"])].copy()
        df["label"] = df["pathogenicity"].map({"BENIGN": 0, "PATHOGENIC": 1})
    else:
        print("   ⚠️  No pathogenicity column found; cannot filter B/P")
        sys.exit(1)
    print(f"   After B/P filter: {len(df):,}")

    # 6. Clean coordinates
    df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", case=False)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    df = df.dropna(subset=["chrom", "pos", "ref", "alt", "gene"]).copy()

    # 7. Liftover GRCh37/hg19 → GRCh38
    df = liftover_hg19_to_hg38(df)

    # 8. Load training holdout sets
    train_vars, train_genes = load_training_genes_and_variants()

    # 8. Gene holdout
    before = len(df)
    df = df[~df["gene"].isin(train_genes)].copy()
    print(f"   After gene holdout:   {len(df):,} (removed {before - len(df):,})")

    # 9. Variant holdout
    before = len(df)
    df["_key"] = list(zip(df["chrom"], df["pos"].astype(str), df["ref"], df["alt"]))
    df = df[~df["_key"].isin(train_vars)].drop(columns=["_key"]).copy()
    print(f"   After variant holdout: {len(df):,} (removed {before - len(df):,})")

    # 10. Final formatting
    df = df.reset_index(drop=True)
    out_df = pd.DataFrame({
        "chrom": df["chrom"],
        "pos": df["pos"].astype(int),
        "ref": df["ref"],
        "alt": df["alt"],
        "gene": df["gene"],
        "label": df["label"].astype(int),
        "source": df.get("source", ""),
        "dataset": "ddd_clinical",
    })

    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n✅ Benchmark saved: {OUTPUT_CSV}")
    print(f"   Total variants: {len(out_df):,}")
    print(f"   Pathogenic (1): {(out_df['label'] == 1).sum():,}")
    print(f"   Benign     (0): {(out_df['label'] == 0).sum():,}")
    print(f"   Genes: {out_df['gene'].nunique():,}")


if __name__ == "__main__":
    main()
