#!/usr/bin/env python3
"""
Prepare PathPredictor S1 Table as an independent external validation benchmark.

Paper: Evans et al. Genome Research 2019
       "Genetic variant pathogenicity prediction trained using disease-specific
        clinical sequencing data sets"
PMC ID: PMC6633260

S1 Table contains missense variants from clinical labs (GeneDx + LMM) with labels:
  B = Benign / Likely Benign      -> label = 0
  P = Pathogenic / Likely Pathogenic -> label = 1
  V = VUS (discarded)

Coordinates are in hg19; we liftover to GRCh38.
"""

import os
import sys
import zipfile
from pathlib import Path
from urllib.request import urlretrieve
from urllib.error import HTTPError

import pandas as pd
from pyliftover import LiftOver

# ── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]          # external_validation/
PROJECT_ROOT = ROOT.parent                           # project root (train.csv)
RAW_DIR = ROOT / "raw_data"
BENCHMARK_DIR = ROOT / "benchmarks"
RAW_DIR.mkdir(parents=True, exist_ok=True)
BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)

ZIP_PATH = RAW_DIR / "pathpredictor_s1.zip"
CSV_PATH = RAW_DIR / "pathpredictor_s1.csv"
OUTPUT_CSV = BENCHMARK_DIR / "benchmark_pathpredictor.csv"

TRAIN_GENES_PATH = ROOT / "train_genes_for_attachment.txt"

# ── Download URLs ────────────────────────────────────────────────────────────
# NCBI PMC often blocks curl; we try a few mirrors / endpoints.
URLS = [
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6633260/bin/supp_gr.240994.118_Supplemental_Table_S1_missenseDiseaseVariants_hg19.csv.zip",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC6633260/bin/supp_gr.240994.118_Supplemental_Table_S1_missenseDiseaseVariants_hg19.csv.zip",
    "https://genome.cshlp.org/content/suppl/2019/06/24/gr.240994.118.DC1/supp_gr.240994.118_Supplemental_Table_S1_missenseDiseaseVariants_hg19.csv.zip",
]

MANUAL_INSTRUCTIONS = """
╔══════════════════════════════════════════════════════════════════════════════╗
║  AUTOMATIC DOWNLOAD FAILED                                                   ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Please download the PathPredictor S1 Table manually:                        ║
║                                                                              ║
║  1. Go to: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6633260/             ║
║  2. Scroll to "Supplementary Material"                                       ║
║  3. Download:                                                                ║
║     "Supplemental Table S1 missenseDiseaseVariants hg19 (CSV)"               ║
║                                                                              ║
║  4. Place the downloaded ZIP file at:                                        ║
║     {zip_path}                                                               ║
║                                                                              ║
║  5. Re-run this script.                                                      ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

# ── Helpers ──────────────────────────────────────────────────────────────────

def annotate_genes_from_refgene(df: pd.DataFrame) -> pd.DataFrame:
    """Map variant positions to gene symbols using refGene exon coordinates via bedtools."""
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
    # refGene columns:
    # 0:bin, 1:name, 2:chrom, 3:strand, 4:txStart, 5:txEnd, 6:cdsStart, 7:cdsEnd,
    # 8:exonCount, 9:exonStarts, 10:exonEnds, 11:score, 12:name2(gene)
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

        # Parse output: var_chrom, var_start, var_end, idx, exon_chrom, exon_start, exon_end, gene
        idx_to_gene = {}
        with open(out_path) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 8:
                    idx = int(cols[3])
                    gene = cols[7]
                    # Keep first gene hit (or we could collect all and pick best)
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


def download() -> bool:
    """Try to download automatically; return True on success."""
    if ZIP_PATH.exists():
        print(f"✅ Found existing zip: {ZIP_PATH}")
        return True

    for url in URLS:
        print(f"⬇️  Trying: {url}")
        try:
            urlretrieve(url, ZIP_PATH)
            # Sanity check: must be a zip (first 2 bytes = PK)
            with open(ZIP_PATH, "rb") as f:
                header = f.read(2)
            if header == b"PK":
                print(f"✅ Downloaded successfully")
                return True
            else:
                print(f"   ⚠️  Returned HTML instead of zip (likely blocked)")
                ZIP_PATH.unlink()
        except HTTPError as e:
            print(f"   ❌ HTTP {e.code}")
        except Exception as e:
            print(f"   ❌ {e}")
    return False


def extract_zip():
    """Extract the CSV from the zip."""
    with zipfile.ZipFile(ZIP_PATH, "r") as zf:
        # Find the CSV inside
        csv_members = [m for m in zf.namelist() if m.lower().endswith(".csv")]
        if not csv_members:
            raise ValueError("No CSV found inside the zip!")
        member = csv_members[0]
        print(f"📦 Extracting {member} ...")
        zf.extract(member, RAW_DIR)
        extracted = RAW_DIR / member
        # Rename to our standard name
        extracted.rename(CSV_PATH)
    print(f"✅ Extracted to {CSV_PATH}")


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
    # 1. Download / extract (or use existing CSV)
    if CSV_PATH.exists():
        print(f"✅ Using existing CSV: {CSV_PATH}")
    else:
        if not download():
            print(MANUAL_INSTRUCTIONS.format(zip_path=ZIP_PATH))
            sys.exit(1)
        extract_zip()

    # 3. Load S1 Table
    print(f"📖 Reading S1 Table...")
    df = pd.read_csv(CSV_PATH, dtype=str, comment="#")
    print(f"   Total rows: {len(df):,}")
    print(f"   Columns: {list(df.columns)}")

    # 4. Standardize column names
    # Expected columns: chrom, pos, ref, alt, class, Disease
    # where class = B/P/V
    rename = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ["chrom", "chr", "chromosome"]:
            rename[c] = "chrom"
        elif cl == "pos":
            rename[c] = "pos"
        elif cl == "ref":
            rename[c] = "ref"
        elif cl == "alt":
            rename[c] = "alt"
        elif cl in ["gene", "genesymbol", "gene_symbol"]:
            rename[c] = "gene"
        elif cl in ["class", "pathogenicity", "label", "category"]:
            rename[c] = "pathogenicity"
        elif cl in ["disease", "panel"]:
            rename[c] = "disease"
    df = df.rename(columns=rename)

    # If gene column is missing, annotate from refGene
    if "gene" not in df.columns or df["gene"].eq("UNKNOWN").all():
        df = annotate_genes_from_refgene(df)

    # 5. Filter to B / P only (drop VUS)
    df = df[df["pathogenicity"].isin(["B", "P"])].copy()
    df["label"] = df["pathogenicity"].map({"B": 0, "P": 1})
    print(f"   After B/P filter: {len(df):,}")

    # 6. Clean coordinates
    df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", case=False)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["ref"] = df["ref"].str.upper()
    df["alt"] = df["alt"].str.upper()
    df = df.dropna(subset=["chrom", "pos", "ref", "alt", "gene"]).copy()

    # 7. Liftover hg19 → GRCh38
    df = liftover_hg19_to_hg38(df)

    # 8. Load training holdout sets
    train_vars, train_genes = load_training_genes_and_variants()

    # 9. Gene holdout
    before = len(df)
    df = df[~df["gene"].isin(train_genes)].copy()
    print(f"   After gene holdout:   {len(df):,} (removed {before - len(df):,})")

    # 10. Variant holdout
    before = len(df)
    df["_key"] = list(zip(df["chrom"], df["pos"].astype(str), df["ref"], df["alt"]))
    df = df[~df["_key"].isin(train_vars)].drop(columns=["_key"]).copy()
    print(f"   After variant holdout: {len(df):,} (removed {before - len(df):,})")

    # 11. Final formatting
    df = df.reset_index(drop=True)
    out_df = pd.DataFrame({
        "chrom": df["chrom"],
        "pos": df["pos"].astype(int),
        "ref": df["ref"],
        "alt": df["alt"],
        "gene": df["gene"],
        "label": df["label"].astype(int),
        "disease": df.get("disease", ""),
        "dataset": "pathpredictor",
    })

    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n✅ Benchmark saved: {OUTPUT_CSV}")
    print(f"   Total variants: {len(out_df):,}")
    print(f"   Pathogenic (1): {(out_df['label'] == 1).sum():,}")
    print(f"   Benign     (0): {(out_df['label'] == 0).sum():,}")
    print(f"   Genes: {out_df['gene'].nunique():,}")


if __name__ == "__main__":
    main()
