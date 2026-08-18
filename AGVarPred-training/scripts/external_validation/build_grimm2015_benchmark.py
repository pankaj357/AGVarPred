#!/usr/bin/env python3
"""
Build an independent benchmark from Grimm et al. 2015 Selected datasets.

Source:
    Grimm, D.G. et al. "The evaluation of tools used to predict the impact of
    missense variants is hindered by two types of circularity." Human Mutation
    36.5 (2015): 513-523.

Datasets used:
    - varibench_selected_tool_scores.csv
    - predictSNP_selected_tool_scores.csv
    - swissvar_selected_tool_scores.csv

These datasets were explicitly designed to avoid circularity with predictor
training sets. We additionally filter them to ensure:
    1. Zero gene overlap with the AlphaGenome training set
    2. Zero variant overlap with train/cal/test
    3. GRCh38 coordinates (lifted from GRCh37)

Output:
    external_validation/benchmarks/benchmark_grimm2015.csv
"""

import os
import sys
import json
import time
import pandas as pd
import pyarrow.parquet as pq
from pyliftover import LiftOver
import requests
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SOURCE_DIR = PROJECT_ROOT / "external_validation" / "source_datasets" / "grimm2015"
OUT_DIR = PROJECT_ROOT / "external_validation" / "benchmarks"
OUT_BENCHMARK = OUT_DIR / "benchmark_grimm2015.csv"

TRAIN_GENES_FILE = PROJECT_ROOT / "external_validation" / "train_genes_upper.json"
TRAIN_VARIANTS_CACHE = PROJECT_ROOT / "external_validation" / "processing" / "grimm2015_train_variants.pkl"
GENE_MAP_CACHE = SOURCE_DIR / "ensembl_to_hgnc.json"

SOURCE_FILES = {
    "varibench_selected": "varibench_selected_tool_scores.csv",
    "predictSNP_selected": "predictSNP_selected_tool_scores.csv",
    "swissvar_selected": "swissvar_selected_tool_scores.csv",
}


def collect_training_variants():
    """Collect all variant IDs from train/cal/test parquet parts using pyarrow.dataset."""
    if TRAIN_VARIANTS_CACHE.exists():
        print(f"Loading cached training variants from {TRAIN_VARIANTS_CACHE}")
        return pd.read_pickle(TRAIN_VARIANTS_CACHE)

    import pyarrow.dataset as ds

    variants = set()
    for split in ["train", "cal", "test"]:
        parts_dir = PROJECT_ROOT / f"final_dataset_parts_{split}"
        if not parts_dir.exists():
            print(f"Warning: {parts_dir} not found, skipping")
            continue
        print(f"Scanning {split} parquet dataset...")
        dataset = ds.dataset(parts_dir, format="parquet")
        table = dataset.to_table(columns=["variant_id"])
        split_variants = set(table.column("variant_id").to_pylist())
        variants.update(split_variants)
        print(f"  {split}: {len(split_variants)} unique variants")

    TRAIN_VARIANTS_CACHE.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(variants, TRAIN_VARIANTS_CACHE)
    print(f"Cached {len(variants)} total training variants")
    return variants


def fetch_ensembl_to_hgnc(ensembl_ids):
    """Map Ensembl Gene IDs to HGNC symbols via Ensembl REST API."""
    if GENE_MAP_CACHE.exists():
        print(f"Loading cached Ensembl->HGNC mapping from {GENE_MAP_CACHE}")
        with open(GENE_MAP_CACHE) as f:
            return json.load(f)

    mapping = {}
    ids = sorted(set(eid for eid in ensembl_ids if pd.notna(eid) and eid != ""))
    print(f"Querying Ensembl REST API for {len(ids)} unique gene IDs...")

    batch_size = 1000
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    for start in range(0, len(ids), batch_size):
        batch = ids[start:start + batch_size]
        payload = json.dumps({"ids": batch})
        for attempt in range(3):
            try:
                r = requests.post(url, headers=headers, data=payload, timeout=60)
                if r.status_code == 200:
                    data = r.json()
                    for eid in batch:
                        info = data.get(eid)
                        if info and "display_name" in info:
                            mapping[eid] = info["display_name"]
                        else:
                            mapping[eid] = None
                    break
                else:
                    print(f"  Batch {start}-{start+len(batch)} HTTP {r.status_code}: {r.text[:200]}")
                    time.sleep(2 ** attempt)
            except Exception as e:
                print(f"  Batch {start}-{start+len(batch)} error: {e}")
                time.sleep(2 ** attempt)
        else:
            for eid in batch:
                mapping[eid] = None
        print(f"  ... mapped {min(start + batch_size, len(ids))}/{len(ids)}")
        time.sleep(0.1)

    with open(GENE_MAP_CACHE, "w") as f:
        json.dump(mapping, f, indent=2)
    print(f"Cached Ensembl->HGNC mapping ({sum(1 for v in mapping.values() if v)} mapped)")
    return mapping


def liftover_row(row, lo):
    """Lift over a single variant from hg19 to hg38."""
    chrom = str(row["CHR"])
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    pos = int(row["Nuc-Pos"])
    result = lo.convert_coordinate(chrom, pos - 1)  # 0-based
    if result and len(result) > 0:
        new_chrom, new_pos_0based, strand, _ = result[0]
        return new_chrom, new_pos_0based + 1
    return None, None


def build_benchmark():
    print("=" * 60)
    print("Building Grimm2015 independent benchmark")
    print("=" * 60)

    # Load training genes
    with open(TRAIN_GENES_FILE) as f:
        train_genes = set(json.load(f))
    train_genes.discard("-")
    print(f"Training genes: {len(train_genes)}")

    # Load training variants
    train_variants = collect_training_variants()
    print(f"Training variants: {len(train_variants)}")

    # Load source datasets
    all_dfs = []
    for name, fname in SOURCE_FILES.items():
        path = SOURCE_DIR / fname
        print(f"\nLoading {name} from {path}")
        df = pd.read_csv(path)
        df["source"] = name
        print(f"  Rows: {len(df)}, Pathogenic: {(df['True Label'] == 1).sum()}, Neutral: {(df['True Label'] == -1).sum()}")
        all_dfs.append(df)

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\nCombined rows before filtering: {len(combined)}")

    # Map Ensembl IDs to HGNC symbols
    gene_map = fetch_ensembl_to_hgnc(combined["Ensembl-Gene-ID"].unique())
    combined["gene"] = combined["Ensembl-Gene-ID"].map(gene_map)

    # Drop rows without a gene symbol
    before = len(combined)
    combined = combined[combined["gene"].notna()].copy()
    print(f"\nDropped {before - len(combined)} rows with unmapped Ensembl gene ID")

    # Filter by gene overlap with training
    before = len(combined)
    combined = combined[~combined["gene"].str.upper().isin({g.upper() for g in train_genes})].copy()
    print(f"Dropped {before - len(combined)} rows in training genes")
    print(f"Remaining rows after gene filter: {len(combined)}")

    # Lift over coordinates
    print("\nLifting over coordinates from GRCh37 to GRCh38...")
    lo = LiftOver("hg19", "hg38")
    lifted = combined.apply(lambda row: liftover_row(row, lo), axis=1)
    combined["chrom_grch38"] = lifted.apply(lambda x: x[0])
    combined["pos_grch38"] = lifted.apply(lambda x: x[1])

    before = len(combined)
    combined = combined[combined["chrom_grch38"].notna()].copy()
    print(f"Dropped {before - len(combined)} rows that failed liftover")
    print(f"Remaining rows after liftover: {len(combined)}")

    # Build variant_id and filter by variant overlap
    combined["variant_id"] = (
        combined["chrom_grch38"].astype(str) + "_" +
        combined["pos_grch38"].astype(int).astype(str) + "_" +
        combined["REF-Nuc"].astype(str) + "_" +
        combined["ALT-Nuc"].astype(str)
    )

    before = len(combined)
    combined = combined[~combined["variant_id"].isin(train_variants)].copy()
    print(f"Dropped {before - len(combined)} rows overlapping train/cal/test variants")
    print(f"Remaining rows after variant filter: {len(combined)}")

    # Convert labels: -1 -> 0 (benign), 1 -> 1 (pathogenic)
    combined["label"] = combined["True Label"].map({1: 1, -1: 0})

    # Build protein change
    combined["protein_change"] = combined.apply(
        lambda row: f"p.{row['REF-AA']}{int(row['AA-Pos'])}{row['ALT-AA']}", axis=1
    )

    # Build output dataframe in benchmark format
    out = pd.DataFrame({
        "chrom": combined["chrom_grch38"].str.replace("chr", ""),
        "pos": combined["pos_grch38"].astype(int),
        "ref": combined["REF-Nuc"],
        "alt": combined["ALT-Nuc"],
        "gene": combined["gene"],
        "rsid": combined["#RS-ID"].replace("?", "").fillna(""),
        "uniprot_ac": combined["Uniprot-Accession"],
        "protein_change": combined["protein_change"],
        "label": combined["label"],
        "category": combined["label"].map({1: "P/LP", 0: "B/LB"}),
        "disease": "-",
        "source": combined["source"],
        "ensembl_gene_id": combined["Ensembl-Gene-ID"],
    })

    # Final deduplication within the benchmark itself
    before = len(out)
    out = out.drop_duplicates(subset=["chrom", "pos", "ref", "alt"]).copy()
    print(f"Dropped {before - len(out)} duplicate variants within Grimm2015")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_BENCHMARK, index=False)
    print(f"\nSaved benchmark to {OUT_BENCHMARK}")
    print(f"Final benchmark size: {len(out)}")
    print(f"  Pathogenic (P/LP): {(out['label'] == 1).sum()}")
    print(f"  Benign (B/LB): {(out['label'] == 0).sum()}")
    print(f"  Genes: {out['gene'].nunique()}")
    print(f"  Sources:")
    print(out["source"].value_counts().to_string())

    # Save summary
    summary_path = OUT_BENCHMARK.with_suffix(".summary.txt")
    with open(summary_path, "w") as f:
        f.write("Grimm2015 Independent Benchmark Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total variants: {len(out)}\n")
        f.write(f"Pathogenic (P/LP): {(out['label'] == 1).sum()}\n")
        f.write(f"Benign (B/LB): {(out['label'] == 0).sum()}\n")
        f.write(f"Unique genes: {out['gene'].nunique()}\n")
        f.write("\nPer source:\n")
        f.write(out["source"].value_counts().to_string())
    print(f"Saved summary to {summary_path}")


if __name__ == "__main__":
    build_benchmark()
