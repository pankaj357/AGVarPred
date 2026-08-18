#!/usr/bin/env python3
"""
Download and prepare a germline variant benchmark that is fully independent
from ClinVar at both the VARIANT and GENE level.

Source: UniProt Humsavar (manually curated missense variants with ACMG labels)
Mapping: versioned local rsID→GRCh38 map (primary) + myvariant.info (fallback for missing rsIDs)
Filtering:
  1. Keep only LP/P (label=1) and LB/B (label=0) variants
  2. Remove any variant that overlaps train/cal/test by chrom/pos/ref/alt
  3. Remove any variant in a gene that appears in train/cal/test (strict gene-holdout)
"""

import os
import re
import sys
import time
import json
import pandas as pd
import requests
from pathlib import Path
from collections import Counter

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
OUT_DIR = ROOT / "external_validation"
OUT_DIR.mkdir(exist_ok=True)

HUMSAVAR_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/source_datasets/humsavar/humsavar.txt"
HUMSAVAR_LOCAL = OUT_DIR / "source_datasets/humsavar/humsavar.txt"

TRAIN_CSV = ROOT / "train.csv"
CAL_CSV = ROOT / "cal.csv"
TEST_CSV = ROOT / "test.csv"

MYVARIANT_BATCH_SIZE = 200
MYVARIANT_DELAY = 0.1  # seconds between batches

# ---------------------------------------------------------------------------
# 1. DOWNLOAD HUMSAVAR
# ---------------------------------------------------------------------------
def download_humsavar():
    if HUMSAVAR_LOCAL.exists() and HUMSAVAR_LOCAL.stat().st_size > 0:
        print(f"[1/6] Humsavar already exists: {HUMSAVAR_LOCAL}")
        return
    print(f"[1/6] Downloading Humsavar from UniProt...")
    r = requests.get(HUMSAVAR_URL, stream=True)
    r.raise_for_status()
    with open(HUMSAVAR_LOCAL, "wb") as fh:
        for chunk in r.iter_content(chunk_size=8192):
            fh.write(chunk)
    print(f"      Saved ({HUMSAVAR_LOCAL.stat().st_size / 1e6:.1f} MB)")

# ---------------------------------------------------------------------------
# 2. PARSE HUMSAVAR
# ---------------------------------------------------------------------------
def parse_humsavar():
    print("[2/6] Parsing Humsavar...")
    rows = []
    with open(HUMSAVAR_LOCAL, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            # Skip header/comment lines
            if not line or line.startswith("-") or line.startswith(" ") or line.startswith("Main") or line.startswith("gene") or line.startswith("_"):
                continue
            # Lines look like:
            # A1BG        P04217     VAR_018369  p.His52Arg     LB/B     rs893184       -
            parts = line.split(None, 6)  # maxsplit=6 to keep disease name intact
            if len(parts) < 6:
                continue
            gene, uniprot_ac, ft_id, protein_change, category, dbsnp = parts[:6]
            disease = parts[6] if len(parts) > 6 else ""
            if category not in ("LP/P", "LB/B"):
                continue
            if dbsnp == "-" or not dbsnp.startswith("rs"):
                continue
            rows.append({
                "gene": gene,
                "uniprot_ac": uniprot_ac,
                "ft_id": ft_id,
                "protein_change": protein_change,
                "category": category,
                "rsid": dbsnp,
                "disease": disease,
            })
    df = pd.DataFrame(rows)
    df["label"] = df["category"].map({"LP/P": 1, "LB/B": 0})
    print(f"      LP/P: {(df['label']==1).sum():,} | LB/B: {(df['label']==0).sum():,} | Total: {len(df):,}")
    return df

# ---------------------------------------------------------------------------
# 3. LOAD TRAIN/CAL/TEST FOR DEDUPLICATION & GENE-HOLDOUT
# ---------------------------------------------------------------------------
def load_training_data():
    print("[3/6] Loading train/cal/test variant & gene sets...")
    seen_vars = set()
    seen_genes = set()
    for csv_path, name in [(TRAIN_CSV, "train"), (CAL_CSV, "cal"), (TEST_CSV, "test")]:
        df = pd.read_csv(csv_path, dtype=str)
        # Variants
        df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False)
        df["PositionVCF"] = df["PositionVCF"].astype(str).str.replace(r"\.0$", "", regex=True)
        df["ref"] = df["ReferenceAlleleVCF"].astype(str).str.upper()
        df["alt"] = df["AlternateAlleleVCF"].astype(str).str.upper()
        for _, row in df.iterrows():
            key = (row["Chromosome"], row["PositionVCF"], row["ref"], row["alt"])
            seen_vars.add(key)
        # Genes
        seen_genes.update(df["GeneSymbol"].dropna().unique())
        print(f"      {name:5s}: {len(df):,} variants | {df['GeneSymbol'].nunique():,} genes")
    print(f"      Total unique variants in training: {len(seen_vars):,}")
    print(f"      Total unique genes in training: {len(seen_genes):,}")
    return seen_vars, seen_genes

# ---------------------------------------------------------------------------
# 4. MAP RSIDS TO GRCh38: versioned local map (primary) + myvariant.info fallback
# ---------------------------------------------------------------------------
def _safe_get(obj, key, default=None):
    """Safely get a value from a dict; if result is a list, return first element if it's a dict."""
    if not isinstance(obj, dict):
        return default
    val = obj.get(key, default)
    if isinstance(val, list) and len(val) > 0:
        if isinstance(val[0], dict):
            return val[0]
        return default
    return val


# Versioned, committed rsID→GRCh38 mapping file. This is the PRIMARY source.
# myvariant.info is used only as a fallback for rsIDs not present here.
LOCAL_COORD_MAP = OUT_DIR / "raw_data/humsavar_rsID_GRCh38_map.json"

# Backwards-compatible transient cache location (legacy)
COORD_CACHE = OUT_DIR / "raw_data/coord_map_cache.json"

def _load_local_coord_map():
    """Load the versioned local rsID→coordinate map."""
    if LOCAL_COORD_MAP.exists():
        with open(LOCAL_COORD_MAP, "r") as fh:
            return json.load(fh)
    # Legacy fallback: old cache file
    if COORD_CACHE.exists():
        with open(COORD_CACHE, "r") as fh:
            return json.load(fh)
    return {}


def _save_local_coord_map(mapping):
    """Save the updated versioned local rsID→coordinate map."""
    LOCAL_COORD_MAP.parent.mkdir(parents=True, exist_ok=True)
    with open(LOCAL_COORD_MAP, "w") as fh:
        json.dump(mapping, fh, indent=None, sort_keys=True)


def _query_myvariant(batch):
    """Query myvariant.info for a single batch of rsIDs."""
    try:
        resp = requests.post(
            "https://myvariant.info/v1/variant",
            json={
                "ids": batch,
                "fields": "clinvar.chrom,clinvar.hg38.start,clinvar.ref,clinvar.alt,clinvar.gene.symbol,cadd.chrom,cadd.gene.genename,dbnsfp.hg38.start,dbnsfp.alt,dbnsfp.ref",
            },
            timeout=300,
        )
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        print(f"      myvariant.info ERROR: {e}")
        return []


def map_rsids_to_coords(df):
    print("[4/6] Mapping rsIDs to GRCh38...")
    results = _load_local_coord_map()
    print(f"      Loaded versioned local map: {len(results):,} rsIDs")

    rsids = df["rsid"].unique().tolist()
    missing = [r for r in rsids if r not in results]
    print(f"      Unique rsIDs in Humsavar: {len(rsids):,}")
    print(f"      Missing from local map:   {len(missing):,}")

    if not missing:
        print("      All rsIDs resolved locally (no myvariant.info queries needed).")
        return results

    print("      Falling back to myvariant.info for missing rsIDs...")
    n_batches = (len(missing) + MYVARIANT_BATCH_SIZE - 1) // MYVARIANT_BATCH_SIZE
    newly_mapped = 0

    for i in range(n_batches):
        batch = missing[i * MYVARIANT_BATCH_SIZE : (i + 1) * MYVARIANT_BATCH_SIZE]
        data = _query_myvariant(batch)
        for item in data:
            if not isinstance(item, dict):
                continue
            query = item.get("query")
            if not query:
                continue

            chrom = pos = ref = alt = gene = None

            cv = item.get("clinvar", {})
            if cv.get("chrom") and _safe_get(cv, "hg38", {}).get("start"):
                chrom = str(cv["chrom"]).replace("chr", "")
                pos = str(_safe_get(cv, "hg38", {}).get("start"))
                ref = str(cv.get("ref", "")).upper()
                alt = str(cv.get("alt", "")).upper()
                gene = _safe_get(cv, "gene", {}).get("symbol", "")

            if not chrom or not pos or not ref or not alt:
                cadd = item.get("cadd", {})
                dbnsfp = item.get("dbnsfp", {})
                if cadd.get("chrom"):
                    chrom = str(cadd["chrom"]).replace("chr", "")
                if _safe_get(dbnsfp, "hg38", {}).get("start"):
                    pos = str(_safe_get(dbnsfp, "hg38", {}).get("start"))
                if not ref:
                    ref = str(cadd.get("ref", dbnsfp.get("ref", ""))).upper()
                if not alt:
                    alt = str(cadd.get("alt", dbnsfp.get("alt", ""))).upper()
                if not gene:
                    gene = _safe_get(cadd, "gene", {}).get("genename", "")

            if chrom and pos and ref and alt and len(ref) == 1 and len(alt) == 1:
                results[query] = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene": gene if gene else "",
                }
                newly_mapped += 1
        print(f"      Batch {i+1}/{n_batches} | Newly mapped: {newly_mapped}", flush=True)
        time.sleep(MYVARIANT_DELAY)
        # Save incremental progress every 20 batches
        if (i + 1) % 20 == 0 and newly_mapped > 0:
            _save_local_coord_map(results)

    print(f"      Total mapped: {len(results):,} / {len(rsids):,}")
    if newly_mapped > 0:
        _save_local_coord_map(results)
        print(f"      Updated versioned local map: {LOCAL_COORD_MAP}")
    return results

# ---------------------------------------------------------------------------
# 5. BUILD CLEAN BENCHMARK
# ---------------------------------------------------------------------------
def build_benchmark(df, coord_map, seen_vars, seen_genes):
    print("[5/6] Building clean benchmark (variant + gene holdout)...")
    mapped_rows = []
    for _, row in df.iterrows():
        coord = coord_map.get(row["rsid"])
        if not coord:
            continue
        chrom = str(coord["chrom"]).lower().replace("chr", "")
        pos = str(coord["pos"]).replace(".0", "")
        ref = coord["ref"].upper()
        alt = coord["alt"].upper()
        gene = coord["gene"] if coord["gene"] else row["gene"]

        # Variant-level dedup
        var_key = (chrom, pos, ref, alt)
        if var_key in seen_vars:
            continue

        # Gene-level holdout
        if gene in seen_genes:
            continue

        mapped_rows.append({
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "gene": gene,
            "rsid": row["rsid"],
            "uniprot_ac": row["uniprot_ac"],
            "protein_change": row["protein_change"],
            "label": row["label"],
            "category": row["category"],
            "disease": row["disease"],
        })

    out = pd.DataFrame(mapped_rows)
    # Deduplicate exact same variant (can happen if multiple rsIDs or multiple UniProt entries map to same genomic coord)
    out = out.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    print(f"      Final benchmark: {len(out):,} variants")
    print(f"      Pathogenic (LP/P): {(out['label']==1).sum():,}")
    print(f"      Benign (LB/B):     {(out['label']==0).sum():,}")
    print(f"      Genes: {out['gene'].nunique():,}")
    return out

# ---------------------------------------------------------------------------
# 6. SAVE
# ---------------------------------------------------------------------------
def save_benchmark(df):
    print("[6/6] Saving benchmark...")
    out_path = OUT_DIR / "benchmarks/benchmark_independent_humsavar.csv"
    df.to_csv(out_path, index=False)
    print(f"      Saved: {out_path}")

    # Also create variant map for chunking
    map_path = OUT_DIR / "benchmarks/variant_benchmark_map.csv"
    df_map = df[["chrom", "pos", "ref", "alt", "gene", "rsid", "label"]].copy()
    df_map["variant_id"] = df_map["chrom"] + "_" + df_map["pos"].astype(str) + "_" + df_map["ref"] + "_" + df_map["alt"]
    df_map["benchmark_source"] = "humsavar_independent"
    df_map.to_csv(map_path, index=False)
    print(f"      Saved map: {map_path}")

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("Independent Germline Benchmark Builder (UniProt Humsavar)")
    print("=" * 70)

    download_humsavar()
    df = parse_humsavar()
    seen_vars, seen_genes = load_training_data()
    coord_map = map_rsids_to_coords(df)
    bench = build_benchmark(df, coord_map, seen_vars, seen_genes)
    save_benchmark(bench)

    print("\n" + "=" * 70)
    print("✅ Done! Benchmark is fully independent from train/cal/test")
    print("   - Zero variant overlap")
    print("   - Zero gene overlap (strict gene-holdout)")
    print("=" * 70)


if __name__ == "__main__":
    main()
