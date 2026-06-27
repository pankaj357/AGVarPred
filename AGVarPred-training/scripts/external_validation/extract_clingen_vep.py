#!/usr/bin/env python3
"""
Fetch VEP annotations for the ClinGen benchmark via Ensembl REST API.

Saves:
  - external_validation/vep_preprocessed/clingen_vep.parquet

This script fetches VEP for all ClinGen variants in a single pass and applies
preprocessing (dtype fixes, consequence one-hot encoding) before saving.
"""

import json
import time
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import requests

PROJECT_ROOT = Path(__file__).resolve().parents[3]
EXT_DIR = PROJECT_ROOT / "external_validation"
BENCH_PATH = EXT_DIR / "benchmarks" / "benchmark_clingen.csv"
OUT_PATH = EXT_DIR / "vep_preprocessed" / "clingen_vep.parquet"

BATCH_SIZE = 200
SLEEP_BETWEEN_BATCHES = 0.5
API_URL = "https://rest.ensembl.org/vep/homo_sapiens/region?canonical=1"


def variant_to_vcf_string(variant_id):
    """Convert chr1_12345_A_G to VCF-style string for Ensembl API."""
    parts = variant_id.split("_")
    if len(parts) != 4:
        return None
    chrom, pos, ref, alt = parts
    chrom = chrom.replace("chr", "")
    return f"{chrom} {pos} . {ref} {alt} . . ."


def parse_api_response(data):
    """Parse Ensembl VEP API response into feature format."""
    results = {}
    for entry in data:
        inp = entry.get("input", "")
        parts = inp.split()
        if len(parts) < 5:
            continue
        chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        variant_id = f"chr{chrom}_{pos}_{ref}_{alt}"

        most_severe = entry.get("most_severe_consequence", "")

        canonical_tc = None
        for tc in entry.get("transcript_consequences", []):
            if tc.get("canonical") == 1:
                canonical_tc = tc
                break

        if canonical_tc is None:
            for tc in entry.get("transcript_consequences", []):
                if tc.get("sift_score") is not None or tc.get("polyphen_score") is not None:
                    canonical_tc = tc
                    break

        if canonical_tc is None:
            tcs = entry.get("transcript_consequences", [])
            if tcs:
                canonical_tc = tcs[0]

        if canonical_tc is None:
            continue

        impact = canonical_tc.get("impact", "")
        consequence_terms = canonical_tc.get("consequence_terms", [])
        if not consequence_terms and most_severe:
            consequence_terms = [most_severe]

        impact_map = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
        impact_score = impact_map.get(impact, 0) if impact else 0

        sift_score = canonical_tc.get("sift_score")
        sift_pred = canonical_tc.get("sift_prediction")
        polyphen_score = canonical_tc.get("polyphen_score")
        polyphen_pred = canonical_tc.get("polyphen_prediction")
        protein_pos = canonical_tc.get("protein_start")

        cons_list = consequence_terms
        result = {
            "vep_SIFT_score": sift_score,
            "vep_SIFT_pred": sift_pred,
            "vep_PolyPhen_score": polyphen_score,
            "vep_PolyPhen_pred": polyphen_pred,
            "vep_IMPACT": impact,
            "vep_IMPACT_score": impact_score,
            "vep_Consequence": "&".join(cons_list) if cons_list else most_severe,
            "vep_is_missense": 1 if "missense_variant" in cons_list else 0,
            "vep_is_synonymous": 1 if "synonymous_variant" in cons_list else 0,
            "vep_is_stop_gained": 1 if "stop_gained" in cons_list else 0,
            "vep_is_frameshift": 1 if "frameshift_variant" in cons_list else 0,
            "vep_is_splice": 1 if any(c in cons_list for c in ["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"]) else 0,
            "vep_LoF": None,
            "vep_is_LoF_HC": 0,
            "vep_Protein_position": protein_pos,
        }
        results[variant_id] = result
    return results


def fetch_batch(batch_variants):
    """Fetch VEP for a batch of variants."""
    vcf_strings = []
    for vid in batch_variants:
        vcf = variant_to_vcf_string(vid)
        if vcf:
            vcf_strings.append(vcf)

    if not vcf_strings:
        return {}

    payload = {"variants": vcf_strings}
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    try:
        resp = requests.post(API_URL, headers=headers, json=payload, timeout=120)
        if resp.status_code == 200:
            return parse_api_response(resp.json())
        elif resp.status_code == 429:
            print("   Rate limited, sleeping 5s...")
            time.sleep(5)
            return fetch_batch(batch_variants)
        else:
            print(f"   API error {resp.status_code}: {resp.text[:200]}")
            return {}
    except Exception as e:
        print(f"   Request error: {e}")
        return {}


def preprocess_vep(df):
    """Apply preprocessing to match training pipeline expectations."""
    df = df.copy()

    # dtype fixes
    for col in df.columns:
        if col in ["vep_SIFT_score", "vep_PolyPhen_score", "vep_IMPACT_score", "vep_Protein_position"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        elif col.startswith("vep_is_"):
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype("int8")

    # Ordinal mappings
    impact_map = {"MODIFIER": 1, "LOW": 2, "MODERATE": 3, "HIGH": 4}
    if "vep_IMPACT" in df.columns:
        df["vep_IMPACT_score"] = df["vep_IMPACT"].map(impact_map).fillna(0).astype(int)

    sift_map = {"deleterious": 2, "deleterious_low_confidence": 1, "tolerated_low_confidence": 0, "tolerated": -1}
    if "vep_SIFT_pred" in df.columns:
        df["vep_SIFT_pred_code"] = df["vep_SIFT_pred"].map(sift_map).fillna(-2).astype(int)

    polyphen_map = {"probably_damaging": 2, "possibly_damaging": 1, "benign": 0, "unknown": -1}
    if "vep_PolyPhen_pred" in df.columns:
        df["vep_PolyPhen_pred_code"] = df["vep_PolyPhen_pred"].map(polyphen_map).fillna(-2).astype(int)

    lof_map = {"HC": 2, "LC": 1}
    if "vep_LoF" in df.columns:
        df["vep_LoF_code"] = df["vep_LoF"].map(lof_map).fillna(0).astype(int)

    # One-hot encode Consequence
    if "vep_Consequence" in df.columns:
        for cat_val in df["vep_Consequence"].dropna().unique():
            df[f"vep_Consequence_{cat_val}"] = (df["vep_Consequence"] == cat_val).astype("int8")

    return df


def main():
    print("=" * 70)
    print("Fetching VEP annotations for ClinGen benchmark")
    print("=" * 70)

    bench = pd.read_csv(BENCH_PATH)
    variant_ids = bench["variant_id"].tolist()
    print(f"Total ClinGen variants: {len(variant_ids):,}")

    n_batches = (len(variant_ids) + BATCH_SIZE - 1) // BATCH_SIZE
    all_results = {}

    for i in tqdm(range(n_batches), desc="Fetching VEP batches"):
        start = i * BATCH_SIZE
        end = min((i + 1) * BATCH_SIZE, len(variant_ids))
        batch = variant_ids[start:end]
        results = fetch_batch(batch)
        all_results.update(results)
        if i < n_batches - 1:
            time.sleep(SLEEP_BETWEEN_BATCHES)

    print(f"\nAnnotations retrieved: {len(all_results):,} / {len(variant_ids):,}")

    # Build dataframe
    rows = []
    for vid in variant_ids:
        row = {"variant_id": vid}
        row.update(all_results.get(vid, {}))
        rows.append(row)

    df = pd.DataFrame(rows).set_index("variant_id")
    df = preprocess_vep(df)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT_PATH)
    print(f"Saved preprocessed VEP to {OUT_PATH}")

    # Coverage report
    print("\nVEP coverage:")
    for col in ["vep_SIFT_score", "vep_PolyPhen_score", "vep_IMPACT_score"]:
        if col in df.columns:
            cov = df[col].notna().mean() * 100
            print(f"  {col}: {cov:.1f}%")


if __name__ == "__main__":
    main()
