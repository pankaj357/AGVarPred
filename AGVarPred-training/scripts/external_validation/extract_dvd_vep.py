#!/usr/bin/env python3
"""
Extract VEP annotations for DVD benchmark via Ensembl REST API.
Saves to external_validation/vep_preprocessed/dvd_vep.parquet
"""
import json
import time
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import requests

ROOT = Path("external_validation")
BENCH_PATH = ROOT / "benchmarks/benchmark_dvd.csv"
VEP_OUT = ROOT / "vep_preprocessed/dvd_vep.parquet"

BATCH_SIZE = 200
SLEEP_BETWEEN_BATCHES = 0.5
API_URL = "https://rest.ensembl.org/vep/homo_sapiens/region?canonical=1"


def variant_to_vcf_string(chrom, pos, ref, alt):
    """Convert to VCF-style string for Ensembl API."""
    chrom = str(chrom).replace("chr", "")
    return f"{chrom} {pos} . {ref} {alt} . . ."


def parse_api_response(data):
    """Parse Ensembl VEP API response into our feature format."""
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
    for vid, chrom, pos, ref, alt in batch_variants:
        vcf = variant_to_vcf_string(chrom, pos, ref, alt)
        if vcf:
            vcf_strings.append(vcf)
    
    if not vcf_strings:
        return {}
    
    payload = {"variants": vcf_strings}
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    try:
        resp = requests.post(API_URL, headers=headers, json=payload, timeout=60)
        if resp.status_code == 200:
            return parse_api_response(resp.json())
        elif resp.status_code == 429:
            print(f"   Rate limited, sleeping 5s...")
            time.sleep(5)
            return fetch_batch(batch_variants)
        else:
            print(f"   API error {resp.status_code}: {resp.text[:200]}")
            return {}
    except Exception as e:
        print(f"   Request error: {e}")
        return {}


def main():
    print("="*70)
    print("DVD VEP EXTRACTION — Ensembl REST API")
    print("="*70)
    
    df = pd.read_csv(BENCH_PATH, dtype=str)
    print(f"Benchmark variants: {len(df):,}")
    
    df["variant_id"] = (
        "chr" + df["chrom"].astype(str).str.replace("chr", "", case=False) + "_" +
        df["pos"].astype(str).str.replace(r"\.0$", "", regex=True) + "_" +
        df["ref"].astype(str).str.upper() + "_" +
        df["alt"].astype(str).str.upper()
    )
    
    variants = list(zip(df["variant_id"], df["chrom"], df["pos"], df["ref"], df["alt"]))
    
    all_results = {}
    n_batches = (len(variants) + BATCH_SIZE - 1) // BATCH_SIZE
    
    for i in tqdm(range(n_batches), desc="Fetching VEP"):
        batch = variants[i*BATCH_SIZE:(i+1)*BATCH_SIZE]
        results = fetch_batch(batch)
        all_results.update(results)
        time.sleep(SLEEP_BETWEEN_BATCHES)
    
    print(f"\nVEP annotations retrieved for {len(all_results):,} / {len(variants):,} variants")
    
    rows = []
    for vid in df["variant_id"]:
        row = {"variant_id": vid}
        row.update(all_results.get(vid, {}))
        rows.append(row)
    
    result_df = pd.DataFrame(rows).set_index("variant_id")
    
    # Preprocess to match training format
    numeric_cols = ["vep_SIFT_score", "vep_PolyPhen_score", "vep_Protein_position"]
    for col in numeric_cols:
        if col in result_df.columns:
            result_df[col] = pd.to_numeric(result_df[col], errors="coerce").astype("float32")
    
    if "vep_IMPACT_score" in result_df.columns:
        result_df["vep_IMPACT_score"] = result_df["vep_IMPACT_score"].fillna(0).astype("int8")
    
    binary_cols = ["vep_is_missense", "vep_is_synonymous", "vep_is_stop_gained",
                   "vep_is_frameshift", "vep_is_splice", "vep_is_LoF_HC"]
    for col in binary_cols:
        if col not in result_df.columns:
            result_df[col] = 0
        result_df[col] = result_df[col].fillna(0).astype("int8")
    
    result_df["vep_has_SIFT"] = result_df["vep_SIFT_score"].notna().astype("int8")
    result_df["vep_has_PolyPhen"] = result_df["vep_PolyPhen_score"].notna().astype("int8")
    
    if "vep_LoF" in result_df.columns:
        lof_map = {"LC": 0, "HC": 1}
        result_df["vep_LoF"] = result_df["vep_LoF"].map(lof_map).fillna(-1).astype("int16")
    
    if "vep_IMPACT" in result_df.columns:
        impact_map = {"MODIFIER": 0, "LOW": 1, "MODERATE": 2, "HIGH": 3}
        result_df["vep_IMPACT"] = result_df["vep_IMPACT"].map(impact_map).fillna(-1).astype("int16")
    
    if "vep_SIFT_pred" in result_df.columns:
        sift_map = {"tolerated": 0, "tolerated_low_confidence": 1, "deleterious_low_confidence": 2, "deleterious": 3}
        result_df["vep_SIFT_pred"] = result_df["vep_SIFT_pred"].map(sift_map).fillna(-1).astype("int16")
    
    if "vep_PolyPhen_pred" in result_df.columns:
        pp_map = {"benign": 0, "possibly_damaging": 1, "probably_damaging": 2, "unknown": 3}
        result_df["vep_PolyPhen_pred"] = result_df["vep_PolyPhen_pred"].map(pp_map).fillna(-1).astype("int16")
    
    if "vep_Consequence" in result_df.columns:
        for cat_val in result_df["vep_Consequence"].dropna().unique():
            col_name = f"vep_Consequence_{cat_val}"
            result_df[col_name] = (result_df["vep_Consequence"] == cat_val).astype("int8")
    
    result_df.to_parquet(VEP_OUT)
    print(f"\n✅ Saved: {VEP_OUT}")
    print(f"   Shape: {result_df.shape}")


if __name__ == "__main__":
    main()
