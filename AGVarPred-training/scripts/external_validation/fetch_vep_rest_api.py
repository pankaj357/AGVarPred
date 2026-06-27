#!/usr/bin/env python3
"""
Fetch missing VEP annotations via Ensembl REST API for external validation benchmarks.

This fills gaps where variants are not present in gnomAD VCF (e.g., MAVE).
Uses the VEP region endpoint with POST batches of 200 variants.
"""

import json
import os
import time
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import requests

ROOT = Path(os.environ.get("PROJECT_ROOT", ".")) / "external_validation"
VEP_DIR = ROOT / "vep"
BATCH_SIZE = 200
SLEEP_BETWEEN_BATCHES = 0.5  # be polite to Ensembl servers
API_URL = "https://rest.ensembl.org/vep/homo_sapiens/region?canonical=1"

KEY_VEP_COLS = ["vep_IMPACT_score", "vep_SIFT_score", "vep_PolyPhen_score"]


def variant_to_vcf_string(variant_id):
    """Convert chr1_12345_A_G to VCF-style string for Ensembl API."""
    parts = variant_id.split("_")
    if len(parts) != 4:
        return None
    chrom, pos, ref, alt = parts
    chrom = chrom.replace("chr", "")
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

        # Use most severe consequence as overall
        most_severe = entry.get("most_severe_consequence", "")
        
        # Find canonical transcript consequences
        canonical_tc = None
        for tc in entry.get("transcript_consequences", []):
            if tc.get("canonical") == 1:
                canonical_tc = tc
                break
        
        # Fallback to first transcript with sift/polyphen if no canonical
        if canonical_tc is None:
            for tc in entry.get("transcript_consequences", []):
                if tc.get("sift_score") is not None or tc.get("polyphen_score") is not None:
                    canonical_tc = tc
                    break
        
        # Ultimate fallback: first transcript
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
            "vep_LoF": None,  # Not available via REST API
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
        resp = requests.post(API_URL, headers=headers, json=payload, timeout=60)
        if resp.status_code == 200:
            return parse_api_response(resp.json())
        elif resp.status_code == 429:
            print(f"   Rate limited, sleeping 5s...")
            time.sleep(5)
            return fetch_batch(batch_variants)  # retry
        else:
            print(f"   API error {resp.status_code}: {resp.text[:200]}")
            return {}
    except Exception as e:
        print(f"   Request error: {e}")
        return {}


def process_benchmark(vep_file):
    """Fetch missing VEP annotations for a single benchmark."""
    vep_path = VEP_DIR / vep_file
    print(f"\n{'='*60}")
    print(f"Processing: {vep_file}")
    print(f"{'='*60}")
    
    df = pd.read_parquet(vep_path)
    print(f"Total variants: {len(df):,}")
    
    # Find variants missing key VEP columns
    missing_mask = pd.Series(True, index=df.index)
    for col in KEY_VEP_COLS:
        if col in df.columns:
            missing_mask &= df[col].isna()
    
    missing_variants = df.index[missing_mask].tolist()
    print(f"Missing key VEP annotations: {len(missing_variants):,}")
    
    if len(missing_variants) == 0:
        print("No missing annotations to fetch.")
        return
    
    # Fetch in batches
    all_results = {}
    n_batches = (len(missing_variants) + BATCH_SIZE - 1) // BATCH_SIZE
    
    for i in tqdm(range(n_batches), desc="Fetching VEP batches"):
        start = i * BATCH_SIZE
        end = min((i + 1) * BATCH_SIZE, len(missing_variants))
        batch = missing_variants[start:end]
        
        results = fetch_batch(batch)
        all_results.update(results)
        
        if i < n_batches - 1:
            time.sleep(SLEEP_BETWEEN_BATCHES)
    
    print(f"Annotations retrieved from API: {len(all_results):,}")
    
    if not all_results:
        print("No annotations retrieved.")
        return
    
    # Merge with existing dataframe
    for vid, ann in all_results.items():
        if vid in df.index:
            for col, val in ann.items():
                if col in df.columns and pd.isna(df.at[vid, col]):
                    df.at[vid, col] = val
    
    # Save updated VEP file (overwrite original)
    df.to_parquet(vep_path)
    print(f"Saved updated VEP to: {vep_path}")
    
    # Report coverage improvement
    for col in KEY_VEP_COLS:
        if col in df.columns:
            cov_after = df[col].notna().mean() * 100
            print(f"  {col}: {cov_after:.1f}% covered")


def main():
    benchmarks = [
        "humsavar_vep.parquet",
        "mave_independent_vep.parquet",
        "gnomad_benign_vep.parquet",
    ]
    
    for fname in benchmarks:
        process_benchmark(fname)
    
    print("\n✅ All benchmarks processed!")


if __name__ == "__main__":
    main()
