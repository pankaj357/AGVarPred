#!/usr/bin/env python3
"""
Fetch complete VEP annotations via Ensembl REST API for ClinVar 3-star subset.
This is needed because the raw VEP file was accidentally deleted during cleanup.
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
BENCH_PATH = ROOT / "benchmarks" / "benchmark_clinvar_3star_subset.csv"
VEP_OUT = ROOT / "clinvar_3star_subset_vep.parquet"
PREPROCESSED_OUT = ROOT / "vep_preprocessed" / "clinvar_3star_subset_vep.parquet"

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

        sift_score = canonical_tc.get("sift_score")
        sift_pred = canonical_tc.get("sift_prediction")
        polyphen_score = canonical_tc.get("polyphen_score")
        polyphen_pred = canonical_tc.get("polyphen_prediction")
        protein_pos = canonical_tc.get("protein_start")

        cons_list = consequence_terms
        impact_map = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
        impact_score = impact_map.get(impact, 0) if impact else 0

        result = {
            "vep_SIFT_score": sift_score,
            "vep_SIFT_pred": sift_pred,
            "vep_PolyPhen_score": polyphen_score,
            "vep_PolyPhen_pred": polyphen_pred,
            "vep_IMPACT": impact,
            "vep_IMPACT_score": impact_score,
            "vep_Consequence": "&".join(cons_list) if cons_list else most_severe,
            "vep_LoF": None,  # Not available via REST API
            "vep_is_missense": 1 if "missense_variant" in cons_list else 0,
            "vep_is_synonymous": 1 if "synonymous_variant" in cons_list else 0,
            "vep_is_stop_gained": 1 if "stop_gained" in cons_list else 0,
            "vep_is_frameshift": 1 if "frameshift_variant" in cons_list else 0,
            "vep_is_splice": 1 if any(c in cons_list for c in ["splice_donor_variant", "splice_acceptor_variant", "splice_region_variant"]) else 0,
            "vep_is_LoF_HC": 0,
            "vep_has_SIFT": 1 if sift_score is not None else 0,
            "vep_has_PolyPhen": 1 if polyphen_score is not None else 0,
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
            print(f"   Rate limited, sleeping 10s...")
            time.sleep(10)
            return fetch_batch(batch_variants)  # retry
        else:
            print(f"   API error {resp.status_code}: {resp.text[:200]}")
            return {}
    except Exception as e:
        print(f"   Request error: {e}")
        return {}


def main():
    # Load benchmark to get variant IDs
    bench = pd.read_csv(BENCH_PATH)
    bench["variant_id"] = (
        "chr" + bench["chrom"].astype(str).str.replace("chr", "", case=False)
        + "_" + bench["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
        + "_" + bench["ref"].astype(str).str.upper()
        + "_" + bench["alt"].astype(str).str.upper()
    )
    all_variants = bench["variant_id"].tolist()
    print(f"Total ClinVar variants to fetch: {len(all_variants):,}")

    # Fetch in batches
    all_results = {}
    n_batches = (len(all_variants) + BATCH_SIZE - 1) // BATCH_SIZE

    for i in tqdm(range(n_batches), desc="Fetching VEP batches"):
        start = i * BATCH_SIZE
        end = min((i + 1) * BATCH_SIZE, len(all_variants))
        batch = all_variants[start:end]

        results = fetch_batch(batch)
        all_results.update(results)

        if i < n_batches - 1:
            time.sleep(SLEEP_BETWEEN_BATCHES)

    print(f"Annotations retrieved from API: {len(all_results):,} / {len(all_variants):,}")

    # Build DataFrame
    rows = []
    for vid in all_variants:
        row = {"variant_id": vid}
        row.update(all_results.get(vid, {}))
        rows.append(row)

    df = pd.DataFrame(rows).set_index("variant_id")

    # Save raw VEP
    VEP_OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(VEP_OUT)
    print(f"Saved raw VEP to: {VEP_OUT} ({len(df):,} variants, {len(df.columns)} columns)")

    # Report coverage
    for col in ["vep_SIFT_score", "vep_PolyPhen_score", "vep_IMPACT"]:
        if col in df.columns:
            cov = df[col].notna().mean() * 100
            print(f"  {col}: {cov:.1f}% covered")

    # Now preprocess
    print("\nPreprocessing VEP...")
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from preprocess_vep_for_scoring import preprocess_vep

    df_pre = preprocess_vep(df.copy())
    PREPROCESSED_OUT.parent.mkdir(parents=True, exist_ok=True)
    df_pre.to_parquet(PREPROCESSED_OUT)
    print(f"Saved preprocessed VEP to: {PREPROCESSED_OUT}")
    print(f"Columns: {len(df_pre.columns)}")

    # Check for required model features
    required = [
        "vep_IMPACT_score", "vep_LoF", "vep_PolyPhen_score", "vep_SIFT_score",
        "vep_is_synonymous", "vep_Consequence_synonymous_variant",
        "vep_Consequence_intron_variant", "vep_Consequence_downstream_gene_variant"
    ]
    for col in required:
        if col in df_pre.columns:
            cov = df_pre[col].notna().mean() * 100
            print(f"  {col}: {cov:.1f}% coverage")
        else:
            print(f"  {col}: MISSING")


if __name__ == "__main__":
    main()
