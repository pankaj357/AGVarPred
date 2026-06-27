#!/usr/bin/env python3
"""
Retry VEP annotation for variants that failed in the initial Grimm2015 VEP run.
Loads existing grimm2015_vep.parquet, identifies missing annotations, and retries.
"""

import json
import time
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import requests

ROOT = Path("external_validation")
BENCH_PATH = ROOT / "benchmarks" / "benchmark_grimm2015.csv"
VEP_OUT = ROOT / "vep_preprocessed" / "grimm2015_vep.parquet"

BATCH_SIZE = 50  # Smaller batches for retry
MAX_RETRIES = 5
API_URL = "https://rest.ensembl.org/vep/homo_sapiens/region?canonical=1"


def variant_to_vcf_string(chrom, pos, ref, alt):
    chrom = str(chrom).replace("chr", "")
    return f"{chrom} {pos} . {ref} {alt} . . ."


def parse_api_response(data):
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
        results[variant_id] = {
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
    return results


def fetch_batch(batch_variants, attempt=0):
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
            sleep_time = min(2 ** attempt, 60)
            print(f"   Rate limited, sleeping {sleep_time}s...")
            time.sleep(sleep_time)
            return fetch_batch(batch_variants, attempt + 1)
        else:
            if attempt < MAX_RETRIES:
                sleep_time = min(2 ** attempt + 5, 120)
                print(f"   API error {resp.status_code}, retrying in {sleep_time}s...")
                time.sleep(sleep_time)
                return fetch_batch(batch_variants, attempt + 1)
            print(f"   API error {resp.status_code} after retries: {resp.text[:100]}")
            return {}
    except Exception as e:
        if attempt < MAX_RETRIES:
            sleep_time = min(2 ** attempt + 5, 120)
            print(f"   Request error: {e}, retrying in {sleep_time}s...")
            time.sleep(sleep_time)
            return fetch_batch(batch_variants, attempt + 1)
        print(f"   Request error after retries: {e}")
        return {}


def main():
    print("=" * 60)
    print("Grimm2015 VEP Retry — Filling missing annotations")
    print("=" * 60)

    df = pd.read_csv(BENCH_PATH)
    variants = []
    for _, row in df.iterrows():
        chrom = str(row["chrom"]).replace("chr", "")
        chrom = "chr" + chrom
        pos = str(row["pos"]).replace(".0", "")
        ref = str(row["ref"]).upper()
        alt = str(row["alt"]).upper()
        vid = f"{chrom}_{pos}_{ref}_{alt}"
        variants.append((vid, chrom, pos, ref, alt))

    # Load existing VEP
    existing_vep = pd.read_parquet(VEP_OUT)
    print(f"Existing VEP rows: {len(existing_vep):,}")
    print(f"Existing with annotation: {existing_vep['vep_Consequence'].notna().sum():,}")

    missing = [v for v in variants if v[0] not in existing_vep.index or pd.isna(existing_vep.loc[v[0], "vep_Consequence"])]
    print(f"Missing annotations: {len(missing):,}")

    if not missing:
        print("No missing annotations. Done!")
        return

    all_results = {}
    n_batches = (len(missing) + BATCH_SIZE - 1) // BATCH_SIZE

    for i in tqdm(range(n_batches), desc="Retrying missing VEP batches"):
        start = i * BATCH_SIZE
        end = min((i + 1) * BATCH_SIZE, len(missing))
        batch = missing[start:end]
        results = fetch_batch(batch)
        all_results.update(results)
        time.sleep(1.0)

    print(f"\nNew annotations retrieved: {len(all_results):,} / {len(missing):,}")

    # Update existing VEP
    for vid, result in all_results.items():
        for col, val in result.items():
            existing_vep.loc[vid, col] = val

    existing_vep.to_parquet(VEP_OUT)
    print(f"Saved updated VEP: {VEP_OUT}")
    print(f"Total with annotation: {existing_vep['vep_Consequence'].notna().sum():,} / {len(existing_vep):,}")
    for col in ["vep_SIFT_score", "vep_PolyPhen_score", "vep_IMPACT_score"]:
        if col in existing_vep.columns:
            cov = existing_vep[col].notna().mean() * 100
            print(f"  {col}: {cov:.1f}% covered")

    print("\n✅ Retry complete!")


if __name__ == "__main__":
    main()
