#!/usr/bin/env python3
import pandas as pd
import numpy as np
import requests
import time
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, f1_score, matthews_corrcoef

import os

ROOT = Path(os.environ.get("PROJECT_ROOT", "."))
EXT_DIR = ROOT / "external_validation"
BENCH_DIR = EXT_DIR / "benchmarks"
OUTDIR = EXT_DIR / "results/benchmark_tools"
OUTDIR.mkdir(parents=True, exist_ok=True)

BATCH_SIZE = 1000
RATE_LIMIT_DELAY = 0.3

BENCHMARKS = {
    "humsavar": "benchmark_independent_humsavar.csv",
    "mave_independent": "benchmark_mave_independent.csv",
    "gnomad_benign": "benchmark_gnomad_benign.csv",
    "vip": "benchmark_vip.csv",
}

def row_to_hgvs(row):
    chrom = str(row['chrom'])
    if not chrom.startswith('chr'):
        chrom = f"chr{chrom}"
    return f"{chrom}:g.{int(row['pos'])}{row['ref']}>{row['alt']}"

def fetch_scores_batch(hgvs_list):
    url = "https://myvariant.info/v1/variant"
    try:
        resp = requests.post(url, json={"ids": hgvs_list}, timeout=60)
        if resp.status_code == 200:
            return resp.json()
        elif resp.status_code == 429:
            time.sleep(5)
            resp = requests.post(url, json={"ids": hgvs_list}, timeout=60)
            return resp.json() if resp.status_code == 200 else []
        return []
    except Exception as e:
        print(f"    API error: {e}")
        return []

def parse_scores(data_list):
    results = {}
    for entry in data_list:
        query = entry.get("query")
        if not query:
            continue
        if "notfound" in entry:
            results[query] = {"cadd": np.nan, "revel": np.nan, "alphamissense": np.nan}
            continue
        
        cadd_phred = entry.get("cadd", {}).get("phred") if isinstance(entry.get("cadd"), dict) else None
        
        revel_score = None
        db = entry.get("dbnsfp")
        if isinstance(db, dict) and "revel" in db:
            revel = db["revel"]
            if isinstance(revel, dict) and "score" in revel:
                scores = revel["score"]
                if isinstance(scores, list) and scores:
                    revel_score = np.mean([s for s in scores if s is not None])
                elif isinstance(scores, (int, float)):
                    revel_score = scores
        
        am_score = None
        if isinstance(db, dict) and "alphamissense" in db:
            am = db["alphamissense"]
            if isinstance(am, dict) and "score" in am:
                scores = am["score"]
                if isinstance(scores, list) and scores:
                    am_score = np.mean([s for s in scores if s is not None])
                elif isinstance(scores, (int, float)):
                    am_score = scores
        
        results[query] = {"cadd": cadd_phred, "revel": revel_score, "alphamissense": am_score}
    return results

def get_all_scores(hgvs_list):
    batches = [hgvs_list[i:i+BATCH_SIZE] for i in range(0, len(hgvs_list), BATCH_SIZE)]
    all_scores = {}
    print(f"Fetching scores for {len(hgvs_list)} variants in {len(batches)} batches...")
    for i, batch in enumerate(batches):
        data = fetch_scores_batch(batch)
        scores = parse_scores(data)
        all_scores.update(scores)
        if (i + 1) % 5 == 0 or i == len(batches) - 1:
            print(f"  Batch {i+1}/{len(batches)}: {len(all_scores)} scores retrieved")
        time.sleep(RATE_LIMIT_DELAY)
    return all_scores

def compute_metrics(y_true, y_scores, threshold=0.5):
    mask = ~np.isnan(y_scores)
    y_true_f = y_true[mask]
    y_scores_f = y_scores[mask]
    if len(y_true_f) == 0 or len(np.unique(y_true_f)) < 2:
        return {"n_scored": int(mask.sum()), "n_missing": int((~mask).sum())}
    
    y_pred = (y_scores_f >= threshold).astype(int)
    return {
        "n_scored": int(len(y_true_f)),
        "n_missing": int(len(y_true) - len(y_true_f)),
        "roc_auc": roc_auc_score(y_true_f, y_scores_f),
        "pr_auc": average_precision_score(y_true_f, y_scores_f),
        "accuracy": accuracy_score(y_true_f, y_pred),
        "f1": f1_score(y_true_f, y_pred, zero_division=0),
        "mcc": matthews_corrcoef(y_true_f, y_pred),
        "sensitivity": float(np.sum((y_pred == 1) & (y_true_f == 1)) / max(1, np.sum(y_true_f == 1))),
        "specificity": float(np.sum((y_pred == 0) & (y_true_f == 0)) / max(1, np.sum(y_true_f == 0))),
    }

def fmt(val):
    return f"{val:.4f}" if isinstance(val, (int, float)) and not np.isnan(val) else "N/A"

def process_benchmark(benchmark_name, bench_file):
    print(f"\n{'='*70}")
    print(f"Processing: {benchmark_name}")
    print(f"{'='*70}")
    
    df = pd.read_csv(bench_file)
    print(f"Loaded {len(df)} variants")
    if len(df) == 0 or "label" not in df.columns:
        print("Empty or no labels, skipping")
        return None
    
    df["hgvs"] = df.apply(row_to_hgvs, axis=1)
    hgvs_list = df["hgvs"].tolist()
    scores = get_all_scores(hgvs_list)
    
    df["cadd"] = df["hgvs"].map(lambda x: scores.get(x, {}).get("cadd", np.nan))
    df["revel"] = df["hgvs"].map(lambda x: scores.get(x, {}).get("revel", np.nan))
    df["alphamissense"] = df["hgvs"].map(lambda x: scores.get(x, {}).get("alphamissense", np.nan))
    
    score_file = OUTDIR / f"{benchmark_name}_tool_scores.csv"
    df.to_csv(score_file, index=False)
    print(f"Saved scores to {score_file}")
    print(f"Coverage: CADD={df['cadd'].notna().sum()}, REVEL={df['revel'].notna().sum()}, AM={df['alphamissense'].notna().sum()}")
    
    y_true = df["label"].values
    results = {"benchmark": benchmark_name, "n_total": len(df)}
    
    cadd_m = compute_metrics(y_true, df["cadd"].values, threshold=10)
    for k, v in cadd_m.items():
        results[f"cadd_{k}"] = v
    
    revel_m = compute_metrics(y_true, df["revel"].values, threshold=0.5)
    for k, v in revel_m.items():
        results[f"revel_{k}"] = v
    
    am_m = compute_metrics(y_true, df["alphamissense"].values, threshold=0.564)
    for k, v in am_m.items():
        results[f"alphamissense_{k}"] = v
    
    print(f"\nResults for {benchmark_name}:")
    print(f"  CADD:      ROC-AUC={fmt(results.get('cadd_roc_auc'))}, PR-AUC={fmt(results.get('cadd_pr_auc'))}, n_scored={results.get('cadd_n_scored', 'N/A')}")
    print(f"  REVEL:     ROC-AUC={fmt(results.get('revel_roc_auc'))}, PR-AUC={fmt(results.get('revel_pr_auc'))}, n_scored={results.get('revel_n_scored', 'N/A')}")
    print(f"  AlphaMiss: ROC-AUC={fmt(results.get('alphamissense_roc_auc'))}, PR-AUC={fmt(results.get('alphamissense_pr_auc'))}, n_scored={results.get('alphamissense_n_scored', 'N/A')}")
    
    return results

def main():
    all_results = []
    for bench_name, bench_file in BENCHMARKS.items():
        bench_path = BENCH_DIR / bench_file
        if not bench_path.exists():
            print(f"Benchmark not found: {bench_path}")
            continue
        result = process_benchmark(bench_name, bench_path)
        if result:
            all_results.append(result)
    
    summary_df = pd.DataFrame(all_results)
    summary_file = OUTDIR / "benchmark_tools_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"\n{'='*70}")
    print(f"Summary saved to {summary_file}")
    print(f"{'='*70}")
    
    print("\n## SUMMARY TABLE")
    print(f"{'Benchmark':<20} {'CADD AUC':<12} {'REVEL AUC':<12} {'AlphaMiss AUC':<15} {'CADD PR':<12} {'REVEL PR':<12} {'AlphaMiss PR':<15}")
    print("-" * 100)
    for _, row in summary_df.iterrows():
        print(f"{row['benchmark']:<20} {fmt(row.get('cadd_roc_auc')):<12} {fmt(row.get('revel_roc_auc')):<12} {fmt(row.get('alphamissense_roc_auc')):<15} {fmt(row.get('cadd_pr_auc')):<12} {fmt(row.get('revel_pr_auc')):<12} {fmt(row.get('alphamissense_pr_auc')):<15}")

if __name__ == "__main__":
    main()
