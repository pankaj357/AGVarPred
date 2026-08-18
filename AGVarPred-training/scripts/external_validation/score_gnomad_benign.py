#!/usr/bin/env python3
"""Score gnomAD Common Benign benchmark (benign-only, all labels = 0).

Since this benchmark contains only genuinely benign variants (common in
population, AF > 0.1%), standard classification metrics (ROC-AUC, MCC) are
not applicable. Instead we report:
- False Positive Rate (FPR): fraction incorrectly called pathogenic
- Average predicted probability
- Calibration (mean predicted prob vs true prob = 0)
- Distribution of scores across genes
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from scoring_core import (
    ROOT, load_extracted_features, load_pipeline, clean_name, bootstrap_ci
)

EXT_DIR = ROOT / "external_validation"
FEATURE_DIR = EXT_DIR / "processing/features/gnomad_benign"
BENCHMARK_PATH = EXT_DIR / "benchmarks/benchmark_gnomad_benign.csv"
MODEL_PATH = ROOT / "final_model_output_regularized/final_pipeline.pkl"
VEP_PATH = EXT_DIR / "vep/gnomad_benign_vep.parquet"
OUTDIR = EXT_DIR / "results/gnomad_benign/regularized"
BENCHMARK_NAME = "gnomad_benign"


def load_benign_labels(benchmark_path):
    print("\n🔄 Loading benchmark labels...")
    bench = pd.read_csv(benchmark_path)
    bench["variant_id"] = (
        "chr" + bench["chrom"].astype(str).str.replace("chr", "", case=False)
        + "_" + bench["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
        + "_" + bench["ref"].astype(str).str.upper()
        + "_" + bench["alt"].astype(str).str.upper()
    )
    print(f"   Benchmark variants: {len(bench):,}")
    print(f"   Genes: {bench['gene'].nunique():,}")
    print(f"   All labels = 0 (benign)")
    print(f"   Mean AF: {bench['af'].mean():.4f}")
    return bench


def score_benign(feature_df, pipe, labels_df, outdir):
    print("\n🔄 Scoring...")
    variant_ids = list(feature_df.index)

    # Add gnomAD AF from benchmark labels
    print("   Adding gnomAD AF from benchmark labels...")
    af_map = dict(zip(labels_df["variant_id"], labels_df["af"]))
    feature_df["gnomAD_AF"] = feature_df.index.map(lambda v: af_map.get(v, 0.0))
    feature_df["log10_gnomAD_AF"] = np.log10(feature_df["gnomAD_AF"] + 1e-8)
    feature_df["AF_missing"] = (feature_df["gnomAD_AF"] == 0).astype(int)
    feature_df["is_ultra_rare"] = (feature_df["gnomAD_AF"] < 0.0001).astype(int)

    if VEP_PATH.exists():
        print(f"   Merging VEP features from {VEP_PATH}")
        vep_df = pd.read_parquet(VEP_PATH)
        vep_df.columns = [clean_name(c) for c in vep_df.columns]
        n_vep_dup = vep_df.index.duplicated().sum()
        if n_vep_dup > 0:
            vep_df = vep_df[~vep_df.index.duplicated(keep="first")]
            print(f"   ⚠ Deduplicated VEP: removed {n_vep_dup} duplicate rows")
        feature_df = feature_df.join(vep_df, how="left")

    features = list(pipe["features"])
    missing = [f for f in features if f not in feature_df.columns]
    if missing:
        print(f"   ⚠ Missing features ({len(missing)}): {missing[:10]}...")
        for f in missing:
            feature_df[f] = np.nan

    X = feature_df[features].copy()
    imputer = pipe["imputer"]
    for col in X.columns:
        if col in imputer.index:
            X[col] = X[col].fillna(imputer[col])
        else:
            X[col] = X[col].fillna(0)

    X_scaled = pipe["scaler"].transform(X)
    probs_raw = pipe["model"].predict_proba(X_scaled)[:, 1]
    probs_cal = pipe["calibrator"].predict(probs_raw)
    threshold = float(pipe["threshold"])
    preds = (probs_cal >= threshold).astype(int)

    pred_df = pd.DataFrame({
        "variant_id": feature_df.index,
        "prob_raw": probs_raw,
        "prob_calibrated": probs_cal,
        "prediction": preds,
    })

    label_map = dict(zip(labels_df["variant_id"], labels_df["label"]))
    gene_map = dict(zip(labels_df["variant_id"], labels_df["gene"]))
    af_map = dict(zip(labels_df["variant_id"], labels_df["af"]))

    pred_df["true_label"] = pred_df["variant_id"].map(label_map)
    pred_df["gene"] = pred_df["variant_id"].map(gene_map)
    pred_df["af"] = pred_df["variant_id"].map(af_map)

    n_unlabeled = pred_df["true_label"].isna().sum()
    if n_unlabeled > 0:
        print(f"   ⚠ {n_unlabeled:,} variants could not be matched to labels")
    pred_df = pred_df.dropna(subset=["true_label"]).copy()
    pred_df["true_label"] = pred_df["true_label"].astype(int)

    n_before = len(pred_df)
    pred_df = pred_df.drop_duplicates(subset=["variant_id"], keep="first")
    n_deduped = n_before - len(pred_df)
    if n_deduped > 0:
        print(f"   ⚠ Deduplicated predictions: removed {n_deduped:,} duplicate variant rows")

    y_true = pred_df["true_label"].values
    y_prob = pred_df["prob_calibrated"].values
    y_pred = pred_df["prediction"].values

    print(f"\n📊 SCORED VARIANTS: {len(pred_df):,}")
    print(f"   All benign (label=0)")

    # Benign-only metrics
    n_total = len(pred_df)
    n_fp = int(y_pred.sum())  # predicted as pathogenic
    fpr = n_fp / n_total if n_total > 0 else 0.0
    mean_prob = float(y_prob.mean())
    median_prob = float(np.median(y_prob))

    all_metrics = {
        "overall_n_variants": n_total,
        "overall_n_benign": n_total,
        "overall_n_fp": n_fp,
        "overall_fpr": fpr,
        "overall_mean_prob": mean_prob,
        "overall_median_prob": median_prob,
        "overall_accuracy": float((y_pred == 0).sum() / n_total) if n_total > 0 else 0.0,
    }

    # Bootstrap 95% CI for FPR (requires both classes flag off because all labels are 0)
    fpr_lo, fpr_hi = bootstrap_ci(
        y_true, y_prob, y_pred,
        metric_fn=lambda yt, yp, yd: float(yd.sum()) / len(yd) if len(yd) > 0 else 0.0,
        n_bootstraps=1000,
        random_state=42,
        requires_both_classes=False,
    )
    if fpr_lo is not None:
        all_metrics["overall_fpr_ci_low"] = fpr_lo
        all_metrics["overall_fpr_ci_high"] = fpr_hi

    # Gene-level FPR
    gene_results = []
    for gene, gdf in pred_df.groupby("gene"):
        if len(gdf) < 5:
            continue
        g_fp = int(gdf["prediction"].sum())
        g_fpr = g_fp / len(gdf)
        gene_results.append({
            "gene": gene,
            "n_variants": len(gdf),
            "n_fp": g_fp,
            "fpr": g_fpr,
            "mean_prob": float(gdf["prob_calibrated"].mean()),
        })

    characterization = {
        "benchmark": "gnomAD Common Benign (non-training genes)",
        "n_total": n_total,
        "n_genes": int(pred_df["gene"].nunique()),
        "n_fp": n_fp,
        "fpr": fpr,
        "mean_prob": mean_prob,
        "median_prob": median_prob,
        "mean_af": float(pred_df["af"].mean()),
    }

    return pred_df, all_metrics, gene_results, characterization


def save_benign_results(pred_df, metrics, gene_results, characterization, outdir):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    pred_path = outdir / f"{BENCHMARK_NAME}_predictions.csv"
    pred_df.to_csv(pred_path, index=False)
    print(f"\n💾 Predictions saved: {pred_path}")

    metrics_path = outdir / f"{BENCHMARK_NAME}_metrics.json"
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    print(f"💾 Metrics saved: {metrics_path}")

    summary = pd.DataFrame([metrics])
    summary_path = outdir / f"{BENCHMARK_NAME}_metrics.csv"
    summary.to_csv(summary_path, index=False)
    print(f"💾 Summary CSV saved: {summary_path}")

    if gene_results:
        gene_df = pd.DataFrame(gene_results).sort_values("fpr", ascending=False)
        gene_path = outdir / f"{BENCHMARK_NAME}_gene_fpr.csv"
        gene_df.to_csv(gene_path, index=False)
        print(f"💾 Gene-level FPR saved: {gene_path}")

    char_path = outdir / "benchmark_characterization.json"
    with open(char_path, "w") as fh:
        json.dump(characterization, fh, indent=2)
    print(f"💾 Characterization saved: {char_path}")

    print("\n" + "="*70)
    print("📊 BENCHMARK CHARACTERIZATION")
    print("="*70)
    for k, v in characterization.items():
        if isinstance(v, float):
            print(f"  {k:25s}: {v:.4f}")
        else:
            print(f"  {k:25s}: {v}")
    print("="*70)

    print("\n" + "="*70)
    print("📊 KEY METRICS (Benign-Only Benchmark)")
    print("="*70)
    print(f"  N variants:     {metrics['overall_n_variants']:,}")
    print(f"  False Positives:{metrics['overall_n_fp']:,}")
    fpr_str = f"{metrics['overall_fpr']:.4f}"
    if 'overall_fpr_ci_low' in metrics and 'overall_fpr_ci_high' in metrics:
        fpr_str += f" (95% CI: {metrics['overall_fpr_ci_low']:.4f}–{metrics['overall_fpr_ci_high']:.4f})"
    print(f"  FPR:            {fpr_str}")
    print(f"  Mean Prob:      {metrics['overall_mean_prob']:.4f}")
    print(f"  Median Prob:    {metrics['overall_median_prob']:.4f}")
    print("="*70)


if __name__ == "__main__":
    print("="*70)
    print(f"GNOMAD COMMON BENIGN BENCHMARK SCORING — Model: regularized")
    print("="*70)

    feature_df = load_extracted_features(FEATURE_DIR)
    pipe = load_pipeline(MODEL_PATH)
    labels_df = load_benign_labels(BENCHMARK_PATH)

    pred_df, metrics, gene_results, characterization = score_benign(
        feature_df, pipe, labels_df, OUTDIR
    )
    save_benign_results(pred_df, metrics, gene_results, characterization, OUTDIR)

    print("\n✅ Scoring complete!")
    print(f"   Output directory: {OUTDIR}/")
