#!/usr/bin/env python3
"""Score PanelApp benign benchmark (benign-only, all labels = 0) across models."""

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
FEATURE_DIR = EXT_DIR / "processing/features/panelapp_benign"
BENCHMARK_PATH = EXT_DIR / "benchmarks/benchmark_panelapp_benign_conf3.csv"
VEP_PATH = EXT_DIR / "vep/panelapp_benign_vep.parquet"
BENCHMARK_NAME = "panelapp_benign"

MODELS = [
    ("regularized", ROOT / "final_model_output_regularized/final_pipeline.pkl", EXT_DIR / "results/panelapp_benign/regularized"),
    ("no_af", ROOT / "model_6_minus_af_output/Model_1_no_AF_pipeline.pkl", EXT_DIR / "results/panelapp_benign/no_af"),
    ("af_only", ROOT / "ablation_feature_groups_output/AF_only_pipeline.pkl", EXT_DIR / "results/panelapp_benign/af_only"),
    ("vep_only", ROOT / "ablation_feature_groups_output/VEP_only_pipeline.pkl", EXT_DIR / "results/panelapp_benign/vep_only"),
    ("af_plus_vep", ROOT / "ablation_feature_groups_output/AF_plus_VEP_pipeline.pkl", EXT_DIR / "results/panelapp_benign/af_plus_vep"),
    ("alphagenome_only", ROOT / "ablation_feature_groups_output/AlphaGenome_only_pipeline.pkl", EXT_DIR / "results/panelapp_benign/alphagenome_only"),
]


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


def score_benign(feature_df, pipe, labels_df, outdir, model_name):
    print(f"\n🔄 Scoring with {model_name}...")

    # Add gnomAD AF from benchmark labels
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

    n_total = len(pred_df)
    n_fp = int(y_pred.sum())
    fpr = n_fp / n_total if n_total > 0 else 0.0
    mean_prob = float(y_prob.mean())
    median_prob = float(np.median(y_prob))

    metrics = {
        "model": model_name,
        "n_variants": n_total,
        "n_fp": n_fp,
        "fpr": fpr,
        "mean_prob": mean_prob,
        "median_prob": median_prob,
        "accuracy": float((y_pred == 0).sum() / n_total) if n_total > 0 else 0.0,
    }

    fpr_lo, fpr_hi = bootstrap_ci(
        y_true, y_prob, y_pred,
        metric_fn=lambda yt, yp, yd: float(yd.sum()) / len(yd) if len(yd) > 0 else 0.0,
        n_bootstraps=1000,
        random_state=42,
        requires_both_classes=False,
    )
    if fpr_lo is not None:
        metrics["fpr_ci_low"] = fpr_lo
        metrics["fpr_ci_high"] = fpr_hi

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
        "benchmark": "PanelApp Benign (non-training disease genes)",
        "model": model_name,
        "n_total": n_total,
        "n_genes": int(pred_df["gene"].nunique()),
        "n_fp": n_fp,
        "fpr": fpr,
        "mean_prob": mean_prob,
        "median_prob": median_prob,
        "mean_af": float(pred_df["af"].mean()),
    }

    return pred_df, metrics, gene_results, characterization


def save_results(pred_df, metrics, gene_results, characterization, outdir):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    pred_path = outdir / f"{BENCHMARK_NAME}_predictions.csv"
    pred_df.to_csv(pred_path, index=False)

    metrics_path = outdir / f"{BENCHMARK_NAME}_metrics.json"
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)

    summary = pd.DataFrame([metrics])
    summary_path = outdir / f"{BENCHMARK_NAME}_metrics.csv"
    summary.to_csv(summary_path, index=False)

    if gene_results:
        gene_df = pd.DataFrame(gene_results).sort_values("fpr", ascending=False)
        gene_path = outdir / f"{BENCHMARK_NAME}_gene_fpr.csv"
        gene_df.to_csv(gene_path, index=False)

    char_path = outdir / "benchmark_characterization.json"
    with open(char_path, "w") as fh:
        json.dump(characterization, fh, indent=2)

    print(f"\n📊 {characterization['model']} RESULTS")
    print(f"   N variants:      {metrics['n_variants']:,}")
    print(f"   False Positives: {metrics['n_fp']:,}")
    fpr_str = f"{metrics['fpr']:.4f}"
    if "fpr_ci_low" in metrics and "fpr_ci_high" in metrics:
        fpr_str += f" (95% CI: {metrics['fpr_ci_low']:.4f}–{metrics['fpr_ci_high']:.4f})"
    print(f"   FPR:             {fpr_str}")
    print(f"   Mean Prob:       {metrics['mean_prob']:.4f}")
    print(f"   Median Prob:     {metrics['median_prob']:.4f}")
    print(f"   Output:          {outdir}")


def main():
    print("=" * 70)
    print("PANELAPP BENIGN BENCHMARK SCORING")
    print("=" * 70)

    feature_df = load_extracted_features(FEATURE_DIR)
    labels_df = load_benign_labels(BENCHMARK_PATH)

    all_metrics = []
    for model_name, model_path, outdir in MODELS:
        if not model_path.exists():
            print(f"\n⚠ Model not found: {model_path}, skipping {model_name}")
            continue
        try:
            pipe = load_pipeline(model_path)
            pred_df, metrics, gene_results, characterization = score_benign(
                feature_df.copy(), pipe, labels_df, outdir, model_name
            )
            save_results(pred_df, metrics, gene_results, characterization, outdir)
            all_metrics.append(metrics)
        except Exception as e:
            print(f"\n❌ Error scoring {model_name}: {e}")
            import traceback
            traceback.print_exc()

    # Save combined summary
    summary_dir = EXT_DIR / "results/panelapp_benign"
    summary_dir.mkdir(exist_ok=True, parents=True)
    summary_df = pd.DataFrame(all_metrics)
    summary_df.to_csv(summary_dir / "all_models_summary.csv", index=False)
    print("\n" + "=" * 70)
    print("ALL MODELS SUMMARY")
    print("=" * 70)
    print(summary_df.to_string(index=False))
    print("=" * 70)
    print(f"\n✅ Scoring complete! Output: {summary_dir}/")


if __name__ == "__main__":
    main()
