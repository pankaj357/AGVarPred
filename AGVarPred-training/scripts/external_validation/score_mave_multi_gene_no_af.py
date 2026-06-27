#!/usr/bin/env python3
"""Score MAVE Multi-Gene benchmark with the no-AF regularized model."""

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
sys.path.insert(0, str(Path(__file__).parent))

from scoring_core import (
    ROOT, load_extracted_features, load_pipeline, compute_metrics,
    clean_name
)

EXT_DIR = ROOT / "external_validation"
FEATURE_DIR = EXT_DIR / "processing/features/mave_multi_gene"
BENCHMARK_PATH = EXT_DIR / "benchmarks/benchmark_mave_multi_gene.csv"
MODEL_PATH = ROOT / "model_6_minus_af_output/Model_1_no_AF_pipeline.pkl"
VEP_PATH = EXT_DIR / "vep/mave_multi_gene_vep.parquet"
OUTDIR = EXT_DIR / "results/mave_multi_gene/regularized_no_af"
BENCHMARK_NAME = "mave_multi_gene"


def load_mave_labels(benchmark_path):
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
    print(f"   Pathogenic (label=1): {bench['label'].sum():,}")
    print(f"   Benign (label=0): {(bench['label'] == 0).sum():,}")
    print(f"   Score range: [{bench['score'].min():.3f}, {bench['score'].max():.3f}]")
    return bench


def score_mave_no_af(feature_df, pipe, labels_df, outdir):
    print("\n🔄 Scoring...")

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
    score_map = dict(zip(labels_df["variant_id"], labels_df["score"]))

    pred_df["true_label"] = pred_df["variant_id"].map(label_map)
    pred_df["gene"] = pred_df["variant_id"].map(gene_map)
    pred_df["mave_score"] = pred_df["variant_id"].map(score_map)

    n_unlabeled = pred_df["true_label"].isna().sum()
    if n_unlabeled > 0:
        print(f"   ⚠ {n_unlabeled:,} variants could not be matched to labels")
    pred_df = pred_df.dropna(subset=["true_label"]).copy()
    pred_df["true_label"] = pred_df["true_label"].astype(int)

    y_true = pred_df["true_label"].values
    y_prob = pred_df["prob_calibrated"].values
    y_pred = pred_df["prediction"].values

    print(f"\n📊 SCORED VARIANTS: {len(pred_df):,}")
    print(f"   Pathogenic: {y_true.sum():,}")
    print(f"   Benign:     {(y_true == 0).sum():,}")

    all_metrics = compute_metrics(y_true, y_prob, y_pred, prefix="overall_")

    # Correlation with continuous MAVE scores
    print("\n🔄 Computing correlations with MAVE scores...")
    valid_score = pred_df.dropna(subset=["mave_score"])
    corr_results = {}
    if len(valid_score) > 1:
        rho, rho_p = spearmanr(valid_score["prob_calibrated"], valid_score["mave_score"])
        r, r_p = pearsonr(valid_score["prob_calibrated"], valid_score["mave_score"])
        corr_results["score_spearman_rho"] = float(rho)
        corr_results["score_spearman_p"] = float(rho_p)
        corr_results["score_pearson_r"] = float(r)
        corr_results["score_pearson_p"] = float(r_p)
        print(f"   Spearman ρ = {rho:.4f} (p={rho_p:.2e})")
        print(f"   Pearson r  = {r:.4f} (p={r_p:.2e})")

    # Classification subgroups (no AF subgroups)
    subgroups = {}
    for gene in sorted(pred_df["gene"].dropna().unique()):
        subgroups[f"gene_{gene}"] = pred_df["gene"] == gene

    subgroup_results = []
    for name, mask in subgroups.items():
        sub = pred_df[mask]
        if len(sub) == 0:
            continue
        sm = compute_metrics(
            sub["true_label"].values,
            sub["prob_calibrated"].values,
            sub["prediction"].values,
            prefix=""
        )
        sm["subgroup"] = name
        sm["n_variants"] = len(sub)
        subgroup_results.append(sm)

    characterization = {
        "benchmark": "MAVE Multi-Gene (no AF)",
        "n_total": int(len(pred_df)),
        "n_pathogenic": int(y_true.sum()),
        "n_benign": int((y_true == 0).sum()),
        "n_genes": int(pred_df["gene"].nunique()),
    }
    characterization.update(corr_results)

    return pred_df, all_metrics, subgroup_results, characterization


def save_mave_results(pred_df, metrics, subgroup_results, characterization, outdir):
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

    if subgroup_results:
        sub_df = pd.DataFrame(subgroup_results)
        sub_path = outdir / f"{BENCHMARK_NAME}_subgroup_metrics.csv"
        sub_df.to_csv(sub_path, index=False)
        print(f"💾 Subgroup metrics saved: {sub_path}")

        print("\n" + "="*70)
        print("📊 SUBGROUP METRICS")
        print("="*70)
        display_cols = ["subgroup", "n_variants", "roc_auc", "pr_auc", "accuracy", "f1", "mcc", "sensitivity", "specificity"]
        display_cols = [c for c in display_cols if c in sub_df.columns]
        print(sub_df[display_cols].to_string(index=False))

    char_path = outdir / "benchmark_characterization.json"
    with open(char_path, "w") as fh:
        json.dump(characterization, fh, indent=2)
    print(f"\n💾 Characterization saved: {char_path}")

    print("\n" + "="*70)
    print("📊 BENCHMARK CHARACTERIZATION")
    print("="*70)
    print(f"  Total variants:        {characterization['n_total']:,}")
    print(f"  Pathogenic:            {characterization['n_pathogenic']:,}")
    print(f"  Benign:                {characterization['n_benign']:,}")
    print(f"  Genes:                 {characterization['n_genes']:,}")
    if "score_spearman_rho" in characterization:
        print(f"  Spearman ρ (score):    {characterization['score_spearman_rho']:.4f}")
    print("="*70)

    print("\n" + "="*70)
    print("📊 OVERALL METRICS")
    print("="*70)
    for k in ["overall_roc_auc", "overall_pr_auc", "overall_accuracy", "overall_f1", "overall_mcc", "overall_brier_score", "overall_sensitivity", "overall_specificity"]:
        if k in metrics:
            val = metrics[k]
            val_str = f"{val:.4f}" if val is not None else "N/A"
            print(f"  {k.replace('overall_', '').replace('_', ' ').title():15s}: {val_str}")
    print("="*70)


if __name__ == "__main__":
    print("="*70)
    print(f"MAVE MULTI-GENE BENCHMARK SCORING — Model: regularized_no_af")
    print("="*70)

    feature_df = load_extracted_features(FEATURE_DIR)
    pipe = load_pipeline(MODEL_PATH)
    labels_df = load_mave_labels(BENCHMARK_PATH)

    pred_df, metrics, subgroup_results, characterization = score_mave_no_af(
        feature_df, pipe, labels_df, OUTDIR
    )
    save_mave_results(pred_df, metrics, subgroup_results, characterization, OUTDIR)

    print("\n✅ Scoring complete!")
    print(f"   Output directory: {OUTDIR}/")
