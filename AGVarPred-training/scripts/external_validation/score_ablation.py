#!/usr/bin/env python3
"""Generic ablation model scoring script for any model + benchmark combination."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

from scoring_core import (
    ROOT, load_extracted_features, load_pipeline, load_labels,
    load_train_genes, score_benchmark, save_results, compute_metrics, clean_name,
    bootstrap_ci
)

EXT_DIR = ROOT / "external_validation"

BENCHMARKS = {
    "humsavar": {
        "feature_dir": EXT_DIR / "processing/features/humsavar",
        "benchmark_path": EXT_DIR / "benchmarks/benchmark_independent_humsavar.csv",
        "vep_path": EXT_DIR / "vep/humsavar_vep.parquet",
        "outdir_name": "humsavar",
        "result_name": "humsavar",
    },
    "mave_independent": {
        "feature_dir": EXT_DIR / "processing/features/mave_independent",
        "benchmark_path": EXT_DIR / "benchmarks/benchmark_mave_independent.csv",
        "vep_path": EXT_DIR / "vep/mave_independent_vep.parquet",
        "outdir_name": "mave_independent",
        "result_name": "mave_independent",
    },
    "gnomad_benign": {
        "feature_dir": EXT_DIR / "processing/features/gnomad_benign",
        "benchmark_path": EXT_DIR / "benchmarks/benchmark_gnomad_benign.csv",
        "vep_path": EXT_DIR / "vep/gnomad_benign_vep.parquet",
        "outdir_name": "gnomad_benign",
        "result_name": "gnomad_benign",
    },
    "vip": {
        "feature_dir": EXT_DIR / "processing/features/vip",
        "benchmark_path": EXT_DIR / "benchmarks/benchmark_vip.csv",
        "vep_path": EXT_DIR / "vep_preprocessed/vip_vep.parquet",
        "outdir_name": "vip",
        "result_name": "vip",
    },
    "grimm2015": {
        "feature_dir": EXT_DIR / "processing/features/grimm2015",
        "benchmark_path": EXT_DIR / "benchmarks/benchmark_grimm2015.csv",
        "vep_path": EXT_DIR / "vep_preprocessed/grimm2015_vep.parquet",
        "outdir_name": "grimm2015",
        "result_name": "grimm2015",
    },
    "dvd": {
        "feature_dir": EXT_DIR / "processing/features/dvd",
        "benchmark_path": EXT_DIR / "benchmarks/benchmark_dvd.csv",
        "vep_path": EXT_DIR / "vep_preprocessed/dvd_vep.parquet",
        "outdir_name": "dvd",
        "result_name": "dvd",
    },
}


def _load_mave_labels(benchmark_path):
    bench = pd.read_csv(benchmark_path)
    bench["variant_id"] = (
        "chr" + bench["chrom"].astype(str).str.replace("chr", "", case=False)
        + "_" + bench["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
        + "_" + bench["ref"].astype(str).str.upper()
        + "_" + bench["alt"].astype(str).str.upper()
    )
    return bench


def _load_benign_labels(benchmark_path):
    bench = pd.read_csv(benchmark_path)
    bench["variant_id"] = (
        "chr" + bench["chrom"].astype(str).str.replace("chr", "", case=False)
        + "_" + bench["pos"].astype(str).str.replace(r"\.0$", "", regex=True)
        + "_" + bench["ref"].astype(str).str.upper()
        + "_" + bench["alt"].astype(str).str.upper()
    )
    return bench


def score_mave_ablation(feature_df, pipe, labels_df, train_genes, outdir, has_af=True, vep_path=None):
    """MAVE scoring via shared score_benchmark plus continuous-score correlations."""
    score_map = dict(zip(labels_df["variant_id"], labels_df["score"]))

    pred_df, all_metrics, subgroup_results, characterization = score_benchmark(
        feature_df, pipe, labels_df, train_genes, outdir,
        has_af=has_af, vep_path=vep_path
    )

    pred_df["mave_score"] = pred_df["variant_id"].map(score_map)

    print("\n🔄 Computing correlations with MAVE scores...")
    valid_score = pred_df.dropna(subset=["mave_score"])
    corr_results = {}
    if len(valid_score) > 1:
        rho, rho_p = spearmanr(valid_score["prob_calibrated"], valid_score["mave_score"])
        r, r_p = pearsonr(valid_score["prob_calibrated"], valid_score["mave_score"])
        corr_results = {
            "score_spearman_rho": float(rho),
            "score_spearman_p": float(rho_p),
            "score_pearson_r": float(r),
            "score_pearson_p": float(r_p),
            "n_with_score": int(len(valid_score)),
        }
        print(f"   Spearman ρ = {rho:.4f} (p={rho_p:.2e})")
        print(f"   Pearson r  = {r:.4f} (p={r_p:.2e})")
    else:
        print("   ⚠ Not enough variants with MAVE scores for correlation")

    characterization["benchmark"] = "MAVE Independent (non-training genes)"
    characterization.update(corr_results)

    return pred_df, all_metrics, subgroup_results, characterization


def score_benign_ablation(feature_df, pipe, labels_df, outdir, vep_path=None):
    print("\n🔄 Scoring...")

    if vep_path and Path(vep_path).exists():
        print(f"   Merging VEP features from {vep_path}")
        vep_df = pd.read_parquet(vep_path)
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

    n_total = len(pred_df)
    n_fp = int(y_pred.sum())
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


def save_benign_results_ablation(pred_df, metrics, gene_results, characterization, outdir, benchmark_name):
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    pred_path = outdir / f"{benchmark_name}_predictions.csv"
    pred_df.to_csv(pred_path, index=False)
    print(f"\n💾 Predictions saved: {pred_path}")

    metrics_path = outdir / f"{benchmark_name}_metrics.json"
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    print(f"💾 Metrics saved: {metrics_path}")

    summary = pd.DataFrame([metrics])
    summary_path = outdir / f"{benchmark_name}_metrics.csv"
    summary.to_csv(summary_path, index=False)
    print(f"💾 Summary CSV saved: {summary_path}")

    if gene_results:
        gene_df = pd.DataFrame(gene_results).sort_values("fpr", ascending=False)
        gene_path = outdir / f"{benchmark_name}_gene_fpr.csv"
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, help="Path to pipeline .pkl")
    parser.add_argument("--benchmark", required=True, choices=list(BENCHMARKS.keys()))
    parser.add_argument("--model-name", required=True, help="Model name for output dir")
    parser.add_argument("--has-af", action="store_true", help="Query gnomAD AF")
    args = parser.parse_args()

    cfg = BENCHMARKS[args.benchmark]
    feature_df = load_extracted_features(cfg["feature_dir"])
    pipe = load_pipeline(Path(args.model))
    train_genes = load_train_genes()
    outdir = EXT_DIR / f"results/{cfg['outdir_name']}/{args.model_name}"

    print("="*70)
    print(f"{args.benchmark.upper()} BENCHMARK SCORING — Model: {args.model_name}")
    print("="*70)

    if args.benchmark == "mave_independent":
        labels_df = _load_mave_labels(cfg["benchmark_path"])
        pred_df, metrics, subgroup_results, characterization = score_mave_ablation(
            feature_df, pipe, labels_df, train_genes, outdir, has_af=args.has_af, vep_path=str(cfg["vep_path"])
        )
        save_results(pred_df, metrics, subgroup_results, characterization, outdir, cfg["result_name"])
        corr_results = {k: v for k, v in characterization.items()
                        if k.startswith("score_") or k == "n_with_score"}
        if corr_results:
            corr_path = Path(outdir) / f"{cfg['result_name']}_mave_correlation.json"
            with open(corr_path, "w") as fh:
                json.dump(corr_results, fh, indent=2)
            print(f"💾 MAVE correlation saved: {corr_path}")
    elif args.benchmark == "gnomad_benign":
        labels_df = _load_benign_labels(cfg["benchmark_path"])
        pred_df, metrics, gene_results, characterization = score_benign_ablation(
            feature_df, pipe, labels_df, outdir, vep_path=str(cfg["vep_path"])
        )
        save_benign_results_ablation(pred_df, metrics, gene_results, characterization, outdir, cfg["result_name"])
    else:
        labels_df = load_labels(cfg["benchmark_path"])
        pred_df, metrics, subgroup_results, characterization = score_benchmark(
            feature_df, pipe, labels_df, train_genes, outdir,
            has_af=args.has_af, vep_path=str(cfg["vep_path"])
        )
        save_results(pred_df, metrics, subgroup_results, characterization, outdir, cfg["result_name"])

    print("\n✅ Scoring complete!")
    print(f"   Output directory: {outdir}/")


if __name__ == "__main__":
    main()
