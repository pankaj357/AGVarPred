"""
Shared scoring utilities for external validation benchmarks.
"""

import os
import re
import json
import subprocess
import joblib
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from sklearn.metrics import (
    roc_auc_score, average_precision_score, accuracy_score,
    precision_score, recall_score, f1_score, matthews_corrcoef,
    confusion_matrix, brier_score_loss
)

warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
GNOMAD_VCF = os.environ.get(
    "GNOMAD_VCF",
    "external_data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
)
TRAIN_GENES_PATH = ROOT / "train_gene_set.json"


def clean_name(s):
    s = str(s)
    s = s.replace(":", "_").replace("-", "_").replace(" ", "_")
    s = re.sub(r'[^A-Za-z0-9_]', '', s)
    s = re.sub(r'_+', '_', s)
    return s


def bootstrap_ci(y_true, y_prob, y_pred, metric_fn, n_bootstraps=1000, random_state=42, requires_both_classes=True):
    """Compute 95% bootstrap confidence interval for a metric."""
    rng = np.random.RandomState(random_state)
    n = len(y_true)
    if n < 10:
        return None, None
    scores = []
    for _ in range(n_bootstraps):
        idx = rng.choice(n, size=n, replace=True)
        # Require at least one sample from each class and non-trivial predictions
        if requires_both_classes and len(set(y_true[idx])) < 2:
            continue
        try:
            scores.append(metric_fn(y_true[idx], y_prob[idx], y_pred[idx]))
        except Exception:
            continue
    if len(scores) < 100:
        return None, None
    lo = float(np.percentile(scores, 2.5))
    hi = float(np.percentile(scores, 97.5))
    return lo, hi


def compute_metrics(y_true, y_prob, y_pred, prefix="", include_ci=False):
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    y_pred = np.asarray(y_pred)

    if len(y_true) == 0:
        return {}

    cm = confusion_matrix(y_true, y_pred)
    if cm.size == 1:
        if len(set(y_true)) == 1:
            if y_true[0] == 1:
                tn, fp, fn, tp = 0, 0, 0, len(y_true)
            else:
                tn, fp, fn, tp = len(y_true), 0, 0, 0
        else:
            tn = fp = fn = tp = 0
    else:
        tn, fp, fn, tp = cm.ravel()

    metrics = {
        f"{prefix}n_variants": int(len(y_true)),
        f"{prefix}n_pathogenic": int(y_true.sum()),
        f"{prefix}n_benign": int((y_true == 0).sum()),
        f"{prefix}roc_auc": float(roc_auc_score(y_true, y_prob)) if len(set(y_true)) > 1 else None,
        f"{prefix}pr_auc": float(average_precision_score(y_true, y_prob)) if y_true.sum() > 0 else None,
        f"{prefix}accuracy": float(accuracy_score(y_true, y_pred)),
        f"{prefix}precision": float(precision_score(y_true, y_pred, zero_division=0)),
        f"{prefix}recall": float(recall_score(y_true, y_pred, zero_division=0)),
        f"{prefix}f1": float(f1_score(y_true, y_pred, zero_division=0)),
        f"{prefix}mcc": float(matthews_corrcoef(y_true, y_pred)) if len(set(y_true)) > 1 else None,
        f"{prefix}brier_score": float(brier_score_loss(y_true, y_prob)),
        f"{prefix}sensitivity": float(tp / (tp + fn)) if (tp + fn) > 0 else 0.0,
        f"{prefix}specificity": float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0,
        f"{prefix}ppv": float(tp / (tp + fp)) if (tp + fp) > 0 else 0.0,
        f"{prefix}npv": float(tn / (tn + fn)) if (tn + fn) > 0 else 0.0,
    }

    if include_ci and len(y_true) >= 20:
        def _auc_fn(yt, yp, _):
            return roc_auc_score(yt, yp) if len(set(yt)) > 1 else np.nan
        def _mcc_fn(yt, _, yd):
            return matthews_corrcoef(yt, yd) if len(set(yt)) > 1 else np.nan

        auc_lo, auc_hi = bootstrap_ci(y_true, y_prob, y_pred, _auc_fn)
        if auc_lo is not None:
            metrics[f"{prefix}roc_auc_ci_low"] = auc_lo
            metrics[f"{prefix}roc_auc_ci_high"] = auc_hi

        mcc_lo, mcc_hi = bootstrap_ci(y_true, y_prob, y_pred, _mcc_fn)
        if mcc_lo is not None:
            metrics[f"{prefix}mcc_ci_low"] = mcc_lo
            metrics[f"{prefix}mcc_ci_high"] = mcc_hi

    return metrics


def load_extracted_features(feature_dir):
    print("🔄 Loading extracted features...")
    feature_dir = Path(feature_dir)

    # Fast path 1: use a single pre-combined file if available
    combined_path = feature_dir / "combined.parquet"
    if combined_path.exists():
        df = pd.read_parquet(combined_path)
        df.columns = [clean_name(c) for c in df.columns]
        if df.index.name is None:
            df.index.name = "variant_id"
        print(f"   Total variants with features: {len(df):,}")
        print(f"   Features: {len(df.columns):,}")
        return df

    # Fast path 2: use per-run pre-combined files if available
    run_combined = sorted(feature_dir.glob("output_test_run_*_combined.parquet"))
    if run_combined:
        print(f"   Found {len(run_combined):,} per-run combined files")
        dfs = []
        for f in run_combined:
            try:
                dfs.append(pd.read_parquet(f))
            except Exception as e:
                print(f"   ⚠ Error reading {f}: {e}")
        df = pd.concat(dfs, ignore_index=False)
        df.columns = [clean_name(c) for c in df.columns]
        if df.index.name is None:
            df.index.name = "variant_id"
        n_before = len(df)
        df = df[~df.index.duplicated(keep="first")]
        n_deduped = n_before - len(df)
        if n_deduped > 0:
            print(f"   ⚠ Deduplicated: removed {n_deduped:,} duplicate variant rows")
        print(f"   Total variants with features: {len(df):,}")
        print(f"   Features: {len(df.columns):,}")
        return df

    parquet_files = list(feature_dir.rglob("*.parquet"))
    print(f"   Found {len(parquet_files):,} parquet files")
    if not parquet_files:
        raise FileNotFoundError(f"No parquet files found in {feature_dir}")

    # Fallback: read individual files with pandas
    dfs = []
    for f in tqdm(parquet_files, desc="Loading parquets"):
        try:
            dfs.append(pd.read_parquet(f))
        except Exception as e:
            print(f"   ⚠ Error reading {f}: {e}")
    df = pd.concat(dfs, ignore_index=False)

    df.columns = [clean_name(c) for c in df.columns]
    if df.index.name is None:
        df.index.name = "variant_id"
    n_before = len(df)
    df = df[~df.index.duplicated(keep="first")]
    n_deduped = n_before - len(df)
    if n_deduped > 0:
        print(f"   ⚠ Deduplicated: removed {n_deduped:,} duplicate variant rows")
    print(f"   Total variants with features: {len(df):,}")
    print(f"   Features: {len(df.columns):,}")
    return df


def query_gnomad_af(variant_ids, vcf_path, outdir):
    print("\n🔄 Querying gnomAD AF via tabix...")
    Path(outdir).mkdir(exist_ok=True, parents=True)
    bed_entries = []
    vid_lookup = {}
    for vid in variant_ids:
        parts = vid.split("_")
        if len(parts) != 4:
            continue
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
        bed_entries.append((chrom, pos - 1, pos, vid))
        vid_lookup[(chrom, pos, ref, alt)] = vid

    bed_entries.sort(key=lambda x: (x[0], x[1]))
    bed_path = Path(outdir) / "query.bed"
    with open(bed_path, "w") as fh:
        for chrom, start, end, vid in bed_entries:
            fh.write(f"{chrom}\t{start}\t{end}\n")

    cmd = ["tabix", "-R", str(bed_path), str(vcf_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"tabix failed: {proc.stderr}")

    af_map = {}
    for line in proc.stdout.strip().split("\n"):
        cols = line.split("\t")
        if len(cols) < 8:
            continue
        chrom, pos, _, ref, alt = cols[0], int(cols[1]), cols[2], cols[3], cols[4]
        info = cols[7]
        alts = alt.split(",")
        for i, a in enumerate(alts):
            key = (chrom, pos, ref, a)
            if key in vid_lookup:
                vid = vid_lookup[key]
                m = re.search(r'AF=([^;]+)', info)
                if m:
                    af_str = m.group(1)
                    af_vals = af_str.split(",")
                    if i < len(af_vals):
                        try:
                            af_map[vid] = float(af_vals[i])
                        except ValueError:
                            af_map[vid] = 0.0
                    else:
                        af_map[vid] = 0.0
                else:
                    af_map[vid] = 0.0

    n_found = sum(1 for vid in variant_ids if vid in af_map)
    print(f"   Found in gnomAD: {n_found:,} / {len(variant_ids):,}")
    print(f"   Not found (AF=0): {len(variant_ids) - n_found:,}")
    return af_map


def load_pipeline(model_path):
    print("\n🔄 Loading frozen pipeline...")
    pipe = joblib.load(model_path)
    print(f"   Model: {type(pipe['model']).__name__}")
    print(f"   Calibrator: {type(pipe['calibrator']).__name__}")
    print(f"   Scaler: {type(pipe['scaler']).__name__}")
    print(f"   Features: {len(pipe['features']):,}")
    print(f"   Threshold: {pipe['threshold']:.4f}")
    return pipe


def load_labels(benchmark_path):
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
    print(f"   Pathogenic: {bench['label'].sum():,}")
    print(f"   Benign: {(bench['label'] == 0).sum():,}")
    return bench


def load_train_genes():
    if TRAIN_GENES_PATH.exists():
        with open(TRAIN_GENES_PATH) as f:
            return set(json.load(f))
    train = pd.read_csv(ROOT / "train.csv", usecols=["GeneSymbol"])
    cal = pd.read_csv(ROOT / "cal.csv", usecols=["GeneSymbol"])
    test = pd.read_csv(ROOT / "test.csv", usecols=["GeneSymbol"])
    genes = set(train["GeneSymbol"].dropna().unique()) | set(cal["GeneSymbol"].dropna().unique()) | set(test["GeneSymbol"].dropna().unique())
    with open(TRAIN_GENES_PATH, "w") as f:
        json.dump(list(genes), f)
    return genes


def score_benchmark(feature_df, pipe, labels_df, train_genes, outdir, has_af=True, vep_path=None):
    """Generic classification scoring with subgroup analysis."""
    print("\n🔄 Scoring...")
    variant_ids = list(feature_df.index)

    if has_af:
        af_map = query_gnomad_af(variant_ids, GNOMAD_VCF, outdir)
        feature_df["gnomAD_AF"] = feature_df.index.map(lambda v: af_map.get(v, 0.0))
        feature_df["gnomAD_AF"] = feature_df["gnomAD_AF"].astype(float)

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
    if has_af:
        pred_df["gnomAD_AF"] = feature_df["gnomAD_AF"].values

    label_map = dict(zip(labels_df["variant_id"], labels_df["label"]))
    gene_map = dict(zip(labels_df["variant_id"], labels_df["gene"]))
    pred_df["true_label"] = pred_df["variant_id"].map(label_map)
    pred_df["gene"] = pred_df["variant_id"].map(gene_map)

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

    if has_af:
        pred_df["AF_zero"] = (pred_df["gnomAD_AF"] == 0).astype(int)
        pred_df["ultra_rare"] = (pred_df["gnomAD_AF"] < 0.0001).astype(int)
    pred_df["unseen_gene"] = (~pred_df["gene"].isin(train_genes)).astype(int)

    y_true = pred_df["true_label"].values
    y_prob = pred_df["prob_calibrated"].values
    y_pred = pred_df["prediction"].values

    print(f"\n📊 SCORED VARIANTS: {len(pred_df):,}")
    print(f"   Pathogenic: {y_true.sum():,}")
    print(f"   Benign:     {(y_true == 0).sum():,}")

    all_metrics = compute_metrics(y_true, y_prob, y_pred, prefix="overall_", include_ci=True)

    subgroup_results = []
    if has_af:
        subgroups = {
            "AF_zero": pred_df["AF_zero"] == 1,
            "AF_present": pred_df["AF_zero"] == 0,
            "ultra_rare": pred_df["ultra_rare"] == 1,
            "common": pred_df["ultra_rare"] == 0,
            "unseen_genes": pred_df["unseen_gene"] == 1,
            "seen_genes": pred_df["unseen_gene"] == 0,
        }
    else:
        subgroups = {
            "unseen_genes": pred_df["unseen_gene"] == 1,
            "seen_genes": pred_df["unseen_gene"] == 0,
        }

    for name, mask in subgroups.items():
        sub = pred_df[mask]
        if len(sub) == 0:
            continue
        sm = compute_metrics(
            sub["true_label"].values,
            sub["prob_calibrated"].values,
            sub["prediction"].values,
            prefix="", include_ci=True
        )
        sm["subgroup"] = name
        sm["n_variants"] = len(sub)
        subgroup_results.append(sm)

    af_nonzero = pred_df.loc[pred_df["gnomAD_AF"] > 0, "gnomAD_AF"] if has_af else pd.Series([])
    novel_gene_set = set(pred_df.loc[pred_df["unseen_gene"] == 1, "gene"])
    seen_gene_set = set(pred_df.loc[pred_df["unseen_gene"] == 0, "gene"])

    characterization = {
        "n_total": int(len(pred_df)),
        "n_pathogenic": int(y_true.sum()),
        "n_benign": int((y_true == 0).sum()),
        "n_genes": int(pred_df["gene"].nunique()),
        "novel_genes_percent": float(pred_df["unseen_gene"].mean() * 100),
        "n_novel_genes": len(novel_gene_set),
        "n_seen_genes": len(seen_gene_set),
    }
    if has_af:
        characterization["AF_zero_percent"] = float((pred_df["gnomAD_AF"] == 0).mean() * 100)
        characterization["median_AF"] = float(af_nonzero.median()) if len(af_nonzero) > 0 else 0.0
        characterization["mean_AF"] = float(pred_df["gnomAD_AF"].mean())
        characterization["ultra_rare_percent"] = float(pred_df["ultra_rare"].mean() * 100)

    return pred_df, all_metrics, subgroup_results, characterization


def save_results(pred_df, metrics, subgroup_results, characterization, outdir, benchmark_name):
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

    # Print key metrics with CIs if available
    print()
    print("📈 KEY METRICS")
    n = metrics.get("overall_n_variants", 0)
    print(f"   N = {n:,}")
    auc = metrics.get("overall_roc_auc")
    if auc is not None:
        auc_lo = metrics.get("overall_roc_auc_ci_low")
        auc_hi = metrics.get("overall_roc_auc_ci_high")
        if auc_lo is not None and auc_hi is not None:
            print(f"   ROC-AUC: {auc:.3f} (95% CI: {auc_lo:.3f}–{auc_hi:.3f})")
        else:
            print(f"   ROC-AUC: {auc:.3f}")
    mcc = metrics.get("overall_mcc")
    if mcc is not None:
        mcc_lo = metrics.get("overall_mcc_ci_low")
        mcc_hi = metrics.get("overall_mcc_ci_high")
        if mcc_lo is not None and mcc_hi is not None:
            print(f"   MCC:     {mcc:.3f} (95% CI: {mcc_lo:.3f}–{mcc_hi:.3f})")
        else:
            print(f"   MCC:     {mcc:.3f}")

    if subgroup_results:
        sub_df = pd.DataFrame(subgroup_results)
        sub_path = outdir / f"{benchmark_name}_subgroup_metrics.csv"
        sub_df.to_csv(sub_path, index=False)
        print(f"💾 Subgroup metrics saved: {sub_path}")

        print("\n" + "="*70)
        print("📊 SUBGROUP METRICS")
        print("="*70)
        display_cols = ["subgroup", "n_variants", "roc_auc", "pr_auc", "accuracy", "f1", "mcc", "sensitivity", "specificity"]
        for col in display_cols:
            if col not in sub_df.columns:
                display_cols = [c for c in display_cols if c in sub_df.columns]
                break
        print(sub_df[display_cols].to_string(index=False))

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
    print("📊 OVERALL METRICS")
    print("="*70)
    for k in ["overall_roc_auc", "overall_pr_auc", "overall_accuracy", "overall_f1", "overall_mcc", "overall_brier_score", "overall_sensitivity", "overall_specificity"]:
        if k in metrics:
            val = metrics[k]
            val_str = f"{val:.4f}" if val is not None else "N/A"
            print(f"  {k.replace('overall_', '').replace('_', ' ').title():15s}: {val_str}")
    print("="*70)
