#!/usr/bin/env python3
"""
Variant-type analysis for AGVarPred.

Uses only existing predictions and VEP annotations.
Does NOT retrain models or re-extract features.

Outputs:
  - internal_test_variant_type_metrics.csv
  - external_variant_type_metrics.csv
  - variant_type_coverage.csv
  - non_missense_comparison.csv
  - variant_type_auc.pdf/png
  - variant_type_coverage.pdf/png
  - non_missense_roc.pdf/png
"""

import os
import json
import glob
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score,
    matthews_corrcoef, accuracy_score, precision_score, recall_score,
    roc_curve
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(OUT_DIR)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BENCHMARKS = {
    "internal_test": {
        "predictions": "final_model_output_regularized/test_predictions.csv",
        "features": "final_dataset_parts_test/part_*.parquet",
        "vep_parquet": None,
        "label_col": "y_true",
        "score_col": "cal_prob",
    },
    "humsavar": {
        "predictions": "external_validation/results/humsavar/regularized/humsavar_predictions.csv",
        "vep_parquet": "external_validation/humsavar_vep.parquet",
        "label_col": "true_label",
        "score_col": "prob_calibrated",
    },
    "mave_independent": {
        "predictions": "external_validation/results/mave_independent/regularized/mave_independent_predictions.csv",
        "vep_parquet": "external_validation/mave_independent_vep.parquet",
        "label_col": "true_label",
        "score_col": "prob_calibrated",
    },
    "gnomad_benign": {
        "predictions": "external_validation/results/gnomad_benign/regularized/gnomad_benign_predictions.csv",
        "vep_parquet": "external_validation/gnomad_benign_vep.parquet",
        "label_col": "true_label",
        "score_col": "prob_calibrated",
    },
    # VIP is excluded from the independent external-benchmark comparison because
    # all pathogenic labels are ClinVar-derived; it is reported separately as a
    # ClinVar-held-out-gene consistency analysis.
    # "vip": { ... },
    "grimm2015": {
        "predictions": "external_validation/results/grimm2015/regularized/grimm2015_predictions.csv",
        "vep_parquet": "external_validation/vep_preprocessed/grimm2015_vep.parquet",
        "label_col": "true_label",
        "score_col": "prob_calibrated",
    },
    "dvd": {
        "predictions": "external_validation/results/dvd/regularized/dvd_predictions.csv",
        "vep_parquet": "external_validation/vep_preprocessed/dvd_vep.parquet",
        "label_col": "true_label",
        "score_col": "prob_calibrated",
    },
}

COMPETITOR_DATA = {
    "humsavar": {
        "variants": "external_validation/benchmarks_comparative/datasets/humsavar/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/humsavar/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/humsavar/revel_scores.csv",
        "json_key_has_chr": False,
    },
    "mave_independent": {
        "variants": "external_validation/benchmarks_comparative/datasets/mave_independent/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/mave_independent/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/mave_independent/revel_scores.csv",
        "json_key_has_chr": False,
    },
    "gnomad_benign": {
        "variants": "external_validation/benchmarks_comparative/datasets/gnomad_benign/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/gnomad_benign/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/gnomad_benign/revel_scores.csv",
        "json_key_has_chr": False,
    },
    # VIP excluded from main comparison (ClinVar-derived labels).
    # "vip": { ... },
    "grimm2015": {
        "variants": "external_validation/benchmarks_comparative/datasets/grimm2015/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/grimm2015/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/grimm2015/revel_scores.csv",
        "json_key_has_chr": True,
    },
    "dvd": {
        "variants": "external_validation/benchmarks_comparative/datasets/dvd/variants_scored.tsv",
        "json": "external_validation/benchmarks_comparative/datasets/dvd/vep_cadd_alphamissense.json",
        "revel": "external_validation/benchmarks_comparative/datasets/dvd/revel_scores.csv",
        "json_key_has_chr": False,
    },
}

BOOTSTRAP_N = 1000
MIN_N = 10

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def classify_variant_type(df, consequence_col="vep_Consequence"):
    """
    Assign each variant to a single high-level type.
    Priority: stop_gained > frameshift > splice > missense > synonymous > UTR > intron > other.
    """
    ctype = pd.Series("other", index=df.index)

    cons = df[consequence_col].fillna("").astype(str).str.lower()

    # Use explicit flags when available
    if "vep_is_stop_gained" in df.columns:
        mask = df["vep_is_stop_gained"].fillna(0).astype(int) == 1
        ctype.loc[mask] = "stop_gained"

    if "vep_is_frameshift" in df.columns:
        mask = (df["vep_is_frameshift"].fillna(0).astype(int) == 1) & (ctype == "other")
        ctype.loc[mask] = "frameshift"

    if "vep_is_splice" in df.columns:
        mask = (df["vep_is_splice"].fillna(0).astype(int) == 1) & (ctype == "other")
        ctype.loc[mask] = "splice"

    if "vep_is_missense" in df.columns:
        mask = (df["vep_is_missense"].fillna(0).astype(int) == 1) & (ctype == "other")
        ctype.loc[mask] = "missense"

    if "vep_is_synonymous" in df.columns:
        mask = (df["vep_is_synonymous"].fillna(0).astype(int) == 1) & (ctype == "other")
        ctype.loc[mask] = "synonymous"

    # Fallback / refinement from consequence string
    if consequence_col in df.columns:
        def set_if(cond, val):
            mask = cond & (ctype == "other")
            ctype.loc[mask] = val

        set_if(cons.str.contains("stop_gained"), "stop_gained")
        set_if(cons.str.contains("frameshift"), "frameshift")
        set_if(cons.str.contains("splice_donor|splice_acceptor|splice_region"), "splice")
        set_if(cons.str.contains("missense"), "missense")
        set_if(cons.str.contains("synonymous"), "synonymous")
        set_if(cons.str.contains("utr_variant"), "UTR")
        set_if(cons.str.contains("intron_variant"), "intron")
        # Non-coding exon distinct from intron
        set_if(cons.str.contains("non_coding_transcript_exon"), "non_coding_exon")

    return ctype


def merge_small_classes(df, min_n=MIN_N, target="other"):
    """
    Re-classify variants in classes with fewer than min_n members into target class.
    Ensures downstream per-class tables sum to the full benchmark size.
    """
    counts = df["variant_type"].value_counts()
    small = counts[counts < min_n].index.tolist()
    if target in small:
        # ensure target itself is not reclassified; if target is small, keep it
        small.remove(target)
    if small:
        df = df.copy()
        df.loc[df["variant_type"].isin(small), "variant_type"] = target
    return df


def compute_metrics(y_true, y_score, y_pred=None, threshold=0.5):
    """Compute classification metrics. Returns dict."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)

    if y_pred is None:
        y_pred = (y_score >= threshold).astype(int)

    res = {"N": len(y_true)}
    res["N_positive"] = int(np.sum(y_true == 1))
    res["N_negative"] = int(np.sum(y_true == 0))

    # AUC only if both classes present
    if len(np.unique(y_true)) > 1 and len(y_score) >= MIN_N:
        try:
            res["AUC"] = roc_auc_score(y_true, y_score)
        except Exception:
            res["AUC"] = np.nan
        try:
            res["AUPRC"] = average_precision_score(y_true, y_score)
        except Exception:
            res["AUPRC"] = np.nan
    else:
        res["AUC"] = np.nan
        res["AUPRC"] = np.nan

    # FPR for all-negative sets (e.g. gnomAD benign)
    if res["N_negative"] > 0:
        res["FPR"] = float(np.sum((y_pred == 1) & (y_true == 0)) / res["N_negative"])
    else:
        res["FPR"] = np.nan

    # Metrics requiring both classes
    if len(np.unique(y_true)) > 1:
        res["Accuracy"] = accuracy_score(y_true, y_pred)
        res["Precision"] = precision_score(y_true, y_pred, zero_division=0)
        res["Recall"] = recall_score(y_true, y_pred, zero_division=0)
        res["F1"] = f1_score(y_true, y_pred, zero_division=0)
        res["MCC"] = matthews_corrcoef(y_true, y_pred)
        res["Specificity"] = float(np.sum((y_pred == 0) & (y_true == 0)) / res["N_negative"]) if res["N_negative"] > 0 else np.nan
    else:
        res["Accuracy"] = accuracy_score(y_true, y_pred)
        res["Precision"] = np.nan
        res["Recall"] = np.nan
        res["F1"] = np.nan
        res["MCC"] = np.nan
        res["Specificity"] = np.nan

    return res


def bootstrap_metrics(y_true, y_score, threshold=0.5, n_boot=BOOTSTRAP_N, seed=42):
    """Return 95% CI for AUC and AUPRC via stratified bootstrap."""
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)

    aucs, auprcs = [], []
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]

    for _ in range(n_boot):
        if len(pos_idx) > 0 and len(neg_idx) > 0:
            b_pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
            b_neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
            idx = np.concatenate([b_pos, b_neg])
        else:
            idx = rng.choice(len(y_true), size=len(y_true), replace=True)
        yt = y_true[idx]
        ys = y_score[idx]
        if len(np.unique(yt)) > 1:
            aucs.append(roc_auc_score(yt, ys))
            auprcs.append(average_precision_score(yt, ys))

    def ci(x):
        return (np.nanpercentile(x, 2.5), np.nanpercentile(x, 97.5)) if x else (np.nan, np.nan)

    return {"AUC_CI_low": ci(aucs)[0], "AUC_CI_high": ci(aucs)[1],
            "AUPRC_CI_low": ci(auprcs)[0], "AUPRC_CI_high": ci(auprcs)[1]}


# ---------------------------------------------------------------------------
# Internal test analysis
# ---------------------------------------------------------------------------


def analyze_internal_test():
    print("\n[Internal test]")
    cfg = BENCHMARKS["internal_test"]

    pred = pd.read_csv(os.path.join(ROOT_DIR, cfg["predictions"]))

    # Load feature matrix in sorted order and align by variant_id
    part_files = sorted(glob.glob(os.path.join(ROOT_DIR, cfg["features"])))
    feats = []
    for f in part_files:
        feats.append(pd.read_parquet(f))
    feats = pd.concat(feats, ignore_index=True)

    if "variant_id" not in pred.columns:
        raise ValueError("Internal test predictions must contain a 'variant_id' column")
    if "variant_id" not in feats.columns:
        raise ValueError("Internal test feature matrix must contain a 'variant_id' column")

    df = pred.merge(feats, on="variant_id", how="inner")

    if len(df) != len(pred):
        raise ValueError(
            f"Merge mismatch: predictions={len(pred)}, features={len(feats)}, merged={len(df)}"
        )

    # Ensure the label used for evaluation is the one from the predictions file
    if "label" in df.columns:
        df = df.drop(columns=["label"])

    # Reconstruct VEP consequence string from one-hot columns (internal features lack vep_Consequence)
    if "vep_Consequence" not in df.columns:
        cons_cols = [c for c in df.columns if c.startswith("vep_Consequence_")]
        def reconstruct(row):
            active = [c.replace("vep_Consequence_", "") for c in cons_cols if row.get(c, 0) == 1]
            return "&".join(active) if active else "unknown"
        df["vep_Consequence"] = df.apply(reconstruct, axis=1)

    df["variant_type"] = classify_variant_type(df)
    df = merge_small_classes(df)

    records = []
    for vtype in df["variant_type"].unique():
        sub = df[df["variant_type"] == vtype]
        if len(sub) < MIN_N:
            continue
        y_true = sub[cfg["label_col"]].values
        y_score = sub[cfg["score_col"]].values
        metrics = compute_metrics(y_true, y_score)
        ci = bootstrap_metrics(y_true, y_score)
        metrics.update(ci)
        metrics["benchmark"] = "internal_test"
        metrics["variant_type"] = vtype
        records.append(metrics)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# External validation analysis
# ---------------------------------------------------------------------------


def normalize_variant_id(vid, has_chr):
    """Normalize to chrX_Y_Z_A format."""
    vid = str(vid).strip()
    if has_chr and not vid.startswith("chr"):
        return "chr" + vid
    if not has_chr and vid.startswith("chr"):
        return vid[3:]
    return vid


def load_competitor_scores(bench, add_vep=True):
    info = COMPETITOR_DATA[bench]
    variants = pd.read_csv(os.path.join(ROOT_DIR, info["variants"]), sep="\t")
    variants = variants.rename(columns={"label": "true_label"})

    # CADD and AlphaMissense from JSON
    with open(os.path.join(ROOT_DIR, info["json"])) as f:
        cadd_alpha = json.load(f)

    rows = []
    for _, row in variants.iterrows():
        vid = row["variant_id"]
        key_no_chr = vid.replace("chr", "") if vid.startswith("chr") else vid
        key_with_chr = vid if vid.startswith("chr") else "chr" + vid
        key = key_with_chr if info["json_key_has_chr"] else key_no_chr
        ca = cadd_alpha.get(key, {})
        rows.append({
            "variant_id": vid,
            "cadd_phred": ca.get("cadd_phred"),
            "alphamissense_score": ca.get("alphamissense_score"),
        })
    comp = pd.DataFrame(rows)

    # REVEL
    revel = pd.read_csv(os.path.join(ROOT_DIR, info["revel"]))
    revel = revel.rename(columns={"revel_score": "revel_score"})
    comp = comp.merge(revel, on="variant_id", how="left")
    comp = comp.merge(variants[["variant_id", "true_label"]], on="variant_id", how="left")

    # Attach VEP consequence for variant-type classification
    if add_vep:
        cfg = BENCHMARKS[bench]
        vep = pd.read_parquet(os.path.join(ROOT_DIR, cfg["vep_parquet"]))
        if vep.index.name == "variant_id":
            vep = vep.reset_index()
        # Keep only VEP columns needed for classification
        keep = ["variant_id"] + [c for c in vep.columns if c in [
            "vep_Consequence", "vep_is_missense", "vep_is_synonymous",
            "vep_is_stop_gained", "vep_is_frameshift", "vep_is_splice",
            "vep_LoF", "vep_is_LoF_HC"]]
        comp = comp.merge(vep[keep], on="variant_id", how="left")

    return comp


def analyze_external(bench):
    print(f"\n[{bench}]")
    cfg = BENCHMARKS[bench]

    pred = pd.read_csv(os.path.join(ROOT_DIR, cfg["predictions"]))
    vep = pd.read_parquet(os.path.join(ROOT_DIR, cfg["vep_parquet"]))

    # Ensure variant_id is column in vep (some files use it as index)
    if vep.index.name == "variant_id":
        vep = vep.reset_index()

    df = pred.merge(vep, on="variant_id", how="left")
    df["variant_type"] = classify_variant_type(df)
    df = merge_small_classes(df)

    records = []
    for vtype in df["variant_type"].unique():
        sub = df[df["variant_type"] == vtype]
        if len(sub) < MIN_N:
            continue
        y_true = sub[cfg["label_col"]].values
        y_score = sub[cfg["score_col"]].values
        metrics = compute_metrics(y_true, y_score)
        if len(np.unique(y_true)) > 1:
            ci = bootstrap_metrics(y_true, y_score)
            metrics.update(ci)
        else:
            metrics.update({"AUC_CI_low": np.nan, "AUC_CI_high": np.nan,
                            "AUPRC_CI_low": np.nan, "AUPRC_CI_high": np.nan})
        metrics["benchmark"] = bench
        metrics["variant_type"] = vtype
        records.append(metrics)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Competitor coverage by variant type
# ---------------------------------------------------------------------------


def compute_coverage_by_type():
    print("\n[Competitor coverage by variant type]")
    records = []

    for bench in COMPETITOR_DATA:
        cfg = BENCHMARKS[bench]
        pred = pd.read_csv(os.path.join(ROOT_DIR, cfg["predictions"]))
        vep = pd.read_parquet(os.path.join(ROOT_DIR, cfg["vep_parquet"]))
        if vep.index.name == "variant_id":
            vep = vep.reset_index()
        df = pred.merge(vep, on="variant_id", how="left")
        df["variant_type"] = classify_variant_type(df)
        df = merge_small_classes(df)

        # Load competitor scores and merge onto AGVarPred-scored variants
        comp = load_competitor_scores(bench, add_vep=False)
        comp = comp.drop(columns=["true_label"], errors="ignore")
        df = df.merge(comp, on="variant_id", how="left")

        for vtype in df["variant_type"].unique():
            sub = df[df["variant_type"] == vtype]
            n_total = len(sub)

            # AGVarPred coverage = 100% by design
            records.append({
                "benchmark": bench,
                "variant_type": vtype,
                "method": "AGVarPred",
                "N_scored": n_total,
                "N_total": n_total,
                "coverage_pct": 100.0,
            })

            for method, col in [("CADD", "cadd_phred"),
                                ("AlphaMissense", "alphamissense_score"),
                                ("REVEL", "revel_score")]:
                n_scored = int(sub[col].notna().sum())
                records.append({
                    "benchmark": bench,
                    "variant_type": vtype,
                    "method": method,
                    "N_scored": n_scored,
                    "N_total": n_total,
                    "coverage_pct": 100.0 * n_scored / n_total if n_total > 0 else np.nan,
                })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Non-missense-only comparison
# ---------------------------------------------------------------------------


def compute_non_missense_comparison():
    print("\n[Non-missense comparison]")
    records = []

    for bench in COMPETITOR_DATA:
        cfg = BENCHMARKS[bench]
        pred = pd.read_csv(os.path.join(ROOT_DIR, cfg["predictions"]))
        pred = pred.rename(columns={cfg["score_col"]: "agvarpred_score", cfg["label_col"]: "true_label"})

        vep = pd.read_parquet(os.path.join(ROOT_DIR, cfg["vep_parquet"]))
        if vep.index.name == "variant_id":
            vep = vep.reset_index()

        comp = load_competitor_scores(bench, add_vep=False)
        comp = comp.drop(columns=["true_label"], errors="ignore")
        df = pred.merge(vep, on="variant_id", how="left").merge(comp, on="variant_id", how="left")
        df["variant_type"] = classify_variant_type(df)
        df = merge_small_classes(df)
        sub = df[df["variant_type"] != "missense"]
        if len(sub) < MIN_N:
            continue

        y_true = sub["true_label"].values
        for method, col in [
            ("AGVarPred", "agvarpred_score"),
            ("CADD", "cadd_phred"),
            ("AlphaMissense", "alphamissense_score"),
            ("REVEL", "revel_score"),
        ]:
            scores = sub[col].values
            valid = ~np.isnan(scores)
            if valid.sum() < MIN_N:
                continue
            yt = y_true[valid]
            ys = scores[valid]
            if len(np.unique(yt)) > 1:
                auc = roc_auc_score(yt, ys)
                auprc = average_precision_score(yt, ys)
                ci = bootstrap_metrics(yt, ys)
            else:
                auc = np.nan
                auprc = np.nan
                ci = {"AUC_CI_low": np.nan, "AUC_CI_high": np.nan,
                      "AUPRC_CI_low": np.nan, "AUPRC_CI_high": np.nan}
            records.append({
                "benchmark": bench,
                "method": method,
                "N": int(valid.sum()),
                "N_positive": int(np.sum(yt == 1)),
                "AUC": auc,
                "AUPRC": auprc,
                **ci,
            })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def plot_variant_type_auc(metrics_df, out_path):
    """Bar plot of AUC by variant type across benchmarks."""
    df = metrics_df.copy()
    df = df[df["N"] >= MIN_N]
    # Use AUC where available, FPR for all-negative gnomad benign
    df["metric_value"] = df["AUC"].fillna(-df["FPR"])
    df["metric_label"] = df.apply(lambda r: "FPR" if pd.isna(r["AUC"]) else "AUC", axis=1)

    benchmarks = df["benchmark"].unique()
    types = sorted(df["variant_type"].unique())

    fig, ax = plt.subplots(figsize=(16, 7))
    x = np.arange(len(types))
    width = 0.12
    colors = plt.cm.tab10(np.linspace(0, 1, len(benchmarks)))

    for i, bench in enumerate(benchmarks):
        sub = df[df["benchmark"] == bench].set_index("variant_type").reindex(types)
        vals = sub["metric_value"].fillna(0).values
        ax.bar(x + (i - len(benchmarks)/2) * width, vals, width, label=bench, color=colors[i])

    ax.set_xticks(x)
    ax.set_xticklabels(types, rotation=45, ha="right")
    ax.set_ylabel("AUC / -FPR")
    ax.set_title("AGVarPred performance by variant type")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False)
    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.5)
    plt.tight_layout(rect=[0, 0, 0.82, 1])
    fig.savefig(out_path, dpi=300)
    fig.savefig(out_path.replace(".pdf", ".png"), dpi=300)
    plt.close(fig)


def plot_coverage(coverage_df, out_path):
    """Coverage by variant type for AGVarPred vs competitors."""
    # Summarize across benchmarks: mean coverage per method per variant type
    summary = coverage_df.groupby(["variant_type", "method"])["coverage_pct"].mean().reset_index()
    summary = summary.pivot(index="variant_type", columns="method", values="coverage_pct").fillna(0)

    methods = ["AGVarPred", "CADD", "AlphaMissense", "REVEL"]
    summary = summary[[m for m in methods if m in summary.columns]]

    fig, ax = plt.subplots(figsize=(12, 7))
    summary.plot(kind="bar", ax=ax, color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"])
    ax.set_ylabel("Mean coverage (%)")
    ax.set_title("Variant coverage by type")
    ax.set_ylim(0, 105)
    ax.legend(title="Method", loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout(rect=[0, 0, 0.78, 1])
    fig.savefig(out_path, dpi=300)
    fig.savefig(out_path.replace(".pdf", ".png"), dpi=300)
    plt.close(fig)


def plot_non_missense_roc(nonmiss_df, out_path):
    """ROC curves for non-missense variants by benchmark and method."""
    # This requires per-variant scores; nonmiss_df only has aggregate AUC.
    # Skip detailed ROC here; will be generated if per-variant data is saved.
    # Instead, plot AUPRC/AUC summary.
    df = nonmiss_df.copy()
    df = df[df["N"] >= MIN_N]
    benchmarks = df["benchmark"].unique()
    methods = ["AGVarPred", "CADD", "AlphaMissense", "REVEL"]

    fig, ax = plt.subplots(figsize=(14, 7))
    x = np.arange(len(benchmarks))
    width = 0.2
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    for i, method in enumerate(methods):
        sub = df[df["method"] == method].set_index("benchmark").reindex(benchmarks)
        ax.bar(x + (i - 1.5) * width, sub["AUC"].fillna(0).values, width, label=method, color=colors[i])

    ax.set_xticks(x)
    ax.set_xticklabels(benchmarks, rotation=45, ha="right")
    ax.set_ylabel("AUC")
    ax.set_title("Non-missense variant performance")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.5)
    plt.tight_layout(rect=[0, 0, 0.78, 1])
    fig.savefig(out_path, dpi=300)
    fig.savefig(out_path.replace(".pdf", ".png"), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Internal test
    internal_df = analyze_internal_test()
    internal_df.to_csv(os.path.join(OUT_DIR, "internal_test_variant_type_metrics.csv"), index=False)
    print("Saved internal_test_variant_type_metrics.csv")

    # 2. External validation
    external_dfs = []
    for bench in ["humsavar", "mave_independent", "gnomad_benign", "grimm2015", "dvd"]:
        external_dfs.append(analyze_external(bench))
    external_df = pd.concat(external_dfs, ignore_index=True)
    external_df.to_csv(os.path.join(OUT_DIR, "external_variant_type_metrics.csv"), index=False)
    print("Saved external_variant_type_metrics.csv")

    # 3. Coverage
    coverage_df = compute_coverage_by_type()
    coverage_df.to_csv(os.path.join(OUT_DIR, "variant_type_coverage.csv"), index=False)
    print("Saved variant_type_coverage.csv")

    # 4. Non-missense comparison
    nonmiss_df = compute_non_missense_comparison()
    nonmiss_df.to_csv(os.path.join(OUT_DIR, "non_missense_comparison.csv"), index=False)
    print("Saved non_missense_comparison.csv")

    # 5. Figures
    all_metrics = pd.concat([internal_df, external_df], ignore_index=True)
    plot_variant_type_auc(all_metrics, os.path.join(OUT_DIR, "variant_type_auc.pdf"))
    print("Saved variant_type_auc.pdf/png")

    plot_coverage(coverage_df, os.path.join(OUT_DIR, "variant_type_coverage.pdf"))
    print("Saved variant_type_coverage.pdf/png")

    plot_non_missense_roc(nonmiss_df, os.path.join(OUT_DIR, "non_missense_auc_summary.pdf"))
    print("Saved non_missense_auc_summary.pdf/png")

    print("\nVariant-type analysis complete.")


if __name__ == "__main__":
    main()
