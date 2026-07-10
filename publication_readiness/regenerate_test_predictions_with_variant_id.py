#!/usr/bin/env python3
"""
Regenerate final_model_output_regularized/test_predictions.csv with variant_id
and recompute the test-set bootstrap CI using 2,000 stratified resamples.

This script loads the saved final pipeline and re-scores the test set in a
deterministic (sorted parquet) order so that predictions can be merged with
feature matrices and downstream analyses by variant_id.
"""

import os
import re
import glob
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score,
    precision_score, recall_score, matthews_corrcoef,
    brier_score_loss, confusion_matrix
)

OUTDIR = "final_model_output_regularized"


def clean_name(s):
    s = str(s)
    s = s.replace(":", "_").replace("-", "_").replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9_]", "", s)
    s = re.sub(r"_+", "_", s)
    return s


def load_dataset_sorted(path, selected_features=None):
    """Load parquet parts in deterministic sorted order."""
    files = sorted(glob.glob(f"{path}/*.parquet"))
    if not files:
        raise FileNotFoundError(f"No parquet files found in {path}")
    dfs = []
    for i, f in enumerate(files):
        dfs.append(pd.read_parquet(f))
        if i % 500 == 0:
            print(f"{path}: {i}/{len(files)}")
    df = pd.concat(dfs, ignore_index=True)
    df.columns = [clean_name(c) for c in df.columns]
    if selected_features is not None:
        keep = ["variant_id", "label", "gene"] + selected_features
        keep = [c for c in keep if c in df.columns]
        df = df[keep]
    return df


def stratified_bootstrap_auc(y_true, y_prob, n_resamples=2000, random_state=42):
    """Bootstrap AUC preserving class prevalence."""
    rng = np.random.RandomState(random_state)
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    scores = []
    for _ in range(n_resamples):
        b_pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        b_neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        idx = np.concatenate([b_pos, b_neg])
        scores.append(roc_auc_score(y_true[idx], y_prob[idx]))
    return np.percentile(scores, 2.5), np.percentile(scores, 97.5)


def main():
    print("Loading saved pipeline...")
    pipeline = joblib.load(f"{OUTDIR}/final_pipeline.pkl")
    model = pipeline["model"]
    calibrator = pipeline["calibrator"]
    imputer = pipeline["imputer"]
    scaler = pipeline["scaler"]
    selected_features = pipeline["features"]
    threshold = pipeline["threshold"]

    print("Loading test set in deterministic sorted order...")
    test_df = load_dataset_sorted("final_dataset_parts_test", selected_features)

    if "variant_id" not in test_df.columns:
        raise ValueError("variant_id missing from test feature matrix")

    X_test = test_df[selected_features]
    y_test = test_df["label"]
    variant_ids = test_df["variant_id"]

    print("Scoring test set...")
    X_test_imp = X_test.fillna(imputer)
    X_test_f = pd.DataFrame(
        scaler.transform(X_test_imp),
        columns=X_test_imp.columns,
        index=X_test_imp.index,
    )
    probs_test_raw = model.predict_proba(X_test_f)[:, 1]
    probs_test = calibrator.transform(probs_test_raw)
    preds = (probs_test > threshold).astype(int)

    roc_auc = roc_auc_score(y_test, probs_test)
    auprc = average_precision_score(y_test, probs_test)
    f1 = f1_score(y_test, preds)
    precision = precision_score(y_test, preds)
    recall = recall_score(y_test, preds)
    mcc = matthews_corrcoef(y_test, preds)
    brier = brier_score_loss(y_test, probs_test)

    print(f"ROC-AUC: {roc_auc:.6f}")
    print(f"AUPRC:   {auprc:.6f}")

    print("Computing 2,000 stratified bootstrap CI...")
    ci_low, ci_high = stratified_bootstrap_auc(y_test.values, probs_test)
    print(f"95% CI: {ci_low:.6f} - {ci_high:.6f}")

    # Save predictions with variant_id
    pred_df = pd.DataFrame({
        "variant_id": variant_ids.values,
        "y_true": y_test.values,
        "raw_prob": probs_test_raw,
        "cal_prob": probs_test,
    })
    pred_path = f"{OUTDIR}/test_predictions.csv"
    pred_df.to_csv(pred_path, index=False)
    print(f"Saved predictions with variant_id to {pred_path}")

    # Update final_metrics.csv
    metrics_df = pd.DataFrame({
        "metric": [
            "AUC", "AUPRC", "F1", "Precision", "Recall", "MCC",
            "Brier", "AUC_CI_low", "AUC_CI_high"
        ],
        "value": [
            roc_auc, auprc, f1, precision, recall, mcc,
            brier, ci_low, ci_high
        ]
    })
    # Preserve existing ECE and Baseline_AUC if present
    old_metrics_path = f"{OUTDIR}/final_metrics.csv"
    if os.path.exists(old_metrics_path):
        old = pd.read_csv(old_metrics_path)
        extra = old[old["metric"].isin(["ECE", "Baseline_AUC"])]
        if not extra.empty:
            metrics_df = pd.concat([metrics_df, extra], ignore_index=True)
    metrics_path = f"{OUTDIR}/final_metrics.csv"
    metrics_df.to_csv(metrics_path, index=False)
    print(f"Updated metrics saved to {metrics_path}")

    # Save confusion matrix
    cm = confusion_matrix(y_test, preds)
    np.savetxt(f"{OUTDIR}/confusion_matrix.txt", cm, fmt="%d")


if __name__ == "__main__":
    main()
