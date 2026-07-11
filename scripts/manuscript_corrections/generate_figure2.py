#!/usr/bin/env python3
"""
Generate a single 3-panel Figure 2 (ROC, PR, calibration) with 95% stratified
bootstrap confidence intervals for the AGVarPred internal test set.

Reads: final_model_output_regularized/test_predictions.csv
Writes: final_model_output_regularized/figure2_combined.pdf
        final_model_output_regularized/figure2_combined.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve, precision_recall_curve, roc_auc_score,
    average_precision_score, brier_score_loss
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRED_PATH = os.path.join(ROOT, "final_model_output_regularized", "test_predictions.csv")
OUTDIR = os.path.join(ROOT, "final_model_output_regularized")

N_BOOTSTRAPS = 2000
RANDOM_STATE = 42


def stratified_bootstrap_curves(y_true, y_score, curve_func, n_bootstraps=2000, random_state=42):
    """Compute bootstrap curves while preserving class prevalence."""
    rng = np.random.RandomState(random_state)
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]

    sampled_curves = []
    for _ in range(n_bootstraps):
        b_pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        b_neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        idx = np.concatenate([b_pos, b_neg])
        try:
            x, y, _ = curve_func(y_true[idx], y_score[idx])
            sampled_curves.append((x, y))
        except Exception:
            continue
    return sampled_curves


def percentile_band(sampled_curves, x_grid):
    """Interpolate sampled curves onto a common grid and return percentiles."""
    values = []
    for x, y in sampled_curves:
        if len(x) < 2:
            continue
        # Ensure monotonic x for interpolation
        x = np.asarray(x)
        y = np.asarray(y)
        if np.any(np.diff(x) <= 0):
            # For PR curves, recall can be non-monotonic; sort by x
            order = np.argsort(x)
            x, y = x[order], y[order]
            # remove duplicates
            uniq = np.concatenate(([True], np.diff(x) > 0))
            x, y = x[uniq], y[uniq]
        interp_y = np.interp(x_grid, x, y, left=0.0, right=y[-1] if len(y) else 0.0)
        values.append(interp_y)
    if not values:
        return np.full_like(x_grid, np.nan), np.full_like(x_grid, np.nan)
    values = np.array(values)
    return np.percentile(values, 2.5, axis=0), np.percentile(values, 97.5, axis=0)


def main():
    pred = pd.read_csv(PRED_PATH)
    y_true = pred["y_true"].values
    y_score = pred["cal_prob"].values

    auc = roc_auc_score(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # --- Panel a: ROC ---
    ax = axes[0]
    fpr, tpr, _ = roc_curve(y_true, y_score)
    ax.plot(fpr, tpr, lw=2, color="#1f77b4", label=f"AGVarPred (AUC = {auc:.4f})")
    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)

    roc_curves = stratified_bootstrap_curves(y_true, y_score, roc_curve)
    fpr_grid = np.linspace(0, 1, 100)
    low, high = percentile_band(roc_curves, fpr_grid)
    ax.fill_between(fpr_grid, low, high, color="#1f77b4", alpha=0.2)

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("(a) ROC curve")
    ax.legend(loc="lower right", frameon=False)

    # --- Panel b: Precision-Recall ---
    ax = axes[1]
    prec, rec, _ = precision_recall_curve(y_true, y_score)
    ax.plot(rec, prec, lw=2, color="#ff7f0e", label=f"AGVarPred (AUPRC = {auprc:.4f})")

    pr_curves = stratified_bootstrap_curves(y_true, y_score, precision_recall_curve)
    rec_grid = np.linspace(0, 1, 100)
    low, high = percentile_band(pr_curves, rec_grid)
    ax.fill_between(rec_grid, low, high, color="#ff7f0e", alpha=0.2)

    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("(b) Precision-recall curve")
    ax.legend(loc="lower left", frameon=False)

    # --- Panel c: Calibration ---
    ax = axes[2]
    n_bins = 10
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    bin_centers = (bin_boundaries[:-1] + bin_boundaries[1:]) / 2
    bin_ids = np.digitize(y_score, bin_boundaries) - 1
    bin_ids = np.clip(bin_ids, 0, n_bins - 1)
    frac_pos = np.array([y_true[bin_ids == i].mean() if (bin_ids == i).any() else np.nan for i in range(n_bins)])
    mean_pred = np.array([y_score[bin_ids == i].mean() if (bin_ids == i).any() else np.nan for i in range(n_bins)])
    mask = ~np.isnan(frac_pos)
    ax.plot(mean_pred[mask], frac_pos[mask], "o-", color="#2ca02c", label="AGVarPred")
    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
    ax.set_xlabel("Mean predicted probability")
    ax.set_ylabel("Fraction of positives")
    ax.set_title("(c) Calibration")
    ax.legend(loc="upper left", frameon=False)

    # ECE annotation
    ece = sum(
        np.abs(frac_pos[i] - mean_pred[i]) * (bin_ids == i).sum() / len(y_score)
        for i in range(n_bins) if mask[i]
    )
    ax.text(0.05, 0.75, f"ECE = {ece:.4f}\nN = {len(y_true):,}", transform=ax.transAxes, fontsize=10, verticalalignment="top")

    plt.suptitle("AGVarPred performance on the held-out ClinVar test set", fontsize=14, y=1.02)
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        out_path = os.path.join(OUTDIR, f"figure2_combined.{ext}")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Saved {out_path}")

    plt.close()


if __name__ == "__main__":
    main()
