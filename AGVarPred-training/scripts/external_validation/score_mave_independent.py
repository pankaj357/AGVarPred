#!/usr/bin/env python3
"""Score MAVE Independent benchmark with the regularized model.

Uses the shared score_benchmark pipeline so classification metrics are
methodologically identical to the other clinical benchmarks. MAVE-specific
correlation analysis between predicted probabilities and continuous MAVE
scores is computed and saved as an additional output.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, pearsonr

sys.path.insert(0, str(Path(__file__).parent))

from scoring_core import (
    ROOT, load_extracted_features, load_pipeline, load_labels,
    load_train_genes, score_benchmark, save_results
)

EXT_DIR = ROOT / "external_validation"
FEATURE_DIR = EXT_DIR / "processing/features/mave_independent"
BENCHMARK_PATH = EXT_DIR / "benchmarks/benchmark_mave_independent.csv"
VEP_PATH = EXT_DIR / "vep/mave_independent_vep.parquet"
BENCHMARK_NAME = "mave_independent"


def compute_mave_correlations(pred_df):
    """Compute correlation between calibrated probabilities and MAVE scores."""
    valid = pred_df.dropna(subset=["mave_score"])
    corr_results = {}
    if len(valid) > 1:
        rho, rho_p = spearmanr(valid["prob_calibrated"], valid["mave_score"])
        r, r_p = pearsonr(valid["prob_calibrated"], valid["mave_score"])
        corr_results = {
            "score_spearman_rho": float(rho),
            "score_spearman_p": float(rho_p),
            "score_pearson_r": float(r),
            "score_pearson_p": float(r_p),
            "n_with_score": int(len(valid)),
        }
        print(f"   Spearman ρ = {rho:.4f} (p={rho_p:.2e})")
        print(f"   Pearson r  = {r:.4f} (p={r_p:.2e})")
    else:
        print("   ⚠ Not enough variants with MAVE scores for correlation")
    return corr_results


def score_mave(model_path, outdir):
    print("="*70)
    print(f"MAVE INDEPENDENT BENCHMARK SCORING — Model: {Path(model_path).stem}")
    print("="*70)

    feature_df = load_extracted_features(FEATURE_DIR)
    pipe = load_pipeline(model_path)
    labels_df = load_labels(BENCHMARK_PATH)
    train_genes = load_train_genes()

    # Preserve continuous MAVE score for correlation analysis
    score_map = dict(zip(labels_df["variant_id"], labels_df["score"]))

    pred_df, metrics, subgroup_results, characterization = score_benchmark(
        feature_df, pipe, labels_df, train_genes, outdir,
        has_af=True, vep_path=VEP_PATH
    )

    pred_df["mave_score"] = pred_df["variant_id"].map(score_map)

    print("\n🔄 Computing correlations with MAVE scores...")
    corr_results = compute_mave_correlations(pred_df)

    # Add correlation results to characterization
    characterization["benchmark"] = "MAVE Independent (non-training genes)"
    characterization.update(corr_results)

    save_results(pred_df, metrics, subgroup_results, characterization, outdir, BENCHMARK_NAME)

    # Save correlation results separately for convenience
    corr_path = Path(outdir) / f"{BENCHMARK_NAME}_mave_correlation.json"
    with open(corr_path, "w") as fh:
        json.dump(corr_results, fh, indent=2)
    print(f"💾 MAVE correlation saved: {corr_path}")

    print("\n✅ Scoring complete!")
    print(f"   Output directory: {outdir}/")


if __name__ == "__main__":
    MODEL_PATH = ROOT / "final_model_output_regularized/final_pipeline.pkl"
    OUTDIR = EXT_DIR / "results/mave_independent/regularized"
    score_mave(MODEL_PATH, OUTDIR)
