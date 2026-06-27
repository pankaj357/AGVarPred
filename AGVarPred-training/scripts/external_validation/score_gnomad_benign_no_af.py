#!/usr/bin/env python3
"""Score gnomAD Common Benign benchmark with the no-AF regularized model."""

import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from score_gnomad_benign import (
    load_extracted_features, load_benign_labels, score_benign, save_benign_results
)
from scoring_core import load_pipeline

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[3]))
EXT_DIR = ROOT / "external_validation"
FEATURE_DIR = EXT_DIR / "processing/features/gnomad_benign"
BENCHMARK_PATH = EXT_DIR / "benchmarks/benchmark_gnomad_benign.csv"
MODEL_PATH = ROOT / "model_6_minus_af_output/Model_1_no_AF_pipeline.pkl"
OUTDIR = EXT_DIR / "results/gnomad_benign/regularized_no_af"

if __name__ == "__main__":
    print("="*70)
    print(f"GNOMAD COMMON BENIGN BENCHMARK SCORING — Model: regularized_no_af")
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
