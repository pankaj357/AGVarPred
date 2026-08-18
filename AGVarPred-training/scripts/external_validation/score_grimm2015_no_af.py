#!/usr/bin/env python3
"""Score Grimm2015 benchmark with the non-regularized (original) model."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from scoring_core import (
    ROOT, load_extracted_features, load_pipeline, load_labels,
    load_train_genes, score_benchmark, save_results
)

EXT_DIR = ROOT / "external_validation"
FEATURE_DIR = EXT_DIR / "processing/features/grimm2015"
BENCHMARK_PATH = EXT_DIR / "benchmarks/benchmark_grimm2015.csv"
MODEL_PATH = ROOT / "model_6_minus_af_output/Model_1_no_AF_pipeline.pkl"
VEP_PATH = EXT_DIR / "vep_preprocessed/grimm2015_vep.parquet"
OUTDIR = EXT_DIR / "results/grimm2015/regularized_no_af"
BENCHMARK_NAME = "grimm2015"

if __name__ == "__main__":
    print("="*70)
    print(f"GRIMM2015 BENCHMARK SCORING — Model: regularized_no_af")
    print("="*70)

    feature_df = load_extracted_features(FEATURE_DIR)
    pipe = load_pipeline(MODEL_PATH)
    labels_df = load_labels(BENCHMARK_PATH)
    train_genes = load_train_genes()

    pred_df, metrics, subgroup_results, characterization = score_benchmark(
        feature_df, pipe, labels_df, train_genes, OUTDIR,
        has_af=True, vep_path=VEP_PATH
    )
    save_results(pred_df, metrics, subgroup_results, characterization, OUTDIR, BENCHMARK_NAME)

    print("\n✅ Scoring complete!")
    print(f"   Output directory: {OUTDIR}/")
