#!/usr/bin/env python3
"""Score ClinVar 3+ Star benchmark with the regularized model."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from scoring_core import (
    ROOT, load_extracted_features, load_pipeline, load_labels,
    load_train_genes, score_benchmark, save_results
)

EXT_DIR = ROOT / "external_validation"
FEATURE_DIR = EXT_DIR / "processing/features/clinvar_3star"
BENCHMARK_PATH = EXT_DIR / "benchmarks/benchmark_clinvar_3star.csv"
MODEL_PATH = ROOT / "final_model_output_regularized/final_pipeline.pkl"
VEP_PATH = EXT_DIR / "vep/clinvar_3star_vep.parquet"
OUTDIR = EXT_DIR / "results/clinvar_3star/regularized"
BENCHMARK_NAME = "clinvar_3star"

if __name__ == "__main__":
    print("="*70)
    print(f"CLINVAR 3+ STAR BENCHMARK SCORING — Model: regularized")
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
