#!/usr/bin/env python3
"""
Generate publication readiness report and audit tables.
No retraining, rescoring, or new biological analyses.
"""

import os
import pandas as pd
from pathlib import Path

OUT_DIR = Path(__file__).parent
ROOT_DIR = OUT_DIR.parent


def file_exists(path):
    return (ROOT_DIR / path).exists()


# ---------------------------------------------------------------------------
# Figure audit
# ---------------------------------------------------------------------------
figures = [
    # (figure_number, caption, source_script, input_files, output_file, status)
    ("Fig 1", "AGVarPred workflow overview", "TBD (drawn externally)", "N/A", "figures/fig1_workflow.pdf", "MISSING"),
    ("Fig 2", "Internal test performance: ROC, PR, calibration", "regularized_final_model.py", "final_dataset_parts_train/*.parquet, final_dataset_parts_cal/*.parquet, final_dataset_parts_test/*.parquet", "final_model_output_regularized/roc_curve.png, pr_curve.png, calibration_curve.png", "READY"),
    ("Fig 3", "External benchmark comparison with DeLong significance", "delong_analysis/compute_delong_tests.py", "external_validation/results/*/regularized/*_predictions.csv, external_validation/benchmarks_comparative/datasets/*/*", "delong_analysis/figures/delong_significance.pdf", "READY"),
    ("Fig 4", "Variant-type performance and coverage", "variant_type_analysis/compute_variant_type_analysis.py", "external_validation/results/*/regularized/*_predictions.csv, external_validation/*_vep.parquet", "variant_type_analysis/variant_type_auc.pdf, variant_type_coverage.pdf", "READY"),
    ("Fig 5", "Non-missense benchmark comparison", "variant_type_analysis/compute_variant_type_analysis.py", "external_validation/results/*/regularized/*_predictions.csv, external_validation/benchmarks_comparative/datasets/*/*", "variant_type_analysis/non_missense_auc_summary.pdf", "READY"),
    ("Fig 6", "Tissue and assay-type SHAP importance", "tissue_shap_analysis/compute_tissue_shap_analysis.py", "final_model_output_regularized/shap_values_sample.csv, feature_importance_shap.csv", "tissue_shap_analysis/tissue_assay_heatmap.pdf, cumulative_tissue_contribution.pdf", "READY"),
    ("Fig 7", "Decision Curve Analysis", "dca_analysis/compute_decision_curve_analysis.py", "final_model_output_regularized/test_predictions.csv, external_validation/results/*/regularized/*_predictions.csv", "dca_analysis/figures/decision_curve_analysis.pdf, decision_curve_internal_test.pdf", "READY"),
]
fig_df = pd.DataFrame(figures, columns=["figure_number", "caption", "source_script", "input_files", "output_file", "status"])
fig_df.to_csv(OUT_DIR / "figure_audit.csv", index=False)

# ---------------------------------------------------------------------------
# Table audit
# ---------------------------------------------------------------------------
tables = [
    # (table_number, caption, source_file, status, missing_values)
    ("Table 1", "ClinVar gold-standard dataset and splits", "clinvar_gold_grch38_clean.csv, train.csv, cal.csv, test.csv", "READY", "Exact ClinVar release date not recorded in repository"),
    ("Table 2", "External benchmark summary", "external_validation/benchmarks_comparative/COMPREHENSIVE_COMPARISON_REPORT.md", "READY", "None"),
    ("Table 3", "DeLong statistical significance results", "delong_analysis/results/delong_results_fdr.csv", "READY", "None"),
    ("Table 4", "Variant-type performance metrics", "variant_type_analysis/external_variant_type_metrics.csv", "READY", "None"),
    ("Table 5", "Variant-type coverage comparison", "variant_type_analysis/variant_type_coverage.csv", "READY", "None"),
    ("Table 6", "Non-missense comparison", "variant_type_analysis/non_missense_comparison.csv", "READY", "None"),
    ("Table 7", "Tissue/cell-type SHAP importance", "tissue_shap_analysis/tissue_shap_importance.csv", "READY", "None"),
    ("Table 8", "Decision Curve Analysis summary", "dca_analysis/results/decision_curve_summary.csv", "READY", "None"),
]
tab_df = pd.DataFrame(tables, columns=["table_number", "caption", "source_file", "status", "missing_values"])
tab_df.to_csv(OUT_DIR / "table_audit.csv", index=False)

# ---------------------------------------------------------------------------
# Supplementary checklist
# ---------------------------------------------------------------------------
supp = [
    # (item, type, status, notes)
    ("Supplementary Figure S1", "Supplementary Figure", "READY", "ablation_output/ablation_curve.png"),
    ("Supplementary Figure S2", "Supplementary Figure", "READY", "ablation_feature_groups_output/* feature-group ablation"),
    ("Supplementary Figure S3", "Supplementary Figure", "READY", "final_model_output_regularized/shap_summary.png"),
    ("Supplementary Figure S4", "Supplementary Figure", "READY", "tissue_shap_analysis/top_biosamples_barplot.pdf"),
    ("Supplementary Figure S5", "Supplementary Figure", "READY", "dca_analysis/figures/decision_curve_internal_test.pdf"),
    ("Supplementary Table S1", "Supplementary Table", "READY", "Full selected feature list: feature_selection_output_nested/selected_features.txt"),
    ("Supplementary Table S2", "Supplementary Table", "READY", "Full SHAP values: final_model_output_regularized/shap_values_sample.csv"),
    ("Supplementary Table S3", "Supplementary Table", "READY", "DeLong raw results: delong_analysis/results/delong_results.csv"),
    ("Supplementary Table S4", "Supplementary Table", "READY", "Per-biosample SHAP: tissue_shap_analysis/biosample_shap_importance.csv"),
    ("Supplementary Methods", "Supplementary Methods", "NEEDS PREPARATION", "Detailed pipeline commands and environment setup"),
    ("Supplementary Note 1", "Supplementary Note", "NEEDS PREPARATION", "Leakage-prevention statistics table"),
    ("Supplementary Note 2", "Supplementary Note", "NEEDS PREPARATION", "Benchmark independence audit"),
]
supp_df = pd.DataFrame(supp, columns=["item", "type", "status", "notes"])
supp_df.to_csv(OUT_DIR / "supplementary_checklist.csv", index=False)

print("Audit CSVs saved.")
