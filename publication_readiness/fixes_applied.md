# Applied Fixes — Forensic Audit Issues

## 1. Internal variant-type alignment (CRITICAL)

**Problem:** `variant_type_analysis/compute_variant_type_analysis.py` concatenated internal test predictions and features by row index, but predictions were generated from an unsorted `glob` while the analysis used a sorted `glob`. Empirical label agreement was only ~51.7%.

**Fixes:**
- Created `regenerate_test_predictions_with_variant_id.py`, which loads the saved final pipeline, scores the test set in deterministic sorted-parquet order, and writes `final_model_output_regularized/test_predictions.csv` with a `variant_id` column.
- Updated `final_model_output_regularized/final_metrics.csv` with the recomputed 2,000-stratified-bootstrap CI.
- Modified `variant_type_analysis/compute_variant_type_analysis.py` to merge predictions and features on `variant_id` instead of row index.
- Updated `AGVarPred/AGVarPred-training/scripts/model_training/regularized_final_model.py` and the top-level `regularized_final_model.py` to:
  - load parquet parts in sorted order;
  - include `variant_id` in `test_predictions.csv` for future runs.
- Re-ran `compute_variant_type_analysis.py`; regenerated `internal_test_variant_type_metrics.csv`, `external_variant_type_metrics.csv`, `variant_type_coverage.csv`, `non_missense_comparison.csv`, and Figures 4a/4b/5.

**Result:** Internal per-class AUCs are now correctly aligned. The corrected values differ substantially from the prior misaligned report (e.g., "other" AUC = 0.9413 on 7,352 variants; synonymous class has no pathogenic cases so AUC is undefined).

## 2. VIP removed from independent external benchmarks (CRITICAL)

**Problem:** The VIP benchmark contains 1,075 pathogenic labels derived from ClinVar (`VIP_review+ClinVar_*`), creating a label-source circularity risk because AGVarPred is trained on ClinVar.

**Decision:** VIP was **removed from the main independent external-benchmark set** and is now reported separately as a **ClinVar-held-out-gene consistency analysis**.

Exclusion of the ClinVar-derived labels was attempted first, but it removes **all** pathogenic VIP variants, leaving only 2,805 benign variants and making AUC undefined. Therefore, the benchmark is retained only as a consistency/robustness analysis, not as independent validation.

**Fixes:**
- Removed `vip` from the `BENCHMARKS` dictionaries in:
  - `delong_analysis/compute_delong_tests.py`
  - `variant_type_analysis/compute_variant_type_analysis.py`
  - `dca_analysis/compute_decision_curve_analysis.py`
  - `external_validation/benchmarks_comparative/generate_comprehensive_comparison_report.py`
- Removed hardcoded `vip` iterations from `compute_delong_tests.py` plotting code.
- Re-ran `compute_delong_tests.py`, `compute_variant_type_analysis.py`, `compute_decision_curve_analysis.py`, and `generate_comprehensive_comparison_report.py`.
- Updated `external_validation/benchmarks_comparative/COMPREHENSIVE_COMPARISON_REPORT.md` to remove VIP from the summary table and interpretation bullets.
- Created `publication_readiness/vip_consistency_analysis.md` documenting VIP as a separate consistency analysis.
- Updated `external_validation/scripts/build_vip_benchmark.py` to print a warning at build time.
- Updated manuscript text (Results, Methods, Abstract, figure captions) to reflect five independent benchmarks and VIP as a separate consistency analysis.

**Result:** Main external-benchmark figures and statistics now include only Humsavar, MAVE Independent, gnomAD Common Benign, Grimm2015, and DVD. DeLong FDR correction now spans 13 comparisons instead of 14.

## 3. Bootstrap CI procedure (MAJOR)

**Problem:** Methods claimed 2,000 stratified bootstrap resamples; code used 200 unstratified resamples.

**Fixes:**
- Regenerated the internal test CI using 2,000 stratified resamples: **0.9731–0.9775**.
- Updated `AGVarPred/AGVarPred-training/scripts/model_training/regularized_final_model.py`, top-level `regularized_final_model.py`, and `ablation_feature_groups_minus_af.py` to implement 2,000 stratified bootstrap resamples for future runs.
- Updated CI references in:
  - `publication_readiness/abstract.md`
  - `publication_readiness/results_section.md`
  - `publication_readiness/methods_section.md`
  - `final_model_output_regularized/final_metrics.csv`

**Result:** Code and manuscript now match: 2,000 stratified resamples.

## 4. Figure captions reconciled with actual plots

**Problem:** Several captions did not describe the actual figure content.

**Fixes:**
- Created `publication_readiness/generate_figure2.py`, which produces a single 3-panel Figure 2 with ROC, PR, and calibration panels and shaded 95% bootstrap confidence intervals. Output saved to `final_model_output_regularized/figure2_combined.pdf/png`.
- Updated `publication_readiness/figure_plan.md` captions for Figures 2–7 to match the actual plotted content.
- Updated `publication_readiness/figure_audit.csv` to point Figure 2 to the new combined file.

**Result:** Figure captions now accurately describe what is plotted.

## Additional corrections applied during this pass

- `publication_readiness/results_section.md` and `methods_section.md`: corrected per-class internal AUCs, SHAP claims (VEP IMPACT > top AlphaGenome; splicing 30.3%; regulatory 67.5%; H3K9me3 not in top 50), and VIP framing.
- `publication_readiness/abstract.md`: updated CI, changed "six independent benchmarks" to "five independent benchmarks", added VIP consistency note, and softened "broad variant coverage" to "broad SNV consequence-class coverage".

## Files changed

- `regenerate_test_predictions_with_variant_id.py` (new)
- `publication_readiness/generate_figure2.py` (new)
- `publication_readiness/vip_consistency_analysis.md` (new)
- `external_validation/scripts/recompute_vip_metrics_without_clinvar.py` (new)
- `AGVarPred/AGVarPred-training/scripts/model_training/regularized_final_model.py`
- `regularized_final_model.py`
- `ablation_feature_groups_minus_af.py`
- `variant_type_analysis/compute_variant_type_analysis.py`
- `delong_analysis/compute_delong_tests.py`
- `dca_analysis/compute_decision_curve_analysis.py`
- `external_validation/benchmarks_comparative/generate_comprehensive_comparison_report.py`
- `external_validation/benchmarks_comparative/COMPREHENSIVE_COMPARISON_REPORT.md`
- `external_validation/scripts/build_vip_benchmark.py`
- `publication_readiness/figure_plan.md`
- `publication_readiness/figure_audit.csv`
- `publication_readiness/abstract.md`
- `publication_readiness/methods_section.md`
- `publication_readiness/results_section.md`
- `final_model_output_regularized/test_predictions.csv`
- `final_model_output_regularized/final_metrics.csv`
- `final_model_output_regularized/figure2_combined.pdf/png` (new)
- `delong_analysis/results/*.csv`
- `delong_analysis/figures/delong_significance.pdf/png`
- `variant_type_analysis/internal_test_variant_type_metrics.csv`
- `variant_type_analysis/external_variant_type_metrics.csv`
- `variant_type_analysis/variant_type_coverage.csv`
- `variant_type_analysis/non_missense_comparison.csv`
- `variant_type_analysis/variant_type_auc.pdf/png`
- `variant_type_analysis/variant_type_coverage.pdf/png`
- `variant_type_analysis/non_missense_auc_summary.pdf/png`
- `dca_analysis/results/decision_curve_data.csv`
- `dca_analysis/results/decision_curve_summary.csv`
- `dca_analysis/figures/decision_curve_analysis.pdf/png`
- `dca_analysis/figures/decision_curve_internal_test.pdf/png`
