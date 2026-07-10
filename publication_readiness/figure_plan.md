# Figure Plan — AGVarPred Manuscript

**Date:** 2026-06-29  
**Scope:** Main-text and supplementary figures for the corrected AGVarPred manuscript.

This plan is aligned with `figure_audit.csv`, which records the source script, input files, and output files for each figure. The numbering below follows the order in which figures are cited in the manuscript.

---

## Main-text figures

### Figure 1 — AGVarPred workflow overview

- **Status:** Needs to be created.
- **Placement:** Methods §1.
- **Source:** `figure_audit.csv` row 1.
- **Suggested content:**
  1. ClinVar `variant_summary.txt` → filter → 134,002-variant gold standard.
  2. Gene-level 80/10/10 split.
  3. Feature extraction: AlphaGenome SDK, gnomAD AF, Ensembl VEP.
  4. Feature selection: nested elimination → 120 features.
  5. Model: LightGBM → isotonic calibration → threshold 0.42.
  6. Evaluation: internal ClinVar test subset + five independent external benchmarks + separate VIP consistency analysis.
- **Citation text:** "Figure 1 provides an overview of the AGVarPred pipeline."

### Figure 2 — Internal ClinVar test-set performance

- **Status:** Ready; combine three existing plots into a single 3-panel figure.
- **Source:** `figure_audit.csv` row 2.
- **Source file:** `final_model_output_regularized/figure2_combined.pdf/png`
- **Placement:** Results §1.
- **Caption suggestion:** "Performance of AGVarPred on the 12,385 gene-held-out ClinVar test variants represented in the final feature matrices. (a) Receiver-operating-characteristic curve (ROC-AUC = 0.9753, 95% CI 0.9731–0.9775). (b) Precision-recall curve (PR-AUC = 0.9606). (c) Calibration curve with expected calibration error (ECE = 0.0225). Shaded regions in (a) and (b) show 95% bootstrap confidence intervals from 2,000 stratified resamples."
- **Citation text:** "Figure 2 shows the ROC, precision-recall, and calibration curves on the scored internal test set."

### Figure 3 — External benchmark comparison with DeLong significance

- **Status:** Ready.
- **Source:** `figure_audit.csv` row 3.
- **Source file:** `delong_analysis/figures/delong_significance.pdf`
- **Placement:** Results §2.
- **Caption suggestion:** "External benchmark performance and pairwise DeLong comparisons between AGVarPred and CADD, AlphaMissense, and REVEL on five independent benchmarks. Bars show ROC-AUC on the pairwise common subset of variants scored by both AGVarPred and the competitor; stars above each competitor bar indicate DeLong test significance on that common subset (q < 0.05 after Benjamini–Hochberg correction across 11 valid pairwise tests). 'ns' = not significant. The gnomAD Common Benign benchmark is omitted because it contains no pathogenic variants and AUC cannot be computed; the VIP Database is omitted because its pathogenic labels are ClinVar-derived and is reported separately as a ClinVar-held-out-gene consistency analysis."
- **Citation text:** "Figure 3 summarizes external benchmark performance and DeLong significance comparisons."

### Figure 4 — Variant-type performance and coverage

- **Status:** Ready.
- **Source:** `figure_audit.csv` row 4.
- **Source files:**
  - `variant_type_analysis/variant_type_auc.pdf`
  - `variant_type_analysis/variant_type_coverage.pdf`
- **Placement:** Results §3.
- **Caption suggestion:** "Variant-type performance and coverage across benchmarks. (a) AGVarPred performance by variant consequence class across the internal ClinVar test set and external benchmarks. For the all-benign gnomAD Common Benign benchmark, bars show 1 − false positive rate rather than ROC-AUC. Estimates from small subgroups are flagged as exploratory. (b) Percentage of variants scored by each tool within each consequence class."
- **Citation text:** "Figure 4 shows variant-type performance and tool coverage."

### Figure 5 — Non-missense benchmark comparison

- **Status:** Ready.
- **Source:** `figure_audit.csv` row 5.
- **Source file:** `variant_type_analysis/non_missense_auc_summary.pdf`
- **Placement:** Results §4.
- **Caption suggestion:** "Non-missense ROC-AUC comparison on the AGVarPred-scored universe across five independent benchmarks. Comparisons are restricted to variants not classified as missense and evaluated on the subset each competitor could score, avoiding coverage-biased artifacts. The gnomAD Common Benign group is empty because it contains no pathogenic non-missense variants."
- **Citation text:** "Figure 5 compares non-missense performance across benchmarks."

### Figure 6 — Decision-curve analysis

- **Status:** Ready.
- **Source:** `figure_audit.csv` row 6.
- **Source files:**
  - `dca_analysis/figures/decision_curve_internal_test.pdf`
  - `dca_analysis/figures/decision_curve_analysis.pdf`
- **Placement:** Results §6.
- **Caption suggestion:** "Decision-curve analysis of AGVarPred on the scored internal ClinVar test set and five independent external benchmarks. Net benefit is shown relative to treat-all and treat-none reference strategies. (a) Internal test set. (b) External benchmarks."
- **Citation text:** "Figure 6 shows decision-curve analysis on the internal test set and external benchmarks."

### Figure 7 — Tissue and assay-type SHAP importance

- **Status:** Ready.
- **Source:** `figure_audit.csv` row 7.
- **Source files:**
  - `tissue_shap_analysis/tissue_assay_heatmap.pdf`
  - `tissue_shap_analysis/cumulative_tissue_contribution.pdf`
- **Placement:** Results §8.
- **Caption suggestion:** "SHAP-based interpretation of AlphaGenome feature usage. (a) Heatmap of SHAP importance by tissue group and assay type for the top eight tissue groups. (b) Cumulative SHAP contribution across broad tissue groups ranked by total contribution."
- **Citation text:** "Figure 7 summarizes tissue and assay-type SHAP importance."

---

## Supplementary figures

| Supplementary figure | Source file(s) | Content |
|----------------------|----------------|---------|
| Supplementary Figure 1 | `feature_selection_output_nested/feature_selection_curve_auc.png`<br>`feature_selection_output_nested/feature_selection_curve_auprc.png` | Feature-selection curves showing cross-validated AUC and AUPRC across panel sizes. |
| Supplementary Figure 2 | `ablation_output/ablation_curve.png` | Feature-ablation curve. |
| Supplementary Figure 3 | `final_model_output_regularized/shap_summary.png` | Global SHAP summary of top individual features. |
| Supplementary Figure 4 | `final_model_output_regularized/shap_category_barplot.png`<br>`final_model_output_regularized/shap_category_piechart.png` | SHAP contribution by assay category. |
| Supplementary Figure 5 | `ablation_feature_groups_output/Model_1_roc_curve.png`<br>`ablation_feature_groups_output/Model_1_pr_curve.png`<br>`ablation_feature_groups_output/Model_1_calibration_curve.png` | Ablation-model performance curves. |
| Supplementary Figure 6 | `model_6_minus_af_output/Model_1_no_AF_roc_curve.png`<br>`model_6_minus_af_output/Model_1_no_AF_pr_curve.png` | No-allele-frequency model performance curves. |

---

## Figure-insertion checklist

- [ ] Create Figure 1 (workflow overview).
- [ ] Combine `roc_curve.png`, `pr_curve.png`, and `calibration_curve.png` into Figure 2.
- [ ] Use `delong_significance.pdf` for Figure 3.
- [ ] Combine `variant_type_auc.pdf` and `variant_type_coverage.pdf` into Figure 4.
- [ ] Use `non_missense_auc_summary.pdf` for Figure 5.
- [ ] Combine internal-test and external-benchmark DCA plots into Figure 6.
- [ ] Combine tissue-assay heatmap and cumulative contribution into Figure 7.
- [ ] Add figure references in `methods_section.md` and `results_section.md` (done for Figures 1–7).
- [ ] Add captions to the manuscript file or supplementary caption list.
- [ ] Ensure font sizes and line weights are journal-compliant.
- [ ] Convert PNGs to vector/PDF where possible for submission.

---

## Suggested citation locations in the text (7-figure plan)

| Figure | Section | Suggested sentence |
|--------|---------|--------------------|
| Figure 1 | Methods §1 | "Figure 1 provides an overview of the AGVarPred pipeline." |
| Figure 2 | Results §1 | "Figure 2 shows the ROC, precision-recall, and calibration curves on the scored internal test set." |
| Figure 3 | Results §2 | "Figure 3 summarizes external benchmark performance and DeLong significance comparisons." |
| Figure 4 | Results §3 | "Figure 4 shows variant-type performance and tool coverage." |
| Figure 5 | Results §4 | "Figure 5 compares non-missense performance across benchmarks." |
| Figure 6 | Results §6 | "Figure 6 shows decision-curve analysis on the internal test set and external benchmarks." |
| Figure 7 | Results §8 | "Figure 7 summarizes tissue and assay-type SHAP importance." |
