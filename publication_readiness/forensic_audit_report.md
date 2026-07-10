# AGVarPred Forensic Audit Report

**Auditor:** Independent scientific auditor (Genome Biology / Nature Genetics / Bioinformatics / Briefings in Bioinformatics reviewer perspective).  
**Scope:** Full repository at `/data/kanaka/pankaj/feature_extration/new`, manuscript `AGVarPred.docx`, and supporting `publication_readiness/` files.  
**Date:** 2026-07-10.  
**Mandate:** Verify every manuscript claim against code, scripts, outputs, figures, and tables. Do not trust the manuscript. Report any unsupported or false claim, numerical inconsistency, methodological error, or reproducibility problem.

---

## 1. Executive Summary

The AGVarPred project implements a well-engineered LightGBM pathogenicity predictor with strong safeguards against data leakage (gene-level train/cal/test split, GroupKFold, calibration-only threshold selection). Most headline numbers in the manuscript are reproducible from repository outputs.

However, the audit identified **three critical issues** that must be resolved before submission to a top-tier journal:

1. **23.4% of the test set is silently discarded before feature extraction** by a hard-coded 2,000-variant-per-gene cap and a `<10 variants/gene` filter in `build_alphagenome_input_split.py`. The manuscript now discloses the resulting 12,385 scored test variants, but the causal filtering is not transparently described.
2. **Internal variant-type analysis is unreliable.** Predictions and features are concatenated by row index, but predictions were generated from an unsorted `glob` while the analysis uses a sorted `glob`. Empirically, only ~51.7% of labels align between the two orderings, so the internal per-class AUCs in Figure 4a / `internal_test_variant_type_metrics.csv` are effectively random.
3. **VIP benchmark contains 1,075 ClinVar-derived pathogenic labels (27.7% of VIP variants).** Because AGVarPred is trained on ClinVar, this creates a label-source circularity risk that undermines the benchmark's independence claims.

In addition, several methods claims are materially inaccurate: the bootstrap CI uses 200 unstratified resamples (not 2,000 stratified), VEP columns are treated as numeric (not categorical), the feature-selection "1-SE rule" actually uses the standard deviation, and the raw feature count is 4,702 (not ~4,631).

Multiple figure captions do not match the plotted content (Figure 2 bootstrap CIs, Figure 4a competitors, Figure 5 error bars). Figure 1 is missing.

**Bottom line:** The manuscript is not yet submission-ready. The core modeling results appear broadly correct, but the variant-type analysis is broken, the benchmark circularity must be addressed or transparently disclosed, and all methods/figure inconsistencies require correction.

---

## 2. Verified Correct Items

| Claim | Evidence | Status |
|---|---|---|
| ClinVar gold standard: 134,002 variants (40,584 P, 93,418 B) | `build_clean_clinvar_gold.py` output; `clinvar_gold_grch38_clean.csv` | ✅ |
| Gene-level 80/10/10 split, `random_state=42` | `split_clinvar_clean.py` L21–32 | ✅ |
| Pre-feature-extraction split sizes: train 101,604; cal 16,223; test 16,175 | `AGVarPred-zenodo/splits/*.csv` | ✅ |
| Zero gene overlap between splits | `split_clinvar_clean.py` L54–61 | ✅ |
| Scored matrix sizes: train 86,378; cal 13,201; test 12,385 | `feature_labbling.log` L60,70,80; parquet counts | ✅ |
| AlphaGenome SDK v0.6.1 | `manifest.yaml`; `.report_env` `pip show` | ✅ |
| gnomAD exomes r2.1.1 liftover GRCh38 | `feature_labbling.py` L38–41; `data_manifest.json` | ✅ |
| Selected feature panel: 120 features (111 AlphaGenome + 1 `gnomAD_AF` + 8 VEP) | `feature_selection_output_nested/selected_features.txt` | ✅ |
| No-AF model: 119 features, threshold 0.45 | `ablation_feature_groups_minus_af.py`; `model_6_minus_af_output/Model_1_no_AF_threshold.txt` | ✅ |
| Final model threshold 0.42 | `final_model_output_regularized/cv_thresholds.txt` | ✅ |
| Optuna best hyperparameters (n_estimators=485, learning_rate=0.0653, etc.) | `regularized_final_model.log` L115–116 | ✅ |
| Isotonic calibration on cal set | `regularized_final_model.py` L229–247 | ✅ |
| Internal test metrics: AUC 0.9753, AUPRC 0.9606, F1 0.8945, MCC 0.8327, Brier 0.0568, ECE 0.0225 | `final_model_output_regularized/final_metrics.csv` | ✅ |
| Test prevalence ~37% (4,576/12,385) | `final_dataset_parts_test` parquet counts | ✅ |
| Calibration prevalence ~41.8% (5,522/13,201) | `final_dataset_parts_cal` parquet counts | ✅ |
| External benchmark AUC/PR-AUC/MCC values and DeLong q-values | `COMPREHENSIVE_COMPARISON_REPORT.md`; `delong_analysis/results/delong_results_fdr.csv` | ✅ |
| gnomAD AF dominates SHAP (4.648) | `final_model_output_regularized/feature_importance_shap.csv` | ✅ |
| No test-set leakage into training/hyperparameter selection/calibration | Pipeline review | ✅ |

---

## 3. Critical Errors

### CR-1: Silent discard of 23.4% of the test set before feature extraction
- **Severity:** Critical
- **Manuscript location:** Methods §2.2–2.3; Results §3.1; Discussion/Limitations
- **Repository files:** `AGVarPred/AGVarPred-training/scripts/data_preparation/build_alphagenome_input_split.py` L50, L56
- **Evidence:**
  ```python
  if len(gdf) < 10:
      continue                              # drops entire gene
  gdf_sampled = gdf.sample(n=min(len(gdf), 2000), random_state=42)  # caps large genes
  ```
  For the test split: 16,175 CSV variants → 12,385 AlphaGenome inputs. Decomposition verified: 1,907 variants in 625 genes with `<10` variants are dropped; 1,883 BRCA2 variants are capped to 2,000. Total = 3,790.
- **Explanation:** The manuscript states the test set contains 12,385 scored variants but does not clearly state that this reduction is caused by arbitrary filters (a 2,000-variant/gene cap and a `<10 variants/gene` minimum) rather than API failures. The filter removes 68% of test genes and disproportionately affects small genes and BRCA2. This is a severe generalizability and reproducibility issue.
- **Recommended correction:** Either (a) remove the `<10` and `2000` filters and rerun feature extraction, or (b) fully disclose the filters, report the list of discarded genes/variants, and justify why the remaining 12,385 variants are representative. The current Limitations paragraph is insufficient because it does not name the cause.

### CR-2: Internal variant-type analysis is invalid due to row-order misalignment
- **Severity:** Critical
- **Manuscript location:** Results §3.4 / §3.8; Figure 4a; `variant_type_analysis/internal_test_variant_type_metrics.csv`
- **Repository files:** `variant_type_analysis/compute_variant_type_analysis.py` L278–288; `regularized_final_model.py` prediction export
- **Evidence:**
  - Predictions generated by `regularized_final_model.py` load parquets with `glob(f"{path}/*.parquet")` (unsorted).
  - `variant_type_analysis.py` loads the same parquets with `sorted(glob.glob(...))` and concatenates predictions/features by row index: `df = pd.concat([pred.reset_index(drop=True), feats.reset_index(drop=True)], axis=1)`.
  - Empirical label agreement between sorted features and predictions: **51.7%** (effectively random).
- **Explanation:** The internal per-class AUCs reported in Results and plotted in Figure 4a are meaningless because `variant_id` is never used to align predictions to features. The external benchmark portions of the figure are aligned by `variant_id` and remain valid, but the internal_test bars are not.
- **Recommended correction:** Add `variant_id` to `test_predictions.csv` and merge `pred` with `feats` on `variant_id`. Recompute `internal_test_variant_type_metrics.csv` and regenerate Figure 4a.

### CR-3: VIP benchmark contains 1,075 ClinVar-derived pathogenic labels
- **Severity:** Critical
- **Manuscript location:** Results §3.3 (external benchmark comparison); Methods §2.10
- **Repository files:** `external_validation/scripts/build_vip_benchmark.py` L124–130; `external_validation/benchmarks/benchmark_vip.csv`
- **Evidence:**
  - VIP CSV contains 3,880 variants (1,075 P, 2,805 B).
  - `label_source` values include `VIP_review+ClinVar_Pathogenic` (667), `VIP_review+ClinVar_Likely pathogenic` (348), `VIP_review+ClinVar_Pathogenic/Likely pathogenic` (53), etc. — **1,075 pathogenic labels derived from ClinVar**.
  - Because training/calibration/test labels are ClinVar-derived, VIP shares the label source with the training data.
- **Explanation:** Although the specific VIP variants/genes are held out, the pathogenicity annotations for 27.7% of VIP come from the same database used to train AGVarPred. This is a circularity risk that a reviewer will flag as undermining independence claims.
- **Recommended correction:** Either (a) exclude the 1,075 ClinVar-derived VIP variants from the benchmark and recompute all VIP metrics, or (b) explicitly disclose the label-source overlap and report a sensitivity analysis with and without these variants.

---

## 4. Major Errors

### MAJ-1: Bootstrap CI procedure misreported
- **Severity:** Major
- **Manuscript location:** Methods §2.9
- **Repository files:** `AGVarPred/AGVarPred-training/scripts/model_training/regularized_final_model.py` L368–376
- **Evidence:** Manuscript claims "2,000 stratified resamples"; code uses `for _ in range(200)` with simple random resampling (`rng.choice(len(y_test), len(y_test), replace=True)`), no stratification.
- **Explanation:** The reported 95% CI [0.9730–0.9774] comes from 200 unstratified bootstraps, not the described procedure. Stratification is especially important given the ~37% prevalence.
- **Recommended correction:** Either implement 2,000 stratified bootstraps and update the CI, or change the Methods to state "200 unstratified bootstraps."

### MAJ-2: Figure 4a caption claims competitors are plotted, but only AGVarPred is shown
- **Severity:** Major
- **Manuscript location:** Figure 4 caption (panel a)
- **Repository files:** `variant_type_analysis/compute_variant_type_analysis.py`; `variant_type_analysis/variant_type_auc.pdf/png`
- **Evidence:** Figure title is "AGVarPred performance by variant type"; bars are colored by benchmark, not method. No CADD/AlphaMissense/REVEL per-class AUCs are plotted.
- **Explanation:** This is a direct contradiction between caption and figure content. It was already flagged in the prior audit brief.
- **Recommended correction:** Revise the caption to "AGVarPred performance by variant consequence class across benchmarks" or regenerate the figure with competitor bars.

### MAJ-3: AlphaGenome feature extraction silently drops variants on API/RPC errors
- **Severity:** Major
- **Manuscript location:** Methods §2.4
- **Repository files:** `AGVarPred/AGVarPred-training/scripts/feature_extraction/code.py` L120–122; `code_external_generic.py`; `code_external_validation.py`
- **Evidence:** `except Exception as e: print(...); continue` with no retry, no fallback, and no failed-variant manifest.
- **Explanation:** Any AlphaGenome failure removes the variant from the output parquet. Logs (`logs_*/run_*.log`) show `_MultiThreadedRendezvous` errors. This contributes to the train/cal gap beyond the input filters.
- **Recommended correction:** Add retry logic and emit a `failed_variants.csv` manifest so losses are auditable.

### MAJ-4: Multi-allelic gnomAD records use only the first ALT and first VEP entry
- **Severity:** Major
- **Manuscript location:** Methods §2.3 (gnomAD/VEP annotation)
- **Repository files:** `feature_labbling.py` L265–269, L283–286; `AGVarPred/src/agvarpred_core/gnomad.py`
- **Evidence:** `alt = rec.alts[0]`; `parse_vep(vep[0])`.
- **Explanation:** For multi-allelic gnomAD sites, the ClinVar ALT may not be `alts[0]`, and the VEP annotation may correspond to the wrong allele. AF values may be assigned to incorrect variants.
- **Recommended correction:** Iterate over all alts and match the VEP CSQ entry by allele.

### MAJ-5: Feature selection reuses the same folds for model selection and evaluation
- **Severity:** Major
- **Manuscript location:** Methods §2.5
- **Repository files:** `AGVarPred/AGVarPred-training/scripts/feature_selection/feature_selection.py` L174–234
- **Evidence:** For each feature-step size `n`, an importance model is fit, top-`n` features selected, a reduced model is trained, and evaluated on the same GroupKFold split. The 1-SE rule then picks `optimal_n` from those same folds.
- **Explanation:** This creates optimistic selection bias. There is no nested CV (outer loop for unbiased evaluation, inner loop for feature-set size selection).
- **Recommended correction:** Implement nested CV for feature-set size selection, or at minimum hold out a validation fold that is never used during feature-set selection.

### MAJ-6: Methods markdown file gives wrong final-feature composition
- **Severity:** Major
- **Manuscript location:** `publication_readiness/methods_section.md` L49
- **Repository files:** `feature_selection_output_nested/selected_features.txt`
- **Evidence:** Markdown states "111 AlphaGenome + 2 gnomAD + 7 VEP" (`gnomAD_AF` and `log10_gnomAD_AF`; `vep_is_LoF_HC`, `vep_has_SIFT`, `vep_has_PolyPhen`). Actual panel is 111 AlphaGenome + 1 `gnomAD_AF` + 8 VEP features.
- **Explanation:** The compiled `AGVarPred.docx` appears correct, but the source markdown is inconsistent with the trained model and the compiled manuscript.
- **Recommended correction:** Update `methods_section.md` to match the actual 111+1+8 panel.

---

## 5. Moderate Errors

### MOD-1: Raw feature count is 4,702, not ~4,631
- **Severity:** Moderate
- **Manuscript location:** Methods §2.3
- **Repository files:** `feature_selection.log` L3–4; `final_dataset_parts_train/*.parquet`
- **Evidence:** `feature_selection.log` reports `Dataset: (86378, 4706)` and `Numeric columns: 4702`. Direct inspection confirms 4,619 AlphaGenome assay tracks + 79 VEP-derived columns + 4 AF-derived columns = 4,702 raw features.
- **Recommended correction:** Change "approximately 4,631" to "4,702" and explain the breakdown.

### MOD-2: "1-standard-error rule" actually uses fold standard deviation
- **Severity:** Moderate
- **Manuscript location:** Methods §2.5
- **Repository files:** `feature_selection.py` L226–234
- **Evidence:** `best_std = std_auprc[best_idx]` (raw std), threshold = best_auprc − best_std. True 1-SE rule would use `std_auprc / sqrt(n_folds)`.
- **Recommended correction:** Rename to "1-standard-deviation rule" or implement the correct 1-SE rule. Note that `n=120` satisfies both, so the practical outcome may not change.

### MOD-3: VEP columns are numeric, not categorical
- **Severity:** Moderate
- **Manuscript location:** Methods §2.5
- **Repository files:** `feature_selection.py` L141–145; `feature_labbling.py`
- **Evidence:** `feature_selection.log` reports `Categorical columns: 0`. VEP-derived columns are stored as `int8`/`float32`, and `feature_selection.py` detects zero categorical columns.
- **Recommended correction:** Remove the claim that LightGBM treats VEP ordinals as categorical, or encode them as `pd.Categorical` if that was intended.

### MOD-4: Figure 2 caption promises a single 3-panel figure with bootstrap CIs; actual output is three separate PNGs without CIs
- **Severity:** Moderate
- **Manuscript location:** Figure 2 caption
- **Repository files:** `final_model_output_regularized/roc_curve.png`, `pr_curve.png`, `calibration_curve.png`; `regularized_final_model.py` L305–331
- **Evidence:** Three standalone 640×480 PNGs; no shaded confidence bands; no panel labels; AUC/AUPRC/ECE values are not annotated on the plots.
- **Recommended correction:** Combine into a single 3-panel figure with bootstrap confidence bands, or revise the caption.

### MOD-5: Figure 5 caption claims 95% bootstrap CI error bars; none are plotted
- **Severity:** Moderate
- **Manuscript location:** Figure 5 caption
- **Repository files:** `variant_type_analysis/non_missense_auc_summary.pdf/png`; `compute_variant_type_analysis.py` `plot_non_missense_roc()`
- **Evidence:** Bar plot has no error bars.
- **Recommended correction:** Add bootstrap CI error bars or revise the caption.

### MOD-6: gnomAD AF missing is conflated with AF = 0
- **Severity:** Moderate
- **Manuscript location:** Methods §2.3
- **Repository files:** `feature_labbling.py` L276–279, L458–462
- **Evidence:** `AF_missing = (gnomAD_AF == 0).astype('int8')`; variants not in gnomAD and variants with true AF=0 are encoded identically.
- **Recommended correction:** Distinguish "not in gnomAD" from "AF=0" using a separate indicator.

### MOD-7: `variant_benchmark_map.csv` still contains 35,299 contaminated ClinVar entries
- **Severity:** Moderate
- **Manuscript location:** External validation Methods
- **Repository files:** `external_validation/benchmarks/variant_benchmark_map.csv`; `external_validation/AUDIT_REPORT.md`; `external_validation/backup_contaminated_clinvar/`
- **Evidence:** The ClinVar 3★ benchmark was removed due to 88.4% gene overlap / 51.1% variant overlap, but the map file retains `clinvar` entries.
- **Recommended correction:** Rebuild `variant_benchmark_map.csv` from active benchmarks only; remove or quarantine stale backups.

### MOD-8: Tissue SHAP bootstrap quantifies uncertainty across features, not samples
- **Severity:** Moderate
- **Manuscript location:** Results §3.8
- **Repository files:** `tissue_shap_analysis/compute_tissue_shap_analysis.py` L236–246
- **Evidence:** Bootstrap resamples feature indices, not observations.
- **Recommended correction:** Clarify the interpretation or resample variants.

### MOD-9: Tissue/cell-type mapping is a hand-curated heuristic
- **Severity:** Moderate
- **Manuscript location:** Results §3.8
- **Repository files:** `tissue_shap_analysis/compute_tissue_shap_analysis.py` L108–229
- **Evidence:** Dozens of keyword rules (e.g., HCT116 → liver/digestive) without ontology validation.
- **Recommended correction:** Document the mapping and perform sensitivity analysis.

### MOD-10: External validation build script deduplicates against parquet parts, not authoritative CSVs
- **Severity:** Moderate
- **Manuscript location:** Methods §2.10
- **Repository files:** `external_validation/scripts/build_grimm2015_benchmark.py` L61–68
- **Evidence:** The script removes training variants by intersecting with `final_dataset_parts_{split}/` (111,964 scored variants) rather than `train.csv`/`cal.csv`/`test.csv` (134,002 variants). A direct check of Grimm2015 vs. authoritative CSVs found **zero actual overlap**, so no leakage occurred in this case, but the procedure is fragile.
- **Recommended correction:** Standardize all benchmark builders to deduplicate against the union of the authoritative split CSVs.

### MOD-11: Working execution directory is outside version control
- **Severity:** Moderate
- **Manuscript location:** Reproducibility statement
- **Repository files:** Top-level `/data/kanaka/pankaj/feature_extration/new/` vs. git-tracked `AGVarPred/`
- **Evidence:** 16 scripts are duplicated at top level; `external_validation/scripts/` at top level differs from `AGVarPred/AGVarPred-training/scripts/external_validation/` in 45+ files.
- **Recommended correction:** Move execution into the Git repo, reconcile divergent copies, commit or revert uncommitted changes.

---

## 6. Minor Errors

### MIN-1: External benchmark count typos
- **Severity:** Minor
- **Manuscript location:** Results §3.3
- **Repository files:** `external_validation/benchmarks/benchmark_independent_humsavar.csv`; `benchmark_grimm2015.csv`; `benchmark_mave_independent.csv`
- **Evidence:**
  - Humsavar scored: 14,718 (1,317 P / 13,401 B) — manuscript says 1,316/13,402.
  - Grimm2015 scored: 5,258 (512 P / 4,746 B) — manuscript says 509/4,749.
  - MAVE: 7,557 (3,751 P / 3,806 B, 25 genes) — manuscript says 3,757/3,800, 24 genes.
- **Recommended correction:** Update counts in the manuscript.

### MIN-2: DVD AlphaMissense coverage misattributed
- **Severity:** Minor
- **Manuscript location:** Results §3.5
- **Repository files:** `variant_type_analysis/variant_type_coverage.csv`
- **Evidence:** "1.4% of ‘other’ variants" is incorrect. DVD AlphaMissense overall coverage = 143/8,388 = 1.7%; "other" class coverage = 47.0%.
- **Recommended correction:** Replace with "1.7% overall coverage" or report the correct per-class figure.

### MIN-3: gnomAD competitor comparison Ns are inconsistent
- **Severity:** Minor
- **Manuscript location:** Results §3.3
- **Repository files:** `external_validation/benchmarks_comparative/COMPREHENSIVE_COMPARISON_REPORT.md`
- **Evidence:** Reported CADD/AlphaMissense/REVEL Ns (4,205 / 209 / 463) equal the non-missense subset, but the FPRs are full-coverage FPRs.
- **Recommended correction:** Report full-coverage Ns (4,989 / 516 / 1,235) or clarify the subset.

### MIN-4: SHAP claim overstates top AlphaGenome feature
- **Severity:** Minor
- **Manuscript location:** Results §3.7
- **Repository files:** `final_model_output_regularized/feature_importance_shap.csv`
- **Evidence:** "Top AlphaGenome features individually contributed more than any other single functional-annotation feature" — but `vep_IMPACT_score` SHAP = 0.437 vs. top AlphaGenome `CAGE_camera_type_eye` = 0.290.
- **Recommended correction:** Remove or qualify the claim.

### MIN-5: Splicing/regulatory SHAP percentages off by ~1 percentage point
- **Severity:** Minor
- **Manuscript location:** Results §3.8
- **Repository files:** `final_model_output_regularized/feature_importance_shap.csv`; `tissue_shap_analysis/parsed_features.csv`
- **Evidence:** Splicing = 30.3% (manuscript 31.3%); regulatory = 67.5% (manuscript 68.5%).
- **Recommended correction:** Recalculate or revise wording.

### MIN-6: H3K9me3 claimed in top 50 SHAP features; actual minimum rank is 57
- **Severity:** Minor
- **Manuscript location:** Results §3.8
- **Repository files:** `final_model_output_regularized/feature_importance_shap.csv`
- **Evidence:** Lowest H3K9me3 rank among AlphaGenome features = 57; none are in the top 50.
- **Recommended correction:** Remove H3K9me3 from the top-50 claim.

### MIN-7: Internal per-class AUC list omits large "other" class
- **Severity:** Minor
- **Manuscript location:** Results §3.4
- **Repository files:** `variant_type_analysis/internal_test_variant_type_metrics.csv`
- **Evidence:** File contains an `other` class with N=7,352 and AUC=0.9760, not reported in the per-class list.
- **Recommended correction:** Add the `other` class to the Results or explain its exclusion.

### MIN-8: Threshold grid excludes 0.95
- **Severity:** Minor
- **Manuscript location:** Methods §2.9
- **Repository files:** `regularized_final_model.py` L238
- **Evidence:** `np.arange(0.05, 0.95, 0.01)` stops at 0.94.
- **Recommended correction:** Change to `np.arange(0.05, 0.96, 0.01)` or state "0.05 to 0.94."

### MIN-9: Two final-model output directories exist
- **Severity:** Minor
- **Manuscript location:** N/A
- **Repository files:** `final_model_output/` and `final_model_output_regularized/`
- **Evidence:** `final_model_output/` (threshold 0.45, AUC 0.9737) and `final_model_output_regularized/` (threshold 0.42, AUC 0.9753).
- **Recommended correction:** Remove the superseded directory or clearly label it as obsolete.

### MIN-10: Figure 3 caption describes two panels; figure is a single composite
- **Severity:** Minor
- **Manuscript location:** Figure 3 caption
- **Repository files:** `delong_analysis/figures/delong_significance.pdf/png`
- **Evidence:** Single combined bar plot with AUC bars and significance stars.
- **Recommended correction:** Revise caption to describe the actual layout and note gnomAD Common Benign omission.

### MIN-11: Figure 4a y-axis mixes AUC and −FPR
- **Severity:** Minor
- **Manuscript location:** Figure 4 caption
- **Repository files:** `variant_type_analysis/compute_variant_type_analysis.py` `plot_variant_type_auc()`
- **Evidence:** Y-axis label "AUC / -FPR"; gnomAD benign (all-negative) bars plot −FPR, not AUC.
- **Recommended correction:** Clarify in caption or separate gnomAD benign into a different metric panel.

### MIN-12: Figure 5 empty gnomad_benign group unexplained
- **Severity:** Minor
- **Manuscript location:** Figure 5 caption
- **Repository files:** `variant_type_analysis/non_missense_auc_summary.pdf/png`
- **Evidence:** gnomAD benign has no pathogenic non-missense variants, so AUC is undefined and the group is empty.
- **Recommended correction:** Remove or annotate the empty group.

### MIN-13: Figure 7b lacks tissue labels on x-axis
- **Severity:** Minor
- **Manuscript location:** Figure 7 caption
- **Repository files:** `tissue_shap_analysis/cumulative_tissue_contribution.pdf/png`
- **Evidence:** X-axis is "Tissue rank" without labels.
- **Recommended correction:** Add tissue-group labels.

### MIN-14: Figure 1 missing
- **Severity:** Minor (but critical for submission)
- **Manuscript location:** Figure plan
- **Repository files:** `figure_audit.csv`; `manuscript_assembly.md`
- **Evidence:** No `figures/fig1_workflow.pdf` exists.
- **Recommended correction:** Create Figure 1.

### MIN-15: Hard-coded paths and thresholds
- **Severity:** Minor
- **Manuscript location:** Methods §2.10
- **Repository files:** `external_validation/scripts/code_external_validation.py` L12; `build_gnomad_benign_benchmark.py` L29; `feature_labbling.py` L38–40
- **Evidence:** Absolute paths and MIN_AF/MAX_PER_GENE constants that are not documented.
- **Recommended correction:** Move to configuration files or environment variables.

### MIN-16: `build_gnomad_benign_benchmark.py` docstring vs. code threshold mismatch
- **Severity:** Minor
- **Manuscript location:** Methods §2.10
- **Repository files:** `external_validation/scripts/build_gnomad_benign_benchmark.py` L29, L36
- **Evidence:** Docstring says "AF > 1%"; code uses `MIN_AF = 0.001` (0.1%).
- **Recommended correction:** Clarify which threshold is correct and update the other.

---

## 7. False Claims

| # | Claim | Location | Evidence | Correction |
|---|---|---|---|---|
| F-1 | Bootstrap CI from "2,000 stratified resamples" | Methods §2.9 | `regularized_final_model.py` uses 200 unstratified resamples | Correct to 200 unstratified or implement 2,000 stratified |
| F-2 | "LightGBM was allowed to treat the ordinal VEP columns as categorical" | Methods §2.5 | `feature_selection.py` detects 0 categorical columns | Remove categorical claim |
| F-3 | Feature-selection "one-standard-error rule" | Methods §2.5 | Code uses raw fold standard deviation, not SE | Rename or implement true 1-SE rule |
| F-4 | Raw feature matrix "approximately 4,631" columns | Methods §2.3 | Actual raw numeric columns = 4,702 | Correct to 4,702 |
| F-5 | Figure 4a shows "AGVarPred and competitors" | Figure 4 caption | Only AGVarPred is plotted | Revise caption or regenerate figure |
| F-6 | Figure 2 shows "95% bootstrap confidence intervals" as shaded regions | Figure 2 caption | No CIs are plotted | Add CIs or revise caption |
| F-7 | Figure 5 has "95% bootstrap confidence intervals" error bars | Figure 5 caption | No error bars | Add error bars or revise caption |
| F-8 | Top AlphaGenome features exceed "any other single functional-annotation feature" | Results §3.7 | `vep_IMPACT_score` SHAP (0.437) > top AlphaGenome (0.290) | Remove or qualify |
| F-9 | H3K9me3 "appeared in the top 50" SHAP features | Results §3.8 | Lowest H3K9me3 rank = 57 | Remove |
| F-10 | DVD AlphaMissense covered "1.4% of ‘other’ variants" | Results §3.5 | "Other" coverage = 47.0%; overall coverage = 1.7% | Correct |

---

## 8. Unsupported Claims

| # | Claim | Location | Status | Evidence Needed |
|---|---|---|---|---|
| U-1 | "Broad variant coverage" / "broad SNV consequence-class coverage" | Abstract | Unverifiable without fixed variant-type analysis | Recompute internal per-class AUCs after aligning predictions |
| U-2 | "Complementary functional signal" from AlphaGenome | Abstract | Interpretive; supported only if internal variant-type analysis is fixed | Recompute or soften wording |
| U-3 | "AGVarPred produced scores for all variants in every external benchmark file" | Results §3.3 | Partially supported; VIP has ClinVar-derived labels | Re-evaluate after handling VIP circularity |
| U-4 | Clinical utility / deployment readiness | Discussion | Unsupported; no prospective validation | Tone down or remove |

---

## 9. Code–Manuscript Mismatches

| Topic | Manuscript | Code | Severity |
|---|---|---|---|
| Bootstrap CI | 2,000 stratified | 200 unstratified | Major |
| VEP categorical | Treated as categorical | 0 categorical columns | Moderate |
| Feature selection rule | 1-SE | 1-SD | Moderate |
| Raw feature count | ~4,631 | 4,702 | Moderate |
| Final feature panel (markdown) | 111+2+7 | 111+1+8 | Major |
| scale_pos_weight | Per-fold | Global | Minor |
| Threshold grid | 0.05–0.95 | 0.05–0.94 | Minor |
| gnomAD AF missing | Distinct categories | Conflated with AF=0 | Moderate |
| Multi-allelic gnomAD | Correct per-allele handling | First ALT/VEP only | Major |
| Variant-type alignment | By `variant_id` | By row index with mismatched sort order | Critical |
| Benchmark deduplication | Against authoritative splits | Some scripts use parquet parts | Moderate |

---

## 10. Statistical Issues

1. **Non-stratified bootstrap** (Major): Both internal (`regularized_final_model.py`) and external (`scoring_core.py`) bootstraps are unstratified, risking uninformative resamples in small benchmarks.
2. **Feature-selection optimism** (Major): Same folds used for feature-set selection and performance estimation.
3. **VIP label circularity** (Critical): 1,075/3,880 VIP labels are ClinVar-derived.
4. **Internal variant-type AUCs unreliable** (Critical): Misaligned predictions/features invalidate subgroup analysis.
5. **Tissue SHAP bootstrap across features** (Moderate): Unclear frequentist interpretation.
6. **DCA competitor comparisons on different subsets** (Minor): "Treat all" line uses AGVarPred subset prevalence; competitor net-benefit curves use common subsets.

---

## 11. Figure Issues

| Figure | Issue | Severity |
|---|---|---|
| Figure 1 | Missing | High |
| Figure 2 | Caption promises 3-panel layout + bootstrap CIs; actual is 3 separate PNGs with no CIs | High |
| Figure 3 | Caption describes two panels; figure is single composite; gnomAD omitted | Medium |
| Figure 4a | Caption says competitors plotted; only AGVarPred shown; y-axis mixes AUC/−FPR; internal bars misaligned | Critical |
| Figure 4b | OK | — |
| Figure 5 | Caption claims error bars; none plotted; empty gnomad_benign group unexplained | High |
| Figure 6 | OK | — |
| Figure 7b | X-axis lacks tissue labels | Low |

---

## 12. Table Issues

| Table | Issue | Severity |
|---|---|---|
| `internal_test_variant_type_metrics.csv` | Internal metrics unreliable due to misalignment | Critical |
| `variant_type_coverage.csv` | Values correct; DVD "other" percentage misattributed in text | Minor |
| `non_missense_comparison.csv` | Values correct; gnomad_benign group empty but present | Minor |
| `delong_results_fdr.csv` | Verified correct | — |
| `decision_curve_summary.csv` | Verified correct | — |
| `feature_importance_shap.csv` | Verified correct | — |
| `tissue_shap_importance.csv` / `biosample_shap_importance.csv` | Verified correct | — |
| `variant_benchmark_map.csv` | Still contains 35,299 contaminated ClinVar entries | Moderate |

---

## 13. Reproducibility Issues

1. **Execution workspace outside Git** (Moderate): Top-level scripts and `external_validation/scripts/` are not tracked; differ from canonical `AGVarPred/` copies.
2. **Uncommitted changes** (Minor): `git status` shows 6 modified files in the canonical repo.
3. **Dead references** (Minor): `data_manifest.json` references removed scripts.
4. **Hard-coded paths** (Minor): Absolute paths in `code_external_validation.py`, `build_gnomad_benign_benchmark.py`.
5. **No failed-variant manifest** (Major): Silent API drops are not logged in a machine-readable way.
6. **No nested CV for feature selection** (Major): Exact feature-set selection may not reproduce the reported curve.
7. **Figure 1 missing** (High): Workflow figure not in repository.
8. **Top-level duplicate scripts** (Minor): 16 byte-identical copies create confusion.

---

## 14. Reviewer Concerns

1. **Why were 3,790 test variants dropped?** The cause (2,000/gene cap + `<10`/gene filter) must be fully explained or removed.
2. **How can internal per-class AUCs be trusted?** The row-index misalignment means they cannot; they must be recomputed.
3. **Is VIP truly independent?** 1,075 ClinVar-derived labels strongly undermine independence.
4. **Is the bootstrap CI valid?** 200 unstratified resamples is not what was described.
5. **Why are VEP columns treated as numeric despite the categorical claim?** LightGBM handles numeric and categorical features differently.
6. **Can another researcher reproduce the pipeline?** The divergent `external_validation/scripts/` copies and hard-coded paths are red flags.
7. **Is the feature-selection curve optimistically biased?** Same-fold selection/evaluation is a standard concern.
8. **Why are several figure captions inconsistent with the figures?** This suggests insufficient proofreading and risks desk rejection.

---

## 15. Required Corrections Before Submission

### Must-fix (block submission)
1. **Fix internal variant-type alignment:** add `variant_id` to `test_predictions.csv` and merge on `variant_id`; recompute `internal_test_variant_type_metrics.csv` and Figure 4a.
2. **Address VIP circularity:** exclude or disclose the 1,075 ClinVar-derived VIP labels; recompute VIP metrics.
3. **Correct bootstrap CI:** implement the described 2,000 stratified resamples or revise Methods to match the 200 unstratified code.
4. **Fix Figure 4a caption/figure mismatch:** either plot competitor per-class AUCs or revise the caption.
5. **Create Figure 1.**
6. **Reconcile Figure 2 and Figure 5 captions** with actual figure content.

### Should-fix (strongly recommended)
7. Remove or justify the 2,000/gene cap and `<10`/gene filter; produce a discarded-variant manifest.
8. Fix multi-allelic gnomAD handling.
9. Implement nested CV for feature-set size selection.
10. Standardize all benchmark builders to deduplicate against authoritative `train.csv`/`cal.csv`/`test.csv`.
11. Rebuild `variant_benchmark_map.csv` without ClinVar entries.
12. Move execution into the Git-tracked repo and reconcile divergent `external_validation/scripts/` copies.
13. Add retry logic and a failed-variant manifest for AlphaGenome API calls.
14. Update `methods_section.md` to the correct 111+1+8 feature panel.

### Nice-to-fix
15. Distinguish "gnomAD AF missing" from "AF=0".
16. Fix minor count typos in external benchmark descriptions.
17. Correct SHAP percentage claims and H3K9me3 top-50 claim.
18. Add tissue labels to Figure 7b.
19. Remove obsolete `final_model_output/` directory.

---

## 16. Final Verdict

**The manuscript is not scientifically correct in its current form and should not be submitted.**

While the core model training pipeline is sound and most headline metrics are reproducible, the following make submission premature:

- **Critical:** Internal variant-type analysis is broken due to row-order misalignment; the reported per-class AUCs for the internal test set are unreliable.
- **Critical:** VIP benchmark has substantial ClinVar-derived label overlap.
- **Critical/Major:** 23.4% of the test set is silently discarded by undocumented filters.
- **Major:** Bootstrap CI, VEP categorical handling, and feature-selection rule are misdescribed.
- **Major/Medium:** Several figure captions are factually inconsistent with the plotted content.

**I would not sign as co-author confirming the manuscript accurately represents the repository until:**
1. The internal variant-type analysis is recomputed correctly.
2. VIP circularity is resolved or fully disclosed.
3. Methods descriptions match the code exactly.
4. All figure captions match the figures.
5. The missing workflow figure is created.

After these corrections, the manuscript would present a credible, though still limited, computational pathology predictor. External benchmark independence and generalizability would remain the primary scientific concerns for reviewers.

---

*End of forensic audit report.*
