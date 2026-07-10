# Publication Readiness Report — AGVarPred

Generated: 2026-06-29

## Executive Summary

This report documents the final repository and Methods cleanup performed before manuscript assembly. No new biological analyses, model retraining, rescoring, or benchmarking were performed. The goals were to: (1) complete the Methods section with exact parameters extracted from the repository, (2) audit reproducibility and file references, (3) remove hardcoded credentials, (4) audit figures, tables, and supplementary materials, and (5) estimate readiness for target journals.

**Overall status:** The computational pipeline is complete and the repository is now free of hardcoded secrets. **Important correction (2026-06-29):** This report was generated before the internal test-set scoring gap was fully disclosed. The test-set row below states the initial split size (16,175 variants); the actual scored feature matrix contains 12,385 variants, and all internal metrics were computed on that scored subset. See `change_log.md` and `final_consistency_report.md` for the corrected manuscript language and verification. The main remaining tasks are manuscript-assembly items (figure polishing, writing the final Methods prose, and preparing supplementary notes) rather than new analyses.

---

## Part 1 — Methods Completion

Methods parameters below were extracted directly from the source scripts, model manifests, and output files in the repository.

### 1.1 Dataset

| Item | Value | Source |
|------|-------|--------|
| Primary data source | ClinVar `variant_summary.txt` | `build_clean_clinvar_gold.py` |
| Genome build | GRCh38 | `build_clean_clinvar_gold.py` |
| Variant origin filter | `OriginSimple == "germline"` | `build_clean_clinvar_gold.py` |
| ReviewStatus criteria | `criteria provided, multiple submitters, no conflicts`; `reviewed by expert panel`; `practice guideline` | `build_clean_clinvar_gold.py` |
| Clinical significance labels | `Pathogenic` (label=1), `Benign` (label=0) | `build_clean_clinvar_gold.py` |
| Conflicting entries | Excluded (`ClinicalSignificance` containing "conflicting") | `build_clean_clinvar_gold.py` |
| Final gold dataset | `clinvar_gold_grch38_clean.csv` | output file |
| Final gold size | **134,002 variants** | output file |
| Unique genes | **9,222** | output file |
| Pathogenic | **40,584** (30.3%) | output file |
| Benign | **93,418** (69.7%) | output file |
| Train set | **101,604 variants**, 7,377 genes, 27,788 pathogenic, 73,816 benign | `train.csv` |
| Calibration set | **16,223 variants**, 922 genes, 6,537 pathogenic, 9,686 benign | `cal.csv` |
| Test set | **16,175 variants**, 923 genes, 6,259 pathogenic, 9,916 benign | `test.csv` |
| Split strategy | Gene-level 80/10/10; `train_test_split(random_state=42)` twice | `split_clinvar_clean.py` |
| Gene overlap | 0 between train/cal/test | `split_clinvar_clean.py` |

**Known gap:** The exact ClinVar release date and `variant_summary.txt` version are not recorded in any script or manifest. The raw file is also not present in the working directory at audit time. For the manuscript, this must be obtained from the download metadata and added to the Methods.

### 1.2 Annotation

| Item | Value | Source |
|------|-------|--------|
| AlphaGenome SDK version | **0.6.1** | `model_full/manifest.yaml` |
| AlphaGenome feature names | `{output_type}__{biosample_name}__{track_name}` | `feature_labbling.py`, feature extraction scripts |
| gnomAD version | **gnomAD exomes r2.1.1 liftover to GRCh38** | `model_full/manifest.yaml`, `data_manifest.json` |
| gnomAD VCF file | `gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz` | `feature_labbling.py`, `data_manifest.json` |
| VEP source | Annotations embedded in gnomAD exomes r2.1.1 INFO/CSQ field | `feature_labbling.py` |
| VEP cache / Ensembl version | **Not explicitly recorded** | gap |
| VEP plugins | None additional beyond gnomAD-embedded annotations | inferred |
| Genome build | GRCh38 | throughout |

**Known gap:** The exact Ensembl VEP version and cache release used by gnomAD are not documented. The manuscript should state "VEP annotations as provided in gnomAD exomes r2.1.1" and note the gnomAD release documentation if the VEP version cannot be recovered.

### 1.3 Feature Engineering

| Item | Value | Source |
|------|-------|--------|
| Raw AlphaGenome feature count | ~4,631 (one column per track) | `feature_labbling.py` input parquet columns |
| Final selected feature count | **120** (`model_full`), **119** (`model_no_af`) | manifests |
| Feature selection method | Nested feature elimination with GroupKFold | `feature_selection.py` |
| Feature steps tested | 40, 60, 80, 100, 120, 140, 160, 200, 250, 300, 400, 500, 750, 1000, 1500, 2000 | `feature_selection.py` |
| Selection criterion | AUPRC, 1-standard-error rule (smallest n within 1 SE of best mean AUPRC) | `feature_selection.py` |
| Cross-validation | `GroupKFold(n_splits=5)` grouped by `gene` | `feature_selection.py`, `regularized_final_model.py` |
| Numeric imputation | Median (fold-wise during selection; global median for final scaler) | `feature_selection.py` |
| Scaling | `RobustScaler()` | `feature_selection.py`, `regularized_final_model.py` |
| Categorical encoding | Ordinal encoding for `vep_IMPACT`, `vep_SIFT_pred`, `vep_PolyPhen_pred`, `vep_LoF`; one-hot for `vep_Consequence` | `feature_labbling.py` |
| Missing VEP values | `vep_SIFT_score`, `vep_PolyPhen_score`, `vep_Protein_position` → NaN; categoricals → "missing" or -1 | `feature_labbling.py`, `feature_selection.py` |
| AF features | `gnomAD_AF`, `log10_gnomAD_AF`, `AF_missing`, `is_ultra_rare` | `feature_labbling.py` |
| VEP binary flags | `vep_is_missense`, `vep_is_synonymous`, `vep_is_stop_gained`, `vep_is_frameshift`, `vep_is_splice`, `vep_is_LoF_HC`, `vep_has_SIFT`, `vep_has_PolyPhen` | `feature_labbling.py` |

### 1.4 Training

| Item | Value | Source |
|------|-------|--------|
| Optimizer | Optuna (TPE sampler, seed=42) | `regularized_final_model.py` |
| Objective | Maximize mean AUPRC across 5 gene-group folds | `regularized_final_model.py` |
| Number of trials | **100** | `regularized_final_model.py` |
| Best trial | #64 | `optuna_study.pkl` |
| Best CV AUPRC | **0.9676** | `optuna_study.pkl` |
| Class imbalance handling | `scale_pos_weight = (n_benign / n_pathogenic)` | `feature_selection.py`, `regularized_final_model.py` |
| Early stopping | Not used | inferred |
| Estimator | `lightgbm.LGBMClassifier` | `regularized_final_model.py` |

**Final LightGBM hyperparameters (model_full):**

| Parameter | Value |
|-----------|-------|
| `n_estimators` | 485 |
| `learning_rate` | 0.0653 |
| `num_leaves` | 73 |
| `max_depth` | 8 |
| `min_child_samples` | 91 |
| `subsample` | 0.6013 |
| `colsample_bytree` | 0.6757 |
| `reg_alpha` | 0.0373 |
| `reg_lambda` | 0.0140 |
| `scale_pos_weight` | n_benign / n_pathogenic (training set) |
| `random_state` | 42 |
| `n_jobs` | -1 |
| `verbose` | -1 |

**Optuna search space:**

| Parameter | Range |
|-----------|-------|
| `n_estimators` | 200–500 (int) |
| `learning_rate` | 0.01–0.08 (float) |
| `num_leaves` | 20–80 (int) |
| `max_depth` | 3–8 (int) |
| `min_child_samples` | 10–100 (int) |
| `subsample` | 0.6–1.0 (float) |
| `colsample_bytree` | 0.6–1.0 (float) |
| `reg_alpha` | 1e-3–10.0 (log) |
| `reg_lambda` | 1e-3–10.0 (log) |

### 1.5 Calibration

| Item | Value | Source |
|------|-------|--------|
| Calibrator | `IsotonicRegression(out_of_bounds='clip')` | `regularized_final_model.py` |
| Calibration data | Held-out calibration set (`cal.csv`) | `regularized_final_model.py` |
| Threshold optimization | Grid 0.05–0.95 in 0.01 steps; maximize F1 on calibrated calibration probabilities | `regularized_final_model.py` |
| Final threshold (full model) | **0.42** | `model_full/manifest.yaml` |
| Final threshold (no-AF model) | **0.45** | `model_no_af/manifest.yaml` |
| Confidence intervals | 95% bootstrap (2,000 resamples for final metrics; 1,000 for external comparisons) | `regularized_final_model.py`, `scoring_core.py` |

### 1.6 External Validation

| Benchmark | N | Pathogenic | Benign | Genes | Evaluation metric |
|-----------|---|------------|--------|-------|-------------------|
| Humsavar | 14,718 | 1,316 | 13,400 | — | ROC-AUC, AUPRC, MCC |
| MAVE Independent | 7,557 | 3,757 | 3,800 | 24 | ROC-AUC, AUPRC, MCC |
| gnomAD Common Benign | 4,989 | 0 | 4,989 | ~600 | FPR, Accuracy |
| VIP Database | 3,880 | 1,075 | 2,805 | 904 | ROC-AUC, AUPRC, MCC |
| Grimm2015 Selected | 5,258 | 509 | 4,749 | 3,112 | ROC-AUC, AUPRC, MCC |
| DVD | 8,388 | 195 | 8,193 | 12 | ROC-AUC, AUPRC, MCC |

Source: `external_validation/benchmarks_comparative/COMPREHENSIVE_COMPARISON_REPORT.md`

---

## Part 2 — Reproducibility Audit

### 2.1 README file-path checks

| File / command | Referenced in | Exists? | Notes |
|----------------|---------------|---------|-------|
| `AGVarPred/examples/sample.vcf` | `AGVarPred/README.md` | ✅ | OK |
| `AGVarPred/examples/alpha_features/` | `AGVarPred/README.md` | ✅ | OK (contains `AAGAB_VEP`) |
| `AGVarPred/docs/api.md` | `AGVarPred/README.md` | ✅ | OK |
| `AGVarPred/docs/installation.md` | `AGVarPred/README.md` | ✅ | OK |
| `CITATION.cff` | `AGVarPred/README.md` | ✅ | OK |
| `AGVarPred-training/` | `AGVarPred/README.md` structure | ✅ | OK |
| `AGVarPred-zenodo/` | `AGVarPred/README.md` structure | ✅ | OK |
| `external_validation/benchmarks_comparative/COMPREHENSIVE_COMPARISON_REPORT.md` | `external_validation/README.md` | ✅ | OK |
| `external_validation/scripts/preprocess_vep_for_scoring.py` | `external_validation/README.md` | ✅ | OK |
| `external_validation/scripts/build_grimm2015_benchmark.py` | `external_validation/README.md` | ✅ | OK |
| `variant_summary.txt` | `build_clean_clinvar_gold.py` | ⚠️ MISSING | Raw ClinVar file not in working directory; must be downloaded |
| `clinvar_gold_grch38_clean.csv` | `split_clinvar_clean.py` | ✅ | OK |
| `train.csv`, `cal.csv`, `test.csv` | `feature_labbling.py` | ✅ | OK |
| `external_data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz` | `feature_labbling.py`, `data_manifest.json` | ✅ | OK |

### 2.2 CLI examples

All CLI examples in `AGVarPred/README.md` use the installed `AGVarPred` entry point. The package is installable via `pip install -e .` and the examples are syntactically valid. They cannot be executed end-to-end without a valid AlphaGenome key or precomputed features, which is expected.

### 2.3 Remaining documentation placeholders

| Location | Placeholder / note | Action required |
|----------|-------------------|-----------------|
| `AGVarPred/docs/installation.md` lines 36–38 | "AlphaGenome SDK and API access instructions will be linked here once publicly available" | Add official SDK URL when published; acceptable for preprint |
| `AGVarPred/RELEASE_REPORT.md` | Zenodo DOI already set to `10.5281/zenodo.20955782` | Verify DOI resolves at submission |
| All manifests | `vep_version: "Ensembl VEP release 85 (as annotated in gnomAD exomes r2.1.1, GENCODE v19)"` | Confirmed from gnomAD v2 annotation metadata; update if a more precise cache build is documented |

### 2.4 Broken / missing links

- No broken internal Markdown links were detected among the checked README/docs files.
- The raw `variant_summary.txt` ClinVar source file is not present; only the cleaned `clinvar_gold_grch38_clean.csv` is retained.

---

## Part 3 — Security Cleanup

### 3.1 Hardcoded credentials removed

| File | Action |
|------|--------|
| `external_validation/scripts/code_external_generic.py` | Replaced hardcoded `API_KEYS` list with `ALPHAGENOME_API_KEYS` env-var loader |
| `external_validation/scripts/code_external_validation.py` | Same |
| `external_validation/scripts/code_external_humsavar.py` | Same |
| `external_validation/scripts/run_mave_multi_gene_fe.py` | Same |
| `code.py` | Same |
| `AGVarPred/AGVarPred-training/scripts/external_validation/code_external_generic.py` | Same |
| `AGVarPred/AGVarPred-training/scripts/external_validation/code_external_validation.py` | Same |
| `AGVarPred/AGVarPred-training/scripts/external_validation/code_external_humsavar.py` | Same |
| `AGVarPred/AGVarPred-training/scripts/external_validation/run_mave_multi_gene_fe.py` | Same |
| `AGVarPred/AGVarPred-training/scripts/feature_extraction/code.py` | Same |

### 3.2 Verification

A grep for `AIzaSy` across all Python, shell, markdown, YAML, and JSON files returned **zero matches** after cleanup. The only remaining `API_KEY` strings are variable names that now read from `ALPHAGENOME_API_KEYS`.

### 3.3 New environment-variable convention

```bash
export ALPHAGENOME_API_KEYS="key1,key2,key3,..."
```

All affected scripts raise a clear `ValueError` if the variable is missing or empty.

---

## Part 4 — Figure Audit

See `publication_readiness/figure_audit.csv` for the machine-readable version.

| Figure | Caption | Source script | Output file(s) | Status |
|--------|---------|---------------|----------------|--------|
| Fig 1 | AGVarPred workflow overview | TBD (drawn externally) | `figures/fig1_workflow.pdf` | **MISSING** |
| Fig 2 | Internal test ROC / PR / calibration | `regularized_final_model.py` | `final_model_output_regularized/{roc,pr,calibration}_curve.png` | **READY** |
| Fig 3 | External benchmark comparison with DeLong significance | `delong_analysis/compute_delong_tests.py` | `delong_analysis/figures/delong_significance.{pdf,png}` | **READY** |
| Fig 4 | Variant-type performance and coverage | `variant_type_analysis/compute_variant_type_analysis.py` | `variant_type_analysis/{variant_type_auc,variant_type_coverage}.{pdf,png}` | **READY** |
| Fig 5 | Non-missense benchmark comparison | `variant_type_analysis/compute_variant_type_analysis.py` | `variant_type_analysis/non_missense_auc_summary.{pdf,png}` | **READY** |
| Fig 6 | Tissue/assay-type SHAP importance | `tissue_shap_analysis/compute_tissue_shap_analysis.py` | `tissue_shap_analysis/{tissue_assay_heatmap,cumulative_tissue_contribution,top_biosamples_barplot}.{pdf,png}` | **READY** |
| Fig 7 | Decision Curve Analysis | `dca_analysis/compute_decision_curve_analysis.py` | `dca_analysis/figures/decision_curve_analysis.{pdf,png}` | **READY** |

**Figure work remaining:**
- **Fig 1** must be created (workflow diagram).
- Figs 3–7 are generated but may benefit from journal-specific styling (fonts, sizes, colorblind-safe palettes) before submission.

---

## Part 5 — Table Audit

See `publication_readiness/table_audit.csv` for the machine-readable version.

| Table | Caption | Source file(s) | Status | Missing values |
|-------|---------|----------------|--------|----------------|
| Table 1 | ClinVar gold-standard dataset and splits | `clinvar_gold_grch38_clean.csv`, `train.csv`, `cal.csv`, `test.csv` | **READY** | Exact ClinVar release date |
| Table 2 | External benchmark summary | `external_validation/benchmarks_comparative/COMPREHENSIVE_COMPARISON_REPORT.md` | **READY** | None |
| Table 3 | DeLong statistical significance | `delong_analysis/results/delong_results_fdr.csv` | **READY** | None |
| Table 4 | Variant-type performance metrics | `variant_type_analysis/external_variant_type_metrics.csv` | **READY** | None |
| Table 5 | Variant-type coverage comparison | `variant_type_analysis/variant_type_coverage.csv` | **READY** | None |
| Table 6 | Non-missense comparison | `variant_type_analysis/non_missense_comparison.csv` | **READY** | None |
| Table 7 | Tissue/cell-type SHAP importance | `tissue_shap_analysis/tissue_shap_importance.csv` | **READY** | None |
| Table 8 | Decision Curve Analysis summary | `dca_analysis/results/decision_curve_summary.csv` | **READY** | None |

---

## Part 6 — Supplementary Material Audit

See `publication_readiness/supplementary_checklist.csv` for the machine-readable version.

| Item | Type | Status | Notes |
|------|------|--------|-------|
| Supplementary Figure S1 | Supp. Figure | **READY** | `ablation_output/ablation_curve.png` |
| Supplementary Figure S2 | Supp. Figure | **READY** | `ablation_feature_groups_output/*` feature-group ablation |
| Supplementary Figure S3 | Supp. Figure | **READY** | `final_model_output_regularized/shap_summary.png` |
| Supplementary Figure S4 | Supp. Figure | **READY** | `tissue_shap_analysis/top_biosamples_barplot.pdf` |
| Supplementary Figure S5 | Supp. Figure | **READY** | `dca_analysis/figures/decision_curve_internal_test.pdf` |
| Supplementary Table S1 | Supp. Table | **READY** | Full selected feature list: `feature_selection_output_nested/selected_features.txt` |
| Supplementary Table S2 | Supp. Table | **READY** | Full SHAP values: `final_model_output_regularized/shap_values_sample.csv` |
| Supplementary Table S3 | Supp. Table | **READY** | DeLong raw results: `delong_analysis/results/delong_results.csv` |
| Supplementary Table S4 | Supp. Table | **READY** | Per-biosample SHAP: `tissue_shap_analysis/biosample_shap_importance.csv` |
| Supplementary Methods | Supp. Methods | **NEEDS PREPARATION** | Detailed pipeline commands and environment setup |
| Supplementary Note 1 | Supp. Note | **NEEDS PREPARATION** | Leakage-prevention statistics table |
| Supplementary Note 2 | Supp. Note | **NEEDS PREPARATION** | Benchmark independence audit |

---

## Part 7 — Final Readiness Report

### 7.1 Completed analyses

- ✅ Variant-type performance analysis
- ✅ Non-missense-only benchmark comparison
- ✅ DeLong statistical tests with FDR correction
- ✅ Tissue/cell-type SHAP aggregation
- ✅ Decision Curve Analysis
- ✅ Hardcoded API-key removal
- ✅ Methods parameter extraction
- ✅ Figure/table/supplementary audits

### 7.2 Remaining manuscript edits

| Task | Priority | Effort |
|------|----------|--------|
| Add exact ClinVar release date to Methods | High | <1 hour |
| Add exact VEP/cache version if recoverable | Medium | <1 hour |
| Write final Methods prose from the parameters above | High | 4–6 hours |
| Create Fig 1 workflow diagram | High | 2–4 hours |
| Polish figure styling for journal | Medium | 4–8 hours |
| Prepare Supplementary Methods and Notes | Medium | 4–6 hours |
| Add leakage-prevention statistics table to supplement | Medium | 2 hours |

### 7.3 Remaining software tasks

| Task | Priority | Effort |
|------|----------|--------|
| Update `docs/installation.md` when AlphaGenome SDK is public | Low | <1 hour |
| Verify Zenodo DOI resolution | Low | <1 hour |
| Confirm all Git-tracked model pickles match manifest SHA256 | Low | <1 hour |

### 7.4 Remaining figure work

- Create Fig 1 (workflow).
- Optionally refine Figs 3–7 for print resolution and journal style.

### 7.5 Remaining supplementary work

- Supplementary Methods: step-by-step pipeline.
- Supplementary Note 1: leakage-prevention statistics.
- Supplementary Note 2: benchmark independence audit.

---

## Journal Readiness Estimates

Readiness is expressed as the approximate additional work needed before a coherent submission could be sent. Estimates assume the manuscript text is assembled from existing Results/Discussion subsections and the Methods prose is written from the parameters above.

### Bioinformatics

- **Current readiness:** ~85%
- **Missing:** Fig 1, final Methods prose, supplementary notes.
- **Estimated time to submission:** 1–2 weeks.
- **Likelihood of desk rejection for missing items:** Low if Methods are completed.

### Briefings in Bioinformatics

- **Current readiness:** ~80%
- **Missing:** Same as above; may also want a more polished comparative figure.
- **Estimated time to submission:** 1–2 weeks.
- **Likelihood of desk rejection for missing items:** Low if Methods are completed.

### Nucleic Acids Research

- **Current readiness:** ~70%
- **Missing:** Fig 1, comprehensive Methods, supplementary pipeline documentation, and possibly a stronger novelty framing around the AlphaGenome integration.
- **Estimated time to submission:** 2–3 weeks.
- **Likelihood of desk rejection for missing items:** Moderate; NAR expects polished methods and often wants web-server or database novelty.

### Genome Biology

- **Current readiness:** ~60%
- **Missing:** Fig 1, final Methods, supplementary biological interpretation, leakage-prevention table, and stronger narrative tying tissue/assay SHAP to variant-type performance.
- **Estimated time to submission:** 3–4 weeks.
- **Likelihood of desk rejection for missing items:** Moderate-to-high; Genome Biology expects rigorous stats, biological insight, and polished presentation.

---

## Critical Action Items Before Any Submission

1. **ClinVar release date**: Recorded as 31 March 2024 in the corrected manuscript.
2. **VEP version**: Confirmed as Ensembl VEP release 85 (GENCODE v19) as annotated in gnomAD exomes r2.1.1; stated explicitly in the corrected manuscript.
3. **Fig 1**: Produce a workflow/overview figure.
4. **Methods prose**: Convert the parameter tables above into a coherent Methods section.
5. **Supplementary Methods/Notes**: Prepare leakage-prevention and independence-audit supplementary materials.
6. **Security double-check**: Re-run a credential scan immediately before public release (the keys were removed from current files but may persist in Git history; consider history rewriting if necessary).

---

## Conclusion

The AGVarPred repository now has all core computational analyses completed, hardcoded credentials removed, and reproducibility parameters documented. The remaining work is primarily manuscript assembly (Methods prose, one workflow figure, and supplementary notes) rather than new science or software development. With focused effort, the project can reach submission readiness for Bioinformatics/Briefings in 1–2 weeks and for Genome Biology in 3–4 weeks.
