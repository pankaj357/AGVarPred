# Final Consistency Report — AGVarPred Manuscript

**Date:** 2026-06-29  
**Repository:** `/data/kanaka/pankaj/feature_extration/new`  
**Scope:** `publication_readiness/methods_section.md`, `publication_readiness/results_section.md`, `publication_readiness/discussion_section.md`, `publication_readiness/manuscript_assembly.md`, and `publication_readiness/change_log.md`.

---

## 1. Purpose

This report verifies that the corrected AGVarPred manuscript sections contain only claims supported by repository outputs. It confirms that the previously identified unsupported statements—(i) that the internal test set contained 16,175 evaluated variants, and (ii) that only single-nucleotide variants were retained—have been removed or corrected, and that no new unsupported claims were introduced.

---

## 2. Verification methodology

- Each numerical claim in the corrected manuscript sections was matched to a repository output file, CSV, Parquet dataset, or previously verified manifest.
- The `build_clean_clinvar_gold.py` script was inspected to confirm that no variant-type filter is applied.
- Parquet row counts were obtained from file metadata without loading full datasets.
- The final model predictions file was checked against the scored test-set size.
- Language about coverage, variant types, and evaluation subsets was reviewed for consistency across sections.

---

## 3. Verified numerical claims

| Claim in manuscript | Reported value | Source | Status |
|---------------------|----------------|--------|--------|
| ClinVar gold-standard dataset size | 134,002 variants | `clinvar_gold_grch38_clean.csv` | ✅ Verified |
| Gold-standard pathogenic prevalence | 40,584 / 134,002 (30.3%) | `clinvar_gold_grch38_clean.csv` | ✅ Verified |
| Gold-standard unique genes | 9,222 | `clinvar_gold_grch38_clean.csv` | ✅ Verified |
| Initial training split | 101,604 variants, 7,377 genes | `train.csv` | ✅ Verified |
| Initial calibration split | 16,223 variants, 922 genes | `cal.csv` | ✅ Verified |
| Initial test split | 16,175 variants, 923 genes | `test.csv` | ✅ Verified |
| Scored training subset | 86,378 variants | `final_dataset_parts_train/*.parquet` | ✅ Verified |
| Scored calibration subset | 13,201 variants | `final_dataset_parts_cal/*.parquet` | ✅ Verified |
| Scored test subset | 12,385 variants | `final_dataset_parts_test/*.parquet` | ✅ Verified |
| Test predictions file | 12,385 variants + header | `final_model_output_regularized/test_predictions.csv` | ✅ Verified |
| Scored test prevalence | 37.0% pathogenic (4,576 / 12,385) | `final_model_output_regularized/test_predictions.csv` | ✅ Verified |
| Scored calibration prevalence | 41.8% pathogenic (5,522 / 13,201) | `final_dataset_parts_cal/*.parquet` | ✅ Verified |
| Internal ROC-AUC | 0.9753 (95% CI 0.9730–0.9774) | `final_model_output_regularized/final_metrics.csv` | ✅ Verified |
| Internal PR-AUC | 0.9606 | `final_model_output_regularized/final_metrics.csv` | ✅ Verified |
| Internal F1 | 0.8945 | `final_model_output_regularized/final_metrics.csv` | ✅ Verified |
| Internal MCC | 0.8327 | `final_model_output_regularized/final_metrics.csv` | ✅ Verified |
| Internal Brier score | 0.0568 | `final_model_output_regularized/final_metrics.csv` | ✅ Verified |
| Internal ECE | 0.0225 | `final_model_output_regularized/final_metrics.csv` | ✅ Verified |
| Selected feature count | 120 (full model), 119 (no-AF model) | `feature_selection_output_nested/selected_features.txt`, model manifests | ✅ Verified |
| Best Optuna trial / CV AUPRC | Trial 64 / 0.9676 | `final_model_output_regularized/optuna_study.pkl` | ✅ Verified |
| Final threshold | 0.42 (full model), 0.45 (no-AF model) | Model manifests | ✅ Verified |

---

## 4. Corrected unsupported claims

The following claims were present in earlier drafts and have been removed or corrected:

| Unsupported claim | Correction | Status |
|-------------------|------------|--------|
| "The internal test set contained 16,175 variants" (implying all were evaluated). | The manuscript now distinguishes the 16,175-variant initial split from the 12,385-variant scored subset and states that all internal metrics are computed on the scored subset. | ✅ Corrected |
| "Only single-nucleotide variants were retained." | Removed. The Methods now state that the ClinVar dataset contained multiple variant types and that model development was performed on the subset represented in the final feature matrices. | ✅ Corrected |
| Unqualified "genome-wide coverage." | Replaced with "coverage of evaluated variants" and related qualifying language. | ✅ Corrected |
| "The calibrator / threshold was selected on the 16,223-variant calibration set" (without distinguishing scored subset). | The Methods now state that the threshold was selected on the 13,201-variant scored calibration subset. | ✅ Corrected |
| No mention of the test-set scoring gap. | A new Discussion limitation explicitly states that the final evaluation dataset differed from the initial split and that the cause is unknown. | ✅ Corrected |

---

## 5. Consistency across sections

| Topic | Methods | Results | Discussion | Status |
|-------|---------|---------|------------|--------|
| Test-set size | "12,385 test-set variants represented in the final feature matrices" | "initial split contained 16,175 variants ... 12,385 were represented and evaluated" | "initial test split contained 16,175 variants, but internal evaluation was performed on 12,385" | ✅ Consistent |
| Calibration-set size | "16,223 variants in 922 genes ... scored calibration set (13,201 variants)" | — | — | ✅ Consistent |
| Variant types | "dataset contained multiple variant types" | — | "ClinVar gold-standard dataset ... contained multiple variant types" | ✅ Consistent |
| Coverage framing | "coverage of evaluated variants" | "100% coverage of AGVarPred on the evaluated benchmarks" | "coverage of evaluated variants" | ✅ Consistent |
| Cause of scoring gap | — | "3,790 ... were not represented" | "exact source of this discrepancy is currently unknown" | ✅ Consistent |

---

## 6. Remaining items that are not unsupported claims

The following items are explicitly acknowledged as limitations, gaps, or future work. They do not represent unsupported claims.

| Item | Where noted | Status |
|------|-------------|--------|
| Exact ClinVar release date and `variant_summary.txt` version not archived. | Methods §2, Limitations | ✅ Disclosed |
| Exact Ensembl VEP version used by gnomAD not archived. | Methods §4, Limitations | ✅ Disclosed |
| AlphaGenome platform is not publicly available. | Methods §4, Limitations | ✅ Disclosed |
| Cause of the 3,790-variant test-set scoring gap is unknown. | Limitations | ✅ Disclosed |
| Figure 1 (workflow overview) is missing. | `manuscript_assembly.md` | ✅ Listed as task |
| Supplementary Methods and Notes 1–2 are placeholders. | `manuscript_assembly.md` | ✅ Listed as task |
| No figure or table citations in current text. | `manuscript_assembly.md` | ✅ Listed as task |
| Potential secrets in Git history. | Limitations | ✅ Disclosed |

---

## 7. Claims that remain unverified or require future confirmation

No remaining claims in the corrected manuscript sections were found to be unsupported by the repository outputs reviewed. The following are not claims but require action before submission:

1. Insert figure and table references.
2. Generate Figure 1.
3. Complete supplementary materials.
4. Recover or explicitly acknowledge the ClinVar release date and Ensembl VEP version.
5. Investigate and explain the 3,790-variant test-set scoring gap.
6. Purge any secrets from Git history before public release.

---

## 8. Conclusion

The corrected manuscript sections are internally consistent and contain no unsupported numerical or factual claims. The critical test-set size inconsistency has been disclosed, the false SNV-only filtering claim has been removed, coverage language has been qualified, and a new limitation about the unknown cause of the scoring gap has been added. The manuscript is now ready for the remaining assembly tasks (figure generation, supplementary materials, references, and copy-editing) before journal submission.

---

## 9. Sign-off

| Role | Finding |
|------|---------|
| Numerical claims | All verified against repository outputs |
| Internal consistency | No contradictions detected |
| Unsupported claims | None remaining in corrected sections |
| Limitations / gaps | All explicitly disclosed |
| Submission readiness | Cleared for assembly tasks; not yet ready for submission because Figure 1 and supplementary materials are incomplete |
