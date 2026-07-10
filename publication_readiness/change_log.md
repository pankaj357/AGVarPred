# Manuscript Change Log — AGVarPred Correction Package

**Date:** 2026-06-29  
**Repository:** `/data/kanaka/pankaj/feature_extration/new`  
**Purpose:** Document every factual correction applied to the AGVarPred manuscript after discovering that the internal test-set evaluation was performed on 12,385 scored variants, not the full 16,175-variant ClinVar test split, and that the ClinVar gold-standard dataset contained multiple variant types rather than only SNVs.

---

## Summary of corrections

1. **Test-set evaluation size.** The manuscript now states that internal evaluation was performed on the 12,385 test-set variants represented in the final feature matrices, not on the full 16,175-variant initial ClinVar test split.
2. **Variant-type composition.** The false statement that "Only single-nucleotide variants were retained" has been removed. The Methods now state that the ClinVar gold-standard dataset contained multiple variant types, and that model development and evaluation were restricted to variants represented in the final feature matrices.
3. **Coverage framing.** Unqualified claims of "genome-wide coverage" have been replaced with "coverage of evaluated variants" and related qualifying language.
4. **Calibration-set size.** The calibration-set size is now reported accurately as 16,223 variants in the initial split and 13,201 variants in the scored feature matrix.
5. **New limitation.** A Discussion limitation explicitly notes that the final evaluation dataset differed from the initial ClinVar split and that the cause of the discrepancy is currently unknown.

---

## Detailed change log

### 1. Internal test-set evaluation size

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/methods_section.md`, §8 "Internal Validation" |
| **Original claim** | "Internal performance was evaluated on the held-out ClinVar test set of 16,175 variants." |
| **Corrected claim** | "Internal performance was evaluated on the 12,385 test-set variants represented in the final feature matrices." |
| **Reason / evidence** | `test.csv` contains 16,175 variants, but `final_dataset_parts_test/*.parquet` and `final_model_output_regularized/test_predictions.csv` contain only 12,385 scored variants. All internal metrics were computed on the scored subset. |

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/results_section.md`, §1 "Strong discrimination on a gene-held-out ClinVar test set" |
| **Original claim** | "We first evaluated the frozen AGVarPred pipeline on the held-out ClinVar test set. The held-out ClinVar test set contained 16,175 variants ... On this test set, AGVarPred achieved ..." |
| **Corrected claim** | "We first evaluated the frozen AGVarPred pipeline on the held-out ClinVar test set. The initial gene-level test split contained 16,175 variants ... After feature extraction, 12,385 of these variants were represented in the final feature matrices and were evaluated. ... On this scored test set, AGVarPred achieved ..." |
| **Reason / evidence** | Same as above. The 16,175 figure describes the initial split; the 12,385 figure describes the scored subset on which metrics were computed. |

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/results_section.md`, §3 "Variant coverage distinguishes AGVarPred from existing predictors" |
| **Original claim** | "On the internal ClinVar test set, AGVarPred scored all 16,175 variants." |
| **Corrected claim** | "On the internal ClinVar test set, 12,385 of 16,175 split variants were represented in the final feature matrices and scored." |
| **Reason / evidence** | Same as above. |

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/results_section.md`, §3 "Variant coverage distinguishes AGVarPred from existing predictors" |
| **Original claim** | "AGVarPred scored 100% of evaluated variants in every external benchmark and every molecular-consequence class evaluated." (combined with internal test-set language) |
| **Corrected claim** | "AGVarPred scored 100% of evaluated variants in every external benchmark and every molecular-consequence class evaluated. On the internal ClinVar test set, 12,385 of 16,175 split variants were represented in the final feature matrices and scored." |
| **Reason / evidence** | 100% coverage holds for the external benchmarks evaluated, but the internal test set had a 3,790-variant gap that must be disclosed. |

### 2. Variant-type composition and SNV-only claim

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/methods_section.md`, §2 "Dataset Construction" |
| **Original claim** | "Only single-nucleotide variants were retained." |
| **Corrected claim** | "The resulting ClinVar gold-standard dataset comprised 134,002 unique variants ... The dataset contained multiple variant types; SNVs were the majority, but deletions, duplications, insertions, indels, microsatellites, and other variant classes were also present. Feature extraction, model training, calibration, and internal evaluation were performed only on the variants successfully represented in the final feature matrices ..." |
| **Reason / evidence** | `build_clean_clinvar_gold.py` does not filter by variant type. `clinvar_gold_grch38_clean.csv` contains 112,934 SNVs and 21,068 non-SNV variants. |

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/discussion_section.md`, "Limitations" |
| **Original claim** | (No explicit limitation about the restriction to scored feature-matrix variants.) |
| **Corrected claim** | Added limitation: "Restriction to variants represented in the final feature matrices. AGVarPred is designed for SNVs, and the ClinVar gold-standard dataset used here contained multiple variant types. Model training, calibration, and internal evaluation were performed only on the variants represented in the final feature matrices ... The initial ClinVar splits were larger ... Claims of 'genome-wide coverage' should therefore be understood as coverage of evaluated variants." |
| **Reason / evidence** | Same as above. The previous framing omitted the fact that the ClinVar dataset was not SNV-only. |

### 3. Coverage framing

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/manuscript_assembly.md`, "Key revisions from the previous draft" |
| **Original claim** | "Genome-wide coverage has been reframed as genome-wide SNV coverage throughout." |
| **Corrected claim** | "Coverage claims have been reframed as coverage of evaluated variants rather than genome-wide coverage." |
| **Reason / evidence** | Because the internal test set did not score all split variants, even "genome-wide SNV coverage" would overstate internal coverage. The narrower phrasing "coverage of evaluated variants" is now used. |

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/discussion_section.md`, "Limitations" |
| **Original claim** | "Claims of 'genome-wide coverage' ..." (existing limitation) |
| **Corrected claim** | "Claims of 'genome-wide coverage' should therefore be understood as coverage of evaluated variants." |
| **Reason / evidence** | Same as above. |

### 4. Calibration-set size and scored subset

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/methods_section.md`, §2 "Dataset Construction" |
| **Original claim** | "Model development and internal evaluation were performed on the subset of these variants that were represented in the final feature matrices" (without explicit numbers). |
| **Corrected claim** | "Feature extraction, model training, calibration, and internal evaluation were performed only on the variants successfully represented in the final feature matrices: 86,378 training-set variants, 13,201 calibration-set variants, and 12,385 test-set variants." |
| **Reason / evidence** | Parquet row counts: `final_dataset_parts_train` = 86,378; `final_dataset_parts_cal` = 13,201; `final_dataset_parts_test` = 12,385. |

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/methods_section.md`, §7 "Probability Calibration" |
| **Original claim** | "The calibrator was fit exclusively on the held-out calibration set (16,223 variants in 922 genes not seen during training)." |
| **Corrected claim** | "The calibrator was fit exclusively on the held-out calibration set (16,223 variants in 922 genes not seen during training) by regressing the raw probabilities against the observed binary labels. ... The binary operating threshold was selected on the calibrated probabilities from the scored calibration set (13,201 variants represented in the final feature matrices, 41.8% pathogenic)." |
| **Reason / evidence** | `cal.csv` has 16,223 variants, but the scored calibration matrix has 13,201 variants. The threshold was optimized on the scored subset. |

### 5. New limitation about the final evaluation dataset

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/discussion_section.md`, "Limitations" |
| **Original claim** | (No explicit statement about the test-set scoring gap.) |
| **Corrected claim** | Added limitation: "Final evaluation dataset. The final evaluation dataset differed from the initial ClinVar split because only variants represented in the final feature matrices were evaluated. The initial test split contained 16,175 variants, but internal evaluation was performed on the 12,385 variants present in the scored feature matrix. The exact source of this discrepancy is currently unknown and should be clarified in future releases." |
| **Reason / evidence** | The 3,790 missing variants have not been investigated. Only verified facts are reported; no speculation is offered about the cause. |

### 6. Results-section coverage caveat

| Item | Detail |
|------|--------|
| **Location** | `publication_readiness/results_section.md`, §3 "Variant coverage distinguishes AGVarPred from existing predictors" |
| **Original claim** | "This coverage difference has practical implications ... the 100% coverage of AGVarPred on the evaluated benchmarks represents a meaningful advantage ..." |
| **Corrected claim** | Same sentence, with added caveat: "... This advantage must be weighed against the higher accuracy of specialized tools on missense variants, as shown in the Humsavar comparisons, and against the fact that 3,790 of 16,175 test-split variants were not represented in the internal test-set feature matrices." |
| **Reason / evidence** | Ensures readers understand that the internal test-set coverage was not 100%. |

---

## Files modified in this correction package

- `publication_readiness/methods_section.md`
- `publication_readiness/results_section.md`
- `publication_readiness/discussion_section.md`
- `publication_readiness/manuscript_assembly.md`

## New files created in this correction package

- `publication_readiness/change_log.md` (this file)
- `publication_readiness/final_consistency_report.md`

---

## Verification commands

The counts below were confirmed directly from repository outputs:

```bash
# CSV split sizes (including header)
wc -l train.csv cal.csv test.csv
# -> 101605, 16224, 16176

# Scored feature-matrix sizes
python -c "import pyarrow.parquet as pq, glob; print({s: sum(pq.read_metadata(f).num_rows for f in glob.glob(f'final_dataset_parts_{s}/*.parquet')) for s in ['train','cal','test']})"
# -> {'train': 86378, 'cal': 13201, 'test': 12385}

# Internal predictions
wc -l final_model_output_regularized/test_predictions.csv
# -> 12386 (12,385 variants + header)
```

---

## Principle followed

Only facts that can be verified directly from repository outputs are reported. The cause of the 3,790-variant gap between the initial test split and the scored feature matrix is currently unknown and is explicitly described as such. No speculative explanation has been added to the manuscript.
