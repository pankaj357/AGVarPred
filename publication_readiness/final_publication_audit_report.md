# Final Publication Audit Report — AGVarPred Manuscript

**Audit date:** 2026-06-29  
**Auditor role:** Senior editor, Bioinformatics / Briefings in Bioinformatics / Genome Biology  
**Scope:** Abstract, Introduction, Methods, Results, Discussion/Limitations, figure/table plans, supplementary checklist, and all supporting repository outputs.  
**Constraint:** No new analyses, model training, rescoring, or statistical tests were performed. Only existing files were inspected.

---

## Executive Summary

The revised manuscript is substantially more conservative and internally consistent than the previous draft. However, a critical numerical inconsistency regarding the internal test-set evaluation size (12,385 scored variants vs. 16,175 stated variants) is currently undisclosed, and the manuscript text contains no figure or table references. These issues, together with missing Figure 1 and incomplete supplementary materials, prevent submission in its current state.

**Verdict: Major Revision required before submission to any of the three journals.**

| Journal | Readiness (as-is) | Likely path |
|---|---|---|
| **Bioinformatics** | ~70% | Major revision → acceptable after fixing critical issues |
| **Briefings in Bioinformatics** | ~65% | Major revision → acceptable if figure/table issues resolved |
| **Genome Biology** | ~45% | Major revision; conceptual concerns (AlphaGenome access, limited novelty) remain even after fixes |

---

## 1. Numerical Consistency Audit

### 1.1 Internal ClinVar test performance

| Metric | Manuscript value | Source value | Status |
|---|---|---|---|
| Test-set variants stated | 16,175 | `test.csv` has 16,175 | ✓ |
| **Scored variants for metrics** | **Not disclosed** | **12,385** (`test_predictions.csv` header+12,385 rows; confusion matrix sums to 12,385) | **✗ Critical** |
| ROC-AUC | 0.9753 | 0.975308 (`final_metrics.csv`) | ✓ |
| 95% CI | 0.9730–0.9774 | 0.973016–0.977423 (`final_metrics.csv`) | ✓ |
| PR-AUC | 0.9606 | 0.960624 | ✓ |
| F1 | 0.8945 | 0.894478 | ✓ |
| Precision | 0.8952 | 0.895163 | ✓ |
| Recall | 0.8938 | 0.893794 | ✓ |
| MCC | 0.8327 | 0.832719 | ✓ |
| Brier | 0.0568 | 0.056784 | ✓ |
| ECE | 0.0225 | 0.022465 | ✓ |
| Dummy AUC | 0.4959 | 0.495894 | ✓ |
| Specificity at 0.42 | "approximately 0.95" | 0.9387 (7,330 TN / 7,809 TN+FP from confusion matrix) | ✗ Minor |

**Finding:** The manuscript repeatedly states that the internal test set contains 16,175 variants and reports metrics as if they were computed on that full set. In fact, `final_model_output_regularized/test_predictions.csv` and the confusion matrix contain only 12,385 scored variants. The 3,790 unscored variants (23.4% of the stated test set) are not explained. This is the single most serious issue in the manuscript.

### 1.2 External benchmark performance

| Benchmark | Metric | Manuscript | Source | Status |
|---|---|---|---|---|
| Humsavar | N | 14,718 (1,316 P / 13,402 B) | `humsavar_full_coverage.csv`: N=14,718 | ✓ |
| Humsavar | ROC-AUC | 0.826 (0.817–0.835) | 0.8262 (0.8169–0.8351) | ✓ |
| Humsavar | PR-AUC | 0.275 | 0.2747 | ✓ |
| Humsavar | MCC | 0.295 (0.277–0.311) | 0.2949 (0.2769–0.3108) | ✓ |
| VIP | N | 3,880 (1,075 P / 2,805 B) | `vip_full_coverage.csv`: N=3,880 | ✓ |
| VIP | ROC-AUC | 0.892 (0.881–0.903) | 0.8921 (0.8807–0.9034) | ✓ |
| VIP | PR-AUC | 0.753 | 0.7529 | ✓ |
| VIP | MCC | 0.507 (0.484–0.529) | 0.5069 (0.4845–0.5290) | ✓ |
| Grimm2015 | N | 5,258 (509 P / 4,749 B) | `grimm2015_full_coverage.csv`: N=5,258 | ✓ |
| Grimm2015 | ROC-AUC | 0.745 (0.726–0.763) | 0.7449 (0.7255–0.7631) | ✓ |
| Grimm2015 | PR-AUC | 0.188 | 0.1875 | ✓ |
| Grimm2015 | MCC | 0.233 (0.209–0.258) | 0.2328 (0.2089–0.2583) | ✓ |
| DVD | N | 8,388 (195 P / 8,193 B) | `dvd_full_coverage.csv`: N=8,388 | ✓ |
| DVD | ROC-AUC | 0.977 (0.967–0.985) | 0.9767 (0.9674–0.9854) | ✓ |
| DVD | PR-AUC | 0.788 | 0.7878 | ✓ |
| DVD | MCC | 0.258 (0.239–0.275) | 0.2581 (0.2395–0.2747) | ✓ |
| MAVE Independent | N | 7,557 (3,757 P / 3,800 B) | `mave_independent_full_coverage.csv`: N=7,557 | ✓ |
| MAVE Independent | ROC-AUC | 0.528 (0.514–0.541) | 0.5285 (0.5143–0.5410) | ✓ |
| MAVE Independent | PR-AUC | 0.511 | 0.5113 | ✓ |
| MAVE Independent | MCC | 0.051 (0.026–0.073) | 0.0511 (0.0264–0.0728) | ✓ |
| gnomAD Common Benign | N | 4,989 | `gnomad_benign_full_coverage.csv`: N=4,989 | ✓ |
| gnomAD Common Benign | FPR | 0.0024 (0.001–0.004) | 0.0024 (0.0012–0.0038) | ✓ |
| gnomAD Common Benign | Accuracy | 0.998 | 0.9976 | ⚠ Close but not identical; manuscript rounds up |

**Note on gnomAD accuracy:** The source value is 0.9976; the manuscript reports 0.998. This is a minor rounding inconsistency. Either report 0.9976 or explicitly round to three decimals.

### 1.3 Competitor values on external benchmarks

| Benchmark | Competitor | Manuscript AUC | Source AUC | Common subset N | Status |
|---|---|---|---|---|---|
| Humsavar | CADD | 0.871 (pairwise) | 0.8712 | 14,716 | ✓ |
| Humsavar | AlphaMissense | 0.882 (pairwise) | 0.8823 | 4,707 | ✓ |
| Humsavar | REVEL | 0.942 (pairwise) | 0.9423 | 14,379 | ✓ |
| Humsavar | CADD (all-tools) | 0.900 | Cannot verify from available files | 4,631 (per manuscript) | ⚠ Source not inspected |
| Humsavar | AlphaMissense (all-tools) | 0.883 | Cannot verify | 4,631 | ⚠ Source not inspected |
| Humsavar | REVEL (all-tools) | 0.943 | Cannot verify | 4,631 | ⚠ Source not inspected |
| VIP | CADD | 0.927 | 0.9270 | 3,862 | ✓ |
| VIP | AlphaMissense | 0.803 | 0.8031 | 153 | ✓ |
| VIP | REVEL | 0.876 | 0.8763 | 1,715 | ✓ |
| DVD | CADD | 0.997 | 0.9974 | 8,279 | ✓ |
| DVD | REVEL | 0.979 | 0.9787 | 855 | ✓ |
| Grimm2015 | CADD | 0.695 (pairwise) | 0.6951 | 3,819 | ✓ |
| Grimm2015 | REVEL | 0.641 (pairwise) | 0.6408 | 5,188 | ✓ |
| MAVE Independent | CADD | 0.562 | 0.5620 | 6,663 | ✓ |
| MAVE Independent | AlphaMissense | 0.697 | 0.6968 | 483 | ✓ |
| MAVE Independent | REVEL | 0.699 | 0.6991 | 2,299 | ✓ |

**Actionable item:** The all-tools common-subset AUCs for Humsavar (4,631 variants) cannot be verified from the DeLong CSV or full-coverage CSVs alone. The source file for the all-tools common subset must be identified and archived.

### 1.4 Non-missense comparison

| Benchmark | Method | Manuscript AUC | Source AUC | N | Status |
|---|---|---|---|---|---|
| Humsavar | AGVarPred | 0.683 (0.668–0.698) | 0.6833 (0.6677–0.6982) | 5,901 | ✓ |
| Humsavar | CADD | 0.849 | 0.8488 | 5,899 | ✓ |
| Humsavar | AlphaMissense | 0.882 | 0.8815 | 1,867 | ✓ |
| Humsavar | REVEL | 0.935 | 0.9348 | 5,583 | ✓ |
| VIP | AGVarPred | 0.921 (0.908–0.934) | 0.9213 (0.9080–0.9339) | 2,608 | ✓ |
| VIP | CADD | 0.959 | 0.9594 | 2,596 | ✓ |
| VIP | REVEL | 0.878 | 0.8783 | 451 | ✓ |
| VIP | AlphaMissense | 0.811 | 0.8106 | 62 | ✓ |
| DVD | AGVarPred | 0.981 (0.970–0.990) | 0.9812 (0.9698–0.9900) | 7,663 | ✓ |
| DVD | CADD | 0.998 | 0.9983 | 7,555 | ✓ |
| DVD | REVEL | 0.992 | 0.9921 | 131 | ✓ |
| DVD | AlphaMissense | 0.979 | 0.9789 | 73 | ✓ |
| Grimm2015 | AGVarPred | 0.752 (0.708–0.791) | 0.7522 (0.7076–0.7909) | 873 | ✓ |
| Grimm2015 | CADD | 0.731 | 0.7314 | 593 | ✓ |
| Grimm2015 | REVEL | 0.700 | 0.6998 | 823 | ✓ |

### 1.5 DeLong tests

All reported DeLong AUC differences, q-values, and sample sizes match `delong_analysis/results/delong_results_fdr.csv`.

### 1.6 Decision curve analysis

| Benchmark | Metric | Manuscript | Source (`decision_curve_summary.csv`) | Status |
|---|---|---|---|---|
| Internal test | Max NB | 0.366 at 0.01 | 0.3662 at 0.01 | ✓ |
| Internal test | Integrated NB | 0.280 | 0.2802 | ✓ |
| Humsavar | Max NB | 0.082 at 0.01 | 0.0817 at 0.01 | ✓ |
| Humsavar | Integrated NB | −0.239 | −0.2386 | ✓ |
| VIP | Max NB | 0.266 at 0.01 | 0.2660 at 0.01 | ✓ |
| VIP | Integrated NB | −0.072 | −0.0719 | ✓ |
| Grimm2015 | Max NB | 0.086 at 0.01 | 0.0862 at 0.01 | ✓ |
| Grimm2015 | Integrated NB | −0.436 | −0.4363 | ✓ |
| DVD | Max NB | 0.015 at 0.01 | 0.0152 at 0.01 | ✓ |
| DVD | Integrated NB | −0.146 | −0.1456 | ✓ |
| MAVE Independent | Exceeds both refs | 0.49–0.53 | 0.49–0.53 | ✓ |

### 1.7 SHAP values

| Claim | Manuscript | Source (`feature_importance_shap.csv`) | Status |
|---|---|---|---|
| gnomAD AF mean |SHAP| | 4.648 | 4.6476 | ✓ |
| VEP IMPACT mean |SHAP| | 0.437 | 0.4374 | ✓ |
| Top AlphaGenome feature (CAGE eye) | 0.290 | 0.2896 | ✓ |
| Second (splice ascending aorta) | 0.278 | 0.2780 | ✓ |
| Third (splice Purkinje) | 0.225 | 0.2249 | ✓ |

### 1.8 Tissue SHAP aggregation

| Tissue group | Manuscript % | Source % (`tissue_shap_importance.csv`) | Status |
|---|---|---|---|
| Brain / CNS | 22.0% | 21.9944% | ✓ |
| Liver / digestive | 19.7% | 19.6764% | ✓ |
| Heart / vasculature | 9.8% | 9.7888% | ✓ |
| Blood / immune | 8.4% | 8.3610% | ✓ |
| Reproductive | 8.1% | 8.1019% | ✓ |
| Embryonic / stem | 7.9% | 7.9007% | ✓ |

### 1.9 Model and feature parameters

| Parameter | Manuscript | Source | Status |
|---|---|---|---|
| Total selected features | 120 | `feature_selection_output_nested/selected_features.txt` has 120 lines | ✓ |
| AlphaGenome features | 111 | Counted in selected features list | ✓ |
| gnomAD AF features | 2 | `gnomAD_AF`, `log10_gnomAD_AF` | ✓ |
| VEP features | 7 | Counted | ✓ |
| No-AF feature count | 119 | `model_6_minus_af_output/ablation_summary.csv` | ✓ |
| Full-model threshold | 0.42 | `model_full/manifest.yaml` | ✓ |
| No-AF threshold | 0.45 | `model_no_af/manifest.yaml` | ✓ |
| No-AF internal AUC | 0.9640 | 0.963972 | ✓ |
| No-AF gnomAD FPR | 0.037 | `gnomad_benign_full_coverage.csv`: 0.0375 | ✓ |

### 1.10 Dataset sizes

| Split | Manuscript | Source | Status |
|---|---|---|---|
| Total ClinVar gold | 134,002 | `clinvar_gold_grch38_clean.csv` | ✓ |
| Genes | 9,222 | Source | ✓ |
| Pathogenic | 40,584 (30.3%) | Source | ✓ |
| Benign | 93,418 (69.7%) | Source | ✓ |
| Training set | 101,604 / 7,377 genes | `train.csv` | ✓ |
| Calibration set | 16,223 / 922 genes | `cal.csv` | ✓ |
| Test set | 16,175 / 923 genes | `test.csv` | ✓ |

---

## 2. Claim Support Audit

### 2.1 Fully supported claims
- Internal ClinVar test AUC, PR-AUC, F1, MCC, Brier, and ECE.
- External benchmark AUCs and sample sizes.
- 100% SNV coverage on external benchmarks.
- gnomAD AF dominance by SHAP.
- Descriptive nature of tissue/assay SHAP rankings.
- Near-random MAVE performance and its interpretation as function-vs-pathogenicity distinction.
- No-AF model retains high AUC but higher FPR.

### 2.2 Unsupported or misleading claims

**Critical:**
1. **"The held-out ClinVar test set … contained 16,175 variants" combined with reported metrics.** The metrics are computed on 12,385 scored variants. The manuscript implies evaluation on the full 16,175.
2. **"AGVarPred scored 100% of evaluated SNVs" in external benchmarks** is true, but the internal test set had 24% unscored variants, creating an unstated contrast.

**High:**
3. **"At the 0.42 threshold, the model achieved … a specificity of approximately 0.95."** Actual specificity on the scored subset is 0.9387 (~0.94).
4. **"AGVarPred achieved an accuracy of 0.998" on gnomAD Common Benign.** Source is 0.9976; rounding to 0.998 without explanation overstates.

**Medium:**
5. **"On the common subset scored by all four tools (4,631 variants), AGVarPred AUC was 0.833 …"** This all-tools common subset is not present in the inspected files; its source must be verified.
6. **"Calibration quality was also good on the internal test set"** is supported, but the claim is restricted to the scored subset and to ClinVar-like data; this limitation is not explicit in Results.

### 2.3 Overstated claims (now mostly fixed)
The revised manuscript has removed most earlier overstatements. Remaining soft overstated phrases:
- "strong discrimination" and "good calibration" in the Discussion are acceptable but could be quantified inline.
- "competitive performance on clinically curated benchmarks" is fair for DVD/VIP but should acknowledge CADD/REVEL superiority on common subsets.

---

## 3. Logical Consistency Audit

### 3.1 Abstract vs. Results
- Abstract numbers match Results.
- Abstract states "AGVarPred scored 100% of evaluated SNVs" — true for external benchmarks but the internal test-set scoring gap is not mentioned. This is acceptable in an abstract but the main text must disclose the gap.
- Abstract framing is appropriately conservative.

### 3.2 Discussion vs. Results
- Discussion correctly summarizes the benchmark-dependent performance pattern.
- Discussion states the no-AF model FPR of 0.037 versus 0.0024 — supported.
- Discussion notes that external calibration was not evaluated — consistent with Methods.
- Discussion states AGVarPred does not outperform best missense predictors on missense-rich benchmarks — consistent with Results.

### 3.3 Methods vs. Results
- Methods describe the same splits, features, model, calibration, and evaluation as Results reports.
- Methods state calibration was assessed only on internal ClinVar test; Results do not contradict this.
- Methods note the 0.42 threshold was selected on 40.3% prevalence; Results apply the same threshold externally, which is methodologically consistent but should be flagged as a limitation more prominently in Results.

### 3.4 Figure legends vs. manuscript
- **The manuscript text contains zero figure or table citations.** This is a critical omission. The figure audit (Fig 2–7) and table audit (Table 1–8) are not referenced anywhere in the manuscript sections.
- The previous Results draft contained citations to Figure 5A–E and Figure 6A–D; those numbering assumptions conflict with the current figure audit (Fig 4 = variant-type, Fig 5 = non-missense, Fig 6 = SHAP, Fig 7 = DCA).

### 3.5 Internal inconsistency regarding test-set size
- Methods: "Internal performance was evaluated on the held-out ClinVar test set of 16,175 variants."
- Results: "On this held-out test set, AGVarPred achieved …"
- Confusion matrix / predictions: 12,385 scored variants.
- This is the most damaging internal inconsistency.

---

## 4. Figure and Table Reference Audit

### 4.1 Figures
| Figure | Audit status | Cited in manuscript? | Issue |
|---|---|---|---|
| Fig 1 | MISSING | No | Must be created |
| Fig 2 | READY | No | Add citation in Results §1 |
| Fig 3 | READY | No | Add citation in Results §5 |
| Fig 4 | READY | No | Add citation in Results §3/§9 |
| Fig 5 | READY | No | Add citation in Results §4 |
| Fig 6 | READY | No | Add citation in Results §7/§8 |
| Fig 7 | READY | No | Add citation in Results §6 |

**Critical finding:** No figure is cited in the manuscript text. Every figure needs a citation at first mention and, where appropriate, at the end of the relevant subsection.

### 4.2 Tables
| Table | Audit status | Cited in manuscript? | Issue |
|---|---|---|---|
| Table 1 | READY | No | Add citation in Methods §2 or Results §1 |
| Table 2 | READY | No | Add citation in Results §2 |
| Table 3 | READY | No | Add citation in Results §5 |
| Table 4 | READY | No | Add citation in Results §9 |
| Table 5 | READY | No | Add citation in Results §3 |
| Table 6 | READY | No | Add citation in Results §4 |
| Table 7 | READY | No | Add citation in Results §7 |
| Table 8 | READY | No | Add citation in Results §6 |

**Critical finding:** No table is cited in the manuscript text.

### 4.3 Supplementary items
- Supplementary Methods, Note 1, and Note 2 are still placeholders and must be prepared.
- Supplementary figures S1–S5 and tables S1–S4 are listed as ready but are not cited in the main text.

---

## 5. Terminology Consistency

| Term | Status | Notes |
|---|---|---|
| "genome-wide SNV coverage" | ✓ Consistent | Used throughout; no reversion to unqualified "genome-wide coverage" |
| "gene-level split" | ✓ Consistent | Methods and Results use identical phrasing |
| "AlphaGenome" | ✓ Consistent | Capitalization and usage consistent |
| "gnomAD AF" / "gnomAD allele frequency" | ⚠ Mostly consistent | Prose uses "gnomAD allele frequency"; feature names use "gnomAD AF". Acceptable but a short glossary might help. |
| "external benchmark" | ✓ Consistent | |
| "common subset" | ✓ Consistent | |
| "non-missense" | ✓ Consistent | |
| "ClinVar test set" vs. "scored ClinVar test subset" | ✗ Inconsistent | The manuscript must distinguish the 16,175-variant file from the 12,385-variant scored subset. |

---

## 6. Duplicated or Redundant Text

1. **Benchmark descriptions** appear in both Methods §9 and Results §2. This is acceptable but lengthy; consider moving full benchmark provenance to Methods and keeping only key numbers and interpretation in Results.
2. **MAVE interpretation** is stated in Results §2 and reinforced in Results §9 and Discussion. This repetition is justified because it is a subtle conceptual point.
3. **AF-dominance caveat** appears in Results §7, Results §8, Methods §11, and Discussion. Some consolidation may improve flow, but the repetition is not excessive for a high-impact journal where reviewers may read sections independently.
4. **No figure/table references** means there is no duplication issue for citations, but the absence itself is a problem.

---

## 7. Writing Quality

### 7.1 Strengths
- Conservative, cautious tone throughout.
- Clear section headings.
- Technical terms are defined when first introduced.
- Complex statistical concepts (DeLong bias, coverage bias, prevalence dependence) are explained well.

### 7.2 Issues
- **Sentence length:** Some sentences in Methods §5 (feature selection) and Results §4 (non-missense comparison) are long and could be split.
- **Tense consistency:** Mostly present tense for results, past tense for what was done. Acceptable.
- **Passive voice:** Methods uses passive voice appropriately; Results could occasionally use more active constructions.
- **Ambiguous phrasing:** "These results provide an estimate of AGVarPred's ability … within the limits of ClinVar's label distribution and ascertainment bias" is good, but the earlier "16,175 variants" claim undermines it.
- **Missing transitions:** The manuscript sections do not currently flow into one another because they are separate files; ensure transitions when assembled.

### 7.3 Recommended edits
- Break up the feature-selection paragraph in Methods §5.
- Add a one-sentence transition between Introduction and Methods when assembling.
- Standardize whether "e.g.," or "for example" is used (currently mixed).

---

## 8. Limitations Correspondence

| Limitation in Discussion | Corresponds to actual study limitation? | Evidence |
|---|---|---|
| Dependence on AlphaGenome | ✓ Yes | SDK v0.6.1, API keys required, not publicly available at time of writing |
| Restriction to SNVs | ✓ Yes | Only SNVs retained in dataset construction; no indel/SV evaluation |
| Training/calibration on ClinVar | ✓ Yes | Labels from ClinVar; calibration only on internal test |
| External calibration not evaluated | ✓ Yes | No external calibration curves or ECE reported |
| Benchmark independence only vs. ClinVar | ✓ Yes | DeLong/common-subset comparisons do not control for competitor training data |
| SHAP interpretation caveats | ✓ Yes | Correlated features, heuristic tissue mapping |
| Clinical utility limitations | ✓ Yes | DCA integrated NB negative on most external benchmarks |
| Competitor set limited | ✓ Yes | Only CADD, AlphaMissense, REVEL compared |
| Small subgroup power | ✓ Yes | Many per-class subgroups have <10 pathogenic variants |
| Git history / metadata gaps | ✓ Yes | ClinVar date and VEP version not recorded |

All stated limitations are real and proportionally discussed.

---

## 9. Prioritized Issues

### Critical (must fix before any submission)

1. **Disclose the internal test-set scoring gap.** The manuscript reports metrics on 12,385 scored variants but describes the test set as 16,175. Add a sentence in Methods §8/Results §1 explaining that 12,385 of 16,175 test variants had complete feature data and that all reported internal metrics are on this scored subset. Investigate and report why 3,790 variants were excluded.
2. **Add figure citations throughout.** Every figure (Figs 1–7) must be cited in the text.
3. **Add table citations throughout.** Every table (Tables 1–8) must be cited in the text.
4. **Create Figure 1.** A workflow/overview figure is required and currently missing.
5. **Prepare Supplementary Methods and Notes 1–2.** These are placeholders.

### High (strongly recommended)

6. **Correct specificity wording** from "approximately 0.95" to "approximately 0.94" or report the exact value 0.9387.
7. **Correct gnomAD Common Benign accuracy** to 0.9976 or explicitly round to three decimals.
8. **Verify the Humsavar all-tools common-subset** (4,631 variants) and ensure its source file is archived.
9. **Resolve figure-numbering alignment.** Decide whether coverage = Fig 4 / non-missense = Fig 5 / SHAP = Fig 6 / DCA = Fig 7, and update all citations accordingly. Do not use Figure 5A–E / 6A–D unless the figure files are actually organized that way.
10. **Add per-benchmark prevalence discussion in Results** to explain why the fixed 0.42 threshold is applied externally.

### Medium (should fix)

11. **Recover or explicitly note missing metadata:** ClinVar release date and Ensembl VEP version.
12. **Add a supplementary table or note** showing the train/cal/test gene and variant counts.
13. **Shorten the Methods §5 feature-selection paragraph** and split long sentences.
14. **Add a transition paragraph/sentence** between sections when assembling the full manuscript.
15. **Cite supplementary figures/tables** in the main text where appropriate.

### Minor (polish)

16. Standardize "e.g." vs. "for example."
17. Ensure all journal-specific formatting requirements (word count, reference style) are met at submission.
18. Verify Zenodo DOI resolves.
19. Purge Git history of any remaining secrets.
20. Run a final spell-check and grammar check.

---

## 10. Journal-Specific Assessment

### Bioinformatics
- **Fit:** Good. The manuscript is a methods/benchmarking paper with solid engineering.
- **Concerns:** The undisclosed test-set scoring gap would be flagged by reviewers. Missing Figure 1 and supplementary notes are desk-rejection risks if not fixed.
- **Readiness after critical fixes:** ~85%.

### Briefings in Bioinformatics
- **Fit:** Good, though the biological interpretation is deliberately conservative, which may be seen as a strength or weakness depending on the editor.
- **Concerns:** Same as Bioinformatics, plus the need for clearer figure/table integration.
- **Readiness after critical fixes:** ~80%.

### Genome Biology
- **Fit:** Moderate. The work is more methods-focused than typical Genome Biology papers.
- **Concerns:**
  - AlphaGenome dependency remains a fundamental reproducibility barrier.
  - Novelty is incremental (a LightGBM model on a new but inaccessible feature set).
  - External calibration is not evaluated.
  - Competitor comparison is limited.
- **Readiness after critical fixes:** ~60%. Even with all fixes, the paper may be seen as better suited to a methods journal unless a stronger biological or clinical narrative is developed.

---

## 11. Final Recommendation

### **Major Revision**

The manuscript cannot be submitted as-is. The most urgent issues are:

1. The internal test-set evaluation size is misrepresented (12,385 scored vs. 16,175 stated).
2. No figures or tables are cited in the text.
3. Figure 1 is missing.
4. Supplementary Methods and Notes are placeholders.

After these issues are resolved, the manuscript will be suitable for submission to **Bioinformatics** or **Briefings in Bioinformatics** with a realistic chance of positive review. Submission to **Genome Biology** would require additional conceptual strengthening (broader competitor comparison, external calibration analysis, and ideally a pathway to public AlphaGenome access) beyond the scope of a purely editorial revision.

The revised framing is appropriate and the numerical results are largely accurate, but the manuscript currently lacks the basic structural elements required for peer review.
