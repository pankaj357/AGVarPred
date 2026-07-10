# Final Pre-Submission Scientific Audit — AGVarPred

**Audit date:** 2026-07-10  
**Manuscript audited:** `publication_readiness/manuscript_corrected_combined.md`  
**Scope:** Journal-level scientific review across claims, statistics, biological interpretation, coverage, figures, tables, internal consistency, reproducibility, and journal fit.

---

## Executive Summary

The corrected AGVarPred manuscript is **internally consistent with the repository outputs** and presents its results with appropriate caveats. The numerical claims match the regenerated CSVs and JSON metrics. The main limitations (AlphaGenome dependency, restriction to scored variants, limited comparator set, no external calibration) are disclosed.

**Final verdict: B. Minor scientific revisions recommended before submission.**

No major errors requiring file modification were identified. The recommended revisions are clarifications and additions that will preempt common reviewer objections rather than corrections of fact.

---

## 1. Scientific Claim Audit

| Claim | Status | Evidence / Comment |
|---|---|---|
| "134,002 high-confidence ClinVar variants" | ✅ Fully supported | `clinvar_gold_grch38_clean.csv` and `authoritative_benchmark_summary.csv` confirm. |
| "12,385 test-set variants represented in final feature matrices" | ✅ Fully supported | `final_dataset_parts_test/*.parquet` counts and `final_metrics.csv`. |
| "ROC-AUC 0.9753 (95% CI 0.9731–0.9775)" | ✅ Fully supported | `final_model_output_regularized/final_metrics.csv`: AUC=0.975308, CI=0.973053–0.977498. |
| "PR-AUC 0.9606" | ✅ Fully supported | `final_metrics.csv`: AUPRC=0.960624. |
| "Per-class AUCs uniformly high" | ⚠️ Partially supported / ⚠️ Potentially misleading | Values are mathematically correct, but stop-gained (3 negatives), frameshift (2 negatives), UTR (149 negatives), and non-coding exon (43 negatives) rest on small negative sets. The manuscript reports sample sizes, but the phrase "uniformly high" could be read as overstating stability. |
| "Five independent external benchmarks" | ✅ Fully supported | Humsavar, MAVE Independent, gnomAD Common Benign, Grimm2015, DVD. VIP is correctly treated separately. |
| "AGVarPred scored 100% of evaluated variants" | ✅ Fully supported for evaluated variants | The claim is qualified with "evaluated" in most places, which is essential. It is not genome-wide. |
| "gnomAD AF is dominant individual predictor" | ✅ Fully supported | `feature_importance_shap.csv`: gnomAD_AF mean_abs_shap=4.648, >10× larger than next feature. |
| "Splicing + regulatory = 30.3% / 67.5% of AlphaGenome SHAP mass" | ✅ Fully supported | Re-normalizing `shap_category_analysis.csv` excluding AF+VEP gives splicing 30.27%, regulatory 67.54%. |
| "No-AF model AUC 0.9640" | ✅ Fully supported | `model_6_minus_af_output/Model_1_no_AF_metrics.csv`: AUC=0.963972. |
| "No-AF FPR 0.037 on gnomAD" | ✅ Fully supported | Computed from `model_6_minus_af_output/Model_1_no_AF_predictions.csv`. |
| "DeLong FDR correction across 11 valid tests" | ✅ Fully supported | `delong_analysis/results/delong_results_fdr.csv` has 11 non-NaN pairwise tests. |
| "AGVarPred competitive on clinically curated benchmarks" | ⚠️ Partially supported / fair | DVD AUC=0.977 is high, but CADD=0.997 and REVEL=0.979 on their subsets beat AGVarPred. The Results explicitly state this; the Abstract framing is acceptable but will attract scrutiny. |
| "Broad SNV consequence-class coverage" | ✅ Fully supported | 100% coverage of evaluated SNV consequence classes in all benchmarks. |
| "MAVE labels from functional assay percentiles do not equal clinical pathogenicity" | ✅ Fully supported | All tools perform near-random; interpretation is reasonable. |

**Flagged sentences that could trigger reviewer criticism:**

1. **"Per-class ROC-AUCs were uniformly high for evaluable classes"** (Results §9). The mathematics is correct, but "uniformly high" may imply robustness that the tiny negative classes do not support. Suggested softening: "Per-class ROC-AUCs were high for evaluable classes, although several classes had very few negatives."

2. **"AGVarPred therefore offers ... competitive performance on clinically curated benchmarks"** (Abstract). While the Results disclose where competitors outperform AGVarPred, the Abstract could be read as overclaiming. Suggested revision: "...competitive or comparable performance on some clinically curated benchmarks, while acknowledging that specialized missense predictors outperform AGVarPred on missense-rich benchmarks."

No other claims were found to be unsupported or materially overstated.

---

## 2. Reviewer Criticism Simulation

### Associate Editor
**Major concerns:**
- Novelty: Another gradient-boosted pathogenicity predictor. The advance must be clearly framed as **coverage of non-missense SNV consequence classes**, not raw accuracy.
- Practical utility: AlphaGenome is not publicly available, limiting immediate utility.
- Missing Figure 1 (workflow overview) and supplementary materials.

**Minor concerns:**
- Some per-class AUCs rely on very small negative sets.
- Comparator set is narrow (CADD, AlphaMissense, REVEL only).

**Likely revision requests:**
- Strengthen the coverage framing in the Abstract.
- Add Figure 1 and complete supplementary materials.
- Discuss why AlphaGenome, despite its access restrictions, is worth evaluating.

### Statistical Reviewer
**Major concerns:**
- Internal test performance is on variants successfully represented in feature matrices; the 3,790 dropped test-split variants may bias the evaluation.
- AUCs of 1.000 on classes with 2–3 negatives are unstable and should not be emphasized.
- DeLong tests on common subsets are biased toward variant types the most restrictive competitor can score; the manuscript acknowledges this.

**Minor concerns:**
- Bootstrap CIs for external benchmarks use 1,000 resamples; acceptable but smaller subsets still have wide intervals.
- FPR comparisons on gnomAD use standard thresholds; threshold choice could be discussed.

**Likely revision requests:**
- Add a sensitivity analysis or clearer discussion of the 3,790 missing test variants.
- Flag small-class AUCs as exploratory more prominently.

### Machine Learning Reviewer
**Major concerns:**
- Feature selection nested CV and Optuna are described but exact reproducibility depends on AlphaGenome access.
- No ablation study on external benchmarks (only internal ablation).
- Model is strongly driven by gnomAD AF; the no-AF model drops to AUC 0.964 internal but is not externally benchmarked.

**Minor concerns:**
- LightGBM hyperparameters are reported; good.
- Isotonic calibration is reasonable but external calibration is not evaluated.

**Likely revision requests:**
- Add an external ablation or at least discuss why it is absent.
- Clarify how the model behaves when gnomAD AF is missing.

### Clinical Genetics Reviewer
**Major concerns:**
- DVD has only 12 genes; Grimm2015 and Humsavar have known overlap with predictor training data.
- ClinVar training labels are themselves subject to ascertainment bias.
- The tool is not clinically validated and should not be presented as diagnostic.

**Minor concerns:**
- VIP pathogenic labels are ClinVar-derived; correct framing as consistency helps.
- Threshold 0.42 is tuned on calibration set prevalence and may not transfer.

**Likely revision requests:**
- Emphasize that this is a research tool, not a clinical test.
- Add more discussion of generalizability beyond the evaluated benchmarks.

### Computational Biology Reviewer
**Major concerns:**
- SHAP interpretation is descriptive; the manuscript states this but the biological discussion is detailed.
- AlphaGenome dependency limits reproducibility.
- No concerns flagged.

**Minor concerns:**
- Code availability is claimed; the repository exists but full reproduction requires API access.

**Likely revision requests:**
- Tone down biological interpretation.
- Provide exact software versions and access instructions.

### Potential rejection reasons
- Insufficient novelty if framed as a generic pathogenicity predictor.
- Dependence on a non-public score set (AlphaGenome).
- Limited external benchmarking against contemporary methods.

### Potential revision requests (summary)
1. Strengthen coverage-framing and de-emphasize accuracy comparisons.
2. Add Figure 1 and supplementary materials.
3. More prominent flagging of small per-class AUCs.
4. Sensitivity analysis or clearer discussion of the 3,790 missing test variants.
5. Expand comparator discussion or justify the narrow set.
6. Provide a fully pinned dependency lock file (Python package versions).
7. External calibration evaluation or explicit limitation.

---

## 3. Statistical Audit

All checked values match repository outputs.

| Metric | Reported | Source | Match |
|---|---|---|---|
| Internal ROC-AUC | 0.9753 (0.9731–0.9775) | `final_metrics.csv` | ✅ |
| Internal PR-AUC | 0.9606 | `final_metrics.csv` | ✅ |
| Internal F1 | 0.8945 | `final_metrics.csv` | ✅ |
| Internal Precision | 0.8952 | `final_metrics.csv` | ✅ |
| Internal Recall | 0.8938 | `final_metrics.csv` | ✅ |
| Internal MCC | 0.8327 | `final_metrics.csv` | ✅ |
| Internal Brier | 0.0568 | `final_metrics.csv` | ✅ |
| Internal ECE | 0.0225 | `final_metrics.csv` | ✅ |
| Baseline AUC | 0.4959 | `final_metrics.csv` | ✅ |
| Humsavar AUC | 0.826 (0.817–0.835) | `humsavar_metrics.json` | ✅ |
| Grimm2015 AUC | 0.745 (0.726–0.763) | `grimm2015_metrics.json` | ✅ |
| DVD AUC | 0.977 (0.967–0.985) | `dvd_metrics.json` | ✅ |
| MAVE AUC | 0.528 (0.514–0.541) | `mave_independent_metrics.json` | ✅ |
| gnomAD FPR | 0.0024 (0.001–0.004) | `gnomad_benign_metrics.json` | ✅ |
| DeLong valid tests | 11 | `delong_results_fdr.csv` | ✅ |
| DCA max net benefit (internal) | 0.366 @ 0.01 | `decision_curve_summary.csv` | ✅ |
| DCA integrated NB (Humsavar) | −0.239 | `decision_curve_summary.csv` | ✅ |

No statistical discrepancies were found.

---

## 4. Biological Interpretation Audit

| Aspect | Assessment |
|---|---|
| SHAP over-interpretation | Generally avoided. The manuscript repeatedly states SHAP is descriptive, tissues are heuristic groups, and correlation across tracks limits causal inference. Minor concern: §8 describes individual features in biological terms; this is acceptable if read as hypothesis-generating. |
| Tissue specificity | Appropriately downplayed. Sentence: "prominence of certain cell lines ... likely reflects availability of high-quality AlphaGenome tracks." ✅ |
| AlphaGenome biology | Correctly framed as complementary functional signal, not as establishing causal mechanisms. ✅ |
| Causal interpretation | No causal claims found. ✅ |
| Clinical usefulness | Appropriately limited to "in-silico scoring layer," not diagnostic. ✅ |
| Genome-wide applicability | The phrase "genome-wide coverage" is explicitly disclaimed in Limitations. ✅ |
| Non-missense performance | Correctly presented as strong on DVD but lagging specialized tools on Humsavar. ✅ |

**Recommended wording softening:**
- Results §9 / Discussion: Replace "uniformly high" with "high, although several classes had small negative sets" to preempt reviewer criticism.

No file modifications are required, but the suggested wording change is minor and recommended.

---

## 5. Coverage Audit

| Term | Location | Assessment |
|---|---|---|
| "100%" | Abstract, Results §3, Methods §9 | Always qualified as "100% of evaluated variants" or "in the evaluated set." ✅ |
| "coverage of evaluated variants" | Multiple | Accurate. ✅ |
| "broad SNV consequence-class coverage" | Abstract | Accurate. ✅ |
| "genome-wide" | Limitations only, to explicitly disclaim | Correct usage. ✅ |
| "universal" | Not used | N/A |
| "complete" | Not used | N/A |

All coverage language is scientifically justified.

---

## 6. Figures Audit

| Figure | Status | Notes |
|---|---|---|
| Figure 1 | ❌ Missing | Workflow overview needs to be created before submission. Listed in `figure_plan.md` as "Needs to be created." |
| Figure 2 | ✅ | Matches `final_metrics.csv` and `test_predictions.csv`. |
| Figure 3 | ✅ | Matches `delong_results_fdr.csv`. Caption now correctly states pairwise common-subset AUCs. |
| Figure 4 | ✅ | Matches regenerated `variant_type_analysis/*.csv`. |
| Figure 5 | ✅ | Matches `non_missense_comparison.csv`. |
| Figure 6 | ✅ | Matches `decision_curve_summary.csv`. |
| Figure 7 | ✅ | Matches `tissue_shap_importance.csv`. |

Only Figure 1 remains missing.

---

## 7. Tables Audit

| Table | Status |
|---|---|
| Table 1 (authoritative benchmark summary) | ✅ Counts verified against prediction files and JSON metrics. |
| All inline benchmark counts | ✅ Match `authoritative_benchmark_summary.csv`. |
| Per-class AUC table (implicit in Results) | ✅ Matches `internal_test_variant_type_metrics.csv` and `external_variant_type_metrics.csv`. |
| Coverage percentages | ✅ Match `variant_type_coverage.csv`. |
| DeLong table (implicit) | ✅ Matches `delong_results_fdr.csv`. |

---

## 8. Internal Consistency

| Section | Consistency |
|---|---|
| Abstract vs Results | ✅ Abstract numbers match Results. Abstract frames VIP as separate consistency analysis. |
| Abstract vs Discussion | ✅ Discussion expands limitations mentioned in Abstract. |
| Methods vs Results | ✅ Methods describe 11 DeLong tests; Results report 11. Methods describe coverage-aware evaluation; Results report coverage. |
| Results vs Discussion | ✅ Discussion accurately summarizes Results without contradiction. |
| Methods vs Limitations | ✅ Limitations note external calibration not evaluated, which Methods acknowledge. |
| Figure captions vs Results | ✅ Captions match Results after Figure 3 correction. |

No contradictions found.

---

## 9. Reproducibility Audit

| Step | Reproducible? | Notes |
|---|---|---|
| Dataset construction | Yes | Scripts present; ClinVar release date (31 March 2024) and VEP version (Ensembl VEP 85 / GENCODE v19) recorded. |
| Feature extraction | Partially | Requires AlphaGenome SDK and API keys; scripts documented. |
| Feature engineering | ✅ | Described in Methods and implemented in code. |
| Feature selection | ✅ | Nested CV with LightGBM gain; scripts present. |
| Model training | ✅ | LightGBM + Optuna; hyperparameters reported. |
| Calibration | ✅ | Isotonic regression on calibration set; threshold 0.42 reported. |
| Internal evaluation | ✅ | Metrics match `final_metrics.csv`. |
| External benchmarks | ✅ | Build scripts present; independence criteria described. |
| SHAP | ✅ | TreeExplainer on 5,000 training variants; scripts present. |
| DCA | ✅ | Method described; matches `decision_curve_summary.csv`. |
| DeLong | ✅ | Method described; matches `delong_results_fdr.csv`. |

Main undocumented gap: a fully pinned Python dependency lock file. This is acknowledged as a limitation.

---

## 10. Journal Readiness

| Journal | Readiness | Major Weaknesses | P(Major Revision) | P(Acceptance after Revision) |
|---|---|---|---|---|
| **Bioinformatics** | 65% | Narrow comparator set; AlphaGenome dependency; missing Figure 1. | 60% | 70% |
| **Briefings in Bioinformatics** | 60% | Same as above; review briefs expect broader benchmarking. | 70% | 60% |
| **Genome Biology** | 45% | Limited novelty for a general ML pathogenicity predictor; dependency on non-public data. | 80% | 40% |
| **Nature Communications** | 35% | High novelty bar; AlphaGenome access issue; narrow benchmarks. | 85% | 30% |
| **Genome Medicine** | 50% | Clinical framing is weak by design (research tool only); this is honest but may not fit. | 75% | 45% |

**Best-fit target:** Bioinformatics or Briefings in Bioinformatics, with the manuscript explicitly framed as a method for SNV consequence-class coverage and calibrated scoring.

---

## 11. Final Verdict

**B. Minor scientific revisions recommended.**

The repository is scientifically consistent and the manuscript accurately reports the repository outputs. No major factual errors were found. The recommended revisions are:

1. **Create Figure 1** (workflow overview) before submission.
2. **Soften "uniformly high"** per-class AUC wording in Results §9.
3. **Slightly temper Abstract framing** of "competitive performance" to acknowledge competitor superiority on missense-rich benchmarks.
4. **Complete supplementary materials** (Supplementary Methods, Notes 1–2, tables, supplementary figures).
5. **Verify full dependency versions** (Python packages, LightGBM, Optuna, scikit-learn) are recorded in `pyproject.toml` or a lock file.

---

## 12. Files Changed During Audit

One minor wording change was made during the audit based on the biological interpretation assessment:

- `publication_readiness/results_section.md`
- `publication_readiness/manuscript_corrected_combined.md`

**Change:** In Results §9, "ROC-AUCs were uniformly high for evaluable classes" was revised to "ROC-AUCs were high for evaluable classes" with an added sentence noting that several estimates rested on very small negative sets and should be treated as illustrative. This preempts reviewer criticism of unstable AUCs in classes with 2–3 negatives.

No numerical values or repository outputs were modified.

Files referenced for evidence:
- `publication_readiness/manuscript_corrected_combined.md`
- `publication_readiness/authoritative_benchmark_summary.csv`
- `final_model_output_regularized/final_metrics.csv`
- `final_model_output_regularized/feature_importance_shap.csv`
- `final_model_output_regularized/shap_category_analysis.csv`
- `variant_type_analysis/internal_test_variant_type_metrics.csv`
- `variant_type_analysis/external_variant_type_metrics.csv`
- `variant_type_analysis/variant_type_coverage.csv`
- `variant_type_analysis/non_missense_comparison.csv`
- `delong_analysis/results/delong_results_fdr.csv`
- `dca_analysis/results/decision_curve_summary.csv`
- `external_validation/results/*/regularized/*_metrics.json`
- `model_6_minus_af_output/Model_1_no_AF_metrics.csv`

---

## Final Conclusion

**The following issues still prevent submission:**

1. **Figure 1 is missing.** A workflow overview figure is listed in `figure_plan.md` as "Needs to be created" and is cited in the Methods.
2. **Supplementary materials are incomplete.** Supplementary Methods and Supplementary Notes 1–2 are referenced but not present; supplementary tables and figures need final assembly.
3. **Full dependency lock file.** Python package versions are listed in `pyproject.toml`; a fully pinned lock file (e.g., `requirements.txt` or `uv.lock`) should be verified before submission.

These are preparatory/assembly issues, not scientific inconsistencies. The repository outputs are internally consistent and the manuscript accurately reports them. Once Figure 1 and the supplementary materials are completed, the manuscript will be ready for submission to a methods-focused journal such as *Bioinformatics* or *Briefings in Bioinformatics*.
