# AGVarPred Manuscript Assembly

This directory contains the revised manuscript sections for the AGVarPred Genome Biology submission. The sections are written as standalone Markdown files and can be concatenated in the order below to produce the full manuscript text.

## Section order

1. `abstract.md` — Abstract
2. `introduction.md` — Introduction
3. `methods_section.md` — Methods
4. `results_section.md` — Results
5. `discussion_section.md` — Discussion and Limitations

## Key revisions from the previous draft

- **Coverage claims** have been reframed as **coverage of evaluated variants** rather than genome-wide coverage.
- Claims about biological mechanisms, clinical utility, and superiority over competitors have been substantially toned down.
- The manuscript now accurately reflects that the ClinVar gold-standard dataset contained **multiple variant types**; the false claim that "only single-nucleotide variants were retained" has been removed.
- Internal evaluation is now explicitly described as performed on the **12,385 scored test-set variants** represented in the final feature matrices, not the full 16,175-variant initial split.
- The dominant role of **gnomAD allele frequency** is now emphasized; AlphaGenome features are described as providing complementary functional signal.
- The **MAVE Independent** benchmark is interpreted as a distinction between molecular function and clinical pathogenicity, not as a flawed benchmark.
- **Small subgroup analyses** are explicitly flagged as exploratory or uninformative.
- A balanced **Limitations** subsection has been added covering AlphaGenome dependence, restriction to variants represented in feature matrices, SNV restriction, SHAP interpretation, benchmark composition, external calibration, and clinical utility.
- Only verified numerical results are reported; speculative claims about missing variants have been replaced with explicit statements that the cause is unknown.

## Figure plan

A proposed set of main-text and supplementary figures, with source files and citation locations, is documented in `figure_plan.md`. The figure numbering follows `figure_audit.csv` and has been reconciled with the manuscript section order. Figure references for Figures 1–7 have been inserted into `methods_section.md` and `results_section.md`.

## Correction package

The following files document the corrections applied after the test-set size discrepancy was discovered:

- `change_log.md` — Line-by-line record of every corrected claim, with original wording, corrected wording, and repository evidence.
- `final_consistency_report.md` — Verification that all numerical claims in the corrected sections match repository outputs and that no unsupported claims remain.

Earlier audit and readiness reports (`final_publication_audit_report.md`, `publication_readiness_report.md`) reflect the state before this correction package and should be read in conjunction with `change_log.md`.

## Remaining manuscript-assembly tasks

- Add title page, author list, affiliations, and correspondence.
- Generate Figure 1 (workflow overview), which is currently missing.
- Add table references and any additional tables (e.g., benchmark summary, feature list).
- Prepare supplementary materials (Supplementary Methods, Notes 1–2, tables, and figures).
- Add references.
- Verify exact ClinVar release date and Ensembl VEP version; add or explicitly note as missing.
- Investigate why 3,790 of the 16,175 test-split variants are absent from the final scored feature matrices.
- Final copy-edit for journal style and word count.
