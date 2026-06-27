# Model Card: AGVarPred v1.0.0

## Model summary

- **Software version:** 1.0.0
- **Model versions:** `model_full` (primary), `model_no_af` (fallback)
- **Task:** Binary classification of germline variants as Pathogenic (1) or Benign (0)
- **Architecture:** Regularized LightGBM classifier with isotonic calibration
- **Input:** 120 features (`model_full`) or 119 features (`model_no_af`)
- **Output:** Calibrated probability, predicted class, and prediction metadata

## Intended use

AGVarPred is designed to support research and clinical variant interpretation
by providing a quantitative pathogenicity score for single-nucleotide variants
and small indels in protein-coding and splice regions. It is intended to
**complement**, not replace, expert review and orthogonal evidence.

Suitable use cases include:

- Prioritizing variants in rare-disease research.
- Benchmarking against existing pathogenicity predictors.
- Generating hypotheses for further experimental or clinical investigation.

## Out-of-scope use

Do not use AGVarPred for:

- Direct clinical decision-making without expert review.
- Somatic variant classification.
- Structural variants, copy-number variants, or large indels.
- Variants in non-coding regions with sparse AlphaGenome functional coverage.

## Training data

- **Source:** ClinVar germline variants (GRCh38)
- **Labels:** Pathogenic vs. Benign, high-review-status only
- **Conflicting records:** Excluded
- **Split:** 80% train / 10% calibration / 10% test, split by gene
- **External allele frequency:** gnomAD exomes r2.1.1 liftover to GRCh38
- **Functional genomics:** AlphaGenome v0.6.1 recommended variant scorers

## Model details

### Primary model (`model_full`)

- **Features:** 120 selected features
- **Threshold:** 0.42 (F1-optimized on the calibration set)
- **Calibration:** Isotonic regression

### Fallback model (`model_no_af`)

- **Features:** 119 selected features (`gnomAD_AF` excluded)
- **Threshold:** 0.45
- **Calibration:** Isotonic regression
- **Use:** Automatically selected when no local gnomAD VCF or online AF/VEP
  provider is available.

## Feature categories (`model_full`)

| Category              | Count | Examples                                      |
|-----------------------|-------|-----------------------------------------------|
| Histone marks         | 33    | H3K79me2, H3K27ac, H3K9me3                   |
| Splicing              | 30    | Splice-site usage, splice junctions          |
| TF binding            | 18    | ChIP-seq TF tracks                           |
| CAGE                  | 10    | Tissue-specific CAGE signals                 |
| Contact maps          | 9     | 3D genome contact maps                       |
| VEP-derived           | 8     | IMPACT score, PolyPhen, SIFT, LoF, consequence flags |
| DNase                 | 7     | DNase-seq accessibility                      |
| PRO-cap               | 3     | Transcription start-site mapping             |
| Allele frequency      | 1     | gnomAD_AF                                    |
| ATAC                  | 1     | ATAC-seq accessibility                       |

## Performance summary (`model_full` on held-out test set)

| Metric | Value |
|--------|-------|
| AUC    | 0.9753 |
| AUPRC  | 0.9606 |
| F1     | 0.8945 |
| Precision | 0.8952 |
| Recall | 0.8938 |
| MCC    | 0.8327 |
| Brier  | 0.0568 |
| ECE    | 0.0225 |

AUC 95% confidence interval: [0.9730, 0.9774].

The `model_no_af` model is expected to perform similarly on variants where
allele frequency is not the dominant signal, and somewhat worse on variants
where gnomAD_AF provides decisive evidence.

## External validation

The primary model was evaluated on independent benchmarks including ClinVar
3-star holdout, Humsavar, MAVE, VIP, Grimm 2015, gnomAD benign, PanelApp
benign, DVD, and DDD clinical datasets. Full results are available in
`AGVarPred-zenodo/metrics/` and the reproducibility scripts in
`AGVarPred-training/scripts/external_validation/`.

## Known limitations

- Trained on ClinVar labels, which reflect current clinical knowledge and may
  contain biases or evolve over time.
- The full model requires a valid gnomAD VCF for allele frequency and VEP
  annotations. The no-AF fallback is provided for convenience but may be less
  accurate.
- AlphaGenome functional scores require either the AlphaGenome SDK/API or
  precomputed matrices.
- Mitochondrial variants are represented but were a small fraction of training
  data; interpret with caution.
- Performance may vary across populations due to representation biases in
  ClinVar and gnomAD.

## Ethical considerations

- Do not use the model to make deterministic clinical diagnoses.
- Be aware of population biases in ClinVar and gnomAD when applying the model
  to underrepresented groups.
- Predictions should be combined with expert review and orthogonal evidence.

## Software and reproducibility

- Source code: <https://github.com/pankaj357/AGVarPred>
- Zenodo archive: [`10.5281/zenodo.20955782`](https://doi.org/10.5281/zenodo.20955782) (concept DOI; resolves to the latest version)
- License: Apache-2.0

## References

- Landrum, M. J., et al. (2018). ClinVar: improving access to variant
  interpretations and supporting evidence. *Nucleic Acids Research*, 46(D1),
  D1062–D1067. https://doi.org/10.1093/nar/gkx1153
- Karczewski, K. J., et al. (2020). The mutational constraint spectrum quantified
  from variation in 141,456 humans. *Nature*, 581, 434–443.
  https://doi.org/10.1038/s41586-020-2308-7
- McLaren, W., et al. (2016). The Ensembl Variant Effect Predictor. *Genome
  Biology*, 17, 122. https://doi.org/10.1186/s13059-016-0974-4
- Ke, G., et al. (2017). LightGBM: A highly efficient gradient boosting decision
  tree. *NeurIPS*. https://papers.nips.cc/paper/6907-lightgbm
- AlphaGenome functional genomics platform (see `docs/installation.md` for
  access instructions).

## Acknowledgements

We acknowledge the ICAR-Indian Institute of Agricultural Biotechnology,
Ranchi for supporting this work.
