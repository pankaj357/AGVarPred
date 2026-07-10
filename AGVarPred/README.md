# AGVarPred

<!-- Badges -->
[![Tests](https://github.com/pankaj357/AGVarPred/actions/workflows/tests.yml/badge.svg)](https://github.com/pankaj357/AGVarPred/actions/workflows/tests.yml)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20955782.svg)](https://doi.org/10.5281/zenodo.20955782)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**AGVarPred** is a command-line tool and Python library for predicting the pathogenicity of germline variants from a VCF file. It integrates:

- **AlphaGenome** functional genomic scores (epigenomic, transcriptomic, and splicing signals),
- **gnomAD** allele frequency,
- **VEP** annotations (SIFT, PolyPhen, IMPACT, LoF, consequence),

into a calibrated LightGBM classifier trained on high-confidence ClinVar-labeled variants.

## Research motivation

Variant interpretation remains one of the central challenges in clinical
genomics. Existing in-silico predictors often rely on a small set of
conservation and biochemical features. AGVarPred adds a large-scale functional
genomics layer by leveraging AlphaGenome, which scores variants across hundreds
of cell-type-specific regulatory and transcriptomic assays. The resulting model
achieves strong discrimination between pathogenic and benign germline variants
and is fully reproducible.

## Key features

- **One-command prediction** from a VCF file.
- **Automatic model selection**: uses the full model when a gnomAD VCF is
  available and transparently falls back to a bundled no-AF model otherwise.
- **Hybrid AlphaGenome support**: query the AlphaGenome API or load precomputed
  feature matrices.
- **Reproducible research pipeline**: training, feature selection, ablation, and
  external validation scripts are included.
- **Model versioning**: multiple bundled models with manifest files and SHA256
  checksums.
- **Prediction metadata**: every output row records model version, AF source,
  and timestamp.

## Installation

```bash
git clone https://github.com/pankaj357/AGVarPred.git
cd AGVarPred
pip install -e .
```

For development:

```bash
pip install -e ".[test]"
pytest
```

See [`docs/installation.md`](docs/installation.md) for detailed setup, including
AlphaGenome and gnomAD configuration.

## Quick start

```bash
# Predict with precomputed AlphaGenome features.
# If GNOMAD_VCF is set, the full model runs automatically.
# Otherwise AGVarPred falls back to the bundled no-AF model.
AGVarPred predict sample.vcf \
  --output predictions.csv \
  --alpha-mode precomputed \
  --alpha-dir /path/to/alpha_features
```

A minimal example is included in the repository:

```bash
AGVarPred predict examples/sample.vcf \
  --output /tmp/predictions.csv \
  --alpha-mode precomputed \
  --alpha-dir examples/alpha_features/
```

## Automatic model selection

AGVarPred chooses the best available model for each run:

1. **Full model** (`model_full`) — used when a local gnomAD VCF is available
   via `--gnomad-vcf` or the `GNOMAD_VCF` environment variable.
2. **Online annotation** — reserved for a future provider interface; not yet
   implemented.
3. **No-AF fallback** (`model_no_af`) — used automatically when no AF/VEP
   source is available. A warning is printed and the output metadata records
   `model_type=no_AF` and `af_source=none`.

Advanced users can override automatic selection with `--model full` or
`--model no_af`.

## CLI examples

```bash
# Auto-select model (recommended)
AGVarPred predict input.vcf -o predictions.csv \
  --alpha-mode precomputed --alpha-dir alpha_features/

# Provide a local gnomAD VCF to ensure the full model is used
AGVarPred predict input.vcf -o predictions.csv \
  --gnomad-vcf gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz

# Force a specific model
AGVarPred predict input.vcf -o predictions.csv --model no_af \
  --alpha-mode precomputed --alpha-dir alpha_features/

# List installed models
AGVarPred list-models

# Show version
AGVarPred --version
```

## Python API example

```python
from agvarpred_core.feature_generator import FeatureGenerator
from AGVarPred import AGVarPredAutoPredictor

auto = AGVarPredAutoPredictor(gnomad_vcf="gnomad.vcf.bgz")
generator = FeatureGenerator(
    af_source=auto.af_source,
    alpha_mode="precomputed",
    alpha_dir="alpha_features/",
)
features = generator.from_vcf("sample.vcf")
predictions = auto.predict(features)
print(predictions[["variant_id", "probability", "predicted_class"]])
```

See [`docs/api.md`](docs/api.md) for the full API reference.

## Repository structure

```
AGVarPred/
├── src/
│   ├── AGVarPred/          # Public installable package (CLI + predictor)
│   └── agvarpred_core/     # Internal reusable feature-engineering library
├── model/                  # Trained model artifacts (separate from source)
│   ├── active_model.json
│   ├── model_full/         # Primary production model (with AF)
│   └── model_no_af/        # Fallback model (without AF)
├── tests/                  # Unit and integration tests
├── examples/               # Example VCF and precomputed features
├── docs/                   # Documentation
├── AGVarPred-training/     # Full research reproducibility pipeline
├── AGVarPred-zenodo/       # Manuscript archive (model, splits, benchmarks)
└── community files         # LICENSE, CITATION.cff, CONTRIBUTING.md, etc.
```

## Requirements

- Python >= 3.11
- Either the AlphaGenome Python SDK and a valid API key, or a directory of
  precomputed AlphaGenome feature matrices.
- A local gnomAD VCF is **recommended** but not required:
  - Download gnomAD exomes r2.1.1 liftover to GRCh38.
  - Set `GNOMAD_VCF=/path/to/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz`.
  - If absent, the bundled no-AF model runs instead.

## Citation

If you use AGVarPred in your research, please cite the Zenodo archive and the
GitHub repository.

- **Concept DOI** (always resolves to the latest version):  
  <https://doi.org/10.5281/zenodo.20955782>
- **Version DOI for v1.0.3**:  
  <https://doi.org/10.5281/zenodo.20955783>

```text
Kumar, P., & Kanaka, K. K. (2026). AGVarPred: germline variant pathogenicity
prediction using AlphaGenome functional genomics (Version 1.0.3).
Zenodo. https://doi.org/10.5281/zenodo.20955782
GitHub: https://github.com/pankaj357/AGVarPred
```

See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

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
- AlphaGenome functional genomics platform (see [`docs/installation.md`](docs/installation.md)
  for access instructions).

## License

Apache-2.0. See [`LICENSE`](LICENSE).

## Contact

For questions, bug reports, or feature requests, please open a GitHub issue at
<https://github.com/pankaj357/AGVarPred/issues> or contact the maintainers at
<ft.pank@gmail.com> or <kkokay07@gmail.com>.

## Acknowledgements

This work was built on top of ClinVar, gnomAD, Ensembl VEP, and the AlphaGenome
functional genomics platform. We thank the maintainers of these resources and
the broader variant-interpretation community.

We acknowledge the ICAR-Indian Institute of Agricultural Biotechnology,
Ranchi for supporting this work.
