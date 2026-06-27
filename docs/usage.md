# Usage

## Automatic model selection

By default `AGVarPred predict` selects the best available model:

1. If a local gnomAD VCF is available → `model_full`.
2. If an online AF/VEP provider is enabled → `model_full`.
3. Otherwise → `model_no_af` with a warning.

You can inspect the selected model and AF source in the output CSV metadata.

## Command-line interface

### Predict pathogenicity for a VCF

```bash
# Auto-select model (recommended)
AGVarPred predict input.vcf -o predictions.csv --alpha-mode precomputed --alpha-dir features/

# Provide a local gnomAD VCF to ensure the full model is used
AGVarPred predict input.vcf -o predictions.csv --gnomad-vcf gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz

# Force a specific model
AGVarPred predict input.vcf -o predictions.csv --model full
AGVarPred predict input.vcf -o predictions.csv --model no_af

# Use the AlphaGenome SDK
AGVarPred predict input.vcf -o predictions.csv --alpha-mode sdk

# Show the installed version
AGVarPred --version
```

### List installed models

```bash
AGVarPred list-models
```

## Output format

`predictions.csv` contains one row per input variant:

| Column                | Description                                      |
|-----------------------|--------------------------------------------------|
| variant_id            | `chr_pos_ref_alt` identifier                     |
| probability           | Calibrated probability of pathogenicity          |
| predicted_class       | 0 = benign, 1 = pathogenic                       |
| model_version         | Model version used                               |
| model_type            | `full` or `no_AF`                                |
| af_source             | `local_gnomad`, `online`, or `none`              |
| feature_version       | Feature-engineering version                      |
| alphagenome_version   | AlphaGenome SDK version                          |
| gnomAD_version        | gnomAD version used for AF/VEP                   |
| VEP_version           | VEP version used for annotations                 |
| prediction_timestamp  | ISO 8601 timestamp                               |

## Python API

```python
from agvarpred_core.feature_generator import FeatureGenerator
from AGVarPred import AGVarPredAutoPredictor

auto = AGVarPredAutoPredictor(gnomad_vcf="gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz")
generator = FeatureGenerator(
    af_source=auto.af_source,
    alpha_mode="precomputed",
    alpha_dir="alpha_features/",
)
features = generator.from_vcf("sample.vcf")
predictions = auto.predict(features)
```
