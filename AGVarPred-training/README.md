# AGVarPred-training

This directory contains the complete research pipeline used to train and evaluate AGVarPred.

## Purpose

- Reproduce the manuscript from raw ClinVar data through model training.
- Document exact data-cleaning, feature-extraction, and model-selection steps.
- Host external-validation benchmark scripts.

## Requirements

Install the package and training dependencies from the repository root:

```bash
cd ..
pip install -e ".[test]"
pip install -r AGVarPred-training/requirements.txt
```

You also need:

- ClinVar `variant_summary.txt` (place at the repository root).
- gnomAD exomes r2.1.1 liftover to GRCh38 VCF (set `GNOMAD_VCF`).
- AlphaGenome Python SDK and API keys (for feature extraction).

## Reproduction workflow

```bash
cd AGVarPred-training
make features
make train
make evaluate
make benchmark
```

## Directory layout

```
scripts/
├── data_preparation/      # Build ClinVar gold set and splits
├── feature_extraction/    # Chunking, AlphaGenome scoring, feature assembly
├── feature_selection/     # Recursive feature elimination / importance selection
├── model_training/        # Baselines, XGBoost, final LightGBM models
├── ablation/              # Ablation studies
└── external_validation/   # Benchmark scoring and comparison scripts
```

## Notes

- The feature-extraction step is the slowest; it is designed to run in 20
  parallel jobs, one per `RUN_ID=0..19`.
- Intermediate folders (`feature_extraction_*`, `chunks_input_*`, `logs_*`,
  `final_dataset_parts_*`) are excluded from Git and can be regenerated with
  `make features`.
