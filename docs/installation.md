# Installation

## Quick install

```bash
git clone https://github.com/pankaj357/AGVarPred.git
cd AGVarPred
pip install -e .
```

This installs the `AGVarPred` command-line tool and the `agvarpred_core`
feature-engineering library.

## Dependencies

Core runtime dependencies are installed automatically:

- pandas, numpy, pyyaml, joblib
- scikit-learn >= 1.8.0, < 1.10.0
- lightgbm
- pysam

## External data

### AlphaGenome features (choose one)

#### Option A: AlphaGenome SDK (recommended for large batches)

Install the AlphaGenome SDK from its official source and obtain an API key:

```bash
export ALPHAGENOME_API_KEY=your_key_here
AGVarPred predict input.vcf -o predictions.csv
```

> **Publication note:** The AlphaGenome SDK and API access instructions will be
> linked here once publicly available. Until then, users who do not have SDK
> access should use precomputed feature matrices (Option B).

#### Option B: Precomputed feature matrices

If you already have AlphaGenome scores for your variants, organize them as:

```
alpha_features/
├── GENE1_VEP/
│   └── GENE1_ALL_VEP_RAW_SCORE_MATRIX.parquet
└── GENE2_VEP/
    └── GENE2_ALL_VEP_RAW_SCORE_MATRIX.parquet
```

Then run:

```bash
AGVarPred predict input.vcf -o predictions.csv --alpha-mode precomputed --alpha-dir alpha_features/
```

### gnomAD VCF (recommended but optional)

A local gnomAD VCF lets AGVarPred run the full production model (`model_full`).
If it is not available, the software automatically falls back to the bundled
`model_no_af` model.

Download the gnomAD exomes r2.1.1 liftover to GRCh38 VCF:

```bash
mkdir -p external_data
cd external_data
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi
```

Set the environment variable:

```bash
export GNOMAD_VCF=/path/to/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
```

## Development install

```bash
pip install -e ".[test]"
pytest
```
