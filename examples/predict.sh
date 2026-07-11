#!/usr/bin/env bash
# Example: predict pathogenicity for sample.vcf using precomputed AlphaGenome features.
# Requires the gnomAD VCF to be available locally.

set -euo pipefail

GNOMAD_VCF="${GNOMAD_VCF:-external_data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz}"

AGVarPred predict sample.vcf \
  --output predictions.csv \
  --alpha-mode precomputed \
  --alpha-dir alpha_features/ \
  --gnomad-vcf "$GNOMAD_VCF"

echo "Predictions written to predictions.csv"
