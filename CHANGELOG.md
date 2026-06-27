# Changelog

All notable changes to AGVarPred will be documented in this file.

## [1.0.0] - 2026-06-26

### Added
- Initial public release of AGVarPred.
- Command-line interface: `AGVarPred predict` and `AGVarPred list-models`.
- Manifest-driven, versioned model directory with SHA256 checksums.
- Internal `agvarpred_core` library for reusable feature generation.
- Hybrid AlphaGenome support: SDK or precomputed feature matrices.
- Example VCF and precomputed features for offline testing.
- GitHub Actions CI/CD workflows for tests and releases.
- Reproducibility package (`AGVarPred-training/`) and Zenodo archive
  (`AGVarPred-zenodo/`).

### Model
- `model_full`: primary production model (120 features, threshold 0.42).
- `model_no_af`: fallback model without gnomAD_AF (119 features, threshold 0.45).
- Both models are regularized LightGBM classifiers trained on ClinVar germline
  variants (GRCh38), calibrated with isotonic regression.
