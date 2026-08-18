# Changelog

All notable changes to AGVarPred will be documented in this file.

## [1.0.9] - 2026-08-18

### Fixed
- Regenerated `AGVarPred-zenodo/checksums.sha256` so the frozen archive passes
  self-verification again after the manifest metadata updates.
- Restored tracking of `AGVarPred-training/scripts/external_validation/`
  (benchmark construction and scoring scripts), `examples/sample.vcf`,
  `examples/alpha_features/`, per-benchmark result JSONs under
  `external_validation/results/`, and analysis figure sources; these were
  unintentionally excluded by overly broad `.gitignore` rules, which broke the
  README quick-start and the CI fallback test on fresh clones.
- Synchronized the root `data_manifest.json` with the Zenodo copy.

### Changed
- Bumped version to `1.0.9` across the package, model manifests, and citation
  files. No functional code changes.

## [1.0.8] - 2026-07-10

### Changed
- Synchronized all version references to `1.0.8` (`pyproject.toml`, `__init__.py`,
  `CITATION.cff`, model manifests, README citation).
- Removed manuscript drafting materials (`publication_readiness/` markdowns and
  `manuscript_*.txt` files) from the public software repository.
- Preserved reusable correction scripts by moving them to
  `scripts/manuscript_corrections/`.
- Restored proper Python packaging: releases now include wheel + source
  distribution assets built with `python -m build`.

### Fixed
- GitHub Releases v1.0.5–v1.0.7 no longer represent the latest installable
  artifacts; v1.0.8 provides the canonical wheel and sdist.

## [1.0.4] - 2026-06-27

### Fixed
- Use `importlib.resources` for bundled model discovery so `AGVarPred list-models`
  and `AGVarPred predict` work after a non-editable `pip install` or wheel install.
- `sha256_file()` and manifest helpers now accept `Traversable` package resources.
- Tests updated to discover the bundled model directory via `importlib.resources`.

## [1.0.3] - 2026-06-27

### Fixed
- Bundle model files (`model_full/` and `model_no_af/`) inside the `AGVarPred`
  package so they are included in wheel/sdist installs.
- Declare `src/AGVarPred/model/**/*` as package data in `pyproject.toml`.

## [1.0.2] - 2026-06-27

### Changed
- Replaced the Zenodo DOI placeholder in documentation and badges with the
  official DOI: `10.5281/zenodo.20955782`.

## [1.0.1] - 2026-06-27

### Fixed
- Include Zenodo model pickle files in Git tracking (LFS pointer resolution).
- Add `pyarrow>=10.0.0` as a runtime dependency for parquet support.
- Require Python `>=3.11` and update CI compatibility matrix.

### Changed
- Updated release checklist and reproducibility report to reflect v1.0.1 status.

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
