# AGVarPred v1.0.1 Release Readiness Report

Generated: 2026-06-27

## Summary

The AGVarPred repository has been successfully published to GitHub at
`https://github.com/pankaj357/AGVarPred`. The `main` branch and tags
(`v1.0.0`, `v1.0.1`) are live. GitHub Actions now passes on Python 3.11,
3.12, and 3.13.

Because the initial `v1.0.0` tag pointed to a commit with failing CI
(Python 3.10 incompatibility and missing Zenodo pickle files in Git tracking),
a `v1.0.1` tag was created on the fixed commit. The recommended release is
**v1.0.1**.

No scientific methodology, machine learning models, feature engineering, or
prediction logic were modified.

## Repository state

| Item | Value |
|---|---|
| Repository URL | https://github.com/pankaj357/AGVarPred |
| Default branch | `main` |
| Latest commit (main) | `9e8774e2...` |
| Latest tag | `v1.0.1` |
| v1.0.0 tag | `b19c8e50...` (historical; CI fails) |
| v1.0.1 tag | `9e8774e2...` (recommended for release) |
| Working tree | Clean |

## What was published

- All source code (`src/AGVarPred/`, `src/agvarpred_core/`)
- Model artifacts (`model/`, `AGVarPred-zenodo/model/`)
- Documentation (`README.md`, `MODEL_CARD.md`, `docs/`)
- Tests (`tests/`)
- Examples (`examples/`)
- Training reproducibility package (`AGVarPred-training/`)
- Zenodo archive package (`AGVarPred-zenodo/`)
- GitHub Actions workflows (`.github/workflows/`)
- Community files (`LICENSE`, `CITATION.cff`, `CONTRIBUTING.md`, etc.)

## CI compatibility fixes applied after initial push

| Fix | Reason |
|---|---|
| `requires-python` raised to `>=3.11` | scikit-learn 1.8.0+ requires Python ≥3.11 |
| CI matrix updated to 3.11/3.12/3.13 | Matches supported Python versions |
| `pyarrow>=10.0.0` added to runtime deps | Tests read Parquet feature matrices |
| `AGVarPred-zenodo/model/**/*.pkl` added to Git | Zenodo pickle files were excluded by `.gitignore` |

## Validation results

| Check | Result |
|---|---|
| `pytest tests/` (local) | 21 passed, 1 skipped |
| `python -m build` (local) | Success |
| Model pipeline SHA256 | Both models match manifests |
| Zenodo `sha256sum -c checksums.sha256` | All OK |
| GitHub Actions `tests` workflow | Success on 3.11, 3.12, 3.13 |
| GitHub repository size | ~132 MB working tree, ~26 MB `.git` |

## Remaining placeholders requiring owner action

1. **Zenodo DOI** (`10.5281/zenodo.TODO`) appears in:
   - `CITATION.cff`
   - `AGVarPred-training/CITATION.cff`
   - `AGVarPred-zenodo/CITATION.cff`
   - `README.md` (badge + citation)
   - `MODEL_CARD.md`

   Action: Create Zenodo record from GitHub Release v1.0.1, obtain DOI,
   replace all occurrences, regenerate `AGVarPred-zenodo/checksums.sha256`.

2. **AlphaGenome SDK public instructions**
   - `docs/installation.md` has a placeholder note pending public SDK docs.

   Action: Add official AlphaGenome SDK URL/install command when available.

## GitHub Release preparation

### Recommended release tag

`v1.0.1` (commit `9e8774e2...`)

### Release title

AGVarPred v1.0.1

### Release notes

```markdown
## AGVarPred v1.0.1

Initial public release of AGVarPred, a command-line tool and Python library for
germline variant pathogenicity prediction using AlphaGenome functional genomics
scores, gnomAD allele frequency, and VEP annotations.

### Included
- `model_full`: primary production model (120 features)
- `model_no_af`: fallback model without gnomAD AF (119 features)
- Command-line interface (`AGVarPred predict`, `AGVarPred list-models`)
- Python API (`AGVarPredAutoPredictor`, `AGVarPredPredictor`)
- Example VCF and precomputed AlphaGenome features
- Reproducibility scripts (`AGVarPred-training/`)
- Zenodo archive package (`AGVarPred-zenodo/`)

### Requirements
- Python >= 3.11
- Optional: local gnomAD exomes r2.1.1 liftover to GRCh38 VCF
- Optional: AlphaGenome SDK/API key or precomputed feature matrices

### Notes
- This release corrects CI compatibility (Python 3.11+) and ensures all model
  pickle files are present in the Zenodo archive.
- The Zenodo DOI badge and `CITATION.cff` files currently contain a placeholder
  (`10.5281/zenodo.TODO`) that will be updated after the Zenodo record is created.
```

### Assets to upload

- `dist/agvarpred-1.0.0-py3-none-any.whl`
- `dist/agvarpred-1.0.0.tar.gz`

Generate fresh artifacts from the `v1.0.1` tag with:

```bash
git checkout v1.0.1
python -m build
```

## Next steps

1. On GitHub, create a Release from tag `v1.0.1`.
2. Upload the wheel and sdist from `dist/`.
3. Link the repository to Zenodo and publish the Zenodo record.
4. Replace `10.5281/zenodo.TODO` with the real DOI everywhere.
5. Regenerate `AGVarPred-zenodo/checksums.sha256`.
6. Optionally update `docs/installation.md` with official AlphaGenome SDK docs.
