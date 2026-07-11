# AGVarPred v1.0.8 Release Summary

**Release date:** 2026-07-11  
**GitHub Release:** https://github.com/pankaj357/AGVarPred/releases/tag/v1.0.8  
**Git tag:** `v1.0.8`  
**Package version:** `1.0.8`

---

## What was published

### Source code
- Git commit `7dff05b` on `main`.
- Git tag `v1.0.8` (annotated).

### Release assets
| Asset | URL | Size (local) |
|---|---|---|
| Wheel | https://github.com/pankaj357/AGVarPred/releases/download/v1.0.8/agvarpred-1.0.8-py3-none-any.whl | ~3.2 MB |
| Source distribution | https://github.com/pankaj357/AGVarPred/releases/download/v1.0.8/agvarpred-1.0.8.tar.gz | ~3.2 MB |

Both assets were built with `python -m build` in a clean environment.

---

## Changes from v1.0.7

### Repository cleanup
- Removed all manuscript drafting materials from the public software repository:
  - Entire `publication_readiness/` directory (markdown manuscripts, audit reports, planning files).
  - `dca_analysis/manuscript_discussion.txt`
  - `dca_analysis/manuscript_results.txt`
  - `delong_analysis/manuscript_results.txt`
  - `tissue_shap_analysis/manuscript_discussion.txt`
  - `tissue_shap_analysis/manuscript_results.txt`
- Moved reusable correction scripts to `AGVarPred/scripts/manuscript_corrections/`:
  - `generate_figure2.py`
  - `generate_readiness_report.py`
  - `regenerate_test_predictions_with_variant_id.py`
- Added `.gitignore` patterns to prevent re-commit of manuscript materials.
- A local backup of the removed manuscript files is kept in `private_manuscript_backup/` (gitignored).

### Version synchronization
All version references were updated to `1.0.8`:

- `AGVarPred/pyproject.toml`
- `AGVarPred/src/AGVarPred/__init__.py`
- `AGVarPred/CITATION.cff` (including release date)
- `AGVarPred/MODEL_CARD.md`
- `AGVarPred/README.md` citation
- `AGVarPred/src/AGVarPred/model/model_full/manifest.yaml`
- `AGVarPred/src/AGVarPred/model/model_no_af/manifest.yaml`
- `AGVarPred/AGVarPred-zenodo/model/model_full/manifest.yaml`
- `AGVarPred/AGVarPred-zenodo/model/model_no_af/manifest.yaml`
- `AGVarPred/CHANGELOG.md`

### Packaging restoration
- Built clean wheel and source distribution with `python -m build`.
- Uploaded both artifacts to the GitHub Release.

---

## Verification performed

### 1. Build verification
```bash
cd AGVarPred
rm -rf build dist *.egg-info
python -m build
```
Produced:
- `dist/agvarpred-1.0.8-py3-none-any.whl`
- `dist/agvarpred-1.0.8.tar.gz`

### 2. Content inspection
Both the wheel and the source distribution were inspected. Neither contains:
- `publication_readiness/`
- audit/forensic/review reports
- manuscript markdowns
- `manuscript_results.txt` / `manuscript_discussion.txt`
- supplementary drafting files

Both contain only:
- `AGVarPred/` package source and bundled model files
- `agvarpred_core/` package source
- `tests/`
- `README.md`, `LICENSE`, `pyproject.toml`, metadata

### 3. GitHub source archive inspection
Downloaded `https://github.com/pankaj357/AGVarPred/archive/refs/tags/v1.0.8.tar.gz`.  
Confirmed: **no manuscript materials** are present in the source archive.

### 4. Installation tests
- Fresh venv + `pip install dist/agvarpred-1.0.8-py3-none-any.whl` → OK.
- Fresh venv + `pip install dist/agvarpred-1.0.8.tar.gz` → OK.
- Fresh venv + `pip install https://github.com/pankaj357/AGVarPred/releases/download/v1.0.8/agvarpred-1.0.8-py3-none-any.whl` → OK.
- `AGVarPred list-models` → reports `model_full` and `model_no_af`.
- `python -c "import AGVarPred; print(AGVarPred.__version__)"` → prints `1.0.8`.

### 5. Git verification
- Branch: `main`
- Commit: `7dff05b`
- Tag: `v1.0.8`
- Push: successful
- GitHub Release: created with both assets uploaded.

---

## Zenodo status

GitHub→Zenodo auto-sync was left enabled. Because the repository no longer contains manuscript materials, the Zenodo archive generated from the `v1.0.8` release should be software-only. The Zenodo record had not yet appeared at the time of this summary; verify it after the auto-sync completes.

---

## Known limitations / notes

1. **CLI `--version` flag is not implemented.** `AGVarPred --version` currently errors because subcommands are required. Version can be checked via `python -c "import AGVarPred; print(AGVarPred.__version__)"` or `AGVarPred list-models`.
2. **Package name is capitalized.** `from agvarpred import *` does not work; use `from AGVarPred import *` or `import AGVarPred`.
3. **AlphaGenome SDK** remains a separate, optional dependency and is not included in the wheel.
4. **README Version DOI placeholder:** the line "Version DOI for v1.0.8" is marked as "to be added after Zenodo publication." Update it once the Zenodo DOI for v1.0.8 is known.

---

## Installation command

```bash
pip install https://github.com/pankaj357/AGVarPred/releases/download/v1.0.8/agvarpred-1.0.8-py3-none-any.whl
```

---

## Files changed in this release

- `AGVarPred/pyproject.toml`
- `AGVarPred/src/AGVarPred/__init__.py`
- `AGVarPred/CITATION.cff`
- `AGVarPred/MODEL_CARD.md`
- `AGVarPred/README.md`
- `AGVarPred/CHANGELOG.md`
- `AGVarPred/src/AGVarPred/model/model_full/manifest.yaml`
- `AGVarPred/src/AGVarPred/model/model_no_af/manifest.yaml`
- `AGVarPred/AGVarPred-zenodo/model/model_full/manifest.yaml`
- `AGVarPred/AGVarPred-zenodo/model/model_no_af/manifest.yaml`
- `.gitignore`
- Removed: `publication_readiness/` and analysis `manuscript_*.txt` files
- Added: `AGVarPred/scripts/manuscript_corrections/`
- Added: `AGVarPred/RELEASE_SUMMARY.md`
