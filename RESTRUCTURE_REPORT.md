# Repository Restructure Report

## Objective

Restructure the `AGVarPred` repository so that the installable Python package
lives at the repository root instead of under a nested `AGVarPred/` directory.
This brings the project in line with the standard setuptools `src`-layout and
lets GitHub display `README.md`, run Actions, and build the package from the
repository root.

No functional source code was changed, the version remains **1.0.8**, and no
new release, tag, or Zenodo entry was created.

## Layout before and after

### Before

```
/data/kanaka/pankaj/feature_extration/new/   <-- repository root
├── .gitignore                               <-- workspace ignore rules
├── AGVarPred/                               <-- nested package directory
│   ├── .github/workflows/
│   ├── AGVarPred-training/
│   ├── AGVarPred-zenodo/
│   ├── CHANGELOG.md
│   ├── CITATION.cff
│   ├── CODE_OF_CONDUCT.md
│   ├── CONTRIBUTING.md
│   ├── LICENSE
│   ├── MODEL_CARD.md
│   ├── PACKAGING_AUDIT.md
│   ├── README.md
│   ├── RELEASE_*.md
│   ├── SECURITY.md
│   ├── data_manifest.json
│   ├── docs/
│   ├── examples/
│   ├── pyproject.toml
│   ├── requirements.txt
│   ├── scripts/
│   ├── src/
│   ├── tests/
│   └── train_gene_set.json
└── dca_analysis/ ...
```

### After

```
/data/kanaka/pankaj/feature_extration/new/   <-- repository root
├── .github/workflows/
├── AGVarPred-training/
├── AGVarPred-zenodo/
├── CHANGELOG.md
├── CITATION.cff
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
├── LICENSE
├── MODEL_CARD.md
├── PACKAGING_AUDIT.md
├── README.md
├── RELEASE_*.md
├── SECURITY.md
├── data_manifest.json
├── docs/
├── examples/
├── pyproject.toml
├── requirements.txt
├── scripts/
├── src/
├── tests/
├── train_gene_set.json
└── dca_analysis/ ...
```

## What changed

### 1. Moved tracked files to the repository root

All tracked contents of the nested `AGVarPred/` directory were promoted one
level using `git mv` so that Git records them as renames rather than
delete/add operations. This preserves file-level history.

Moved items:

- `.github/`
- `AGVarPred-training/`
- `AGVarPred-zenodo/`
- `docs/`
- `examples/`
- `scripts/`
- `src/`
- `tests/`
- Community/markdown files: `CHANGELOG.md`, `CITATION.cff`,
  `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, `LICENSE`, `MODEL_CARD.md`,
  `PACKAGING_AUDIT.md`, `README.md`, `RELEASE_AUDIT.md`,
  `RELEASE_CHECKLIST.md`, `RELEASE_CHECKLIST_v1.0.8.md`, `RELEASE_REPORT.md`,
  `RELEASE_SUMMARY.md`, `SECURITY.md`
- Project metadata: `data_manifest.json`, `pyproject.toml`,
  `requirements.txt`, `train_gene_set.json`

The now-empty `AGVarPred/` directory was removed.

### 2. Updated `.gitignore`

The root `.gitignore` was rewritten to:

- Track the new top-level project files and directories.
- Keep the package-oriented ignore rules from the old `AGVarPred/.gitignore`
  (e.g., model `.pkl` unignore rules, build artifacts, virtual environments).
- Continue ignoring large workspace artifacts (e.g., `alphagenome_input/`,
  `chunks_input_*/`, `feature_extraction_*_dataset/`, `final_dataset_parts_*/`,
  `external_data/`, `output/`, `external_validation/`, `final_model_output/`,
  `feature_selection_output_nested/`, `xgboost_output/`).
- Ignore top-level scratch duplicates that are not part of the package
  (`/*.py`, `/*.csv`, `/*.sh`, `/*.docx`, `/*.log`).

### 3. Updated documentation paths

- `README.md`: refreshed the repository structure diagram so it no longer
  shows a wrapping `AGVarPred/` directory.
- `CHANGELOG.md`: changed the manuscript-corrections path reference from
  `AGVarPred/scripts/manuscript_corrections/` to
  `scripts/manuscript_corrections/`.

### 4. Fixed CI workflow

`.github/workflows/tests.yml` previously used a hard-coded `Path("model")`
that did not point to the bundled models. The manifest-verification step now
uses `AGVarPred.pipeline.get_model_root()` so it correctly resolves the
bundled model directory in both editable and installed contexts.

## Verification

| Check | Command | Result |
|---|---|---|
| Build produces wheel + sdist | `python -m build` | ✅ `agvarpred-1.0.8-py3-none-any.whl` and `agvarpred-1.0.8.tar.gz` |
| Install from wheel | `pip install dist/agvarpred-1.0.8-py3-none-any.whl` | ✅ Installed successfully |
| Install editable | `pip install -e ".[test]"` | ✅ Installed successfully |
| CLI model discovery | `AGVarPred list-models` | ✅ Shows `model_full` and `model_no_af` |
| Manifest validation | `python -c "from AGVarPred.pipeline import ..."` | ✅ Default/full/no_af mapping correct; both manifests validate |
| Unit tests | `pytest tests/ -v` | ✅ 20 passed, 1 skipped |

## Git history

The restructure was committed as a single Git commit with `git mv` renames so
that `git log --follow` continues to work for each moved file. No history was
rewritten or lost.

## Notes

- Version remains `1.0.8` in `pyproject.toml`, `src/AGVarPred/__init__.py`,
  `CITATION.cff`, `MODEL_CARD.md`, and other manifests.
- No new Git tag or GitHub Release was created, so the existing Zenodo
  integration was not triggered.
- Functional source code was not modified; only paths, the `.gitignore`, the
  README tree, a CHANGELOG reference, and the CI manifest-check step were
  updated.
