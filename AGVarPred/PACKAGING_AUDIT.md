# AGVarPred Packaging Audit

**Audit date:** 2026-07-10  
**Scope:** Python packaging configuration, artifact contents, and manuscript-material leakage risk.

---

## 1. Executive Summary

AGVarPred is configured as a standard `src`-layout setuptools project driven by `pyproject.toml`. The **wheel** and **sdist** produced by `python -m build` are clean: they contain only source code, bundled model files, tests, `LICENSE`, `README.md`, and `pyproject.toml`. They do **not** include `publication_readiness/`, audit reports, or manuscript text files.

The **GitHub source archive** (and therefore the Zenodo archive) is a plain snapshot of the repository at the release tag. Because manuscript drafting directories are committed to the repository, the source archives currently ship manuscript materials. This is the root cause of the Zenodo leakage.

**No `MANIFEST.in`, `setup.py`, or `.gitattributes` file exists**, so artifact contents are controlled entirely by `pyproject.toml`, `tool.setuptools.packages.find`, `tool.setuptools.package-data`, and `.gitignore`.

---

## 2. Build Configuration Audit

### 2.1 `pyproject.toml`

- **Build backend:** `setuptools.build_meta`
- **Build-system requires:** `setuptools>=77.0.0`, `wheel`
- **Package name:** `AGVarPred`
- **Current version:** `1.0.4` (mismatch with GitHub release tag `v1.0.7`)
- **Layout:** `src/` layout (`src/AGVarPred/`, `src/agvarpred_core/`)
- **Package discovery:** `tool.setuptools.packages.find.where = ["src"]`
- **Package data:** `AGVarPred = ["py.typed", "model/**/*"]`; `agvarpred_core = ["py.typed"]`
- **Entry point:** `AGVarPred = "AGVarPred.cli:main"`
- **Dependencies:** pinned only by lower bounds; no lock file referenced in packaging metadata

### 2.2 `MANIFEST.in`

**Not present.** Setuptools falls back to default sdist inclusion rules plus `pyproject.toml` metadata. Default behavior includes:

- `README.md` (because `readme = "README.md"`)
- `LICENSE` (because `license = "Apache-2.0"`)
- `pyproject.toml`
- All Python packages discovered under `src/`
- Tests directory (`tests/`) is included because it is at the repository root and setuptools’ default sdist collector picks it up

### 2.3 `setup.py`

**Not present.** Build is fully PEP 518/PEP 621 via `pyproject.toml`.

### 2.4 `.gitignore`

- **Top-level `.gitignore`:** excludes large data/model directories and top-level scratch files. It does **not** exclude `publication_readiness/`, `dca_analysis/`, `delong_analysis/`, or `tissue_shap_analysis/`.
- **`AGVarPred/.gitignore`:** excludes build artifacts, virtual environments, logs, `*.pkl` (with exceptions for model `.pkl` files), and large external data files. It does **not** exclude manuscript text files.

### 2.5 `.gitattributes`

**Not present.** GitHub source archives are produced with default line endings and no export-subst/exclusion configuration.

### 2.6 Package Layout

```
AGVarPred/
├── pyproject.toml
├── LICENSE
├── README.md
├── CHANGELOG.md
├── CITATION.cff
├── RELEASE_REPORT.md
├── RELEASE_CHECKLIST.md
├── MODEL_CARD.md
├── docs/
├── src/
│   ├── AGVarPred/          # main package
│   └── agvarpred_core/     # core library
├── tests/
└── AGVarPred-zenodo/       # zenodo payload (not part of wheel/sdist)
```

`AGVarPred-zenodo/` is **not** included in wheel/sdist because it is outside the discovered `src/` packages. It is, however, included in the GitHub source archive.

---

## 3. Artifact Comparison

| Artifact | What it is | How it is produced | What it currently contains | Manuscript leakage? |
|---|---|---|---|---|
| **Wheel (`.whl`)** | Binary/installable distribution | `python -m build` → `bdist_wheel` | `AGVarPred/`, `agvarpred_core/`, model data, `dist-info` metadata, `LICENSE` | **No** |
| **Source distribution (`.tar.gz`)** | Source distribution for `pip install` | `python -m build` → `sdist` | `src/AGVarPred/`, `src/agvarpred_core/`, `tests/`, `README.md`, `LICENSE`, `pyproject.toml`, `setup.cfg`, egg-info | **No** |
| **GitHub source archive** | Plain git snapshot at tag | GitHub auto-generates from tag | Everything tracked in the repo at the tag: code, docs, `AGVarPred-zenodo/`, `publication_readiness/`, analysis folders, manuscript `.txt` files | **Yes** |
| **Zenodo archive** | Mirror of GitHub release | Zenodo auto-imports GitHub release | Same as GitHub source archive | **Yes** |

### 3.1 Wheel contents (verified build of current `main`)

Built with `python -m build --outdir /tmp/agvarpred_dist`:

```
agvarpred-1.0.4-py3-none-any.whl
├── AGVarPred/                  # package source
│   ├── model/                  # bundled models (JSON, YAML, CSV, .pkl)
│   └── ...
├── agvarpred_core/             # core library
├── agvarpred-1.0.4.dist-info/  # metadata + LICENSE
└── (no publication_readiness/, no manuscript files)
```

### 3.2 Source distribution contents (verified build)

```
agvarpred-1.0.4.tar.gz
├── README.md
├── LICENSE
├── pyproject.toml
├── setup.cfg
├── src/AGVarPred/
├── src/agvarpred_core/
└── tests/
```

### 3.3 GitHub source archive contents (verified download of `v1.0.7`)

```
AGVarPred-1.0.7.tar.gz
├── AGVarPred/                  # package code
├── publication_readiness/      # ❌ manuscript markdowns & audit reports
├── dca_analysis/manuscript_*.txt
├── delong_analysis/manuscript_results.txt
├── tissue_shap_analysis/manuscript_*.txt
└── ...
```

---

## 4. Findings

1. **Wheel and sdist are already clean.** The build system correctly excludes top-level manuscript directories because they are outside the discovered packages.
2. **GitHub/Zenodo leakage is caused by repository content, not packaging configuration.** The source archives are faithful git snapshots.
3. **Version drift:** `pyproject.toml` and `src/AGVarPred/__init__.py` report `1.0.4`; `CITATION.cff` reports `1.0.3`; GitHub release tags are at `v1.0.7`.
4. **Missing proper release assets:** `v1.0.5`, `v1.0.6`, and `v1.0.7` were created without uploading the wheel/sdist artifacts produced by `python -m build`.
5. **No `MANIFEST.in` is required** for the current clean wheel/sdist, but adding one would make inclusion/exclusion explicit and reviewer-friendly.
6. **`AGVarPred-zenodo/` is shipped in the GitHub source archive** but is not part of the installable package. This is acceptable if Zenodo is meant to host the full reproducibility payload, but it should be reviewed to ensure it contains no manuscript raw data.

---

## 5. Recommendations

### Short-term (for v1.0.8)

1. **Remove manuscript drafting material from the git-tracked repository** (or move it to a private/non-release branch).
2. **Preserve reusable scripts** by moving them out of `publication_readiness/` into an appropriate code directory before excluding the folder.
3. **Add explicit `.gitignore` patterns** for remaining manuscript text files:
   - `publication_readiness/` (after moving scripts)
   - `*/manuscript_results.txt`
   - `*/manuscript_discussion.txt`
4. **Bump all version strings** to `1.0.8` (`pyproject.toml`, `src/AGVarPred/__init__.py`, `CITATION.cff`).
5. **Upload wheel + sdist to the GitHub release**; do not rely solely on the auto-generated source archive.
6. **For Zenodo:** either (a) disable GitHub auto-sync and manually upload the clean wheel + sdist, or (b) keep auto-sync after the repository is cleaned and verify the archived source archive contains no manuscript materials.

### Long-term

1. Add a `MANIFEST.in` to document intentional inclusions, e.g.:
   ```
   include LICENSE
   include README.md
   include pyproject.toml
   recursive-include src *
   recursive-include tests *
   ```
2. Automate release asset upload via GitHub Actions.
3. Publish to PyPI so users can `pip install agvarpred`.
4. Keep manuscript materials in a separate, private repository; link to it from the software README if needed.
