# Clean v1.0.8 Release Checklist — AGVarPred

**Goal:** Publish a software-only GitHub release and Zenodo archive that contains no manuscript drafting materials, includes proper wheel + sdist assets, and aligns all version strings.

---

## Phase 1 — Repository cleanup (no release yet)

- [ ] **Back up manuscript materials**
  - Copy the entire `publication_readiness/` directory to a safe location outside the git repository (e.g., local backup or a private `AGVarPred-manuscript` repo).
  - Copy `dca_analysis/manuscript_discussion.txt` and `manuscript_results.txt` to the backup.
  - Copy `delong_analysis/manuscript_results.txt` to the backup.
  - Copy `tissue_shap_analysis/manuscript_discussion.txt` and `manuscript_results.txt` to the backup.

- [ ] **Preserve reusable scripts currently in `publication_readiness/`**
  Move the following code files to a software-appropriate location, e.g. `AGVarPred/scripts/manuscript_corrections/` or `AGVarPred/AGVarPred-training/scripts/corrections/`:
  - `generate_figure2.py`
  - `generate_readiness_report.py`
  - `regenerate_test_predictions_with_variant_id.py`

- [ ] **Remove manuscript files from the working tree**
  - `rm -rf publication_readiness/`
  - `rm dca_analysis/manuscript_discussion.txt`
  - `rm dca_analysis/manuscript_results.txt`
  - `rm delong_analysis/manuscript_results.txt`
  - `rm tissue_shap_analysis/manuscript_discussion.txt`
  - `rm tissue_shap_analysis/manuscript_results.txt`

- [ ] **Update `.gitignore` to prevent accidental re-commit**
  Add to the top-level `.gitignore`:
  ```gitignore
  # Manuscript drafting materials (not part of software releases)
  /publication_readiness/
  */manuscript_results.txt
  */manuscript_discussion.txt
  ```

- [ ] **Verify no manuscript files remain in the git index**
  ```bash
  git rm -r --cached publication_readiness/ dca_analysis/manuscript_*.txt delong_analysis/manuscript_results.txt tissue_shap_analysis/manuscript_*.txt
  git status --short
  ```
  Expected: deletions staged; no new manuscript files tracked.

- [ ] **Decide on `AGVarPred-zenodo/` content**
  - Review `AGVarPred/AGVarPred-zenodo/` to confirm it contains only software/model/benchmark artifacts and no manuscript markdowns.
  - If it contains manuscript materials, remove or relocate them.

---

## Phase 2 — Version alignment

- [ ] **Bump version in `AGVarPred/pyproject.toml`**
  ```toml
  version = "1.0.8"
  ```

- [ ] **Bump version in `AGVarPred/src/AGVarPred/__init__.py`**
  ```python
  __version__ = "1.0.8"
  ```

- [ ] **Bump version in `AGVarPred/CITATION.cff`**
  ```yaml
  version: 1.0.8
  ```

- [ ] **Update `AGVarPred/CHANGELOG.md`**
  Add a `v1.0.8` section summarizing:
  - Repository cleaned of manuscript drafting materials
  - Packaging restored with wheel + sdist assets
  - Version strings aligned

- [ ] **Update `AGVarPred/RELEASE_REPORT.md`** (optional)
  - Update Zenodo DOI reference or remove version-specific DOI claims.

---

## Phase 3 — Build and local verification

- [ ] **Create a clean build environment**
  ```bash
  cd /tmp
  python -m venv build_env
  source build_env/bin/activate
  pip install build
  ```

- [ ] **Build wheel and sdist**
  ```bash
  cd /data/kanaka/pankaj/feature_extration/new/AGVarPred
  python -m build --outdir /tmp/agvarpred_dist
  ```
  Expected output:
  ```
  /tmp/agvarpred_dist/agvarpred-1.0.8-py3-none-any.whl
  /tmp/agvarpred_dist/agvarpred-1.0.8.tar.gz
  ```

- [ ] **Inspect wheel contents**
  ```bash
  unzip -l /tmp/agvarpred_dist/agvarpred-1.0.8-py3-none-any.whl
  ```
  Confirm presence of:
  - `AGVarPred/` package files
  - `agvarpred_core/` package files
  - `AGVarPred/model/` bundled models
  - `agvarpred-1.0.8.dist-info/` metadata + `LICENSE`

  Confirm absence of:
  - `publication_readiness/`
  - Any `manuscript_*.txt` files
  - Audit/forensic/review markdowns

- [ ] **Inspect sdist contents**
  ```bash
  tar -tzf /tmp/agvarpred_dist/agvarpred-1.0.8.tar.gz
  ```
  Confirm presence of:
  - `src/AGVarPred/`
  - `src/agvarpred_core/`
  - `tests/`
  - `README.md`, `LICENSE`, `pyproject.toml`

  Confirm absence of manuscript materials.

- [ ] **Install wheel in a fresh environment and run tests**
  ```bash
  python -m venv /tmp/wheel_test_env
  source /tmp/wheel_test_env/bin/activate
  pip install /tmp/agvarpred_dist/agvarpred-1.0.8-py3-none-any.whl
  pip install pytest pytest-cov
  python -c "import AGVarPred; print(AGVarPred.__version__)"  # should print 1.0.8
  AGVarPred --version
  pytest --pyargs AGVarPred  # or pytest /path/to/tests
  ```

- [ ] **Install sdist in a fresh environment and run tests**
  ```bash
  python -m venv /tmp/sdist_test_env
  source /tmp/sdist_test_env/bin/activate
  pip install /tmp/agvarpred_dist/agvarpred-1.0.8.tar.gz
  pip install pytest pytest-cov
  python -c "import AGVarPred; print(AGVarPred.__version__)"  # should print 1.0.8
  AGVarPred --version
  pytest --pyargs AGVarPred
  ```

- [ ] **Run `twine check`** (optional but recommended)
  ```bash
  pip install twine
  twine check /tmp/agvarpred_dist/agvarpred-1.0.8-py3-none-any.whl /tmp/agvarpred_dist/agvarpred-1.0.8.tar.gz
  ```

---

## Phase 4 — Git tag and push

- [ ] **Commit all cleanup and version-bump changes**
  ```bash
  git add -A
  git commit -m "chore: clean manuscript materials and bump version to 1.0.8

  - Remove publication_readiness/ and manuscript text files from tracked tree.
  - Add .gitignore patterns to prevent re-commit.
  - Align version strings in pyproject.toml, __init__.py, and CITATION.cff.
  - Update CHANGELOG."
  ```

- [ ] **Tag v1.0.8**
  ```bash
  git tag -a v1.0.8 -m "v1.0.8 - software-only release with clean packaging"
  ```

- [ ] **Push commit and tag**
  ```bash
  git push origin main
  git push origin v1.0.8
  ```

---

## Phase 5 — GitHub release

- [ ] **Create the GitHub release**
  - Use tag `v1.0.8`.
  - Title: `v1.0.8 - Clean software release`.
  - Body: summarize that this release removes manuscript drafting materials and restores wheel/sdist assets.

- [ ] **Upload build artifacts**
  Attach both:
  - `agvarpred-1.0.8-py3-none-any.whl`
  - `agvarpred-1.0.8.tar.gz`

- [ ] **Verify release page**
  Ensure the release shows four assets:
  - `agvarpred-1.0.8-py3-none-any.whl`
  - `agvarpred-1.0.8.tar.gz`
  - `Source code.zip`
  - `Source code.tar.gz`

---

## Phase 6 — Zenodo strategy

Choose **Option A** or **Option B**:

### Option A — Manual Zenodo upload (recommended for full control)

- [ ] Disable GitHub→Zenodo auto-sync at https://zenodo.org/account/settings/github/ (toggle off `pankaj357/AGVarPred`).
- [ ] Delete or retract the existing manuscript-containing Zenodo records (requires owner action).
- [ ] Create a new Zenodo upload manually.
- [ ] Upload only:
  - `agvarpred-1.0.8-py3-none-any.whl`
  - `agvarpred-1.0.8.tar.gz`
- [ ] Add metadata (title, authors, license Apache-2.0, keywords).
- [ ] Publish and record the new DOI.
- [ ] Update `AGVarPred/README.md` and `CITATION.cff` with the new Zenodo DOI.

### Option B — Clean repo + auto-sync

- [ ] Keep GitHub→Zenodo auto-sync enabled.
- [ ] After Phase 5, wait for Zenodo to create a new version from the `v1.0.8` release.
- [ ] Download the Zenodo archive and verify it contains no `publication_readiness/` or manuscript text files.
- [ ] If Zenodo allows, upload the wheel + sdist as additional files to the same record.

---

## Phase 7 — Post-release verification

- [ ] `pip install https://github.com/pankaj357/AGVarPred/releases/download/v1.0.8/agvarpred-1.0.8-py3-none-any.whl` works.
- [ ] `pip install https://github.com/pankaj357/AGVarPred/releases/download/v1.0.8/agvarpred-1.0.8.tar.gz` works.
- [ ] `AGVarPred --version` prints `1.0.8`.
- [ ] Zenodo record resolves and contains no manuscript markdowns.
- [ ] README badge/link points to the correct release and DOI.

---

## Optional next steps

- [ ] Set up a GitHub Actions workflow to build and upload wheel + sdist automatically on tagged releases.
- [ ] Publish `v1.0.8` to PyPI (`twine upload dist/*`).
- [ ] Create a `MANIFEST.in` for explicit sdist control.
