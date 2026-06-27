# AGVarPred v1.0.0 Release Readiness Report

Generated: 2026-06-27

## Summary

The AGVarPred repository has been audited, metadata placeholders replaced, and
publication artifacts prepared for the v1.0.0 release. All validation checks
pass. The only remaining item requiring external action is the Zenodo DOI.

## Files modified

### Metadata and packaging
- `pyproject.toml` — real authors/maintainers, emails, GitHub URLs, SPDX license.
- `CITATION.cff` — real authors, ORCIDs, affiliations, GitHub URL.
- `AGVarPred-training/CITATION.cff` — same metadata.
- `AGVarPred-zenodo/CITATION.cff` — same metadata.
- `AGVarPred-zenodo/code_reference.txt` — real GitHub release URL.
- `AGVarPred-zenodo/data_manifest.json` — added gnomAD file sizes.
- `AGVarPred-zenodo/checksums.sha256` — regenerated after edits.

### Documentation
- `README.md` — real GitHub URLs, contact emails, acknowledgements, references.
- `MODEL_CARD.md` — real GitHub URL, acknowledgements, references.
- `docs/installation.md` — real GitHub URL, AlphaGenome note.
- `docs/usage.md` — added `--version` example.
- `docs/api.md` — clarified `AGVarPredPredictor` constructor.
- `docs/reproducibility.md` — ClinVar download URL, Makefile note.
- `CONTRIBUTING.md` — added security disclosure note.
- `CODE_OF_CONDUCT.md` — real contact emails.
- `SECURITY.md` — real contact emails, supported-versions table.
- `RELEASE_CHECKLIST.md` — updated with completed steps.

### Training reproducibility
- `AGVarPred-training/Makefile` — improved `evaluate` and `benchmark` targets.
- `AGVarPred-training/scripts/external_validation/benchmark_cadd_revel_alphamissense.py` — replaced hardcoded path with `PROJECT_ROOT`, HTTPS for myvariant.info.
- `AGVarPred-training/scripts/external_validation/fetch_vep_rest_api.py` — replaced hardcoded path with `PROJECT_ROOT`.
- `AGVarPred-training/scripts/external_validation/fetch_clinvar_vep_rest_api.py` — replaced hardcoded path with `PROJECT_ROOT`.
- `AGVarPred-training/scripts/external_validation/create_benchmark_comparison_figure.py` — replaced hardcoded path with `PROJECT_ROOT`.
- `AGVarPred-training/scripts/external_validation/build_clinvar_benchmark.py` — relative default VCF path.
- `AGVarPred-training/scripts/external_validation/run_grimm2015_fe_parallel.sh` — use existing `PROJECT_ROOT` env var.

### Source code
- `src/AGVarPred/predictor.py` — corrected docstring for `model_name` parameter.

### Model manifests
- `model/model_full/manifest.yaml` — filled `git_commit`.
- `model/model_no_af/manifest.yaml` — filled `git_commit`.
- `AGVarPred-zenodo/model/model_full/manifest.yaml` — filled `git_commit`.
- `AGVarPred-zenodo/model/model_no_af/manifest.yaml` — filled `git_commit`.

### Repository hygiene
- Removed committed generated artifacts (`dist/`, `src/AGVarPred.egg-info/`, `tests/__pycache__/`, `.pytest_cache/`).
- Initialized git repository, made release commit, and tagged `v1.0.0`.
- Renamed default branch to `main` to match GitHub Actions workflows.

## Git state

- Branch: `main`
- Tag: `v1.0.0` → `afc474aa1d91067aef03650e9fff85a468936542`
- Model release commit (in manifest `git_commit`): `52a8ede5e1c7b437d7bdf760016325a35f9848f7`
- Working tree: clean

## Remaining placeholders requiring owner input

1. **Zenodo DOI**
   - Currently `10.5281/zenodo.TODO` in:
     - `CITATION.cff`
     - `AGVarPred-training/CITATION.cff`
     - `AGVarPred-zenodo/CITATION.cff`
     - `README.md` (badge and citation)
     - `MODEL_CARD.md`
   - Action required: Create Zenodo record, obtain DOI, replace all occurrences.

2. **AlphaGenome SDK public instructions**
   - `docs/installation.md` currently notes that SDK instructions will be added
     when publicly available.
   - Action required: Add official AlphaGenome SDK installation URL/command.

## Inconsistencies found and resolved

| Issue | Resolution |
|-------|------------|
| Placeholder author names/ORCID in CITATION.cff files | Replaced with real author metadata. |
| Placeholder GitHub URLs (`<GITHUB_ORG>`) | Replaced with `pankaj357`. |
| Placeholder contact emails | Replaced with `ft.pank@gmail.com` and `kkokay07@gmail.com`. |
| Placeholder Zenodo DOI | Intentionally retained; DOI must be obtained from Zenodo. |
| Hardcoded local paths in training scripts | Replaced with `PROJECT_ROOT` environment variable. |
| `http://myvariant.info` (insecure) | Changed to `https://myvariant.info`. |
| Missing gnomAD file sizes in data manifest | Added real sizes (~85.3 GiB VCF, ~913 KiB index). |
| Generated artifacts committed | Removed; `.gitignore` already covers them. |
| Incomplete Makefile `evaluate`/`benchmark` targets | Improved with useful output and best-effort benchmark runs. |
| `project.license` deprecation warning | Updated to SPDX string `license = "Apache-2.0"` and bumped setuptools. |
| git_commit placeholders in model manifests | Filled with release commit SHA. |

## Validation results

| Check | Result |
|-------|--------|
| `pytest tests/` | 21 passed, 1 skipped (gnomAD VCF not available) |
| `python -m build` | Success (wheel + sdist built) |
| CLI on `examples/sample.vcf` (no gnomAD) | Success, 5 predictions |
| `sha256sum -c AGVarPred-zenodo/checksums.sha256` | All OK |
| Git working tree | Clean |

## Recommendations before creating the GitHub repository

1. **Create the public GitHub repository** at `https://github.com/pankaj357/AGVarPred`.
2. **Push the local repository**:
   ```bash
   cd /data/kanaka/pankaj/feature_extration/new/AGVarPred
   git remote add origin https://github.com/pankaj357/AGVarPred.git
   git push -u origin main
   git push origin v1.0.0
   ```
3. **Create a GitHub Release** from tag `v1.0.0` and attach the built wheel/sdist
   from `dist/`.
4. **Enable GitHub Actions** so the test and release workflows run.
5. **Link the repository to Zenodo**, create the v1.0.0 Zenodo record, and obtain
   the DOI.
6. **Replace `10.5281/zenodo.TODO`** everywhere with the real DOI, commit, and
   (optionally) move the `v1.0.0` tag to include that final DOI update.
7. **Regenerate `AGVarPred-zenodo/checksums.sha256`** after the DOI replacement.
8. **Add repository topics** on GitHub (e.g., `variant-pathogenicity`,
   `machine-learning`, `genomics`, `AlphaGenome`).
9. **Add the official AlphaGenome SDK installation instructions** to
   `docs/installation.md` when they become available.

## Conclusion

The repository is **ready for a public v1.0.0 release** once the Zenodo DOI is
assigned and the final DOI replacement step is completed. No scientific
methodology, machine learning models, feature engineering, or prediction logic
were changed.
