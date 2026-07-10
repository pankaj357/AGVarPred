# AGVarPred Release Audit

**Audit date:** 2026-07-10  
**Scope:** GitHub releases v1.0.4–v1.0.7, Zenodo archiving, and asset hygiene.

---

## 1. Release History

| Release | Tag date | GitHub assets | Wheel present? | SDist present? | Source archive contains manuscript? | Zenodo archived? |
|---|---|---|---|---|---|---|
| **v1.0.4** | earlier | `agvarpred-1.0.4-py3-none-any.whl`, `agvarpred-1.0.4.tar.gz`, `Source code.zip`, `Source code.tar.gz` | ✅ Yes | ✅ Yes | ❓ Likely yes (repository contained `publication_readiness/` at that time) | ✅ Yes (DOI `10.5281/zenodo.20955782`) |
| **v1.0.5** | 2026-07-10 | Only `Source code.zip`, `Source code.tar.gz` | ❌ No | ❌ No | ✅ Yes | ✅ Yes (DOI `10.5281/zenodo.21300056`) |
| **v1.0.6** | 2026-07-10 | Only `Source code.zip`, `Source code.tar.gz` | ❌ No | ❌ No | ✅ Yes | ⏳ Pending sync |
| **v1.0.7** | 2026-07-10 | Only `Source code.zip`, `Source code.tar.gz` | ❌ No | ❌ No | ✅ Yes | ⏳ Pending sync |

> Note: `v1.0.5`, `v1.0.6`, and `v1.0.7` were created with the GitHub Releases API but no wheel/sdist files were uploaded.

---

## 2. Verified Contents of GitHub Source Archive (v1.0.7)

Downloaded `https://github.com/pankaj357/AGVarPred/archive/refs/tags/v1.0.7.tar.gz` and inspected:

**Manuscript materials present:**

- `publication_readiness/abstract.md`
- `publication_readiness/abstract_corrected.md`
- `publication_readiness/authoritative_benchmark_summary.csv`
- `publication_readiness/change_log.md`
- `publication_readiness/discussion_section.md`
- `publication_readiness/final_consistency_report.md`
- `publication_readiness/final_publication_audit_report.md`
- `publication_readiness/final_submission_audit.md`
- `publication_readiness/forensic_audit_report.md`
- `publication_readiness/manuscript_assembly.md`
- `publication_readiness/manuscript_corrected_combined.md`
- `publication_readiness/methods_section.md`
- `publication_readiness/publication_readiness_report.md`
- `publication_readiness/results_section.md`
- `publication_readiness/supplementary_materials_combined.md`
- `publication_readiness/vip_consistency_analysis.md`
- `dca_analysis/manuscript_discussion.txt`
- `dca_analysis/manuscript_results.txt`
- `delong_analysis/manuscript_results.txt`
- `tissue_shap_analysis/manuscript_discussion.txt`
- `tissue_shap_analysis/manuscript_results.txt`

**Also present (non-package, reproducibility payload):**

- `AGVarPred-zenodo/` — external benchmark CSVs, model `.pkl`, splits
- `AGVarPred/AGVarPred-training/` — full training pipeline scripts
- `dca_analysis/`, `delong_analysis/`, `tissue_shap_analysis/`, `variant_type_analysis/` — analysis scripts and small CSVs

These are software/reproducibility artifacts and are acceptable in a software repository, but the manuscript text files should not be part of a software release.

---

## 3. Issues Identified

### Issue 1 — Manuscript material in source archives
**Severity:** High  
**Impact:** Zenodo archives and GitHub source archives contain manuscript drafts, audit reports, and review reports. This contradicts the goal of keeping GitHub as a software repository.

### Issue 2 — Missing installable assets on v1.0.5–v1.0.7
**Severity:** High  
**Impact:** Users cannot `pip install` from the GitHub release page. Only the source archive is available, which is not a PyPI-style distribution.

### Issue 3 — Version string mismatch
**Severity:** Medium  
**Impact:**

- GitHub tag: `v1.0.7`
- `pyproject.toml`: `1.0.4`
- `src/AGVarPred/__init__.py`: `1.0.4`
- `CITATION.cff`: `1.0.3`
- `AGVarPred-zenodo/model/*/manifest.yaml`: reports model versions but package version is not aligned

This is confusing and breaks `pip show AGVarPred` expectations.

### Issue 4 — Zenodo mirrors GitHub source archive exactly
**Severity:** Medium  
**Impact:** Even if wheel/sdist are uploaded to GitHub, Zenodo will still archive the full source snapshot. The repository must be cleaned **before** tagging, or Zenodo integration must be disabled in favor of manual uploads.

### Issue 5 — `RELEASE_REPORT.md` references stale DOI
**Severity:** Low  
**Impact:** `AGVarPred/RELEASE_REPORT.md` points to `10.5281/zenodo.20955782` (v1.0.4 concept DOI). It should be updated or removed from release-specific claims.

---

## 4. Recommendations

1. **Delete or hide manuscript materials from the main branch** before creating `v1.0.8`. Options:
   - Move them to a private repository.
   - Move them to a local directory outside git.
   - Keep them on a `manuscript` branch that is never tagged for release.
2. **Build wheel + sdist with `python -m build`** and upload both to the GitHub release.
3. **Align all version strings** to `1.0.8`.
4. **For Zenodo:**
   - **Preferred:** Disable GitHub auto-sync and manually create a Zenodo upload containing only the wheel and sdist.
   - **Alternative:** Clean the repo, tag `v1.0.8`, let Zenodo auto-archive the source snapshot, and manually upload the wheel/sdist as additional files.
5. **Do not delete existing GitHub releases** unless necessary; instead, publish `v1.0.8` as the clean release and update the README to point users to it.
6. **Consider deleting the Zenodo records for v1.0.5–v1.0.7** if they remain unpublished/pending, so they do not become the default resolved versions.
