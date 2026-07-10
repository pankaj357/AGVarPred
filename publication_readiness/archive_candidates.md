# Archive Candidates — AGVarPred Repository Synchronization

**Generated:** 2026-07-10 during final manuscript synchronization pass.

This list identifies files that are obsolete, duplicated, or superseded after the corrected analyses were regenerated. **No files have been deleted automatically.** The repository owner should review this list and decide whether to archive or remove each item.

| File | Reason | Safe to archive? |
|---|---|---|
| `publication_readiness/abstract_corrected.md` | Exact duplicate of `publication_readiness/abstract.md` after synchronization. Originally intended as a corrected draft but now redundant. | Yes |
| `publication_readiness/introduction_corrected.md` | Exact duplicate of `publication_readiness/introduction.md`. No longer needed once `manuscript_corrected_combined.md` becomes the single source of truth. | Yes |
| `publication_readiness/final_publication_audit_report.md` | Superseded by `forensic_audit_report.md` and `final_consistency_report.md`, which reflect the post-correction state. | Yes (after confirming no unique content) |
| `publication_readiness/publication_readiness_report.md` | Older readiness report generated before the variant-type alignment fix and VIP removal. Superseded by `manuscript_corrected_combined.md` and the final consistency report. | Yes |
| `publication_readiness/change_log.md` | Historical log of claim corrections. Useful for transparency but no longer an active manuscript source. Consider archiving or converting to appendix. | Yes (historical reference) |
| `publication_readiness/table_audit.csv` | Appears to be an empty/stale audit template. Not referenced by any active script or manuscript section. | Review contents first; likely safe |
| `publication_readiness/supplementary_checklist.csv` | Contains placeholder supplementary entries. Should be replaced by the final supplementary materials or removed if unused. | Review before archiving |
| `publication_readiness/manuscript_assembly.md` | Outdated assembly instructions. Will be superseded by `manuscript_corrected_combined.md`. | Yes, once combined manuscript is finalized |
| `external_validation/scripts/` (git-tracked copy in `AGVarPred/` subtree) | The working `external_validation/scripts/` directory has diverged from the git-tracked version. The git-tracked copy may be stale. Resolve divergence before archiving. | No — requires manual reconciliation |
| `test (3.11)/` and `test (3.13)/` directories | Appear to contain single text files with environment setup notes. Directory names with spaces and parentheses are awkward for scripting. | Review; likely safe to archive or merge into `docs/` |

## Notes

- **Do not delete** any file without confirming it is not referenced by an active script, figure, or manuscript section.
- The authoritative outputs after synchronization are:
  - `publication_readiness/abstract.md`
  - `publication_readiness/introduction.md`
  - `publication_readiness/methods_section.md`
  - `publication_readiness/results_section.md`
  - `publication_readiness/discussion_section.md`
  - `publication_readiness/figure_plan.md`
  - `publication_readiness/manuscript_corrected_combined.md` (created in this pass)
  - `publication_readiness/authoritative_benchmark_summary.csv`
  - `variant_type_analysis/*.csv`
  - `delong_analysis/results/*.csv`
  - `dca_analysis/results/*.csv`
