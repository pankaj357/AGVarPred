# VIP Database — ClinVar-Label Consistency Analysis

## Status

The VIP Database is **not included in the main independent external-benchmark comparison** because all 1,075 of its pathogenic labels are derived from overlapping ClinVar annotations. It is reported here as a separate **ClinVar-held-out-gene consistency analysis**.

## Why VIP is not an independent benchmark

| Property | VIP | Independent benchmarks |
|---|---|---|
| Variant/gene holdout from train/cal/test | ✅ Yes | ✅ Yes |
| Label source | ClinVar-derived (`review` + ClinVar Path/LP) | Humsavar, MAVE assays, Grimm2015, DVD, gnomAD AF |
| Tests generalization to new genes | ✅ Yes | ✅ Yes |
| Tests generalization to independent curation | ❌ No | ✅ Yes |

Because AGVarPred is trained on ClinVar labels, using ClinVar-derived labels for external validation creates a label-source circularity risk. The VIP analysis is therefore interpreted as consistency of ClinVar classifications in genes that were not seen during training, not as independent proof of generalization.

## Benchmark composition

| Metric | Value |
|---|---|
| Total variants | 3,880 |
| Pathogenic / Likely pathogenic | 1,075 |
| Benign / Likely benign | 2,805 |
| Genes | 904 |

Breakdown of ClinVar-derived pathogenic labels:

| `label_source` | Count |
|---|---|
| `VIP_review+ClinVar_Pathogenic` | 667 |
| `VIP_review+ClinVar_Likely pathogenic` | 348 |
| `VIP_review+ClinVar_Pathogenic/Likely pathogenic` | 53 |
| Other composite ClinVar categories | 7 |
| **Total ClinVar-derived pathogenic** | **1,075** |

## AGVarPred performance on VIP

| Metric | Value |
|---|---|
| ROC-AUC (95% CI) | 0.892 (0.881–0.903) |
| PR-AUC | 0.753 |
| MCC (95% CI) | 0.507 (0.484–0.529) |
| Full coverage | 3,880 / 3,880 (100%) |

## Competitor performance on common subsets

| Tool | N | AUC |
|---|---|---|
| AGVarPred | 3,862 | 0.892 |
| CADD | 3,862 | 0.927 |
| REVEL | 1,715 | 0.876 |
| AlphaMissense | 153 | 0.803 |

## Interpretation

AGVarPred reproduces ClinVar-derived pathogenic classifications in held-out genes with good discrimination (AUC ~0.89). The benchmark also illustrates the coverage advantage of AGVarPred: it scored all 3,880 VIP variants, whereas AlphaMissense scored only 153 (3.9%) and REVEL scored 1,715 (44.2%). However, because the pathogenic labels originate from ClinVar, this result cannot be claimed as generalization beyond the ClinVar curation source.

## Manuscript placement

VIP is discussed in Results §2.1 ("ClinVar-label consistency in held-out genes") and is excluded from Figures 3–6 and from the independent external-benchmark count.
