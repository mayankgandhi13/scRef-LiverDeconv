# scRef-LiverDeconv

## Overview
This project applies RNA deconvolution to bulk RNA-seq data from human liver tissue to characterize how cell type composition varies with age and sex in healthy individuals.

We replicate the methodology of Conning-Rowland et al. (2025) — originally applied to heart and skeletal muscle — and extend it to human liver using two deconvolution tools: **CIBERSORTx** and **BayesPrism**.

## Biological Question
Does the cellular composition of healthy human liver differ by age or sex?

## Data Sources
| Data | Source | Description |
|---|---|---|
| Single-cell reference | MacParland et al. 2018, GEO: GSE115469 | 8,444 liver cells, 20 annotated populations |
| Bulk RNA-seq | GTEx v8 liver | ~226 healthy human donors |

## Cell Lineages
| Lineage | Source cell types | Cells in reference |
|---|---|---|
| Hepatocyte | Hepatocyte_1 to _6 | 200 |
| Kupffer | Inflammatory + Non-inflammatory Macrophage | 200 |
| LSEC | Central venous + Periportal LSECs | 200 |
| Lymphoid | T cells, NK cells, Plasma cells | 200 |

## Pipeline
- **Phase 1** — Build single-cell reference from MacParland atlas
- **Phase 2** — In silico validation using 100 synthetic tissues
- **Phase 3** — Apply to GTEx liver and test age/sex associations

## Tools
- R 4.4.2
- Seurat 5.4.0
- CIBERSORTx (web portal)
- BayesPrism 2.2.2

## Template
Conning-Rowland et al. (2025) *Heliyon* — Application of CIBERSORTx and BayesPrism to deconvolution of bulk RNA-seq data from human myocardium and skeletal muscle.
