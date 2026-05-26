# Peripheral Blood Transcriptomics Identifies a Six-Gene Diagnostics Signature with Network-Aware Biomarker Prioritization and Uncovers Transcriptomic Overlap with Motor Neuron Disease Mimics

This repository contains the complete R pipeline for the manuscript:

**"Peripheral Blood Transcriptomics Identifies a Six-Gene Diagnostics
Signature with Network-Aware Biomarker Prioritization and Uncovers
Transcriptomic Overlap with Motor Neuron Disease Mimics"**

---

## Overview

This project presents a systems-integrated, leakage-free bioinformatics framework for
peripheral blood biomarker discovery in ALS. Using GSE112676 (n = 702) as the sole
discovery cohort, the pipeline applies differential expression analysis, WGCNA, and a
triple-consensus machine learning framework (Random Forest, LASSO, SVM-RFE) to derive
a six-gene diagnostic signature — **BRI3, ABCA1, QPCT, PPP2R5A, ETFRF1, and SLC37A3**.
A novel Biomarker Priority Score (BPS) integrating cross-cohort diagnostic stability
with network centrality is introduced for principled biomarker ranking. The signature is
validated across five independent arms spanning three external cohorts, including a
dedicated ALS-versus-mimic discrimination analysis.

---

## Data Availability

All datasets are publicly available through the NCBI Gene Expression Omnibus (GEO):

| Dataset | Platform | Role | n | Groups |
|---------|----------|------|---|--------|
| GSE112676 | Illumina HumanHT-12 V3.0 | Discovery | 702 | ALS, CON |
| GSE112680 | Illumina HumanHT-12 V4.0 | Validation (held-out) | 376 | ALS, CON, MIM |
| GSE28253 | Agilent 4×44K GPL4133 | Validation (external) | 22 | ALS, CON |
| GSE234297 | Illumina RNA-seq | Validation (cross-platform) | 144 | ALS, CON |

**Note:** GSE112680 was held out entirely and never accessed during discovery, feature
selection, or model training.

---

## Prerequisites

The analysis was performed in **R version 4.5.2** with a global seed (`set.seed(1122)`).
To install all required packages, run:

```r
source("00_Install_Dependencies.R")
```

---

## Pipeline Structure

Run scripts in sequential order. Each script loads intermediate objects saved by the
preceding script from the `Processed_Data/` directory.

| Script | Description |
|--------|-------------|
| `00_Install_Dependencies.R` | Installs all required CRAN and Bioconductor packages |
| `01_Data_Prep_and_DEG.R` | Data preprocessing, quality control, and differential expression analysis |
| `02_WGCNA_Analysis.R` | Co-expression network construction, module-trait correlation, and hub gene extraction |
| `03_Machine_Learning_Consensus.R` | Triple-consensus feature selection and logistic regression model training; discovery-cohort gene-wise scaling parameters (mean and SD) saved alongside the trained model for deployment-ready normalization |
| `04_Clinical_Validation_and_BPS.R` | Five-arm independent validation using discovery-locked z-scoring, BPS computation, and ALS-mimic discrimination analysis |
| `05_Immune_Analysis.R` | EPIC, MCP-counter, and ssGSEA immune deconvolution and ALS immune subtype stratification |
| `06_Drug_Repurposing.R` | DSigDB drug repurposing via Enrichr |
| `07_Cohort_Table_and_Covariate_Analysis.R` | Cohort characteristics table summarising sex, age at onset, site of onset, and survival across all groups in GSE112676 and GSE112680; sex-adjusted limma sensitivity analysis |
| `08_Threshold_CI_Analysis.R` | Youden-optimal decision threshold with 2000-bootstrap 95% confidence intervals for sensitivity and specificity across GSE112680 Arms 1–3 |
| `09_Immune_Subtype_Validation.R` | Clinical anchoring of ALS immune subtypes against survival (Kaplan-Meier, log-rank), age at onset (Mann-Whitney U), site of onset, and sex (Fisher's exact); external validation by EPIC deconvolution and centroid projection in GSE112680 (n = 164 ALS) |

### Note on Functional Enrichment and PPI Analysis

Functional enrichment (GO Biological Process, GO Cellular Component, Reactome) and
PPI network construction were performed using the **STRING database web server (v12.0)**
(https://string-db.org/). The 240 ML hub candidates were used as enrichment input and
the six consensus genes with first-degree interactors were used for the PPI network.
These steps are not included as R scripts as they use the STRING graphical interface.

---

## Reproducibility Notes

- **Deterministic Seed Control:** A global seed (`set.seed(1122)`) is strictly enforced
  across all scripts. This guarantees that all stochastic processes — Random Forest tree
  generation, LASSO cross-validation folds, SVM recursive feature elimination, and WGCNA
  module clustering — yield identical, perfectly reproducible results on every run.

- **Sequential Execution:** The pipeline relies on a continuous computational environment.
  Scripts `04` through `09` explicitly load filtered matrices, trained models, and
  signature lists saved to the `Processed_Data/` directory by the preceding scripts.
  Please run scripts in exact numerical order.

- **Discovery-Locked Normalization:** Gene-wise means and
  standard deviations are computed once from the GSE112676 discovery cohort (n = 702) and
  stored in `Step4_ML_Results.rds`. All validation cohorts are scaled using these fixed
  parameters (`z = (x − μ_discovery) / σ_discovery`) rather than their own cohort
  statistics, making the model applicable to a single future patient without requiring
  cohort-level data.

- **Automated Directories:** The scripts automatically create `Processed_Data/` and
  `Manuscript_Figures/` directories in the working folder to store intermediate data
  objects (`.rds`, `.csv`) and high-resolution figures (`.pdf`, `.png`).

- **GEO Data Retrieval (Timeouts):** Scripts `01`, `04`, `07`, and `09` download
  transcriptomic matrices and phenotype records directly from NCBI GEO via `GEOquery`.
  If you encounter timeout errors, increase the R timeout limit before running:
  `options(timeout = 600)`. Local series matrix caches (e.g.,
  `GSE112676_series_matrix.txt.gz`) are detected automatically and loaded in preference
  to live downloads where available.

- **Hardware Recommendations:** WGCNA soft-threshold topology fitting and triple-consensus
  machine learning cross-validation are memory-intensive. A workstation with at least
  8–16 GB of RAM is recommended for seamless execution.

---

## Citation

If you use this pipeline in your research, please cite:

> [Citation to be added upon publication]

---

## Contact

- **Academic & Collaborative Inquiries:** For questions regarding the manuscript,
  methodology, or potential collaborations, please contact the corresponding authors
  via email as listed in the publication.
- **Code & Pipeline Issues:** For technical questions, troubleshooting, or code-related
  bugs, please [open an Issue](https://github.com/RabiaSultan8/ALS_Project/issues)
  directly in this GitHub repository.
