# Multi-Cohort Peripheral Blood Transcriptomics Identifies a Six-Gene ALS Diagnostic Signature with Network-Aware Biomarker Prioritization

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
|---|---|---|---|---|
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
|---|---|
| `00_Install_Dependencies.R` | Installs all required CRAN and Bioconductor packages |
| `01_Data_Prep_and_DEG.R` | Data preprocessing, quality control, and differential expression analysis |
| `02_WGCNA_Analysis.R` | Co-expression network construction, module-trait correlation, and hub gene extraction |
| `03_Machine_Learning_Consensus.R` | Triple-consensus feature selection and logistic regression model training |
| `04_Clinical_Validation.R` | Five-arm validation, BPS computation, and ALS-mimic discrimination analysis |
| `05_Immune_Analysis.R` | EPIC, MCP-counter, and ssGSEA immune deconvolution and subtype stratification |
| `06_Drug_Repurposing.R` | DSigDB drug repurposing via Enrichr |

### Note on Functional Enrichment and PPI Analysis

Functional enrichment (GO Biological Process, GO Cellular Component, Reactome) and 
PPI network construction were performed using the **STRING database web server (v12.0)** 
(https://string-db.org/). The 240 ML hub candidates were used as enrichment input and 
the six consensus genes with first-degree interactors were used for the PPI network. 
These steps are not included as R scripts as they use the STRING graphical interface.

---

## Reproducibility Notes

- **Deterministic Seed Control:** A global seed (`set.seed(1122)`) is strictly enforced across all scripts. This is critical to guarantee that all stochastic processes, such as Random Forest tree generation, LASSO cross-validation folds, SVM recursive feature elimination, and WGCNA module clustering - yield identical, perfectly reproducible results on every run.
- **Sequential Execution:** The pipeline relies on a continuous computational environment. Downstream scripts (`04` through `07`) explicitly load the filtered matrices, trained models, and signature lists saved to the `Processed_Data/` directory by scripts `01` through `03`. Please ensure the scripts are run in exact numerical order.
- **Automated Directories:** The scripts are designed to automatically create `Processed_Data/` and `Manuscript_Figures/` directories in your root working folder to securely store intermediate data objects (.rds, .csv) and high-resolution plots (.pdf, .png) without cluttering your environment.
- **GEO Data Retrieval (Timeouts):** Scripts `01` and `04` use `GEOquery` to download large transcriptomic matrices directly from the NCBI servers. If you experience "timeout" errors due to internet connectivity, simply increase your R timeout limit before running the script by executing: `options(timeout = 600)`.
- **Hardware Recommendations:** While the pipeline is highly streamlined, WGCNA topology calculation and triple-consensus machine learning cross-validation are memory-intensive. A standard modern workstation with at least 8GB to 16GB of RAM is recommended for seamless execution.

---

## Citation

If you use this pipeline in your research, please cite:

> [Citation to be added upon publication]

---

## Contact

- **Academic & Collaborative Inquiries:** For questions regarding the manuscript, methodology, or potential collaborations, please contact the corresponding authors via email as listed in the publication.
- **Code & Pipeline issues:**For technical questions, troubleshooting, or code-related bugs, please [open an Issue](https://github.com/RabiaSultan8/ALS_Project/issues) directly in this GitHub repository.
