# Peripheral Blood Transcriptomics Identifies XPO1-Led Six-Gene Diagnostic Signature for ALS

This repository contains the complete R scripts and computational pipeline for the manuscript: 
**"Peripheral Blood Transcriptomics Identifies XPO1-Led Six-Gene Diagnostic Signature for ALS via Consensus Machine Learning and Network-Based Biomarker Prioritization."**

## Overview
This project provides a systems-integrated bioinformatics framework to identify and rigorously prioritize peripheral blood diagnostic biomarkers for Amyotrophic Lateral Sclerosis (ALS). The pipeline harmonizes large-scale microarray data, performs co-expression network modeling, and utilizes a strict triple-consensus machine learning framework (Random Forest, LASSO, SVM-RFE) to isolate a robust 6-gene diagnostic signature led by *XPO1*.

## Data Availability
All raw datasets used in this pipeline are publicly available via the NCBI Gene Expression Omnibus (GEO):
- **Discovery Cohorts:** GSE112676, GSE112680 (Illumina HumanHT-12 V3.0)
- **Validation Cohort 1:** GSE28253 (Microarray)
- **Validation Cohort 2:** GSE234297 (RNA-seq)

## Prerequisites and Installation
The analysis was performed in **R (version 4.5.2)** under deterministic seed control (`set.seed(1122)`). 
To automatically install all required CRAN and Bioconductor packages, run the initial setup script:
- `00_Install_Dependencies.R` 

## Pipeline Structure
_(Note: Please run the scripts in the following sequential order to replicate the study's findings)_
- `01_Data_Prep_and_DEG.R` Downloads GEO datasets, maps probes to symbols, merges the discovery cohorts, applies `ComBat` batch correction, and performs `limma` differential expression analysis. (Generates Figure 2).
- `02_WGCNA_Analysis.R` Constructs the Weighted Gene Co-expression Network, identifies the ALS-associated turquoise module, and extracts high-connectivity hub genes. (Generates Figure 3).
- `03_Machine_Learning_Consensus.R` Runs Random Forest, LASSO, and SVM-RFE. Calculates the strict intersection of all three algorithms to extract the 24-gene consensus signature and plots the discovery ROC. (Generates Figure 5).
- `04_Clinical_Validation_and_BPS.R` Constructs the logistic regression model, blindly validates the 24-gene signature on the independent GSE28253 and GSE234297 cohorts, and calculates the novel Biomarker Priority Score (BPS) to derive the final elite 6-gene panel. (Generates Figure 6A, 6C).
- `05_Single_Gene_ROC_Validation.R` A script that isolates the top 6 BPS biomarkers and calculates their individual diagnostic performance (AUC) across all discovery and validation cohorts. (Generates Figues 6B).
- `06_Immune_Infiltration_ssGSEA.R` Executes ssGSEA using a curated marker dictionary to estimate peripheral immune cell fractions and correlates them with the elite gene signature. (Generates Figure 6D, 6E).
- `07_Drug_Repurposing.R` Queries the consensus signature against the DSigDB database via `enrichR` to identify highly significant therapeutic repurposing candidates. (Generates Figure 6F).

### Note on Functional Enrichment and PPI Analysis
Functional enrichment and Protein-Protein Interaction (PPI) network construction were performed using the **STRING database web server (v12.0)** (https://string-db.org/). 
- **Pathway Enrichment (GO/Cell Process/Reactome):** To characterize the broad biological processes of the ALS-associated network, the highly significant core genes of the WGCNA turquoise module were used as the input (`Turquoise_Core_Enrichment_Genes.csv`). 
- **PPI Network:** To map the exact topological interactions of the diagnostic features, the 24-gene machine learning consensus signature was used as the input (`Final_Gene_Signature.csv`).

*(Because these steps utilize the STRING graphical web interface, they are not included as standalone R scripts in this repository).*

## Notes on Reproducibility
- **Deterministic Seed Control:** A global seed (`set.seed(1122)`) is strictly enforced across all scripts. This is critical to guarantee that all stochastic processes, such as Random Forest tree generation, LASSO cross-validation folds, SVM recursive feature elimination, and WGCNA module clustering - yield identical, perfectly reproducible results on every run.
- **Sequential Execution:** The pipeline relies on a continuous computational environment. Downstream scripts (`04` through `07`) explicitly load the filtered matrices, trained models, and signature lists saved to the `Processed_Data/` directory by scripts `01` through `03`. Please ensure the scripts are run in exact numerical order.
- **Automated Directories:** The scripts are designed to automatically create `Processed_Data/` and `Manuscript_Figures/` directories in your root working folder to securely store intermediate data objects (.rds, .csv) and high-resolution plots (.pdf, .png) without cluttering your environment.
- **GEO Data Retrieval (Timeouts):** Scripts `01` and `04` use `GEOquery` to download large transcriptomic matrices directly from the NCBI servers. If you experience "timeout" errors due to internet connectivity, simply increase your R timeout limit before running the script by executing: `options(timeout = 600)`.
- **Hardware Recommendations:** While the pipeline is highly streamlined, WGCNA topology calculation and triple-consensus machine learning cross-validation are memory-intensive. A standard modern workstation with at least 8GB to 16GB of RAM is recommended for seamless execution.

## Contact
- **Academic & Collaborative Inquiries:** For questions regarding the manuscript, methodology, or potential collaborations, please contact the corresponding authors via email as listed in the publication.
- **Code & Pipeline Issues:** For technical questions, troubleshooting, or code-related bugs, please [open an Issue](https://github.com/RabiaSultan8/ALS_Project/issues) directly in this GitHub repository.
