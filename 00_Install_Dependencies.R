# ==============================================================================
# 00_Install_Dependencies.R
#
# Master package installer for the ALS peripheral blood transcriptomics
# biomarker discovery pipeline. Run this script once before executing any
# other pipeline scripts to ensure all required dependencies are available.
#
# Dependencies are organized by source:
#   - CRAN packages (installed via install.packages)
#   - Bioconductor packages (installed via BiocManager)
#   - GitHub packages (installed via remotes)
#
# Tested on R >= 4.4.0 with Bioconductor >= 3.19
# ==============================================================================

message("Starting dependency installation check...")
message("================================================================")

# ------------------------------------------------------------------------------
# BiocManager — required before any Bioconductor installation
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  install.packages("BiocManager")
}

# ------------------------------------------------------------------------------
# CRAN Packages
# ------------------------------------------------------------------------------
cran_pkgs <- c(
  # Core data manipulation and visualization
  "ggplot2",        # Grammar of graphics plotting
  "dplyr",          # Data manipulation
  "tidyr",          # Data tidying and reshaping
  "stringr",        # String operations
  "forcats",        # Factor handling
  "tibble",         # Modern data frames
  "patchwork",      # Multi-panel figure composition
  "scales",         # Axis and color scale utilities

  # Specialized visualization
  "ggsci",          # Scientific journal color palettes
  "ggrepel",        # Non-overlapping text labels
  "ggpubr",         # Publication-ready figure utilities
  "ggforce",        # Additional ggplot2 extensions
  "ggVennDiagram",  # Venn diagram visualization
  "pheatmap",       # Heatmap generation
  "RColorBrewer",   # Color palette library
  "viridis",        # Perceptually uniform color scales
  "circlize",       # Circular visualization (ComplexHeatmap dependency)
  "corrr",          # Correlation data frames

  # Network and co-expression analysis
  "igraph",         # Network construction and visualization

  # Machine learning
  "randomForest",   # Random forest classification
  "glmnet",         # LASSO and elastic net regression
  "caret",          # Machine learning training framework
  "e1071",          # SVM implementation (required by caret for SVM-RFE)

  # Statistical analysis
  "pROC",           # ROC curve analysis and AUC computation
  "cluster",        # Cluster analysis and silhouette scoring

  # Survival analysis
  "survival",       # Kaplan-Meier and Cox regression
  "survminer",      # Survival curve visualization

  # Drug repurposing and enrichment
  "enrichR",        # Enrichr API interface for gene set enrichment

  # Utility
  "remotes"         # GitHub package installation
)

message("Checking and installing CRAN packages...")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  Installing: %s", pkg))
    install.packages(pkg, quiet = TRUE)
  } else {
    message(sprintf("  [OK] %s", pkg))
  }
}

# ------------------------------------------------------------------------------
# Bioconductor Packages
# ------------------------------------------------------------------------------
bioc_pkgs <- c(
  "GEOquery",           # GEO data retrieval
  "limma",              # Linear models for microarray differential expression
  "sva",                # Surrogate variable analysis (batch correction utilities)
  "illuminaHumanv3.db", # Illumina HumanHT-12 V3.0 probe annotation (GPL6947)
  "illuminaHumanv4.db", # Illumina HumanHT-12 V4.0 probe annotation (GPL10558)
  "ComplexHeatmap",     # Advanced heatmap visualization
  "clusterProfiler",    # GO and pathway enrichment analysis
  "org.Hs.eg.db",       # Human gene annotation database
  "GSVA",               # Gene set variation analysis including ssGSEA
  "WGCNA"               # Weighted gene co-expression network analysis
)

message("Checking and installing Bioconductor packages...")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  Installing: %s", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    message(sprintf("  [OK] %s", pkg))
  }
}

# ------------------------------------------------------------------------------
# GitHub Packages
# Immune deconvolution tools not available on CRAN or Bioconductor.
# Requires an active internet connection and the remotes package.
# ------------------------------------------------------------------------------
message("Checking and installing GitHub packages...")

if (!requireNamespace("EPIC", quietly = TRUE)) {
  message("  Installing: EPIC (Racle et al., eLife 2017)")
  remotes::install_github("GfellerLab/EPIC", build_vignettes = FALSE)
} else {
  message("  [OK] EPIC")
}

if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  message("  Installing: MCPcounter (Becht et al., Genome Biology 2016)")
  remotes::install_github("ebecht/MCPcounter",
                           ref     = "master",
                           subdir  = "Source")
} else {
  message("  [OK] MCPcounter")
}

# ------------------------------------------------------------------------------
# Installation Verification
# ------------------------------------------------------------------------------
message("================================================================")
message("Verifying all installations...")

all_pkgs <- c(cran_pkgs, bioc_pkgs, "EPIC", "MCPcounter")
failed   <- all_pkgs[!sapply(all_pkgs, requireNamespace, quietly = TRUE)]

if (length(failed) == 0) {
  message("All packages successfully installed and verified.")
} else {
  message("WARNING: The following packages could not be installed:")
  for (pkg in failed) message(sprintf("  FAILED: %s", pkg))
  message("Please install failed packages manually before running the pipeline.")
}

message("================================================================")
