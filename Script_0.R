# ==============================================================================
# MASTER PACKAGE INSTALLER FOR ALS BIOMARKER PIPELINE
# ==============================================================================
message("Starting master installation...")

# 1. Install BiocManager first if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Define all CRAN Packages
cran_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "stringr", "forcats", "patchwork", 
  "ggsci", "ggrepel", "circlize", "WGCNA", "igraph", "randomForest", 
  "glmnet", "caret", "ggVennDiagram", "pROC", "scales", "viridis", 
  "ggpubr", "ggforce", "pheatmap", "RColorBrewer", "corrr", "enrichR"
)

# 3. Define all Bioconductor Packages
bioc_pkgs <- c(
  "GEOquery", "sva", "illuminaHumanv3.db", "limma", 
  "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db", "GSVA"
)

# 4. Install missing CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing CRAN package:", pkg))
    install.packages(pkg)
  } else {
    message(paste("✔", pkg, "is already installed."))
  }
}

# 5. Install missing Bioconductor packages
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing Bioconductor package:", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    message(paste("✔", pkg, "is already installed."))
  }
}

message("================================================================")
message("All required packages are successfully installed and ready!")
message("================================================================")
