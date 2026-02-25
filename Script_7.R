# ==============================================================================
# ADDENDUM SCRIPT: Single-Gene ROC Analysis for Top 6 BPS Biomarkers
# ==============================================================================

# ── 0. Load Required Packages & Data ──────────────────────────────────────────
library(pROC); library(dplyr); library(ggplot2)
set.seed(1122)
dir.create("Manuscript_Figures/Validation/Single_Gene_ROCs", recursive = TRUE, showWarnings = FALSE)

message("Loading previous step data for single-gene analysis...")
# Load Discovery Data
step1_data <- readRDS("Processed_Data/Step1_ComBat_Data.rds")
disc_expr <- as.data.frame(t(step1_data$expr_combat))
disc_labels <- ifelse(step1_data$group_factor == "ALS", 1, 0)

# Load Validation Data (These must exist from a successful run of Script 4)
# We need to re-load the processed validation matrices from Step 7 environment. 
# NOTE: If you closed R, you need to re-run the data loading parts of Script 4 first.
# Assuming val_expr_28253, labels_28253, val_expr_234297_final, labels_234297 exist in memory.

# Load Final Scores to get Top 6
final_scores <- read.csv("Processed_Data/Step8_Final_Composite_Scores.csv")
top_6_genes <- head(final_scores$Gene, 6)
message(paste("Top 6 genes selected for individual ROC analysis:", paste(top_6_genes, collapse = ", ")))

# ── 1. Define Generic ROC Plotting Function ───────────────────────────────────
plot_single_roc <- function(gene, expr_matrix, labels_vec, cohort_name, file_suffix) {
  
  # Ensure gene exists in data
  if (!gene %in% colnames(expr_matrix)) {
    message(sprintf("Skipping %s in %s (gene not found)", gene, cohort_name))
    return(NULL)
  }
  
  # Prepare data & ensure numeric labels
  gene_values <- as.numeric(expr_matrix[, gene])
  clean_labels <- as.numeric(ifelse(labels_vec %in% c("ALS", "1"), 1, 0))
  
  # Calculate ROC
  roc_obj <- roc(clean_labels, gene_values, quiet = TRUE, ci = TRUE)
  auc_val <- round(as.numeric(roc_obj$auc), 3)
  
  # Define file names
  pdf_file <- paste0("Manuscript_Figures/Validation/Single_Gene_ROCs/", gene, "_", file_suffix, ".pdf")
  png_file <- paste0("Manuscript_Figures/Validation/Single_Gene_ROCs/", gene, "_", file_suffix, ".png")
  
  # Plot Title based on AUC performance
  plot_color <- ifelse(auc_val >= 0.7, "#BC3C29FF", ifelse(auc_val >= 0.6, "#E18727FF", "#0072B5FF"))
  
  # Save PDF
  pdf(pdf_file, width = 5, height = 5)
  plot(roc_obj, col = plot_color, lwd = 3, 
       main = paste0(gene, "\n", cohort_name),
       cex.main = 1.2,
       print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.15, print.auc.cex = 1.3)
  dev.off()
  
  # Save PNG
  png(png_file, width = 5, height = 5, units = "in", res = 600)
  plot(roc_obj, col = plot_color, lwd = 3, 
       main = paste0(gene, "\n", cohort_name),
       cex.main = 1.2,
       print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.15, print.auc.cex = 1.3)
  dev.off()
  
  message(sprintf("Generated ROC for %s in %s (AUC = %s)", gene, cohort_name, auc_val))
}

# ── 2. Execute Loop for Top 6 Genes across All Cohorts ────────────────────────
message("\n--- Generating Single-Gene ROC Curves ---")

for (gene in top_6_genes) {
  # 1. Discovery Cohort
  plot_single_roc(gene, disc_expr, disc_labels, "Discovery Cohort", "Discovery")
  
  # 2. Validation Cohort 1 (GSE28253)
  # Ensure these variables exist from Script 4 run
  if(exists("val_expr_28253") & exists("labels_28253")) {
      val1_mat_t <- as.data.frame(t(val_expr_28253))
      plot_single_roc(gene, val1_mat_t, labels_28253, "Validation 1 (Microarray)", "GSE28253")
  } else { message("Skipping Validation 1 (Data not found in memory). Run Script 4 first.") }

  # 3. Validation Cohort 2 (GSE234297)
  # Ensure these variables exist from Script 4 run
  if(exists("val_expr_234297_final") & exists("labels_234297")) {
      val2_mat_t <- as.data.frame(t(val_expr_234297_final))
      plot_single_roc(gene, val2_mat_t, labels_234297, "Validation 2 (RNA-Seq)", "GSE234297")
  } else { message("Skipping Validation 2 (Data not found in memory). Run Script 4 first.") }
}

message("================================================================")
message("Single-gene ROC analysis complete. Check 'Validation/Single_Gene_ROCs' folder.")
message("================================================================")
