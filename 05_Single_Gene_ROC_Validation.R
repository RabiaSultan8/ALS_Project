# ==============================================================================
# 05_Single_Gene_ROC_Validation.R
# Single-Gene ROC Analysis for Top 6 BPS Biomarkers
# Generates Manuscript Figure 6B components
# ==============================================================================
library(pROC); library(dplyr); library(ggplot2)
set.seed(1122)
dir.create("Manuscript_Figures/Validation/Single_Gene_ROCs", recursive = TRUE, showWarnings = FALSE)
step1_data <- readRDS("Processed_Data/Step1_ComBat_Data.rds")
disc_expr <- as.data.frame(t(step1_data$expr_combat))
disc_labels <- ifelse(step1_data$group_factor == "ALS", 1, 0)
final_scores <- read.csv("Processed_Data/Step8_Final_Composite_Scores.csv")
top_6_genes <- head(final_scores$Gene, 6)
message(paste("Top 6 genes selected for individual ROC analysis:", paste(top_6_genes, collapse = ", ")))
plot_single_roc <- function(gene, expr_matrix, labels_vec, cohort_name, file_suffix) {
  if (!gene %in% colnames(expr_matrix)) { return(NULL) }
  gene_values <- as.numeric(expr_matrix[, gene])
  clean_labels <- as.numeric(ifelse(labels_vec %in% c("ALS", "1"), 1, 0))
  
  roc_obj <- roc(clean_labels, gene_values, quiet = TRUE, ci = TRUE)
  auc_val <- round(as.numeric(roc_obj$auc), 3)
  
  pdf_file <- paste0("Manuscript_Figures/Validation/Single_Gene_ROCs/Figure6B_", gene, "_", file_suffix, ".pdf")
  png_file <- paste0("Manuscript_Figures/Validation/Single_Gene_ROCs/Figure6B_", gene, "_", file_suffix, ".png")
  plot_color <- ifelse(auc_val >= 0.7, "#BC3C29FF", ifelse(auc_val >= 0.6, "#E18727FF", "#0072B5FF"))
  
  pdf(pdf_file, width = 5, height = 5)
  plot(roc_obj, col = plot_color, lwd = 3, main = paste0("6B. ", gene, "\n", cohort_name), cex.main = 1.2, print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.15, print.auc.cex = 1.3)
  dev.off()
  png(png_file, width = 5, height = 5, units = "in", res = 600)
  plot(roc_obj, col = plot_color, lwd = 3, main = paste0("6B. ", gene, "\n", cohort_name), cex.main = 1.2, print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.15, print.auc.cex = 1.3)
  dev.off()
  message(sprintf("Generated ROC for %s in %s (AUC = %s)", gene, cohort_name, auc_val))
}
for (gene in top_6_genes) {
  plot_single_roc(gene, disc_expr, disc_labels, "Discovery Cohort", "Discovery")
  if(exists("val_expr_28253") & exists("labels_28253")) {
      val1_mat_t <- as.data.frame(t(val_expr_28253))
      plot_single_roc(gene, val1_mat_t, labels_28253, "Validation 1 (Microarray)", "GSE28253")
  } 
  if(exists("val_expr_234297_final") & exists("labels_234297")) {
      val2_mat_t <- as.data.frame(t(val_expr_234297_final))
      plot_single_roc(gene, val2_mat_t, labels_234297, "Validation 2 (RNA-Seq)", "GSE234297")
  } 
}

message("Step 9 Single-gene ROC analysis complete. Figure 6B components generated.")
