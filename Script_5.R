# ==============================================================================
# MASTER SCRIPT 5: Immune Infiltration Analysis (Step 9) - CORRECTED
# ==============================================================================

# ── 0. Load Required Packages ────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GSVA", quietly = TRUE)) BiocManager::install("GSVA", update = FALSE, ask = FALSE)

packages <- c("GSVA", "pheatmap", "ggplot2", "dplyr", "tidyr", "ggpubr", "RColorBrewer")
for (pkg in packages) { 
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE)) 
}

set.seed(1122)
dir.create("Manuscript_Figures/Immune_Analysis", recursive = TRUE, showWarnings = FALSE)
message("Starting Step 9: Immune Infiltration Analysis (ssGSEA)...")

# 1. Load Discovery Data & The Core 6-Gene Signature
expr_data <- readRDS("Processed_Data/Step1_ComBat_Data.rds")
core_genes <- c("XPO1", "PPP1CC", "DHX15", "EIF4G2", "ZFAND5", "FBXO11")

expr_mat <- expr_data$expr_combat
group_factor <- expr_data$group_factor

# 2. Define a Curated Immune Cell Marker Dictionary
immune_signatures <- list(
  `B cells` = c("CD19", "CD79A", "CD79B", "MS4A1", "BANK1"),
  `CD4+ T cells` = c("CD4", "CD40LG", "CXCR3", "CCR4", "IL2RA"),
  `CD8+ T cells` = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG"),
  `NK cells` = c("NCAM1", "KLRB1", "KLRC1", "KLRD1", "KLRF1", "NKG7"),
  `Macrophages` = c("CD68", "CD163", "MRC1", "MSR1", "CHIT1"),
  `Neutrophils` = c("CEACAM8", "FCGR3B", "CSF3R", "S100A12", "CXCR2"),
  `Dendritic cells` = c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "CD1C"),
  `Monocytes` = c("CD14", "FCGR3A", "CD33", "ITGAM", "CCR2"),
  `Tregs` = c("FOXP3", "IL2RA", "CTLA4", "CCR8", "IKZF2"),
  `Mast cells` = c("TPSAB1", "TPSB2", "CPA3", "HDC", "MS4A2")
)

# 3. Run ssGSEA to estimate Immune Cell Abundance
message("Running ssGSEA to estimate immune cell fractions...")
tryCatch({
  gsva_param <- gsvaParam(expr_mat, immune_signatures, maxDiff = TRUE)
  immune_scores <- gsva(gsva_param)
}, error = function(e) {
  immune_scores <<- gsva(expr_mat, immune_signatures, method = "ssgsea", verbose = FALSE)
})

immune_df <- as.data.frame(t(immune_scores))
immune_df$Diagnosis <- factor(ifelse(group_factor == "ALS", "ALS", "Control"), levels = c("Control", "ALS"))

# 4. Figure 9A: Boxplots of Immune Cell Abundance
immune_long <- immune_df %>%
  pivot_longer(cols = -Diagnosis, names_to = "Cell_Type", values_to = "ssGSEA_Score")

p_immune_box <- ggplot(immune_long, aes(x = reorder(Cell_Type, ssGSEA_Score, FUN = median), y = ssGSEA_Score, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  scale_fill_manual(values = c("ALS" = "#BC3C29FF", "Control" = "#0072B5FF")) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
        legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5)) +
  stat_compare_means(aes(group = Diagnosis), label = "p.signif", method = "wilcox.test", hide.ns = TRUE) +
  labs(title = "A. Immune Cell Infiltration Profile (ssGSEA)", x = "", y = "Enrichment Score")

ggsave("Manuscript_Figures/Immune_Analysis/Figure9A_Immune_Boxplots.png", plot = p_immune_box, width = 10, height = 6, dpi = 600)
ggsave("Manuscript_Figures/Immune_Analysis/Figure9A_Immune_Boxplots.pdf", plot = p_immune_box, width = 10, height = 6)

# 5. Figure 9B: Correlation Heatmap between Core Genes and Immune Cells
message("Calculating correlations between 6 Core Genes and Immune Cells...")

core_expr <- t(expr_mat[core_genes, ])

# CRITICAL FIX: Use base R subsetting to prevent package masking errors
immune_numeric <- immune_df[, colnames(immune_df) != "Diagnosis"]
cor_mat <- cor(x = core_expr, y = immune_numeric, method = "spearman")

col_fun <- colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100)

# CRITICAL FIX: Split long lines to prevent console copy-paste truncation
png_file <- "Manuscript_Figures/Immune_Analysis/Figure9B_Correlation.png"
png(png_file, width = 8, height = 6, units = "in", res = 600)
pheatmap(cor_mat, 
         color = col_fun,
         display_numbers = round(cor_mat, 2),
         fontsize_number = 10,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         angle_col = 45,
         main = "B. Correlation: Core Biomarkers vs. Immune Cells")
dev.off()

pdf_file <- "Manuscript_Figures/Immune_Analysis/Figure9B_Correlation.pdf"
pdf(pdf_file, width = 8, height = 6)
pheatmap(cor_mat, 
         color = col_fun,
         display_numbers = round(cor_mat, 2),
         fontsize_number = 10,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         angle_col = 45,
         main = "B. Correlation: Core Biomarkers vs. Immune Cells")
dev.off()

message("================================================================")
message("Step 9 Complete! Immune boxplots and correlation heatmap saved.")
message("================================================================")
