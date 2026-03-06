# ==============================================================================
# 04_Clinical_Validation_and_BPS.R
# Clinical Validation & Biomarker Scoring (BPS)
# Generates components for Manuscript Figure 6
# ==============================================================================

# ── 0. Load Required Libraries ────────────────────────────────────────────────
library(GEOquery); library(dplyr); library(pROC); library(ggplot2); library(ggforce); library(ggpubr); library(clusterProfiler); library(org.Hs.eg.db)

dir.create("Manuscript_Figures/Validation", recursive = TRUE, showWarnings = FALSE)
set.seed(1122)

# ==============================================================================
# Step 6: ML Model Construction & Discovery Cohort Evaluation
# ==============================================================================
expr_data <- readRDS("Processed_Data/Step1_ComBat_Data.rds")
consensus_genes <- read.csv("Processed_Data/Final_Gene_Signature.csv")$Gene

message(sprintf("Building model using %d consensus signature genes...", length(consensus_genes)))

disc_data_raw <- as.data.frame(t(expr_data$expr_combat[consensus_genes, ]))
disc_data <- as.data.frame(scale(disc_data_raw)) 
disc_data$Diagnosis <- ifelse(expr_data$group_factor == "ALS", 1, 0)

logit_model <- glm(Diagnosis ~ ., data = disc_data, family = "binomial")
disc_data$RiskScore <- predict(logit_model, type = "response")

# FIGURE 6A (Part 1): Discovery ROC Curve
roc_disc <- roc(disc_data$Diagnosis, disc_data$RiskScore, ci = TRUE)
png("Manuscript_Figures/Validation/Figure6A_Discovery_ROC.png", width = 6, height = 6, units = "in", res = 600)
plot(roc_disc, col = "#BC3C29FF", lwd = 3, main = "6A. ROC Curve: Discovery Cohort (n=1,042)", print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.2, print.auc.cex = 1.2)
dev.off()
pdf("Manuscript_Figures/Validation/Figure6A_Discovery_ROC.pdf", width = 6, height = 6)
plot(roc_disc, col = "#BC3C29FF", lwd = 3, main = "6A. ROC Curve: Discovery Cohort (n=1,042)", print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.2, print.auc.cex = 1.2)
dev.off()

# FIGURE 6C: Statistical Violin Plot of Risk Scores
disc_data$Group <- factor(ifelse(disc_data$Diagnosis == 1, "ALS", "Control"), levels = c("Control", "ALS"))

p_violin <- ggplot(disc_data, aes(x = Group, y = RiskScore, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  scale_fill_manual(values = c("ALS" = "#BC3C29FF", "Control" = "#0072B5FF")) +
  theme_classic(base_size = 15) +
  stat_compare_means(method = "wilcox.test", label.y = 1.1, size = 5, fontface = "bold") + 
  labs(title = "6C. Diagnostic Risk Score Distribution", y = "Predicted Probability of ALS", x = "") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("Manuscript_Figures/Validation/Figure6C_RiskScore_Violin.png", plot = p_violin, width = 5, height = 6, dpi = 600)
ggsave("Manuscript_Figures/Validation/Figure6C_RiskScore_Violin.pdf", plot = p_violin, width = 5, height = 6)

evaluate_validation <- function(val_expr, val_labels, cohort_name, plot_prefix, model, sig_genes) {
  available_genes <- intersect(sig_genes, rownames(val_expr))
  val_df_raw <- as.data.frame(t(val_expr[available_genes, , drop=FALSE]))
  val_df <- as.data.frame(scale(val_df_raw))
  missing_genes <- setdiff(sig_genes, available_genes)
  for (gene in missing_genes) { val_df[[gene]] <- 0 }
  val_df$Diagnosis <- ifelse(val_labels == "ALS", 1, 0)
  val_df$RiskScore <- predict(model, newdata = val_df, type = "response")
  val_roc <- roc(val_df$Diagnosis, val_df$RiskScore, ci = TRUE, quiet = TRUE)
  
  pdf(paste0("Manuscript_Figures/Validation/", plot_prefix, ".pdf"), width = 6, height = 6)
  plot(val_roc, col = "#0072B5FF", lwd = 3, main = paste("ROC Curve:", cohort_name), print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.2, print.auc.cex = 1.2)
  dev.off()
  png(paste0("Manuscript_Figures/Validation/", plot_prefix, ".png"), width = 6, height = 6, units = "in", res = 600)
  plot(val_roc, col = "#0072B5FF", lwd = 3, main = paste("ROC Curve:", cohort_name), print.auc = TRUE, print.auc.x = 0.4, print.auc.y = 0.2, print.auc.cex = 1.2)
  dev.off()
  message(sprintf("✔ Successfully validated %s! AUC = %.3f", cohort_name, val_roc$auc))
  return(val_roc)
}

saveRDS(logit_model, "Processed_Data/Step6_Logit_Model.rds")
saveRDS(evaluate_validation, "Processed_Data/Step6_Validation_Function.rds")

message("Step 6 Complete. Discovery model locked and saved.")

# ==============================================================================
# Step 7: Independent Cohort Validation (Generates Figure 6A parts 2 & 3)
# ==============================================================================
message("Downloading Validation Cohort 1: GSE28253 (Microarray)...")
gse28253_list <- getGEO("GSE28253", getGPL = TRUE)
gse28253 <- gse28253_list[[1]]
expr_28253 <- exprs(gse28253); pdata_28253 <- pData(gse28253)
all_clinical_text <- paste(pdata_28253$title, pdata_28253$source_name_ch1, pdata_28253$characteristics_ch1)
labels_28253 <- ifelse(grepl("ALS|Amyotrophic", all_clinical_text, ignore.case = TRUE), "ALS", ifelse(grepl("Control|Healthy", all_clinical_text, ignore.case = TRUE), "Control", NA))
message(sprintf("Found %d ALS and %d Control patients in GSE28253.", sum(labels_28253 == "ALS", na.rm=TRUE), sum(labels_28253 == "Control", na.rm=TRUE)))

fd_28253 <- fData(gse28253)
sym_col_name <- grep("GENE_SYMBOL|GeneSymbol|Symbol", colnames(fd_28253), ignore.case = TRUE, value = TRUE)[1]
if(!is.na(sym_col_name)) { gene_symbols <- fd_28253[[sym_col_name]] } else { gene_symbols <- rownames(fd_28253) }
gene_symbols_clean <- sapply(strsplit(as.character(gene_symbols), " |/"), `[`, 1)

valid_probes <- !is.na(gene_symbols_clean) & gene_symbols_clean != "" & gene_symbols_clean != "NA"
expr_28253_clean <- expr_28253[valid_probes, ]; rownames(expr_28253_clean) <- gene_symbols_clean[valid_probes]
df_28253 <- as.data.frame(expr_28253_clean); df_28253$gene <- rownames(df_28253)
df_28253_agg <- df_28253 %>% group_by(gene) %>% summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
val_expr_28253 <- as.matrix(df_28253_agg[, -1]); rownames(val_expr_28253) <- df_28253_agg$gene

message(sprintf("GSE28253 successfully matched %d out of %d signature genes!", length(intersect(consensus_genes, rownames(val_expr_28253))), length(consensus_genes)))
roc_28253 <- evaluate_validation(val_expr = val_expr_28253, val_labels = labels_28253, cohort_name = "GSE28253 (Independent Microarray)", plot_prefix = "Figure6A_Validation_ROC_GSE28253", model = logit_model, sig_genes = consensus_genes)

message("\nDownloading Validation Cohort 2: GSE234297 (RNA-Seq Clinical Data)...")
gse234297_list <- getGEO("GSE234297", getGPL = FALSE)
gse234297 <- gse234297_list[[1]]; pdata_234297 <- pData(gse234297)

if("disease state:ch1" %in% colnames(pdata_234297)) {
  labels_raw <- pdata_234297$`disease state:ch1`
  labels_234297 <- ifelse(grepl("ALS|Amyotrophic", labels_raw, ignore.case = TRUE), "ALS", ifelse(grepl("Control|Healthy", labels_raw, ignore.case = TRUE), "Control", NA))
} else {
  all_clinical_text_23 <- paste(pdata_234297$title, pdata_234297$source_name_ch1, pdata_234297$characteristics_ch1)
  labels_234297 <- ifelse(grepl("ALS|sALS|Amyotrophic", all_clinical_text_23, ignore.case = TRUE), "ALS", ifelse(grepl("Control|Healthy", all_clinical_text_23, ignore.case = TRUE), "Control", NA))
}
message(sprintf("Found %d ALS and %d Control patients in GSE234297.", sum(labels_234297 == "ALS", na.rm=TRUE), sum(labels_234297 == "Control", na.rm=TRUE)))

message("Locating/Downloading RNA-Seq Count Matrix...")
if (!dir.exists("GSE234297") || length(list.files("GSE234297", pattern = "raw_counts")) == 0) { suppressMessages(getGEOSuppFiles("GSE234297", makeDirectory = TRUE)) }
file_paths <- list.files("GSE234297", full.names = TRUE)
counts_file <- grep("raw_counts", file_paths, value = TRUE, ignore.case = TRUE)[1]

if(!is.na(counts_file)) {
  val_expr_234297 <- read.table(gzfile(counts_file), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  raw_rownames <- as.character(rownames(val_expr_234297))
  if (grepl("^[0-9]+$", raw_rownames[1])) {
    message("Detected ENTREZ IDs. Translating to standard Gene Symbols...")
    gene_trans <- suppressMessages(bitr(raw_rownames, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db))
    val_expr_234297$ENTREZID <- raw_rownames
    val_expr_merged <- merge(gene_trans, val_expr_234297, by = "ENTREZID")
    val_expr_merged <- val_expr_merged %>% mutate(across(where(is.numeric), ~ log2(.x + 1)))
    df_23 <- as.data.frame(val_expr_merged)
    df_23_agg <- df_23 %>% group_by(SYMBOL) %>% summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
    val_expr_234297_final <- as.matrix(df_23_agg[, -1]); rownames(val_expr_234297_final) <- df_23_agg$SYMBOL
  } else { val_expr_234297_final <- log2(as.matrix(val_expr_234297) + 1) }

  message(sprintf("GSE234297 successfully matched %d out of %d signature genes!", length(intersect(consensus_genes, rownames(val_expr_234297_final))), length(consensus_genes)))
  roc_234297 <- evaluate_validation(val_expr = val_expr_234297_final, val_labels = labels_234297, cohort_name = "GSE234297 (Independent RNA-Seq)", plot_prefix = "Figure6A_Validation_ROC_GSE234297", model = logit_model, sig_genes = consensus_genes)
}

message("Step 7 Complete. External cohort validation done.")

# ==============================================================================
# Step 8: Composite Biomarker Priority Scoring (BPS)
# ==============================================================================
disc_expr_t <- as.data.frame(t(expr_data$expr_combat[consensus_genes, ]))
val1_mat <- as.data.frame(t(val_expr_28253))
val2_mat <- as.data.frame(t(val_expr_234297_final))
disc_labels <- ifelse(expr_data$group_factor == "ALS", 1, 0)

get_safe_auc <- function(gene, expr_matrix, labels) {
  if (gene %in% colnames(expr_matrix)) {
    vals <- as.numeric(expr_matrix[, gene])
    tryCatch({ r <- roc(labels, vals, quiet = TRUE); return(as.numeric(r$auc)) }, error = function(e) return(NA))
  } else { return(NA) }
}

wgcna_data <- readRDS("Processed_Data/Step3_WGCNA_Data.rds")
top_mod <- wgcna_data$top_module_colors[1]
kME_matrix <- as.data.frame(cor(wgcna_data$datExpr, wgcna_data$MEs, use = "p"))
kME_col <- paste0("ME", top_mod)
results_list <- list()

for (gene in consensus_genes) {
  auc_disc <- get_safe_auc(gene, disc_expr_t, disc_labels)
  auc_val1 <- get_safe_auc(gene, val1_mat, ifelse(labels_28253 == "ALS", 1, 0))
  auc_val2 <- get_safe_auc(gene, val2_mat, ifelse(labels_234297 == "ALS", 1, 0))
  mean_auc <- mean(c(auc_disc, auc_val1, auc_val2), na.rm = TRUE)
  kme_val <- NA
  if (gene %in% rownames(kME_matrix)) { kme_val <- abs(kME_matrix[gene, kME_col]) }
  results_list[[gene]] <- data.frame(Gene = gene, AUC_Microarray = round(auc_val1, 3), AUC_RNASeq = round(auc_val2, 3), Mean_AUC = mean_auc, WGCNA_kME = kme_val)
}

final_df <- do.call(rbind, results_list)
min_max_scale <- function(x) { (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) }

final_df <- final_df %>%
  filter(!is.na(WGCNA_kME) & !is.na(Mean_AUC)) %>%
  mutate(Scaled_AUC = min_max_scale(Mean_AUC), Scaled_kME = min_max_scale(WGCNA_kME), Composite_Score = (Scaled_AUC * 0.5) + (Scaled_kME * 0.5), Final_BPS_100 = round(Composite_Score * 100, 1)) %>% arrange(desc(Final_BPS_100))

message("\n=======================================================")
message("🏆 TOP BIOMARKERS BY COMPOSITE PRIORITY SCORE (BPS) 🏆")
print(head(final_df[, c("Gene", "AUC_Microarray", "AUC_RNASeq", "WGCNA_kME", "Final_BPS_100")], 10))
message("=======================================================\n")

write.csv(final_df, "Processed_Data/Step8_Final_Composite_Scores.csv", row.names = FALSE)

message("Step 8 Complete. Figure 6 components generated.")
