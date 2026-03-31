# ==============================================================================
# 01_Data_Prep_and_DEG.R
#
# Data preparation and differential expression analysis.
# Discovery cohort: GSE112676 only (Illumina HumanHT-12 V3.0, GPL6947).
# GSE112680 is reserved entirely as the held-out validation cohort and is
# never loaded or processed in this script.
#
# Key methodological decisions documented inline:
#   - Best-probe selection by IQR (superior to multi-probe averaging)
#   - Adaptive normalization with QC-driven decision
#   - Sample outlier detection before DEG analysis
#   - Conservative |logFC| >= 0.25 threshold justified for blood transcriptomics
#
# Outputs:
#   Processed_Data/Step1_Discovery_Data.rds
#   Processed_Data/Step2_DEG_all.csv
#   Processed_Data/Step2_DEG_sig.csv
#   Manuscript_Figures/Figure2A_Discovery_PCA.{pdf,png}
#   Manuscript_Figures/Figure2B_Volcano.{pdf,png}
#   Manuscript_Figures/Figure2C_Heatmap.{pdf,png}
#   Manuscript_Figures/Figure2D_MA_Plot.{pdf,png}
#   Manuscript_Figures/QC/Sample_Outlier_Detection.{pdf,png}
#   Manuscript_Figures/QC/Sample_Distribution_Boxplot.{pdf,png}
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(GEOquery)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(illuminaHumanv3.db)
library(patchwork)
library(ggsci)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

dir.create("Manuscript_Figures/QC",        recursive = TRUE,  showWarnings = FALSE)
dir.create("Manuscript_Figures",            showWarnings = FALSE)
dir.create("Processed_Data",               showWarnings = FALSE)

# ==============================================================================
# STEP 1: Data Acquisition
# ==============================================================================
message("Loading GSE112676 (discovery cohort only)...")

# Load from local cache for reproducibility. To download and cache:
#   gse <- getGEO("GSE112676", getGPL = FALSE)[[1]]
#   saveRDS(gse, "GSE112676_cached.rds")
# or use:
#   getGEO(filename = "GSE112676_series_matrix.txt.gz", getGPL = FALSE)

if (file.exists("GSE112676_series_matrix.txt.gz")) {
  gse112676 <- getGEO(filename = "GSE112676_series_matrix.txt.gz",
                      getGPL   = FALSE)
} else {
  message("Local cache not found. Downloading GSE112676 from GEO...")
  gse112676 <- getGEO("GSE112676", getGPL = FALSE)[[1]]
}

expr_raw  <- exprs(gse112676)
pheno     <- pData(gse112676)
group_raw <- toupper(trimws(as.character(pheno$`diagnosis:ch1`)))

# GSE112676 contains only ALS and CON, but the filter is explicit for safety
keep      <- group_raw %in% c("ALS", "CON")
expr_raw  <- expr_raw[, keep]
group_raw <- group_raw[keep]

message(sprintf("Samples retained: %d total | ALS: %d | CON: %d",
                sum(keep),
                sum(group_raw == "ALS"),
                sum(group_raw == "CON")))

# ==============================================================================
# STEP 2: Probe-to-Gene Mapping — Best Probe by IQR
#
# Standard practice (used in both previous pipelines) averages all probes
# mapping to the same gene. This is suboptimal: it dilutes the signal from
# the most informative probe with noise from less-variable probes.
#
# Here we select the single probe per gene with the highest inter-quartile
# range (IQR) across samples. The highest-IQR probe captures the maximum
# biological variance for that gene and is the standard recommended by the
# WGCNA authors (Langfelder & Horvath, 2008) for BeadChip preprocessing.
# ==============================================================================
message("Mapping probes to gene symbols — selecting best probe per gene by IQR...")

probes  <- rownames(expr_raw)
symbols <- mapIds(
  illuminaHumanv3.db,
  keys      = probes,
  column    = "SYMBOL",
  keytype   = "PROBEID",
  multiVals = "first"
)

valid        <- !is.na(symbols) & symbols != ""
expr_valid   <- expr_raw[valid, ]
syms_valid   <- symbols[valid]

message(sprintf("Probe mapping: %d / %d probes successfully mapped to HGNC symbols.",
                sum(valid), length(valid)))

# Apply log2 transform BEFORE IQR calculation — IQR on raw intensities
# is dominated by high-expression genes and does not reflect relative variance
if (max(expr_valid, na.rm = TRUE) > 100) {
  message("Linear-scale intensities detected — applying log2(x + 1) transformation.")
  expr_valid <- log2(expr_valid + 1)
}

# Select best probe per gene by IQR
probe_iqr  <- apply(expr_valid, 1, IQR, na.rm = TRUE)
probe_df   <- data.frame(
  probe  = rownames(expr_valid),
  symbol = syms_valid,
  iqr    = probe_iqr,
  stringsAsFactors = FALSE
)

best_probes <- probe_df %>%
  group_by(symbol) %>%
  slice_max(order_by = iqr, n = 1, with_ties = FALSE) %>%
  pull(probe)

expr_discovery           <- expr_valid[best_probes, ]
rownames(expr_discovery) <- probe_df$symbol[match(best_probes, probe_df$probe)]

message(sprintf("Best-probe selection complete: %d unique genes retained.",
                nrow(expr_discovery)))

# ==============================================================================
# STEP 3: Adaptive Normalization with QC Decision
#
# GEO-deposited Illumina BeadChip data is frequently already quantile-
# normalized by BeadStudio/GenomeStudio prior to submission. Blindly applying
# quantile normalization on pre-normalized data is methodologically unjustified.
#
# Decision logic:
#   Compute the coefficient of variation (CV) of per-sample medians.
#   If CV > 0.02 (>2% relative variability in sample medians), the data
#   is heterogeneous enough to warrant between-array normalization.
#   Otherwise, document that the data is already in a normalized state.
#
# In both cases, a sample distribution boxplot is generated as a QC record.
# ==============================================================================
message("Assessing sample distribution for normalization decision...")

sample_medians <- apply(expr_discovery, 2, median, na.rm = TRUE)
cv_medians     <- sd(sample_medians) / mean(sample_medians)

message(sprintf(
  "Sample median CV = %.4f (threshold = 0.02 | %s normalization)",
  cv_medians,
  ifelse(cv_medians > 0.02, "APPLYING", "SKIPPING")
))

# Generate QC boxplot BEFORE any normalization decision
box_df_before <- as.data.frame(expr_discovery) %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(-gene, names_to = "Sample", values_to = "Expression") %>%
  mutate(Group = ifelse(group_raw[match(Sample, colnames(expr_discovery))] == "ALS",
                        "ALS", "Control"))

# Subsample for speed — boxplot with 500+ samples is unreadable anyway
set.seed(1122)
sample_subset <- sample(colnames(expr_discovery), min(60, ncol(expr_discovery)))
box_df_sub <- box_df_before %>% filter(Sample %in% sample_subset)

p_box_before <- ggplot(box_df_sub,
                        aes(x = reorder(Sample, Expression, FUN = median),
                            y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  scale_fill_manual(values = c("ALS" = "#DC143C", "Control" = "#008B8B")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title    = element_text(face = "bold", hjust = 0.5)) +
  labs(title    = "Sample Distribution (Pre-Normalization Check)",
       subtitle = sprintf("n = %d samples shown | CV of medians = %.4f",
                          length(sample_subset), cv_medians),
       x = "Samples", y = "log2 Expression")

if (cv_medians > 0.02) {
  expr_discovery <- normalizeBetweenArrays(expr_discovery, method = "quantile")
  message("Quantile normalization applied.")

  box_df_after <- as.data.frame(expr_discovery[, sample_subset]) %>%
    mutate(gene = rownames(.)) %>%
    pivot_longer(-gene, names_to = "Sample", values_to = "Expression") %>%
    mutate(Group = ifelse(group_raw[match(Sample, colnames(expr_discovery))] == "ALS",
                          "ALS", "Control"))

  p_box_after <- ggplot(box_df_after,
                         aes(x = reorder(Sample, Expression, FUN = median),
                             y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
    scale_fill_manual(values = c("ALS" = "#DC143C", "Control" = "#008B8B")) +
    theme_classic(base_size = 11) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title    = element_text(face = "bold", hjust = 0.5)) +
    labs(title = "Sample Distribution (Post-Normalization)",
         x = "Samples", y = "log2 Expression")

  p_box_combined <- p_box_before / p_box_after
} else {
  message("Data already normalized — quantile normalization skipped.")
  p_box_combined <- p_box_before +
    labs(subtitle = sprintf(
      "n = %d samples | CV = %.4f — Data is pre-normalized, no adjustment applied.",
      length(sample_subset), cv_medians))
}

ggsave("Manuscript_Figures/QC/Sample_Distribution_Boxplot.pdf",
       plot = p_box_combined, width = 14, height = ifelse(cv_medians > 0.02, 10, 6))
ggsave("Manuscript_Figures/QC/Sample_Distribution_Boxplot.png",
       plot = p_box_combined, width = 14, height = ifelse(cv_medians > 0.02, 10, 6), dpi = 300)
message("QC boxplot saved.")

# ==============================================================================
# STEP 4: Sample Outlier Detection
#
# Neither previous pipeline performs this step. A single outlier sample can
# disproportionately influence limma's empirical Bayes variance shrinkage,
# artificially inflating or deflating test statistics for hundreds of genes.
#
# Method: Compute inter-sample Pearson correlation matrix on the full
# expression matrix. Calculate each sample's mean correlation with all others
# (its "connectivity"). Samples with connectivity < mean - 2*SD are flagged
# as outliers. Removal is automatic but reported explicitly in the output.
#
# A connectivity plot is generated as a QC figure regardless of whether
# outliers are found.
# ==============================================================================
message("Running sample outlier detection...")

cor_mat   <- cor(expr_discovery, method = "pearson", use = "pairwise.complete.obs")
mean_conn <- rowMeans(cor_mat, na.rm = TRUE)
conn_mean <- mean(mean_conn)
conn_sd   <- sd(mean_conn)
threshold <- conn_mean - 2 * conn_sd
outliers  <- names(mean_conn[mean_conn < threshold])

conn_df <- data.frame(
  Sample      = names(mean_conn),
  Connectivity = mean_conn,
  Group       = ifelse(group_raw == "ALS", "ALS", "Control"),
  Outlier     = names(mean_conn) %in% outliers
)

p_conn <- ggplot(conn_df,
                  aes(x = seq_along(Sample), y = Connectivity,
                      color = Group, shape = Outlier)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_hline(yintercept = threshold,
             linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c("ALS" = "#DC143C", "Control" = "#008B8B")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4),
                     labels = c("Normal", "Outlier")) +
  annotate("text", x = 1, y = threshold - 0.003,
           label = sprintf("Threshold = %.4f (mean - 2SD)", threshold),
           hjust = 0, color = "red", size = 3.5) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title    = "Sample Connectivity Outlier Detection",
       subtitle = sprintf("n = %d samples | Outliers flagged: %d",
                          ncol(expr_discovery), length(outliers)),
       x = "Sample Index", y = "Mean Inter-Sample Correlation")

ggsave("Manuscript_Figures/QC/Sample_Outlier_Detection.pdf",
       plot = p_conn, width = 11, height = 5)
ggsave("Manuscript_Figures/QC/Sample_Outlier_Detection.png",
       plot = p_conn, width = 11, height = 5, dpi = 600)

if (length(outliers) > 0) {
  message(sprintf("WARNING: %d outlier sample(s) detected and removed:", length(outliers)))
  for (s in outliers) {
    message(sprintf("  Removed: %s (connectivity = %.4f)", s, mean_conn[s]))
  }
  keep_samples   <- !colnames(expr_discovery) %in% outliers
  expr_discovery <- expr_discovery[, keep_samples]
  group_raw      <- group_raw[keep_samples]
  message(sprintf("Samples after outlier removal: %d", ncol(expr_discovery)))
} else {
  message("No outlier samples detected. All samples retained.")
}

group_factor <- factor(
  ifelse(group_raw == "ALS", "ALS", "Control"),
  levels = c("Control", "ALS")
)

message(sprintf("Final discovery matrix: %d genes x %d samples | ALS: %d | CON: %d",
                nrow(expr_discovery), ncol(expr_discovery),
                sum(group_factor == "ALS"),
                sum(group_factor == "Control")))

# ==============================================================================
# STEP 5: Figure 2A — Discovery Cohort PCA
# ==============================================================================
message("Generating Figure 2A (PCA)...")

pca_res <- prcomp(t(expr_discovery), scale. = TRUE)
pct_var <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

df_pca <- data.frame(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  Diagnosis = group_factor
)

pub_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text       = element_text(color = "black"),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 11)
  )

fig2A <- ggplot(df_pca, aes(x = PC1, y = PC2,
                              color = Diagnosis, fill = Diagnosis)) +
  geom_point(alpha = 0.8, size = 2.2) +
  scale_color_manual(values = c("ALS" = "#DC143C", "Control" = "#008B8B")) +
  scale_fill_manual( values = c("ALS" = "#DC143C", "Control" = "#008B8B")) +
  pub_theme +
  labs(
    title    = "Discovery Cohort — Principal Component Analysis",
    subtitle = sprintf("GSE112676 | n = %d samples (ALS: %d, CON: %d)",
                       ncol(expr_discovery),
                       sum(group_factor == "ALS"),
                       sum(group_factor == "Control")),
    x = sprintf("PC1 (%s%%)", pct_var[1]),
    y = sprintf("PC2 (%s%%)", pct_var[2])
  )

# Confidence ellipses — drawn only when n is sufficient per group
if (min(table(group_factor)) >= 4) {
  fig2A <- fig2A +
    stat_ellipse(geom = "polygon", alpha = 0.08, type = "norm", linewidth = 0.5)
}

ggsave("Manuscript_Figures/Figure2A_Discovery_PCA.pdf",
       plot = fig2A, width = 7, height = 6)
ggsave("Manuscript_Figures/Figure2A_Discovery_PCA.png",
       plot = fig2A, width = 7, height = 6, dpi = 600)
message("Figure 2A saved.")

# ==============================================================================
# STEP 6: Save Discovery Data
# ==============================================================================
saveRDS(
  list(
    expr_discovery = expr_discovery,
    group_factor   = group_factor,
    dataset_id     = "GSE112676",
    n_samples      = ncol(expr_discovery),
    n_als          = sum(group_factor == "ALS"),
    n_con          = sum(group_factor == "Control"),
    outliers_removed = outliers,
    normalized     = cv_medians > 0.02
  ),
  "Processed_Data/Step1_Discovery_Data.rds"
)
message("Step 1 complete. Discovery data saved.")

# ==============================================================================
# STEP 7: Differential Expression Analysis — Figure 2B and 2C
#
# limma empirical Bayes framework with a binary ALS vs Control contrast.
# BH FDR correction across all tested genes.
#
# |logFC| >= 0.25 threshold is intentionally conservative.
# Blood transcriptomic signals in neurodegenerative disease are characteristically
# modest — leukocytes are indirect reporters of CNS-centred pathology.
# Stringent fold-change cutoffs (e.g., >= 1.0) are appropriate for tissue-level
# studies but systematically discard biologically relevant peripheral signals.
# The 0.25 threshold is consistent with published ALS blood transcriptomics
# (van Rheenen et al., 2018; Zhao et al., 2025).
# ==============================================================================
message("Running limma differential expression analysis...")

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Control", "ALS")

fit   <- lmFit(expr_discovery, design)
contr <- makeContrasts(ALS_vs_Control = ALS - Control, levels = design)
fit2  <- eBayes(contrasts.fit(fit, contr))

deg_all      <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
deg_all$gene <- rownames(deg_all)
deg_sig      <- deg_all %>% filter(abs(logFC) >= 0.25 & adj.P.Val < 0.05)

message(sprintf("DEGs identified: %d total | Up: %d | Down: %d",
                nrow(deg_sig),
                sum(deg_sig$logFC > 0),
                sum(deg_sig$logFC < 0)))

# Sensitivity note — how many genes survive a stricter threshold?
n_strict <- sum(abs(deg_all$logFC) >= 0.5 & deg_all$adj.P.Val < 0.05)
message(sprintf("Sensitivity check: %d DEGs survive stricter |logFC| >= 0.5 threshold.",
                n_strict))

write.csv(deg_all, "Processed_Data/Step2_DEG_all.csv", row.names = FALSE)
write.csv(deg_sig, "Processed_Data/Step2_DEG_sig.csv", row.names = FALSE)

# ── Figure 2B: Volcano Plot ────────────────────────────────────────────────────
vol <- deg_all %>%
  mutate(
    negLogFDR = -log10(adj.P.Val),
    Direction = case_when(
      adj.P.Val < 0.05 & logFC >=  0.25 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC <= -0.25 ~ "Downregulated",
      TRUE                               ~ "Not Significant"
    )
  )

top_labeled <- vol %>%
  filter(Direction != "Not Significant") %>%
  arrange(adj.P.Val) %>%
  head(20)

fig2B <- ggplot(vol, aes(x = logFC, y = negLogFDR, color = Direction)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c(
    "Downregulated"   = "#008B8B",
    "Upregulated"     = "#DC143C",
    "Not Significant" = "#E0E0E0"
  )) +
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_text_repel(
    data          = top_labeled,
    aes(label     = gene),
    size          = 4.2,
    fontface      = "italic",
    box.padding   = 0.5,
    max.overlaps  = Inf,
    color         = "black",
    segment.color = "grey60",
    segment.size  = 0.3
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "top",
    legend.title    = element_blank(),
    plot.title      = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40", size = 11)
  ) +
  labs(
    title    = "Differentially Expressed Genes",
    subtitle = sprintf("GSE112676 | %d DEGs (|logFC| >= 0.25, FDR < 0.05)",
                       nrow(deg_sig)),
    x = bquote(log[2]("Fold Change")),
    y = bquote(-log[10]("FDR"))
  )

ggsave("Manuscript_Figures/Figure2B_Volcano.pdf",
       plot = fig2B, width = 7.5, height = 7.5)
ggsave("Manuscript_Figures/Figure2B_Volcano.png",
       plot = fig2B, width = 7.5, height = 7.5, dpi = 600)
message("Figure 2B saved.")

# ── Figure 2C: Top DEG Heatmap ─────────────────────────────────────────────────
n_top     <- min(100, nrow(deg_sig))
top_genes <- deg_sig %>% arrange(adj.P.Val) %>% head(n_top) %>% pull(gene)

hm_scaled <- t(scale(t(expr_discovery[top_genes, ])))

ha_col <- HeatmapAnnotation(
  Diagnosis = group_factor,
  col = list(Diagnosis = c("ALS" = "#DC143C", "Control" = "#008B8B")),
  annotation_name_gp = gpar(fontface = "bold", fontsize = 11)
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#6495ED", "#FFFACD", "#FF6347"))

ht <- Heatmap(
  hm_scaled,
  name              = "Z-score",
  top_annotation    = ha_col,
  col               = col_fun,
  show_row_names    = FALSE,
  show_column_names = FALSE,
  cluster_columns   = TRUE,
  cluster_rows      = TRUE,
  show_row_dend     = TRUE,
  show_column_dend  = FALSE,
  column_title      = sprintf("Top %d DEGs — GSE112676 Discovery Cohort", n_top),
  column_title_gp   = gpar(fontsize = 14, fontface = "bold"),
  use_raster        = TRUE
)

pdf("Manuscript_Figures/Figure2C_Heatmap.pdf", width = 8.5, height = 7.5)
draw(ht, merge_legend = TRUE)
dev.off()
png("Manuscript_Figures/Figure2C_Heatmap.png",
    width = 8.5, height = 7.5, units = "in", res = 600)
draw(ht, merge_legend = TRUE)
dev.off()
message(sprintf("Figure 2C saved (top %d DEGs).", n_top))

# ── Figure 2D: MA Plot ─────────────────────────────────────────────────────────
# MA plot shows fold change as a function of mean expression — this reveals
# whether differential expression is artifactually enriched at either extreme
# of the expression range (a common problem in poorly normalized data).
fig2D <- ggplot(vol, aes(x = AveExpr, y = logFC, color = Direction)) +
  geom_point(alpha = 0.6, size = 1.4) +
  scale_color_manual(values = c(
    "Downregulated"   = "#008B8B",
    "Upregulated"     = "#DC143C",
    "Not Significant" = "#E0E0E0"
  ), name = "Direction") +
  geom_hline(yintercept = c(-0.25, 0.25),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_smooth(data = subset(vol, Direction != "Not Significant"),
              method = "loess", se = FALSE,
              color = "navy", linewidth = 0.8, linetype = "dotted") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "MA Plot — Mean Expression vs Fold Change (GSE112676)",
    x     = "Average log2 Expression",
    y     = bquote(log[2]("Fold Change"))
  )

ggsave("Manuscript_Figures/Figure2D_MA_Plot.pdf",
       plot = fig2D, width = 6.5, height = 6)
ggsave("Manuscript_Figures/Figure2D_MA_Plot.png",
       plot = fig2D, width = 6.5, height = 6, dpi = 600)
message("Figure 2D saved.")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
message("================================================================")
message("Script 01 complete.")
message(sprintf("  Dataset            : GSE112676 (discovery only)"))
message(sprintf("  Samples            : %d (ALS: %d | CON: %d)",
                ncol(expr_discovery),
                sum(group_factor == "ALS"),
                sum(group_factor == "Control")))
message(sprintf("  Outliers removed   : %d", length(outliers)))
message(sprintf("  Normalized         : %s (CV = %.4f)",
                ifelse(cv_medians > 0.02, "YES", "NO — pre-normalized"),
                cv_medians))
message(sprintf("  Genes in matrix    : %d", nrow(expr_discovery)))
message(sprintf("  Significant DEGs   : %d (Up: %d | Down: %d)",
                nrow(deg_sig),
                sum(deg_sig$logFC > 0),
                sum(deg_sig$logFC < 0)))
message("  GSE112680          : NOT LOADED — reserved for validation")
message("================================================================")
