# ==============================================================================
# 05_Immune_Analysis.R
#
# Three-layer immune deconvolution of GSE112676 discovery cohort.
#   Layer 1A: EPIC (primary — quantitative cell fractions)
#   Layer 1B: MCP-counter (orthogonal validation)
#   Layer 2:  ssGSEA (complementary enrichment scoring)
#   Layer 3:  ALS immune subtype stratification
#
# Inputs:
#   Processed_Data/Step1_Discovery_Data.rds
#   Processed_Data/Step5_PerGene_AUCs.csv
#
# Outputs:
#   Manuscript_Figures/Immune_Analysis/Figure6D_EPIC_Deconvolution.{pdf,png}
#   Manuscript_Figures/Immune_Analysis/Figure6E_ssGSEA_Boxplots.{pdf,png}
#   Manuscript_Figures/Immune_Analysis/Figure6F_Correlation_Heatmap.{pdf,png}
#   Manuscript_Figures/Immune_Analysis/Figure6G_ALS_Immune_Subtypes.{pdf,png}
#   Manuscript_Figures/Immune_Analysis/FigureSupp_MCPcounter.{pdf,png}
#   Manuscript_Figures/Immune_Analysis/FigureSupp_MultiMethod_Consensus.png
#   Manuscript_Figures/Immune_Analysis/FigureSupp_Subtype_Expression.{pdf,png}
#   Processed_Data/Supp_EPIC_Cell_Fractions.csv
#   Processed_Data/Supp_MCPcounter_Scores.csv
#   Processed_Data/ALS_Immune_Subtypes.csv
#   Processed_Data/Step6_Immune_Results.rds
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(EPIC)
library(MCPcounter)
library(GSVA)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(stringr)

dir.create("Manuscript_Figures/Immune_Analysis", recursive = TRUE,
           showWarnings = FALSE)
dir.create("Processed_Data", showWarnings = FALSE)

# ==============================================================================
# STEP 1: Load Data
# ==============================================================================
message("Loading expression data and signature genes...")

expr_data    <- readRDS("Processed_Data/Step1_Discovery_Data.rds")
expr_mat     <- expr_data$expr_discovery
group_factor <- expr_data$group_factor

if (is.null(names(group_factor))) {
  names(group_factor) <- colnames(expr_mat)
}

diagnosis <- factor(
  ifelse(group_factor == "ALS", "ALS", "Control"),
  levels = c("Control", "ALS")
)

# Load top 6 BPS genes — use correct file from this pipeline
bps_file <- "Processed_Data/Step5_PerGene_AUCs.csv"
if (!file.exists(bps_file)) {
  stop(sprintf("BPS file not found: %s\nRun Script 04 first.", bps_file))
}

bps_df     <- read.csv(bps_file)
core_genes <- head(bps_df$Gene[order(bps_df$BPS, decreasing = TRUE)], 6)

message(sprintf("Matrix: %d genes x %d samples | ALS: %d | CON: %d",
                nrow(expr_mat), ncol(expr_mat),
                sum(diagnosis == "ALS"),
                sum(diagnosis == "Control")))
message(sprintf("Core genes for correlation: %s",
                paste(core_genes, collapse = ", ")))

# EPIC and MCP-counter require linear (non-log) scale
# Data is already log2 — back-transform
if (max(expr_mat, na.rm = TRUE) < 50) {
  message("Log-scale detected — back-transforming to linear for EPIC/MCP-counter...")
  expr_linear <- 2^expr_mat - 1
} else {
  expr_linear <- expr_mat
}
expr_linear[expr_linear < 0] <- 0

# ==============================================================================
# STEP 2: EPIC Deconvolution (Primary Method)
# Racle et al., eLife 2017
# ==============================================================================
message("Running EPIC immune deconvolution (primary method)...")

epic_result <- tryCatch({
  EPIC::EPIC(bulk = expr_linear)
}, error = function(e) {
  message(sprintf("EPIC error: %s", e$message))
  NULL
})

if (is.null(epic_result)) {
  stop("EPIC failed. Verify expression matrix contains sufficient gene coverage.")
}

epic_fracs           <- as.data.frame(epic_result$cellFractions)
epic_fracs$Sample    <- rownames(epic_fracs)
epic_fracs$Diagnosis <- diagnosis[match(epic_fracs$Sample, colnames(expr_mat))]

cell_cols_epic <- colnames(epic_fracs)[
  !colnames(epic_fracs) %in% c("Sample", "Diagnosis")
]

message(sprintf("EPIC: %d cell types | %d samples",
                length(cell_cols_epic), nrow(epic_fracs)))
message(sprintf("Cell types: %s", paste(cell_cols_epic, collapse = ", ")))

write.csv(epic_fracs,
          "Processed_Data/Supp_EPIC_Cell_Fractions.csv",
          row.names = FALSE)

# Keep cell types with mean fraction > 0.5%
epic_means <- colMeans(epic_fracs[, cell_cols_epic], na.rm = TRUE)
keep_epic  <- names(epic_means[epic_means > 0.005])

message(sprintf("Retaining %d / %d cell types (mean fraction > 0.5%%)",
                length(keep_epic), length(cell_cols_epic)))

epic_long <- epic_fracs %>%
  filter(!is.na(Diagnosis)) %>%
  select(Sample, Diagnosis, all_of(keep_epic)) %>%
  pivot_longer(
    cols      = -c(Sample, Diagnosis),
    names_to  = "Cell_Type",
    values_to = "Fraction"
  )

# ── Figure 6D: EPIC Violin Plots ───────────────────────────────────────────────
p_epic <- ggplot(
  epic_long,
  aes(x = reorder(Cell_Type, Fraction, FUN = median),
      y = Fraction, fill = Diagnosis)
) +
  geom_violin(trim = FALSE, alpha = 0.7,
              position = position_dodge(0.85)) +
  geom_boxplot(width = 0.12, fill = "white", color = "black",
               outlier.shape = NA,
               position = position_dodge(0.85)) +
  scale_fill_manual(
    values = c("ALS" = "#BC3C29FF", "Control" = "#0072B5FF")
  ) +
  stat_compare_means(
    aes(group = Diagnosis),
    label   = "p.signif",
    method  = "wilcox.test",
    hide.ns = TRUE,
    size    = 4
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1,
                                    color = "black", face = "bold"),
    axis.text.y     = element_text(color = "black"),
    legend.position = "top",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "6D. Immune Cell Composition (EPIC Deconvolution)",
    subtitle = "Racle et al., eLife 2017 | Quantitative cell fractions",
    x        = "Immune Cell Type",
    y        = "Estimated Cell Fraction"
  )

ggsave("Manuscript_Figures/Immune_Analysis/Figure6D_EPIC_Deconvolution.pdf",
       plot = p_epic, width = 13, height = 7)
ggsave("Manuscript_Figures/Immune_Analysis/Figure6D_EPIC_Deconvolution.png",
       plot = p_epic, width = 13, height = 7, dpi = 600)
message("Figure 6D saved.")

# ==============================================================================
# STEP 3: MCP-counter (Orthogonal Validation)
# Becht et al., Genome Biology 2016
# ==============================================================================
message("Running MCP-counter (orthogonal validation)...")

mcp_result <- tryCatch({
  MCPcounter::MCPcounter.estimate(
    expression   = expr_linear,
    featuresType = "HUGO_symbols"
  )
}, error = function(e) {
  message(sprintf("MCP-counter error: %s", e$message))
  NULL
})

mcp_df        <- NULL
cell_cols_mcp <- NULL

if (!is.null(mcp_result)) {
  mcp_df           <- as.data.frame(t(mcp_result))
  mcp_df$Sample    <- rownames(mcp_df)
  mcp_df$Diagnosis <- diagnosis[match(mcp_df$Sample, colnames(expr_mat))]
  cell_cols_mcp    <- colnames(mcp_df)[
    !colnames(mcp_df) %in% c("Sample", "Diagnosis")
  ]

  write.csv(mcp_df,
            "Processed_Data/Supp_MCPcounter_Scores.csv",
            row.names = FALSE)
  message(sprintf("MCP-counter: %d cell types estimated.",
                  length(cell_cols_mcp)))

  mcp_long <- mcp_df %>%
    filter(!is.na(Diagnosis)) %>%
    select(Sample, Diagnosis, all_of(cell_cols_mcp)) %>%
    pivot_longer(
      cols      = -c(Sample, Diagnosis),
      names_to  = "Cell_Type",
      values_to = "Score"
    )

  # Supplementary Figure: MCP-counter
  p_mcp <- ggplot(
    mcp_long,
    aes(x = reorder(Cell_Type, Score, FUN = median),
        y = Score, fill = Diagnosis)
  ) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    scale_fill_manual(
      values = c("ALS" = "#BC3C29FF", "Control" = "#0072B5FF")
    ) +
    stat_compare_means(
      aes(group = Diagnosis),
      label   = "p.signif",
      method  = "wilcox.test",
      hide.ns = TRUE,
      size    = 4
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1,
                                      color = "black", face = "bold"),
      legend.position = "top",
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(hjust = 0.5, color = "grey40")
    ) +
    labs(
      title    = "Supp. Immune Cell Abundance (MCP-counter)",
      subtitle = "Becht et al., Genome Biology 2016 | Orthogonal validation",
      x        = "Immune Cell Type",
      y        = "MCP-counter Score"
    )

  ggsave("Manuscript_Figures/Immune_Analysis/FigureSupp_MCPcounter.pdf",
         plot = p_mcp, width = 11, height = 6)
  ggsave("Manuscript_Figures/Immune_Analysis/FigureSupp_MCPcounter.png",
         plot = p_mcp, width = 11, height = 6, dpi = 600)
  message("Supplementary Figure (MCP-counter) saved.")

  # Multi-method consensus log2FC heatmap
  common_cells <- intersect(keep_epic, cell_cols_mcp)
  if (length(common_cells) >= 3) {
    diag_vec <- as.character(diagnosis)

    calc_log2fc <- function(wide_df, diag_v, cells) {
      mat      <- as.matrix(wide_df[, cells, drop = FALSE])
      als_mn   <- colMeans(mat[diag_v == "ALS",     , drop = FALSE],
                           na.rm = TRUE)
      ctrl_mn  <- colMeans(mat[diag_v == "Control", , drop = FALSE],
                           na.rm = TRUE)
      log2((als_mn + 0.001) / (ctrl_mn + 0.001))
    }

    epic_fc <- calc_log2fc(epic_fracs, diag_vec, common_cells)
    mcp_fc  <- calc_log2fc(mcp_df,    diag_vec, common_cells)

    consensus_mat <- cbind(EPIC = epic_fc, `MCP-counter` = mcp_fc)
    col_consensus <- colorRamp2(
      c(-0.5, 0, 0.5),
      c("#0072B5FF", "white", "#BC3C29FF")
    )

    ht_consensus <- Heatmap(
      consensus_mat,
      name              = "log2FC\n(ALS/Ctrl)",
      col               = col_consensus,
      cluster_rows      = TRUE,
      cluster_columns   = FALSE,
      show_row_names    = TRUE,
      show_column_names = TRUE,
      row_names_gp      = gpar(fontsize = 10),
      column_names_gp   = gpar(fontsize = 12, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(sprintf("%.2f", consensus_mat[i, j]),
                        x, y, gp = gpar(fontsize = 9))
      },
      column_title    = "Supp. Multi-Method Immune Consensus (log2FC: ALS vs Control)",
      column_title_gp = gpar(fontsize = 12, fontface = "bold")
    )

    png("Manuscript_Figures/Immune_Analysis/FigureSupp_MultiMethod_Consensus.png",
        width = 6, height = 7, units = "in", res = 600)
    draw(ht_consensus)
    dev.off()
    message("Multi-method consensus heatmap saved.")
  }
} else {
  message("MCP-counter failed — continuing without orthogonal validation.")
}

# ==============================================================================
# STEP 4: ssGSEA Immune Activity Scoring
# ==============================================================================
message("Running ssGSEA immune activity scoring...")

immune_signatures <- list(
  `B cells`         = c("CD19",     "CD79A",    "MS4A1",    "BANK1",    "CD79B"),
  `CD4+ T cells`    = c("CD4",      "CD40LG",   "IL7R",     "TCF7",     "CD27",   "MAL"),
  `CD8+ T cells`    = c("CD8A",     "CD8B",     "GZMA",     "GZMB",     "PRF1"),
  `NK cells`        = c("NCAM1",    "KLRB1",    "KLRD1",    "NKG7",     "KLRC1"),
  `Macrophages`     = c("CD68",     "CD163",    "MRC1",     "MSR1",     "CHIT1"),
  `Neutrophils`     = c("CEACAM8",  "FCGR3B",   "CSF3R",    "S100A12",  "CXCR2"),
  `Dendritic cells` = c("HLA-DPA1", "HLA-DPB1", "CD1C",     "HLA-DRA",  "HLA-DQA1"),
  `Monocytes`       = c("CD14",     "FCGR3A",   "CD33",     "ITGAM",    "CCR2"),
  `Tregs`           = c("FOXP3",    "IL2RA",    "CTLA4",    "CCR8",     "IKZF2"),
  `Mast cells`      = c("TPSAB1",   "TPSB2",    "CPA3",     "HDC",      "MS4A2")
)

message("Marker gene recovery per cell type:")
for (sig_name in names(immune_signatures)) {
  found <- intersect(immune_signatures[[sig_name]], rownames(expr_mat))
  message(sprintf("  %-18s : %d / %d markers found",
                  sig_name,
                  length(found),
                  length(immune_signatures[[sig_name]])))
}

# GSVA version-adaptive call
immune_scores <- tryCatch({
  gsva(ssgseaParam(expr_mat, immune_signatures, minSize = 2),
       verbose = FALSE)
}, error = function(e1) {
  message(sprintf("ssgseaParam() failed (%s) -- trying legacy...", e1$message))
  tryCatch(
    gsva(expr_mat, immune_signatures, method = "ssgsea", verbose = FALSE),
    error = function(e2) {
      message(sprintf("Legacy failed (%s) -- gsvaParam fallback...", e2$message))
      gsva(gsvaParam(expr_mat, immune_signatures, minSize = 2),
           verbose = FALSE)
    }
  )
})

immune_df <- as.data.frame(t(immune_scores))
immune_df$Diagnosis <- diagnosis

immune_long <- immune_df %>%
  pivot_longer(
    cols      = -Diagnosis,
    names_to  = "Cell_Type",
    values_to = "ssGSEA_Score"
  )

# ── Figure 6E: ssGSEA Boxplots ─────────────────────────────────────────────────
p_ssgsea <- ggplot(
  immune_long,
  aes(x = reorder(Cell_Type, ssGSEA_Score, FUN = median),
      y = ssGSEA_Score, fill = Diagnosis)
) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(
    values = c("ALS" = "#BC3C29FF", "Control" = "#0072B5FF")
  ) +
  stat_compare_means(
    aes(group = Diagnosis),
    label   = "p.signif",
    method  = "wilcox.test",
    hide.ns = TRUE,
    size    = 4
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1,
                                    color = "black", face = "bold"),
    axis.text.y     = element_text(color = "black"),
    legend.position = "top",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "6E. Immune Cell Activity Profile (ssGSEA)",
    subtitle = "Curated marker gene sets | 10 immune cell types",
    x        = "Immune Cell Type",
    y        = "ssGSEA Enrichment Score"
  )

ggsave("Manuscript_Figures/Immune_Analysis/Figure6E_ssGSEA_Boxplots.pdf",
       plot = p_ssgsea, width = 11, height = 6)
ggsave("Manuscript_Figures/Immune_Analysis/Figure6E_ssGSEA_Boxplots.png",
       plot = p_ssgsea, width = 11, height = 6, dpi = 600)
message("Figure 6E saved.")

# ── Figure 6F: Spearman Correlation Heatmap ────────────────────────────────────
available_core <- intersect(core_genes, rownames(expr_mat))

if (length(available_core) < length(core_genes)) {
  message(sprintf("Missing core genes from expression matrix: %s",
                  paste(setdiff(core_genes, rownames(expr_mat)),
                        collapse = ", ")))
}

shared_samples <- intersect(colnames(expr_mat), epic_fracs$Sample)
core_for_cor   <- t(expr_mat[available_core, shared_samples, drop = FALSE])
epic_for_cor   <- epic_fracs[
  match(shared_samples, epic_fracs$Sample),
  keep_epic, drop = FALSE
]

cor_mat  <- cor(x = core_for_cor, y = epic_for_cor, method = "spearman")
col_heat <- colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF"))(100)

heatmap_args <- list(
  mat             = cor_mat,
  color           = col_heat,
  display_numbers = round(cor_mat, 2),
  fontsize_number = 9,
  fontsize_row    = 11,
  fontsize_col    = 10,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  angle_col       = "45",
  border_color    = "grey85",
  main            = "6F. Spearman Correlation: ALS Signature vs Immune Cells (EPIC)"
)

pdf("Manuscript_Figures/Immune_Analysis/Figure6F_Correlation_Heatmap.pdf",
    width = 10, height = 6)
do.call(pheatmap, heatmap_args)
dev.off()

png("Manuscript_Figures/Immune_Analysis/Figure6F_Correlation_Heatmap.png",
    width = 10, height = 6, units = "in", res = 600)
do.call(pheatmap, heatmap_args)
dev.off()
message("Figure 6F saved.")

# ==============================================================================
# STEP 5: ALS Immune Subtype Stratification
# ==============================================================================
message("Running ALS immune subtype stratification...")

als_samples <- epic_fracs$Sample[
  epic_fracs$Diagnosis == "ALS" & !is.na(epic_fracs$Diagnosis)
]

als_frac_mat           <- as.matrix(
  epic_fracs[match(als_samples, epic_fracs$Sample), keep_epic, drop = FALSE]
)
rownames(als_frac_mat) <- als_samples

# Silhouette-optimal k (2 to 5)
sil_scores <- sapply(2:5, function(k) {
  set.seed(1122)
  km  <- kmeans(als_frac_mat, centers = k, nstart = 25, iter.max = 100)
  sil <- silhouette(km$cluster, dist(als_frac_mat))
  mean(sil[, 3])
})
names(sil_scores) <- paste0("k=", 2:5)

message("Silhouette scores:")
print(round(sil_scores, 3))

optimal_k <- which.max(sil_scores) + 1
message(sprintf("Optimal k: %d immune subtypes", optimal_k))

set.seed(1122)
km_final       <- kmeans(als_frac_mat, centers = optimal_k,
                          nstart = 25, iter.max = 100)
subtype_labels <- paste0("Immune-", LETTERS[km_final$cluster])
names(subtype_labels) <- als_samples

message("ALS immune subtype distribution:")
print(table(subtype_labels))

write.csv(
  data.frame(Sample  = als_samples,
             Subtype = subtype_labels,
             stringsAsFactors = FALSE),
  "Processed_Data/ALS_Immune_Subtypes.csv",
  row.names = FALSE
)

# ── Figure 6G: Stratification Heatmap ─────────────────────────────────────────
subtype_colors <- setNames(
  RColorBrewer::brewer.pal(max(3, optimal_k), "Set2")[1:optimal_k],
  paste0("Immune-", LETTERS[1:optimal_k])
)

ha_top <- HeatmapAnnotation(
  `Immune Subtype`     = subtype_labels[als_samples],
  col                  = list(`Immune Subtype` = subtype_colors),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 11, fontface = "bold")
)

col_strat <- colorRamp2(
  c(0, 0.10, 0.35),
  c("#4DBBD5FF", "#FFFACD", "#E64B35FF")
)

ht_strat <- Heatmap(
  t(als_frac_mat),
  name              = "Cell\nFraction",
  col               = col_strat,
  top_annotation    = ha_top,
  show_column_names = FALSE,
  show_row_names    = TRUE,
  cluster_columns   = TRUE,
  cluster_rows      = TRUE,
  row_names_gp      = gpar(fontsize = 10),
  column_split      = subtype_labels[als_samples],
  column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
  column_gap        = unit(3, "mm"),
  use_raster        = TRUE
)

png("Manuscript_Figures/Immune_Analysis/Figure6G_ALS_Immune_Subtypes.png",
    width = 12, height = 7, units = "in", res = 600)
draw(ht_strat,
     column_title    = "6G. ALS Patient Immune Subtype Stratification (EPIC)",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

pdf("Manuscript_Figures/Immune_Analysis/Figure6G_ALS_Immune_Subtypes.pdf",
    width = 12, height = 7)
draw(ht_strat,
     column_title    = "6G. ALS Patient Immune Subtype Stratification (EPIC)",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
message("Figure 6G saved.")

# Supplementary: biomarker expression per immune subtype
als_expr_t <- as.data.frame(
  t(expr_mat[available_core,
              intersect(als_samples, colnames(expr_mat)),
              drop = FALSE])
)
als_expr_t$Subtype <- factor(subtype_labels[rownames(als_expr_t)])

als_expr_long <- als_expr_t %>%
  pivot_longer(cols = -Subtype,
               names_to  = "Gene",
               values_to = "Expression")

p_subtype <- ggplot(
  als_expr_long,
  aes(x = Subtype, y = Expression, fill = Subtype)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, fill = "white",
               outlier.shape = NA, color = "black") +
  scale_fill_manual(values = subtype_colors) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  stat_compare_means(method = "kruskal.test",
                     label  = "p.format",
                     size   = 3.5) +
  theme_classic(base_size = 12) +
  theme(
    strip.text      = element_text(face = "bold.italic", size = 11),
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Supp. ALS Signature Gene Expression Across Immune Subtypes",
    x     = "Immune Subtype",
    y     = "Normalised Expression"
  )

ggsave("Manuscript_Figures/Immune_Analysis/FigureSupp_Subtype_Expression.pdf",
       plot = p_subtype, width = 12, height = 7)
ggsave("Manuscript_Figures/Immune_Analysis/FigureSupp_Subtype_Expression.png",
       plot = p_subtype, width = 12, height = 7, dpi = 600)
message("Supplementary subtype expression figure saved.")

# ==============================================================================
# STEP 6: Save All Results
# ==============================================================================
saveRDS(
  list(
    epic_fracs     = epic_fracs,
    keep_epic      = keep_epic,
    mcp_df         = mcp_df,
    immune_df      = immune_df,
    subtype_labels = subtype_labels,
    optimal_k      = optimal_k,
    core_genes     = core_genes
  ),
  "Processed_Data/Step6_Immune_Results.rds"
)

message("================================================================")
message("Script 05 complete.")
message(sprintf("  EPIC cell types          : %d (retained: %d)",
                length(cell_cols_epic), length(keep_epic)))
message(sprintf("  MCP-counter              : %s",
                ifelse(!is.null(mcp_df), "completed", "failed - skipped")))
message(sprintf("  ALS immune subtypes      : %d", optimal_k))
message(sprintf("  Core genes used          : %s",
                paste(available_core, collapse = ", ")))
message("  Output: Manuscript_Figures/Immune_Analysis/")
message("================================================================")
