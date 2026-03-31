# ==============================================================================
# 02_WGCNA_Analysis.R
#
# Weighted Gene Co-expression Network Analysis on the GSE112676 discovery cohort.
#
# Key methodological decisions:
#   - Top 5,000 genes by MAD (robust to outliers; superior to variance)
#   - Soft-thresholding power: lowest power achieving R2 >= 0.80 with
#     acceptable mean connectivity (>= 5). NOT pickSoftThreshold's default
#     estimate, which optimizes R2 at the cost of connectivity and produces
#     over-sparse networks with too few modules.
#   - deepSplit = 4 for fine-grained module detection
#   - signedKME() for kME computation (correct for signed hybrid networks)
#   - Hub candidates extracted from ALL significant modules
#   - traitData aligned by sample name matching, not positional assumption
#
# Inputs:
#   Processed_Data/Step1_Discovery_Data.rds
#
# Outputs:
#   Processed_Data/Step3_WGCNA_Data.rds
#   Processed_Data/Candidate_Hubs_for_ML.csv
#   Processed_Data/Step3_Top200_for_STRING.csv
#   Manuscript_Figures/WGCNA/Figure3A_ScaleFree_Fit.{pdf,png}
#   Manuscript_Figures/WGCNA/Figure3B_Mean_Connectivity.{pdf,png}
#   Manuscript_Figures/WGCNA/Figure3C_Dendrogram.{pdf,png}
#   Manuscript_Figures/WGCNA/Figure3D_Module_Trait.{pdf,png}
#   Manuscript_Figures/WGCNA/Figure3E_ME_Correlation.{pdf,png}
#   Manuscript_Figures/WGCNA/Figure3F_GS_vs_MM_{module}.{pdf,png}
#   Manuscript_Figures/WGCNA/Figure3G_ME_Violin.{pdf,png}
#   Manuscript_Figures/WGCNA/Figure3H_Gene_Count_Bar.{pdf,png}
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(WGCNA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)
library(ggpubr)

allowWGCNAThreads()

dir.create("Manuscript_Figures/WGCNA", recursive = TRUE, showWarnings = FALSE)
dir.create("Processed_Data",           showWarnings = FALSE)

# ==============================================================================
# STEP 1: Load Discovery Data
# ==============================================================================
message("Loading discovery cohort data...")

step1        <- readRDS("Processed_Data/Step1_Discovery_Data.rds")
expr_mat     <- step1$expr_discovery     # genes x samples
group_factor <- step1$group_factor       # factor, levels = Control/ALS

if (is.null(names(group_factor))) {
  names(group_factor) <- colnames(expr_mat)
}

message(sprintf("Discovery matrix: %d genes x %d samples | ALS: %d | CON: %d",
                nrow(expr_mat), ncol(expr_mat),
                sum(group_factor == "ALS"),
                sum(group_factor == "Control")))

# ==============================================================================
# STEP 2: Select Top 5,000 Genes by MAD
# ==============================================================================
gene_mad  <- apply(expr_mat, 1, mad, na.rm = TRUE)
n_genes   <- min(5000, length(gene_mad))
top_genes <- names(sort(gene_mad, decreasing = TRUE))[1:n_genes]
datExpr   <- t(expr_mat[top_genes, ])   # WGCNA format: samples x genes

message(sprintf("Top %d genes by MAD selected for network construction.",
                ncol(datExpr)))

# ==============================================================================
# STEP 3: WGCNA Quality Check
# ==============================================================================
message("Running WGCNA goodSamplesGenes check...")

gsg <- goodSamplesGenes(datExpr, verbose = 0)

if (!gsg$allOK) {
  n_g <- sum(!gsg$goodGenes)
  n_s <- sum(!gsg$goodSamples)
  message(sprintf("  Removing %d genes and %d samples failing WGCNA QC.", n_g, n_s))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
} else {
  message("  All genes and samples passed WGCNA QC.")
}

message(sprintf("WGCNA input: %d samples x %d genes",
                nrow(datExpr), ncol(datExpr)))

# ==============================================================================
# STEP 4: Build Trait Matrix
# ==============================================================================
sample_ids    <- rownames(datExpr)
expr_ids      <- names(group_factor)
idx           <- match(sample_ids, expr_ids)

if (any(is.na(idx))) {
  stop(sprintf(
    "%d datExpr samples have no match in group_factor. Check preprocessing.",
    sum(is.na(idx))
  ))
}

group_aligned <- as.character(group_factor)[idx]

traitData <- data.frame(
  ALS       = as.numeric(group_aligned == "ALS"),
  row.names = sample_ids
)

message(sprintf("Trait alignment: ALS = %d | CON = %d | NA = %d",
                sum(traitData$ALS == 1),
                sum(traitData$ALS == 0),
                sum(is.na(traitData$ALS))))

if (any(is.na(traitData$ALS))) {
  stop("NA values in traitData after alignment. Investigate sample ID mismatch.")
}

# ==============================================================================
# STEP 5: Soft-Thresholding Power Selection
# ==============================================================================
message("Selecting soft-thresholding power...")

powers <- c(1:20, seq(22, 30, by = 2))
sft    <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed hybrid",
  verbose     = 0
)

sft_df <- data.frame(
  Power = sft$fitIndices$Power,
  R2    = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
  MeanK = sft$fitIndices$mean.k.
)

message("Scale-free topology summary:")
print(sft_df[sft_df$R2 > 0.5, ])

# ==============================================================================
# POWER SELECTION -- computed here so p3A can use chosen_power for the highlight
# ==============================================================================
good_r2   <- sft_df$Power[sft_df$R2 >= 0.80]
good_conn <- sft_df$Power[sft_df$MeanK >= 5]
both_ok   <- intersect(good_r2, good_conn)

if (length(both_ok) > 0) {
  chosen_power <- min(both_ok)
  message(sprintf("Selected power = %d (R2 = %.3f, mean.k = %.2f) - both criteria met.",
                  chosen_power,
                  sft_df$R2[sft_df$Power == chosen_power],
                  sft_df$MeanK[sft_df$Power == chosen_power]))
} else if (length(good_r2) > 0) {
  chosen_power <- min(good_r2)
  message("WARNING: No power satisfies both R2 >= 0.80 AND mean.k >= 5.")
  message(sprintf("Using power = %d (R2 = %.3f, mean.k = %.2f) - R2 criterion only.",
                  chosen_power,
                  sft_df$R2[sft_df$Power == chosen_power],
                  sft_df$MeanK[sft_df$Power == chosen_power]))
} else {
  stop(paste(
    "No power achieves R2 >= 0.80.",
    "The dataset may not fit scale-free topology.",
    "Review Figure 3A before proceeding."
  ))
}

achieved_r2    <- sft_df$R2[sft_df$Power == chosen_power]
achieved_meank <- sft_df$MeanK[sft_df$Power == chosen_power]

message(sprintf("Final: power = %d | R2 = %.3f | mean.k = %.2f",
                chosen_power, achieved_r2, achieved_meank))

# ==============================================================================
# FIGURE 3A: Scale-Free Topology Fit
# AESTHETIC CHANGE: large filled circles with white bold power labels inside;
#                   red hollow circle highlights chosen_power;
#                   "Selected power = X" annotation below highlight;
#                   x-axis label updated to "Soft Threshold (Power)"
# ==============================================================================
p3A <- ggplot(sft_df, aes(x = Power, y = R2)) +
  geom_hline(yintercept = 0.80,
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_point(size = 8, color = "#BC3C29FF", fill = "#BC3C29FF", shape = 21) +
  geom_text(aes(label = Power), color = "white", size = 4, fontface = "bold") +
  geom_point(data = sft_df[sft_df$Power == chosen_power, ],
             aes(x = Power, y = R2),
             size = 8, shape = 21, color = "black", fill = NA, stroke = 1) +
  annotate("label",
           x         = chosen_power + 3,
           y         = min(sft_df$R2) + 0.08,
           label     = paste0("Selected power = ", chosen_power),
           size      = 3.8, color = "black", fontface = "bold",
           fill      = "white", label.size = 0.4) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "3A. Scale-Free Topology Fit",
       x     = "Soft Threshold (Power)",
       y     = "Scale-Free Topology R2")



# ==============================================================================
# FIGURE 3B: Mean Connectivity
# AESTHETIC CHANGE: large filled circles with white bold power labels inside;
#                   x-axis label updated to "Soft Threshold (Power)"
# ==============================================================================
p3B <- ggplot(sft_df, aes(x = Power, y = MeanK)) +
  geom_hline(yintercept = 5,
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_point(size = 8, color = "#0072B5FF", fill = "#0072B5FF", shape = 21) +
  geom_text(aes(label = Power), color = "white", size = 4, fontface = "bold") +
  coord_cartesian(ylim = c(min(sft_df$MeanK) - 50, max(sft_df$MeanK) + 50)) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "3B. Mean Connectivity",
       x     = "Soft Threshold (Power)",
       y     = "Mean Connectivity")

ggsave("Manuscript_Figures/WGCNA/Figure3A_ScaleFree_Fit.pdf",
       plot = p3A, width = 6, height = 5)
ggsave("Manuscript_Figures/WGCNA/Figure3A_ScaleFree_Fit.png",
       plot = p3A, width = 6, height = 5, dpi = 600)
ggsave("Manuscript_Figures/WGCNA/Figure3B_Mean_Connectivity.pdf",
       plot = p3B, width = 6, height = 5)
ggsave("Manuscript_Figures/WGCNA/Figure3B_Mean_Connectivity.png",
       plot = p3B, width = 6, height = 5, dpi = 600)
message("Figures 3A and 3B saved.")


# ==============================================================================
# STEP 6: Network Construction and Module Detection
# ==============================================================================
message(sprintf(
  "Building signed hybrid network (power = %d, deepSplit = 4)...",
  chosen_power))
message("This step may take 5-20 minutes depending on hardware...")

net <- blockwiseModules(
  datExpr,
  power             = chosen_power,
  networkType       = "signed hybrid",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  deepSplit         = 4,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 0
)

moduleColors <- net$colors
MEs          <- net$MEs
n_modules    <- length(unique(moduleColors)) - 1
n_grey       <- sum(moduleColors == "grey")
pct_grey     <- round(n_grey / length(moduleColors) * 100, 1)

message(sprintf("Network built: %d modules | grey (unassigned): %d genes (%.1f%%)",
                n_modules, n_grey, pct_grey))

if (pct_grey > 30) {
  message(sprintf("WARNING: %.1f%% of genes are unassigned (grey).", pct_grey))
  message("Consider running with lower power or minModuleSize = 20 if this persists.")
}

mod_table <- sort(table(moduleColors), decreasing = TRUE)
message("Module sizes (top 10):")
print(head(mod_table, 10))

# ==============================================================================
# FIGURE 3C: Dendrogram (UNCHANGED)
# ==============================================================================
pdf("Manuscript_Figures/WGCNA/Figure3C_Dendrogram.pdf", width = 14, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "3C. Gene Dendrogram and Module Colors"
)
dev.off()

png("Manuscript_Figures/WGCNA/Figure3C_Dendrogram.png",
    width = 14, height = 6, units = "in", res = 600)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "3C. Gene Dendrogram and Module Colors"
)
dev.off()
message("Figure 3C saved.")

# ==============================================================================
# STEP 7: Module-Trait Correlation
# ==============================================================================
message("Computing module-trait correlations...")

MEs_clean         <- removeGreyME(MEs, greyMEName = "MEgrey")
moduleTraitCor    <- cor(MEs_clean, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

message("Module-trait correlations:")
result_df <- data.frame(
  Module    = gsub("ME", "", rownames(moduleTraitCor)),
  r         = round(moduleTraitCor[, "ALS"], 3),
  p         = signif(moduleTraitPvalue[, "ALS"], 3),
  direction = ifelse(moduleTraitCor[, "ALS"] > 0, "POS", "NEG"),
  n_genes   = as.integer(table(moduleColors)[gsub("ME", "",
                                                    rownames(moduleTraitCor))])
)
result_df <- result_df[order(abs(result_df$r), decreasing = TRUE), ]
print(result_df)

# ==============================================================================
# FIGURE 3D: Module-Trait Heatmap (UNCHANGED)
# ==============================================================================
n_mods   <- nrow(moduleTraitCor)
text_mat <- matrix(
  paste0(round(moduleTraitCor, 2), "\n(",
         signif(moduleTraitPvalue, 1), ")"),
  nrow     = n_mods,
  dimnames = dimnames(moduleTraitCor)
)

for (ext in c("pdf", "png")) {
  if (ext == "pdf") {
    pdf("Manuscript_Figures/WGCNA/Figure3D_Module_Trait.pdf",
        width = 5, height = max(6, n_mods * 0.45))
  } else {
    png("Manuscript_Figures/WGCNA/Figure3D_Module_Trait.png",
        width = 5, height = max(6, n_mods * 0.45),
        units = "in", res = 600)
  }
  par(mar = c(4, 10, 3, 3))
  labeledHeatmap(
    Matrix        = moduleTraitCor,
    xLabels       = colnames(moduleTraitCor),
    yLabels       = rownames(moduleTraitCor),
    ySymbols      = rownames(moduleTraitCor),
    colorLabels   = FALSE,
    colors        = blueWhiteRed(50),
    textMatrix    = text_mat,
    setStdMargins = FALSE,
    cex.text      = 0.85,
    zlim          = c(-1, 1),
    main          = "3D. Module-Trait Relationships"
  )
  dev.off()
}
message("Figure 3D saved.")

# ==============================================================================
# FIGURE 3E: Module Eigengene Correlation Heatmap
# AESTHETIC CHANGE: each cell now shows r value AND p-value below it as
#                   "r\n(p)" with proper <1e-100 handling, matching reference
# ==============================================================================
modCor  <- cor(MEs_clean, use = "p")
modPval <- corPvalueStudent(modCor, nSamples = nrow(datExpr))

format_p <- function(p) {
  ifelse(p < 1e-100, "<1e-100", formatC(p, format = "e", digits = 1))
}
formatted_p_vals <- apply(modPval, c(1, 2), format_p)

textMatCor      <- paste(signif(modCor, 2), "\n(", formatted_p_vals, ")", sep = "")
dim(textMatCor) <- dim(modCor)

fig_dim <- max(8, n_mods * 0.7)
for (ext in c("pdf", "png")) {
  if (ext == "pdf") {
    pdf("Manuscript_Figures/WGCNA/Figure3E_ME_Correlation.pdf",
        width = fig_dim, height = fig_dim)
  } else {
    png("Manuscript_Figures/WGCNA/Figure3E_ME_Correlation.png",
        width = fig_dim, height = fig_dim, units = "in", res = 600)
  }
  par(mar = c(8, 8, 4, 3))
  labeledHeatmap(
    Matrix        = modCor,
    xLabels       = names(MEs_clean),
    yLabels       = names(MEs_clean),
    ySymbols      = names(MEs_clean),
    xSymbols      = names(MEs_clean),
    colorLabels   = FALSE,
    colors        = blueWhiteRed(50),
    textMatrix    = textMatCor,
    setStdMargins = FALSE,
    cex.text      = 0.6,
    zlim          = c(-1, 1),
    main          = "3E. Module Eigengene Correlation Heatmap",
    xLabelsAngle  = 45
  )
  dev.off()
}
message("Figure 3E saved.")

# ==============================================================================
# STEP 8: Identify Significant Modules (UNCHANGED)
# ==============================================================================
message("Identifying significantly correlated modules...")

valid_mask <- rownames(moduleTraitCor) != "MEgrey"
valid_cor  <- moduleTraitCor[valid_mask, "ALS", drop = FALSE]
valid_pval <- moduleTraitPvalue[valid_mask, "ALS", drop = FALSE]

sig_mask <- !is.na(valid_pval[, "ALS"]) & valid_pval[, "ALS"] < 0.05

if (sum(sig_mask) < 1) {
  message("No modules at p < 0.05. Relaxing to p < 0.10...")
  sig_mask <- !is.na(valid_pval[, "ALS"]) & valid_pval[, "ALS"] < 0.10
}

if (sum(sig_mask) < 1) {
  message("Still none. Using all non-NA modules...")
  sig_mask <- !is.na(valid_pval[, "ALS"])
}

sig_module_names <- rownames(valid_cor)[sig_mask]
sig_module_cols  <- gsub("ME", "", sig_module_names)

message(sprintf("Significant modules (%d total):", length(sig_module_cols)))
for (mn in sig_module_names) {
  mc <- gsub("ME", "", mn)
  message(sprintf("  %s: r = %.3f | p = %.2e | n = %d | %s",
                  mc, valid_cor[mn, "ALS"], valid_pval[mn, "ALS"],
                  sum(moduleColors == mc),
                  ifelse(valid_cor[mn, "ALS"] > 0, "POSITIVE", "NEGATIVE")))
}

pos_modules <- sig_module_names[
  !is.na(valid_cor[sig_module_names, "ALS"]) &
    valid_cor[sig_module_names, "ALS"] > 0
]

if (length(pos_modules) == 0) {
  pos_modules <- sig_module_names[
    which.max(abs(valid_cor[sig_module_names, "ALS"]))
  ]
  message("No positive modules found - using highest |r| module as primary.")
}

primary_module     <- pos_modules[which.max(abs(valid_cor[pos_modules, "ALS"]))]
primary_module_col <- gsub("ME", "", primary_module)

message(sprintf("Primary module (BPS reference): %s (r = %.3f, p = %.2e)",
                primary_module_col,
                valid_cor[primary_module, "ALS"],
                valid_pval[primary_module, "ALS"]))

# ==============================================================================
# STEP 9: Compute signedKME (UNCHANGED)
# ==============================================================================
message("Computing signed kME values...")

kME_all <- signedKME(datExpr, MEs_clean, outputColumnName = "kME")

message(sprintf("kME matrix: %d genes x %d columns",
                nrow(kME_all), ncol(kME_all)))

primary_kME_col <- paste0("kME", primary_module_col)

if (!primary_kME_col %in% colnames(kME_all)) {
  idx <- grep(primary_module_col, colnames(kME_all), ignore.case = TRUE)
  if (length(idx) > 0) {
    primary_kME_col <- colnames(kME_all)[idx[1]]
    message(sprintf("Primary kME column matched: %s", primary_kME_col))
  } else {
    primary_kME_col <- colnames(kME_all)[1]
    message(sprintf("WARNING: Primary kME not found - using %s as fallback.",
                    primary_kME_col))
  }
}

# ==============================================================================
# STEP 10: Gene Significance vs Module Membership -- Figure 3F
# AESTHETIC CHANGE: nudge_x = 1 - top_hub$MM anchors all labels to x = 1,
#                   direction = "y" stacks them in a clean vertical column,
#                   x-axis label updated to "Module Membership in X module"
# ==============================================================================
GS_ALS        <- as.numeric(cor(datExpr, traitData$ALS, use = "p"))
names(GS_ALS) <- colnames(datExpr)

neg_modules <- sig_module_names[
  !is.na(valid_cor[sig_module_names, "ALS"]) &
    valid_cor[sig_module_names, "ALS"] < 0
]

top_neg_col <- if (length(neg_modules) > 0) {
  gsub("ME", "", neg_modules[which.max(abs(valid_cor[neg_modules, "ALS"]))])
} else NULL

plot_modules <- unique(c(primary_module_col, top_neg_col))

for (mod_col in plot_modules) {
  kME_col <- paste0("kME", mod_col)
  if (!kME_col %in% colnames(kME_all)) {
    message(sprintf("Skipping Figure 3F for %s - kME column not found.", mod_col))
    next
  }

  gene_name_vec <- names(moduleColors)
  if (is.null(gene_name_vec)) gene_name_vec <- colnames(datExpr)

  inModule_genes <- gene_name_vec[moduleColors == mod_col]
  avail_genes    <- intersect(inModule_genes, rownames(kME_all))

  if (length(avail_genes) < 10) {
    message(sprintf("Skipping Figure 3F for %s - too few genes (%d).",
                    mod_col, length(avail_genes)))
    next
  }

  df_gsm <- data.frame(
    MM   = kME_all[avail_genes, kME_col],
    GS   = GS_ALS[avail_genes],
    Gene = avail_genes
  )

  top_hub <- df_gsm %>% arrange(desc(abs(MM))) %>% head(10)

  p_gsm <- ggplot(df_gsm, aes(x = MM, y = GS)) +
    geom_point(fill = mod_col, color = "black",
               alpha = 0.8, size = 2.5, shape = 21, stroke = 0.3) +
    geom_smooth(method = "lm", se = TRUE,
                color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_label_repel(
      data          = top_hub,
      aes(label     = Gene),
      size          = 4.5,
      fontface      = "bold.italic",
      color         = "black",
      fill          = alpha("white", 0.8),
      label.size    = 0,
      nudge_x       = 1 - top_hub$MM,
      direction     = "y",
      hjust         = 0,
      segment.color = "grey60",
      segment.size  = 0.4,
      max.overlaps  = Inf
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    theme_classic(base_size = 15) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(
      title = sprintf("3F. GS vs MM - %s Module", toupper(mod_col)),
      x     = paste("Module Membership in", mod_col, "module"),
      y     = "Gene Significance (ALS Trait)"
    )

  ggsave(sprintf("Manuscript_Figures/WGCNA/Figure3F_GS_vs_MM_%s.pdf", mod_col),
         plot = p_gsm, width = 8, height = 6.5)
  ggsave(sprintf("Manuscript_Figures/WGCNA/Figure3F_GS_vs_MM_%s.png", mod_col),
         plot = p_gsm, width = 8, height = 6.5, dpi = 600)
  message(sprintf("Figure 3F (%s) saved.", mod_col))
}

# ==============================================================================
# STEP 11: Module Eigengene Violin Plot -- Figure 3G
# AESTHETIC CHANGE: trimmed to top 4 modules by |r|; nrow = 2 forces a clean
#                   2x2 grid; fixed 10x10 dimensions for balanced panels
# ==============================================================================
top4_idx     <- order(abs(valid_cor[sig_module_names, "ALS"]),
                      decreasing = TRUE)[1:min(4, length(sig_module_names))]
top4_modules <- sig_module_names[top4_idx]

message(sprintf("Figure 3G: top %d modules by |r|: %s",
                length(top4_modules),
                paste(gsub("ME", "", top4_modules), collapse = ", ")))

me_df <- as.data.frame(MEs_clean[, top4_modules, drop = FALSE])
me_df$Diagnosis <- factor(
  ifelse(traitData$ALS == 1, "ALS", "Control"),
  levels = c("Control", "ALS")
)

me_long <- me_df %>%
  pivot_longer(cols      = -Diagnosis,
               names_to  = "Module",
               values_to = "Eigengene") %>%
  mutate(Module = factor(Module, levels = top4_modules))

p_violin <- ggplot(me_long,
                   aes(x = Diagnosis, y = Eigengene, fill = Diagnosis)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, fill = "white",
               color = "black", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  stat_compare_means(
    method      = "wilcox.test",
    label       = "p.format",
    comparisons = list(c("Control", "ALS")),
    size        = 3.5
  ) +
  scale_fill_manual(values = c("ALS" = "#BC3C29FF", "Control" = "#0072B5FF")) +
  facet_wrap(~ Module, scales = "free_y", nrow = 2,
             labeller = labeller(Module = setNames(
               sapply(top4_modules, function(m) {
                 mc    <- gsub("ME", "", m)
                 r_val <- valid_cor[m, "ALS"]
                 p_val <- valid_pval[m, "ALS"]
                 if (is.na(r_val)) return(mc)
                 sprintf("%s\nr=%.2f, p=%.1e", mc, r_val, p_val)
               }),
               top4_modules
             ))) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    legend.position  = "top",
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "3G. Module Eigengene Expression (ALS vs Control)",
    subtitle = sprintf("Top %d modules by |r|", length(top4_modules)),
    x        = "",
    y        = "Module Eigengene Value"
  )

ggsave("Manuscript_Figures/WGCNA/Figure3G_ME_Violin.pdf",
       plot = p_violin, width = 10, height = 10)
ggsave("Manuscript_Figures/WGCNA/Figure3G_ME_Violin.png",
       plot = p_violin, width = 10, height = 10, dpi = 600)
message("Figure 3G saved.")

# ==============================================================================
# STEP 12: Gene Count Per Module -- Figure 3H (UNCHANGED)
# ==============================================================================
mc_table <- as.data.frame(table(moduleColors), stringsAsFactors = FALSE)
colnames(mc_table) <- c("Module", "GeneCount")
mc_table <- mc_table[mc_table$Module != "grey", ]
mc_table <- mc_table[order(mc_table$GeneCount, decreasing = TRUE), ]

p_bar <- ggplot(mc_table,
                aes(x = reorder(Module, -GeneCount),
                    y = GeneCount,
                    fill = Module)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = setNames(as.character(mc_table$Module),
                      as.character(mc_table$Module))
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "3H. Number of Genes Per Module",
       x = "Module", y = "Gene Count")

ggsave("Manuscript_Figures/WGCNA/Figure3H_Gene_Count_Bar.pdf",
       plot = p_bar, width = 10, height = 6)
ggsave("Manuscript_Figures/WGCNA/Figure3H_Gene_Count_Bar.png",
       plot = p_bar, width = 10, height = 6, dpi = 600)
message("Figure 3H saved.")

# ==============================================================================
# STEP 13: Hub Gene Extraction from ALL Significant Modules (UNCHANGED)
# ==============================================================================
message("Extracting hub genes from all significant modules...")

gene_name_vec <- names(moduleColors)
if (is.null(gene_name_vec)) gene_name_vec <- colnames(datExpr)

all_hub_genes <- character(0)

for (mod_col in sig_module_cols) {
  kME_col <- paste0("kME", mod_col)
  if (!kME_col %in% colnames(kME_all)) {
    message(sprintf("  Skipping %s - kME column not found.", mod_col))
    next
  }

  inModule_genes <- gene_name_vec[moduleColors == mod_col]
  avail_genes    <- intersect(inModule_genes, rownames(kME_all))

  if (length(avail_genes) == 0) {
    message(sprintf("  Module %s: no genes in kME_all.", mod_col))
    next
  }

  kme_vals   <- kME_all[avail_genes, kME_col]
  hub_ranked <- avail_genes[order(abs(kme_vals), decreasing = TRUE)]
  top_hubs   <- hub_ranked[1:min(30, length(hub_ranked))]

  all_hub_genes <- unique(c(all_hub_genes, top_hubs))

  message(sprintf("  Module %s: %d genes | top hub = %s (kME = %.3f)",
                  mod_col,
                  length(avail_genes),
                  hub_ranked[1],
                  max(abs(kme_vals))))
}

message(sprintf("Total hub candidates for ML: %d", length(all_hub_genes)))

primary_genes_all   <- gene_name_vec[moduleColors == primary_module_col]
primary_genes_avail <- intersect(primary_genes_all, rownames(kME_all))
primary_kme_vals    <- kME_all[primary_genes_avail, primary_kME_col]
primary_ranked      <- primary_genes_avail[
  order(abs(primary_kme_vals), decreasing = TRUE)
]
top200_string <- primary_ranked[1:min(200, length(primary_ranked))]

write.csv(data.frame(Gene = all_hub_genes),
          "Processed_Data/Candidate_Hubs_for_ML.csv",
          row.names = FALSE)
write.csv(data.frame(Gene = top200_string),
          "Processed_Data/Step3_Top200_for_STRING.csv",
          row.names = FALSE)

message(sprintf("Saved: %d hub genes for ML | %d for STRING.",
                length(all_hub_genes), length(top200_string)))

# ==============================================================================
# STEP 14: Save Complete WGCNA Results (UNCHANGED)
# ==============================================================================
saveRDS(
  list(
    datExpr            = datExpr,
    moduleColors       = moduleColors,
    gene_name_vec      = gene_name_vec,
    MEs                = MEs_clean,
    traitData          = traitData,
    moduleTraitCor     = moduleTraitCor,
    moduleTraitPvalue  = moduleTraitPvalue,
    sig_module_names   = sig_module_names,
    sig_module_cols    = sig_module_cols,
    primary_module_col = primary_module_col,
    primary_kME_col    = primary_kME_col,
    kME_all            = kME_all,
    GS_ALS             = GS_ALS,
    chosen_power       = chosen_power,
    achieved_r2        = achieved_r2,
    achieved_meank     = achieved_meank,
    n_modules          = n_modules,
    top200_string      = top200_string,
    valid_cor          = valid_cor,
    valid_pval         = valid_pval
  ),
  "Processed_Data/Step3_WGCNA_Data.rds"
)

message("================================================================")
message("Script 02 complete.")
message(sprintf("  Power            : %d (R2=%.3f, mean.k=%.2f)",
                chosen_power, achieved_r2, achieved_meank))
message(sprintf("  Modules total    : %d", n_modules))
message(sprintf("  Grey (unassigned): %d genes (%.1f%%)", n_grey, pct_grey))
message(sprintf("  Significant      : %d modules", length(sig_module_cols)))
message(sprintf("  Primary module   : %s (r=%.3f)",
                primary_module_col,
                valid_cor[primary_module, "ALS"]))
message(sprintf("  Hub candidates   : %d genes for ML", length(all_hub_genes)))
message("================================================================")
