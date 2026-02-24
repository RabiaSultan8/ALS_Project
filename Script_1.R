# ==============================================================================
# MASTER SCRIPT 1: Data Prep & Differential Expression (Steps 1 & 2)
# ==============================================================================

# ── 0. Reproducibility & Dependency Management ────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("sva", "illuminaHumanv3.db", "limma", "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db")
cran_pkgs <- c("GEOquery", "ggplot2", "dplyr", "patchwork", "ggsci", "ggrepel", "circlize")

for (pkg in bioc_pkgs) { if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE) }
for (pkg in cran_pkgs) { if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg) }

library(GEOquery); library(sva); library(ggplot2); library(dplyr); library(illuminaHumanv3.db)
library(patchwork); library(ggsci); library(limma); library(ggrepel); library(ComplexHeatmap); library(circlize)

set.seed(1122)
dir.create("Manuscript_Figures", showWarnings = FALSE)
dir.create("Processed_Data", showWarnings = FALSE)

# ==============================================================================
# Step 1: Data Prep & Publication-Ready Batch Correction (PCA)
# ==============================================================================
message("Downloading/Loading GEO data...")
gse112676 <- getGEO("GSE112676", getGPL = FALSE)[[1]]
gse112680 <- getGEO("GSE112680", getGPL = FALSE)[[1]]

expr_112676 <- exprs(gse112676)
expr_112680 <- exprs(gse112680)

group_112676 <- toupper(trimws(as.character(pData(gse112676)$`diagnosis:ch1`)))
group_112680 <- toupper(trimws(as.character(pData(gse112680)$`diagnosis:ch1`)))

keep_112676 <- group_112676 %in% c("ALS", "CON")
keep_112680 <- group_112680 %in% c("ALS", "CON")

expr_112676 <- expr_112676[, keep_112676]; group_112676 <- group_112676[keep_112676]
expr_112680 <- expr_112680[, keep_112680]; group_112680 <- group_112680[keep_112680]

map_using_bioc <- function(expr) {
  probes <- rownames(expr)
  symbols <- mapIds(illuminaHumanv3.db, keys = probes, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  keep <- !is.na(symbols) & symbols != ""
  df <- data.frame(expr[keep, ]); df$gene <- symbols[keep]
  df_agg <- df %>% group_by(gene) %>% summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
  mat <- as.matrix(df_agg[, -1]); rownames(mat) <- df_agg$gene
  return(mat)
}

expr_112676_g <- map_using_bioc(expr_112676)
expr_112680_g <- map_using_bioc(expr_112680)

common_genes <- intersect(rownames(expr_112676_g), rownames(expr_112680_g))
expr_combined <- cbind(expr_112676_g[common_genes, ], expr_112680_g[common_genes, ])

batch <- c(rep("GSE112676", sum(keep_112676)), rep("GSE112680", sum(keep_112680)))
group_factor <- factor(ifelse(c(group_112676, group_112680) == "ALS", "ALS", "Control"))

if (max(expr_combined, na.rm = TRUE) > 100) expr_combined <- log2(expr_combined + 1)

# Pre-Correction PCA
pca_pre <- prcomp(t(expr_combined), scale. = TRUE)
df_pca_pre <- data.frame(PC1 = pca_pre$x[,1], PC2 = pca_pre$x[,2], Batch = batch, Diagnosis = group_factor)

pub_theme <- theme_classic(base_size = 14) + theme(axis.text = element_text(color = "black"), legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5))

plot_pre <- ggplot(df_pca_pre, aes(x = PC1, y = PC2, color = Batch, shape = Diagnosis)) + geom_point(alpha = 0.7, size = 2) + scale_color_jco() + pub_theme + labs(title = "Before Batch Correction", x = "PC1", y = "PC2")

# ComBat Correction
mod <- model.matrix(~ group_factor)
expr_combat <- ComBat(dat = expr_combined, batch = batch, mod = mod, par.prior = TRUE)

# Post-Correction PCA
pca_post <- prcomp(t(expr_combat), scale. = TRUE)
df_pca_post <- data.frame(PC1 = pca_post$x[,1], PC2 = pca_post$x[,2], Batch = batch, Diagnosis = group_factor)

plot_post <- ggplot(df_pca_post, aes(x = PC1, y = PC2, color = Batch, shape = Diagnosis)) + geom_point(alpha = 0.7, size = 2) + scale_color_jco() + pub_theme + labs(title = "After ComBat Correction", x = "PC1", y = "PC2")

fig1 <- plot_pre + plot_post + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("Manuscript_Figures/Figure1_Batch_PCA.pdf", plot = fig1, width = 12, height = 6)
ggsave("Manuscript_Figures/Figure1_Batch_PCA.png", plot = fig1, width = 12, height = 6, dpi = 600)

saveRDS(list(expr_combat = expr_combat, group_factor = group_factor, batch = batch), "Processed_Data/Step1_ComBat_Data.rds")
message("Step 1 Complete.")

# ==============================================================================
# Step 2: Differential Expression Analysis & Publication-Ready Figures
# ==============================================================================
design <- model.matrix(~0 + group_factor)
colnames(design) <- c("Control", "ALS")

fit <- lmFit(expr_combat, design)
contr <- makeContrasts(ALS_vs_Control = ALS - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))

deg_all <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
deg_all$gene <- rownames(deg_all)
deg_sig <- deg_all %>% filter(abs(logFC) >= 0.25 & adj.P.Val < 0.05)

write.csv(deg_all, "Processed_Data/Step2_DEG_all.csv", row.names = FALSE)
write.csv(deg_sig, "Processed_Data/Step2_DEG_sig.csv", row.names = FALSE)

vol <- deg_all %>% mutate(negLogFDR = -log10(adj.P.Val), Sig = case_when(adj.P.Val < 0.05 & logFC >= 0.25 ~ "Upregulated", adj.P.Val < 0.05 & logFC <= -0.25 ~ "Downregulated", TRUE ~ "Not Significant"))
top_genes <- vol %>% filter(Sig != "Not Significant") %>% arrange(adj.P.Val) %>% head(20)

plot_volcano <- ggplot(vol, aes(x = logFC, y = negLogFDR, color = Sig)) + 
  geom_point(alpha = 0.8, size = 1.8) + 
  scale_color_manual(values = c("Downregulated" = "#008B8B", "Upregulated" = "#DC143C", "Not Significant" = "#E0E0E0")) + 
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey30", linewidth = 0.5) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", linewidth = 0.5) + 
  geom_text_repel(data = top_genes, aes(label = gene), size = 4.5, fontface = "italic", box.padding = 0.5, max.overlaps = Inf, color = "black", segment.color = NA) + 
  theme_classic(base_size = 15) + 
  theme(legend.position = "top", legend.title = element_blank(), plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
  labs(title = "Differentially Expressed Genes", x = bquote(log[2]("Fold Change")), y = bquote(-log[10]("FDR")))

ggsave("Manuscript_Figures/Figure2A_Volcano.pdf", plot = plot_volcano, width = 7.5, height = 7.5)
ggsave("Manuscript_Figures/Figure2A_Volcano.png", plot = plot_volcano, width = 7.5, height = 7.5, dpi = 600)

top100 <- deg_sig %>% arrange(adj.P.Val) %>% head(100) %>% pull(gene)
hm_data <- expr_combat[top100, ]
hm_scaled <- t(scale(t(hm_data)))

ha_col <- HeatmapAnnotation(Diagnosis = group_factor, Batch = batch, col = list(Diagnosis = c("ALS" = "#DC143C", "Control" = "#008B8B"), Batch = c("GSE112676" = "#4A4A4A", "GSE112680" = "#B3B3B3")), show_annotation_name = TRUE, annotation_name_gp = gpar(fontface = "bold"))
col_fun <- colorRamp2(c(-2, 0, 2), c("#6495ED", "#FFFACD", "#FF6347"))

ht <- Heatmap(hm_scaled, name = "Z-score", top_annotation = ha_col, col = col_fun, show_row_names = FALSE, show_column_names = FALSE, cluster_columns = TRUE, cluster_rows = TRUE, show_row_dend = TRUE, show_column_dend = FALSE, column_title = "Top 100 DEGs Expression Profile", column_title_gp = gpar(fontsize = 16, fontface = "bold"), use_raster = TRUE)

pdf("Manuscript_Figures/Figure2B_Heatmap.pdf", width = 8.5, height = 7.5)
draw(ht, merge_legend = TRUE)
dev.off()
png("Manuscript_Figures/Figure2B_Heatmap.png", width = 8.5, height = 7.5, units = "in", res = 600)
draw(ht, merge_legend = TRUE)
dev.off()
message("Step 2 Complete.")
