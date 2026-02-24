# ==============================================================================
# MASTER SCRIPT 2: WGCNA (Step 3)
# ==============================================================================
library(WGCNA); library(dplyr); library(ggplot2); library(ggrepel); library(igraph)

set.seed(1122)
tryCatch({ allowWGCNAThreads() }, error = function(e) message("Running in single-thread mode for stability."))
dir.create("Manuscript_Figures/WGCNA", recursive = TRUE, showWarnings = FALSE)

step1_data <- readRDS("Processed_Data/Step1_ComBat_Data.rds")
expr_combat <- step1_data$expr_combat
group_factor <- step1_data$group_factor

datExpr <- t(expr_combat) 
gene_vars <- apply(datExpr, 2, var)
datExpr_filtered <- datExpr[, names(sort(gene_vars, decreasing = TRUE))[1:5000]]

# PLOT 1 & 2: Soft-Thresholding
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr_filtered, powerVector = powers, verbose = 0)
df_sft <- data.frame(Power = sft$fitIndices[,1], Fit = -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], Connectivity = sft$fitIndices[,5])

p1 <- ggplot(df_sft, aes(x = Power, y = Fit)) + geom_hline(yintercept = 0.85, linetype = "dashed", color = "red", linewidth = 0.8) + geom_point(size = 8, color = "#BC3C29FF", fill = "#BC3C29FF", shape = 21) + geom_text(aes(label = Power), color = "white", size = 4, fontface = "bold") + theme_classic(base_size = 14) + labs(title = "Plot 1. Scale-Free Topology Fit", x = "Soft Threshold (Power)", y = expression("Scale Free Topology R"^2)) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("Manuscript_Figures/WGCNA/Plot1_ScaleFree_Topology.pdf", plot = p1, width = 6, height = 5)
ggsave("Manuscript_Figures/WGCNA/Plot1_ScaleFree_Topology.png", plot = p1, width = 6, height = 5, dpi = 600)

p2 <- ggplot(df_sft, aes(x = Power, y = Connectivity)) + geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) + geom_point(size = 8, color = "#0072B5FF", fill = "#0072B5FF", shape = 21) + geom_text(aes(label = Power), color = "white", size = 4, fontface = "bold") + theme_classic(base_size = 14) + labs(title = "Plot 2. Mean Connectivity", x = "Soft Threshold (Power)", y = "Mean Connectivity") + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + coord_cartesian(ylim = c(min(df_sft$Connectivity) - 100, max(df_sft$Connectivity) + 100)) 
ggsave("Manuscript_Figures/WGCNA/Plot2_Mean_Connectivity.pdf", plot = p2, width = 6, height = 5)
ggsave("Manuscript_Figures/WGCNA/Plot2_Mean_Connectivity.png", plot = p2, width = 6, height = 5, dpi = 600)

# PLOT 3: Sample Clustering
sampleTree <- hclust(dist(datExpr_filtered), method = "average")
pdf("Manuscript_Figures/WGCNA/Plot3_Sample_Clustering.pdf", width = 12, height = 6)
par(mar = c(0, 5, 4, 2), cex = 0.8)
plot(sampleTree, main = "Plot 3. Sample Clustering to Detect Outliers", sub = "", xlab = "", cex.main = 1.4, labels = FALSE)
abline(h = 200, col = "red", lty = 2, lwd = 2)
dev.off()
png("Manuscript_Figures/WGCNA/Plot3_Sample_Clustering.png", width = 12, height = 6, units = "in", res = 600)
par(mar = c(0, 5, 4, 2), cex = 0.8)
plot(sampleTree, main = "Plot 3. Sample Clustering to Detect Outliers", sub = "", xlab = "", cex.main = 1.4, labels = FALSE)
abline(h = 200, col = "red", lty = 2, lwd = 2)
dev.off()

# CONSTRUCT NETWORK
net <- blockwiseModules(datExpr_filtered, power = 8, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "Processed_Data/TOM", verbose = 0)
moduleColors <- labels2colors(net$colors)
MEs <- orderMEs(moduleEigengenes(datExpr_filtered, moduleColors)$eigengenes)

# PLOT 4: Gene Dendrogram
pdf("Manuscript_Figures/WGCNA/Plot4_Gene_Dendrogram.pdf", width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module Colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Plot 4. Gene Dendrogram and Module Colors")
dev.off()
png("Manuscript_Figures/WGCNA/Plot4_Gene_Dendrogram.png", width = 12, height = 6, units = "in", res = 600)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module Colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Plot 4. Gene Dendrogram and Module Colors")
dev.off()

# PLOT 5: Module-Trait Heatmap
traitData <- data.frame(ALS = ifelse(group_factor == "ALS", 1, 0))
rownames(traitData) <- rownames(datExpr_filtered)
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr_filtered))
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf("Manuscript_Figures/WGCNA/Plot5_ModuleTrait_Heatmap.pdf", width = 6, height = 10)
par(mar = c(6, 9, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "ALS Diagnosis", yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1, 1), main = "Plot 5. Module-Trait Relationships")
dev.off()
png("Manuscript_Figures/WGCNA/Plot5_ModuleTrait_Heatmap.png", width = 6, height = 10, units = "in", res = 600)
par(mar = c(6, 9, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "ALS Diagnosis", yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1, 1), main = "Plot 5. Module-Trait Relationships")
dev.off()

# PLOT 6: Eigengene Dendrogram
pdf("Manuscript_Figures/WGCNA/Plot6_Eigengene_Dendrogram.pdf", width = 10, height = 6)
plotEigengeneNetworks(MEs, "Plot 6. Eigengene Dendrogram", marDendro = c(0, 4, 2, 0), plotHeatmaps = FALSE)
dev.off()
png("Manuscript_Figures/WGCNA/Plot6_Eigengene_Dendrogram.png", width = 10, height = 6, units = "in", res = 600)
plotEigengeneNetworks(MEs, "Plot 6. Eigengene Dendrogram", marDendro = c(0, 4, 2, 0), plotHeatmaps = FALSE)
dev.off()

# PLOT 7: Correlation Heatmap
modCor <- cor(MEs, use = "p")
modCorP <- corPvalueStudent(modCor, nrow(datExpr_filtered))
format_p <- function(p) { ifelse(p < 1e-100, "<1e-100", formatC(p, format = "e", digits = 1)) }
formatted_p_vals <- apply(modCorP, c(1,2), format_p)
textMatCor <- paste(signif(modCor, 2), "\n(", formatted_p_vals, ")", sep = "")
dim(textMatCor) <- dim(modCor)

pdf("Manuscript_Figures/WGCNA/Plot7_Module_Correlation_Heatmap.pdf", width = 11, height = 10)
par(mar = c(6, 8, 4, 2)) 
labeledHeatmap(Matrix = modCor, xLabels = names(MEs), yLabels = names(MEs), ySymbols = names(MEs), xSymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatCor, setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1, 1), main = "Plot 7. Top Module Eigengene Correlation Heatmap", cex.main = 1.5, cex.lab = 0.8, xLabelsAngle = 45) 
dev.off()
png("Manuscript_Figures/WGCNA/Plot7_Module_Correlation_Heatmap.png", width = 11, height = 10, units = "in", res = 600)
par(mar = c(6, 8, 4, 2)) 
labeledHeatmap(Matrix = modCor, xLabels = names(MEs), yLabels = names(MEs), ySymbols = names(MEs), xSymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatCor, setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1, 1), main = "Plot 7. Top Module Eigengene Correlation Heatmap", cex.main = 1.5, cex.lab = 0.8, xLabelsAngle = 45) 
dev.off()

valid_modules <- moduleTraitCor[rownames(moduleTraitCor) != "MEgrey", , drop = FALSE]
top_modules <- rownames(valid_modules)[order(abs(valid_modules[,"ALS"]), decreasing = TRUE)][1:2]
top_module_colors <- gsub("ME", "", top_modules)

# PLOT 8: GS vs MM
GS_ALS <- as.numeric(cor(datExpr_filtered, traitData$ALS, use = "p"))
names(GS_ALS) <- colnames(datExpr_filtered)

for (mod_col in top_module_colors) {
  inModule <- (moduleColors == mod_col)
  MM <- as.numeric(cor(datExpr_filtered[, inModule], MEs[, paste0("ME", mod_col)], use = "p"))
  df_plot <- data.frame(MM = MM, GS = GS_ALS[inModule], Gene = colnames(datExpr_filtered)[inModule])
  top_hub <- df_plot %>% arrange(desc(abs(MM))) %>% head(10)
  
  p8 <- ggplot(df_plot, aes(x = MM, y = GS)) + geom_point(fill = mod_col, color = "black", alpha = 0.8, size = 2.5, shape = 21, stroke = 0.3) + geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") + geom_label_repel(data = top_hub, aes(label = Gene), size = 4.5, fontface = "bold.italic", color = "black", fill = alpha("white", 0.8), label.size = 0, nudge_x = 1 - top_hub$MM, direction = "y", hjust = 0, segment.color = "grey60", segment.size = 0.4, max.overlaps = Inf) + theme_classic(base_size = 15) + theme(plot.title = element_text(face = "bold", hjust = 0.5)) + scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) + labs(title = paste0("Plot 8. GS vs MM — ", toupper(mod_col), " Module"), x = paste("Module Membership in", mod_col, "module"), y = "Gene Significance (ALS Trait)")
  ggsave(paste0("Manuscript_Figures/WGCNA/Plot8_GS_vs_MM_", mod_col, ".pdf"), plot = p8, width = 8, height = 6.5)
  ggsave(paste0("Manuscript_Figures/WGCNA/Plot8_GS_vs_MM_", mod_col, ".png"), plot = p8, width = 8, height = 6.5, dpi = 600)
}

# PLOT 9: Gene Count Bar Plot
module_counts <- as.data.frame(table(moduleColors)) %>% arrange(desc(Freq))
colnames(module_counts) <- c("Module", "GeneCount")
p9 <- ggplot(module_counts, aes(x = reorder(Module, -GeneCount), y = GeneCount, fill = Module)) + geom_bar(stat = "identity", color = "black", linewidth = 0.3) + scale_fill_manual(values = setNames(as.character(module_counts$Module), as.character(module_counts$Module))) + theme_classic(base_size = 13) + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none") + labs(title = "Plot 9. Number of Genes Per Module", x = "Module", y = "Gene Count")
ggsave("Manuscript_Figures/WGCNA/Plot9_Gene_Count_BarPlot.pdf", plot = p9, width = 10, height = 6)
ggsave("Manuscript_Figures/WGCNA/Plot9_Gene_Count_BarPlot.png", plot = p9, width = 10, height = 6, dpi = 600)

# PLOT 10: Hub Gene Network
TOM_file <- list.files("Processed_Data", pattern = "TOM.*\\.RData$", full.names = TRUE)[1]
load(TOM_file)
TOM_matrix <- as.matrix(TOM); colnames(TOM_matrix) <- rownames(TOM_matrix) <- colnames(datExpr_filtered)
top_mod <- top_module_colors[1]
inModule <- moduleColors == top_mod
hub_genes <- names(sort(rowSums(TOM_matrix[inModule, inModule]), decreasing = TRUE))[1:30]
adj_matrix <- TOM_matrix[hub_genes, hub_genes]; diag(adj_matrix) <- 0; adj_matrix[adj_matrix < 0.05] <- 0  
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

pdf("Manuscript_Figures/WGCNA/Plot10_Hub_Gene_Network.pdf", width = 10, height = 9)
set.seed(42)
plot(g, vertex.label = V(g)$name, vertex.size = 12, vertex.color = "#89CFF0", vertex.label.cex = 0.9, vertex.label.color = "black", vertex.frame.color = "black", edge.color = "grey70", edge.width = E(g)$weight * 15, layout = layout_with_kk(g), main = paste0("Plot 10. Hub Gene Network — ", toupper(top_mod), " Module"))
dev.off()
png("Manuscript_Figures/WGCNA/Plot10_Hub_Gene_Network.png", width = 10, height = 9, units = "in", res = 600)
set.seed(42)
plot(g, vertex.label = V(g)$name, vertex.size = 12, vertex.color = "#89CFF0", vertex.label.cex = 0.9, vertex.label.color = "black", vertex.frame.color = "black", edge.color = "grey70", edge.width = E(g)$weight * 15, layout = layout_with_kk(g), main = paste0("Plot 10. Hub Gene Network — ", toupper(top_mod), " Module"))
dev.off()

saveRDS(list(datExpr = datExpr_filtered, moduleColors = moduleColors, MEs = MEs, traitData = traitData, moduleTraitCor = moduleTraitCor, moduleTraitPvalue = moduleTraitPvalue, top_module_colors = top_module_colors, GS_ALS = GS_ALS, module_genes = net$colors), "Processed_Data/Step3_WGCNA_Data.rds")
write.csv(data.frame(gene = hub_genes), "Processed_Data/Candidate_Hubs_for_ML.csv", row.names = FALSE)
message("Step 3 Complete.")
