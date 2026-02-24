# ==============================================================================
# MASTER SCRIPT 3: Machine Learning Feature Selection (Step 5)
# ==============================================================================
packages <- c("randomForest", "glmnet", "caret", "ggplot2", "dplyr", "tidyr", "ggVennDiagram", "pROC", "scales", "viridis")
for (pkg in packages) { if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg); library(pkg, character.only = TRUE) }
set.seed(2026)
dir.create("Manuscript_Figures/Step5_MachineLearning", recursive = TRUE, showWarnings = FALSE)

theme_publication <- function(base_size = 14) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 2, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      legend.position = "right", legend.title = element_text(face = "bold", size = base_size - 2),
      legend.text = element_text(size = base_size - 2),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3), panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

expr_data <- readRDS("Processed_Data/Step1_ComBat_Data.rds")
hubs <- read.csv("Processed_Data/Candidate_Hubs_for_ML.csv")$gene
ml_data <- t(expr_data$expr_combat[hubs, ])
ml_labels <- factor(expr_data$group_factor)

# A. RANDOM FOREST
rf_model <- randomForest(x = ml_data, y = ml_labels, ntree = 500, importance = TRUE)
oob_long <- data.frame(Trees = 1:500, OOB_Error = rf_model$err.rate[, "OOB"], ALS_Error = rf_model$err.rate[, "ALS"], Control_Error = rf_model$err.rate[, "Control"]) %>% pivot_longer(cols = c(OOB_Error, ALS_Error, Control_Error), names_to = "Error_Type", values_to = "Error_Rate")
oob_long$Error_Type <- factor(oob_long$Error_Type, levels = c("OOB_Error", "ALS_Error", "Control_Error"), labels = c("OOB", "ALS", "Control"))

p_oob <- ggplot(oob_long, aes(x = Trees, y = Error_Rate, color = Error_Type)) + geom_line(linewidth = 1.2, alpha = 0.9) + scale_color_manual(values = c("OOB" = "#000000", "ALS" = "#E64B35", "Control" = "#4DBBD5"), name = "Error Type") + scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.45), expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + labs(title = "Random Forest Out-of-Bag Error Rate", x = "Number of Trees", y = "Error Rate") + theme_publication() + theme(legend.position = "right", legend.justification = "center", legend.background = element_rect(fill = "white", color = NA), legend.key = element_rect(fill = "white"), legend.key.size = unit(0.9, "lines"), legend.spacing.y = unit(0.15, "cm"), legend.margin = ggplot2::margin(5, 5, 5, 5, "pt"))
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5A_RF_OOB.pdf", plot = p_oob, width = 8, height = 6)
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5A_RF_OOB.png", plot = p_oob, width = 8, height = 6, dpi = 600)

rf_imp <- importance(rf_model, type = 1)
rf_imp_df <- data.frame(Gene = rownames(rf_imp), MeanDecreaseAccuracy = rf_imp[, 1]) %>% arrange(desc(MeanDecreaseAccuracy)) %>% head(15)
rf_genes <- data.frame(Gene = rownames(rf_imp), MeanDecreaseAccuracy = rf_imp[, 1]) %>% arrange(desc(MeanDecreaseAccuracy)) %>% head(30) %>% pull(Gene)

p_rf <- ggplot(rf_imp_df, aes(x = MeanDecreaseAccuracy, y = reorder(Gene, MeanDecreaseAccuracy), fill = MeanDecreaseAccuracy)) + geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = reorder(Gene, MeanDecreaseAccuracy), yend = reorder(Gene, MeanDecreaseAccuracy)), color = "gray70", linewidth = 0.5) + geom_point(size = 5, shape = 21, color = "black", stroke = 0.8) + scale_fill_gradient2(low = "#4DBBD5", mid = "#FED439", high = "#E64B35", midpoint = median(rf_imp_df$MeanDecreaseAccuracy), name = "Importance") + labs(title = "Top Genes by Random Forest Importance", x = "Mean Decrease in Accuracy", y = NULL) + theme_publication() + theme(panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3), axis.text.y = element_text(face = "italic"))
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5B_RF_Importance.pdf", plot = p_rf, width = 8, height = 6)
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5B_RF_Importance.png", plot = p_rf, width = 8, height = 6, dpi = 600)

# C. SVM-RFE
svm_ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 5, verbose = FALSE)
svm_rfe <- rfe(x = ml_data, y = ml_labels, sizes = c(1:30, 35, 40, 45, 50), rfeControl = svm_ctrl, method = "svmLinear")
svm_genes <- predictors(svm_rfe)
svm_res <- svm_rfe$results
optimal_idx <- which.max(svm_res$Accuracy)

p_svm <- ggplot(svm_res, aes(x = Variables, y = Accuracy)) + geom_line(color = "#4DBBD5", linewidth = 1.3) + geom_point(size = 3, color = "#4DBBD5", fill = "white", shape = 21, stroke = 1.5) + geom_point(data = svm_res[optimal_idx, ], aes(x = Variables, y = Accuracy), color = "#E64B35", size = 6, shape = 18) + geom_vline(xintercept = svm_res$Variables[optimal_idx], linetype = "dashed", color = "#E64B35", linewidth = 0.6, alpha = 0.7) + annotate("segment", x = svm_res$Variables[optimal_idx], xend = svm_res$Variables[optimal_idx] + 8, y = svm_res$Accuracy[optimal_idx], yend = svm_res$Accuracy[optimal_idx], arrow = arrow(length = unit(0.2, "cm")), color = "#E64B35", linewidth = 0.7) + annotate("text", x = svm_res$Variables[optimal_idx] + 10, y = svm_res$Accuracy[optimal_idx], label = paste0("Optimal\n", svm_res$Variables[optimal_idx], " features"), hjust = 0, size = 4, fontface = "bold", color = "#E64B35") + scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(min(svm_res$Accuracy) - 0.02, max(svm_res$Accuracy) + 0.02)) + scale_x_continuous(breaks = seq(0, max(svm_res$Variables), by = 10)) + labs(title = "SVM-RFE Cross-Validation Performance", x = "Number of Selected Features", y = "Cross-Validation Accuracy") + theme_publication()
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5C_SVM_RFE.pdf", plot = p_svm, width = 8, height = 6)
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5C_SVM_RFE.png", plot = p_svm, width = 8, height = 6, dpi = 600)

# D-E. LASSO
y_num <- ifelse(ml_labels == "ALS", 1, 0)
cv_lasso <- cv.glmnet(ml_data, y_num, family = "binomial", alpha = 1, nfolds = 10)
lasso_full <- glmnet(ml_data, y_num, family = "binomial", alpha = 1)
coef_matrix <- as.matrix(coef(lasso_full)); gene_names <- rownames(coef_matrix)[-1]
coef_df <- data.frame(Lambda = lasso_full$lambda, LogLambda = log(lasso_full$lambda))
for (i in 1:length(gene_names)) { coef_df[[gene_names[i]]] <- coef_matrix[i + 1, ] }
coef_long <- coef_df %>% pivot_longer(cols = -c(Lambda, LogLambda), names_to = "Gene", values_to = "Coefficient")

top_coefs <- coef(cv_lasso, s = "lambda.min")
top_genes_lasso <- names(sort(abs(top_coefs[-1, ]), decreasing = TRUE)[1:15])
coef_long$Gene_Group <- ifelse(coef_long$Gene %in% top_genes_lasso, coef_long$Gene, "Other")
gene_colors <- setNames(c(viridis::viridis(length(top_genes_lasso), option = "D"), "gray70"), c(top_genes_lasso, "Other"))

p_lasso_path <- ggplot(coef_long, aes(x = LogLambda, y = Coefficient, group = Gene, color = Gene_Group)) + geom_line(data = filter(coef_long, Gene_Group == "Other"), linewidth = 0.3, alpha = 0.2) + geom_line(data = filter(coef_long, Gene_Group != "Other"), linewidth = 1.2, alpha = 0.8) + geom_vline(xintercept = log(cv_lasso$lambda.min), linetype = "dashed", color = "#E64B35", linewidth = 1) + geom_vline(xintercept = log(cv_lasso$lambda.1se), linetype = "dotted", color = "#4DBBD5", linewidth = 1) + scale_color_manual(values = gene_colors, guide = "none") + labs(title = "LASSO Coefficient Regularization Path", x = expression(bold(Log(lambda))), y = "Standardized Coefficients") + theme_publication()
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5D_LASSO_Path.pdf", plot = p_lasso_path, width = 8, height = 6)
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5D_LASSO_Path.png", plot = p_lasso_path, width = 8, height = 6, dpi = 600)

cv_data <- data.frame(LogLambda = log(cv_lasso$lambda), Deviance = cv_lasso$cvm, SE_upper = cv_lasso$cvup, SE_lower = cv_lasso$cvlo, NonZero = cv_lasso$nzero)
p_lasso_cv <- ggplot(cv_data, aes(x = LogLambda, y = Deviance)) + geom_errorbar(aes(ymin = SE_lower, ymax = SE_upper), width = 0.1, color = "gray50", linewidth = 0.4, alpha = 0.6) + geom_point(color = "#E64B35", size = 2.5, alpha = 0.8) + geom_vline(xintercept = log(cv_lasso$lambda.min), linetype = "dashed", color = "#4DBBD5", linewidth = 1) + geom_vline(xintercept = log(cv_lasso$lambda.1se), linetype = "dotted", color = "#00A087", linewidth = 1) + annotate("text", x = log(cv_lasso$lambda.min), y = max(cv_data$SE_upper) * 0.97, label = expression(lambda[min]), hjust = 1.2, fontface = "bold", size = 4.5, color = "#4DBBD5") + annotate("text", x = log(cv_lasso$lambda.1se), y = max(cv_data$SE_upper) * 0.97, label = expression(lambda["1SE"]), hjust = -0.2, fontface = "bold", size = 4.5, color = "#00A087") + scale_x_continuous(sec.axis = sec_axis(~ ., breaks = log(cv_lasso$lambda)[round(seq(1, length(cv_lasso$lambda), length.out = 10))], labels = cv_data$NonZero[round(seq(1, nrow(cv_data), length.out = 10))], name = "Number of Non-Zero Coefficients")) + labs(title = "LASSO 10-Fold Cross-Validation", x = expression(bold(Log(lambda))), y = "Binomial Deviance") + theme_publication()
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5E_LASSO_CV.pdf", plot = p_lasso_cv, width = 8, height = 6)
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5E_LASSO_CV.png", plot = p_lasso_cv, width = 8, height = 6, dpi = 600)

lasso_genes <- setdiff(rownames(top_coefs)[top_coefs[, 1] != 0], "(Intercept)")

# F. TRIPLE CONSENSUS
gene_lists <- list("Random Forest" = rf_genes, "LASSO" = lasso_genes, "SVM-RFE" = svm_genes)
consensus_genes <- Reduce(intersect, gene_lists)

p_venn <- ggVennDiagram(gene_lists, label = "count", label_alpha = 0, edge_size = 1.5, set_size = 5) + scale_fill_gradient(low = "#FFF5EB", high = "#E64B35", name = "Gene Count") + scale_color_manual(values = c("#4DBBD5", "#E64B35", "#00A087")) + labs(title = paste0("Triple-Algorithm Consensus (n = ", length(consensus_genes), " genes)")) + theme_void(base_size = 14) + theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), legend.position = "right", legend.title = element_text(face = "bold"))
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5F_Venn.pdf", plot = p_venn, width = 8, height = 6)
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5F_Venn.png", plot = p_venn, width = 8, height = 6, dpi = 600)
write.csv(data.frame(Gene = consensus_genes), "Processed_Data/Final_21Gene_Signature.csv", row.names = FALSE)

# G. CONSENSUS ROC
glm_model <- glm(ml_labels ~ ., data = data.frame(ml_data[, consensus_genes], ml_labels), family = binomial)
prob <- predict(glm_model, type = "response")
roc_obj <- roc(ml_labels, prob, levels = c("Control", "ALS"), quiet = TRUE)

# Calculate TRUE shape Confidence Intervals (Reviewer-proof)
ci_shape <- ci.se(roc_obj, specificities = roc_obj$specificities)

roc_data <- data.frame(
  Specificity = roc_obj$specificities,
  Sensitivity = roc_obj$sensitivities,
  CI_Lower = as.numeric(ci_shape[, 1]), # True lower bound
  CI_Upper = as.numeric(ci_shape[, 3])  # True upper bound
)
ci_obj <- ci.auc(roc_obj, conf.level = 0.95)

p_roc <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  # Using the mathematically true confidence bounds
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#E64B35", alpha = 0.2) + 
  geom_line(color = "#E64B35", linewidth = 2.5) +
  annotate("text", x = 0.6, y = 0.25,
           label = paste0("AUC = ", round(auc(roc_obj), 3), "\n",
                          "95% CI: [", round(ci_obj[1], 3), "-", round(ci_obj[3], 3), "]"),
           hjust = 0, size = 5.5, fontface = "bold", color = "#E64B35") +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "Consensus Gene Signature ROC Curve",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_publication() +
  theme(panel.grid.major = element_line(color = "gray90", linewidth = 0.3), aspect.ratio = 1)

ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5G_Consensus_ROC.pdf", plot = p_roc, width = 7, height = 7)
ggsave("Manuscript_Figures/Step5_MachineLearning/Figure5G_Consensus_ROC.png", plot = p_roc, width = 7, height = 7, dpi = 600)
message("Step 5 Complete.")
