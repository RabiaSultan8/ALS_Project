# ==============================================================================
# 03_Machine_Learning_Consensus.R
#
# Triple-consensus machine learning feature selection on GSE112676 discovery.
# Three independent algorithms: Random Forest, LASSO, SVM-RFE.
# Only genes selected by ALL THREE are retained as the consensus signature.
#
# Key methodological decisions:
#   - Feature pool: union of WGCNA hub candidates AND significant DEGs
#   - SVM-RFE tests every integer from 1-30, then steps of 10 above 30
#   - make.names() sanitization with full reverse-mapping to HGNC symbols
#   - Fallback to RF+LASSO if 3-way intersection is empty
#   - Logistic regression trained on consensus, saved for blind validation
#   - Discovery AUC with 2000-bootstrap CI
#
# Inputs:
#   Processed_Data/Step1_Discovery_Data.rds
#   Processed_Data/Step3_WGCNA_Data.rds
#   Processed_Data/Step2_DEG_sig.csv
#   Processed_Data/Candidate_Hubs_for_ML.csv
#
# Outputs:
#   Processed_Data/Step4_ML_Results.rds
#   Processed_Data/Final_Gene_Signature.csv
#   Manuscript_Figures/ML/Figure5A_RF_OOB_Error.{pdf,png}
#   Manuscript_Figures/ML/Figure5B_RF_Importance.{pdf,png}
#   Manuscript_Figures/ML/Figure5C_LASSO_Path.{pdf,png}
#   Manuscript_Figures/ML/Figure5D_LASSO_CV.{pdf,png}
#   Manuscript_Figures/ML/Figure5E_SVM_RFE.{pdf,png}
#   Manuscript_Figures/ML/Figure5F_Discovery_ROC.{pdf,png}
#   Manuscript_Figures/ML/Figure5G_RiskScore_Violin.{pdf,png}
#   Manuscript_Figures/ML/Figure5H_Venn_Consensus.{pdf,png}
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(randomForest)
library(glmnet)
library(caret)
library(e1071)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggVennDiagram)
library(ggrepel)
library(ggpubr)   # for stat_compare_means in violin
library(scales)   # for number_format in 5E and sec_axis in 5D

dir.create("Manuscript_Figures/ML", recursive = TRUE, showWarnings = FALSE)
dir.create("Processed_Data",        showWarnings = FALSE)

# ==============================================================================
# STEP 1: Load Data
# ==============================================================================
message("Loading discovery data and WGCNA results...")

step1_data   <- readRDS("Processed_Data/Step1_Discovery_Data.rds")
wgcna_data   <- readRDS("Processed_Data/Step3_WGCNA_Data.rds")

expr_mat     <- step1_data$expr_discovery
group_factor <- step1_data$group_factor

if (is.null(names(group_factor))) {
  names(group_factor) <- colnames(expr_mat)
}

diagnosis_bin <- ifelse(group_factor == "ALS", 1, 0)

# ==============================================================================
# STEP 2: Build Feature Pool
# ==============================================================================
hub_candidates <- read.csv("Processed_Data/Candidate_Hubs_for_ML.csv")$Gene
deg_genes      <- read.csv("Processed_Data/Step2_DEG_sig.csv")$gene

candidates <- unique(union(
  intersect(hub_candidates, rownames(expr_mat)),
  intersect(deg_genes,      rownames(expr_mat))
))

message(sprintf(
  "Feature pool: %d genes | Hub candidates: %d | DEGs: %d | Overlap: %d",
  length(candidates),
  sum(hub_candidates %in% rownames(expr_mat)),
  sum(deg_genes      %in% rownames(expr_mat)),
  sum(hub_candidates %in% deg_genes)
))

# ==============================================================================
# STEP 3: Build Feature Matrix
# ==============================================================================
feat_raw    <- as.data.frame(t(expr_mat[candidates, ]))
feat_scaled <- as.data.frame(scale(feat_raw))

original_names        <- colnames(feat_scaled)
clean_names           <- make.names(original_names, unique = TRUE)
colnames(feat_scaled) <- clean_names
colnames(feat_raw)    <- clean_names

name_map <- setNames(original_names, clean_names)

feat_scaled$Diagnosis <- factor(
  ifelse(group_factor == "ALS", "ALS", "Control")
)

message(sprintf("Feature matrix: %d samples x %d genes",
                nrow(feat_scaled), length(candidates)))

# ==============================================================================
# STEP 4: Algorithm 1 -- Random Forest
# ==============================================================================
message("Running Random Forest (ntree = 500)...")

set.seed(1122)
rf_model <- randomForest(
  Diagnosis ~ .,
  data       = feat_scaled,
  ntree      = 500,
  importance = TRUE
)

message(sprintf("Random Forest OOB error: %.2f%%",
                rf_model$err.rate[500, "OOB"] * 100))

rf_importance <- as.data.frame(importance(rf_model))
rf_importance$CleanName <- rownames(rf_importance)
rf_importance <- rf_importance %>%
  arrange(desc(MeanDecreaseAccuracy))

rf_genes_clean <- rf_importance$CleanName[
  1:round(nrow(rf_importance) * 0.5)
]

message(sprintf("Random Forest: %d genes selected (top 50%% by MDA).",
                length(rf_genes_clean)))

# -- Figure 5A: OOB Error Convergence -----------------------------------------
# AESTHETIC CHANGE: light grid background added
oob_df <- data.frame(
  Trees   = seq_len(nrow(rf_model$err.rate)),
  OOB     = rf_model$err.rate[, "OOB"],
  ALS     = rf_model$err.rate[, "ALS"],
  Control = rf_model$err.rate[, "Control"]
) %>%
  pivot_longer(-Trees, names_to = "Error_Type", values_to = "Error_Rate")

p5A <- ggplot(oob_df, aes(x = Trees, y = Error_Rate * 100,
                            color = Error_Type)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("OOB" = "black",
               "ALS" = "#BC3C29FF",
               "Control" = "#0072B5FF"),
    name = "Error Type"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5),
        legend.position  = "right",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2)) +
  labs(title = "5A. Random Forest Out-of-Bag Error Rate",
       x     = "Number of Trees",
       y     = "Error Rate (%)")

ggsave("Manuscript_Figures/ML/Figure5A_RF_OOB_Error.pdf",
       plot = p5A, width = 7, height = 5)
ggsave("Manuscript_Figures/ML/Figure5A_RF_OOB_Error.png",
       plot = p5A, width = 7, height = 5, dpi = 600)
message("Figure 5A saved.")

# -- Figure 5B: RF Importance Dot Plot ----------------------------------------
# AESTHETIC CHANGE: size aesthetic added + color AND size legends shown + grid
top_rf_plot <- rf_importance %>%
  head(20) %>%
  mutate(Gene = name_map[CleanName])

p5B <- ggplot(top_rf_plot,
               aes(x     = MeanDecreaseAccuracy,
                   y     = reorder(Gene, MeanDecreaseAccuracy),
                   color = MeanDecreaseAccuracy,
                   size  = MeanDecreaseAccuracy)) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy,
                   y = Gene, yend = Gene),
               linewidth = 0.6, color = "grey50") +
  geom_point() +
  scale_color_gradient(low  = "#E18727FF", high = "#BC3C29FF",
                       name = "MDA",
                       guide = guide_colorbar(order = 1)) +
  scale_size_continuous(range = c(3, 8),
                        name  = "MDA",
                        guide = guide_legend(order = 2)) +
  theme_classic(base_size = 13) +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5),
        axis.text.y     = element_text(face = "bold.italic",
                                       color = "black", size = 11),
        axis.text.x     = element_text(color = "black"),
        legend.position = "right",
        panel.grid      = element_blank()) +
  labs(title = "5B. Top Genes by Random Forest Importance",
       x     = "Mean Decrease in Accuracy",
       y     = "")


ggsave("Manuscript_Figures/ML/Figure5B_RF_Importance.pdf",
       plot = p5B, width = 7, height = 7)
ggsave("Manuscript_Figures/ML/Figure5B_RF_Importance.png",
       plot = p5B, width = 7, height = 7, dpi = 600)
message("Figure 5B saved.")

# ==============================================================================
# STEP 5: Algorithm 2 -- LASSO Regression
# ==============================================================================
message("Running LASSO regression (10-fold CV)...")

set.seed(1122)
y_bin    <- as.numeric(feat_scaled$Diagnosis == "ALS")
x_mat    <- as.matrix(feat_scaled[, clean_names])

lasso_cv <- cv.glmnet(
  x_mat, y_bin,
  family = "binomial",
  alpha  = 1,
  nfolds = 10
)

lasso_coef        <- coef(lasso_cv, s = "lambda.min")
lasso_genes_clean <- rownames(lasso_coef)[
  which(abs(as.vector(lasso_coef)) > 0 &
        rownames(lasso_coef) != "(Intercept)")
]

message(sprintf("LASSO (lambda.min): %d genes selected.",
                length(lasso_genes_clean)))
message(sprintf("LASSO (lambda.1se): %d genes (reported for reference).",
                sum(abs(as.vector(
                  coef(lasso_cv, s = "lambda.1se"))) > 0) - 1))

# -- Figure 5C: LASSO Coefficient Path ----------------------------------------
# AESTHETIC CHANGE: top 15 genes colored uniquely, others in steelblue,
#                   lambda annotations with gene count, subtitle added
lasso_full <- glmnet(x_mat, y_bin, family = "binomial", alpha = 1)
coef_mat   <- as.matrix(coef(lasso_full))[-1, ]
coef_df    <- as.data.frame(t(coef_mat))
coef_df$LogLambda <- log(lasso_full$lambda)

coef_long <- coef_df %>%
  pivot_longer(-LogLambda,
               names_to  = "Gene",
               values_to = "Coefficient")

top15_genes <- lasso_genes_clean[
  order(abs(lasso_coef[lasso_genes_clean, 1]),
        decreasing = TRUE)
][1:min(15, length(lasso_genes_clean))]

top15_colors <- setNames(
  scales::hue_pal()(length(top15_genes)),
  top15_genes
)

coef_long <- coef_long %>%
  mutate(GeneGroup = ifelse(Gene %in% top15_genes, Gene, "other"))

p5C <- ggplot() +
  geom_line(data    = coef_long[coef_long$GeneGroup == "other", ],
            aes(x   = LogLambda, y = Coefficient, group = Gene),
            alpha   = 0.25, linewidth = 0.35, color = "steelblue") +
  geom_line(data    = coef_long[coef_long$GeneGroup != "other", ],
            aes(x   = LogLambda, y = Coefficient,
                group = Gene, color = GeneGroup),
            linewidth = 0.9) +
  scale_color_manual(values = top15_colors, name = "Top 15 Genes") +
  geom_vline(xintercept = log(lasso_cv$lambda.min),
             linetype = "dashed", color = "#0072B5FF", linewidth = 1) +
  geom_vline(xintercept = log(lasso_cv$lambda.1se),
             linetype = "dotted", color = "#20854EFF", linewidth = 1) +
  annotate("text",
           x     = log(lasso_cv$lambda.min),
           y     = max(abs(coef_long$Coefficient), na.rm = TRUE) * 0.95,
           label = expression(lambda[min]),
           hjust = -0.2, size = 4,
           color = "#0072B5FF", fontface = "bold") +
  annotate("text",
           x     = log(lasso_cv$lambda.1se),
           y     = max(abs(coef_long$Coefficient), na.rm = TRUE) * 0.75,
           label = expression(lambda["1SE"]),
           hjust = -0.2, size = 4,
           color = "#20854EFF", fontface = "bold") +
  theme_classic(base_size = 14) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 11),
        legend.position  = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2)) +
  labs(title    = "5C. LASSO Coefficient Regularization Path",
       subtitle = sprintf("lambda.min retains %d features",
                          length(lasso_genes_clean)),
       x        = "Log(Lambda)  [decreasing penalty ->]",
       y        = "Standardized Coefficient")


ggsave("Manuscript_Figures/ML/Figure5C_LASSO_Path.pdf",
       plot = p5C, width = 7, height = 5)
ggsave("Manuscript_Figures/ML/Figure5C_LASSO_Path.png",
       plot = p5C, width = 7, height = 5, dpi = 600)
message("Figure 5C saved.")

# -- Figure 5D: LASSO Cross-Validation ----------------------------------------
# AESTHETIC CHANGE: theme_bw with grid, scale_x_reverse, sec_axis for
#                   non-zero coefficient counts, expression() lambda labels,
#                   grey errorbars, green vline for lambda.1se
cv_df <- data.frame(
  LogLambda = log(lasso_cv$lambda),
  CVM       = lasso_cv$cvm,
  CVUP      = lasso_cv$cvup,
  CVLO      = lasso_cv$cvlo,
  NZero     = lasso_cv$nzero
)

n_top_axis <- cv_df[seq(1, nrow(cv_df), length.out = 8), ]

p5D <- ggplot(cv_df, aes(x = LogLambda, y = CVM)) +
  geom_errorbar(aes(ymin = CVLO, ymax = CVUP),
                width = 0.04, color = "grey60", alpha = 0.8) +
  geom_point(color = "#BC3C29FF", size = 2.5) +
  geom_vline(xintercept = log(lasso_cv$lambda.min),
             linetype = "dashed", color = "#0072B5FF", linewidth = 1) +
  geom_vline(xintercept = log(lasso_cv$lambda.1se),
             linetype = "dotted", color = "#20854EFF", linewidth = 1) +
  annotate("text",
           x     = log(lasso_cv$lambda.min),
           y     = max(cv_df$CVUP) * 1.02,
           label = expression(lambda[min]),
           hjust = 1.2, size = 4,
           color = "#0072B5FF", fontface = "bold") +
  annotate("text",
           x     = log(lasso_cv$lambda.1se),
           y     = max(cv_df$CVUP) * 1.02,
           label = expression(lambda["1SE"]),
           hjust = -0.2, size = 4,
           color = "#20854EFF", fontface = "bold") +
  scale_x_reverse(
    sec.axis = sec_axis(
      ~ .,
      breaks = n_top_axis$LogLambda,
      labels = n_top_axis$NZero,
      name   = "Number of Non-Zero Coefficients"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
        axis.text        = element_text(color = "black")) +
  labs(title = "5D. LASSO 10-Fold Cross-Validation",
       x     = "Log(Lambda)",
       y     = "Binomial Deviance")

ggsave("Manuscript_Figures/ML/Figure5D_LASSO_CV.pdf",
       plot = p5D, width = 7, height = 5)
ggsave("Manuscript_Figures/ML/Figure5D_LASSO_CV.png",
       plot = p5D, width = 7, height = 5, dpi = 600)
message("Figure 5D saved.")

# ==============================================================================
# STEP 6: Algorithm 3 -- SVM-RFE
# ==============================================================================
message("Running SVM-RFE (5-fold CV, fine-grained 1-30)...")
message("This step typically takes 15-25 minutes...")

ctrl_svm <- rfeControl(
  functions    = caretFuncs,
  method       = "cv",
  number       = 5,
  verbose      = FALSE,
  returnResamp = "final"
)

train_ctrl <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter     = FALSE
)

subsets <- c(1:30, seq(35, min(100, length(candidates)), by = 10))

set.seed(1122)
svm_rfe <- rfe(
  x          = feat_scaled[, clean_names],
  y          = feat_scaled$Diagnosis,
  sizes      = subsets,
  rfeControl = ctrl_svm,
  method     = "svmLinear",
  metric     = "ROC",
  trControl  = train_ctrl
)

svm_perf <- svm_rfe$results
perf_col <- if ("ROC" %in% colnames(svm_perf)) {
  "ROC"
} else if ("Accuracy" %in% colnames(svm_perf)) {
  "Accuracy"
} else {
  num_cols <- colnames(svm_perf)[
    sapply(svm_perf, is.numeric) & colnames(svm_perf) != "Variables"
  ]
  num_cols[1]
}

message(sprintf("SVM-RFE performance metric: %s", perf_col))
message("SVM-RFE results (1-30 range):")
print(svm_perf[svm_perf$Variables <= 30, c("Variables", perf_col)])

svm_perf30 <- svm_perf[svm_perf$Variables <= 30, ]
optimal_n  <- svm_perf30$Variables[which.max(svm_perf30[[perf_col]])]

svm_genes_clean <- predictors(svm_rfe)[
  1:min(optimal_n, length(predictors(svm_rfe)))
]

message(sprintf("SVM-RFE: optimal subset = %d genes.", optimal_n))

# -- Figure 5E: SVM-RFE Performance -------------------------------------------
# AESTHETIC CHANGE: grid background, larger points, number_format on y-axis
plot_df <- svm_perf[svm_perf$Variables <= 30, ]

p5E <- ggplot(plot_df, aes(x = Variables, y = .data[[perf_col]])) +
  geom_line(color = "#0072B5FF", linewidth = 1) +
  geom_point(color = "#0072B5FF", size = 4) +
  geom_point(
    data  = plot_df[plot_df$Variables == optimal_n, ],
    aes(x = Variables, y = .data[[perf_col]]),
    color = "#BC3C29FF", size = 6
  ) +
  annotate("text",
           x     = optimal_n,
           y     = plot_df[[perf_col]][plot_df$Variables == optimal_n] + 0.005,
           label = paste0("Optimal: ", optimal_n, " features"),
           hjust = -0.15, vjust = 0,
           color = "#BC3C29FF", fontface = "bold", size = 4.5) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  theme_classic(base_size = 14) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 11),
        axis.text        = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2)) +
  labs(title    = "5E. SVM-RFE Cross-Validation Performance",
       subtitle = "Feature subset sizes 1-30 (single-gene resolution)",
       x        = "Number of Selected Features",
       y        = paste("Cross-Validation", perf_col))


ggsave("Manuscript_Figures/ML/Figure5E_SVM_RFE.pdf",
       plot = p5E, width = 7, height = 5)
ggsave("Manuscript_Figures/ML/Figure5E_SVM_RFE.png",
       plot = p5E, width = 7, height = 5, dpi = 600)
message("Figure 5E saved.")

# ==============================================================================
# STEP 7: Triple-Algorithm Consensus
# ==============================================================================
message("Computing triple-algorithm consensus...")
message(sprintf("  Random Forest  : %d genes", length(rf_genes_clean)))
message(sprintf("  LASSO          : %d genes", length(lasso_genes_clean)))
message(sprintf("  SVM-RFE        : %d genes", length(svm_genes_clean)))

consensus_clean <- Reduce(intersect,
                           list(rf_genes_clean,
                                lasso_genes_clean,
                                svm_genes_clean))

message(sprintf("  Three-way intersection: %d genes", length(consensus_clean)))

fallback_used <- FALSE
if (length(consensus_clean) == 0) {
  message("WARNING: Three-way empty. Falling back to RF + LASSO.")
  consensus_clean <- intersect(rf_genes_clean, lasso_genes_clean)
  fallback_used   <- TRUE
  message(sprintf("  RF + LASSO intersection: %d genes", length(consensus_clean)))
}

if (length(consensus_clean) == 0) {
  stop("No consensus genes. Review feature pool and algorithm parameters.")
}

consensus_genes <- unique(name_map[consensus_clean])
consensus_genes <- consensus_genes[!is.na(consensus_genes)]

message(sprintf("Final consensus signature: %d genes", length(consensus_genes)))
message(sprintf("Genes: %s", paste(consensus_genes, collapse = ", ")))

original_6 <- c("XPO1", "PPP1CC", "DHX15", "EIF4G2", "ZFAND5", "FBXO11")
recovered  <- intersect(original_6, consensus_genes)
message(sprintf("Original 6 recovered: %d/6 -- %s",
                length(recovered),
                ifelse(length(recovered) > 0,
                       paste(recovered, collapse = ", "),
                       "none")))

# -- Figure 5H: Three-Way Venn Diagram ----------------------------------------
# AESTHETIC CHANGE: legend.position = "right" added; renamed to 5H
venn_list <- list(
  `Random Forest` = name_map[rf_genes_clean],
  `LASSO`         = name_map[lasso_genes_clean],
  `SVM-RFE`       = name_map[svm_genes_clean]
)

p5H <- ggVennDiagram(venn_list, label_alpha = 0, edge_size = 1) +
  scale_fill_gradient(low = "#FDDBC7", high = "#BC3C29FF") +
  scale_color_manual(values = rep("black", 3)) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14)) +
  labs(title = "5H. Triple-Algorithm Consensus")

ggsave("Manuscript_Figures/ML/Figure5H_Venn_Consensus.pdf",
       plot = p5H, width = 7, height = 6)
ggsave("Manuscript_Figures/ML/Figure5H_Venn_Consensus.png",
       plot = p5H, width = 7, height = 6, dpi = 600)
message("Figure 5H saved.")

# ==============================================================================
# STEP 8: Logistic Regression Model and Discovery ROC
# ==============================================================================
message("Training logistic regression model on consensus signature...")

missing_genes <- setdiff(consensus_genes, rownames(expr_mat))
if (length(missing_genes) > 0) {
  message(sprintf("WARNING: %d genes missing from expression matrix.",
                  length(missing_genes)))
  consensus_genes <- intersect(consensus_genes, rownames(expr_mat))
}

disc_data_raw       <- as.data.frame(
  t(expr_mat[consensus_genes, , drop = FALSE])
)
disc_data           <- as.data.frame(scale(disc_data_raw))
colnames(disc_data) <- make.names(colnames(disc_data), unique = TRUE)
model_gene_names    <- colnames(disc_data)
disc_data$Diagnosis <- diagnosis_bin

logit_model         <- glm(Diagnosis ~ ., data = disc_data,
                            family = "binomial")
disc_data$RiskScore <- predict(logit_model, type = "response")

# Group factor needed for violin plot
disc_data$Group <- factor(
  ifelse(disc_data$Diagnosis == 1, "ALS", "Control"),
  levels = c("Control", "ALS")
)

set.seed(1122)
roc_disc <- roc(
  disc_data$Diagnosis,
  disc_data$RiskScore,
  ci        = TRUE,
  ci.method = "bootstrap",
  boot.n    = 2000,
  quiet     = TRUE
)

message(sprintf("Discovery AUC: %.3f (95%% CI: %.3f - %.3f)",
                as.numeric(roc_disc$auc),
                as.numeric(roc_disc$ci[1]),
                as.numeric(roc_disc$ci[3])))

message(sprintf("Risk score -- ALS median: %.3f | CON median: %.3f",
                median(disc_data$RiskScore[diagnosis_bin == 1]),
                median(disc_data$RiskScore[diagnosis_bin == 0])))

# -- Figure 5F: Discovery ROC -------------------------------------------------
# AESTHETIC CHANGE: theme_bw with light grid background
roc_df <- data.frame(
  Specificity = roc_disc$specificities,
  Sensitivity = roc_disc$sensitivities
)

auc_label <- sprintf("AUC = %.3f\n95%% CI: %.3f-%.3f",
                      as.numeric(roc_disc$auc),
                      as.numeric(roc_disc$ci[1]),
                      as.numeric(roc_disc$ci[3]))

p5F <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity)) +
  geom_line(color = "#BC3C29FF", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.65, y = 0.15,
           label = auc_label, size = 5,
           color = "#BC3C29FF", fontface = "bold") +
  theme_bw(base_size = 14) +
  theme(plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, color = "grey40"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2)) +
  labs(
    title    = sprintf("5F. Discovery Cohort ROC (GSE112676, n = %d)",
                       nrow(disc_data)),
    subtitle = sprintf("Consensus signature: %d genes",
                       length(consensus_genes)),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  )

ggsave("Manuscript_Figures/ML/Figure5F_Discovery_ROC.pdf",
       plot = p5F, width = 6, height = 6)
ggsave("Manuscript_Figures/ML/Figure5F_Discovery_ROC.png",
       plot = p5F, width = 6, height = 6, dpi = 600)
message("Figure 5F saved.")

# -- Figure 5G: Diagnostic Risk Score Violin ----------------------------------
# NEW: added from reference script, matching WGCNA violin aesthetics
p5G <- ggplot(disc_data,
               aes(x = Group, y = RiskScore, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.1, fill = "white",
               color = "black", outlier.shape = NA) +
  scale_fill_manual(values = c("ALS" = "#BC3C29FF", "Control" = "#0072B5FF")) +
  stat_compare_means(method   = "wilcox.test",
                     label    = "p.format",
                     size     = 5,
                     fontface = "bold",
                     label.y  = max(disc_data$RiskScore) * 1.15) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none",
        plot.title      = element_text(face = "bold", hjust = 0.5),
        axis.text       = element_text(color = "black")) +
  labs(title = "5G. Diagnostic Risk Score Distribution (Discovery Cohort)",
       x     = "",
       y     = "Predicted Probability of ALS")

ggsave("Manuscript_Figures/ML/Figure5G_RiskScore_Violin.pdf",
       plot = p5G, width = 5, height = 6)
ggsave("Manuscript_Figures/ML/Figure5G_RiskScore_Violin.png",
       plot = p5G, width = 5, height = 6, dpi = 600)
message("Figure 5G (Violin) saved.")

# ==============================================================================
# STEP 9: Per-Gene Individual AUC in Discovery
# ==============================================================================
message("Computing per-gene individual AUC in discovery cohort...")

auc_disc <- sapply(consensus_genes, function(g) {
  if (!g %in% rownames(expr_mat)) return(NA)
  tryCatch({
    r <- roc(as.numeric(group_factor == "ALS"),
             as.numeric(expr_mat[g, ]),
             quiet = TRUE)
    as.numeric(r$auc)
  }, error = function(e) NA)
})

auc_disc_df <- data.frame(
  Gene          = consensus_genes,
  AUC_Discovery = round(auc_disc, 3)
) %>% arrange(desc(AUC_Discovery))

message("Individual AUC in discovery:")
print(auc_disc_df)

write.csv(auc_disc_df,
          "Processed_Data/Step4_PerGene_AUC_Discovery.csv",
          row.names = FALSE)

# ==============================================================================
# STEP 10: Save All Results
# ==============================================================================
write.csv(
  data.frame(Gene = consensus_genes),
  "Processed_Data/Final_Gene_Signature.csv",
  row.names = FALSE
)

saveRDS(
  list(
    consensus_genes  = consensus_genes,
    rf_genes         = name_map[rf_genes_clean],
    lasso_genes      = name_map[lasso_genes_clean],
    svm_genes        = name_map[svm_genes_clean],
    rf_importance    = rf_importance,
    lasso_cv         = lasso_cv,
    logit_model      = logit_model,
    model_gene_names = model_gene_names,
    roc_disc         = roc_disc,
    disc_data        = disc_data,
    auc_disc_df      = auc_disc_df,
    fallback_used    = fallback_used,
    optimal_svm_n    = optimal_n,
    name_map         = name_map
  ),
  "Processed_Data/Step4_ML_Results.rds"
)

message("================================================================")
message("Script 03 complete.")
message(sprintf("  Feature pool          : %d genes", length(candidates)))
message(sprintf("  Random Forest         : %d genes (top 50%% MDA)",
                length(rf_genes_clean)))
message(sprintf("  LASSO (lambda.min)    : %d genes",
                length(lasso_genes_clean)))
message(sprintf("  SVM-RFE               : %d genes (optimal n = %d)",
                length(svm_genes_clean), optimal_n))
message(sprintf("  Consensus signature   : %d genes",
                length(consensus_genes)))
message(sprintf("  Fallback used         : %s",
                ifelse(fallback_used, "YES (RF+LASSO)", "NO (3-way)")))
message(sprintf("  Discovery AUC         : %.3f (CI: %.3f-%.3f)",
                as.numeric(roc_disc$auc),
                as.numeric(roc_disc$ci[1]),
                as.numeric(roc_disc$ci[3])))
message(sprintf("  Consensus genes       : %s",
                paste(consensus_genes, collapse = ", ")))
message("================================================================")
