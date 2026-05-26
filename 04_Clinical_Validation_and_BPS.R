# ==============================================================================
# 04_Clinical_Validation_and_BPS.R
#
# Five-arm independent validation of the consensus diagnostic signature.
# The trained logistic regression model from Script 03 is applied without
# retraining to each validation dataset.
#
# VALIDATION ARMS:
#   Arm 1 — GSE112680 ALS vs CON  : primary replication (Illumina V4.0, n~301)
#   Arm 2 — GSE112680 ALS vs MIM  : clinical specificity — novel finding
#   Arm 3 — GSE112680 MIM vs CON  : specificity control (expected AUC ~0.5)
#   Arm 4 — GSE28253  ALS vs CON  : independent microarray (Agilent GPL4133)
#   Arm 5 — GSE234297 ALS vs CON  : independent RNA-seq cohort
#
# BPS COMPUTATION:
#   Per-gene AUC: Arms 1, 4, 5 + discovery (diagnostic arms only).
#   Arms 2 and 3 address clinical specificity — excluded from BPS by design.
#   kME from each gene's own module (not a single reference module).
#
# Inputs:
#   Processed_Data/Step4_ML_Results.rds
#   Processed_Data/Step1_Discovery_Data.rds
#   Processed_Data/Step3_WGCNA_Data.rds
#   GSE112680_series_matrix.txt.gz (local cache preferred)
#   GSE28253_series_matrix.txt.gz  (local cache preferred)
#   GSE234297 (downloaded via GEOquery + supplementary counts)
#
# Outputs:
#   Processed_Data/Step5_Validation_Results.rds
#   Processed_Data/Step5_PerGene_AUCs.csv
#   Processed_Data/Step5_Validation_Summary.csv
#   Manuscript_Figures/Validation/Figure6A_*.{pdf,png}
#   Manuscript_Figures/Validation/Figure6B_*.{pdf,png}
#   Manuscript_Figures/Validation/Figure6C_RiskScore_Violin.{pdf,png}
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(GEOquery)
library(dplyr)
library(tidyr)
library(pROC)
library(ggplot2)
library(ggpubr)
library(illuminaHumanv4.db)
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("Manuscript_Figures/Validation", recursive = TRUE,
           showWarnings = FALSE)
dir.create("Processed_Data", showWarnings = FALSE)

# ==============================================================================
# STEP 1: Load Trained Model, Consensus Signature, Discovery Data
# ==============================================================================
message("Loading trained model and consensus signature...")

ml_results       <- readRDS("Processed_Data/Step4_ML_Results.rds")
step1_data       <- readRDS("Processed_Data/Step1_Discovery_Data.rds")
wgcna_data       <- readRDS("Processed_Data/Step3_WGCNA_Data.rds")

logit_model      <- ml_results$logit_model
consensus_genes  <- ml_results$consensus_genes
model_gene_names <- ml_results$model_gene_names
disc_center      <- ml_results$disc_center
disc_scale_sd    <- ml_results$disc_scale_sd
disc_expr        <- step1_data$expr_discovery
disc_labels      <- step1_data$group_factor

if (is.null(names(disc_labels))) {
  names(disc_labels) <- colnames(disc_expr)
}

message(sprintf("Consensus signature: %d genes — %s",
                length(consensus_genes),
                paste(consensus_genes, collapse = ", ")))

# ==============================================================================
# STEP 2: Helper Functions
# ==============================================================================

# Prepare validation data for model prediction.
# Z-score scales INDEPENDENTLY of discovery — correct for blind validation.
# Missing genes imputed as 0 after scaling.
prepare_val_data <- function(val_expr, val_labels, sig_genes,
                              model, model_names,
                              disc_center, disc_scale_sd) {
  available <- intersect(sig_genes, rownames(val_expr))
  missing   <- setdiff(sig_genes, rownames(val_expr))

  if (length(available) == 0) {
    warning("No signature genes found in dataset.")
    return(NULL)
  }

  if (length(missing) > 0) {
    message(sprintf("  Missing genes (imputed as 0): %s",
                    paste(missing, collapse = ", ")))
  }

  val_df <- as.data.frame(t(val_expr[available, , drop = FALSE]))
  for (g in available) {
    if (g %in% names(disc_center) && disc_scale_sd[g] > 0) {
	  val_df[[g]] <- (val_df[[g]] - disc_center[g]) / disc_scale_sd[g]
	}
  }

  for (g in missing) val_df[[g]] <- 0

  colnames(val_df) <- make.names(colnames(val_df), unique = TRUE)

  for (nm in model_names) {
    if (!nm %in% colnames(val_df)) val_df[[nm]] <- 0
  }

  val_df <- val_df[, model_names, drop = FALSE]

  val_df$Diagnosis <- as.numeric(val_labels == "ALS")
  val_df$RiskScore <- predict(model, newdata = val_df, type = "response")
  val_df$Group     <- val_labels

  return(val_df)
}

# Compute ROC with 95% CI
compute_roc <- function(val_df) {
  roc(val_df$Diagnosis, val_df$RiskScore,
      ci = TRUE, quiet = TRUE)
}

# Save a single ROC curve plot
save_roc_plot <- function(roc_obj, title_str, subtitle_str,
                           file_prefix, line_color = "#0072B5FF") {
  roc_df <- data.frame(
    Spec = roc_obj$specificities,
    Sens = roc_obj$sensitivities
  )

  auc_label <- sprintf("AUC = %.3f\n95%% CI: %.3f-%.3f",
                        as.numeric(roc_obj$auc),
                        as.numeric(roc_obj$ci[1]),
                        as.numeric(roc_obj$ci[3]))

  p <- ggplot(roc_df, aes(x = 1 - Spec, y = Sens)) +
    geom_line(color = line_color, linewidth = 1.3) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_ribbon(aes(ymin = Sens - 0.015, ymax = Sens + 0.015),
                alpha = 0.1, fill = line_color) +
    annotate("text", x = 0.6, y = 0.18,
             label = auc_label, size = 4.8,
             color = line_color, fontface = "bold") +
    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 11),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2)
    ) +
    labs(title    = title_str,
         subtitle = subtitle_str,
         x = "1 - Specificity (False Positive Rate)",
         y = "Sensitivity (True Positive Rate)")


  ggsave(paste0(file_prefix, ".pdf"), plot = p, width = 6, height = 6)
  ggsave(paste0(file_prefix, ".png"), plot = p, width = 6, height = 6,
         dpi = 600)
  return(invisible(p))
}

# Build combined single-gene ROC panel for one cohort
build_gene_roc_panel <- function(val_expr, labels_vec,
                                  sig_genes, cohort_name, file_prefix) {
  clean_labels <- as.numeric(labels_vec == "ALS")
  gene_colors  <- c("#BC3C29FF", "#E18727FF", "#20854EFF",
                     "#0072B5FF", "#7876B1FF", "#6F99ADFF")
  roc_rows <- list()

  for (i in seq_along(sig_genes)) {
    g <- sig_genes[i]
    if (!g %in% rownames(val_expr)) next
    r <- tryCatch(
      roc(clean_labels, as.numeric(val_expr[g, ]), quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(r)) next
    auc_val <- round(as.numeric(r$auc), 3)
    roc_rows[[g]] <- data.frame(
      Gene        = sprintf("%s (AUC=%.3f)", g, auc_val),
      Specificity = r$specificities,
      Sensitivity = r$sensitivities,
      Color       = gene_colors[i],
      stringsAsFactors = FALSE
    )
  }

  if (length(roc_rows) == 0) {
    message(sprintf("No genes found for %s — skipping.", cohort_name))
    return(NULL)
  }

  roc_df    <- do.call(rbind, roc_rows)
  color_map <- unique(roc_df[, c("Gene", "Color")])
  cv        <- setNames(color_map$Color, color_map$Gene)

  p <- ggplot(roc_df, aes(x = 1 - Specificity,
                            y = Sensitivity, color = Gene)) +
    geom_line(linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_color_manual(values = cv) +
    theme_bw(base_size = 12) +
    theme(
      legend.position   = c(0.62, 0.28),
      legend.text       = element_text(size = 8, face = "italic"),
      legend.title      = element_blank(),
      legend.background = element_rect(fill = alpha("white", 0.8)),
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle     = element_text(hjust = 0.5, color = "grey40", size = 10),
      panel.grid.major  = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor  = element_line(color = "grey95", linewidth = 0.2)
    ) +
    labs(
      title    = sprintf("6B. Individual Gene ROC -- %s", cohort_name),
      subtitle = "Each gene's independent discriminative power",
      x = "1 - Specificity",
      y = "Sensitivity"
    )


  ggsave(paste0(file_prefix, ".pdf"), plot = p, width = 7, height = 6)
  ggsave(paste0(file_prefix, ".png"), plot = p, width = 7, height = 6,
         dpi = 600)
  message(sprintf("Figure 6B (%s) saved.", cohort_name))
  return(invisible(p))
}

# Compute single-gene AUC safely
get_gene_auc <- function(gene, expr_matrix, labels) {
  binary_labels <- as.numeric(labels == "ALS")
  if (!gene %in% rownames(expr_matrix)) return(NA_real_)
  if (length(unique(binary_labels)) < 2) return(NA_real_)
  tryCatch({
    r <- roc(binary_labels, as.numeric(expr_matrix[gene, ]), quiet = TRUE)
    as.numeric(r$auc)
  }, error = function(e) NA_real_)
}

# ==============================================================================
# STEP 3: Load and Process GSE112680
# Illumina HumanHT-12 V4.0 (GPL10558) — NOT V3
# All three groups: ALS, CON, MIM
# This dataset was never touched during discovery — zero leakage
# ==============================================================================
message("
==============================================================
  Loading GSE112680 (Validation Cohort 1)
==============================================================")

if (file.exists("GSE112680_series_matrix.txt.gz")) {
  message("Loading from local cache...")
  gse112680 <- getGEO(filename = "GSE112680_series_matrix.txt.gz",
                      getGPL   = FALSE)
} else {
  message("Downloading GSE112680...")
  gse112680 <- getGEO("GSE112680", getGPL = FALSE)[[1]]
}

expr_112680_raw  <- exprs(gse112680)
pheno_112680     <- pData(gse112680)
group_112680_raw <- toupper(trimws(
  as.character(pheno_112680$`diagnosis:ch1`)
))

keep_112680     <- group_112680_raw %in% c("ALS", "CON", "MIM")
expr_112680_raw <- expr_112680_raw[, keep_112680]
group_112680    <- group_112680_raw[keep_112680]

message(sprintf("GSE112680: ALS=%d | CON=%d | MIM=%d",
                sum(group_112680 == "ALS"),
                sum(group_112680 == "CON"),
                sum(group_112680 == "MIM")))

# Probe mapping — illuminaHumanv4.db (V4.0 chip)
message("Mapping probes using illuminaHumanv4.db...")

probes_v4  <- rownames(expr_112680_raw)
symbols_v4 <- mapIds(
  illuminaHumanv4.db,
  keys      = probes_v4,
  column    = "SYMBOL",
  keytype   = "PROBEID",
  multiVals = "first"
)

valid_v4       <- !is.na(symbols_v4) & symbols_v4 != ""
expr_v4_valid  <- expr_112680_raw[valid_v4, ]
syms_v4_valid  <- symbols_v4[valid_v4]

message(sprintf("GSE112680: %d / %d probes mapped.",
                sum(valid_v4), length(valid_v4)))

if (max(expr_v4_valid, na.rm = TRUE) > 100) {
  expr_v4_valid <- log2(expr_v4_valid + 1)
  message("GSE112680: log2 transformation applied.")
}

# Best probe by IQR — consistent with discovery preprocessing
probe_iqr_v4 <- apply(expr_v4_valid, 1, IQR, na.rm = TRUE)
probe_df_v4  <- data.frame(
  probe  = rownames(expr_v4_valid),
  symbol = syms_v4_valid,
  iqr    = probe_iqr_v4,
  stringsAsFactors = FALSE
)

best_probes_v4 <- probe_df_v4 %>%
  group_by(symbol) %>%
  slice_max(order_by = iqr, n = 1, with_ties = FALSE) %>%
  pull(probe)

val_expr_112680           <- expr_v4_valid[best_probes_v4, ]
rownames(val_expr_112680) <- probe_df_v4$symbol[
  match(best_probes_v4, probe_df_v4$probe)
]

message(sprintf("GSE112680 matrix: %d genes x %d samples",
                nrow(val_expr_112680), ncol(val_expr_112680)))
message(sprintf("Consensus gene coverage: %d / %d",
                sum(consensus_genes %in% rownames(val_expr_112680)),
                length(consensus_genes)))

# ── Arm 1: ALS vs CON ─────────────────────────────────────────────────────────
message("--- Arm 1: GSE112680 ALS vs CON ---")

keep_arm1   <- group_112680 %in% c("ALS", "CON")
labels_arm1 <- ifelse(group_112680[keep_arm1] == "ALS", "ALS", "Control")

val_df_arm1 <- prepare_val_data(
  val_expr_112680[, keep_arm1], labels_arm1,
  consensus_genes, logit_model, model_gene_names,
  disc_center, disc_scale_sd
)
roc_arm1 <- compute_roc(val_df_arm1)

message(sprintf("Arm 1: n=%d | AUC=%.3f (CI: %.3f-%.3f)",
                sum(keep_arm1),
                as.numeric(roc_arm1$auc),
                as.numeric(roc_arm1$ci[1]),
                as.numeric(roc_arm1$ci[3])))

save_roc_plot(
  roc_arm1,
  title_str    = "GSE112680 -- ALS vs Control (Primary Replication)",
  subtitle_str = sprintf("n=%d | Illumina V4.0 | Independent cohort",
                          sum(keep_arm1)),
  file_prefix  = "Manuscript_Figures/Validation/Figure6A_Arm1_GSE112680_ALS_CON",
  line_color   = "#0072B5FF"
)

# ── Arm 2: ALS vs MIM ─────────────────────────────────────────────────────────
message("--- Arm 2: GSE112680 ALS vs MIM ---")

keep_arm2          <- group_112680 %in% c("ALS", "MIM")
labels_arm2_binary <- ifelse(group_112680[keep_arm2] == "ALS",
                              "ALS", "Control")

val_df_arm2 <- prepare_val_data(
  val_expr_112680[, keep_arm2], labels_arm2_binary,
  consensus_genes, logit_model, model_gene_names,
  disc_center, disc_scale_sd
)
roc_arm2 <- compute_roc(val_df_arm2)

message(sprintf("Arm 2: n=%d | AUC=%.3f (CI: %.3f-%.3f)",
                sum(keep_arm2),
                as.numeric(roc_arm2$auc),
                as.numeric(roc_arm2$ci[1]),
                as.numeric(roc_arm2$ci[3])))

save_roc_plot(
  roc_arm2,
  title_str    = "GSE112680 -- ALS vs Mimic (Clinical Specificity)",
  subtitle_str = sprintf("n=%d (ALS=%d, MIM=%d) | Novel finding",
                          sum(keep_arm2),
                          sum(group_112680[keep_arm2] == "ALS"),
                          sum(group_112680[keep_arm2] == "MIM")),
  file_prefix  = "Manuscript_Figures/Validation/Figure6A_Arm2_GSE112680_ALS_MIM",
  line_color   = "#BC3C29FF"
)

# ── Arm 3: MIM vs CON ─────────────────────────────────────────────────────────
message("--- Arm 3: GSE112680 MIM vs CON ---")
message("    MIM coded as positive class for ROC directionality only.")
message("    Expected AUC near 0.5 if signature is ALS-specific.")

keep_arm3   <- group_112680 %in% c("MIM", "CON")
labels_arm3 <- ifelse(group_112680[keep_arm3] == "MIM",
                       "ALS", "Control")

val_df_arm3 <- prepare_val_data(
  val_expr_112680[, keep_arm3], labels_arm3,
  consensus_genes, logit_model, model_gene_names,
  disc_center, disc_scale_sd
)
roc_arm3 <- compute_roc(val_df_arm3)

message(sprintf("Arm 3: n=%d | AUC=%.3f (CI: %.3f-%.3f)",
                sum(keep_arm3),
                as.numeric(roc_arm3$auc),
                as.numeric(roc_arm3$ci[1]),
                as.numeric(roc_arm3$ci[3])))

save_roc_plot(
  roc_arm3,
  title_str    = "GSE112680 -- Mimic vs Control (Specificity Check)",
  subtitle_str = sprintf("n=%d (MIM=%d, CON=%d) | Expected AUC near 0.5",
                          sum(keep_arm3),
                          sum(group_112680[keep_arm3] == "MIM"),
                          sum(group_112680[keep_arm3] == "CON")),
  file_prefix  = "Manuscript_Figures/Validation/Figure6A_Arm3_GSE112680_MIM_CON",
  line_color   = "#7876B1FF"
)

# ── Three-Group Risk Score Violin — Figure 6C ──────────────────────────────────
message("Generating three-group risk score violin (Figure 6C)...")

all_112680_df <- prepare_val_data(
  val_expr_112680, group_112680,
  consensus_genes, logit_model, model_gene_names,
  disc_center, disc_scale_sd
)

all_112680_df$Group <- factor(group_112680, levels = c("CON", "MIM", "ALS"))

p_violin_112680 <- ggplot(
  all_112680_df,
  aes(x = Group, y = RiskScore, fill = Group)
) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, fill = "white",
               color = "black", outlier.shape = NA) +
  stat_compare_means(
    comparisons = list(c("CON", "ALS"), c("MIM", "ALS"), c("CON", "MIM")),
    method      = "wilcox.test",
    label       = "p.format",
    size        = 3.5
  ) +
  scale_fill_manual(values = c(
    "ALS" = "#BC3C29FF",
    "CON" = "#0072B5FF",
    "MIM" = "#E18727FF"
  )) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Diagnostic Risk Score Distribution (GSE112680)",
    subtitle = "ALS vs Control vs Mimic — three-group comparison",
    x = "", y = "Predicted ALS Probability"
  )

ggsave("Manuscript_Figures/Validation/Figure6C_RiskScore_Violin.pdf",
       plot = p_violin_112680, width = 7, height = 7)
ggsave("Manuscript_Figures/Validation/Figure6C_RiskScore_Violin.png",
       plot = p_violin_112680, width = 7, height = 7, dpi = 600)
message("Figure 6C saved.")

# ==============================================================================
# STEP 4: Load and Process GSE28253 (Agilent GPL4133)
# ==============================================================================
message("
==============================================================
  Loading GSE28253 (Validation Cohort 2)
==============================================================")

load_gse28253 <- function() {
  if (file.exists("GSE28253_series_matrix.txt.gz")) {
    message("Loading from local cache...")
    gse <- getGEO(filename = "GSE28253_series_matrix.txt.gz",
                  getGPL   = FALSE)
    fd  <- fData(gse)
    sym_col <- grep("^GENE_SYMBOL$", colnames(fd),
                     ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(sym_col) && nrow(fd) > 0) return(gse)
    message("No annotation in cache — downloading with GPL...")
  }
  tryCatch(
    getGEO("GSE28253", getGPL = TRUE)[[1]],
    error = function(e) {
      message("GPL download failed — series matrix only.")
      getGEO("GSE28253", getGPL = FALSE)[[1]]
    }
  )
}

gse28253    <- load_gse28253()
expr_28253  <- exprs(gse28253)
pdata_28253 <- pData(gse28253)
fd_28253    <- fData(gse28253)

# Identify gene symbol column
sym_col <- NA_character_
for (pattern in c("^GENE_SYMBOL$", "^GeneSymbol$",
                   "^gene_symbol$", "^Symbol$",
                   "GENE_SYMBOL", "GeneSymbol")) {
  matches <- grep(pattern, colnames(fd_28253),
                   ignore.case = FALSE, value = TRUE)
  if (length(matches) > 0) {
    sym_col <- matches[1]
    break
  }
}

if (is.na(sym_col)) {
  stop("GSE28253: Cannot identify gene symbol column.")
}

message(sprintf("GSE28253: Using symbol column '%s'", sym_col))

gene_sym_28253 <- sapply(
  strsplit(as.character(fd_28253[[sym_col]]), " |///|/"), `[`, 1
)

valid_28253    <- !is.na(gene_sym_28253) &
                  gene_sym_28253 != "" &
                  gene_sym_28253 != "---"

df_28253       <- as.data.frame(expr_28253[valid_28253, ])
df_28253$gene  <- gene_sym_28253[valid_28253]

val_expr_28253 <- df_28253 %>%
  group_by(gene) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame()
rownames(val_expr_28253) <- val_expr_28253$gene
val_expr_28253 <- as.matrix(val_expr_28253[, -1])

if (max(val_expr_28253, na.rm = TRUE) > 100) {
  val_expr_28253 <- log2(val_expr_28253 + 1)
  message("GSE28253: log2 transformation applied.")
}

# Label extraction
all_text_28253 <- paste(
  pdata_28253$title,
  pdata_28253$source_name_ch1,
  pdata_28253$characteristics_ch1
)

labels_28253 <- ifelse(
  grepl("ALS|Amyotrophic", all_text_28253, ignore.case = TRUE), "ALS",
  ifelse(grepl("Control|Healthy|Normal",
               all_text_28253, ignore.case = TRUE), "Control", NA)
)

keep_28253 <- !is.na(labels_28253)
message(sprintf("GSE28253: %d samples (ALS=%d | CON=%d)",
                sum(keep_28253),
                sum(labels_28253[keep_28253] == "ALS"),
                sum(labels_28253[keep_28253] == "Control")))
message(sprintf("Consensus gene coverage: %d / %d",
                sum(consensus_genes %in% rownames(val_expr_28253)),
                length(consensus_genes)))

val_df_28253 <- prepare_val_data(
  val_expr_28253[, keep_28253],
  labels_28253[keep_28253],
  consensus_genes, logit_model, model_gene_names,
  disc_center, disc_scale_sd
)
roc_28253 <- compute_roc(val_df_28253)

message(sprintf("Arm 4: n=%d | AUC=%.3f (CI: %.3f-%.3f)",
                sum(keep_28253),
                as.numeric(roc_28253$auc),
                as.numeric(roc_28253$ci[1]),
                as.numeric(roc_28253$ci[3])))

save_roc_plot(
  roc_28253,
  title_str    = "GSE28253 -- ALS vs Control (External Microarray)",
  subtitle_str = sprintf("n=%d | Agilent GPL4133 | Independent cohort",
                          sum(keep_28253)),
  file_prefix  = "Manuscript_Figures/Validation/Figure6A_Arm4_GSE28253",
  line_color   = "#20854EFF"
)

# ==============================================================================
# STEP 5: Load and Process GSE234297 (RNA-seq)
# ==============================================================================
message("
==============================================================
  Loading GSE234297 (Validation Cohort 3 -- RNA-seq)
==============================================================")

gse234297    <- getGEO("GSE234297", getGPL = FALSE)[[1]]
pdata_234297 <- pData(gse234297)

# Label extraction
label_col <- if ("disease state:ch1" %in% colnames(pdata_234297)) {
  "disease state:ch1"
} else {
  grep("characteristics", colnames(pdata_234297),
        ignore.case = TRUE, value = TRUE)[1]
}

labels_raw_234297 <- as.character(pdata_234297[[label_col]])

labels_234297 <- ifelse(
  grepl("ALS|sALS|fALS|Amyotrophic",
        labels_raw_234297, ignore.case = TRUE), "ALS",
  ifelse(grepl("Control|Healthy|Normal",
               labels_raw_234297, ignore.case = TRUE), "Control", NA)
)

keep_234297 <- !is.na(labels_234297)
message(sprintf("GSE234297: %d samples (ALS=%d | CON=%d)",
                sum(keep_234297),
                sum(labels_234297[keep_234297] == "ALS"),
                sum(labels_234297[keep_234297] == "Control")))

# Download raw count matrix
if (!dir.exists("GSE234297") ||
    length(list.files("GSE234297",
                       pattern = "raw_counts|count",
                       ignore.case = TRUE)) == 0) {
  message("Downloading GSE234297 supplementary files...")
  suppressMessages(getGEOSuppFiles("GSE234297", makeDirectory = TRUE))
}

counts_file <- grep(
  "raw_counts|count",
  list.files("GSE234297", full.names = TRUE),
  value = TRUE, ignore.case = TRUE
)[1]

if (is.na(counts_file)) {
  stop("GSE234297 count file not found. Check supplementary files.")
}

message(sprintf("Loading: %s", basename(counts_file)))

val_expr_raw_234297 <- read.table(
  gzfile(counts_file),
  header           = TRUE,
  row.names        = 1,
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

# Map identifiers to HGNC symbols
raw_ids <- as.character(rownames(val_expr_raw_234297))

if (grepl("^[0-9]+$", raw_ids[1])) {
  message("GSE234297: Entrez IDs detected...")
  gene_map <- suppressMessages(
    bitr(raw_ids, fromType = "ENTREZID",
         toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  )
  val_expr_raw_234297$ENTREZID <- raw_ids
  df_234297 <- merge(gene_map, val_expr_raw_234297, by = "ENTREZID") %>%
    mutate(across(where(is.numeric), ~ log2(.x + 1))) %>%
    group_by(SYMBOL) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    as.data.frame()
  rownames(df_234297) <- df_234297$SYMBOL
  val_expr_234297     <- as.matrix(df_234297[, -1])

} else if (grepl("^ENSG", raw_ids[1])) {
  message("GSE234297: Ensembl IDs detected...")
  gene_map <- suppressMessages(
    bitr(raw_ids, fromType = "ENSEMBL",
         toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  )
  val_expr_raw_234297$ENSEMBL <- raw_ids
  df_234297 <- merge(gene_map, val_expr_raw_234297, by = "ENSEMBL") %>%
    mutate(across(where(is.numeric), ~ log2(.x + 1))) %>%
    group_by(SYMBOL) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    as.data.frame()
  rownames(df_234297) <- df_234297$SYMBOL
  val_expr_234297     <- as.matrix(df_234297[, -1])

} else {
  message("GSE234297: Gene symbols detected directly.")
  val_expr_234297 <- log2(as.matrix(val_expr_raw_234297) + 1)
}

message(sprintf("GSE234297 matrix: %d genes x %d samples",
                nrow(val_expr_234297), ncol(val_expr_234297)))
message(sprintf("Consensus gene coverage: %d / %d",
                sum(consensus_genes %in% rownames(val_expr_234297)),
                length(consensus_genes)))

val_df_234297 <- prepare_val_data(
  val_expr_234297[, keep_234297],
  labels_234297[keep_234297],
  consensus_genes, logit_model, model_gene_names,
  disc_center, disc_scale_sd
)
roc_234297 <- compute_roc(val_df_234297)

message(sprintf("Arm 5: n=%d | AUC=%.3f (CI: %.3f-%.3f)",
                sum(keep_234297),
                as.numeric(roc_234297$auc),
                as.numeric(roc_234297$ci[1]),
                as.numeric(roc_234297$ci[3])))

save_roc_plot(
  roc_234297,
  title_str    = "GSE234297 -- ALS vs Control (RNA-seq Validation)",
  subtitle_str = sprintf("n=%d | Illumina RNA-seq | Cross-platform",
                          sum(keep_234297)),
  file_prefix  = "Manuscript_Figures/Validation/Figure6A_Arm5_GSE234297",
  line_color   = "#E18727FF"
)

# ==============================================================================
# STEP 6: Combined Multi-Cohort ROC Panel — Figure 6A
# ==============================================================================
message("Generating combined ROC panel (Figure 6A)...")

disc_label <- sprintf("Discovery GSE112676 (n=%d)", ncol(disc_expr))
arm1_label <- sprintf("GSE112680 ALS/CON (n=%d)",   sum(keep_arm1))
arm4_label <- sprintf("GSE28253 (n=%d)",             sum(keep_28253))
arm5_label <- sprintf("GSE234297 RNA-seq (n=%d)",    sum(keep_234297))

roc_list_diag <- list(
  disc_label = list(roc = ml_results$roc_disc, color = "#BC3C29FF"),
  arm1_label = list(roc = roc_arm1,            color = "#0072B5FF"),
  arm4_label = list(roc = roc_28253,            color = "#20854EFF"),
  arm5_label = list(roc = roc_234297,           color = "#E18727FF")
)
names(roc_list_diag) <- c(disc_label, arm1_label, arm4_label, arm5_label)

combined_roc_df <- do.call(rbind, lapply(names(roc_list_diag), function(nm) {
  r <- roc_list_diag[[nm]]$roc
  data.frame(
    Cohort      = sprintf("%s | AUC=%.3f", nm, as.numeric(r$auc)),
    Specificity = r$specificities,
    Sensitivity = r$sensitivities,
    Color       = roc_list_diag[[nm]]$color,
    stringsAsFactors = FALSE
  )
}))

color_map <- unique(combined_roc_df[, c("Cohort", "Color")])
color_vec <- setNames(color_map$Color, color_map$Cohort)

p6A_combined <- ggplot(
  combined_roc_df,
  aes(x = 1 - Specificity, y = Sensitivity, color = Cohort)
) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_color_manual(values = color_vec) +
  theme_bw(base_size = 13) +
  theme(
    legend.position   = c(0.62, 0.22),
    legend.text       = element_text(size = 9),
    legend.title      = element_blank(),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    plot.title        = element_text(face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(hjust = 0.5, color = "grey40"),
    panel.grid.major  = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey95", linewidth = 0.2)
  ) +
  labs(
    title    = "6A. Multi-Cohort Diagnostic Validation",
    subtitle = "Consensus signature across four independent cohorts",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  )


ggsave("Manuscript_Figures/Validation/Figure6A_Combined_ROC.pdf",
       plot = p6A_combined, width = 8, height = 7)
ggsave("Manuscript_Figures/Validation/Figure6A_Combined_ROC.png",
       plot = p6A_combined, width = 8, height = 7, dpi = 600)
message("Figure 6A (combined ROC) saved.")

# ==============================================================================
# STEP 7: Individual Gene ROC Panels — Figure 6B
# ==============================================================================
message("Generating individual gene ROC panels (Figure 6B)...")

build_gene_roc_panel(
  disc_expr, as.character(disc_labels),
  consensus_genes, "Discovery (GSE112676)",
  "Manuscript_Figures/Validation/Figure6B_Discovery"
)

build_gene_roc_panel(
  val_expr_112680[, keep_arm1], labels_arm1,
  consensus_genes, "GSE112680 (ALS vs CON)",
  "Manuscript_Figures/Validation/Figure6B_GSE112680"
)

build_gene_roc_panel(
  val_expr_28253[, keep_28253], labels_28253[keep_28253],
  consensus_genes, "GSE28253 (Microarray)",
  "Manuscript_Figures/Validation/Figure6B_GSE28253"
)

build_gene_roc_panel(
  val_expr_234297[, keep_234297], labels_234297[keep_234297],
  consensus_genes, "GSE234297 (RNA-seq)",
  "Manuscript_Figures/Validation/Figure6B_GSE234297"
)

# ==============================================================================
# STEP 8: Per-Gene AUC Table for BPS
# Arms 1, 4, 5 + discovery only (diagnostic arms)
# Arms 2 and 3 excluded — address clinical specificity, not diagnostic AUC
# ==============================================================================
message("Computing per-gene AUC table for BPS scoring...")

per_gene_aucs <- data.frame(Gene = consensus_genes,
                              stringsAsFactors = FALSE)

per_gene_aucs$AUC_Discovery <- sapply(consensus_genes, function(g)
  get_gene_auc(g, disc_expr, as.character(disc_labels)))

per_gene_aucs$AUC_GSE112680 <- sapply(consensus_genes, function(g)
  get_gene_auc(g, val_expr_112680[, keep_arm1], labels_arm1))

per_gene_aucs$AUC_GSE28253 <- sapply(consensus_genes, function(g)
  get_gene_auc(g, val_expr_28253[, keep_28253],
               labels_28253[keep_28253]))

per_gene_aucs$AUC_GSE234297 <- sapply(consensus_genes, function(g)
  get_gene_auc(g, val_expr_234297[, keep_234297],
               labels_234297[keep_234297]))

per_gene_aucs$Mean_AUC <- rowMeans(
  per_gene_aucs[, c("AUC_Discovery", "AUC_GSE112680",
                     "AUC_GSE28253",  "AUC_GSE234297")],
  na.rm = TRUE
)

message("Per-gene AUC table:")
print_df <- cbind(
  Gene = per_gene_aucs$Gene,
  as.data.frame(lapply(
    per_gene_aucs[, c("AUC_Discovery", "AUC_GSE112680",
                       "AUC_GSE28253",  "AUC_GSE234297", "Mean_AUC")],
    round, 3
  ))
)
print(print_df)

# ==============================================================================
# STEP 9: Biomarker Priority Score (BPS)
# BPS = (0.5 x Scaled Mean AUC) + (0.5 x Scaled kME)
# kME from each gene's own module — not a single reference module
# ==============================================================================
message("Computing Biomarker Priority Scores (BPS)...")

moduleColors_vec <- wgcna_data$moduleColors
kME_all          <- wgcna_data$kME_all

per_gene_aucs$Module <- NA_character_
per_gene_aucs$kME    <- NA_real_

for (idx in seq_len(nrow(per_gene_aucs))) {
  g <- per_gene_aucs$Gene[idx]
  if (g %in% names(moduleColors_vec)) {
    mod <- moduleColors_vec[g]
    per_gene_aucs$Module[idx] <- mod
    kme_col <- paste0("kME", mod)
    if (!is.na(mod) && mod != "grey" &&
        kme_col %in% colnames(kME_all) &&
        g %in% rownames(kME_all)) {
      per_gene_aucs$kME[idx] <- abs(kME_all[g, kme_col])
    }
  }
}

min_max <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (r[2] == r[1]) return(rep(1, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

per_gene_aucs$Scaled_AUC <- min_max(per_gene_aucs$Mean_AUC)
per_gene_aucs$Scaled_kME <- min_max(per_gene_aucs$kME)
per_gene_aucs$BPS <- round(
  (0.5 * per_gene_aucs$Scaled_AUC +
   0.5 * per_gene_aucs$Scaled_kME) * 100, 1
)
per_gene_aucs <- per_gene_aucs %>% arrange(desc(BPS))

message("Biomarker Priority Scores:")
bps_print <- cbind(
  Gene   = per_gene_aucs$Gene,
  Module = per_gene_aucs$Module,
  as.data.frame(lapply(
    per_gene_aucs[, c("AUC_Discovery", "AUC_GSE112680",
                       "Mean_AUC", "kME", "BPS")],
    round, 3
  ))
)
print(bps_print)

write.csv(per_gene_aucs,
          "Processed_Data/Step5_PerGene_AUCs.csv",
          row.names = FALSE)

top_gene <- per_gene_aucs$Gene[1]
top_bps  <- per_gene_aucs$BPS[1]

# ==============================================================================
# STEP 10: Validation Summary Table
# ==============================================================================
summary_table <- data.frame(
  Arm        = c("Discovery","Arm1","Arm2","Arm3","Arm4","Arm5"),
  Dataset    = c("GSE112676","GSE112680","GSE112680","GSE112680",
                  "GSE28253","GSE234297"),
  Comparison = c("ALS/CON","ALS/CON","ALS/MIM",
                  "MIM/CON","ALS/CON","ALS/CON"),
  Platform   = c("Illumina V3","Illumina V4","Illumina V4",
                  "Illumina V4","Agilent","RNA-seq"),
  n = c(
    ncol(disc_expr), sum(keep_arm1), sum(keep_arm2),
    sum(keep_arm3),  sum(keep_28253), sum(keep_234297)
  ),
  AUC = round(c(
    as.numeric(ml_results$roc_disc$auc),
    as.numeric(roc_arm1$auc),
    as.numeric(roc_arm2$auc),
    as.numeric(roc_arm3$auc),
    as.numeric(roc_28253$auc),
    as.numeric(roc_234297$auc)), 3),
  CI_Lower = round(c(
    as.numeric(ml_results$roc_disc$ci[1]),
    as.numeric(roc_arm1$ci[1]),
    as.numeric(roc_arm2$ci[1]),
    as.numeric(roc_arm3$ci[1]),
    as.numeric(roc_28253$ci[1]),
    as.numeric(roc_234297$ci[1])), 3),
  CI_Upper = round(c(
    as.numeric(ml_results$roc_disc$ci[3]),
    as.numeric(roc_arm1$ci[3]),
    as.numeric(roc_arm2$ci[3]),
    as.numeric(roc_arm3$ci[3]),
    as.numeric(roc_28253$ci[3]),
    as.numeric(roc_234297$ci[3])), 3),
  stringsAsFactors = FALSE
)

message("================================================================")
message("VALIDATION SUMMARY TABLE:")
print(summary_table)
message("================================================================")

write.csv(summary_table,
          "Processed_Data/Step5_Validation_Summary.csv",
          row.names = FALSE)

# ==============================================================================
# STEP 11: Gene-Level MIM Diagnostic
# Reports individual gene behavior in Arms 2 and 3
# ==============================================================================
message("Gene-level MIM arm analysis...")

cat("\n=== Individual gene AUC: MIM vs CON (Arm 3) ===\n")
cat("Expected ~0.5 for each gene if signature is ALS-specific\n\n")

for (g in consensus_genes) {
  if (!g %in% rownames(val_expr_112680)) next
  r <- tryCatch(
    roc(as.numeric(labels_arm3 == "ALS"),
        as.numeric(val_expr_112680[g, keep_arm3]),
        quiet = TRUE),
    error = function(e) NULL
  )
  if (!is.null(r)) {
    auc_val <- as.numeric(r$auc)
    flag    <- if (auc_val > 0.65) "  <<< elevated" else ""
    cat(sprintf("%-12s  AUC = %.3f%s\n", g, auc_val, flag))
  }
}

cat("\n=== Individual gene AUC: ALS vs MIM (Arm 2) ===\n")
cat("Higher AUC = better ALS-specific discrimination\n\n")

for (g in consensus_genes) {
  if (!g %in% rownames(val_expr_112680)) next
  r <- tryCatch(
    roc(as.numeric(labels_arm2_binary == "ALS"),
        as.numeric(val_expr_112680[g, keep_arm2]),
        quiet = TRUE),
    error = function(e) NULL
  )
  if (!is.null(r)) {
    auc_val <- as.numeric(r$auc)
    flag    <- if (auc_val > 0.65) "  <<< ALS-specific" else ""
    cat(sprintf("%-12s  AUC = %.3f%s\n", g, auc_val, flag))
  }
}

# ==============================================================================
# STEP 12: Save All Results
# ==============================================================================
saveRDS(
  list(
    val_expr_112680    = val_expr_112680,
    group_112680       = group_112680,
    labels_arm1        = labels_arm1,
    labels_arm2_binary = labels_arm2_binary,
    labels_arm3        = labels_arm3,
    keep_arm1          = keep_arm1,
    keep_arm2          = keep_arm2,
    keep_arm3          = keep_arm3,
    roc_arm1           = roc_arm1,
    roc_arm2           = roc_arm2,
    roc_arm3           = roc_arm3,
    val_expr_28253     = val_expr_28253,
    labels_28253       = labels_28253,
    keep_28253         = keep_28253,
    roc_28253          = roc_28253,
    val_expr_234297    = val_expr_234297,
    labels_234297      = labels_234297,
    keep_234297        = keep_234297,
    roc_234297         = roc_234297,
    per_gene_aucs      = per_gene_aucs,
    summary_table      = summary_table
  ),
  "Processed_Data/Step5_Validation_Results.rds"
)

message("================================================================")
message("Script 04 complete.")
message(sprintf("  Arm 1 GSE112680 ALS/CON : AUC = %.3f",
                as.numeric(roc_arm1$auc)))
message(sprintf("  Arm 2 GSE112680 ALS/MIM : AUC = %.3f",
                as.numeric(roc_arm2$auc)))
message(sprintf("  Arm 3 GSE112680 MIM/CON : AUC = %.3f",
                as.numeric(roc_arm3$auc)))
message(sprintf("  Arm 4 GSE28253  ALS/CON : AUC = %.3f",
                as.numeric(roc_28253$auc)))
message(sprintf("  Arm 5 GSE234297 ALS/CON : AUC = %.3f",
                as.numeric(roc_234297$auc)))
message(sprintf("  Top BPS gene            : %s (BPS = %.1f)",
                top_gene, top_bps))
message("  All outputs: Processed_Data/ and Manuscript_Figures/Validation/")
message("================================================================")
