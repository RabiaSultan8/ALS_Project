# ==============================================================================
# 09_Immune_Subtype_Validation.R
#
# Clinical anchoring and independent cohort validation of the two ALS immune
# subtypes (Immune-A, Immune-B) identified in Script 05.
#
# Added in response to Reviewer 2, Major Comment 4:
#   "The two immune subtypes are based on inferred cell fractions with moderate
#    separation. Their relevance should be tested against clinical variables
#    such as ALSFRS-R decline, survival, site of onset, disease duration, or
#    at least validated in an independent cohort."
#
# PART 1 — Clinical anchoring (GSE112676 ALS samples, n = 223)
#   Survival         : Kaplan-Meier + log-rank test
#   Site of onset    : Fisher's exact test (bulbar vs spinal)
#   Age at onset     : Mann-Whitney U test
#   Sex distribution : Fisher's exact test
#   Note: ALSFRS-R, disease duration, and treatment status are not deposited
#   in the GEO records for this dataset and could not be retrieved.
#
# PART 2 — External cohort validation (GSE112680 ALS samples, n = 164)
#   EPIC deconvolution on GSE112680 ALS samples
#   Centroid projection: each sample assigned to Immune-A or Immune-B by
#     Euclidean distance to discovery-derived centroids
#   Cell fraction comparison between projected subtypes
#
# Technical note:
#   km_final is not saved in Step6_Immune_Results.rds. Discovery centroids are
#   recomputed here as subtype-wise column means of the saved epic_fracs matrix,
#   which is algebraically equivalent to km$centers.
#   GSE112680 expression is back-transformed (2^x - 1) from log2(x+1) scale
#   before EPIC, consistent with the discovery preprocessing in Script 05.
#
# Inputs:
#   Processed_Data/Step5_Validation_Results.rds
#   Processed_Data/Step6_Immune_Results.rds
#   GSE112676_series_matrix.txt.gz  (local cache or GEO download)
#
# Outputs:
#   Processed_Data/Step7_Subtype_Validation.rds
#   Processed_Data/Step7_Clinical_Summary.csv
#   Manuscript_Figures/Immune_Validation/FigS_KM_Survival.{pdf,png}
#   Manuscript_Figures/Immune_Validation/FigS_Clinical_Variables.{pdf,png}
#   Manuscript_Figures/Immune_Validation/FigS_GSE112680_Validation.{pdf,png}
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(GEOquery)
library(survival)
library(survminer)
library(EPIC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

dir.create("Manuscript_Figures/Immune_Validation", recursive = TRUE,
           showWarnings = FALSE)

# ==============================================================================
# STEP 1: Load Immune Analysis and Validation Results
# ==============================================================================
message("Loading immune analysis results and validation data...")

immune_res <- readRDS("Processed_Data/Step6_Immune_Results.rds")
val_res    <- readRDS("Processed_Data/Step5_Validation_Results.rds")

epic_fracs      <- immune_res$epic_fracs
subtype_labels  <- immune_res$subtype_labels
keep_epic       <- immune_res$keep_epic
val_expr_112680 <- val_res$val_expr_112680
group_112680    <- val_res$group_112680

message(sprintf("Immune subtypes loaded: %d ALS samples | %s",
                length(subtype_labels),
                paste(names(table(subtype_labels)),
                      table(subtype_labels), sep = "=", collapse = " | ")))

# ==============================================================================
# STEP 2: Load GSE112676 Phenotype Data
#
# Phenotype data is not saved in Step1_Discovery_Data.rds. The series matrix
# is reloaded here to access survival_yr:ch1, censored:ch1, site_onset:ch1,
# age_onset:ch1, and Sex:ch1.
# ==============================================================================
message("Loading GSE112676 phenotype data...")

if (file.exists("GSE112676_series_matrix.txt.gz")) {
  gse_pheno <- getGEO(filename = "GSE112676_series_matrix.txt.gz",
                      getGPL   = FALSE)
} else {
  message("Local cache not found. Downloading GSE112676 from GEO...")
  gse_pheno <- getGEO("GSE112676", getGPL = FALSE)[[1]]
}

pheno_raw <- pData(gse_pheno)

# Flexible column detection — handles minor naming differences in GEO records
get_col <- function(df, patterns) {
  for (p in patterns) {
    hits <- grep(p, colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  return(NA_character_)
}

col_survival <- get_col(pheno_raw, c("survival_yr", "survival"))
col_censored <- get_col(pheno_raw, c("censored"))
col_site     <- get_col(pheno_raw, c("site_onset", "site.onset"))
col_age      <- get_col(pheno_raw, c("age_onset", "age.onset"))
col_sex      <- get_col(pheno_raw, c("^Sex", "gender", "sex"))
col_diag     <- get_col(pheno_raw, c("diagnosis"))

message(sprintf(
  "Matched columns: survival=%s | censored=%s | site=%s | age=%s | sex=%s",
  col_survival, col_censored, col_site, col_age, col_sex
))

# ==============================================================================
# STEP 3: Build Clinical Data Frame Aligned to Immune Subtypes
# ==============================================================================
message("Building clinical data frame...")

pheno_als <- pheno_raw %>%
  filter(toupper(trimws(.data[[col_diag]])) == "ALS") %>%
  select(
    SampleID    = geo_accession,
    survival_yr = all_of(col_survival),
    censored    = all_of(col_censored),
    site_onset  = all_of(col_site),
    age_onset   = all_of(col_age),
    sex         = all_of(col_sex)
  ) %>%
  mutate(
    survival_yr = suppressWarnings(as.numeric(as.character(survival_yr))),
    censored    = suppressWarnings(as.numeric(as.character(censored))),
    age_onset   = suppressWarnings(as.numeric(as.character(age_onset))),
    sex         = toupper(trimws(as.character(sex))),
    site_onset  = toupper(trimws(as.character(site_onset)))
  )

subtype_df <- data.frame(
  SampleID = names(subtype_labels),
  Subtype  = subtype_labels,
  stringsAsFactors = FALSE
)

clinical_df <- pheno_als %>%
  inner_join(subtype_df, by = "SampleID")

message(sprintf("Clinical data frame: %d ALS samples with subtype and clinical data",
                nrow(clinical_df)))

na_counts <- sapply(
  c("survival_yr", "censored", "site_onset", "age_onset", "sex"),
  function(x) sum(is.na(clinical_df[[x]]))
)
message("NA counts per clinical field:")
print(na_counts)

# ==============================================================================
# STEP 4: Survival Analysis (Kaplan-Meier + Log-Rank)
#
# Censoring convention: censored = 0 means patient died (event occurred),
# censored = 1 means still alive at last follow-up (right-censored).
# Verified from GSE112676 GEO series record.
# ==============================================================================
message("Running survival analysis...")

surv_df <- clinical_df %>%
  filter(!is.na(survival_yr) & !is.na(censored)) %>%
  mutate(event = ifelse(censored == 0, 1, 0))

message(sprintf("Survival analysis: n = %d | events = %d | censored = %d",
                nrow(surv_df), sum(surv_df$event), sum(surv_df$event == 0)))

surv_obj <- Surv(time = surv_df$survival_yr, event = surv_df$event)
km_fit   <- survfit(surv_obj ~ Subtype, data = surv_df)
log_rank <- survdiff(surv_obj ~ Subtype, data = surv_df)
lr_pval  <- 1 - pchisq(log_rank$chisq, df = length(log_rank$n) - 1)

message(sprintf("Log-rank test: chi-sq = %.3f | p = %.4f", log_rank$chisq, lr_pval))
message("Median survival (years) by subtype:")
print(quantile(km_fit, probs = 0.5))

# -- Figure S[X]: Kaplan-Meier Survival Curve ----------------------------------
km_plot <- ggsurvplot(
  km_fit,
  data               = surv_df,
  pval               = TRUE,
  pval.method        = FALSE,
  conf.int           = TRUE,
  conf.int.alpha     = 0.15,
  risk.table         = TRUE,
  risk.table.height  = 0.28,
  xlab               = "Survival (years)",
  ylab               = "Overall Survival Probability",
  title              = "Kaplan-Meier Survival - ALS Immune Subtypes (GSE112676)",
  legend.title       = "Immune Subtype",
  legend.labs        = names(table(surv_df$Subtype)),
  palette            = c("#BC3C29FF", "#0072B5FF"),
  ggtheme            = theme_bw(base_size = 13),
  fontsize           = 4.5
)

pdf("Manuscript_Figures/Immune_Validation/FigS_KM_Survival.pdf",
    width = 8, height = 8)
print(km_plot)
dev.off()

png("Manuscript_Figures/Immune_Validation/FigS_KM_Survival.png",
    width = 8, height = 8, units = "in", res = 600)
print(km_plot)
dev.off()

message("Kaplan-Meier figure saved.")

# ==============================================================================
# STEP 5: Site of Onset — Fisher's Exact Test
#
# Values in site_onset:ch1 are single-letter codes: B = bulbar, S = spinal.
# ==============================================================================
message("Running site of onset analysis...")

clinical_df <- clinical_df %>%
  mutate(
    onset_bs = ifelse(site_onset == "B", "Bulbar",
               ifelse(site_onset == "S", "Spinal", NA))
  )

message("Site of onset distribution by subtype:")
onset_table <- table(clinical_df$Subtype, clinical_df$onset_bs)
print(onset_table)

fisher_onset <- fisher.test(onset_table)
message(sprintf("Fisher's exact (site of onset): OR = %.3f | p = %.4f",
                fisher_onset$estimate, fisher_onset$p.value))

# ==============================================================================
# STEP 6: Age at Onset — Mann-Whitney U
# ==============================================================================
message("Running age at onset analysis...")

age_df   <- clinical_df %>% filter(!is.na(age_onset))
age_test <- wilcox.test(age_onset ~ Subtype, data = age_df,
                         exact = FALSE, conf.int = TRUE)

message(sprintf("Mann-Whitney U (age at onset): W = %.1f | p = %.4f",
                age_test$statistic, age_test$p.value))

age_df %>%
  group_by(Subtype) %>%
  summarise(
    n      = n(),
    Median = round(median(age_onset, na.rm = TRUE), 1),
    IQR_lo = round(quantile(age_onset, 0.25, na.rm = TRUE), 1),
    IQR_hi = round(quantile(age_onset, 0.75, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  print()

# ==============================================================================
# STEP 7: Sex Distribution — Fisher's Exact Test
# ==============================================================================
message("Running sex distribution analysis...")

sex_df <- clinical_df %>%
  filter(!is.na(sex) & sex != "") %>%
  mutate(
    sex_bin = case_when(
      grepl("^M|^1", sex, ignore.case = TRUE) ~ "Male",
      grepl("^F|^2", sex, ignore.case = TRUE) ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sex_bin))

sex_table  <- table(sex_df$Subtype, sex_df$sex_bin)
message("Sex distribution by subtype:")
print(sex_table)

fisher_sex <- fisher.test(sex_table)
message(sprintf("Fisher's exact (sex): OR = %.3f | p = %.4f",
                fisher_sex$estimate, fisher_sex$p.value))

# ==============================================================================
# STEP 8: Clinical Variables Panel Figure
# ==============================================================================
message("Generating clinical variables panel figure...")

p_age <- ggplot(age_df, aes(x = Subtype, y = age_onset, fill = Subtype)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 16, outlier.size = 1.5,
               width = 0.5) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35, color = "grey30") +
  scale_fill_manual(values = c("Immune-A" = "#BC3C29FF",
                                "Immune-B" = "#0072B5FF")) +
  annotate("text",
           x     = 1.5,
           y     = max(age_df$age_onset, na.rm = TRUE) + 2,
           label = sprintf("p = %.3f", age_test$p.value),
           size  = 4.5, fontface = "bold") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title      = element_text(face = "bold", hjust = 0.5),
        axis.text       = element_text(color = "black")) +
  labs(title = "Age at Onset", x = "", y = "Age (years)")

sex_prop <- sex_df %>%
  count(Subtype, sex_bin) %>%
  group_by(Subtype) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_sex <- ggplot(sex_prop, aes(x = Subtype, y = pct, fill = sex_bin)) +
  geom_col(position = "stack", color = "black", linewidth = 0.3, width = 0.5) +
  scale_fill_manual(values = c("Male" = "#7876B1FF", "Female" = "#6F99ADFF"),
                    name = "Sex") +
  annotate("text",
           x     = 1.5,
           y     = 108,
           label = sprintf("p = %.3f", fisher_sex$p.value),
           size  = 4.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 118)) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "Sex Distribution", x = "", y = "Percentage (%)")

clinical_panel <- ggarrange(p_age, p_sex,
                              ncol   = 2,
                              labels = c("A", "B"),
                              align  = "hv")

ggsave("Manuscript_Figures/Immune_Validation/FigS_Clinical_Variables.pdf",
       plot = clinical_panel, width = 9, height = 5)
ggsave("Manuscript_Figures/Immune_Validation/FigS_Clinical_Variables.png",
       plot = clinical_panel, width = 9, height = 5, dpi = 600)
message("Clinical variables figure saved.")

# ==============================================================================
# STEP 9: Clinical Summary Table
# ==============================================================================
message("Generating clinical summary table...")

clinical_summary <- clinical_df %>%
  group_by(Subtype) %>%
  summarise(
    n                  = n(),
    Age_median         = round(median(age_onset, na.rm = TRUE), 1),
    Age_IQR_lo         = round(quantile(age_onset, 0.25, na.rm = TRUE), 1),
    Age_IQR_hi         = round(quantile(age_onset, 0.75, na.rm = TRUE), 1),
    Survival_median_yr = round(median(survival_yr, na.rm = TRUE), 2),
    Bulbar_pct         = round(mean(onset_bs == "Bulbar", na.rm = TRUE) * 100, 1),
    Male_pct           = round(mean(grepl("^M|^1", sex, ignore.case = TRUE),
                                     na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

message("Clinical summary by subtype:")
print(clinical_summary)

write.csv(clinical_summary,
          "Processed_Data/Step7_Clinical_Summary.csv",
          row.names = FALSE)

# ==============================================================================
# STEP 10: EPIC Deconvolution — GSE112680 ALS Samples
#
# Expression is back-transformed from log2(x+1) to linear scale before EPIC,
# consistent with the discovery preprocessing applied in Script 05.
# ==============================================================================
message("Running EPIC deconvolution on GSE112680 ALS samples...")

keep_als_112680 <- toupper(trimws(group_112680)) == "ALS"
expr_als_112680 <- val_expr_112680[, keep_als_112680, drop = FALSE]

message(sprintf("GSE112680 ALS samples: n = %d", ncol(expr_als_112680)))

expr_als_linear <- 2^expr_als_112680 - 1

epic_112680 <- tryCatch(
  EPIC::EPIC(bulk = expr_als_linear),
  error = function(e) {
    message(sprintf("EPIC error: %s", e$message))
    NULL
  }
)

if (is.null(epic_112680)) {
  stop("EPIC failed on GSE112680. Check gene coverage and expression scale.")
}

cell_fracs_112680           <- as.data.frame(epic_112680$mRNAProportions)
cell_fracs_112680$SampleID  <- rownames(cell_fracs_112680)

message(sprintf("GSE112680 EPIC complete: %d cell types | %d samples",
                ncol(cell_fracs_112680) - 1, nrow(cell_fracs_112680)))

# ==============================================================================
# STEP 11: Discovery Centroid Derivation
#
# km_final is not saved in Step6_Immune_Results.rds. Centroids are recomputed
# as subtype-wise column means of the discovery EPIC fractions, which is
# algebraically equivalent to km$centers.
# ==============================================================================
message("Computing discovery subtype centroids...")

disc_frac_mat <- as.matrix(
  epic_fracs[match(names(subtype_labels), epic_fracs$Sample),
             keep_epic, drop = FALSE]
)
rownames(disc_frac_mat) <- names(subtype_labels)
disc_frac_mat           <- disc_frac_mat[complete.cases(disc_frac_mat), ,
                                          drop = FALSE]

valid_labels <- subtype_labels[rownames(disc_frac_mat)]

centroid_A <- colMeans(
  disc_frac_mat[valid_labels == "Immune-A", , drop = FALSE],
  na.rm = TRUE
)
centroid_B <- colMeans(
  disc_frac_mat[valid_labels == "Immune-B", , drop = FALSE],
  na.rm = TRUE
)

message("Discovery centroids:")
print(round(rbind(`Immune-A` = centroid_A, `Immune-B` = centroid_B), 4))

# ==============================================================================
# STEP 12: Centroid Projection — GSE112680 ALS Samples
#
# Each sample is assigned to Immune-A or Immune-B by Euclidean distance to
# the discovery centroids computed in Step 11.
# ==============================================================================
message("Projecting GSE112680 ALS samples onto discovery subtypes...")

common_cells <- intersect(keep_epic, colnames(cell_fracs_112680))

message(sprintf("Cell types used for projection: %s",
                paste(common_cells, collapse = ", ")))

val_frac_mat           <- as.matrix(
  cell_fracs_112680[, common_cells, drop = FALSE]
)
rownames(val_frac_mat) <- cell_fracs_112680$SampleID

subtypes_112680 <- setNames(
  apply(val_frac_mat, 1, function(r) {
    dA <- sqrt(sum((r[common_cells] - centroid_A[common_cells])^2, na.rm = TRUE))
    dB <- sqrt(sum((r[common_cells] - centroid_B[common_cells])^2, na.rm = TRUE))
    if (dA <= dB) "Immune-A" else "Immune-B"
  }),
  rownames(val_frac_mat)
)

message("GSE112680 subtype distribution:")
print(table(subtypes_112680))

prop_disc <- prop.table(table(subtype_labels))
prop_val  <- prop.table(table(subtypes_112680))

prop_compare <- data.frame(
  Subtype       = names(prop_disc),
  Discovery_pct = round(as.numeric(prop_disc) * 100, 1),
  GSE112680_pct = round(as.numeric(prop_val[names(prop_disc)]) * 100, 1)
)
message("Subtype proportions — Discovery vs GSE112680:")
print(prop_compare)

cohort_table  <- rbind(
  Discovery = table(subtype_labels),
  GSE112680 = table(subtypes_112680)[names(table(subtype_labels))]
)
fisher_cohort <- fisher.test(cohort_table)
message(sprintf("Proportion replication Fisher's: OR = %.3f | p = %.4f",
                fisher_cohort$estimate, fisher_cohort$p.value))

# ==============================================================================
# STEP 13: GSE112680 Validation Figure
# ==============================================================================
message("Generating GSE112680 validation figure...")

val_long <- cell_fracs_112680 %>%
  select(SampleID, all_of(common_cells)) %>%
  mutate(Subtype = subtypes_112680[SampleID]) %>%
  filter(!is.na(Subtype)) %>%
  pivot_longer(cols      = all_of(common_cells),
               names_to  = "CellType",
               values_to = "Fraction")

p_val <- ggplot(val_long, aes(x = Subtype, y = Fraction, fill = Subtype)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 16, outlier.size = 0.8) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Immune-A" = "#BC3C29FF",
                                "Immune-B" = "#0072B5FF")) +
  stat_compare_means(method      = "wilcox.test",
                     label       = "p.format",
                     size        = 3.2,
                     fontface    = "bold",
                     label.y.npc = 0.92) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.text      = element_text(face = "bold", size = 9),
        plot.title      = element_text(face = "bold", hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40")) +
  labs(
    title    = "EPIC Cell Fractions by Immune Subtype -- GSE112680 ALS (n = 164)",
    subtitle = sprintf(
      "Centroid-projected from GSE112676 discovery | Proportion p = %.3f",
      fisher_cohort$p.value
    ),
    x = "",
    y = "mRNA Fraction"
  )

ggsave("Manuscript_Figures/Immune_Validation/FigS_GSE112680_Validation.pdf",
       plot = p_val, width = 13, height = 8)
ggsave("Manuscript_Figures/Immune_Validation/FigS_GSE112680_Validation.png",
       plot = p_val, width = 13, height = 8, dpi = 600)
message("GSE112680 validation figure saved.")

# ==============================================================================
# STEP 14: Save All Results
# ==============================================================================
saveRDS(
  list(
    # Clinical anchoring
    clinical_df      = clinical_df,
    clinical_summary = clinical_summary,
    km_fit           = km_fit,
    log_rank         = log_rank,
    lr_pval          = lr_pval,
    age_test         = age_test,
    fisher_onset     = fisher_onset,
    fisher_sex       = fisher_sex,

    # External validation
    epic_112680       = epic_112680,
    cell_fracs_112680 = cell_fracs_112680,
    subtypes_112680   = subtypes_112680,
    centroid_A        = centroid_A,
    centroid_B        = centroid_B,
    common_cells      = common_cells,
    prop_compare      = prop_compare,
    fisher_cohort     = fisher_cohort
  ),
  "Processed_Data/Step7_Subtype_Validation.rds"
)

message("================================================================")
message("Script 09 complete.")
message(sprintf("  Clinical anchoring (GSE112676, n = %d ALS samples):",
                nrow(clinical_df)))
message(sprintf("    Survival log-rank p     : %.4f", lr_pval))
message(sprintf("    Immune-A median survival: %.2f years",
                median(surv_df$survival_yr[surv_df$Subtype == "Immune-A"],
                       na.rm = TRUE)))
message(sprintf("    Immune-B median survival: %.2f years",
                median(surv_df$survival_yr[surv_df$Subtype == "Immune-B"],
                       na.rm = TRUE)))
message(sprintf("    Age at onset p          : %.4f", age_test$p.value))
message(sprintf("    Site of onset p         : %.4f", fisher_onset$p.value))
message(sprintf("    Sex distribution p      : %.4f", fisher_sex$p.value))
message(sprintf("  External validation (GSE112680, n = %d ALS samples):",
                ncol(expr_als_112680)))
message(sprintf("    Subtype distribution    : %s",
                paste(names(table(subtypes_112680)),
                      table(subtypes_112680), sep = "=", collapse = " | ")))
message(sprintf("    Proportion replication p: %.4f", fisher_cohort$p.value))
message(sprintf("    Cell type projection    : %d common types used",
                length(common_cells)))
message("  Figures : Manuscript_Figures/Immune_Validation/")
message("  Data    : Processed_Data/Step7_Subtype_Validation.rds")
message("================================================================")
