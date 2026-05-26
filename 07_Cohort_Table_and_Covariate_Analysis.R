# ==============================================================================
# 07_Cohort_Table_and_Covariate_Analysis.R
#
# Cohort characteristics table and sex-adjusted differential expression
# sensitivity analysis for the GSE112676 discovery and GSE112680 validation
# cohorts.
#
# Added in response to Reviewer 2, Major Comment 2:
#   "Age, sex, disease duration, ALSFRS-R, site of onset, and treatment status
#    may influence whole-blood expression. The authors should provide a cohort
#    characteristics table and repeat differential expression, WGCNA, and
#    classification analyses with appropriate covariate adjustment where
#    metadata permit."
#
# PART 1 — Cohort characteristics table
#   Summarises sex, age at onset, site of onset, and survival for all groups
#   across GSE112676 (ALS, CON) and GSE112680 (ALS, CON, MIM).
#   Saved as Supplementary Table S3.
#
# PART 2 — Sex-adjusted differential expression sensitivity analysis
#   Sex is the only covariate available for all samples in GSE112676.
#   Age at onset cannot be used: it is absent for control subjects.
#   ALSFRS-R, disease duration, and treatment status are not deposited in
#   the GEO records for this dataset.
#   A sex-adjusted limma model is run and the six consensus genes are
#   compared between unadjusted and sex-adjusted results.
#
# Inputs:
#   Processed_Data/Step1_Discovery_Data.rds
#   Processed_Data/Step2_DEG_all.csv
#   Processed_Data/Step4_ML_Results.rds
#   GSE112676_series_matrix.txt.gz  (local cache or GEO download)
#   GSE112680_series_matrix.txt.gz  (local cache or GEO download)
#
# Outputs:
#   Processed_Data/Cohort_Characteristics_Table.csv   (Supplementary Table S3)
#   Processed_Data/Step2_DEG_adj_all.csv
#   Processed_Data/Step2_DEG_adj_sig.csv
#   Processed_Data/SixGene_SexAdjustment_Check.csv
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(GEOquery)
library(limma)
library(dplyr)

dir.create("Processed_Data", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: Load Saved Pipeline Results
# ==============================================================================
message("Loading saved pipeline results...")

step1_data     <- readRDS("Processed_Data/Step1_Discovery_Data.rds")
ml_results     <- readRDS("Processed_Data/Step4_ML_Results.rds")

expr_discovery <- step1_data$expr_discovery   # genes x samples (post-QC)
group_factor   <- step1_data$group_factor     # factor: ALS / Control

deg_all        <- read.csv("Processed_Data/Step2_DEG_all.csv",
                            row.names = 1, stringsAsFactors = FALSE)

consensus_genes <- ml_results$consensus_genes

message(sprintf("Discovery matrix: %d genes x %d samples",
                nrow(expr_discovery), ncol(expr_discovery)))
message(sprintf("Consensus genes : %s", paste(consensus_genes, collapse = ", ")))

# ==============================================================================
# STEP 2: Load GSE112676 Phenotype Data
#
# Sex, age at onset, site of onset, and survival are extracted from the GEO
# series matrix. pheno is matched to post-QC samples by geo_accession.
# Dutch sex encoding: "V" (vrouw) = Female, "M" = Male.
# ==============================================================================
message("Loading GSE112676 phenotype data...")

if (file.exists("GSE112676_series_matrix.txt.gz")) {
  gse676 <- getGEO(filename = "GSE112676_series_matrix.txt.gz",
                   getGPL   = FALSE)
} else {
  message("Local cache not found. Downloading GSE112676...")
  gse676 <- getGEO("GSE112676", getGPL = FALSE)[[1]]
}

pheno676_raw <- pData(gse676)

# ==============================================================================
# STEP 3: Load GSE112680 Phenotype Data
# ==============================================================================
message("Loading GSE112680 phenotype data...")

if (file.exists("GSE112680_series_matrix.txt.gz")) {
  gse680 <- getGEO(filename = "GSE112680_series_matrix.txt.gz",
                   getGPL   = FALSE)
} else {
  message("Local cache not found. Downloading GSE112680...")
  gse680 <- getGEO("GSE112680", getGPL = FALSE)[[1]]
}

pheno680_raw <- pData(gse680)

# ==============================================================================
# STEP 4: Clean Phenotype Data
#
# Dutch "V" (vrouw) recoded to "F" for consistency.
# Survival and age at onset coerced to numeric; non-parseable entries become NA.
# ==============================================================================
message("Cleaning phenotype data...")

clean_pheno <- function(pd, cohort_name) {
  df <- data.frame(
    Cohort      = cohort_name,
    SampleID    = pd$geo_accession,
    Diagnosis   = toupper(trimws(as.character(pd$`diagnosis:ch1`))),
    Sex         = toupper(trimws(as.character(pd$`Sex:ch1`))),
    Age_onset   = suppressWarnings(
                    as.numeric(as.character(pd$`age_onset:ch1`))),
    Site_onset  = toupper(trimws(as.character(pd$`site_onset:ch1`))),
    Survival_yr = suppressWarnings(
                    as.numeric(as.character(pd$`survival_yr:ch1`))),
    Censored    = suppressWarnings(
                    as.numeric(as.character(pd$`censored:ch1`))),
    stringsAsFactors = FALSE
  )
  df$Sex[df$Sex == "V"] <- "F"   # Dutch: vrouw = female
  df
}

pheno676 <- clean_pheno(pheno676_raw, "GSE112676") %>%
  filter(Diagnosis %in% c("ALS", "CON"))

pheno680 <- clean_pheno(pheno680_raw, "GSE112680")

message(sprintf("GSE112676 cleaned: %d samples | Groups: %s",
                nrow(pheno676),
                paste(names(table(pheno676$Diagnosis)),
                      table(pheno676$Diagnosis), sep = "=", collapse = " | ")))
message(sprintf("GSE112680 cleaned: %d samples | Groups: %s",
                nrow(pheno680),
                paste(names(table(pheno680$Diagnosis)),
                      table(pheno680$Diagnosis), sep = "=", collapse = " | ")))

# ==============================================================================
# STEP 5: Build Cohort Characteristics Table (Supplementary Table S3)
# ==============================================================================
message("Building cohort characteristics table...")

summarise_group <- function(df, cohort, diag) {

  sub <- df %>% filter(Cohort == cohort, Diagnosis == diag)
  n   <- nrow(sub)

  if (n == 0) {
    return(data.frame(
      Cohort          = cohort,
      Group           = diag,
      N               = 0,
      Sex_ratio       = NA,
      Age_onset_mean  = NA,
      Age_onset_sd    = NA,
      Age_onset_range = NA,
      Site_spinal_n   = NA,
      Site_bulbar_n   = NA,
      Site_other_n    = NA,
      Survival_mean   = NA,
      Survival_sd     = NA,
      stringsAsFactors = FALSE
    ))
  }

  male_n    <- sum(sub$Sex == "M", na.rm = TRUE)
  female_n  <- sum(sub$Sex == "F", na.rm = TRUE)
  age_vals  <- sub$Age_onset[!is.na(sub$Age_onset)]
  surv_vals <- sub$Survival_yr[!is.na(sub$Survival_yr)]
  site_s    <- sum(sub$Site_onset == "S", na.rm = TRUE)
  site_b    <- sum(sub$Site_onset == "B", na.rm = TRUE)
  site_o    <- sum(
    !is.na(sub$Site_onset) &
    sub$Site_onset != "NA" &
    sub$Site_onset != "S"  &
    sub$Site_onset != "B",
    na.rm = TRUE
  )

  data.frame(
    Cohort          = cohort,
    Group           = diag,
    N               = n,
    Sex_ratio       = sprintf("%d M / %d F", male_n, female_n),
    Age_onset_mean  = ifelse(length(age_vals)  > 0,
                             round(mean(age_vals),  1), NA),
    Age_onset_sd    = ifelse(length(age_vals)  > 0,
                             round(sd(age_vals),    1), NA),
    Age_onset_range = ifelse(length(age_vals)  > 0,
                             sprintf("%.1f-%.1f",
                                     min(age_vals), max(age_vals)), NA),
    Site_spinal_n   = site_s,
    Site_bulbar_n   = site_b,
    Site_other_n    = site_o,
    Survival_mean   = ifelse(length(surv_vals) > 0,
                             round(mean(surv_vals), 2), NA),
    Survival_sd     = ifelse(length(surv_vals) > 0,
                             round(sd(surv_vals),   2), NA),
    stringsAsFactors = FALSE
  )
}

cohort_table <- bind_rows(
  lapply(c("ALS", "CON"),
         function(g) summarise_group(pheno676, "GSE112676", g)),
  lapply(c("ALS", "CON", "MIM"),
         function(g) summarise_group(pheno680, "GSE112680", g))
)

message("Cohort characteristics table:")
print(cohort_table)

write.csv(cohort_table,
          "Processed_Data/Cohort_Characteristics_Table.csv",
          row.names = FALSE)

message("Cohort characteristics table saved.")

# ==============================================================================
# STEP 6: Extract Sex Variable and Align to Post-QC Samples
#
# pheno676 contains all samples including any outliers removed during QC in
# Script 01. expr_discovery columns represent the post-QC sample IDs. Sex is
# extracted for post-QC samples only by matching on geo_accession.
# ==============================================================================
message("Extracting and aligning sex variable to post-QC samples...")

sex_lookup <- setNames(pheno676$Sex, pheno676$SampleID)
sex_aligned <- sex_lookup[colnames(expr_discovery)]

sex_factor <- factor(
  ifelse(tolower(sex_aligned) == "m", "Male", "Female"),
  levels = c("Female", "Male")
)
names(sex_factor) <- colnames(expr_discovery)

message(sprintf(
  "Sex distribution (post-QC): Male = %d | Female = %d | NA = %d",
  sum(sex_factor == "Male",   na.rm = TRUE),
  sum(sex_factor == "Female", na.rm = TRUE),
  sum(is.na(sex_factor))
))

# ==============================================================================
# STEP 7: Sex-Adjusted Differential Expression Analysis
#
# Sex is the only covariate available for all samples. Age at onset is absent
# for control subjects and cannot be used uniformly. ALSFRS-R, disease
# duration, and treatment status are not deposited in GEO for this dataset.
#
# Design: ~ 0 + Group + Sex
# Contrast: ALS - Control (identical to Script 01 unadjusted analysis)
# ==============================================================================
message("Running sex-adjusted differential expression analysis...")

na_sex    <- is.na(sex_factor)
expr_adj  <- expr_discovery[, !na_sex]
group_adj <- group_factor[!na_sex]
sex_adj   <- droplevels(sex_factor[!na_sex])

if (any(na_sex)) {
  message(sprintf("  %d sample(s) with missing sex excluded from adjusted analysis.",
                  sum(na_sex)))
}

message(sprintf("Adjusted analysis: %d samples | ALS: %d | CON: %d",
                ncol(expr_adj),
                sum(group_adj == "ALS"),
                sum(group_adj == "Control")))

design_adj  <- model.matrix(~ 0 + group_adj + sex_adj)
colnames(design_adj) <- c("Control", "ALS", "SexMale")

fit_adj    <- lmFit(expr_adj, design_adj)
contr_adj  <- makeContrasts(ALS_vs_Control = ALS - Control,
                              levels = design_adj)
fit2_adj   <- eBayes(contrasts.fit(fit_adj, contr_adj))

deg_adj_all      <- topTable(fit2_adj, coef = 1, number = Inf,
                               adjust.method = "BH")
deg_adj_all$gene <- rownames(deg_adj_all)

deg_adj_sig <- deg_adj_all %>%
  filter(abs(logFC) >= 0.25 & adj.P.Val < 0.05)

message(sprintf("Sex-adjusted DEGs: %d significant | Up: %d | Down: %d",
                nrow(deg_adj_sig),
                sum(deg_adj_sig$logFC > 0),
                sum(deg_adj_sig$logFC < 0)))

write.csv(deg_adj_all,
          "Processed_Data/Step2_DEG_adj_all.csv",
          row.names = FALSE)
write.csv(deg_adj_sig,
          "Processed_Data/Step2_DEG_adj_sig.csv",
          row.names = FALSE)

# ==============================================================================
# STEP 8: Consensus Gene Survival Check
#
# Each of the six consensus genes is compared between the unadjusted and
# sex-adjusted models. A gene is considered to survive adjustment if it
# retains FDR < 0.05 and |logFC| >= 0.25 in the adjusted model.
# ==============================================================================
message("Checking consensus gene survival under sex adjustment...")

# Ensure deg_all is indexed by gene name for direct lookup
if (!"gene" %in% colnames(deg_all)) {
  deg_all$gene <- rownames(deg_all)
}
rownames(deg_all)     <- deg_all$gene
rownames(deg_adj_all) <- deg_adj_all$gene

available_genes <- intersect(consensus_genes, rownames(deg_all))
missing_genes   <- setdiff(consensus_genes, rownames(deg_all))

if (length(missing_genes) > 0) {
  message(sprintf("WARNING: %d consensus gene(s) not in DEG table: %s",
                  length(missing_genes),
                  paste(missing_genes, collapse = ", ")))
}

check_df <- data.frame(
  Gene          = available_genes,
  logFC_unadj   = round(deg_all[available_genes, "logFC"],       4),
  FDR_unadj     = signif(deg_all[available_genes, "adj.P.Val"],  3),
  logFC_adj     = round(deg_adj_all[available_genes, "logFC"],     4),
  FDR_adj       = signif(deg_adj_all[available_genes, "adj.P.Val"], 3),
  row.names     = NULL,
  stringsAsFactors = FALSE
)
check_df$Survives_Adjustment <- (
  check_df$FDR_adj   < 0.05 &
  abs(check_df$logFC_adj) >= 0.25
)

message("\n=== Six-gene signature: unadjusted vs sex-adjusted ===")
print(check_df)

n_survive <- sum(check_df$Survives_Adjustment)
message(sprintf("\n%d / %d consensus genes survive sex adjustment.",
                n_survive, nrow(check_df)))

write.csv(check_df,
          "Processed_Data/SixGene_SexAdjustment_Check.csv",
          row.names = FALSE)

# ==============================================================================
# STEP 9: Save All Results
# ==============================================================================
saveRDS(
  list(
    cohort_table    = cohort_table,
    sex_factor      = sex_factor,
    deg_adj_all     = deg_adj_all,
    deg_adj_sig     = deg_adj_sig,
    check_df        = check_df
  ),
  "Processed_Data/Step9_Covariate_Analysis.rds"
)

message("================================================================")
message("Script 07 complete.")
message("  PART 1 — Cohort characteristics table:")
message(sprintf("    Cohorts          : GSE112676 (ALS, CON) | GSE112680 (ALS, CON, MIM)"))
message(sprintf("    Saved to         : Processed_Data/Cohort_Characteristics_Table.csv"))
message("  PART 2 — Sex-adjusted differential expression:")
message(sprintf("    Covariate        : Sex (only covariate available for all samples)"))
message(sprintf("    Adjusted DEGs    : %d significant (|logFC| >= 0.25, FDR < 0.05)",
                nrow(deg_adj_sig)))
message(sprintf("    Consensus genes  : %d / %d survive sex adjustment",
                n_survive, nrow(check_df)))
message(sprintf("    Saved to         : Processed_Data/Step2_DEG_adj_all.csv"))
message(sprintf("                       Processed_Data/SixGene_SexAdjustment_Check.csv"))
message("================================================================")
