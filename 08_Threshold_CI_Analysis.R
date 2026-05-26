# ==============================================================================
# 08_Threshold_CI_Analysis.R
#
# Youden-optimal decision threshold with 2000-bootstrap 95% confidence
# intervals for sensitivity and specificity. Computed for GSE112680 Arms 1,
# 2, and 3, where threshold-dependent operating-point estimation is clinically
# informative.
#
# The Youden index (J = sensitivity + specificity - 1) identifies the decision
# threshold that maximises the sum of sensitivity and specificity simultaneously.
# When multiple thresholds produce identical J, the first (lowest) is retained
# via best.policy = "first". Bootstrap resampling is stratified by class label
# to preserve class proportions across replicates.
#
# Added in response to Reviewer 2, Minor Comment 5:
#   "Report confidence intervals for sensitivity, specificity, and clinically
#    relevant thresholds, not only AUC."
#
# Inputs:
#   Processed_Data/Step4_ML_Results.rds
#   Processed_Data/Step5_Validation_Results.rds
#
# Outputs:
#   Processed_Data/Supplementary_Table_S4_Threshold_CIs.csv
# ==============================================================================

options(timeout = 600)
set.seed(1122)

library(pROC)

dir.create("Processed_Data", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: Load Saved ROC Objects
# ==============================================================================
message("Loading ROC objects...")

ml_res  <- readRDS("Processed_Data/Step4_ML_Results.rds")
val_res <- readRDS("Processed_Data/Step5_Validation_Results.rds")

arms <- list(
  list(name       = "Arm 1",
       comparison = "ALS vs CON",
       n          = 301,
       roc        = val_res$roc_arm1),
  list(name       = "Arm 2",
       comparison = "ALS vs MIM",
       n          = 239,
       roc        = val_res$roc_arm2),
  list(name       = "Arm 3",
       comparison = "MIM vs CON",
       n          = 212,
       roc        = val_res$roc_arm3)
)

message(sprintf("ROC objects loaded: %d arms | Cohort: GSE112680",
                length(arms)))

# ==============================================================================
# STEP 2: Operating-Point Extraction Function
#
# Point estimates via coords() at Youden-optimal threshold.
# 95% CIs via ci.coords() with 2000 stratified bootstrap replicates.
# Ties broken by best.policy = "first" (lowest threshold retained).
# ==============================================================================
extract_op <- function(arm) {

  roc_obj    <- arm$roc
  arm_name   <- arm$name
  comparison <- arm$comparison
  n          <- arm$n

  auc_val <- round(as.numeric(roc_obj$auc), 4)

  # -- Point estimates --
  pt <- coords(
    roc_obj,
    x           = "best",
    best.method = "youden",
    ret         = c("threshold", "sensitivity", "specificity"),
    transpose   = FALSE
  )
  pt   <- as.data.frame(pt)[1, , drop = FALSE]   # first row if ties
  thr  <- as.numeric(pt$threshold)
  sens <- as.numeric(pt$sensitivity)
  spec <- as.numeric(pt$specificity)

  message(sprintf("  [%s] AUC = %.4f | Threshold = %.4f | Sens = %.4f | Spec = %.4f",
                  arm_name, auc_val, thr, sens, spec))

  # -- Bootstrap CIs --
  ci_out <- tryCatch(
    ci.coords(
      roc_obj,
      x               = "best",
      best.method     = "youden",
      best.policy     = "first",
      ret             = c("threshold", "sensitivity", "specificity"),
      boot.n          = 2000,
      boot.stratified = TRUE,
      progress        = "none"
    ),
    error = function(e) {
      message(sprintf("  WARNING [%s] ci.coords failed: %s", arm_name, e$message))
      NULL
    }
  )

  if (!is.null(ci_out)) {
    thr_lo  <- round(ci_out[["threshold"]][1],   4)
    thr_hi  <- round(ci_out[["threshold"]][3],   4)
    sens_lo <- round(ci_out[["sensitivity"]][1], 4)
    sens_hi <- round(ci_out[["sensitivity"]][3], 4)
    spec_lo <- round(ci_out[["specificity"]][1], 4)
    spec_hi <- round(ci_out[["specificity"]][3], 4)
    note    <- NA_character_
    message(sprintf("  [%s] Sens CI: %.4f-%.4f | Spec CI: %.4f-%.4f",
                    arm_name, sens_lo, sens_hi, spec_lo, spec_hi))
  } else {
    thr_lo <- thr_hi <- sens_lo <- sens_hi <- spec_lo <- spec_hi <- NA_real_
    note   <- "ci.coords failed"
  }

  data.frame(
    Cohort            = "GSE112680",
    Arm               = arm_name,
    Comparison        = comparison,
    N                 = n,
    AUC               = auc_val,
    Threshold         = round(thr,  4),
    Threshold_CI_lo   = thr_lo,
    Threshold_CI_hi   = thr_hi,
    Sensitivity       = round(sens, 4),
    Sensitivity_CI_lo = sens_lo,
    Sensitivity_CI_hi = sens_hi,
    Specificity       = round(spec, 4),
    Specificity_CI_lo = spec_lo,
    Specificity_CI_hi = spec_hi,
    Note              = note,
    stringsAsFactors  = FALSE
  )
}

# ==============================================================================
# STEP 3: Run Extraction Across All Arms
# ==============================================================================
message("Extracting Youden-optimal thresholds and bootstrap CIs...")
message("Note: 2000 bootstrap replicates per arm — expect 2-5 minutes total.")

s4_table <- do.call(rbind, lapply(arms, extract_op))

# ==============================================================================
# STEP 4: Print and Save
# ==============================================================================
message("Supplementary Table S4:")
print(s4_table, row.names = FALSE)

write.csv(
  s4_table,
  "Processed_Data/Supplementary_Table_S4_Threshold_CIs.csv",
  row.names = FALSE
)

message("================================================================")
message("Script 08 complete.")
message(sprintf("  Arms processed          : %d (GSE112680 Arms 1-3)",
                nrow(s4_table)))
message(sprintf("  Threshold method        : Youden index (J = Sens + Spec - 1)"))
message(sprintf("  Bootstrap replicates    : 2000 (stratified by class label)"))
message(sprintf("  Tie-breaking policy     : first (lowest threshold retained)"))
message(sprintf("  Output                  : Processed_Data/Supplementary_Table_S4_Threshold_CIs.csv"))
message("================================================================")
