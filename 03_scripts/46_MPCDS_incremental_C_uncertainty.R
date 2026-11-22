# =======================================================================
# Script 46 - Incremental MPCDS C-index: patient-level uncertainty audit
# =======================================================================
# Purpose
#   1. Faithfully reproduce Script 18's 20 x 5 leakage-safe nested CV
#      comparing clinical predictors with clinical predictors + MPCDS.
#   2. Save individual OUT-OF-FOLD predictions (not available in Script 18).
#   3. Resample PATIENTS, with the SAME patient resampling weights across
#      all 20 repeats, to estimate conditional sampling uncertainty in
#      the average of fold-specific paired Harrell C-index differences.
#
# Important limitations
#   - The patient-bootstrap CI is CONDITIONAL on the fitted CV models,
#     selected 47 genes, and observed fold allocations. Model refitting
#     and signature-discovery uncertainty are NOT included.
#   - The 20 repeats reuse patients and are NOT 20 independent datasets.
#   - The 2.5%/97.5% quantiles of repeat-level deltas are split-variability
#     summaries and must not be called patient-sampling confidence intervals.
#   - No previously saved analysis is changed.
#
# Required packages: glmnet, survival
# =======================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
REPEATS <- 20L
OUTER_K <- 5L
INNER_K <- 5L
BASE_SEED <- 12345L
BOOTSTRAP_B <- 500L
BOOTSTRAP_SEED <- 460046L

in_file <- function(...) file.path(ROOT, ...)
out_dir <- in_file("04_results", "mpcds_incremental_uncertainty")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("survival", "glmnet")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing R package: ", pkg)
  }
}

existing_fold_file <- in_file(
  "04_results", "mpcds", "TCGA_LUAD_clinical_MPCDS_repeated_CV.csv"
)
existing_repeat_file <- in_file(
  "04_results", "mpcds", "TCGA_LUAD_clinical_MPCDS_repeat_summary.csv"
)
existing_summary_file <- in_file(
  "04_results", "mpcds", "TCGA_LUAD_MPCDS_incremental_prediction_summary.csv"
)
meta_file <- in_file(
  "02_processed_data", "survival", "TCGA_LUAD_PCD_survival_metadata.csv"
)
expr_file <- in_file(
  "02_processed_data", "survival", "TCGA_LUAD_PCD_survival_expression.rds"
)
gene_file <- in_file(
  "04_results", "differential_expression",
  "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv"
)
all_files <- c(existing_fold_file, existing_repeat_file, existing_summary_file,
               meta_file, expr_file, gene_file)
missing_files <- all_files[!file.exists(all_files)]
if (length(missing_files)) {
  stop("Missing files:\n", paste(missing_files, collapse = "\n"))
}

original_folds <- read.csv(existing_fold_file, stringsAsFactors = FALSE)
original_repeats <- read.csv(existing_repeat_file, stringsAsFactors = FALSE)
original_summary <- read.csv(existing_summary_file, stringsAsFactors = FALSE)
meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
expr <- readRDS(expr_file)
gene_tab <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
genes <- as.character(gene_tab$gene_name)

stopifnot(length(genes) == 47L,
          !anyDuplicated(genes),
          all(genes %in% rownames(expr)),
          identical(as.character(colnames(expr)), as.character(meta$submitter_id)),
          nrow(original_folds) == REPEATS * OUTER_K,
          nrow(original_repeats) == REPEATS,
          nrow(original_summary) == 1L,
          !anyDuplicated(meta$submitter_id))

# Exactly the stage/sex coding from Script 18, without expanded categories.
stage <- as.character(meta$ajcc_pathologic_stage)
stage_group <- rep(NA_character_, length(stage))
stage_group[stage %in% c("Stage I", "Stage IA", "Stage IB")] <- "Stage I"
stage_group[stage %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
stage_group[stage %in% c("Stage IIIA", "Stage IIIB")] <- "Stage III"
stage_group[stage %in% c("Stage IV")] <- "Stage IV"

clinical_all <- data.frame(
  patient_id = as.character(meta$submitter_id),
  OS_time = as.numeric(meta$OS_time),
  OS_status = as.integer(meta$OS_status),
  age_10yr = as.numeric(meta$age_at_diagnosis_years) / 10,
  sex_at_birth = factor(as.character(meta$sex_at_birth),
                        levels = c("female", "male")),
  stage_group = factor(stage_group,
                       levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
  stringsAsFactors = FALSE
)

eligible <- complete.cases(clinical_all) &
  clinical_all$OS_time > 0 &
  clinical_all$OS_status %in% c(0L, 1L)
clinical <- clinical_all[eligible, , drop = FALSE]
x <- t(expr[genes, clinical$patient_id, drop = FALSE])
mode(x) <- "numeric"
stopifnot(ncol(x) == 47L,
          identical(rownames(x), clinical$patient_id),
          all(is.finite(x)),
          nrow(clinical) == sum(original_folds$n_test[original_folds$repeat_id == 1L]),
          sum(clinical$OS_status) == sum(original_folds$events[
            original_folds$repeat_id == 1L]))

make_stratified_folds <- function(status, k = 5L, seed = 12345L) {
  # Identical random-number use to Script 18.
  set.seed(seed)
  foldid <- integer(length(status))
  for (event_value in sort(unique(status))) {
    idx <- which(status == event_value)
    idx <- sample(idx)
    foldid[idx] <- rep(1:k, length.out = length(idx))
  }
  foldid
}

calc_c <- function(time, status, lp, weights = NULL) {
  d <- data.frame(time = time, status = status, lp = lp)
  if (is.null(weights)) {
    return(as.numeric(survival::concordance(
      survival::Surv(time, status) ~ lp,
      data = d, reverse = TRUE, timewt = "n"
    )$concordance))
  }
  d$weight <- as.integer(weights)
  d <- d[d$weight > 0L, , drop = FALSE]
  if (nrow(d) < 3L || sum(d$status) < 2L) return(NA_real_)
  tryCatch(
    as.numeric(survival::concordance(
      survival::Surv(time, status) ~ lp,
      data = d, weights = weight, reverse = TRUE, timewt = "n"
    )$concordance),
    error = function(e) NA_real_
  )
}

cat("\n===== SCRIPT 46: INPUT CHECK =====\n")
cat("Complete-case clinical patients:", nrow(clinical),
    "| events:", sum(clinical$OS_status), "\n")
cat("PCD genes:", ncol(x), "| repeats:", REPEATS,
    "| outer folds:", OUTER_K, "\n")

pred_parts <- vector("list", REPEATS * OUTER_K)
fold_parts <- vector("list", REPEATS * OUTER_K)
ix <- 0L

for (r in seq_len(REPEATS)) {
  repeat_seed <- BASE_SEED + (r - 1L) * 1000L
  foldid <- make_stratified_folds(clinical$OS_status, OUTER_K, repeat_seed)
  for (k in seq_len(OUTER_K)) {
    ix <- ix + 1L
    itest <- which(foldid == k)
    itrain <- which(foldid != k)
    tr <- clinical[itrain, , drop = FALSE]
    te <- clinical[itest, , drop = FALSE]
    xtr <- x[itrain, , drop = FALSE]
    xte <- x[itest, , drop = FALSE]
    inner_folds <- make_stratified_folds(tr$OS_status, INNER_K, repeat_seed + k)

    ridge_cv <- glmnet::cv.glmnet(
      x = xtr,
      y = survival::Surv(tr$OS_time, tr$OS_status),
      family = "cox", alpha = 0,
      foldid = inner_folds,
      type.measure = "C", standardize = TRUE,
      nlambda = 100, cox.ties = "breslow"
    )
    ridge_fit <- glmnet::glmnet(
      x = xtr,
      y = survival::Surv(tr$OS_time, tr$OS_status),
      family = "cox", alpha = 0,
      lambda = ridge_cv$lambda.min,
      standardize = TRUE, cox.ties = "breslow"
    )
    tr$MPCDS_cv <- as.numeric(predict(ridge_fit, newx = xtr, type = "link"))
    te$MPCDS_cv <- as.numeric(predict(ridge_fit, newx = xte, type = "link"))
    train_mean <- mean(tr$MPCDS_cv)
    train_sd <- stats::sd(tr$MPCDS_cv)
    if (!is.finite(train_sd) || train_sd <= 0) stop("Zero SD in CV fold")
    tr$MPCDS_cv_z <- (tr$MPCDS_cv - train_mean) / train_sd
    te$MPCDS_cv_z <- (te$MPCDS_cv - train_mean) / train_sd

    clinical_fit <- survival::coxph(
      survival::Surv(OS_time, OS_status) ~ age_10yr + sex_at_birth + stage_group,
      data = tr, ties = "efron"
    )
    combined_fit <- survival::coxph(
      survival::Surv(OS_time, OS_status) ~
        MPCDS_cv_z + age_10yr + sex_at_birth + stage_group,
      data = tr, ties = "efron"
    )
    cl_lp <- as.numeric(predict(clinical_fit, newdata = te, type = "lp"))
    co_lp <- as.numeric(predict(combined_fit, newdata = te, type = "lp"))
    cl_c <- calc_c(te$OS_time, te$OS_status, cl_lp)
    co_c <- calc_c(te$OS_time, te$OS_status, co_lp)

    pred_parts[[ix]] <- data.frame(
      patient_id = te$patient_id,
      repeat_id = r, outer_fold = k,
      OS_time = te$OS_time, OS_status = te$OS_status,
      clinical_lp = cl_lp, combined_lp = co_lp,
      stringsAsFactors = FALSE
    )
    fold_parts[[ix]] <- data.frame(
      repeat_id = r, outer_fold = k,
      n_test = length(itest), events = sum(te$OS_status),
      lambda = ridge_cv$lambda.min,
      clinical_C = cl_c, combined_C = co_c,
      delta_C = co_c - cl_c
    )
  }
  v <- do.call(rbind, fold_parts[((r-1L)*OUTER_K + 1L):(r*OUTER_K)])
  cat(sprintf("Repeat %2d/%d reproduced: clinical C=%.4f | combined C=%.4f | delta=%.5f\n",
              r, REPEATS, mean(v$clinical_C), mean(v$combined_C), mean(v$delta_C)))
}

pred <- do.call(rbind, pred_parts)
fold <- do.call(rbind, fold_parts)
rownames(pred) <- NULL
rownames(fold) <- NULL
stopifnot(nrow(pred) == REPEATS * nrow(clinical),
          !anyDuplicated(pred[, c("patient_id", "repeat_id")]),
          all(is.finite(pred$clinical_lp)), all(is.finite(pred$combined_lp)))

# Reproduce Script 18 exactly BEFORE proceeding with uncertainty analysis.
key <- paste(fold$repeat_id, fold$outer_fold)
old_key <- paste(original_folds$repeat_id, original_folds$outer_fold)
ref <- original_folds[match(key, old_key), , drop = FALSE]
if (anyNA(ref$clinical_C)) stop("Could not match all original CV folds")

check_cols <- c("n_test", "events", "lambda", "clinical_C", "combined_C", "delta_C")
audit <- data.frame(repeat_id = fold$repeat_id, outer_fold = fold$outer_fold)
for (nm in check_cols) {
  audit[[paste0("original_", nm)]] <- ref[[nm]]
  audit[[paste0("reproduced_", nm)]] <- fold[[nm]]
  audit[[paste0("absolute_difference_", nm)]] <-
    abs(as.numeric(ref[[nm]]) - as.numeric(fold[[nm]]))
}
write.csv(audit, file.path(out_dir, "fold_reproduction_audit.csv"), row.names = FALSE)
for (nm in check_cols) {
  tol <- if (nm %in% c("n_test", "events")) 0 else 1e-6
  if (max(audit[[paste0("absolute_difference_", nm)]], na.rm = TRUE) > tol) {
    stop("REPRODUCTION FAILED for ", nm,
         ". Inspect fold_reproduction_audit.csv. No bootstrap was performed.")
  }
}

repeats <- do.call(rbind, lapply(split(fold, fold$repeat_id), function(a) {
  data.frame(
    repeat_id = a$repeat_id[1],
    mean_clinical_C = mean(a$clinical_C),
    mean_combined_C = mean(a$combined_C),
    mean_delta_C = mean(a$delta_C)
  )
}))
repeats <- repeats[order(repeats$repeat_id), , drop = FALSE]
old_repeats <- original_repeats[match(repeats$repeat_id,
                                       original_repeats$repeat_id), , drop = FALSE]
stopifnot(max(abs(repeats$mean_delta_C - old_repeats$mean_delta_C)) < 1e-6,
          abs(mean(repeats$mean_delta_C) - original_summary$overall_delta_C[1]) < 1e-6)
write.csv(repeats, file.path(out_dir, "reproduced_repeat_summary.csv"),
          row.names = FALSE)
write.csv(pred, file.path(out_dir, "out_of_fold_patient_predictions.csv"),
          row.names = FALSE)

cat("\n===== ORIGINAL SCRIPT 18 REPRODUCED =====\n")
cat(sprintf("Clinical C=%.7f | combined C=%.7f | delta=%.7f\n",
            mean(fold$clinical_C), mean(fold$combined_C), mean(fold$delta_C)))
cat("Split-variability 2.5–97.5 percentiles:",
    sprintf("%.5f to %.5f\n",
            stats::quantile(repeats$mean_delta_C, .025),
            stats::quantile(repeats$mean_delta_C, .975)))
cat("These split percentiles are NOT a patient-sampling CI.\n")

# Patient-cluster bootstrap of held-out predictions, NO model refitting.
# One set of patient multiplicities is used in all 20 repeats. Evaluation
# stays within each original outer test fold (avoids comparing uncalibrated
# linear predictors from different Cox fits across different folds).
patient_index <- seq_len(nrow(clinical))
names(patient_index) <- clinical$patient_id
pred$patient_idx <- unname(patient_index[pred$patient_id])
stopifnot(!anyNA(pred$patient_idx))
fold_groups <- split(pred, paste(pred$repeat_id, pred$outer_fold, sep = ":"))

cat("\n===== PATIENT-LEVEL BOOTSTRAP (CONDITIONAL ON CV FITS) =====\n")
cat("Resamples:", BOOTSTRAP_B, "| patient clusters:", nrow(clinical), "\n")
set.seed(BOOTSTRAP_SEED)
boot <- numeric(BOOTSTRAP_B)
for (b in seq_len(BOOTSTRAP_B)) {
  sampled <- sample.int(nrow(clinical), size = nrow(clinical), replace = TRUE)
  counts <- tabulate(sampled, nbins = nrow(clinical))
  fold_delta <- vapply(fold_groups, function(d) {
    w <- counts[d$patient_idx]
    a <- calc_c(d$OS_time, d$OS_status, d$clinical_lp, w)
    z <- calc_c(d$OS_time, d$OS_status, d$combined_lp, w)
    z - a
  }, numeric(1))
  if (any(!is.finite(fold_delta))) {
    stop("Bootstrap ", b, " has invalid fold C-index. Retry with a different seed.")
  }
  boot[b] <- mean(fold_delta)
  if (b %% 50L == 0L || b == BOOTSTRAP_B) {
    cat("Bootstrap", b, "of", BOOTSTRAP_B, "complete\n")
  }
}

boot_df <- data.frame(bootstrap_id = seq_len(BOOTSTRAP_B), delta_C = boot)
write.csv(boot_df, file.path(out_dir, "patient_bootstrap_delta_C.csv"),
          row.names = FALSE)
ci <- as.numeric(stats::quantile(boot, probs = c(.025, .975), na.rm = TRUE))
summary_df <- data.frame(
  n_patients = nrow(clinical), events = sum(clinical$OS_status),
  repeats = REPEATS, outer_folds = OUTER_K,
  original_clinical_C = original_summary$overall_clinical_C[1],
  original_combined_C = original_summary$overall_combined_C[1],
  original_delta_C = original_summary$overall_delta_C[1],
  repeated_CV_split_Q025 = original_summary$Q025_delta_C[1],
  repeated_CV_split_Q975 = original_summary$Q975_delta_C[1],
  bootstrap_resamples = BOOTSTRAP_B,
  patient_bootstrap_conditional_CI_lower = ci[1],
  patient_bootstrap_conditional_CI_upper = ci[2],
  bootstrap_fraction_positive = mean(boot > 0),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "incremental_uncertainty_summary.csv"),
          row.names = FALSE)
notes <- c(
  "SCRIPT 46: MPCDS incremental C-index patient-bootstrap sensitivity",
  "Exact reproduction required: 20 repeated five-fold nested CV with inner ridge tuning; same Script 18 seeds and predictors.",
  "Target statistic: arithmetic average of the 100 outer-fold paired C-index differences, identical to Script 18.",
  "A patient is resampled once per bootstrap iteration; their multiplicity is reused across ALL CV repeats.",
  "Within each fold, concordance is evaluated with patient case weights on OUT-OF-FOLD predictions.",
  "Bootstrap interval: percentile interval conditional on observed dataset, CV model fits, gene set and fold allocation.",
  "This is NOT a fully model-refitted bootstrap interval, does not reflect signature discovery or refitting uncertainty, and is not external validation.",
  "CV repeat 2.5th/97.5th percentiles describe fold-partition variability and are not patient-sampling CIs.",
  "No previous Script 18 outputs were altered."
)
writeLines(notes, file.path(out_dir, "methodological_notes.txt"))

cat("\n===== SCRIPT 46 FINISHED =====\n")
print(summary_df, row.names = FALSE)
cat("Conditional patient-bootstrap percentile 95% CI:",
    sprintf("%.5f to %.5f\n", ci[1], ci[2]))
cat("Original ΔC:", sprintf("%.5f\n", original_summary$overall_delta_C[1]))
cat("NOTE: bootstrap does not refit models or redo gene selection.\n")
cat("Output folder:", out_dir, "\n")
