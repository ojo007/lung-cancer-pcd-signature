# ============================================================================
# SCRIPT 48: GEO LUAD MPCDS 3-/5-year cumulative/dynamic AUC uncertainty
# ============================================================================
# Purpose: Reproduce Script 22, then quantify fixed-score time-dependent AUC
# variability with 500 patient-level nonparametric bootstrap resamples.
#
# Critical interpretation:
# * This is a conditional bootstrap: the 47-gene MPCDS fit and the external
#   cohort-wise normalization are fixed. We do NOT refit, restandardize,
#   reselect genes, or bootstrap external cohort membership.
# * The original GSE72094 5-year analysis has only 5 known survivors beyond
#   60 months. Bootstrapping cannot repair that poor late follow-up support.
# * No prior outputs are overwritten.
# ============================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
COHORTS <- c("GSE68465", "GSE72094", "GSE31210", "GSE50081")
HORIZONS <- c(36, 60)  # months, 3 and 5 years
TIME_ROC_HORIZONS <- c(12, 36, 60)  # match Script 22 precisely
N_BOOT <- 500L
BASE_SEED <- 480000L

OUT <- file.path(ROOT, "04_results", "mpcds_5yr_AUC_uncertainty")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("timeROC", quietly = TRUE) ||
    !requireNamespace("survival", quietly = TRUE)) {
  stop("Packages timeROC and survival are required. Script 22 used these packages.")
}

original_summary_file <- file.path(
  ROOT, "04_results", "mpcds", "time_dependent_ROC",
  "LUAD_MPCDS_external_GEO_time_dependent_ROC_summary.csv"
)
if (!file.exists(original_summary_file)) {
  stop("Missing Script 22 ROC summary: ", original_summary_file)
}
original <- read.csv(original_summary_file, stringsAsFactors = FALSE,
                     check.names = FALSE)
needed_original <- c("cohort", "time_value", "AUC", "CI_lower", "CI_upper")
if (!all(needed_original %in% names(original))) {
  stop("Script 22 summary is missing: ",
       paste(setdiff(needed_original, names(original)), collapse = ", "))
}

get_scores_file <- function(cohort) {
  in_cohort <- file.path(ROOT, "04_results", "mpcds", "external_validation",
                         cohort, paste0(cohort, "_MPCDS_patient_scores.csv"))
  if (file.exists(in_cohort)) return(in_cohort)
  if (cohort == "GSE68465") {
    legacy <- file.path(ROOT, "04_results", "mpcds", "external_validation",
                        "GSE68465_MPCDS_patient_scores.csv")
    if (file.exists(legacy)) return(legacy)
  }
  stop("Missing Script 21 patient scores for ", cohort, ": ", in_cohort)
}

# Same marginal-censoring cumulative/dynamic ROC as Script 22.
calculate_roc <- function(d) {
  out <- suppressWarnings(timeROC::timeROC(
    T = d$OS_time,
    delta = d$OS_status,
    marker = d$MPCDS_z,
    cause = 1,
    weighting = "marginal",
    times = TIME_ROC_HORIZONS,
    iid = FALSE
  ))
  which_horizon <- match(HORIZONS, as.numeric(out$times))
  if (anyNA(which_horizon)) stop("timeROC did not return expected time points.")
  aucs <- as.numeric(out$AUC[which_horizon])
  if (length(aucs) != length(HORIZONS) || any(!is.finite(aucs)) ||
      any(aucs < 0 | aucs > 1)) {
    stop("Nonfinite/out-of-range AUC from timeROC.")
  }
  aucs
}

horizon_counts <- function(d, t) {
  c(events_by_horizon = sum(d$OS_status == 1L & d$OS_time <= t),
    survivors_beyond_horizon = sum(d$OS_time > t),
    censored_by_horizon = sum(d$OS_status == 0L & d$OS_time <= t))
}

all_results <- vector("list", length(COHORTS))
all_boot <- vector("list", length(COHORTS))

for (i in seq_along(COHORTS)) {
  cohort <- COHORTS[i]
  cat("\n=============== ", cohort, " ===============\n", sep = "")
  d <- read.csv(get_scores_file(cohort), stringsAsFactors = FALSE,
                check.names = FALSE)
  need <- c("sample_id", "OS_time", "OS_status", "MPCDS_z")
  if (!all(need %in% names(d))) {
    stop(cohort, ": missing columns ",
         paste(setdiff(need, names(d)), collapse = ", "))
  }
  d <- d[, need, drop = FALSE]
  d$sample_id <- as.character(d$sample_id)
  d$OS_time <- as.numeric(d$OS_time)
  d$OS_status <- as.integer(as.character(d$OS_status))
  d$MPCDS_z <- as.numeric(d$MPCDS_z)
  if (anyDuplicated(d$sample_id) || anyNA(d$sample_id) ||
      any(!nzchar(d$sample_id)) ||
      any(!is.finite(d$OS_time)) || any(d$OS_time <= 0) ||
      anyNA(d$OS_status) || any(!d$OS_status %in% c(0L, 1L)) ||
      any(!is.finite(d$MPCDS_z))) {
    stop(cohort, ": invalid patient-level survival data")
  }

  baseline <- calculate_roc(d)
  original_rows <- original[original$cohort == cohort &
                            original$time_value %in% HORIZONS, , drop = FALSE]
  if (nrow(original_rows) != 2L ||
      anyDuplicated(original_rows$time_value)) {
    stop(cohort, ": missing/duplicated original 3-/5-year AUC rows")
  }
  original_rows <- original_rows[match(HORIZONS, original_rows$time_value), , drop = FALSE]
  max_diff <- max(abs(baseline - original_rows$AUC))
  if (!is.finite(max_diff) || max_diff > 1e-7) {
    stop(cohort, ": FAILED Script 22 AUC reproduction; max delta=", max_diff)
  }
  cat("Original AUCs reproduced: 3-year=", sprintf("%.6f", baseline[1]),
      " | 5-year=", sprintf("%.6f", baseline[2]), "\n", sep = "")
  n <- nrow(d)
  five_counts <- horizon_counts(d, 60)
  cat("Patients: ", n, " | deaths: ", sum(d$OS_status),
      " | 5-year events: ", five_counts["events_by_horizon"],
      " | known survivors: ", five_counts["survivors_beyond_horizon"],
      " | censored by five years: ", five_counts["censored_by_horizon"], "\n", sep = "")

  set.seed(BASE_SEED + i * 10000L)
  boot <- vector("list", N_BOOT * length(HORIZONS))
  pos <- 0L
  for (b in seq_len(N_BOOT)) {
    ids <- sample.int(n, size = n, replace = TRUE)
    sample_d <- d[ids, , drop = FALSE]
    supports <- lapply(HORIZONS, function(t) horizon_counts(sample_d, t))
    # Fit once per bootstrap sample to retain matched 3- and 5-year estimates.
    if (all(vapply(supports, function(x) x["events_by_horizon"] >= 2L &&
                   x["survivors_beyond_horizon"] >= 2L, logical(1)))) {
      roc_attempt <- tryCatch(calculate_roc(sample_d), error = function(e) e)
      if (inherits(roc_attempt, "error")) {
        vals <- rep(NA_real_, length(HORIZONS))
        error_msg <- conditionMessage(roc_attempt)
      } else {
        vals <- roc_attempt
        error_msg <- ""
      }
    } else {
      vals <- rep(NA_real_, length(HORIZONS))
      error_msg <- "Insufficient cases/known survivors at 3 or 5 years"
    }
    for (j in seq_along(HORIZONS)) {
      pos <- pos + 1L
      supp <- supports[[j]]
      boot[[pos]] <- data.frame(
        cohort = cohort, bootstrap_id = b, horizon_months = HORIZONS[j],
        auc = vals[j],
        events_by_horizon = as.integer(supp["events_by_horizon"]),
        survivors_beyond_horizon = as.integer(supp["survivors_beyond_horizon"]),
        censored_by_horizon = as.integer(supp["censored_by_horizon"]),
        failure_reason = error_msg, stringsAsFactors = FALSE
      )
    }
    if (b %% 100L == 0L) cat("Bootstrap ", b, "/", N_BOOT, " complete\n", sep = "")
  }
  all_boot[[i]] <- do.call(rbind, boot)

  rows <- vector("list", length(HORIZONS))
  for (j in seq_along(HORIZONS)) {
    t <- HORIZONS[j]
    sub <- all_boot[[i]][all_boot[[i]]$horizon_months == t, , drop = FALSE]
    valid <- is.finite(sub$auc)
    n_valid <- sum(valid)
    # Report CIs only when >=90% of bootstrap replicates yielded finite AUCs.
    q <- if (n_valid >= ceiling(0.9 * N_BOOT))
      unname(stats::quantile(sub$auc[valid], c(0.025, 0.975)))
    else c(NA_real_, NA_real_)
    supp <- horizon_counts(d, t)
    # Descriptive sparse-follow-up flag, not a validated statistical threshold.
    followup_flag <- if (supp["survivors_beyond_horizon"] < 10L)
      "VERY_SPARSE_KNOWN_SURVIVORS; DO_NOT_OVERINTERPRET"
    else if (supp["survivors_beyond_horizon"] < 20L)
      "LIMITED_KNOWN_SURVIVORS"
    else "NO_EXTREME_SPARSE_SURVIVOR_FLAG"
    rows[[j]] <- data.frame(
      cohort = cohort,
      horizon_months = t,
      n_patients = n,
      total_deaths = sum(d$OS_status),
      events_by_horizon = as.integer(supp["events_by_horizon"]),
      survivors_beyond_horizon = as.integer(supp["survivors_beyond_horizon"]),
      censored_by_horizon = as.integer(supp["censored_by_horizon"]),
      original_AUC = baseline[j],
      original_timeROC_CI_lower = original_rows$CI_lower[j],
      original_timeROC_CI_upper = original_rows$CI_upper[j],
      n_bootstrap = N_BOOT,
      n_bootstrap_valid = n_valid,
      percent_bootstrap_valid = 100 * n_valid / N_BOOT,
      bootstrap_mean_AUC = if (n_valid > 0) mean(sub$auc[valid]) else NA_real_,
      conditional_bootstrap_CI_lower = q[1],
      conditional_bootstrap_CI_upper = q[2],
      min_bootstrap_known_survivors = min(sub$survivors_beyond_horizon),
      fraction_bootstrap_fewer_than_10_survivors =
        mean(sub$survivors_beyond_horizon < 10L),
      followup_flag = followup_flag,
      stringsAsFactors = FALSE
    )
  }
  all_results[[i]] <- do.call(rbind, rows)
  print(all_results[[i]][, c("cohort", "horizon_months", "original_AUC",
                           "survivors_beyond_horizon", "censored_by_horizon",
                           "n_bootstrap_valid", "conditional_bootstrap_CI_lower",
                           "conditional_bootstrap_CI_upper", "followup_flag")],
        row.names = FALSE)
}

summary <- do.call(rbind, all_results)
boot_results <- do.call(rbind, all_boot)
rownames(summary) <- rownames(boot_results) <- NULL

stopifnot(nrow(summary) == 8L,
          nrow(boot_results) == length(COHORTS) * N_BOOT * length(HORIZONS),
          all(summary$original_AUC >= 0 & summary$original_AUC <= 1))

write.csv(summary, file.path(OUT, "GEO_MPCDS_3yr_5yr_AUC_uncertainty_summary.csv"),
          row.names = FALSE)
write.csv(boot_results,
          file.path(OUT, "GEO_MPCDS_3yr_5yr_AUC_patient_bootstrap.csv"),
          row.names = FALSE)
writeLines(c(
  "SCRIPT 48 — Fixed-MPCDS external GEO cumulative/dynamic AUC sensitivity",
  "Source: original Script 21 patient-level score files and Script 22 timeROC summary.",
  "3/5-year AUCs reproduce original Script 22 (marginal-censoring IPCW, cumulative/dynamic ROC).",
  "Each of the four cohorts: 500 individual patient bootstrap resamples with replacement, fixed scores.",
  "Percentile intervals use finite bootstrap AUCs when at least 90% of resamples are valid.",
  "No model retraining, gene reselection, recalibration, or recalculation of cohort-wide z-score is performed.",
  "Bootstrap intervals are conditional on fixed scores and the observed cohort's follow-up/censoring patterns.",
  "GSE72094 has only FIVE observed survivors beyond 60 months: regardless of",
  "the calculated interval, its 5-year AUC must be described as unstable / poorly supported.",
  "The <10-survivor flag is a descriptive caution, not a statistically validated cutoff.",
  "The bootstrap excludes invalid replicates and reports their frequency; do not impute missing AUCs.",
  "These analyses do not demonstrate clinical utility or independent external prospective deployment.",
  "Existing Script 21 and 22 outputs remain unchanged."
), file.path(OUT, "methods_and_interpretation.txt"))

cat("\n============= SCRIPT 48 FINISHED =============\n")
print(summary[, c("cohort", "horizon_months", "original_AUC", "n_bootstrap_valid",
                  "conditional_bootstrap_CI_lower", "conditional_bootstrap_CI_upper",
                  "survivors_beyond_horizon", "followup_flag")], row.names = FALSE)
cat("Outputs: ", OUT, "\n", sep = "")
cat("Interpretation caution: GSE72094 5-year AUC has only five known survivors.\n")
