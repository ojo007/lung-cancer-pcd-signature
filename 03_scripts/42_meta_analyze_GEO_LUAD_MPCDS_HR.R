# ==============================================================================
# SCRIPT 42 — RANDOM-EFFECTS META-ANALYSIS OF EXTERNAL LUAD MPCDS HAZARD RATIOS
# Kankarofi et al. — Cell Death Discovery scientific corrections
#
# Purpose:
#   Independently audit Script 41's unadjusted estimates against Script 21;
#   pool unadjusted, age/sex-adjusted, and stage-adjusted log-hazard ratios
#   from four external LUAD cohorts, using REML random effects.
#   Report Cochran's Q, Q-based I-squared, REML tau-squared, conventional
#   random-effects CI and conservative ad hoc Hartung–Knapp t-based CI.
#   Perform leave-one-cohort-out analyses for each model specification.
#
# Important:
#   - NO model fitting/retraining and NO patient-level outcome reprocessing.
#   - NO pooling of C-indices; effects are log-HRs per within-cohort SD MPCDS.
#   - GSE68465 uses pT/pN proxies, other cohorts use stage I versus II-IV.
#   - Only four cohorts; random-effects heterogeneity is imprecise.
#   - GSE31210 has just 35 deaths. PH assumptions also warrant review.
#   - Existing Script 41 and Script 21 outputs are never overwritten.
#   - Only base R packages (stats, utils) are used.
#
# Source data:
#   04_results/mpcds_clinically_adjusted_external/
#      LUAD_MPCDS_clinical_Cox_by_cohort.csv
#   04_results/mpcds/external_validation/<cohort>/
#      <cohort>_MPCDS_validation_summary.csv
#
# Outputs:
#   04_results/mpcds_clinically_adjusted_external/meta_analysis/
#      LUAD_MPCDS_REML_meta_summary.csv
#      LUAD_MPCDS_meta_cohort_contributions.csv
#      LUAD_MPCDS_meta_leave_one_out.csv
#      LUAD_MPCDS_script21_unadjusted_audit.csv
#      LUAD_MPCDS_meta_methods_and_cautions.txt
# ==============================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
COHORTS <- c("GSE68465", "GSE72094", "GSE31210", "GSE50081")
MODELS <- c("unadjusted", "age_sex_adjusted", "stage_adjusted")
INPUT <- file.path(ROOT, "04_results", "mpcds_clinically_adjusted_external",
                   "LUAD_MPCDS_clinical_Cox_by_cohort.csv")
OUT <- file.path(ROOT, "04_results", "mpcds_clinically_adjusted_external",
                 "meta_analysis")

if (!file.exists(INPUT)) stop("Script 41 summary is missing: ", INPUT)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

fits <- read.csv(INPUT, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("cohort", "model", "status", "n_analysed", "events_analysed",
            "log_HR", "log_HR_SE", "HR_per_MPCDS_SD", "CI95_lower",
            "CI95_upper")
missing <- setdiff(needed, names(fits))
if (length(missing)) stop("Script 41 summary missing fields: ",
                          paste(missing, collapse = ", "))
if (anyDuplicated(paste(fits$cohort, fits$model, sep = "|"))) {
  stop("Duplicate cohort/model rows in Script 41 summary")
}

selected <- fits[fits$cohort %in% COHORTS & fits$model %in% MODELS, , drop = FALSE]
selected <- selected[match(as.vector(outer(COHORTS, MODELS, paste, sep = "|")),
                           paste(selected$cohort, selected$model, sep = "|")), , drop = FALSE]
# Match can produce NA rows if a combination is absent. Reject rather than pool silently.
if (nrow(selected) != length(COHORTS) * length(MODELS) ||
    anyNA(selected$cohort) || anyNA(selected$model)) {
  stop("All four cohorts must have all three planned model specifications")
}
if (any(selected$status != "OK")) {
  bad <- selected[selected$status != "OK", c("cohort", "model", "status")]
  print(bad)
  stop("A required model failed or was skipped; cannot conduct planned meta-analysis")
}
num_fields <- c("log_HR", "log_HR_SE", "HR_per_MPCDS_SD",
                "CI95_lower", "CI95_upper", "n_analysed", "events_analysed")
for (field in num_fields) {
  selected[[field]] <- suppressWarnings(as.numeric(selected[[field]]))
}
if (any(!is.finite(as.matrix(selected[, num_fields, drop = FALSE]))) ||
    any(selected$log_HR_SE <= 0) ||
    any(selected$events_analysed < 1) ||
    any(selected$n_analysed < selected$events_analysed)) {
  stop("Invalid/nonfinite model estimates, SEs, patients, or event counts")
}
if (any(abs(log(selected$HR_per_MPCDS_SD) - selected$log_HR) > 1e-6)) {
  stop("Hazard ratios are inconsistent with log hazard ratios")
}
if (any(abs(log(selected$CI95_lower) -
            (selected$log_HR - 1.96 * selected$log_HR_SE)) > 1e-5) ||
    any(abs(log(selected$CI95_upper) -
            (selected$log_HR + 1.96 * selected$log_HR_SE)) > 1e-5)) {
  stop("Script 41 confidence limits inconsistent with coefficients and SE")
}

# -------- Cross-check unadjusted HRs against existing Script 21 results --------
check_rows <- lapply(COHORTS, function(cohort) {
  path <- file.path(ROOT, "04_results", "mpcds", "external_validation",
                    cohort, paste0(cohort, "_MPCDS_validation_summary.csv"))
  if (!file.exists(path)) stop("Missing Script 21 validation summary: ", path)
  prior <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  columns <- c("continuous_HR_per_SD", "n_OS_samples", "n_events")
  if (nrow(prior) != 1L || !all(columns %in% names(prior))) {
    stop("Unexpected Script 21 schema for ", cohort,
         "; inspect original summary before pooling")
  }
  cur <- selected[selected$cohort == cohort & selected$model == "unadjusted", ]
  old_hr <- as.numeric(prior$continuous_HR_per_SD)
  new_hr <- cur$HR_per_MPCDS_SD
  old_n <- as.integer(prior$n_OS_samples)
  old_events <- as.integer(prior$n_events)
  ok <- is.finite(old_hr) && abs(old_hr - new_hr) <= 0.001 &&
    old_n == cur$n_analysed && old_events == cur$events_analysed
  data.frame(cohort = cohort, script21_HR = old_hr, script41_HR = new_hr,
             absolute_HR_difference = abs(old_hr - new_hr),
             script21_n = old_n, script41_n = cur$n_analysed,
             script21_events = old_events, script41_events = cur$events_analysed,
             verified = ok, stringsAsFactors = FALSE)
})
audit <- do.call(rbind, check_rows)
row.names(audit) <- NULL
cat("\n===== UNADJUSTED VALIDATION REPRODUCIBILITY =====\n")
print(audit, row.names = FALSE, digits = 6)
if (!all(audit$verified)) {
  stop("The Script 41 unadjusted estimates do not reproduce Script 21. ",
       "Investigate the discrepancy before conducting meta-analysis.")
}
write.csv(audit, file.path(OUT, "LUAD_MPCDS_script21_unadjusted_audit.csv"),
          row.names = FALSE)

# -------- Base R REML estimator / random-effects inference --------
meta_reml <- function(x) {
  y <- as.numeric(x$log_HR)
  se <- as.numeric(x$log_HR_SE)
  vi <- se^2
  k <- length(y)
  if (k < 3L || any(!is.finite(y)) || any(!is.finite(vi)) || any(vi <= 0)) {
    stop("At least three valid log-HRs and positive sampling variances are required")
  }
  # Restricted negative log-likelihood, up to a tau-independent constant.
  # The one-dimensional optimisation avoids any external package dependency.
  nll <- function(tau2) {
    v <- vi + tau2
    w <- 1 / v
    mu <- sum(w * y) / sum(w)
    0.5 * (sum(log(v)) + log(sum(w)) + sum(w * (y - mu)^2))
  }
  upper <- max(1, 10 * diff(range(y))^2, 50 * max(vi))
  sol <- stats::optimize(nll, interval = c(0, upper), tol = 1e-11)
  tau2 <- if (nll(0) <= sol$objective) 0 else sol$minimum
  if (abs(tau2 - upper) / upper < 1e-4) {
    warning("REML heterogeneity estimate near optimisation upper bound")
  }
  w <- 1 / (vi + tau2)
  mu <- sum(w * y) / sum(w)
  se_mu <- sqrt(1 / sum(w))
  z <- stats::qnorm(0.975)
  wald_lo <- mu - z * se_mu
  wald_hi <- mu + z * se_mu
  # Hartung–Knapp method with ad hoc safeguard (never shrink variance).
  # Particularly important with k = 4; no CI may exclude zero by assumption.
  q_re <- sum(w * (y - mu)^2)
  hk_scale <- max(1, q_re / (k - 1))
  se_hk <- se_mu * sqrt(hk_scale)
  tcrit <- stats::qt(0.975, df = k - 1)
  hk_lo <- mu - tcrit * se_hk
  hk_hi <- mu + tcrit * se_hk
  p_hk <- 2 * stats::pt(-abs(mu / se_hk), df = k - 1)
  # Cochran's Q and conventional Q-based I-squared.
  wf <- 1 / vi
  muf <- sum(wf * y) / sum(wf)
  Q <- sum(wf * (y - muf)^2)
  Q_p <- stats::pchisq(Q, df = k - 1, lower.tail = FALSE)
  I2 <- if (Q > 0) 100 * max(0, (Q - (k - 1)) / Q) else 0
  data.frame(
    cohorts = k,
    analysed_patients_sum = sum(x$n_analysed),
    events_sum = sum(x$events_analysed),
    pooled_log_HR = mu,
    pooled_HR = exp(mu),
    tau2_REML = tau2,
    tau_REML = sqrt(tau2),
    Cochran_Q = Q,
    Q_df = k - 1,
    Q_p = Q_p,
    I2_percent_Q_based = I2,
    conventional_REML_HR_CI_lower = exp(wald_lo),
    conventional_REML_HR_CI_upper = exp(wald_hi),
    HK_ad_hoc_HR_CI_lower = exp(hk_lo),
    HK_ad_hoc_HR_CI_upper = exp(hk_hi),
    HK_ad_hoc_p = p_hk,
    HK_scale_multiplier = hk_scale,
    stringsAsFactors = FALSE
  )
}

meta_results <- list()
contrib <- list()
loo <- list()

for (model_name in MODELS) {
  dat <- selected[selected$model == model_name, , drop = FALSE]
  dat <- dat[match(COHORTS, dat$cohort), , drop = FALSE]
  cat("\n============= ", model_name, " =============\n", sep = "")
  print(dat[, c("cohort", "n_analysed", "events_analysed",
                 "HR_per_MPCDS_SD", "CI95_lower", "CI95_upper")],
        row.names = FALSE, digits = 5)
  result <- meta_reml(dat)
  result$model <- model_name
  result <- result[, c("model", setdiff(names(result), "model"))]
  meta_results[[model_name]] <- result
  weights <- 1 / (dat$log_HR_SE^2 + result$tau2_REML)
  dat$random_effects_weight_percent <- 100 * weights / sum(weights)
  dat$tau2_REML <- result$tau2_REML
  contrib[[model_name]] <- dat[, c("model", "cohort", "n_analysed",
                                   "events_analysed", "log_HR", "log_HR_SE",
                                   "HR_per_MPCDS_SD", "CI95_lower", "CI95_upper",
                                   "tau2_REML", "random_effects_weight_percent")]
  for (omit in COHORTS) {
    sub <- dat[dat$cohort != omit, , drop = FALSE]
    tmp <- meta_reml(sub)
    tmp$model <- model_name
    tmp$omitted_cohort <- omit
    tmp <- tmp[, c("model", "omitted_cohort",
                   setdiff(names(tmp), c("model", "omitted_cohort")))]
    loo[[paste(model_name, omit, sep = "__")]] <- tmp
  }
  cat(sprintf("REML pooled HR %.4f | conventional 95%% CI %.4f–%.4f\n",
              result$pooled_HR, result$conventional_REML_HR_CI_lower,
              result$conventional_REML_HR_CI_upper))
  cat(sprintf("Conservative HK 95%% CI %.4f–%.4f | p %.4g\n",
              result$HK_ad_hoc_HR_CI_lower, result$HK_ad_hoc_HR_CI_upper,
              result$HK_ad_hoc_p))
  cat(sprintf("I-squared %.1f%% | REML tau-squared %.5f | Q p %.4g\n",
              result$I2_percent_Q_based, result$tau2_REML, result$Q_p))
}

meta_df <- do.call(rbind, meta_results)
contrib_df <- do.call(rbind, contrib)
loo_df <- do.call(rbind, loo)
row.names(meta_df) <- NULL
row.names(contrib_df) <- NULL
row.names(loo_df) <- NULL
write.csv(meta_df,
          file.path(OUT, "LUAD_MPCDS_REML_meta_summary.csv"), row.names = FALSE)
write.csv(contrib_df,
          file.path(OUT, "LUAD_MPCDS_meta_cohort_contributions.csv"),
          row.names = FALSE)
write.csv(loo_df,
          file.path(OUT, "LUAD_MPCDS_meta_leave_one_out.csv"),
          row.names = FALSE)

notes <- c(
  "Kankarofi et al. — External LUAD MPCDS meta-analysis, Script 42",
  "Cohorts: GSE68465, GSE72094, GSE31210, GSE50081.",
  "Input: validated per-cohort Cox log-HR and model-based SE (Script 41).",
  "Each effect represents the HR per one within-cohort SD of MPCDS.",
  "Estimand: cross-cohort mean of study-specific associations; not a causal effect.",
  "Primary model specification: stage_adjusted, with cohort-specific clinical staging.",
  "GSE68465 adjusts for age, sex, pT2+ and pN1+; other cohorts adjust for",
  "age, sex and stage I versus II-IV. Stage covariates are not fully harmonised.",
  "Unadjusted and age-sex adjusted models are supporting sensitivity summaries.",
  "Normal-normal random-effects inverse-variance meta-analysis, tau2 via REML.",
  "Cochran Q and I2 use within-study SE only (conventional Q-based I2).",
  "Both conventional REML z CIs and conservative ad hoc Hartung-Knapp",
  "t CIs with k-1 degrees of freedom are included; prefer HK for reporting.",
  "At k=4, tau2 and I2 are imprecise and power for heterogeneity tests is low.",
  "Cohort non-independence/duplicate patients has NOT been ruled out yet.",
  "GSE68465 shows evidence of non-proportional hazards in Script 21.",
  "Results should be described as average prognostic associations with",
  "potential time-variation and clinical adjustment limitations.",
  "Do not interpret combined sample sizes as fully independent without overlap audit.",
  "No C-indices, subgroup associations, or individual participant data are pooled."
)
writeLines(notes,
           file.path(OUT, "LUAD_MPCDS_meta_methods_and_cautions.txt"))

cat("\n================ SCRIPT 42 FINISHED ================\n")
cat("Original Script 21 unadjusted estimates: verified in all four cohorts\n")
cat("Output directory:", OUT, "\n")
cat("No existing Script 21 or Script 41 outputs were modified.\n")
cat("IMPORTANT: stage adjustment differs by cohort; only 4 cohorts;",
    "check patient overlap and PH assumptions before manuscript integration.\n")
cat("====================================================\n")
