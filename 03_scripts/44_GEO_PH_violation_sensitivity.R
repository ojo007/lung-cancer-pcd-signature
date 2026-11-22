# ==============================================================================
# SCRIPT 44 — EXTERNAL LUAD PROPORTIONAL-HAZARDS SENSITIVITY ANALYSIS (v2: survSplit fix)
# Kankarofi et al. — Cell Death Discovery manuscript corrections
#
# Reads the verified Script 41 patient-level data and Script 43 PH diagnostics.
# 1. GSE68465: estimate interval-specific MPCDS effects (0–12, 12–36, >36
#    months) in a time-split, stage-adjusted Cox model. A likelihood-ratio test
#    compares interval-specific MPCDS effects with one constant effect.
# 2. GSE31210: sensitivity with sex/stage STRATIFIED baseline hazards, given
#    the observed PH violations in those two variables. Only 35 deaths, so
#    this sensitivity is exploratory and potentially imprecise.
#
# Original Script 41 and Script 42 outputs are never modified. DO NOT directly
# pool these time-specific or stratified estimates with the original pooled HR.
# No coefficients are refitted for the MPCDS gene signature itself.
# ==============================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
BASE <- file.path(ROOT, "04_results", "mpcds_clinically_adjusted_external")
INPUT <- file.path(BASE, "clinical_harmonized_patient_variables.csv")
SUMMARY <- file.path(BASE, "LUAD_MPCDS_clinical_Cox_by_cohort.csv")
PH_INPUT <- file.path(BASE, "proportional_hazards_diagnostics",
                      "Cox_Schoenfeld_PH_tests.csv")
OUT <- file.path(BASE, "PH_violation_sensitivity")

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The survival package is required: install.packages('survival')")
}
# survSplit inspects the formula syntax and specifically requires a bare
# Surv(...) call (rather than survival::Surv(...)). Bind it without attaching
# the package. Also expose strata() for the Cox formula evaluation environment.
Surv <- survival::Surv
strata <- survival::strata
for (f in c(INPUT, SUMMARY, PH_INPUT)) {
  if (!file.exists(f)) stop("Required verified input not found: ", f)
}
patients <- read.csv(INPUT, stringsAsFactors = FALSE, check.names = FALSE)
previous <- read.csv(SUMMARY, stringsAsFactors = FALSE, check.names = FALSE)
ph <- read.csv(PH_INPUT, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("cohort", "sample_id", "OS_time", "OS_status", "MPCDS_z",
              "age_10yr", "male", "t2plus", "npos", "advanced_stage")
if (!all(required %in% names(patients))) {
  stop("Script 41 patient table missing: ",
       paste(setdiff(required, names(patients)), collapse = ", "))
}
if (anyDuplicated(paste(patients$cohort, patients$sample_id))) {
  stop("Duplicate patient/cohort identifiers")
}

make_data <- function(cohort, vars) {
  x <- patients[patients$cohort == cohort, , drop = FALSE]
  fields <- c("OS_time", "OS_status", vars)
  cc <- stats::complete.cases(x[, fields, drop = FALSE]) &
    is.finite(x$OS_time) & x$OS_time > 0 &
    x$OS_status %in% c(0L, 1L)
  dat <- x[cc, c("sample_id", fields), drop = FALSE]
  if (anyDuplicated(dat$sample_id)) stop(cohort, ": duplicate sample IDs")
  record <- previous[previous$cohort == cohort &
                       previous$model == "stage_adjusted" &
                       previous$status == "OK", , drop = FALSE]
  if (nrow(record) != 1L) stop(cohort, ": missing unique successful Script 41 stage model")
  if (nrow(dat) != record$n_analysed ||
      sum(dat$OS_status == 1L) != record$events_analysed) {
    stop(cohort, ": analysis population differs from Script 41")
  }
  list(dat = dat, reference = record)
}

extract_mpcds <- function(fit, term = "MPCDS_z") {
  b <- stats::coef(fit)
  v <- stats::vcov(fit)
  if (!term %in% names(b) || !is.finite(b[[term]]) ||
      !is.finite(v[term, term]) || v[term, term] <= 0) {
    stop("Invalid Cox coefficient or variance for ", term)
  }
  se <- sqrt(v[term, term])
  z <- b[[term]] / se
  data.frame(term = term, log_HR = unname(b[[term]]), SE = unname(se),
             HR = exp(b[[term]]),
             CI95_lower = exp(b[[term]] - 1.96 * se),
             CI95_upper = exp(b[[term]] + 1.96 * se),
             p_value = 2 * stats::pnorm(-abs(z)),
             stringsAsFactors = FALSE)
}

# -------------------------------- GSE68465 ------------------------------------
cat("\n============= GSE68465: MPCDS TIME-VARYING SENSITIVITY =============\n")
vars_a <- c("MPCDS_z", "age_10yr", "male", "t2plus", "npos")
res_a <- make_data("GSE68465", vars_a)
a <- res_a$dat
original_a <- survival::coxph(
  survival::Surv(OS_time, OS_status) ~ MPCDS_z + age_10yr + male + t2plus + npos,
  data = a, ties = "efron", singular.ok = FALSE)
original_hr_a <- exp(unname(stats::coef(original_a)["MPCDS_z"]))
if (abs(original_hr_a - res_a$reference$HR_per_MPCDS_SD) > 1e-6) {
  stop("GSE68465 original stage-adjusted HR does not reproduce Script 41")
}

# These intervals are descriptive, pre-specified for this diagnostic and are
# NOT chosen by optimising model performance or statistical significance.
cut_months <- c(12, 36)
s <- survival::survSplit(
  Surv(OS_time, OS_status) ~ .,
  data = a, cut = cut_months, start = "tstart", end = "tstop",
  event = "event", id = "source_row_id"
)
period_levels <- c("0-12 months", "12-36 months", ">36 months")
s$period <- factor(
  ifelse(s$tstart < cut_months[1], period_levels[1],
         ifelse(s$tstart < cut_months[2], period_levels[2], period_levels[3])),
  levels = period_levels
)
# Check splitting preserves each patient's final outcome exactly once.
if (nrow(s) < nrow(a) || length(unique(s$source_row_id)) != nrow(a) ||
    sum(s$event) != sum(a$OS_status)) {
  stop("GSE68465 survival splitting did not preserve patients/events")
}
interval_counts <- do.call(rbind, lapply(period_levels, function(z) {
  sub <- s[s$period == z, , drop = FALSE]
  data.frame(cohort = "GSE68465", period = z,
             patients_at_risk_during_interval = length(unique(sub$source_row_id)),
             events = sum(sub$event), stringsAsFactors = FALSE)
}))
if (any(interval_counts$events < 5L)) {
  warning("A GSE68465 interval has fewer than five deaths; estimates may be unstable")
}
print(interval_counts, row.names = FALSE)

# Stratifying baseline hazard by time period is equivalent to allowing
# separate unspecified baseline hazards in each of the intervals.
constant <- survival::coxph(
  survival::Surv(tstart, tstop, event) ~
    MPCDS_z + age_10yr + male + t2plus + npos + strata(period),
  data = s, ties = "efron", singular.ok = FALSE)

s$risk_0_12 <- ifelse(s$period == period_levels[1], s$MPCDS_z, 0)
s$risk_12_36 <- ifelse(s$period == period_levels[2], s$MPCDS_z, 0)
s$risk_gt36 <- ifelse(s$period == period_levels[3], s$MPCDS_z, 0)
variable <- survival::coxph(
  survival::Surv(tstart, tstop, event) ~
    risk_0_12 + risk_12_36 + risk_gt36 + age_10yr + male +
    t2plus + npos + strata(period),
  data = s, ties = "efron", singular.ok = FALSE)

term_map <- c("risk_0_12", "risk_12_36", "risk_gt36")
interval_hr <- do.call(rbind, lapply(seq_along(term_map), function(i) {
  x <- extract_mpcds(variable, term_map[i])
  x$cohort <- "GSE68465"
  x$period <- period_levels[i]
  x$n <- nrow(a)
  x$events_total <- sum(a$OS_status)
  x$events_in_period <- interval_counts$events[i]
  x
}))
interval_hr <- interval_hr[, c("cohort", "period", "n", "events_total",
                               "events_in_period", "term", "log_HR", "SE",
                               "HR", "CI95_lower", "CI95_upper", "p_value")]

# The two models use the same split data, risk sets, baseline strata and
# clinical covariates. Only the MPCDS coefficient is allowed to vary.
LR <- max(0, 2 * (as.numeric(variable$loglik[2]) -
                   as.numeric(constant$loglik[2])))
extra_df <- length(stats::coef(variable)) - length(stats::coef(constant))
if (extra_df != 2L) stop("Unexpected GSE68465 time-varying model degrees of freedom")
p_LR <- stats::pchisq(LR, df = extra_df, lower.tail = FALSE)
constant_hr <- extract_mpcds(constant)
if (abs(constant_hr$HR - original_hr_a) > 1e-5) {
  warning("Constant-effect split model differs from original: inspect risk sets")
}
comparison <- data.frame(
  cohort = "GSE68465", n = nrow(a), events = sum(a$OS_status),
  original_stage_adjusted_HR = original_hr_a,
  constant_effect_split_HR = constant_hr$HR,
  likelihood_ratio_chisquare = LR, additional_df = extra_df,
  p_for_MPCDS_effect_varying_by_period = p_LR,
  period_cut_months = "12; 36", stringsAsFactors = FALSE
)
cat("\nGSE68465 interval-specific MPCDS HRs:\n")
print(interval_hr[, c("period", "events_in_period", "HR", "CI95_lower",
                       "CI95_upper", "p_value")], row.names = FALSE, digits = 4)
cat(sprintf("GSE68465 constant vs interval-varying MPCDS effect: LR chi2=%.4f, df=%d, p=%.4g\n",
            LR, extra_df, p_LR))

# ------------------------------- GSE31210 -------------------------------------
cat("\n============= GSE31210: STRATIFIED STAGE/SEX SENSITIVITY =============\n")
vars_b <- c("MPCDS_z", "age_10yr", "male", "advanced_stage")
res_b <- make_data("GSE31210", vars_b)
b <- res_b$dat
original_b <- survival::coxph(
  survival::Surv(OS_time, OS_status) ~
    MPCDS_z + age_10yr + male + advanced_stage,
  data = b, ties = "efron", singular.ok = FALSE)
original_hr_b <- exp(unname(stats::coef(original_b)["MPCDS_z"]))
if (abs(original_hr_b - res_b$reference$HR_per_MPCDS_SD) > 1e-6) {
  stop("GSE31210 original stage-adjusted HR does not reproduce Script 41")
}

# Allow each observed male-by-stage group a distinct baseline hazard. This
# avoids requiring proportional hazards for sex or the grouped stage term.
b$sex_stage <- interaction(b$male, b$advanced_stage, drop = TRUE, sep = "_")
stratum_counts <- do.call(rbind, lapply(levels(b$sex_stage), function(z) {
  x <- b[b$sex_stage == z, , drop = FALSE]
  data.frame(cohort = "GSE31210", sex_stage_stratum = as.character(z),
             n = nrow(x), events = sum(x$OS_status), stringsAsFactors = FALSE)
}))
print(stratum_counts, row.names = FALSE)
if (any(stratum_counts$events < 5L)) {
  warning("Sparse GSE31210 stratum (<5 deaths): interpretation is exploratory")
}
if (sum(stratum_counts$events > 0L) < 2L) {
  warning("Insufficient event-bearing strata for a useful stratified comparison")
  strat_row <- data.frame(cohort = "GSE31210", n = nrow(b),
                          events = sum(b$OS_status),
                          original_stage_adjusted_HR = original_hr_b,
                          stratified_HR = NA_real_, CI95_lower = NA_real_,
                          CI95_upper = NA_real_, p_value = NA_real_,
                          MPCDS_PH_p = NA_real_, status = "SKIPPED_FEW_EVENT_STRATA")
} else {
  strat_fit <- tryCatch(
    survival::coxph(
      survival::Surv(OS_time, OS_status) ~
        MPCDS_z + age_10yr + strata(sex_stage),
      data = b, ties = "efron", singular.ok = FALSE),
    error = function(e) e
  )
  if (inherits(strat_fit, "error")) {
    warning("GSE31210 stratified Cox failed: ", conditionMessage(strat_fit))
    strat_row <- data.frame(cohort = "GSE31210", n = nrow(b),
                            events = sum(b$OS_status),
                            original_stage_adjusted_HR = original_hr_b,
                            stratified_HR = NA_real_, CI95_lower = NA_real_,
                            CI95_upper = NA_real_, p_value = NA_real_,
                            MPCDS_PH_p = NA_real_,
                            status = paste0("FIT_FAILED: ", conditionMessage(strat_fit)))
  } else {
    val <- extract_mpcds(strat_fit)
    ph_res <- tryCatch(survival::cox.zph(strat_fit, transform = "km",
                                         terms = FALSE, global = TRUE),
                       error = function(e) NULL)
    ph_p <- if (!is.null(ph_res) && "MPCDS_z" %in% rownames(ph_res$table))
      ph_res$table["MPCDS_z", "p"] else NA_real_
    strat_row <- data.frame(cohort = "GSE31210", n = nrow(b),
                            events = sum(b$OS_status),
                            original_stage_adjusted_HR = original_hr_b,
                            stratified_HR = val$HR, CI95_lower = val$CI95_lower,
                            CI95_upper = val$CI95_upper, p_value = val$p_value,
                            MPCDS_PH_p = ph_p, status = "OK")
  }
}
cat("\nGSE31210 stratified sensitivity:\n")
print(strat_row, row.names = FALSE, digits = 4)

# --------------------------------- Outputs ------------------------------------
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(interval_counts,
          file.path(OUT, "GSE68465_time_interval_counts.csv"), row.names = FALSE)
write.csv(interval_hr,
          file.path(OUT, "GSE68465_time_interval_MPCDS_HRs.csv"), row.names = FALSE)
write.csv(comparison,
          file.path(OUT, "GSE68465_constant_vs_varying_test.csv"), row.names = FALSE)
write.csv(stratum_counts,
          file.path(OUT, "GSE31210_stage_sex_stratum_counts.csv"), row.names = FALSE)
write.csv(strat_row,
          file.path(OUT, "GSE31210_stratified_MPCDS_sensitivity.csv"), row.names = FALSE)
writeLines(c(
  "Time cutoffs were specified as 12 and 36 months, not outcome-optimised.",
  "GSE68465: interval-specific hazard ratios are descriptive and should NOT",
  "be substituted for the common-HR estimates in Script 42 meta-analysis.",
  "GSE31210: 35 events; stage/sex-stratified sensitivity can be imprecise.",
  "The global PH test for GSE31210 implicated both sex and stage, not MPCDS.",
  "No new treatment-effect, model-superiority or causal claims are warranted.",
  "Patient overlap across GEO studies is still not ruled out.",
  "This script leaves Scripts 21, 41, 42 and 43 outputs unchanged."
), file.path(OUT, "methods_and_interpretation_cautions.txt"))
cat("\n============= SCRIPT 44 FINISHED =============\n")
cat("Output folder: ", OUT, "\n", sep = "")
cat("Sensitivity estimates are separate from the original random-effects meta-analysis.\n")
