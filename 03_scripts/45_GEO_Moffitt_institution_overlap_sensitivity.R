# =============================================================================
# SCRIPT 45 — EXTERNAL GEO INSTITUTION-OVERLAP SENSITIVITY
# Kankarofi et al. | Cell Death Discovery revisions
#
# Rationale:
#   79 HLM (H. Lee Moffitt) samples in GSE68465 share institutional provenance
#   with GSE72094. Distinct GSM IDs do not establish distinct patients.
#   This script EXCLUDES ALL GSE68465 HLM samples as a conservative sensitivity,
#   WITHOUT asserting that any patient is actually duplicated.
#
# Inputs (existing, unmodified):
#   Script 41 clinical_harmonized_patient_variables.csv and Cox estimates
#   Script 42 random-effects meta summary
#   Script 19 GSE68465 selected GEO phenotype annotations
#
# Outputs (new folder only):
#   04_results/mpcds_clinically_adjusted_external/
#      institutional_overlap_sensitivity/
#      GSE68465_site_and_survival_counts.csv
#      GSE68465_HLM_exclusion_cox_comparison.csv
#      GSE68465_HLM_exclusion_PH_tests.csv
#      LUAD_HLM_exclusion_meta_comparison.csv
#      LUAD_HLM_exclusion_meta_cohort_contributions.csv
#      methods_and_limitations.txt
#
# IMPORTANT:
#   Original MPCDS_z values and clinical definitions remain unchanged. The
#   remaining samples are NOT re-standardized, to isolate effect of exclusion.
#   The base Cox fit uses Efron ties, matching Script 41.
#   The REML / ad hoc Hartung-Knapp meta estimator matches Script 42.
#   The results cannot prove non-overlap, and any GSE68465 PH violation remains
#   a methodological limitation of the pooled constant-HR analyses.
# =============================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
BASE <- file.path(ROOT, "04_results", "mpcds_clinically_adjusted_external")
OUT <- file.path(BASE, "institutional_overlap_sensitivity")
COHORTS <- c("GSE68465", "GSE72094", "GSE31210", "GSE50081")
MODELS <- c("unadjusted", "age_sex_adjusted", "stage_adjusted")

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The survival package is required: install.packages('survival')")
}

patient_path <- file.path(BASE, "clinical_harmonized_patient_variables.csv")
cox_path <- file.path(BASE, "LUAD_MPCDS_clinical_Cox_by_cohort.csv")
meta_path <- file.path(BASE, "meta_analysis", "LUAD_MPCDS_REML_meta_summary.csv")
pheno_path <- file.path(ROOT, "02_processed_data", "GEO", "external_validation",
                        "GSE68465", "GSE68465_phenotype_selected.csv")
for (p in c(patient_path, cox_path, meta_path, pheno_path)) {
  if (!file.exists(p)) stop("Required input not found: ", p)
}

patients <- read.csv(patient_path, stringsAsFactors = FALSE, check.names = FALSE)
prior <- read.csv(cox_path, stringsAsFactors = FALSE, check.names = FALSE)
prior_meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
pheno <- read.csv(pheno_path, stringsAsFactors = FALSE, check.names = FALSE)

required_patient <- c("cohort", "sample_id", "OS_time", "OS_status", "MPCDS_z",
                      "age_10yr", "male", "t2plus", "npos")
required_prior <- c("cohort", "model", "status", "n_analysed", "events_analysed",
                    "log_HR", "log_HR_SE", "HR_per_MPCDS_SD")
required_meta <- c("model", "pooled_HR", "HK_ad_hoc_HR_CI_lower",
                   "HK_ad_hoc_HR_CI_upper", "HK_ad_hoc_p", "I2_percent_Q_based")

check_columns <- function(d, cols, what) {
  missing <- setdiff(cols, colnames(d))
  if (length(missing)) stop(what, ": missing columns: ", paste(missing, collapse = ", "))
}
check_columns(patients, required_patient, "Script 41 patient data")
check_columns(prior, required_prior, "Script 41 Cox summary")
check_columns(prior_meta, required_meta, "Script 42 meta summary")
check_columns(pheno, c("geo_accession", "source_name_ch1"), "GSE68465 phenotype")
if (anyDuplicated(paste(patients$cohort, patients$sample_id, sep = "|"))) {
  stop("Duplicate sample identifiers within the Script 41 patient table")
}
if (anyDuplicated(pheno$geo_accession)) {
  stop("Duplicate GSE68465 GEO accessions in phenotype")
}
if (anyDuplicated(paste(prior$cohort, prior$model, sep = "|"))) {
  stop("Duplicate cohort/model records in Script 41 Cox summary")
}

# Site label is assigned strictly by matching GEO accession, never row position.
suball <- patients[patients$cohort == "GSE68465", , drop = FALSE]
match_ix <- match(suball$sample_id, trimws(as.character(pheno$geo_accession)))
if (anyNA(match_ix)) stop("GSE68465 survival samples do not all match phenotype")
suball$site <- toupper(trimws(as.character(pheno$source_name_ch1[match_ix])))
if (anyNA(suball$site) || any(suball$site == "")) stop("Missing institutional code")

all_sites <- toupper(trimws(as.character(pheno$source_name_ch1)))
site_codes <- sort(unique(all_sites))
site_counts <- data.frame(
  site = site_codes,
  n_phenotype_samples = as.integer(table(factor(all_sites, levels = site_codes))),
  n_survival_samples = as.integer(table(factor(suball$site, levels = site_codes))),
  stringsAsFactors = FALSE
)
site_counts$n_excluded_from_survival <- site_counts$n_phenotype_samples -
  site_counts$n_survival_samples
if (!"HLM" %in% site_counts$site) stop("HLM institutional code absent")
n_hlm <- sum(suball$site == "HLM")
if (n_hlm != 79L) {
  warning("Expected 79 HLM survival samples based on previous audit; found ", n_hlm,
          ". Inspect data provenance before interpreting results.")
}
if (n_hlm == 0L || n_hlm == nrow(suball)) stop("No valid HLM exclusion possible")
non_hlm <- suball[suball$site != "HLM", , drop = FALSE]

cat("\n===== GSE68465 INSTITUTIONAL PROVENANCE =====\n")
print(site_counts, row.names = FALSE)
cat("HLM samples excluded:", n_hlm, "\n")
cat("GSE68465 survival patients retained:", nrow(non_hlm), "of", nrow(suball), "\n")

get_covariates <- function(model_name) {
  switch(model_name,
         unadjusted = character(0),
         age_sex_adjusted = c("age_10yr", "male"),
         stage_adjusted = c("age_10yr", "male", "t2plus", "npos"),
         stop("Unknown model specification: ", model_name))
}

fit_cox <- function(d, model_name, label) {
  covars <- get_covariates(model_name)
  vars <- c("MPCDS_z", covars)
  needed <- c("OS_time", "OS_status", vars)
  cc <- complete.cases(d[, needed, drop = FALSE]) &
    is.finite(d$OS_time) & d$OS_time > 0 &
    is.finite(d$MPCDS_z) & d$OS_status %in% c(0L, 1L)
  dat <- d[cc, , drop = FALSE]
  n <- nrow(dat)
  deaths <- sum(dat$OS_status == 1L)
  if (n < 25 || deaths < 15 || deaths / length(vars) < 10) {
    stop(label, "/", model_name, ": insufficient complete-case sample/events")
  }
  for (v in vars) {
    if (length(unique(dat[[v]])) < 2) stop(label, "/", model_name, ": ", v, " has no variation")
  }
  form <- stats::as.formula(paste(
    "survival::Surv(OS_time, OS_status) ~", paste(vars, collapse = " + ")
  ))
  fit <- survival::coxph(form, data = dat, ties = "efron", x = TRUE, y = TRUE,
                         singular.ok = FALSE)
  beta <- unname(stats::coef(fit)["MPCDS_z"])
  se <- sqrt(unname(stats::vcov(fit)["MPCDS_z", "MPCDS_z"]))
  if (!is.finite(beta) || !is.finite(se) || se <= 0) stop("Invalid Cox coefficient or SE")
  tab <- data.frame(
    data_subset = label, model = model_name,
    n_analysed = n, events_analysed = deaths,
    log_HR = beta, log_HR_SE = se, HR_per_MPCDS_SD = exp(beta),
    CI95_lower = exp(beta - 1.96 * se), CI95_upper = exp(beta + 1.96 * se),
    p_value = 2 * stats::pnorm(-abs(beta / se)),
    stringsAsFactors = FALSE
  )
  list(tab = tab, fit = fit)
}

full_fit <- list()
filtered_fit <- list()
comparison <- list()
ph_rows <- list()
for (model_name in MODELS) {
  full <- fit_cox(suball, model_name, "all_GSE68465_sites")
  filtered <- fit_cox(non_hlm, model_name, "excluding_HLM")
  prior_row <- prior[prior$cohort == "GSE68465" & prior$model == model_name &
                       prior$status == "OK", , drop = FALSE]
  if (nrow(prior_row) != 1L) stop("Missing or duplicate Script 41 estimate: ", model_name)
  if (full$tab$n_analysed != prior_row$n_analysed ||
      full$tab$events_analysed != prior_row$events_analysed ||
      abs(full$tab$HR_per_MPCDS_SD - prior_row$HR_per_MPCDS_SD) > 1e-6 ||
      abs(full$tab$log_HR_SE - prior_row$log_HR_SE) > 1e-6) {
    stop("Unfiltered GSE68465 Cox model failed Script 41 reproduction: ", model_name)
  }
  full_fit[[model_name]] <- full$tab
  filtered_fit[[model_name]] <- filtered$tab
  comparison[[model_name]] <- data.frame(
    model = model_name,
    full_n = full$tab$n_analysed, full_events = full$tab$events_analysed,
    full_HR = full$tab$HR_per_MPCDS_SD,
    full_CI_lower = full$tab$CI95_lower, full_CI_upper = full$tab$CI95_upper,
    excluding_HLM_n = filtered$tab$n_analysed,
    excluding_HLM_events = filtered$tab$events_analysed,
    excluding_HLM_HR = filtered$tab$HR_per_MPCDS_SD,
    excluding_HLM_CI_lower = filtered$tab$CI95_lower,
    excluding_HLM_CI_upper = filtered$tab$CI95_upper,
    delta_HR = filtered$tab$HR_per_MPCDS_SD - full$tab$HR_per_MPCDS_SD,
    delta_log_HR = filtered$tab$log_HR - full$tab$log_HR,
    stringsAsFactors = FALSE
  )
  for (which in c("all_GSE68465_sites", "excluding_HLM")) {
    fit <- if (which == "excluding_HLM") filtered$fit else full$fit
    zph <- survival::cox.zph(fit, transform = "km", terms = FALSE, global = TRUE)
    ztable <- as.data.frame(zph$table)
    ph_rows[[paste(model_name, which, sep = "__")]] <- data.frame(
      data_subset = which, model = model_name, term = rownames(ztable),
      chisquare = ztable$chisq, df = ztable$df, p_value = ztable$p,
      stringsAsFactors = FALSE
    )
  }
  cat(sprintf("%-17s | all HR %.4f (n=%d, events=%d) | excluding HLM HR %.4f (n=%d, events=%d)\n",
              model_name, full$tab$HR_per_MPCDS_SD,
              full$tab$n_analysed, full$tab$events_analysed,
              filtered$tab$HR_per_MPCDS_SD,
              filtered$tab$n_analysed, filtered$tab$events_analysed))
}
comparison_df <- do.call(rbind, comparison)
ph_df <- do.call(rbind, ph_rows)
rownames(comparison_df) <- NULL
rownames(ph_df) <- NULL

# Same restricted-likelihood random-effects estimator and ad hoc HK CI as Script 42.
meta_reml_42 <- function(x) {
  y <- as.numeric(x$log_HR)
  vi <- as.numeric(x$log_HR_SE)^2
  k <- length(y)
  if (k < 3L || any(!is.finite(y)) || any(!is.finite(vi)) || any(vi <= 0)) {
    stop("Meta-analysis requires >= 3 valid log-HRs and positive variances")
  }
  nll <- function(tau2) {
    v <- vi + tau2
    w <- 1 / v
    mu <- sum(w * y) / sum(w)
    0.5 * (sum(log(v)) + log(sum(w)) + sum(w * (y - mu)^2))
  }
  upper <- max(1, 10 * diff(range(y))^2, 50 * max(vi))
  sol <- stats::optimize(nll, interval = c(0, upper), tol = 1e-11)
  tau2 <- if (nll(0) <= sol$objective) 0 else sol$minimum
  w <- 1 / (vi + tau2)
  mu <- sum(w * y) / sum(w)
  se_mu <- sqrt(1 / sum(w))
  q_re <- sum(w * (y - mu)^2)
  hk_scale <- max(1, q_re / (k - 1))
  se_hk <- se_mu * sqrt(hk_scale)
  tcrit <- stats::qt(0.975, df = k - 1)
  wf <- 1 / vi
  mu_f <- sum(wf * y) / sum(wf)
  Q <- sum(wf * (y - mu_f)^2)
  I2 <- if (Q > 0) 100 * max(0, (Q - (k - 1)) / Q) else 0
  result <- data.frame(
    cohorts = k, n_analysed_sum = sum(x$n_analysed),
    events_sum = sum(x$events_analysed),
    pooled_log_HR = mu, pooled_HR = exp(mu),
    conventional_CI_lower = exp(mu - stats::qnorm(.975) * se_mu),
    conventional_CI_upper = exp(mu + stats::qnorm(.975) * se_mu),
    HK_CI_lower = exp(mu - tcrit * se_hk),
    HK_CI_upper = exp(mu + tcrit * se_hk),
    HK_p = 2 * stats::pt(-abs(mu / se_hk), df = k - 1),
    I2_percent = I2, tau2_REML = tau2,
    Cochran_Q_p = stats::pchisq(Q, df = k - 1, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  list(summary = result, weights = w)
}

meta_comparison <- list()
meta_contribution <- list()
for (model_name in MODELS) {
  original_rows <- prior[prior$model == model_name & prior$cohort %in% COHORTS &
                           prior$status == "OK", , drop = FALSE]
  original_rows <- original_rows[match(COHORTS, original_rows$cohort), , drop = FALSE]
  if (nrow(original_rows) != 4 || anyNA(original_rows$cohort)) {
    stop("Missing required Script 41 cohort estimates for model: ", model_name)
  }
  original <- meta_reml_42(original_rows)
  old_meta <- prior_meta[prior_meta$model == model_name, , drop = FALSE]
  if (nrow(old_meta) != 1L ||
      abs(original$summary$pooled_HR - old_meta$pooled_HR) > 1e-5 ||
      abs(original$summary$HK_CI_lower - old_meta$HK_ad_hoc_HR_CI_lower) > 1e-5 ||
      abs(original$summary$HK_CI_upper - old_meta$HK_ad_hoc_HR_CI_upper) > 1e-5 ||
      abs(original$summary$I2_percent - old_meta$I2_percent_Q_based) > 1e-4) {
    stop("Meta-analysis reproduction does not match Script 42 for: ", model_name)
  }
  sensitivity_rows <- original_rows
  ix <- match("GSE68465", sensitivity_rows$cohort)
  new <- filtered_fit[[model_name]]
  sensitivity_rows$log_HR[ix] <- new$log_HR
  sensitivity_rows$log_HR_SE[ix] <- new$log_HR_SE
  sensitivity_rows$HR_per_MPCDS_SD[ix] <- new$HR_per_MPCDS_SD
  sensitivity_rows$n_analysed[ix] <- new$n_analysed
  sensitivity_rows$events_analysed[ix] <- new$events_analysed
  sensitivity <- meta_reml_42(sensitivity_rows)
  meta_comparison[[model_name]] <- data.frame(
    model = model_name,
    original_pooled_HR = original$summary$pooled_HR,
    original_HK_CI_lower = original$summary$HK_CI_lower,
    original_HK_CI_upper = original$summary$HK_CI_upper,
    original_HK_p = original$summary$HK_p,
    original_I2_percent = original$summary$I2_percent,
    excluding_HLM_pooled_HR = sensitivity$summary$pooled_HR,
    excluding_HLM_HK_CI_lower = sensitivity$summary$HK_CI_lower,
    excluding_HLM_HK_CI_upper = sensitivity$summary$HK_CI_upper,
    excluding_HLM_HK_p = sensitivity$summary$HK_p,
    excluding_HLM_I2_percent = sensitivity$summary$I2_percent,
    excluding_HLM_n_sum = sensitivity$summary$n_analysed_sum,
    excluding_HLM_events_sum = sensitivity$summary$events_sum,
    change_pooled_HR = sensitivity$summary$pooled_HR - original$summary$pooled_HR,
    stringsAsFactors = FALSE
  )
  sensitivity_rows$revised_weight_percent <- 100 * sensitivity$weights / sum(sensitivity$weights)
  sensitivity_rows$model <- model_name
  meta_contribution[[model_name]] <- sensitivity_rows[, c(
    "model", "cohort", "n_analysed", "events_analysed", "log_HR", "log_HR_SE",
    "HR_per_MPCDS_SD", "revised_weight_percent")]
}
meta_df <- do.call(rbind, meta_comparison)
contrib_df <- do.call(rbind, meta_contribution)
rownames(meta_df) <- NULL
rownames(contrib_df) <- NULL

# Save everything only after all validation and analyses have succeeded.
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(site_counts, file.path(OUT, "GSE68465_site_and_survival_counts.csv"), row.names = FALSE)
write.csv(comparison_df, file.path(OUT, "GSE68465_HLM_exclusion_cox_comparison.csv"), row.names = FALSE)
write.csv(ph_df, file.path(OUT, "GSE68465_HLM_exclusion_PH_tests.csv"), row.names = FALSE)
write.csv(meta_df, file.path(OUT, "LUAD_HLM_exclusion_meta_comparison.csv"), row.names = FALSE)
write.csv(contrib_df, file.path(OUT, "LUAD_HLM_exclusion_meta_cohort_contributions.csv"), row.names = FALSE)
writeLines(c(
  "Script 45 — Institutional-provenance exclusion sensitivity (GSE68465 / GSE72094).",
  "HLM denotes the H. Lee Moffitt institution subset in GSE68465. GSE72094 is Moffitt-derived.",
  "Distinct GEO accessions and distinct-looking sample labels do not establish distinct patients.",
  "Remove all HLM samples from GSE68465, regardless of suspected identity; no duplicates are asserted.",
  "Cox baseline covariates and model specifications reproduce Script 41 (Efron ties).",
  "MPCDS_z is kept at its original full-cohort scale; no re-standardisation after filtering.",
  "Rerun only GSE68465 Cox models on non-HLM subjects; other 3 cohorts unchanged.",
  "Meta-analysis uses Script 42 REML and conservative ad hoc Hartung-Knapp estimator.",
  "PH is re-tested via cox.zph for full and non-HLM GSE68465 samples.",
  "Only exploratory sensitivity. Institutional exclusion may remove genuine independent patients.",
  "Site exclusion does not prove no overlap; direct patient-level identity cannot be verified.",
  "Any Cox PH violations, different stage encodings and n=4 meta cohorts remain limitations.",
  "Existing Script 41/42 results are not overwritten."
), file.path(OUT, "methods_and_limitations.txt"))

cat("\n===== GSE68465 PROPORTIONAL-HAZARDS CHECK =====\n")
print(ph_df[ph_df$model == "stage_adjusted" & ph_df$term %in% c("MPCDS_z", "GLOBAL"),
            c("data_subset", "term", "p_value")], row.names = FALSE, digits = 5)
cat("\n===== POOLED HR SENSITIVITY: EXCLUDE HLM FROM GSE68465 =====\n")
print(meta_df[, c("model", "original_pooled_HR", "excluding_HLM_pooled_HR",
                  "excluding_HLM_HK_CI_lower", "excluding_HLM_HK_CI_upper",
                  "excluding_HLM_HK_p", "excluding_HLM_I2_percent")],
      row.names = FALSE, digits = 5)
cat("\n===== SCRIPT 45 FINISHED =====\n")
cat("All three original GSE68465 models reproduced Script 41.\n")
cat("All three original pooled estimates reproduced Script 42.\n")
cat("Folder:", OUT, "\n")
cat("The exclusion sensitivity does NOT establish patient identity or no overlap.\n")
