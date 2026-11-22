# ==============================================================================
# SCRIPT 41 — CLINICALLY ADJUSTED EXTERNAL LUAD MPCDS COX REGRESSION
# Kankarofi et al. — Cell Death Discovery scientific corrections
#
# INPUTS: Unmodified results of Script 21 (MPCDS patient scores) and
#         Script 19 (GEO selected phenotype tables).
#
# Design:
#   - Reuse the already validated MPCDS_z and overall-survival outcomes.
#   - One patient/sample per existing MPCDS survival record; strict ID matching.
#   - Per-cohort Cox models, NO pooling in this script and NO coefficient refit.
#   - Age/sex model; stage-adjusted model; optional smoking sensitivity.
#   - GSE68465 has pT/pN, NOT a conventional overall stage: use binary pT2+
#     and pN1+ proxies rather than manufacture an AJCC stage.
#   - Other cohorts: harmonize pathological stage I vs stage II–IV to limit
#     degrees of freedom in smaller event-count cohorts.
#   - Smoking analysis skipped where <10 deaths per modeled parameter.
#   - Complete-case analysis; no invented clinical values or imputation.
#   - Save only NEW outputs in 04_results/mpcds_clinically_adjusted_external.
#
# Scientific cautions:
#   Different baseline stages and missingness limit interpretation of the
#   stage-adjusted pooled HR. GSE31210 has only 35 events; report explicitly.
#   These are observational prognostic associations, not treatment effects.
#   Existing within-cohort z-scoring is not a frozen-platform deployment test.
# ==============================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
COHORTS <- c("GSE68465", "GSE72094", "GSE31210", "GSE50081")
OUT <- file.path(ROOT, "04_results", "mpcds_clinically_adjusted_external")

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The survival package is required. Install it before running Script 41.")
}
if (!dir.exists(ROOT)) stop("Project root does not exist: ", ROOT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

score_file <- function(cohort) {
  if (cohort == "GSE68465") {
    return(file.path(ROOT, "04_results", "mpcds", "external_validation",
                     "GSE68465_MPCDS_patient_scores.csv"))
  }
  file.path(ROOT, "04_results", "mpcds", "external_validation", cohort,
            paste0(cohort, "_MPCDS_patient_scores.csv"))
}
phenotype_file <- function(cohort) {
  file.path(ROOT, "02_processed_data", "GEO", "external_validation", cohort,
            paste0(cohort, "_phenotype_selected.csv"))
}

clean_text <- function(x) {
  s <- trimws(as.character(x))
  s <- sub("^[^:]+:\\s*", "", s, perl = TRUE) # e.g. 'gender: M'
  missing <- is.na(x) | toupper(s) %in%
    c("", "NA", "N/A", "NULL", "UNKNOWN", "MISSING", "--", "---", "NONE")
  s[missing] <- NA_character_
  s
}
parse_numeric <- function(x) {
  s <- clean_text(x)
  s <- gsub(",", "", s, fixed = TRUE)
  suppressWarnings(as.numeric(s))
}
parse_male <- function(x) {
  s <- toupper(clean_text(x))
  ans <- rep(NA_integer_, length(s))
  ans[s %in% c("MALE", "M")] <- 1L
  ans[s %in% c("FEMALE", "F")] <- 0L
  ans
}
parse_smoking_ever <- function(x) {
  s <- toupper(clean_text(x))
  ans <- rep(NA_integer_, length(s))
  never <- !is.na(s) & grepl("NEVER|NON.?SMOK", s, perl = TRUE)
  ans[never] <- 0L
  ever <- !is.na(s) & grepl("EVER|CURRENT|PAST|FORMER|EX.?SMOK", s, perl = TRUE)
  ans[ever & !never] <- 1L
  ans
}
parse_stage_advanced <- function(x) {
  s <- toupper(clean_text(x))
  s <- gsub("\\s+", "", s, perl = TRUE)
  s <- sub("^STAGE", "", s)
  ans <- rep(NA_integer_, length(s))
  stage_one <- !is.na(s) &
    grepl("^(1[A-D]?|I[A-D]?)$", s, perl = TRUE)
  stage_later <- !is.na(s) &
    grepl("^(2[A-D]?|3[A-D]?|4[A-D]?|II[A-D]?|III[A-D]?|IV[A-D]?)$",
          s, perl = TRUE)
  ans[stage_one] <- 0L
  ans[stage_later] <- 1L
  ans
}
parse_tn <- function(x, prefix) {
  s <- toupper(clean_text(x))
  # Available GSE68465 values are combined pN/pT strings, e.g. pN0pT2.
  pat <- paste0("P", prefix, "[0-4]")
  hit <- regexpr(pat, s, perl = TRUE)
  ans <- rep(NA_integer_, length(s))
  valid <- !is.na(s) & hit > 0L
  if (any(valid)) {
    matched <- regmatches(s[valid], regexpr(pat, s[valid], perl = TRUE))
    ans[valid] <- as.integer(substr(matched, 3L, 3L))
  }
  ans
}

clinical_fields <- list(
  GSE68465 = list(age = "age:ch1", sex = "Sex:ch1",
                  stage = "disease_stage:ch1", smoke = "smoking_history:ch1"),
  GSE72094 = list(age = "age_at_diagnosis:ch1", sex = "gender:ch1",
                  stage = "Stage:ch1", smoke = "smoking_status:ch1"),
  GSE31210 = list(age = "age (years):ch1", sex = "gender:ch1",
                  stage = "pathological stage:ch1", smoke = "smoking status:ch1"),
  GSE50081 = list(age = "age:ch1", sex = "Sex:ch1",
                  stage = "Stage:ch1", smoke = "smoking:ch1")
)

all_models <- list()
all_qc <- list()
all_terms <- list()
patient_rows <- list()

fit_model <- function(dat, cohort, model_name, covariates,
                      allow_sparse = FALSE) {
  required <- c("OS_time", "OS_status", "MPCDS_z", covariates)
  cc <- complete.cases(dat[, required, drop = FALSE]) &
    is.finite(dat$OS_time) & dat$OS_time > 0 &
    is.finite(dat$MPCDS_z) & dat$OS_status %in% c(0, 1)
  sub <- dat[cc, , drop = FALSE]
  n <- nrow(sub)
  events <- sum(sub$OS_status == 1L)
  n_params <- 1L + length(covariates) # MPCDS + binary/continuous covariates
  status <- "OK"
  note <- ""
  model <- NULL

  if (n < 25L || events < 15L) {
    status <- "SKIPPED_TOO_FEW_PATIENTS_OR_EVENTS"
    note <- paste("n=", n, "; events=", events, sep = "")
  } else if (!allow_sparse && (events / n_params < 10)) {
    status <- "SKIPPED_LOW_EVENTS_PER_PARAMETER"
    note <- sprintf("%.2f events per parameter (<10)", events / n_params)
  } else if (any(vapply(sub[, c("MPCDS_z", covariates), drop = FALSE],
                         function(x) length(unique(x)) < 2L, logical(1)))) {
    status <- "SKIPPED_NO_VARIATION"
    note <- "At least one predictor has no variation among complete cases"
  } else {
    f <- stats::as.formula(paste(
      "survival::Surv(OS_time, OS_status) ~",
      paste(c("MPCDS_z", covariates), collapse = " + ")
    ))
    warnings <- character(0)
    model <- tryCatch(
      withCallingHandlers(
        survival::coxph(f, data = sub, ties = "efron", singular.ok = FALSE,
                        model = FALSE, x = FALSE, y = FALSE),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        note <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(model)) status <- "FAILED_COX_FIT"
    if (length(warnings)) note <- paste(c(note, warnings), collapse = " | ")
  }

  log_hr <- se <- hr <- lo <- hi <- p <- NA_real_
  if (!is.null(model)) {
    beta <- stats::coef(model)
    vv <- stats::vcov(model)
    if (!"MPCDS_z" %in% names(beta) ||
        !is.finite(beta[["MPCDS_z"]]) ||
        !is.finite(vv["MPCDS_z", "MPCDS_z"]) ||
        vv["MPCDS_z", "MPCDS_z"] <= 0) {
      status <- "FAILED_NONFINITE_MPCDS_ESTIMATE"
      note <- paste(note, "Nonfinite MPCDS coefficient/variance")
    } else {
      log_hr <- unname(beta[["MPCDS_z"]])
      se <- sqrt(unname(vv["MPCDS_z", "MPCDS_z"]))
      hr <- exp(log_hr)
      lo <- exp(log_hr - 1.96 * se)
      hi <- exp(log_hr + 1.96 * se)
      p <- 2 * stats::pnorm(-abs(log_hr / se))
    }
  }
  one <- data.frame(
    cohort = cohort, model = model_name, status = status,
    n_total = nrow(dat), n_analysed = n, events_analysed = events,
    events_per_parameter = if (n_params) events / n_params else NA_real_,
    variables = paste(c("MPCDS_z", covariates), collapse = " + "),
    log_HR = log_hr, log_HR_SE = se,
    HR_per_MPCDS_SD = hr, CI95_lower = lo, CI95_upper = hi,
    p_value = p, warning_or_note = note, stringsAsFactors = FALSE
  )
  list(summary = one, fit = model, sample_ids = sub$sample_id)
}

for (cohort in COHORTS) {
  cat("\n=============== ", cohort, " ===============\n", sep = "")
  sf <- score_file(cohort)
  pf <- phenotype_file(cohort)
  if (!file.exists(sf)) stop("Missing MPCDS score file: ", sf)
  if (!file.exists(pf)) stop("Missing selected phenotype file: ", pf)
  scores <- read.csv(sf, stringsAsFactors = FALSE, check.names = FALSE)
  pheno <- read.csv(pf, stringsAsFactors = FALSE, check.names = FALSE)

  score_time <- if (cohort == "GSE68465") "OS_months" else "OS_time"
  sc_required <- c("sample_id", score_time, "OS_status", "MPCDS_z")
  missing_sc <- setdiff(sc_required, names(scores))
  if (length(missing_sc)) stop(cohort, ": score columns absent: ",
                               paste(missing_sc, collapse = ", "))
  if (!"geo_accession" %in% names(pheno)) {
    stop(cohort, ": phenotype lacks geo_accession")
  }
  f <- clinical_fields[[cohort]]
  missing_f <- setdiff(unlist(f, use.names = FALSE), names(pheno))
  if (length(missing_f)) stop(cohort, ": clinical columns absent: ",
                              paste(missing_f, collapse = ", "))
  if (anyNA(scores$sample_id) || anyDuplicated(scores$sample_id) ||
      anyNA(pheno$geo_accession) || anyDuplicated(pheno$geo_accession)) {
    stop(cohort, ": duplicated or missing sample identifiers")
  }
  ix <- match(trimws(scores$sample_id), trimws(pheno$geo_accession))
  if (anyNA(ix)) stop(cohort, ": some survival-score patients lack phenotype matches")

  dat <- data.frame(
    cohort = cohort,
    sample_id = as.character(scores$sample_id),
    OS_time = suppressWarnings(as.numeric(scores[[score_time]])),
    OS_status = suppressWarnings(as.integer(as.character(scores$OS_status))),
    MPCDS_z = suppressWarnings(as.numeric(scores$MPCDS_z)),
    age_10yr = parse_numeric(pheno[[f$age]][ix]) / 10,
    male = parse_male(pheno[[f$sex]][ix]),
    ever_smoked = parse_smoking_ever(pheno[[f$smoke]][ix]),
    original_stage = clean_text(pheno[[f$stage]][ix]),
    stringsAsFactors = FALSE
  )
  if (cohort == "GSE68465") {
    tnum <- parse_tn(dat$original_stage, "T")
    nnum <- parse_tn(dat$original_stage, "N")
    dat$t2plus <- ifelse(is.na(tnum), NA_integer_, as.integer(tnum >= 2L))
    dat$npos <- ifelse(is.na(nnum), NA_integer_, as.integer(nnum >= 1L))
    dat$advanced_stage <- NA_integer_
    stage_terms <- c("t2plus", "npos")
  } else {
    dat$advanced_stage <- parse_stage_advanced(dat$original_stage)
    dat$t2plus <- NA_integer_
    dat$npos <- NA_integer_
    stage_terms <- "advanced_stage"
  }
  if (any(!is.na(dat$OS_status) & !(dat$OS_status %in% c(0L, 1L)))) {
    stop(cohort, ": OS_status must be coded 0/1")
  }
  if (any(!is.na(dat$age_10yr) &
          (dat$age_10yr < 1.8 | dat$age_10yr > 12))) {
    stop(cohort, ": age parsing produced values outside 18–120 years")
  }
  if (any(!is.finite(dat$OS_time) | dat$OS_time <= 0, na.rm = TRUE)) {
    warning(cohort, ": nonpositive/invalid survival times will be excluded")
  }
  stage_ok <- complete.cases(dat[, stage_terms, drop = FALSE])
  if (mean(!is.na(dat$age_10yr)) < 0.90 || mean(!is.na(dat$male)) < 0.90 ||
      mean(stage_ok) < 0.90) {
    stop(cohort, ": under 90% age, sex or staging completeness after parsing. ",
         "Inspect the source fields before fitting Cox models.")
  }
  if (!nrow(dat) || all(is.na(dat$MPCDS_z))) stop(cohort, ": no MPCDS values")

  cat("Score patients:", nrow(dat), "| events:",
      sum(dat$OS_status == 1L, na.rm = TRUE), "\n")
  cat("Clinical completeness: age=", sum(!is.na(dat$age_10yr)),
      "; sex=", sum(!is.na(dat$male)),
      "; stage=", sum(stage_ok),
      "; smoking=", sum(!is.na(dat$ever_smoked)), "\n", sep = "")
  cat("Stage coding (GSE68465: pT2+ and pN1+; other cohorts: stage II–IV):\n")
  for (v in stage_terms) {
    cat("  ", v, ": ", sep = "")
    print(table(dat[[v]], useNA = "ifany"))
  }

  qc <- data.frame(
    cohort = cohort, n_score_patients = nrow(dat),
    n_events = sum(dat$OS_status == 1L, na.rm = TRUE),
    age_valid = sum(!is.na(dat$age_10yr)),
    sex_valid = sum(!is.na(dat$male)),
    stage_valid = sum(stage_ok),
    smoking_valid = sum(!is.na(dat$ever_smoked)),
    stage_strategy = if (cohort == "GSE68465")
      "pT2plus_and_pN1plus" else "stage_I_vs_II_to_IV",
    stringsAsFactors = FALSE
  )
  all_qc[[cohort]] <- qc
  patient_rows[[cohort]] <- dat

  # The unadjusted model is a reproducibility check against Script 21.
  # age+sex and stage-adjusted models are the primary planned analyses.
  # Smoking is a complete-case sensitivity; skip if too few events/parameter.
  specifications <- list(
    unadjusted = character(0),
    age_sex_adjusted = c("age_10yr", "male"),
    stage_adjusted = c("age_10yr", "male", stage_terms),
    stage_smoking_sensitivity = c("age_10yr", "male", stage_terms, "ever_smoked")
  )
  for (model_name in names(specifications)) {
    allow_sparse <- model_name != "stage_smoking_sensitivity"
    res <- fit_model(dat, cohort, model_name, specifications[[model_name]],
                     allow_sparse = allow_sparse)
    all_models[[paste(cohort, model_name, sep = "__")]] <- res$summary
    if (!is.null(res$fit)) {
      fit <- res$fit
      betas <- stats::coef(fit)
      se <- sqrt(diag(stats::vcov(fit)))
      all_terms[[paste(cohort, model_name, sep = "__")]] <- data.frame(
        cohort = cohort, model = model_name, term = names(betas),
        coefficient = as.numeric(betas), se = as.numeric(se),
        HR = exp(as.numeric(betas)),
        stringsAsFactors = FALSE
      )
    }
  }
  tmp <- do.call(rbind, all_models[grepl(paste0("^", cohort, "__"), names(all_models))])
  print(tmp[, c("model", "status", "n_analysed", "events_analysed",
                "HR_per_MPCDS_SD", "CI95_lower", "CI95_upper", "p_value")],
        row.names = FALSE, digits = 4)
}

summary <- do.call(rbind, all_models)
row.names(summary) <- NULL
qc <- do.call(rbind, all_qc)
row.names(qc) <- NULL
all_patient_data <- do.call(rbind, patient_rows)
row.names(all_patient_data) <- NULL

write.csv(summary,
          file.path(OUT, "LUAD_MPCDS_clinical_Cox_by_cohort.csv"),
          row.names = FALSE, na = "")
write.csv(qc,
          file.path(OUT, "clinical_harmonization_QC.csv"),
          row.names = FALSE, na = "")
write.csv(all_patient_data,
          file.path(OUT, "clinical_harmonized_patient_variables.csv"),
          row.names = FALSE, na = "")
if (length(all_terms)) {
  tt <- do.call(rbind, all_terms)
  row.names(tt) <- NULL
  write.csv(tt,
            file.path(OUT, "clinical_Cox_all_coefficients.csv"),
            row.names = FALSE, na = "")
}

cat("\n================ SCRIPT 41 FINISHED ================\n")
cat("Output directory:", OUT, "\n")
cat("No meta-analysis has been conducted yet.\n")
cat("Compare unadjusted HRs with existing Script 21 results before pooling.\n")
cat("Notes: stage proxies differ between cohorts; complete-case models;\n")
cat("GSE31210 (35 deaths) is especially sensitive to overfitting.\n")
cat("====================================================\n")
