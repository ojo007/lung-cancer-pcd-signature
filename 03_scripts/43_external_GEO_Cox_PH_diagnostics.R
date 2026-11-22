# ==============================================================================
# SCRIPT 43 — EXTERNAL LUAD MPCDS COX PROPORTIONAL-HAZARDS DIAGNOSTICS
# Kankarofi et al. — Cell Death Discovery scientific corrections
#
# INPUTS: Script 41 verified patient-level clinical variables and Cox summaries.
# PURPOSE:
#   1. Reproduce the 12 Script 41 Cox model estimates without changing them.
#   2. Apply Schoenfeld-residual proportional-hazards tests (cox.zph) to
#      MPCDS and all model covariates, including a global test.
#   3. Screen for repeated GEO sample-accession IDs across the four datasets.
#
# IMPORTANT: Unique GEO sample accessions cannot exclude patient overlap
# across studies. Confirm study provenance and shared individuals separately.
# This is a diagnostic script; do not interpret PH p-values as guarantees
# that the hazard ratio is constant in time. No old files are overwritten.
# ==============================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
BASE <- file.path(ROOT, "04_results", "mpcds_clinically_adjusted_external")
OUT <- file.path(BASE, "proportional_hazards_diagnostics")
COHORTS <- c("GSE68465", "GSE72094", "GSE31210", "GSE50081")
MODELS <- c("unadjusted", "age_sex_adjusted", "stage_adjusted")

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Install the survival package first: install.packages('survival')")
}
patient_file <- file.path(BASE, "clinical_harmonized_patient_variables.csv")
summary_file <- file.path(BASE, "LUAD_MPCDS_clinical_Cox_by_cohort.csv")
stopifnot(file.exists(patient_file), file.exists(summary_file))
patients <- read.csv(patient_file, stringsAsFactors = FALSE, check.names = FALSE)
previous <- read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
req <- c("cohort", "sample_id", "OS_time", "OS_status", "MPCDS_z",
         "age_10yr", "male", "t2plus", "npos", "advanced_stage")
missing <- setdiff(req, names(patients))
if (length(missing)) stop("Missing patient data fields: ", paste(missing, collapse = ", "))

if (anyNA(patients$sample_id) || anyDuplicated(paste(patients$cohort, patients$sample_id))) {
  stop("Missing or duplicate patient IDs within a cohort")
}
if (!all(COHORTS %in% unique(patients$cohort))) {
  stop("Expected GEO cohorts not all present in Script 41 patient CSV")
}

get_vars <- function(cohort, model) {
  if (model == "unadjusted") return("MPCDS_z")
  if (model == "age_sex_adjusted") return(c("MPCDS_z", "age_10yr", "male"))
  if (cohort == "GSE68465") {
    return(c("MPCDS_z", "age_10yr", "male", "t2plus", "npos"))
  }
  c("MPCDS_z", "age_10yr", "male", "advanced_stage")
}

results <- list()
reproduction <- list()
for (cohort in COHORTS) {
  cat("\n============ ", cohort, " ============\n", sep = "")
  whole <- patients[patients$cohort == cohort, , drop = FALSE]
  for (model in MODELS) {
    vars <- get_vars(cohort, model)
    needed <- c("OS_time", "OS_status", vars)
    complete <- stats::complete.cases(whole[, needed, drop = FALSE]) &
      is.finite(whole$OS_time) & whole$OS_time > 0 &
      whole$OS_status %in% c(0L, 1L)
    dat <- whole[complete, , drop = FALSE]
    n <- nrow(dat)
    events <- sum(dat$OS_status == 1L)

    prev <- previous[previous$cohort == cohort & previous$model == model &
                       previous$status == "OK", , drop = FALSE]
    if (nrow(prev) != 1L) stop("Missing/duplicate successful Script 41 result: ", cohort, "/", model)
    if (n != prev$n_analysed || events != prev$events_analysed) {
      stop("Script 41 sample/event mismatch for ", cohort, "/", model)
    }

    formula <- stats::as.formula(paste(
      "survival::Surv(OS_time, OS_status) ~", paste(vars, collapse = " + ")
    ))
    fit <- survival::coxph(formula, data = dat, ties = "efron", x = TRUE,
                           y = TRUE, singular.ok = FALSE)
    hr <- exp(unname(stats::coef(fit)["MPCDS_z"]))
    difference <- abs(hr - prev$HR_per_MPCDS_SD)
    if (!is.finite(difference) || difference > 1e-6) {
      stop(sprintf("Cox fit does not reproduce Script 41 for %s/%s: delta=%.8f",
                   cohort, model, difference))
    }
    reproduction[[length(reproduction) + 1L]] <- data.frame(
      cohort = cohort, model = model, n = n, events = events,
      previous_MPCDS_HR = prev$HR_per_MPCDS_SD,
      reproduced_MPCDS_HR = hr, absolute_difference = difference,
      stringsAsFactors = FALSE
    )

    ph <- tryCatch(
      survival::cox.zph(fit, transform = "km", terms = FALSE, global = TRUE),
      error = function(e) e
    )
    if (inherits(ph, "error")) {
      warning(cohort, "/", model, ": cox.zph failed: ", conditionMessage(ph))
      results[[length(results) + 1L]] <- data.frame(
        cohort = cohort, model = model, n = n, events = events,
        tested_term = "ERROR", chisquare = NA_real_, df = NA_real_,
        p_value = NA_real_, PH_flag_p_lt_0_05 = NA,
        note = conditionMessage(ph), stringsAsFactors = FALSE
      )
    } else {
      t <- as.data.frame(ph$table)
      term <- rownames(t)
      results[[length(results) + 1L]] <- data.frame(
        cohort = cohort, model = model, n = n, events = events,
        tested_term = term, chisquare = t$chisq, df = t$df,
        p_value = t$p, PH_flag_p_lt_0_05 = t$p < 0.05,
        note = "Schoenfeld test; KM time transformation",
        stringsAsFactors = FALSE
      )
      mpcds_p <- t$p[match("MPCDS_z", term)]
      global_p <- t$p[match("GLOBAL", term)]
      cat(sprintf("%-20s n=%3d events=%3d | MPCDS PH p=%.4g | global PH p=%.4g\n",
                  model, n, events, mpcds_p, global_p))
    }
  }
}

reproduced <- do.call(rbind, reproduction)
ph_table <- do.call(rbind, results)
stopifnot(nrow(reproduced) == length(COHORTS) * length(MODELS))

# A conservative accession-only screen. Distinct GSM identifiers across GEO
# studies are NOT evidence of distinct biological individuals.
ids_by_cohort <- split(unique(patients[, c("cohort", "sample_id")]), patients$cohort)
cohort_combinations <- utils::combn(COHORTS, 2, simplify = FALSE)
accession_screen <- do.call(rbind, lapply(cohort_combinations, function(pair) {
  a <- unique(patients$sample_id[patients$cohort == pair[1]])
  b <- unique(patients$sample_id[patients$cohort == pair[2]])
  hits <- intersect(a, b)
  data.frame(cohort_a = pair[1], cohort_b = pair[2],
             same_GEO_accession_count = length(hits),
             repeated_GEO_accessions = paste(hits, collapse = "; "),
             interpretation = "Accession screen only; patient overlap remains unverified",
             stringsAsFactors = FALSE)
}))

# Save only in a separate output directory.
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(reproduced, file.path(OUT, "Cox_Script41_reproduction.csv"), row.names = FALSE)
write.csv(ph_table, file.path(OUT, "Cox_Schoenfeld_PH_tests.csv"), row.names = FALSE)
write.csv(accession_screen, file.path(OUT, "GEO_sample_accession_overlap_screen.csv"), row.names = FALSE)
writeLines(c(
  "Schoenfeld PH tests are diagnostics, not proof of constant hazards.",
  "Primary model: stage_adjusted; include MPCDS-specific and global tests.",
  "PH violations can make a single Cox hazard ratio an average over follow-up.",
  "Accession-only screen does NOT establish absence of patient overlap.",
  "Cross-cohort identity requires provenance assessment and, when available, shared clinical/sample identifiers.",
  "No prior Script 21, 41 or 42 results were modified."
), file.path(OUT, "interpretation_cautions.txt"))
cat("\n============ SCRIPT 43 FINISHED ============\n")
cat("Stage-adjusted MPCDS and global PH results:\n")
print(ph_table[ph_table$model == "stage_adjusted" &
                 ph_table$tested_term %in% c("MPCDS_z", "GLOBAL"),
               c("cohort", "tested_term", "chisquare", "p_value")], row.names = FALSE)
cat("All 12 Cox estimates reproduced within tolerance.\n")
cat("Sample ID overlap screen is NOT a patient identity audit.\n")
cat("Output folder: ", OUT, "\n", sep = "")
