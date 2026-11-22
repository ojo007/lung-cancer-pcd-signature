# ==============================================================================
# Script 40: Audit clinical covariate availability in external LUAD GEO cohorts
# Project: Kankarofi et al., Cell Death Discovery manuscript corrections
#
# PURPOSE ONLY: Inventory (not model) age/sex/stage/T/N/smoking information for
# clinically adjusted external survival analyses. Does not fit any survival
# model, overwrite earlier results, or edit the source datasets.
#
# Input: phenotype_selected.csv from Script 19 and published MPCDS patient
# scores from Script 21 for GSE68465/GSE72094/GSE31210/GSE50081.
# Output: New audit CSVs under 04_results/mpcds_clinical_covariate_audit.
# ==============================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
OUT <- file.path(ROOT, "04_results", "mpcds_clinical_covariate_audit")
COHORTS <- c("GSE68465", "GSE72094", "GSE31210", "GSE50081")

stopifnot(dir.exists(ROOT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

clinical_pattern <- paste(
  c("age", "sex", "gender", "stage", "tnm", "nodal", "smok",
    "pack", "histolog", "grade", "therap", "treat", "recurr",
    "patholog", "characteristic", "clinical", "adjuvant", "resect"),
  collapse = "|"
)

missing_tokens <- c("", "NA", "N/A", "NULL", "NONE", "UNKNOWN",
                    "NOT AVAILABLE", "NOT REPORTED", "NOT APPLICABLE",
                    "UNAVAILABLE", "MISSING", "--", "---")

get_missing <- function(x) {
  y <- toupper(trimws(as.character(x)))
  is.na(x) | is.na(y) | y %in% missing_tokens
}

source_score_path <- function(cohort) {
  if (cohort == "GSE68465") {
    return(file.path(ROOT, "04_results", "mpcds", "external_validation",
                     "GSE68465_MPCDS_patient_scores.csv"))
  }
  file.path(ROOT, "04_results", "mpcds", "external_validation",
            cohort, paste0(cohort, "_MPCDS_patient_scores.csv"))
}

cohort_rows <- list()
field_rows <- list()

for (cohort in COHORTS) {
  cat("\n============== ", cohort, " ==============\n", sep = "")

  phenotype_file <- file.path(ROOT, "02_processed_data", "GEO",
                              "external_validation", cohort,
                              paste0(cohort, "_phenotype_selected.csv"))
  score_file <- source_score_path(cohort)

  if (!file.exists(phenotype_file) || !file.exists(score_file)) {
    cat("Missing input file(s):\n")
    if (!file.exists(phenotype_file)) cat("  ", phenotype_file, "\n", sep = "")
    if (!file.exists(score_file)) cat("  ", score_file, "\n", sep = "")
    cohort_rows[[cohort]] <- data.frame(
      cohort = cohort, input_status = "MISSING_INPUT", n_scores = NA_integer_,
      n_matched = NA_integer_, n_events = NA_integer_, n_candidate_fields = NA_integer_
    )
    next
  }

  pheno <- read.csv(phenotype_file, check.names = FALSE,
                    stringsAsFactors = FALSE)
  scores <- read.csv(score_file, check.names = FALSE,
                     stringsAsFactors = FALSE)

  stopifnot("geo_accession" %in% names(pheno),
            "sample_id" %in% names(scores),
            "OS_status" %in% names(scores))

  ids <- trimws(as.character(pheno$geo_accession))
  score_ids <- trimws(as.character(scores$sample_id))
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop(cohort, ": ambiguous or missing phenotype GEO identifiers")
  }
  if (anyNA(score_ids) || anyDuplicated(score_ids)) {
    stop(cohort, ": ambiguous or missing score sample identifiers")
  }
  ix <- match(score_ids, ids)
  matched <- !is.na(ix)
  nmatch <- sum(matched)
  cat("Survival-score patients:", nrow(scores),
      "| matched phenotype records:", nmatch, "\n")
  cat("Events in score file:",
      sum(suppressWarnings(as.numeric(scores$OS_status)) == 1, na.rm = TRUE),
      "\n")

  candidate_cols <- names(pheno)[grepl(clinical_pattern, names(pheno),
                                      ignore.case = TRUE)]
  candidate_cols <- setdiff(candidate_cols, "geo_accession")
  cat("Clinical/characteristics candidate columns:",
      length(candidate_cols), "\n")

  cohort_rows[[cohort]] <- data.frame(
    cohort = cohort, input_status = "OK", n_scores = nrow(scores),
    n_matched = nmatch,
    n_events = sum(suppressWarnings(as.numeric(scores$OS_status)) == 1,
                   na.rm = TRUE),
    n_candidate_fields = length(candidate_cols)
  )

  if (!length(candidate_cols)) {
    cat("No clinical columns identified by their names.\n")
    cat("All phenotype column names:\n")
    print(head(names(pheno), 40L))
    next
  }

  details <- lapply(candidate_cols, function(nm) {
    values <- pheno[[nm]][ix[matched]]
    miss <- get_missing(values)
    distinct <- unique(trimws(as.character(values[!miss])))
    distinct <- distinct[!is.na(distinct)]
    samples <- paste(head(substr(distinct, 1L, 65L), 3L), collapse = " | ")
    data.frame(
      cohort = cohort, field = nm,
      n_available = sum(!miss),
      n_missing = sum(miss),
      available_percent = if (nmatch > 0) round(100 * sum(!miss) / nmatch, 1)
                          else NA_real_,
      n_distinct = length(distinct), example_values = samples,
      stringsAsFactors = FALSE
    )
  })
  dat <- do.call(rbind, details)
  field_rows[[cohort]] <- dat

  # Print the most informative fields first, without dumping entire phenotype.
  ord <- order(-dat$n_available, dat$field)
  preview <- dat[head(ord, 35L),
                 c("field", "n_available", "available_percent",
                   "n_distinct", "example_values"), drop = FALSE]
  print(preview, row.names = FALSE, right = FALSE)
  if (nrow(dat) > 35L) {
    cat("...", nrow(dat) - 35L,
        "additional candidate fields are saved in the CSV.\n")
  }
}

summary <- do.call(rbind, cohort_rows)
write.csv(summary, file.path(OUT, "clinical_covariate_cohort_coverage.csv"),
          row.names = FALSE)
if (length(field_rows)) {
  fields <- do.call(rbind, field_rows)
  write.csv(fields, file.path(OUT, "clinical_covariate_field_inventory.csv"),
            row.names = FALSE)
}

cat("\n=============== AUDIT COMPLETED ===============\n")
print(summary, row.names = FALSE)
cat("No clinical models have been fitted.\n")
cat("Output directory: ", OUT, "\n", sep = "")
