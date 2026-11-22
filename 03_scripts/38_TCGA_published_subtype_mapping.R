# =============================================================
# Script 38 — PCD expression groups versus published TCGA subtypes
# Study: LUAD / LUSC PCD expression patterns
# Purpose: Reproduce the original case-level mapping, provide
#          descriptive cross-tabs, Pearson chi-square tests,
#          Cramer's V and missing-label coverage checks.
# Inputs:  04_results/clustering/TCGA_<histology>_PCD_cluster_assignments.csv
#          Published subtype annotations via TCGAbiolinks::TCGAquery_subtype
# Outputs: 04_results/TCGA_published_subtype_mapping/
# Does not overwrite original clustering, model or figure outputs.
# =============================================================

root <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  stop("Install Bioconductor package 'TCGAbiolinks' before continuing.")
}

out_dir <- file.path(root, "04_results", "TCGA_published_subtype_mapping")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

histologies <- c("LUAD", "LUSC")
summary_rows <- list()

for (cancer in histologies) {

  cat("\n========== ", cancer, " ==========\n", sep = "")

  cluster_file <- file.path(
    root, "04_results", "clustering",
    paste0("TCGA_", cancer, "_PCD_cluster_assignments.csv")
  )
  if (!file.exists(cluster_file)) stop("Missing cluster file: ", cluster_file)

  clusters <- utils::read.csv(cluster_file, stringsAsFactors = FALSE,
                              check.names = FALSE)
  stopifnot(all(c("patient_id", "PCD_cluster") %in% colnames(clusters)))
  clusters$patient_id <- toupper(trimws(as.character(clusters$patient_id)))
  clusters$PCD_cluster <- trimws(as.character(clusters$PCD_cluster))
  if (anyNA(clusters$patient_id) || anyDuplicated(clusters$patient_id)) {
    stop("Invalid or duplicate PCD patient IDs for ", cancer)
  }
  if (!all(clusters$PCD_cluster %in% c("PCD_C1", "PCD_C2"))) {
    stop("Unexpected PCD group labels in ", cancer)
  }
  expected_cluster_n <- if (cancer == "LUAD") 517L else 501L
  if (nrow(clusters) != expected_cluster_n) {
    stop("Unexpected original cohort size for ", cancer)
  }

  published <- TCGAbiolinks::TCGAquery_subtype(tumor = cancer)
  subtype_column <- if (cancer == "LUAD") "expression_subtype" else "Expression.Subtype"
  if (!all(c("patient", subtype_column) %in% colnames(published))) {
    stop("Expected published subtype columns unavailable for ", cancer,
         ". Check the TCGAbiolinks result before proceeding.")
  }

  reference <- data.frame(
    patient_id = toupper(substr(as.character(published$patient), 1L, 12L)),
    subtype_original = trimws(as.character(published[[subtype_column]])),
    stringsAsFactors = FALSE
  )
  reference$subtype_original[
    is.na(reference$subtype_original) |
      reference$subtype_original %in% c("", "NA", "N/A", "NULL")
  ] <- NA_character_

  if (anyNA(reference$patient_id) || anyDuplicated(reference$patient_id)) {
    stop("Missing or duplicate reference patient IDs for ", cancer)
  }

  # Retain original published labels in an explicit column; do not
  # infer any label for patients not present in the TCGA reference.
  matched <- clusters
  matched$published_subtype <- reference$subtype_original[
    match(clusters$patient_id, reference$patient_id)
  ]

  # Consistent display labels, without changing the published labels above.
  matched$published_subtype_display <- matched$published_subtype
  if (cancer == "LUAD") {
    matched$published_subtype_display[
      matched$published_subtype_display == "prox.-inflam"
    ] <- "Proximal-inflammatory"
    matched$published_subtype_display[
      matched$published_subtype_display == "prox.-prolif."
    ] <- "Proximal-proliferative"
    expected_subtypes <- c("TRU", "Proximal-inflammatory", "Proximal-proliferative")
    source_doi <- "10.1038/nature13385"
  } else {
    expected_subtypes <- c("basal", "classical", "primitive", "secretory")
    source_doi <- "10.1038/nature11404"
  }
  observed <- unique(stats::na.omit(matched$published_subtype_display))
  if (!all(observed %in% expected_subtypes)) {
    stop("Unexpected published subtype labels for ", cancer, ": ",
         paste(setdiff(observed, expected_subtypes), collapse = ", "))
  }

  has_label <- !is.na(matched$published_subtype_display)
  n_matched <- sum(has_label)
  expected_n <- if (cancer == "LUAD") 230L else 178L
  if (n_matched != expected_n) {
    stop("Expected ", expected_n, " labelled ", cancer, " patients; got ", n_matched,
         ". Check source version and case matching.")
  }

  annotated <- matched[has_label, , drop = FALSE]
  tab <- table(
    PCD_group = factor(annotated$PCD_cluster, levels = c("PCD_C1", "PCD_C2")),
    Published_TCGA_subtype = factor(
      annotated$published_subtype_display, levels = expected_subtypes
    )
  )
  if (any(colSums(tab) == 0L)) stop("A subtype has zero matched cases in ", cancer)

  test <- stats::chisq.test(tab, correct = FALSE)
  if (min(test$expected) < 5) {
    stop("Low expected cell count for ", cancer,
         "; Pearson chi-square approximation should be reviewed.")
  }
  cramer_v <- sqrt(
    as.numeric(test$statistic) /
      (sum(tab) * min(nrow(tab) - 1L, ncol(tab) - 1L))
  )

  # Coverage check: are published labels disproportionately available
  # for one of the two PCD groups?
  availability <- table(
    PCD_group = factor(matched$PCD_cluster, levels = c("PCD_C1", "PCD_C2")),
    Label_available = factor(has_label, levels = c(FALSE, TRUE))
  )
  coverage_test <- stats::chisq.test(availability, correct = FALSE)

  cat("Original PCD patients:", nrow(matched), "\n")
  cat("Patients with published subtype:", n_matched, "\n")
  cat("Label coverage:", round(100 * n_matched/nrow(matched), 1), "%\n\n")
  cat("Cross-tabulation (counts):\n")
  print(addmargins(tab))
  cat("\nWithin-PCD-group percentages:\n")
  print(round(100 * prop.table(tab, margin = 1), 1))
  cat("\nPearson chi-square:", round(unname(test$statistic), 4),
      "  df:", unname(test$parameter),
      "  p:", format.pval(test$p.value, digits = 4), "\n")
  cat("Cramer's V:", round(cramer_v, 4), "\n")
  cat("Minimum expected count:", round(min(test$expected), 3), "\n")
  cat("Availability-by-PCD-group p:",
      format.pval(coverage_test$p.value, digits = 4), "\n")

  prefix <- file.path(out_dir, paste0("TCGA_", cancer))
  utils::write.csv(matched,
                   paste0(prefix, "_patient_level_reference_mapping.csv"),
                   row.names = FALSE, na = "")
  utils::write.csv(as.data.frame.matrix(tab),
                   paste0(prefix, "_cross_tab_counts.csv"))
  utils::write.csv(as.data.frame.matrix(round(100 * prop.table(tab, 1), 3)),
                   paste0(prefix, "_within_PCD_group_percent.csv"))
  utils::write.csv(as.data.frame.matrix(round(100 * prop.table(tab, 2), 3)),
                   paste0(prefix, "_within_published_subtype_percent.csv"))
  utils::write.csv(as.data.frame.matrix(availability),
                   paste0(prefix, "_label_availability_counts.csv"))

  summary_rows[[cancer]] <- data.frame(
    histology = cancer,
    source_doi = source_doi,
    original_PCD_patients = nrow(matched),
    labelled_patients = n_matched,
    unlabelled_patients = sum(!has_label),
    labelled_percent = 100 * n_matched / nrow(matched),
    chi_square = unname(test$statistic),
    degrees_freedom = unname(test$parameter),
    chi_square_p = test$p.value,
    cramers_v = cramer_v,
    minimum_expected_count = min(test$expected),
    label_availability_p = coverage_test$p.value,
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_rows)
utils::write.csv(summary_table,
                 file.path(out_dir, "TCGA_published_subtype_mapping_summary.csv"),
                 row.names = FALSE)

cat("\n===== FINAL SUMMARY =====\n")
print(summary_table, row.names = FALSE)
cat("\nOutput folder:\n", out_dir, "\n", sep = "")
cat("\nInterpretation caution: matched published reference cases are a subset\n",
    "of the original TCGA cohorts. Association does not imply subtype identity.\n", sep = "")
