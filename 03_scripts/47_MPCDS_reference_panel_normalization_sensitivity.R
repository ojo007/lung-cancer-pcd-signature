# =======================================================================
# SCRIPT 47 — MPCDS cohort-wise z-score transportability sensitivity
# =======================================================================
# Goal:
#   Examine dependence of EXTERNAL MPCDS discrimination on the reference
#   population used for per-gene z-standardization (Script 19/21).
#
#   A fixed fraction of expression samples (randomly drawn WITHOUT using
#   outcome information) provides each gene's reference mean and SD.
#   The remaining samples are scored; only eligible survival patients
#   among those remaining samples enter the C-index comparison.
#
#   For a FAIR paired comparison, the published full-cohort z-scored
#   MPCDS is also evaluated on those exact same remaining patients.
#
# This is a reference-panel stability / sensitivity assessment, NOT
# deployment validation, and NOT an independent external cohort.
# No model coefficients are refitted and no original result is altered.
# =======================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
COHORTS <- c("GSE68465", "GSE72094", "GSE31210", "GSE50081")
REFERENCE_FRACTIONS <- c(0.20, 0.40, 0.60)
N_REPEATS <- 50L
BASE_SEED <- 470000L

in_file <- function(...) file.path(ROOT, ...)
OUT <- in_file("04_results", "mpcds_normalization_transportability")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The survival package is required.")
}

weights_file <- in_file(
  "04_results", "mpcds", "external_validation",
  "TCGA_LUAD_MPCDS_standardized_weights_47.csv"
)
if (!file.exists(weights_file)) stop("Missing standardized weights: ", weights_file)
wt_table <- read.csv(weights_file, stringsAsFactors = FALSE, check.names = FALSE)
required_wt <- c("gene", "standardized_weight")
if (!all(required_wt %in% names(wt_table))) {
  stop("Standardized weights file does not have expected columns.")
}
if (anyDuplicated(wt_table$gene) || nrow(wt_table) != 47L) {
  stop("Expected 47 unique gene weights.")
}
weights <- as.numeric(wt_table$standardized_weight)
names(weights) <- as.character(wt_table$gene)
if (any(!is.finite(weights))) stop("Nonfinite TCGA-derived weights.")

score_file <- function(cohort) {
  if (cohort == "GSE68465") {
    in_file("04_results", "mpcds", "external_validation",
            "GSE68465_MPCDS_patient_scores.csv")
  } else {
    in_file("04_results", "mpcds", "external_validation", cohort,
            paste0(cohort, "_MPCDS_patient_scores.csv"))
  }
}

calc_c <- function(time, status, score) {
  d <- data.frame(time = as.numeric(time), event = as.integer(status),
                  risk = as.numeric(score))
  if (nrow(d) < 20L || sum(d$event == 1L) < 5L ||
      sum(d$event == 0L) < 5L) return(NA_real_)
  as.numeric(survival::concordance(
    survival::Surv(time, event) ~ risk,
    data = d, reverse = TRUE
  )$concordance)
}

all_qc <- list()
all_repeats <- list()

for (cohort_idx in seq_along(COHORTS)) {
  cohort <- COHORTS[cohort_idx]
  cat("\n=============== ", cohort, " ===============\n", sep = "")
  prep <- in_file("02_processed_data", "GEO", "external_validation", cohort)
  log_file <- file.path(prep, paste0(cohort, "_candidate47_expression_log2.rds"))
  z_file <- file.path(prep, paste0(cohort, "_candidate47_expression_z.rds"))
  sf <- score_file(cohort)
  for (file in c(log_file, z_file, sf)) {
    if (!file.exists(file)) stop(cohort, ": missing input ", file)
  }

  x <- as.matrix(readRDS(log_file))
  orig_z <- as.matrix(readRDS(z_file))
  d <- read.csv(sf, stringsAsFactors = FALSE, check.names = FALSE)
  time_col <- if (cohort == "GSE68465") "OS_months" else "OS_time"
  needed <- c("sample_id", time_col, "OS_status", "MPCDS")
  if (!all(needed %in% names(d))) {
    stop(cohort, ": missing score columns: ",
         paste(setdiff(needed, names(d)), collapse = ", "))
  }
  d$sample_id <- trimws(as.character(d$sample_id))
  d$time <- as.numeric(d[[time_col]])
  d$event <- as.integer(as.character(d$OS_status))
  if (anyDuplicated(d$sample_id) || anyNA(d$sample_id) ||
      any(!is.finite(d$time)) || any(d$time <= 0) ||
      anyNA(d$event) || any(!d$event %in% c(0L, 1L)) ||
      any(!is.finite(d$MPCDS))) {
    stop(cohort, ": invalid/duplicate survival-score data.")
  }
  if (anyDuplicated(colnames(x)) || anyDuplicated(rownames(x)) ||
      anyDuplicated(colnames(orig_z)) || anyDuplicated(rownames(orig_z))) {
    stop(cohort, ": duplicated expression genes or sample IDs.")
  }
  if (!identical(colnames(x), colnames(orig_z))) {
    stop(cohort, ": log expression and z expression have different sample order.")
  }
  if (!all(d$sample_id %in% colnames(x))) {
    stop(cohort, ": survival patients not found in expression matrix.")
  }

  genes <- intersect(names(weights), rownames(orig_z))
  if (length(genes) < 38L || !all(genes %in% rownames(x))) {
    stop(cohort, ": available model genes below threshold or missing in log data.")
  }
  # Exactly Script 21's available model genes and fixed standardized weights.
  x <- x[genes, , drop = FALSE]
  orig_z <- orig_z[genes, , drop = FALSE]
  w <- weights[genes]
  if (any(!is.finite(x)) || any(!is.finite(orig_z))) {
    stop(cohort, ": missing or nonfinite gene expression.")
  }

  # Reproduce the Script 19 gene-wise z-scoring directly from the
  # log-expression matrix across ALL selected expression samples.
  gene_mu <- rowMeans(x)
  gene_sd <- apply(x, 1, stats::sd)
  if (any(!is.finite(gene_sd)) || any(gene_sd <= 0)) {
    stop(cohort, ": gene SD is zero/nonfinite.")
  }
  z_reconstructed <- sweep(sweep(x, 1, gene_mu, "-"), 1, gene_sd, "/")
  max_z_difference <- max(abs(z_reconstructed - orig_z))
  if (max_z_difference > 1e-7) {
    stop(cohort, ": original GEO z-scores do not reproduce; max discrepancy = ",
         format(max_z_difference, digits = 8))
  }
  score_from_z <- as.numeric(crossprod(w, orig_z))
  names(score_from_z) <- colnames(orig_z)
  max_score_difference <- max(abs(score_from_z[d$sample_id] - d$MPCDS))
  if (max_score_difference > 1e-6) {
    stop(cohort, ": saved MPCDS fails reproduction; max discrepancy = ",
         format(max_score_difference, digits = 8))
  }

  baseline_c <- calc_c(d$time, d$event, d$MPCDS)
  cat("Expression samples:", ncol(x), "| survival patients:", nrow(d),
      "| events:", sum(d$event), "| model genes:", length(genes), "\n")
  cat(sprintf("Full-cohort reference reproduced: max z diff %.3g; max score diff %.3g; C = %.5f\n",
              max_z_difference, max_score_difference, baseline_c))
  all_qc[[cohort]] <- data.frame(
    cohort = cohort, expression_samples = ncol(x),
    survival_patients = nrow(d), events = sum(d$event),
    n_model_genes = length(genes),
    max_z_difference = max_z_difference,
    max_score_difference = max_score_difference,
    existing_C_index = baseline_c,
    stringsAsFactors = FALSE
  )

  n_all <- ncol(x)
  for (frac in REFERENCE_FRACTIONS) {
    n_ref <- round(n_all * frac)
    if (n_ref < 12L || n_ref >= n_all - 20L) {
      stop(cohort, ": reference fraction creates too few samples.")
    }
    results <- vector("list", N_REPEATS)
    for (r in seq_len(N_REPEATS)) {
      set.seed(BASE_SEED + cohort_idx * 100000L +
               as.integer(frac * 1000) * 100L + r)
      ref_idx <- sample.int(n_all, size = n_ref, replace = FALSE)
      ref_ids <- colnames(x)[ref_idx]
      test <- d[!(d$sample_id %in% ref_ids), , drop = FALSE]
      ref_mu <- rowMeans(x[, ref_idx, drop = FALSE])
      ref_sd <- apply(x[, ref_idx, drop = FALSE], 1, stats::sd)
      if (any(!is.finite(ref_sd)) || any(ref_sd <= 0)) {
        stop(cohort, ": reference SD zero/nonfinite at fraction=", frac,
             ", repeat=", r)
      }
      xt <- x[, test$sample_id, drop = FALSE]
      zt <- sweep(sweep(xt, 1, ref_mu, "-"), 1, ref_sd, "/")
      reference_score <- as.numeric(crossprod(w, zt))
      original_score <- score_from_z[test$sample_id]
      c_reference <- calc_c(test$time, test$event, reference_score)
      c_original <- calc_c(test$time, test$event, original_score)
      if (!is.finite(c_reference) || !is.finite(c_original)) {
        stop(cohort, ": too few comparable survival observations at repeat ", r)
      }
      rank_rho <- suppressWarnings(stats::cor(
        reference_score, original_score, method = "spearman"
      ))
      results[[r]] <- data.frame(
        cohort = cohort, reference_fraction = frac, repeat_id = r,
        n_reference_expression = n_ref, n_test_survival = nrow(test),
        n_test_events = sum(test$event),
        C_full_cohort_z_same_test = c_original,
        C_reference_panel_z_same_test = c_reference,
        delta_C_reference_minus_full = c_reference - c_original,
        spearman_scores_same_test = rank_rho,
        stringsAsFactors = FALSE
      )
    }
    tab <- do.call(rbind, results)
    all_repeats[[paste(cohort, frac)]] <- tab
    cat(sprintf("Reference %2.0f%% (n=%d): mean paired delta C %+0.4f; mean rank rho %.4f; test events %d–%d\n",
                100 * frac, n_ref,
                mean(tab$delta_C_reference_minus_full),
                mean(tab$spearman_scores_same_test),
                min(tab$n_test_events), max(tab$n_test_events)))
  }
}

qc <- do.call(rbind, all_qc)
rep_table <- do.call(rbind, all_repeats)
rownames(qc) <- rownames(rep_table) <- NULL

summary_table <- do.call(rbind, lapply(split(
  rep_table, paste(rep_table$cohort, rep_table$reference_fraction)
), function(d) {
  v <- d$delta_C_reference_minus_full
  data.frame(
    cohort = d$cohort[1], reference_fraction = d$reference_fraction[1],
    n_reference_expression = d$n_reference_expression[1],
    n_repeats = nrow(d),
    mean_test_patients = mean(d$n_test_survival),
    mean_test_events = mean(d$n_test_events),
    mean_full_cohort_C_same_test = mean(d$C_full_cohort_z_same_test),
    mean_reference_panel_C_same_test = mean(d$C_reference_panel_z_same_test),
    mean_paired_delta_C = mean(v),
    q025_delta_C = unname(stats::quantile(v, 0.025)),
    q975_delta_C = unname(stats::quantile(v, 0.975)),
    proportion_abs_delta_C_below_0_01 = mean(abs(v) < 0.01),
    mean_spearman_scores = mean(d$spearman_scores_same_test),
    stringsAsFactors = FALSE
  )
}))
rownames(summary_table) <- NULL
summary_table <- summary_table[order(match(summary_table$cohort, COHORTS),
                                     summary_table$reference_fraction), , drop = FALSE]

write.csv(qc, file.path(OUT, "normalization_reproduction_QC.csv"), row.names = FALSE)
write.csv(rep_table, file.path(OUT, "reference_panel_repeat_results.csv"), row.names = FALSE)
write.csv(summary_table, file.path(OUT, "reference_panel_summary.csv"), row.names = FALSE)
writeLines(c(
  "Script 47: cohort-wise normalization/reference-panel sensitivity",
  "Script 19/21 original gene-wise z-scores and fixed TCGA standardized weights are reproduced before sensitivity testing.",
  "Unsupervised expression-only reference panel, without using survival outcomes to select reference samples.",
  "All comparisons of C-index use IDENTICAL held-out eligible survival cases for original vs reference-panel scaling.",
  "Genes/coefficients are never refitted. Reference fractions 20%, 40%, 60%, each assessed over 50 random splits.",
  "Reference samples and evaluation samples are disjoint; reference panels can contain samples without survival information.",
  "The ranges across 50 repeats describe reference-panel and test-composition variability, not bootstrap patient confidence intervals.",
  "This DOES NOT validate transfer across RNA-seq vs microarray or establish calibration for individual-patient clinical deployment.",
  "Original complete-cohort normalization is transductive: gene means/SDs use all expression samples, including evaluated cases.",
  "No original Script 19/21 results are changed."
), file.path(OUT, "methods_and_limitations.txt"))

cat("\n===== SCRIPT 47 COMPLETE =====\n")
print(qc, row.names = FALSE)
print(summary_table[, c("cohort", "reference_fraction", "n_repeats",
                         "mean_paired_delta_C", "q025_delta_C", "q975_delta_C",
                         "mean_spearman_scores")], row.names = FALSE, digits = 4)
cat("Output folder:", OUT, "\n")
