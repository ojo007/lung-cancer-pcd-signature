# ============================================================================
# Script 39: External LUAD comparison with published PCD prognostic indices
# Study: Cell Death Discovery manuscript corrections (Kankarofi et al.)
#
# Published models, coefficients transcribed from:
#  Wei et al., Journal of Cancer Research and Clinical Oncology (2023),
#    DOI: 10.1007/s00432-023-05000-w; PMID: 37378675.
#  Wang et al., npj Precision Oncology (2024),
#    DOI: 10.1038/s41698-024-00538-5; PMID: 38409471.
#
# This script:
#   (1) Reuses locally cached GEO series matrix files and the exact probe maps
#       saved by Script 19; DOES NOT re-download GEO or rerun Script 19.
#   (2) Uses Script 19's probe-median and log2-if-needed preprocessing.
#   (3) Renames the verified historical alias ERO1L to ERO1A.
#   (4) Computes both published signatures with their published coefficients;
#       DOES NOT fit or select new model weights.
#   (5) Compares C-indices on identical survival-eligible patients with the
#       existing, already validated MPCDS scores saved by Script 21.
#   (6) Gives descriptive paired patient bootstrap C-index difference CIs.
#
# IMPORTANT LIMITATIONS:
#   - Platform/protocol differences may affect absolute published-index values.
#   - GSE72094 participated in Wang et al.'s model-selection procedure; it is
#     NOT an entirely independent evaluation for Wang et al.
#   - No direction-flipping or post-hoc coefficient fitting is performed.
#   - No original results are overwritten. Output is in a new folder.
# ============================================================================

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
OUT <- file.path(ROOT, "04_results", "published_PCD_model_benchmarks")
BOOTSTRAP_REPS <- 500L  # increase to 1000L or 2000L after checking the first run
BOOTSTRAP_SEED <- 20261119L

packages <- c("GEOquery", "Biobase", "survival")
not_installed <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed)) stop("Missing R packages: ", paste(not_installed, collapse = ", "))
if (!dir.exists(ROOT)) stop("Project root does not exist: ", ROOT)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

models <- list(
  Wei_2023_PCDI = c(
    ERO1A = 0.2249586, GAPDH = 0.1896544, TFAP2A = 0.1241279,
    CDK5R1 = 0.1281776, CX3CR1 = -0.1197703
  ),
  Wang_2024_PCDI = c(
    CHEK2 = -0.028463286, KRT18 = 0.099357443,
    GAPDH = 0.156214065, MMP1 = 0.01656202,
    CHRNA5 = 0.007265954, TMPRSS4 = -0.003710263,
    ITGB4 = 0.047539452, CD79A = -0.037942721,
    CTLA4 = -0.049874867, RRM2 = 0.09887913
  )
)
stopifnot(length(models[[1]]) == 5L, length(models[[2]]) == 10L)
all_needed <- unique(unlist(lapply(models, names), use.names = FALSE))
write.csv(
  do.call(rbind, lapply(names(models), function(model) {
    data.frame(model = model, gene = names(models[[model]]),
               published_coefficient = unname(models[[model]]))
  })),
  file.path(OUT, "published_model_coefficients_as_reported.csv"), row.names = FALSE
)

check_file <- function(p) {
  if (!file.exists(p)) stop("Required local file missing: ", p)
  p
}

# Reproduce Script 19: apply log2(x+1) only for intensities on a linear scale.
ensure_log2_like_script19 <- function(x, cohort) {
  finite <- x[is.finite(x)]
  if (!length(finite)) stop(cohort, ": empty/nonfinite expression matrix")
  needs_log2 <- min(finite) >= 0 &&
    (as.numeric(stats::quantile(finite, .99, na.rm = TRUE)) > 100 ||
       max(finite) > 1000)
  if (needs_log2) x <- log2(x + 1)
  cat(cohort, "log2(x+1) applied:", needs_log2, "\n")
  attr(x, "log2_transform_applied") <- needs_log2
  x
}

# Parse the ALREADY-DOWNLOADED local series matrix rather than fetching GEO.
read_cached_geo_matrix <- function(cohort) {
  local_path <- check_file(file.path(
    ROOT, "01_raw_data", "GEO", cohort,
    paste0(cohort, "_series_matrix.txt.gz")
  ))
  args <- list(filename = local_path, GSEMatrix = TRUE, getGPL = FALSE)
  # Newer GEOquery defaults to SummarizedExperiment; ask for ExpressionSet.
  if ("returnType" %in% names(formals(GEOquery::getGEO))) {
    args$returnType <- "ExpressionSet"
  }
  parsed <- do.call(GEOquery::getGEO, args)
  if (is.list(parsed) && !inherits(parsed, "ExpressionSet")) {
    if (length(parsed) != 1L) {
      stop(cohort, ": multiple ExpressionSets in cached file; inspect platforms.")
    }
    parsed <- parsed[[1L]]
  }
  if (!inherits(parsed, "ExpressionSet")) {
    stop(cohort, ": GEOquery did not return an ExpressionSet.")
  }
  x <- Biobase::exprs(parsed)
  storage.mode(x) <- "double"
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(cohort, ": duplicate probe or sample IDs in source matrix")
  }
  x
}

# Replicate the Script 19 strategy: distinct probe-symbol pairs; median if a
# gene is represented by multiple probes. Only 14 unique required genes.
extract_model_genes <- function(x, cohort, sample_ids) {
  map_file <- check_file(file.path(
    ROOT, "02_processed_data", "GEO", "external_validation", cohort,
    paste0(cohort, "_probe_gene_mapping.csv")
  ))
  mp <- read.csv(map_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("probe_id", "gene_symbol") %in% names(mp))) {
    stop(cohort, ": expected probe_id and gene_symbol mapping columns missing")
  }
  mp$probe_id <- trimws(as.character(mp$probe_id))
  mp$original_symbol <- trimws(as.character(mp$gene_symbol))
  mp$gene_symbol <- toupper(mp$original_symbol)
  mp$gene_symbol[mp$gene_symbol %in% c("ERO1L", "ERO1-L", "ERO1LA")] <- "ERO1A"
  mp <- unique(mp[mp$gene_symbol %in% all_needed &
                    mp$probe_id %in% rownames(x),
                  c("probe_id", "gene_symbol", "original_symbol")])
  missing_genes <- setdiff(all_needed, unique(mp$gene_symbol))
  if (length(missing_genes)) {
    stop(cohort, ": mapped probes not present for all required genes: ",
         paste(missing_genes, collapse = ", "))
  }
  missing_ids <- setdiff(sample_ids, colnames(x))
  if (length(missing_ids)) {
    stop(cohort, ": clinical sample IDs absent from GEO matrix: ",
         paste(head(missing_ids, 10L), collapse = ", "))
  }
  # Apply the transformation before gene-level probe collapse, exactly as 19.
  # The source-level range determines whether transformation is necessary.
  x <- ensure_log2_like_script19(x, cohort)
  log2_done <- isTRUE(attr(x, "log2_transform_applied"))
  x <- x[, sample_ids, drop = FALSE]
  out <- matrix(NA_real_, nrow = length(all_needed), ncol = length(sample_ids),
                dimnames = list(all_needed, sample_ids))
  audit <- list()
  for (gene in all_needed) {
    these <- unique(mp$probe_id[mp$gene_symbol == gene])
    if (!length(these)) stop(cohort, ": no probes for ", gene)
    v <- x[these, , drop = FALSE]
    out[gene, ] <- if (nrow(v) == 1L) as.numeric(v[1L, ]) else
      apply(v, 2L, stats::median, na.rm = TRUE)
    audit[[gene]] <- data.frame(
      cohort = cohort, gene = gene, n_probes = length(these),
      probe_ids = paste(these, collapse = ";"),
      original_symbols = paste(unique(mp$original_symbol[mp$gene_symbol == gene]),
                               collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  if (any(!is.finite(out))) stop(cohort, ": nonfinite gene expression after collapsing probes")
  list(expression = out, mapping_audit = do.call(rbind, audit),
       log2_applied = log2_done)
}

c_index <- function(d, score_col) {
  survival::concordance(
    stats::as.formula(paste0("survival::Surv(time, event) ~ ", score_col)),
    data = d, reverse = TRUE
  )$concordance
}

# Paired patient bootstrap CI for each score's C-index difference from MPCDS.
# Reuses the SAME resampled patients for both signatures (not independent CIs).
paired_bootstrap <- function(d, comparator, B, seed) {
  set.seed(seed)
  n <- nrow(d)
  observed <- c_index(d, "MPCDS") - c_index(d, comparator)
  vals <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    dd <- d[idx, , drop = FALSE]
    if (length(unique(dd$event)) < 2L) next
    vals[b] <- tryCatch(
      c_index(dd, "MPCDS") - c_index(dd, comparator),
      error = function(e) NA_real_
    )
  }
  good <- vals[is.finite(vals)]
  if (length(good) < max(50L, ceiling(.8 * B))) {
    stop("Insufficient valid bootstrap draws for ", comparator)
  }
  data.frame(
    comparator = comparator, delta_C_MPCDS_minus_comparator = observed,
    bootstrap_CI_low = as.numeric(quantile(good, .025)),
    bootstrap_CI_high = as.numeric(quantile(good, .975)),
    valid_bootstrap_reps = length(good),
    requested_bootstrap_reps = B,
    stringsAsFactors = FALSE
  )
}

score_files <- c(
  GSE68465 = file.path(
    ROOT, "04_results", "mpcds", "external_validation",
    "GSE68465_MPCDS_patient_scores.csv"),
  GSE72094 = file.path(
    ROOT, "04_results", "mpcds", "external_validation",
    "GSE72094", "GSE72094_MPCDS_patient_scores.csv")
)

all_samples <- list()
all_c <- list()
all_pairs <- list()
all_audit <- list()

for (cohort in names(score_files)) {
  cat("\n============ ", cohort, " ============\n", sep = "")
  existing <- read.csv(check_file(score_files[[cohort]]),
                       stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("sample_id", "OS_status", "MPCDS") %in% names(existing))) {
    stop(cohort, ": existing MPCDS scores lack required columns")
  }
  time_col <- if (cohort == "GSE68465") "OS_months" else "OS_time"
  if (!(time_col %in% names(existing))) {
    stop(cohort, ": expected OS follow-up column missing: ", time_col)
  }
  existing$sample_id <- trimws(as.character(existing$sample_id))
  if (anyDuplicated(existing$sample_id)) stop(cohort, ": repeated sample IDs")
  existing$time <- suppressWarnings(as.numeric(existing[[time_col]]))
  existing$event <- suppressWarnings(as.integer(existing$OS_status))
  existing$MPCDS <- suppressWarnings(as.numeric(existing$MPCDS))
  keep <- is.finite(existing$time) & existing$time > 0 &
    !is.na(existing$event) & existing$event %in% 0:1 &
    is.finite(existing$MPCDS)
  if (!all(keep)) {
    cat("Omitting", sum(!keep), "patients with incomplete OS/MPCDS.\n")
    existing <- existing[keep, , drop = FALSE]
  }
  if (nrow(existing) < 30L || sum(existing$event) < 10L) {
    stop(cohort, ": unexpectedly few eligible survival cases/events")
  }
  cat("Existing MPCDS patients:", nrow(existing),
      "| events:", sum(existing$event), "\n")

  raw_x <- read_cached_geo_matrix(cohort)
  cat("Original GEO matrix:", nrow(raw_x), "probes x", ncol(raw_x), "samples\n")
  recovered <- extract_model_genes(raw_x, cohort, existing$sample_id)
  rm(raw_x)
  invisible(gc())
  x <- recovered$expression
  all_audit[[cohort]] <- recovered$mapping_audit
  cat("Model gene expression:", nrow(x), "genes x", ncol(x), "samples\n")

  d <- data.frame(
    cohort = cohort, sample_id = existing$sample_id,
    time = existing$time, event = existing$event,
    MPCDS = existing$MPCDS, stringsAsFactors = FALSE
  )
  for (model in names(models)) {
    beta <- models[[model]]
    if (!all(names(beta) %in% rownames(x))) stop("Missing a gene in ", model)
    # The coefficient vector is in PUBLISHED order and is not retrained.
    d[[model]] <- as.numeric(crossprod(beta, x[names(beta), , drop = FALSE]))
  }
  if (any(!is.finite(as.matrix(d[, c("MPCDS", names(models)), drop = FALSE])))) {
    stop(cohort, ": nonfinite model score")
  }
  all_samples[[cohort]] <- d
  score_cols <- c("MPCDS", names(models))
  for (m in score_cols) {
    cc <- survival::concordance(
      stats::as.formula(paste0("survival::Surv(time, event) ~ ", m)),
      data = d, reverse = TRUE
    )
    all_c[[paste(cohort, m)]] <- data.frame(
      cohort = cohort, score = m, n = nrow(d), events = sum(d$event),
      C_index = unname(cc$concordance), SE = sqrt(unname(cc$var)),
      CI_lower_approx = max(0, unname(cc$concordance - 1.96 * sqrt(cc$var))),
      CI_upper_approx = min(1, unname(cc$concordance + 1.96 * sqrt(cc$var))),
      stringsAsFactors = FALSE
    )
  }
  for (j in seq_along(models)) {
    m <- names(models)[j]
    cat("Bootstrap comparison MPCDS vs", m, "...\n")
    p <- paired_bootstrap(d, m, BOOTSTRAP_REPS,
                          BOOTSTRAP_SEED + 1000L * match(cohort, names(score_files)) + j)
    p$cohort <- cohort
    # Wang 2024 used GSE72094 during model selection. Do not call it
    # independent external validation for that Wang comparison.
    p$model_selection_overlap_warning <-
      cohort == "GSE72094" && m == "Wang_2024_PCDI"
    all_pairs[[paste(cohort, m)]] <- p
  }
  cat("C-index summary for", cohort, ":\n")
  print(do.call(rbind, all_c[startsWith(names(all_c), paste0(cohort, " "))])[
    , c("score", "n", "events", "C_index")], row.names = FALSE)
}

write.csv(do.call(rbind, all_samples),
          file.path(OUT, "published_PCD_models_patient_scores.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, all_c),
          file.path(OUT, "published_PCD_models_C_indices.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, all_pairs),
          file.path(OUT, "published_PCD_models_paired_delta_C_bootstrap.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, all_audit),
          file.path(OUT, "published_PCD_models_probe_mapping_audit.csv"),
          row.names = FALSE)

cat("\n=========== COMPLETED ===========\n")
cat("All formula genes available; no published model retrained.\n")
cat("Note: Wang et al. used GSE72094 in model selection.\n")
cat("NOTE: These are cross-platform, published-formula comparisons;\n")
cat("      differences may depend on expression normalization.\n")
cat("Output folder:\n", OUT, "\n", sep = "")
