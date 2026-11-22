# ============================================================
# 37_MPCDS_random_and_meanZ_benchmark.R
# Professor Zaravinos — major concerns 1 and 3
#
# Benchmarks the 47-gene PCD ridge-Cox model against:
#   (A) an UNWEIGHTED mean of 47 gene-wise z-scores,
#       with a single Cox slope fitted on training patients;
#   (B) equally sized RANDOM 47-gene signatures,
#       each fitted as a ridge Cox model on training patients.
#
# All model fitting, standardisation, lambda selection and
# score-direction fitting occur INSIDE outer training folds.
# Identical stratified outer folds are used for comparisons.
#
# IMPORTANT:
#  - This is a benchmark of a prespecified RIDGE model,
#    whereas original Script 18 searched an alpha grid.
#    Results are NOT guaranteed to reproduce its C=0.604.
#  - Random sets are drawn from sufficiently expressed,
#    variable protein-coding genes (excluding the 47 PCD genes).
#  - PCD candidates were already selected from the full TCGA
#    tumour-vs-normal data before this survival benchmark.
#  - PILOT is a code/data check, NOT final statistical evidence.
#  - Runs only in LUAD, where the published MPCDS was developed.
#  - No existing scripts, model files or figures are modified.
# ============================================================

# Full benchmark enabled after successful pilot input and model checks.
# Set back to "pilot" only if you want to rerun the short test.
RUN_MODE <- "full"  # "pilot" or "full"

ROOT <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"

if (!RUN_MODE %in% c("pilot", "full")) stop("Invalid RUN_MODE")

N_RANDOM <- if (RUN_MODE == "pilot") 3L else 100L
N_REPEATS <- if (RUN_MODE == "pilot") 1L else 20L
OUTER_K <- 5L
INNER_K <- 5L
BASE_SEED <- 12345L

required <- c("dplyr", "survival", "glmnet")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))

path <- function(...) file.path(ROOT, ...)
out_dir <- path("04_results", "mpcds_benchmarks", RUN_MODE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

needed <- c(
  path("02_processed_data", "expression", "TCGA_LUAD_TPM.rds"),
  path("02_processed_data", "expression", "TCGA_LUAD_raw_counts.rds"),
  path("02_processed_data", "expression", "TCGA_LUAD_gene_annotation.rds"),
  path("02_processed_data", "expression", "TCGA_LUAD_sample_metadata.rds"),
  path("02_processed_data", "survival", "TCGA_LUAD_PCD_survival_expression.rds"),
  path("02_processed_data", "survival", "TCGA_LUAD_PCD_survival_metadata.csv"),
  path("04_results", "differential_expression", "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv")
)
if (any(!file.exists(needed))) stop("Missing input files:\n", paste(needed[!file.exists(needed)], collapse = "\n"))

cat("\n=== BENCHMARK INPUT CHECKS (", toupper(RUN_MODE), ") ===\n", sep = "")

ann <- readRDS(needed[3])
meta <- readRDS(needed[4])
tpm <- readRDS(needed[1])
counts <- readRDS(needed[2])
ref_expr <- readRDS(needed[5])
surv_meta <- read.csv(needed[6], stringsAsFactors = FALSE, check.names = FALSE)
pcd <- read.csv(needed[7], stringsAsFactors = FALSE, check.names = FALSE)
pcd_genes <- as.character(pcd$gene_name)

stopifnot(length(pcd_genes) == 47L, !anyDuplicated(pcd_genes))
stopifnot(identical(colnames(tpm), meta$cases))
stopifnot(identical(colnames(counts), meta$cases))
stopifnot(identical(as.character(ann$gene_id), rownames(tpm)))
stopifnot(identical(rownames(counts), rownames(tpm)))
stopifnot(identical(colnames(ref_expr), surv_meta$submitter_id))
stopifnot(all(c("OS_time", "OS_status") %in% names(surv_meta)))
stopifnot(all(is.finite(surv_meta$OS_time)), all(surv_meta$OS_time > 0))
stopifnot(all(surv_meta$OS_status %in% 0:1))

# Exactly reproduce Script 12's primary-tumour/aliquot selection.
meta$library_size <- colSums(counts)
rm(counts)
invisible(gc())

selected <- meta |>
  dplyr::filter(sample_type == "Primary Tumor") |>
  dplyr::group_by(sample.submitter_id) |>
  dplyr::slice_max(order_by = library_size, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::group_by(cases.submitter_id) |>
  dplyr::slice_max(order_by = library_size, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

stopifnot(nrow(selected) == 517L)
selected <- selected[match(surv_meta$submitter_id, selected$cases.submitter_id), , drop = FALSE]
stopifnot(nrow(selected) == 504L, !anyNA(selected$cases))
stopifnot(identical(as.character(selected$cases.submitter_id), surv_meta$submitter_id))

sample_columns <- match(selected$cases, colnames(tpm))
pcd_rows <- match(pcd_genes, ann$gene_name)
stopifnot(!anyNA(sample_columns), !anyNA(pcd_rows))

# This exact comparison catches aliquot/sample or transform mismatches.
x_pcd <- t(log2(tpm[pcd_rows, sample_columns, drop = FALSE] + 1))
colnames(x_pcd) <- pcd_genes
rownames(x_pcd) <- surv_meta$submitter_id
ref <- ref_expr[pcd_genes, surv_meta$submitter_id, drop = FALSE]
max_diff <- max(abs(t(x_pcd) - ref))
cat("Patients:", nrow(x_pcd), "  Events:", sum(surv_meta$OS_status), "\n")
cat("PCD gene count:", ncol(x_pcd), "\n")
cat("Max difference from Script 12 survival expression:", format(max_diff, digits = 8), "\n")
if (!is.finite(max_diff) || max_diff > 1e-7) {
  stop("PCD expression does not reproduce Script 12. Stop and review patient/sample mapping.")
}

# Define a defensible background for random 47-gene panels.
# One row per gene symbol (first occurrence, like Script 12).
mask <- !is.na(ann$gene_name) & nzchar(ann$gene_name) &
  !duplicated(ann$gene_name) &
  ann$gene_type == "protein_coding" &
  !(ann$gene_name %in% pcd_genes)
bg_rows <- which(mask)
if (length(bg_rows) < 47L) stop("Too few annotated protein-coding background genes")

bg_tpm <- tpm[bg_rows, sample_columns, drop = FALSE]
detect_rate <- rowMeans(bg_tpm > 1)
bg_log <- log2(bg_tpm + 1)
bg_sd <- apply(bg_log, 1, stats::sd)
keep <- is.finite(bg_sd) & bg_sd >= 0.2 & detect_rate >= 0.2
bg_symbols <- as.character(ann$gene_name[bg_rows[keep]])
x_background <- t(bg_log[keep, , drop = FALSE])
colnames(x_background) <- bg_symbols
rownames(x_background) <- surv_meta$submitter_id
rm(tpm, bg_tpm, bg_log, ann, meta, selected, ref_expr, ref)
invisible(gc())

stopifnot(ncol(x_background) >= 47L, !anyNA(x_pcd), !anyNA(x_background))
cat("Eligible random-gene background:", ncol(x_background), "genes\n")
cat("No PCD/background overlap:", !any(pcd_genes %in% bg_symbols), "\n")

# Reproduce Script 18's stratified fold generator.
make_folds <- function(status, k, seed) {
  set.seed(seed)
  fold <- integer(length(status))
  for (ev in sort(unique(status))) {
    idx <- sample(which(status == ev))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

# Calculate outer-fold C-index with higher risk implying shorter survival.
c_index <- function(time, status, risk) {
  if (length(unique(risk)) < 2L) return(0.5)
  unname(survival::concordance(
    survival::Surv(time, status) ~ risk, reverse = TRUE
  )$concordance)
}

# Inner CV uses only outer training patients. The SAME inner fold IDs
# are used for PCD and random gene sets within each outer fold.
score_ridge <- function(x_train, x_test, y_train, inner_foldid) {
  cv <- glmnet::cv.glmnet(
    x = x_train,
    y = y_train,
    family = "cox",
    alpha = 0,
    foldid = inner_foldid,
    type.measure = "C",
    standardize = TRUE,
    nlambda = 100,
    cox.ties = "breslow"
  )
  as.numeric(stats::predict(cv, newx = x_test, s = "lambda.min", type = "link"))
}

# The unweighted mean-z score has no gene-specific fitted weights.
# Gene means and SDs are calculated on training patients only;
# a SINGLE Cox slope is learned only on training patients to
# orient the risk direction (avoids unfair fixed sign assumptions).
score_mean_z <- function(x_train, x_test, time_train, status_train) {
  mu <- colMeans(x_train)
  sig <- apply(x_train, 2, stats::sd)
  sig[!is.finite(sig) | sig < 1e-8] <- 1
  train_mean <- rowMeans(sweep(sweep(x_train, 2, mu), 2, sig, "/"))
  test_mean <- rowMeans(sweep(sweep(x_test, 2, mu), 2, sig, "/"))
  fit <- survival::coxph(survival::Surv(time_train, status_train) ~ train_mean)
  beta <- unname(stats::coef(fit)[1])
  if (!is.finite(beta)) stop("Non-finite mean-z Cox slope")
  as.numeric(beta * test_mean)
}

status <- as.integer(surv_meta$OS_status)
time <- as.numeric(surv_meta$OS_time)
y <- survival::Surv(time, status)

# Prepare splits and inner fold IDs once for matched comparisons.
prepare_splits <- function(seed) {
  outer <- make_folds(status, OUTER_K, seed)
  parts <- lapply(seq_len(OUTER_K), function(k) {
    tr <- which(outer != k)
    te <- which(outer == k)
    list(tr = tr, te = te,
         inner = make_folds(status[tr], INNER_K, seed + k))
  })
  parts
}

score_model_on_splits <- function(x, parts, kind = "ridge") {
  fold_c <- numeric(length(parts))
  for (k in seq_along(parts)) {
    pp <- parts[[k]]
    risk <- if (kind == "ridge") {
      score_ridge(x[pp$tr, , drop = FALSE],
                  x[pp$te, , drop = FALSE], y[pp$tr], pp$inner)
    } else {
      score_mean_z(x[pp$tr, , drop = FALSE],
                   x[pp$te, , drop = FALSE], time[pp$tr], status[pp$tr])
    }
    fold_c[k] <- c_index(time[pp$te], status[pp$te], risk)
  }
  fold_c
}

# Fair, paired repeated comparison of PCD ridge versus the simple score.
paired <- vector("list", N_REPEATS)
first_pcd_c <- NA_real_
first_mean_c <- NA_real_
first_parts <- NULL

cat("\n=== PAIRED PCD VS MEAN-Z CROSS-VALIDATION ===\n")
for (r in seq_len(N_REPEATS)) {
  seed <- BASE_SEED + (r - 1L) * 1000L
  parts <- prepare_splits(seed)
  pcd_fold <- score_model_on_splits(x_pcd, parts, "ridge")
  mean_fold <- score_model_on_splits(x_pcd, parts, "mean_z")
  if (r == 1L) {
    first_parts <- parts
    first_pcd_c <- mean(pcd_fold)
    first_mean_c <- mean(mean_fold)
  }
  paired[[r]] <- data.frame(
    repeat_id = r,
    pcd_ridge_C = mean(pcd_fold),
    meanZ_C = mean(mean_fold),
    delta_C = mean(pcd_fold) - mean(mean_fold)
  )
  cat(sprintf("Repeat %d/%d: PCD ridge=%.4f; mean-z=%.4f; difference=%.4f\n",
              r, N_REPEATS, mean(pcd_fold), mean(mean_fold), mean(pcd_fold)-mean(mean_fold)))
}
paired <- do.call(rbind, paired)
write.csv(paired, file.path(out_dir, "PCD_ridge_vs_meanZ_repeat_comparison.csv"), row.names = FALSE)

# Random signatures are sampled independently of survival and
# evaluated on the EXACT SAME first-repeat 5 outer folds.
cat("\n=== RANDOM 47-GENE RIDGE SIGNATURES ===\n")
set.seed(20261119L)
random_idx <- replicate(N_RANDOM, sample.int(ncol(x_background), 47L), simplify = FALSE)
random_tbl <- vector("list", N_RANDOM)
random_genes <- vector("list", N_RANDOM)
for (i in seq_len(N_RANDOM)) {
  chosen <- random_idx[[i]]
  x_rand <- x_background[, chosen, drop = FALSE]
  fold_scores <- score_model_on_splits(x_rand, first_parts, "ridge")
  random_tbl[[i]] <- data.frame(signature_id = i, n_genes = 47L,
                                 mean_outer_C = mean(fold_scores))
  random_genes[[i]] <- data.frame(signature_id = i, gene_symbol = colnames(x_rand))
  cat(sprintf("Random signature %d/%d: held-out C=%.4f\n", i, N_RANDOM, mean(fold_scores)))
}
random_tbl <- do.call(rbind, random_tbl)
random_genes <- do.call(rbind, random_genes)
write.csv(random_tbl, file.path(out_dir, "random_47gene_ridge_C_indices.csv"), row.names = FALSE)
write.csv(random_genes, file.path(out_dir, "random_47gene_signature_membership.csv"), row.names = FALSE)

# Add-one empirical exceedance: meaningful only with sufficient random sets.
p_empirical <- if (RUN_MODE == "pilot") NA_real_ else
  (1 + sum(random_tbl$mean_outer_C >= first_pcd_c)) / (N_RANDOM + 1)
summary_row <- data.frame(
  mode = RUN_MODE,
  patients = length(status),
  events = sum(status),
  pcd_genes = 47L,
  background_genes = ncol(x_background),
  random_signatures = N_RANDOM,
  paired_repeats = N_REPEATS,
  first_repeat_pcd_ridge_C = first_pcd_c,
  first_repeat_meanZ_C = first_mean_c,
  repeated_mean_pcd_ridge_C = mean(paired$pcd_ridge_C),
  repeated_mean_meanZ_C = mean(paired$meanZ_C),
  mean_paired_delta_C = mean(paired$delta_C),
  random_mean_C = mean(random_tbl$mean_outer_C),
  random_max_C = max(random_tbl$mean_outer_C),
  random_empirical_exceedance_p = p_empirical
)
write.csv(summary_row, file.path(out_dir, "benchmark_summary.csv"), row.names = FALSE)

cat("\n=============== BENCHMARK FINISHED ===============\n")
print(summary_row, row.names = FALSE)
cat("Output folder:\n", out_dir, "\n", sep = "")
if (RUN_MODE == "pilot") {
  cat("PILOT ONLY: 3 random signatures give no useful scientific p-value.\n")
  cat("Do not quote these preliminary results in the manuscript.\n")
}
cat("NOTE: fixed-ridge benchmark differs from Script 18's alpha-grid CV.\n")
cat("=================================================\n")
