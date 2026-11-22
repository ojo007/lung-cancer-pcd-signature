# ============================================================
# 35_proliferation_purity_sensitivity.R
# Exploratory sensitivity analyses for PCD expression groups
# and the LUAD MPCDS using proliferation and an ESTIMATE-based
# RELATIVE tumour-purity proxy (not measured tumour purity).
#
# Inputs: existing outputs from scripts 14, 15 and 18, plus
# TCGA full TPM, raw-count, annotation, and sample metadata RDS.
# Outputs: NEW folder 04_results/proliferation_purity_sensitivity
#
# Does not overwrite existing files or refit the MPCDS score.
# The Cox models are exploratory, development-cohort models,
# NOT cross-validated or externally validated performance.
# ============================================================

project_root <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"

if (!dir.exists(project_root)) stop("Project root not found. Edit project_root at the top of this script.")
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("R package 'survival' is required. Run install.packages('survival') separately, then rerun.")
}

expression_dir <- file.path(project_root, "02_processed_data", "expression")
cluster_dir    <- file.path(project_root, "04_results", "clustering")
immune_dir     <- file.path(project_root, "04_results", "immune")
mpcds_dir      <- file.path(project_root, "04_results", "mpcds")
candidate_file <- file.path(project_root, "04_results", "differential_expression",
                            "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv")
output_dir     <- file.path(project_root, "04_results", "proliferation_purity_sensitivity")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Ten pre-specified, exploratory cell-cycle/proliferation markers.
# This is a proxy panel, NOT a validated clinical proliferation assay.
markers <- c("MKI67", "PCNA", "TOP2A", "MCM2", "MCM5",
             "CCNB1", "CDK1", "AURKA", "BIRC5", "UBE2C")

need_file <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  path
}

standardise <- function(x, what) {
  if (any(!is.finite(x))) stop(what, " contains missing/non-finite values")
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) stop(what, " has zero/non-finite variance")
  as.numeric((x - mean(x)) / s)
}

get_lm_row <- function(model, term, cohort, outcome, adjustment) {
  z <- summary(model)$coefficients
  if (!term %in% rownames(z)) stop("Missing model term ", term)
  data.frame(histology = cohort, outcome = outcome, adjustment = adjustment,
             contrast = "PCD_C2 vs PCD_C1", estimate = unname(z[term, "Estimate"]),
             standard_error = unname(z[term, "Std. Error"]),
             p_value = unname(z[term, "Pr(>|t|)"]),
             sample_size = nobs(model), row.names = NULL)
}

# Reads one histology at a time, preserving the original script-14
# selection rule: Primary Tumor; for duplicate aliquots/samples and
# patients, keep the one with the greatest raw-count library size.
process_histology <- function(cancer, pcd_genes) {
  message("\n=== Processing ", cancer, " ===")
  tpm <- readRDS(need_file(file.path(expression_dir, paste0("TCGA_", cancer, "_TPM.rds"))))
  counts <- readRDS(need_file(file.path(expression_dir, paste0("TCGA_", cancer, "_raw_counts.rds"))))
  ann <- readRDS(need_file(file.path(expression_dir, paste0("TCGA_", cancer, "_gene_annotation.rds"))))
  meta <- readRDS(need_file(file.path(expression_dir, paste0("TCGA_", cancer, "_sample_metadata.rds"))))
  cluster <- read.csv(need_file(file.path(cluster_dir,
                           paste0("TCGA_", cancer, "_PCD_cluster_assignments.csv"))),
                      stringsAsFactors = FALSE)
  estimate <- read.csv(need_file(file.path(immune_dir,
                            paste0("TCGA_", cancer, "_ESTIMATE_scores.csv"))),
                       stringsAsFactors = FALSE)

  stopifnot(is.matrix(tpm), is.matrix(counts),
            identical(dim(tpm), dim(counts)),
            identical(rownames(tpm), rownames(counts)),
            identical(colnames(tpm), colnames(counts)),
            identical(colnames(tpm), as.character(meta$cases)),
            identical(rownames(tpm), as.character(ann$gene_id)),
            nrow(ann) == nrow(tpm),
            all(c("gene_name", "gene_id") %in% names(ann)),
            all(c("sample_type", "sample.submitter_id", "cases.submitter_id", "cases") %in% names(meta)),
            all(c("patient_id", "PCD_cluster") %in% names(cluster)),
            all(c("patient_id", "ESTIMATEScore") %in% names(estimate)),
            !anyDuplicated(cluster$patient_id),
            !anyDuplicated(estimate$patient_id))

  expected <- if (cancer == "LUAD") 517L else 501L
  if (nrow(cluster) != expected || nrow(estimate) != expected) {
    stop(cancer, ": expected ", expected, " matched cluster/ESTIMATE patients")
  }
  if (!setequal(cluster$patient_id, estimate$patient_id)) {
    stop(cancer, ": ESTIMATE patient identifiers do not match clustering")
  }

  # If a marker also belongs to the PCD candidate set, omit it to
  # prevent direct gene reuse between the proliferation and MPCDS scores.
  eligible_markers <- setdiff(markers, pcd_genes)
  overlap <- intersect(markers, pcd_genes)
  if (length(overlap)) message("Excluded overlapping PCD markers: ", paste(overlap, collapse = ", "))
  eligible_markers <- intersect(eligible_markers, unique(as.character(ann$gene_name)))
  if (length(eligible_markers) < 8L) stop(cancer, ": <8 nonoverlapping proliferation markers available")
  message("Proliferation markers used (", length(eligible_markers), "): ",
          paste(eligible_markers, collapse = ", "))

  n_primary_tumour_profiles <- sum(meta$sample_type == "Primary Tumor")
  meta$library_size <- as.numeric(colSums(counts))
  meta$original_order <- seq_len(nrow(meta))
  selected <- meta[meta$sample_type == "Primary Tumor" &
                     meta$cases.submitter_id %in% cluster$patient_id, , drop = FALSE]
  if (nrow(selected) == 0L) stop("No matching primary tumour samples")
  # Highest library count first; stable first-in-input tie breaker.
  selected <- selected[order(-selected$library_size, selected$original_order), , drop = FALSE]
  selected <- selected[!duplicated(selected$sample.submitter_id), , drop = FALSE]
  selected <- selected[!duplicated(selected$cases.submitter_id), , drop = FALSE]
  selected <- selected[match(cluster$patient_id, selected$cases.submitter_id), , drop = FALSE]
  if (anyNA(selected$cases)) stop(cancer, ": missing tumour profiles for clustered patients")
  if (!identical(as.character(selected$cases.submitter_id), as.character(cluster$patient_id))) {
    stop("Selected tumour profiles not aligned with cluster patient IDs")
  }
  if (anyDuplicated(selected$cases)) stop("Duplicate selected sample barcodes")

  # Only 10 rows of the matrix are needed from this point onwards.
  # Duplicate symbols, if any, are averaged within each marker.
  marker_log2 <- matrix(NA_real_, nrow = length(eligible_markers),
                        ncol = nrow(selected), dimnames = list(eligible_markers, cluster$patient_id))
  for (j in seq_along(eligible_markers)) {
    ix <- which(as.character(ann$gene_name) == eligible_markers[j])
    v <- tpm[ix, selected$cases, drop = FALSE]
    if (any(!is.finite(v)) || any(v < 0)) stop("Invalid TPM values for ", eligible_markers[j])
    marker_log2[j, ] <- colMeans(log2(v + 1))
  }
  rm(tpm, counts, ann, meta)
  invisible(gc())

  # Z-score each marker over the eligible tumour patients and then
  # average. Finally standardise composite within the histology.
  marker_z <- t(apply(marker_log2, 1, standardise, what = "marker expression"))
  rownames(marker_z) <- rownames(marker_log2)
  proliferation <- colMeans(marker_z)
  proliferation_z <- standardise(proliferation, "composite proliferation")

  es <- estimate$ESTIMATEScore[match(cluster$patient_id, estimate$patient_id)]
  # Larger ESTIMATE score = more immune/stromal admixture = lower
  # relative tumour purity. Minus sign reverses orientation.
  # RNA-seq scores cannot be turned into validated absolute purity.
  purity_proxy_z <- standardise(-es, "negative ESTIMATE score")

  dat <- data.frame(patient_id = cluster$patient_id,
                    PCD_cluster = factor(cluster$PCD_cluster, levels = c("PCD_C1", "PCD_C2")),
                    proliferation_raw = as.numeric(proliferation),
                    proliferation_z = as.numeric(proliferation_z),
                    ESTIMATEScore = as.numeric(es),
                    relative_purity_proxy_z = as.numeric(purity_proxy_z),
                    stringsAsFactors = FALSE)
  if (anyNA(dat$PCD_cluster)) stop("Unrecognised PCD cluster label")
  if (anyNA(dat)) stop("Missing generated covariates")

  write.csv(dat, file.path(output_dir, paste0(cancer, "_proliferation_ESTIMATE_covariates.csv")), row.names = FALSE)
  write.csv(data.frame(marker = eligible_markers, overlaps_47_PCD_genes = FALSE),
            file.path(output_dir, paste0(cancer, "_proliferation_marker_audit.csv")), row.names = FALSE)

  prolif_raw <- lm(proliferation_z ~ PCD_cluster, dat)
  prolif_adj <- lm(proliferation_z ~ PCD_cluster + relative_purity_proxy_z, dat)
  purity_raw <- lm(relative_purity_proxy_z ~ PCD_cluster, dat)
  purity_adj <- lm(relative_purity_proxy_z ~ PCD_cluster + proliferation_z, dat)
  associations <- rbind(
    get_lm_row(prolif_raw, "PCD_clusterPCD_C2", cancer, "proliferation_z", "unadjusted"),
    get_lm_row(prolif_adj, "PCD_clusterPCD_C2", cancer, "proliferation_z", "adjusted for relative purity proxy"),
    get_lm_row(purity_raw, "PCD_clusterPCD_C2", cancer, "relative_purity_proxy_z", "unadjusted"),
    get_lm_row(purity_adj, "PCD_clusterPCD_C2", cancer, "relative_purity_proxy_z", "adjusted for proliferation"))

  list(data = dat, associations = associations,
       audit = data.frame(histology = cancer, full_tumor_profiles = n_primary_tumour_profiles,
                          cluster_patients = nrow(cluster), matched_estimate = sum(cluster$patient_id %in% estimate$patient_id),
                          selected_unique_patients = nrow(selected), marker_count = length(eligible_markers),
                          marker_overlap_removed = paste(overlap, collapse = "; "), row.names = NULL))
}

candidate <- read.csv(need_file(candidate_file), stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene_name" %in% names(candidate) || nrow(candidate) != 47L) {
  stop("Candidate file must contain exactly 47 genes with column gene_name")
}
pcd_genes <- unique(as.character(candidate$gene_name))
if (length(pcd_genes) != 47L) stop("Duplicate PCD gene symbols")

luad <- process_histology("LUAD", pcd_genes)
lusc <- process_histology("LUSC", pcd_genes)

write.csv(rbind(luad$associations, lusc$associations),
          file.path(output_dir, "PCD_cluster_proliferation_ESTIMATE_associations.csv"), row.names = FALSE)

# LUAD-only MPCDS exploratory survival association. 504 patients
# with existing original scores, event and survival-time data.
mpcds <- read.csv(need_file(file.path(mpcds_dir, "TCGA_LUAD_MPCDS_patient_scores.csv")),
                  stringsAsFactors = FALSE)
stopifnot(nrow(mpcds) == 504L,
          all(c("patient_id", "OS_time", "OS_status", "MPCDS_z") %in% names(mpcds)),
          !anyDuplicated(mpcds$patient_id),
          all(mpcds$patient_id %in% luad$data$patient_id))

ix <- match(mpcds$patient_id, luad$data$patient_id)
surv_data <- data.frame(mpcds[, c("patient_id", "OS_time", "OS_status", "MPCDS_z"), drop = FALSE],
                        luad$data[ix, c("PCD_cluster", "proliferation_z", "ESTIMATEScore", "relative_purity_proxy_z"),
                                  drop = FALSE])
if (anyNA(surv_data)) stop("Missing fields in LUAD survival analysis")
if (any(!is.finite(surv_data$OS_time)) || any(surv_data$OS_time <= 0)) {
  stop("Survival times must be finite and >0")
}
if (!all(surv_data$OS_status %in% c(0, 1))) stop("Expected 0/1 event coding")
if (sum(surv_data$OS_status) < 15) stop("Too few survival events to fit adjustment models")
write.csv(surv_data, file.path(output_dir, "LUAD_MPCDS_adjustment_patient_level.csv"), row.names = FALSE)

formulas <- list(
  MPCDS_only = survival::Surv(OS_time, OS_status) ~ MPCDS_z,
  MPCDS_plus_proliferation = survival::Surv(OS_time, OS_status) ~ MPCDS_z + proliferation_z,
  MPCDS_plus_ESTIMATE_proxy = survival::Surv(OS_time, OS_status) ~ MPCDS_z + relative_purity_proxy_z,
  MPCDS_plus_both = survival::Surv(OS_time, OS_status) ~ MPCDS_z + proliferation_z + relative_purity_proxy_z
)
models <- lapply(formulas, function(f) survival::coxph(f, data = surv_data, ties = "efron", x = TRUE))

get_cox_row <- function(model, model_name, variable) {
  s <- summary(model)
  if (!variable %in% rownames(s$coefficients)) stop("Missing Cox variable ", variable)
  data.frame(model = model_name, covariate = variable,
             hazard_ratio = unname(s$conf.int[variable, "exp(coef)"]),
             lower_95_CI = unname(s$conf.int[variable, "lower .95"]),
             upper_95_CI = unname(s$conf.int[variable, "upper .95"]),
             p_value = unname(s$coefficients[variable, "Pr(>|z|)"]),
             N = model$n, events = model$nevent, row.names = NULL)
}

cox_out <- do.call(rbind, lapply(names(models), function(nm) {
  mm <- models[[nm]]
  do.call(rbind, lapply(names(stats::coef(mm)), function(term) get_cox_row(mm, nm, term)))
}))
write.csv(cox_out, file.path(output_dir, "LUAD_MPCDS_Cox_adjustment_models.csv"), row.names = FALSE)

ph <- survival::cox.zph(models$MPCDS_plus_both)
ph_tab <- data.frame(term = rownames(ph$table), ph$table, check.names = FALSE)
write.csv(ph_tab, file.path(output_dir, "LUAD_adjusted_Cox_PH_diagnostics.csv"), row.names = FALSE)

spearman <- lapply(c("proliferation_z", "relative_purity_proxy_z"), function(v) {
  r <- suppressWarnings(stats::cor.test(surv_data$MPCDS_z, surv_data[[v]], method = "spearman", exact = FALSE))
  data.frame(variable = v, rho = unname(r$estimate), p_value = r$p.value,
             N = nrow(surv_data))
})
write.csv(do.call(rbind, spearman),
          file.path(output_dir, "LUAD_MPCDS_covariate_correlations.csv"), row.names = FALSE)

write.csv(rbind(luad$audit, lusc$audit),
          file.path(output_dir, "input_alignment_and_marker_audit.csv"), row.names = FALSE)

notes <- c(
  "EXPLORATORY PROLIFERATION AND ESTIMATE ADJUSTMENTS",
  "================================================",
  "PCD groups: 517 LUAD and 501 LUSC previously defined patients.",
  "Survival models: existing 504-patient LUAD MPCDS scores only.",
  "Proliferation: mean of within-histology per-gene Z-scores of log2(TPM+1)",
  paste("Pre-specified marker panel:", paste(markers, collapse = ", ")),
  "Any markers overlapping the 47 PCD candidates are excluded automatically.",
  "Tumour purity: relative proxy = standardised NEGATIVE ESTIMATEScore.",
  "A higher proxy value indicates comparatively less inferred immune/stromal admixture.",
  "It is NOT a measured or calibrated absolute tumour-purity percentage.",
  "Do NOT use the Affymetrix-based ESTIMATE-to-purity cosine formula with TCGA RNA-seq scores.",
  "Exploratory Cox models are fitted in the existing development cohort.",
  "They do NOT replace the original nested cross-validation or external validation.",
  "Cox proportional hazards diagnostic for the fully adjusted model is included.",
  "Read model outputs and diagnostics before drawing biological or clinical conclusions.",
  "Scripts 00-34 and all original research outputs remain unchanged.",
  "Obtain author review and consider an independent tumour-purity estimate (e.g. ABSOLUTE)",
  "if a genuine purity-adjusted analysis is needed for the revised manuscript."
)
writeLines(notes, con = file.path(output_dir, "ANALYSIS_NOTES.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

cat("\n=============== COMPLETED ===============\n")
cat("LUAD clustering patients:", nrow(luad$data), "\n")
cat("LUSC clustering patients:", nrow(lusc$data), "\n")
cat("LUAD MPCDS survival patients:", nrow(surv_data), "\n")
cat("LUAD survival events:", sum(surv_data$OS_status), "\n")
cat("\nLUAD MPCDS coefficients from survival models:\n")
print(cox_out[cox_out$covariate == "MPCDS_z", , drop = FALSE], row.names = FALSE)
cat("\nFully adjusted Cox PH diagnostic:\n")
print(ph_tab, row.names = FALSE)
cat("\nResults saved in:\n", output_dir, "\n")
cat("=========================================\n")
