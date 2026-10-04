# ============================================================
# Script 20: Validate PCD subtypes across external GEO cohorts
# ============================================================
# Fixed-centroid external projection across eight GEO cohorts.
# No external reclustering. Held-out mechanism analyses exclude
# all genes used for projection to reduce circularity.
# ============================================================

source("03_scripts/00_project_config.R")

required_packages <- c("dplyr", "tidyr", "cluster")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

geo_prepared_dir <- file.path(processed_dir, "GEO", "external_validation")
external_results_dir <- file.path(results_dir, "external_validation", "subtypes")
clustering_dir <- file.path(results_dir, "clustering")
expression_dir <- file.path(processed_dir, "expression")
pcd_dir <- file.path(processed_dir, "PCD_genes")
de_dir <- file.path(results_dir, "differential_expression")
ensure_dir(external_results_dir)

manifest_file <- file.path(geo_prepared_dir, "external_GEO_cohort_manifest.csv")
check_files_exist(manifest_file)
cohort_manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("accession", "histology") %in% colnames(cohort_manifest)), nrow(cohort_manifest) == 8)

candidate_file <- file.path(de_dir, "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv")
pcd_master_file <- file.path(pcd_dir, "PCD_master_unique_genes.csv")
luad_assignment_file <- file.path(clustering_dir, "TCGA_LUAD_PCD_cluster_assignments.csv")
lusc_assignment_file <- file.path(clustering_dir, "TCGA_LUSC_PCD_cluster_assignments.csv")
luad_pcd_tpm_file <- file.path(expression_dir, "TCGA_LUAD_PCD_TPM.rds")
lusc_pcd_tpm_file <- file.path(expression_dir, "TCGA_LUSC_PCD_TPM.rds")
luad_counts_file <- file.path(expression_dir, "TCGA_LUAD_raw_counts.rds")
lusc_counts_file <- file.path(expression_dir, "TCGA_LUSC_raw_counts.rds")
luad_metadata_file <- file.path(expression_dir, "TCGA_LUAD_sample_metadata.rds")
lusc_metadata_file <- file.path(expression_dir, "TCGA_LUSC_sample_metadata.rds")

check_files_exist(c(candidate_file, pcd_master_file, luad_assignment_file, lusc_assignment_file,
                    luad_pcd_tpm_file, lusc_pcd_tpm_file, luad_counts_file, lusc_counts_file,
                    luad_metadata_file, lusc_metadata_file))

candidate_genes <- read.csv(candidate_file, stringsAsFactors = FALSE, check.names = FALSE)
pcd_master <- read.csv(pcd_master_file, stringsAsFactors = FALSE, check.names = FALSE)
luad_assignments <- read.csv(luad_assignment_file, stringsAsFactors = FALSE, check.names = FALSE)
lusc_assignments <- read.csv(lusc_assignment_file, stringsAsFactors = FALSE, check.names = FALSE)
luad_pcd_tpm <- readRDS(luad_pcd_tpm_file)
lusc_pcd_tpm <- readRDS(lusc_pcd_tpm_file)
luad_counts <- readRDS(luad_counts_file)
lusc_counts <- readRDS(lusc_counts_file)
luad_metadata <- readRDS(luad_metadata_file)
lusc_metadata <- readRDS(lusc_metadata_file)

model_genes <- unique(candidate_genes$gene_name)
stopifnot(length(model_genes) == 47, nrow(pcd_master) == 296,
          all(model_genes %in% rownames(luad_pcd_tpm)),
          all(model_genes %in% rownames(lusc_pcd_tpm)))

pcd_membership <- pcd_master |>
  dplyr::select(gene_symbol, PCD_types) |>
  tidyr::separate_rows(PCD_types, sep = ";\\s*") |>
  dplyr::rename(PCD_type = PCD_types) |>
  dplyr::distinct(gene_symbol, PCD_type)

select_one_primary_tumor <- function(counts, metadata, assignments) {
  stopifnot(identical(colnames(counts), metadata$cases))
  metadata |>
    dplyr::mutate(patient_id = substr(cases, 1, 12), library_size = colSums(counts)) |>
    dplyr::filter(sample_type == "Primary Tumor", patient_id %in% assignments$patient_id) |>
    dplyr::arrange(patient_id, dplyr::desc(library_size)) |>
    dplyr::group_by(patient_id) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::left_join(assignments |> dplyr::select(patient_id, PCD_cluster), by = "patient_id")
}

luad_selected <- select_one_primary_tumor(luad_counts, luad_metadata, luad_assignments)
lusc_selected <- select_one_primary_tumor(lusc_counts, lusc_metadata, lusc_assignments)
stopifnot(nrow(luad_selected) == 517, nrow(lusc_selected) == 501,
          all(table(luad_selected$PCD_cluster) == c(PCD_C1 = 284, PCD_C2 = 233)),
          all(table(lusc_selected$PCD_cluster) == c(PCD_C1 = 286, PCD_C2 = 215)))

extract_selected_expression <- function(tpm, selected) {
  out <- tpm[, match(selected$cases, colnames(tpm)), drop = FALSE]
  colnames(out) <- selected$patient_id
  log2(out + 1)
}

zscore_genes <- function(x) {
  mu <- rowMeans(x, na.rm = TRUE)
  s <- apply(x, 1, stats::sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0
  z <- sweep(x[keep, , drop = FALSE], 1, mu[keep], "-")
  sweep(z, 1, s[keep], "/")
}

luad_tcga_pcd_z <- zscore_genes(extract_selected_expression(luad_pcd_tpm, luad_selected))
lusc_tcga_pcd_z <- zscore_genes(extract_selected_expression(lusc_pcd_tpm, lusc_selected))

build_centroids <- function(tcga_z, selected, projection_genes) {
  genes <- intersect(projection_genes, rownames(tcga_z))
  z <- tcga_z[genes, selected$patient_id, drop = FALSE]
  labels <- selected$PCD_cluster[match(colnames(z), selected$patient_id)]
  list(
    genes = genes,
    PCD_C1 = rowMeans(z[, labels == "PCD_C1", drop = FALSE]),
    PCD_C2 = rowMeans(z[, labels == "PCD_C2", drop = FALSE])
  )
}

project_external_subtypes <- function(external_z, centroids) {
  genes <- centroids$genes
  x <- external_z[genes, , drop = FALSE]
  c1 <- centroids$PCD_C1[genes]
  c2 <- centroids$PCD_C2[genes]
  cor_c1 <- apply(x, 2, function(v) stats::cor(v, c1, method = "pearson", use = "complete.obs"))
  cor_c2 <- apply(x, 2, function(v) stats::cor(v, c2, method = "pearson", use = "complete.obs"))
  projected <- ifelse(cor_c1 >= cor_c2, "PCD_C1", "PCD_C2")
  data.frame(
    sample_id = colnames(x),
    correlation_PCD_C1 = cor_c1,
    correlation_PCD_C2 = cor_c2,
    correlation_margin = abs(cor_c1 - cor_c2),
    projected_PCD_cluster = projected,
    stringsAsFactors = FALSE
  )
}

calculate_projection_silhouette <- function(external_z, projection) {
  cc <- table(projection$projected_PCD_cluster)
  if (length(cc) < 2 || any(cc < 2)) return(list(mean_silhouette = NA_real_, silhouette = NULL))
  x <- external_z[, projection$sample_id, drop = FALSE]
  d <- stats::dist(t(x), method = "euclidean")
  cl <- as.integer(factor(projection$projected_PCD_cluster, levels = c("PCD_C1", "PCD_C2")))
  sil <- cluster::silhouette(cl, d)
  sil_df <- data.frame(sample_id = projection$sample_id,
                       projected_PCD_cluster = projection$projected_PCD_cluster,
                       silhouette_width = sil[, "sil_width"], stringsAsFactors = FALSE)
  list(mean_silhouette = mean(sil_df$silhouette_width), silhouette = sil_df)
}

calculate_mechanism_scores <- function(z_matrix, membership_table) {
  out <- lapply(sort(unique(membership_table$PCD_type)), function(mech) {
    genes <- membership_table$gene_symbol[membership_table$PCD_type == mech]
    genes <- intersect(unique(genes), rownames(z_matrix))
    if (length(genes) == 0) return(NULL)
    data.frame(sample_id = colnames(z_matrix), PCD_type = mech,
               score = as.numeric(colMeans(z_matrix[genes, , drop = FALSE], na.rm = TRUE)),
               n_genes = length(genes), stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(out)
}

test_mechanism_scores <- function(scores, cluster_data, cluster_column) {
  cluster_small <- cluster_data |> dplyr::select(sample_id, cluster = dplyr::all_of(cluster_column))
  merged <- scores |> dplyr::left_join(cluster_small, by = "sample_id") |> dplyr::filter(!is.na(cluster))
  merged |>
    dplyr::group_by(PCD_type) |>
    dplyr::group_modify(function(.x, .y) {
      c1 <- .x$score[.x$cluster == "PCD_C1"]
      c2 <- .x$score[.x$cluster == "PCD_C2"]
      if (length(c1) < 2 || length(c2) < 2) {
        return(data.frame(n_C1 = length(c1), n_C2 = length(c2), median_C1 = NA_real_, median_C2 = NA_real_, difference_C2_minus_C1 = NA_real_, pvalue = NA_real_))
      }
      wt <- stats::wilcox.test(c2, c1, exact = FALSE)
      data.frame(n_C1 = length(c1), n_C2 = length(c2), median_C1 = stats::median(c1), median_C2 = stats::median(c2),
                 difference_C2_minus_C1 = stats::median(c2) - stats::median(c1), pvalue = wt$p.value)
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      padj = stats::p.adjust(pvalue, method = "BH"),
      direction = dplyr::case_when(difference_C2_minus_C1 > 0 ~ "PCD_C2",
                                   difference_C2_minus_C1 < 0 ~ "PCD_C1",
                                   TRUE ~ "No difference")
    )
}

luad_cluster_data <- data.frame(sample_id = luad_selected$patient_id, PCD_cluster = luad_selected$PCD_cluster)
lusc_cluster_data <- data.frame(sample_id = lusc_selected$patient_id, PCD_cluster = lusc_selected$PCD_cluster)

validate_external_cohort <- function(accession, histology) {
  cat("\n============================================================\n",
      "VALIDATING ", accession, " (", histology, ")\n",
      "============================================================\n", sep = "")
  
  cohort_dir <- file.path(geo_prepared_dir, accession)
  cohort_result_dir <- file.path(external_results_dir, accession)
  ensure_dir(cohort_result_dir)
  
  candidate_z_file <- file.path(cohort_dir, paste0(accession, "_candidate47_expression_z.rds"))
  pcd_z_file <- file.path(cohort_dir, paste0(accession, "_PCD296_expression_z.rds"))
  check_files_exist(c(candidate_z_file, pcd_z_file))
  external_candidate_z <- readRDS(candidate_z_file)
  external_pcd_z <- readRDS(pcd_z_file)
  
  common_projection_genes <- intersect(model_genes, rownames(external_candidate_z))
  if (length(common_projection_genes) < 38) stop(accession, ": fewer than 38/47 projection genes available.")
  
  if (histology == "LUAD") {
    tcga_z <- luad_tcga_pcd_z; tcga_selected <- luad_selected; tcga_cluster_data <- luad_cluster_data
  } else if (histology == "LUSC") {
    tcga_z <- lusc_tcga_pcd_z; tcga_selected <- lusc_selected; tcga_cluster_data <- lusc_cluster_data
  } else stop(accession, ": unsupported histology ", histology)
  
  centroids <- build_centroids(tcga_z, tcga_selected, common_projection_genes)
  projection <- project_external_subtypes(external_candidate_z[common_projection_genes, , drop = FALSE], centroids)
  projection$accession <- accession
  projection$histology <- histology
  projection$n_projection_genes <- length(common_projection_genes)
  
  write.csv(projection, file.path(cohort_result_dir, paste0(accession, "_PCD_subtype_projection.csv")), row.names = FALSE)
  write.csv(data.frame(gene = common_projection_genes,
                       PCD_C1 = centroids$PCD_C1[common_projection_genes],
                       PCD_C2 = centroids$PCD_C2[common_projection_genes]),
            file.path(cohort_result_dir, paste0(accession, "_matched_TCGA_centroids.csv")), row.names = FALSE)
  
  sil <- calculate_projection_silhouette(external_candidate_z[common_projection_genes, , drop = FALSE], projection)
  if (!is.null(sil$silhouette)) write.csv(sil$silhouette, file.path(cohort_result_dir, paste0(accession, "_projection_silhouette.csv")), row.names = FALSE)
  
  heldout_gene_sets <- lapply(sort(unique(pcd_membership$PCD_type)), function(mech) {
    mech_genes <- unique(pcd_membership$gene_symbol[pcd_membership$PCD_type == mech])
    retained <- setdiff(mech_genes, common_projection_genes)
    common <- Reduce(intersect, list(retained, rownames(tcga_z), rownames(external_pcd_z)))
    data.frame(PCD_type = mech, gene = common, stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()
  
  heldout_coverage <- heldout_gene_sets |> dplyr::count(PCD_type, name = "n_common_heldout_genes")
  write.csv(heldout_coverage, file.path(cohort_result_dir, paste0(accession, "_heldout_mechanism_gene_coverage.csv")), row.names = FALSE)
  if (nrow(heldout_coverage) != 5 || any(heldout_coverage$n_common_heldout_genes < 3)) stop(accession, ": insufficient held-out mechanism-gene coverage.")
  
  common_heldout_genes <- unique(heldout_gene_sets$gene)
  matched_membership <- heldout_gene_sets |> dplyr::transmute(gene_symbol = gene, PCD_type) |> dplyr::distinct()
  
  tcga_scores <- calculate_mechanism_scores(tcga_z[common_heldout_genes, , drop = FALSE], matched_membership)
  external_scores <- calculate_mechanism_scores(external_pcd_z[common_heldout_genes, , drop = FALSE], matched_membership)
  external_cluster_data <- projection |> dplyr::transmute(sample_id, projected_PCD_cluster)
  tcga_tests <- test_mechanism_scores(tcga_scores, tcga_cluster_data, "PCD_cluster")
  external_tests <- test_mechanism_scores(external_scores, external_cluster_data, "projected_PCD_cluster")
  
  write.csv(external_scores, file.path(cohort_result_dir, paste0(accession, "_heldout_PCD_mechanism_scores.csv")), row.names = FALSE)
  write.csv(external_tests, file.path(cohort_result_dir, paste0(accession, "_heldout_PCD_mechanism_tests.csv")), row.names = FALSE)
  
  concordance <- dplyr::inner_join(
    tcga_tests |> dplyr::select(PCD_type, TCGA_difference = difference_C2_minus_C1, TCGA_padj = padj, TCGA_direction = direction),
    external_tests |> dplyr::select(PCD_type, GEO_difference = difference_C2_minus_C1, GEO_padj = padj, GEO_direction = direction),
    by = "PCD_type"
  ) |>
    dplyr::mutate(direction_concordant = TCGA_direction == GEO_direction,
                  strict_FDR_replication = direction_concordant & TCGA_padj < 0.05 & GEO_padj < 0.05)
  
  write.csv(concordance, file.path(cohort_result_dir, paste0(accession, "_heldout_mechanism_concordance.csv")), row.names = FALSE)
  
  n_c1 <- sum(projection$projected_PCD_cluster == "PCD_C1")
  n_c2 <- sum(projection$projected_PCD_cluster == "PCD_C2")
  n_direction <- sum(concordance$direction_concordant)
  n_strict <- sum(concordance$strict_FDR_replication)
  strict_mechanisms <- paste(concordance$PCD_type[concordance$strict_FDR_replication], collapse = ";")
  if (strict_mechanisms == "") strict_mechanisms <- "None"
  
  summary_row <- data.frame(
    accession = accession, histology = histology, n_samples = ncol(external_candidate_z),
    n_projection_genes = length(common_projection_genes), projected_PCD_C1 = n_c1, projected_PCD_C2 = n_c2,
    mean_correlation_margin = mean(projection$correlation_margin), median_correlation_margin = stats::median(projection$correlation_margin),
    mean_silhouette = sil$mean_silhouette, mechanisms_directionally_concordant = n_direction,
    mechanisms_tested = nrow(concordance), directional_concordance_percent = 100 * n_direction / nrow(concordance),
    strict_FDR_replicated_mechanisms = n_strict, replicated_mechanisms = strict_mechanisms,
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_row, file.path(cohort_result_dir, paste0(accession, "_external_subtype_validation_summary.csv")), row.names = FALSE)
  
  cat("Samples: ", summary_row$n_samples, "\n",
      "Projection genes: ", summary_row$n_projection_genes, "/47\n",
      "Projected PCD_C1: ", n_c1, "\n",
      "Projected PCD_C2: ", n_c2, "\n",
      "Mean silhouette: ", round(summary_row$mean_silhouette, 4), "\n",
      "Held-out mechanism directional concordance: ", n_direction, "/", nrow(concordance), "\n",
      "Strict FDR-replicated mechanisms: ", strict_mechanisms, "\n", sep = "")
  
  summary_row
}

validation_results <- vector("list", nrow(cohort_manifest))
for (i in seq_len(nrow(cohort_manifest))) {
  accession <- cohort_manifest$accession[i]
  histology <- cohort_manifest$histology[i]
  validation_results[[i]] <- tryCatch(
    validate_external_cohort(accession, histology),
    error = function(e) {
      warning(accession, " validation failed: ", conditionMessage(e))
      data.frame(accession = accession, histology = histology, n_samples = NA_integer_, n_projection_genes = NA_integer_,
                 projected_PCD_C1 = NA_integer_, projected_PCD_C2 = NA_integer_, mean_correlation_margin = NA_real_,
                 median_correlation_margin = NA_real_, mean_silhouette = NA_real_, mechanisms_directionally_concordant = NA_integer_,
                 mechanisms_tested = NA_integer_, directional_concordance_percent = NA_real_, strict_FDR_replicated_mechanisms = NA_integer_,
                 replicated_mechanisms = paste0("FAILED: ", conditionMessage(e)), stringsAsFactors = FALSE)
    }
  )
}

validation_summary <- dplyr::bind_rows(validation_results)
write.csv(validation_summary, file.path(external_results_dir, "external_GEO_PCD_subtype_validation_summary.csv"), row.names = FALSE)

histology_summary <- validation_summary |>
  dplyr::group_by(histology) |>
  dplyr::summarise(
    n_cohorts = dplyr::n(),
    n_successful = sum(!grepl("^FAILED:", replicated_mechanisms)),
    total_samples = sum(n_samples, na.rm = TRUE),
    median_projection_genes = stats::median(n_projection_genes, na.rm = TRUE),
    median_silhouette = stats::median(mean_silhouette, na.rm = TRUE),
    median_directional_concordance_percent = stats::median(directional_concordance_percent, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(histology_summary, file.path(external_results_dir, "external_GEO_PCD_subtype_histology_summary.csv"), row.names = FALSE)

concordance_files <- list.files(external_results_dir, pattern = "_heldout_mechanism_concordance\\.csv$", recursive = TRUE, full.names = TRUE)
mechanism_replication <- lapply(concordance_files, function(f) {
  accession <- basename(dirname(f))
  histology <- cohort_manifest$histology[match(accession, cohort_manifest$accession)]
  dat <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  dat$accession <- accession
  dat$histology <- histology
  dat
}) |> dplyr::bind_rows()

write.csv(mechanism_replication, file.path(external_results_dir, "external_GEO_heldout_mechanism_all_cohorts.csv"), row.names = FALSE)

mechanism_frequency <- mechanism_replication |>
  dplyr::group_by(histology, PCD_type) |>
  dplyr::summarise(n_cohorts = dplyr::n(),
                   n_direction_concordant = sum(direction_concordant),
                   n_strict_FDR_replication = sum(strict_FDR_replication),
                   .groups = "drop")
write.csv(mechanism_frequency, file.path(external_results_dir, "external_GEO_heldout_mechanism_replication_frequency.csv"), row.names = FALSE)

n_failed <- sum(grepl("^FAILED:", validation_summary$replicated_mechanisms))

cat("\n============================================================\n",
    "SCRIPT 20 EXTERNAL PCD SUBTYPE VALIDATION SUMMARY\n",
    "============================================================\n", sep = "")
print(validation_summary)
cat("\nHistology-level summary:\n")
print(histology_summary)
cat("\nHeld-out mechanism replication frequency:\n")
print(mechanism_frequency)
cat("\nCohorts attempted: ", nrow(cohort_manifest), "\n",
    "Cohorts successfully validated: ", nrow(cohort_manifest) - n_failed, "\n",
    "Failed cohorts: ", n_failed, "\n",
    "Outputs: ", external_results_dir, "\n", sep = "")

if (n_failed > 0) {
  cat("\nIMPORTANT: One or more external subtype validations failed. Resolve them before interpretation.\n")
} else {
  cat("\nAll eight external GEO cohorts were projected using fixed TCGA centroids without external reclustering.\n")
}

cat("============================================================\n",
    "\n========================================\n",
    "SCRIPT 20 COMPLETED SUCCESSFULLY\n",
    "========================================\n", sep = "")