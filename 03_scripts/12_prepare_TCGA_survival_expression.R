# ============================================================
# 12_prepare_TCGA_survival_expression.R
#
# Purpose:
# Prepare patient-level PCD expression matrices for survival
# analysis in TCGA-LUAD and TCGA-LUSC.
#
# Input:
# - Full TPM expression matrices
# - Full raw-count matrices
# - RNA-seq sample metadata
# - Final survival metadata
# - 47 concordant LUAD/LUSC PCD candidates
#
# Output:
# 02_processed_data/survival/
#
# TCGA_LUAD_PCD_survival_expression.rds
# TCGA_LUAD_PCD_survival_metadata.csv
# TCGA_LUSC_PCD_survival_expression.rds
# TCGA_LUSC_PCD_survival_metadata.csv
#
# Expression transformation:
# log2(TPM + 1)
#
# Sample-selection rule:
# 1. Primary Tumor only
# 2. One aliquot per biological sample:
#    highest raw-count library size
# 3. One tumor sample per patient:
#    highest raw-count library size
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

expression_processed_dir <- file.path(
  processed_dir,
  "expression"
)

clinical_processed_dir <- file.path(
  processed_dir,
  "clinical"
)

de_results_dir <- file.path(
  results_dir,
  "differential_expression"
)

survival_processed_dir <- file.path(
  processed_dir,
  "survival"
)

ensure_dir(survival_processed_dir)


# ------------------------------------------------------------
# 2. Required package
# ------------------------------------------------------------

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

suppressPackageStartupMessages(
  library(dplyr)
)

# Explicit dplyr namespaces are used below to avoid masking by
# Bioconductor packages in interactive R sessions.


# ------------------------------------------------------------
# 2A. Required inputs
# ------------------------------------------------------------

required_input_files <- c(
  file.path(
    de_results_dir,
    "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv"
  ),
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_raw_counts.rds"
  ),
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_TPM.rds"
  ),
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_raw_counts.rds"
  ),
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_TPM.rds"
  ),
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_gene_annotation.rds"
  ),
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_gene_annotation.rds"
  ),
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_sample_metadata.rds"
  ),
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_sample_metadata.rds"
  ),
  file.path(
    clinical_processed_dir,
    "TCGA_LUAD_survival_final.rds"
  ),
  file.path(
    clinical_processed_dir,
    "TCGA_LUSC_survival_final.rds"
  )
)

check_files_exist(
  required_input_files,
  label = "Script 12 input file(s)"
)


# ------------------------------------------------------------
# 3. Load 47 concordant PCD candidates
# ------------------------------------------------------------

concordant_pcd_candidates <- read.csv(
  file.path(
    de_results_dir,
    "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv"
  ),
  stringsAsFactors = FALSE
)

survival_candidate_genes <-
  concordant_pcd_candidates$gene_name

stopifnot(
  length(survival_candidate_genes) == 47
)

stopifnot(
  !anyDuplicated(
    survival_candidate_genes
  )
)


# ------------------------------------------------------------
# 4. Load RNA-seq matrices
# ------------------------------------------------------------

luad_counts <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_raw_counts.rds"
  )
)

luad_tpm <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_TPM.rds"
  )
)

lusc_counts <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_raw_counts.rds"
  )
)

lusc_tpm <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_TPM.rds"
  )
)


# ------------------------------------------------------------
# 5. Load gene annotations
# ------------------------------------------------------------

luad_gene_annotation <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_gene_annotation.rds"
  )
)

lusc_gene_annotation <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_gene_annotation.rds"
  )
)


# ------------------------------------------------------------
# 6. Load RNA-seq sample metadata
# ------------------------------------------------------------

luad_sample_metadata <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_sample_metadata.rds"
  )
)

lusc_sample_metadata <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_sample_metadata.rds"
  )
)


# ------------------------------------------------------------
# 7. Validate raw matrix / metadata alignment
# ------------------------------------------------------------

stopifnot(
  identical(
    colnames(luad_counts),
    luad_sample_metadata$cases
  )
)

stopifnot(
  identical(
    colnames(luad_tpm),
    luad_sample_metadata$cases
  )
)

stopifnot(
  identical(
    colnames(lusc_counts),
    lusc_sample_metadata$cases
  )
)

stopifnot(
  identical(
    colnames(lusc_tpm),
    lusc_sample_metadata$cases
  )
)


# ------------------------------------------------------------
# 8. Calculate raw-count library sizes
# ------------------------------------------------------------

luad_sample_metadata$library_size <-
  colSums(luad_counts)

lusc_sample_metadata$library_size <-
  colSums(lusc_counts)


# ============================================================
# PART A — SELECT ONE LUAD TUMOR PROFILE PER PATIENT
# ============================================================


# ------------------------------------------------------------
# 9. LUAD Primary Tumor only
# ------------------------------------------------------------

luad_tumor_metadata <- luad_sample_metadata |>
  dplyr::filter(
    sample_type == "Primary Tumor"
  )


# ------------------------------------------------------------
# 10. One LUAD aliquot per biological sample
# ------------------------------------------------------------

luad_tumor_one_aliquot <- luad_tumor_metadata |>
  dplyr::group_by(
    sample.submitter_id
  ) |>
  dplyr::slice_max(
    order_by = library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()


# ------------------------------------------------------------
# 11. One LUAD tumor sample per patient
# ------------------------------------------------------------

luad_tumor_one_patient <- luad_tumor_one_aliquot |>
  dplyr::group_by(
    cases.submitter_id
  ) |>
  dplyr::slice_max(
    order_by = library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

stopifnot(
  nrow(luad_tumor_one_patient) == 517
)

stopifnot(
  !anyDuplicated(
    luad_tumor_one_patient$cases.submitter_id
  )
)


# ============================================================
# PART B — SELECT ONE LUSC TUMOR PROFILE PER PATIENT
# ============================================================


# ------------------------------------------------------------
# 12. LUSC Primary Tumor only
# ------------------------------------------------------------

lusc_tumor_metadata <- lusc_sample_metadata |>
  dplyr::filter(
    sample_type == "Primary Tumor"
  )


# ------------------------------------------------------------
# 13. One LUSC aliquot per biological sample
# ------------------------------------------------------------

lusc_tumor_one_aliquot <- lusc_tumor_metadata |>
  dplyr::group_by(
    sample.submitter_id
  ) |>
  dplyr::slice_max(
    order_by = library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()


# ------------------------------------------------------------
# 14. One LUSC tumor sample per patient
# ------------------------------------------------------------

lusc_tumor_one_patient <- lusc_tumor_one_aliquot |>
  dplyr::group_by(
    cases.submitter_id
  ) |>
  dplyr::slice_max(
    order_by = library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

stopifnot(
  nrow(lusc_tumor_one_patient) == 501
)

stopifnot(
  !anyDuplicated(
    lusc_tumor_one_patient$cases.submitter_id
  )
)


# ============================================================
# PART C — LOAD SURVIVAL DATA
# ============================================================


# ------------------------------------------------------------
# 15. LUAD/LUSC final survival tables
# ------------------------------------------------------------

luad_survival <- readRDS(
  file.path(
    clinical_processed_dir,
    "TCGA_LUAD_survival_final.rds"
  )
)

lusc_survival <- readRDS(
  file.path(
    clinical_processed_dir,
    "TCGA_LUSC_survival_final.rds"
  )
)

stopifnot(
  nrow(luad_survival) == 509
)

stopifnot(
  nrow(lusc_survival) == 496
)


# ------------------------------------------------------------
# 16. Identify expression-survival overlap
# ------------------------------------------------------------

luad_survival_overlap <- intersect(
  luad_tumor_one_patient$cases.submitter_id,
  luad_survival$submitter_id
)

lusc_survival_overlap <- intersect(
  lusc_tumor_one_patient$cases.submitter_id,
  lusc_survival$submitter_id
)

stopifnot(
  length(luad_survival_overlap) == 504
)

stopifnot(
  length(lusc_survival_overlap) == 493
)


# ============================================================
# PART D — LUAD SURVIVAL EXPRESSION
# ============================================================


# ------------------------------------------------------------
# 17. LUAD survival metadata
# ------------------------------------------------------------

luad_survival_metadata <- luad_survival |>
  dplyr::filter(
    submitter_id %in%
      luad_survival_overlap
  )


# ------------------------------------------------------------
# 18. Map LUAD patients to selected RNA aliquots
# ------------------------------------------------------------

luad_expression_map <- luad_tumor_one_patient |>
  dplyr::filter(
    cases.submitter_id %in%
      luad_survival_metadata$submitter_id
  ) |>
  dplyr::select(
    patient_id = cases.submitter_id,
    sample_id = sample.submitter_id,
    aliquot_id = cases
  )

luad_expression_map <-
  luad_expression_map[
    match(
      luad_survival_metadata$submitter_id,
      luad_expression_map$patient_id
    ),
  ]

stopifnot(
  identical(
    luad_survival_metadata$submitter_id,
    luad_expression_map$patient_id
  )
)


# ------------------------------------------------------------
# 19. Locate LUAD TPM columns
# ------------------------------------------------------------

luad_tpm_columns <- match(
  luad_expression_map$aliquot_id,
  colnames(luad_tpm)
)

stopifnot(
  !anyNA(
    luad_tpm_columns
  )
)


# ------------------------------------------------------------
# 20. Locate LUAD candidate genes
# ------------------------------------------------------------

luad_candidate_rows <- match(
  survival_candidate_genes,
  luad_gene_annotation$gene_name
)

stopifnot(
  !anyNA(
    luad_candidate_rows
  )
)


# ------------------------------------------------------------
# 21. Extract and transform LUAD expression
# ------------------------------------------------------------

luad_survival_expression <-
  luad_tpm[
    luad_candidate_rows,
    luad_tpm_columns,
    drop = FALSE
  ]

rownames(
  luad_survival_expression
) <- survival_candidate_genes

colnames(
  luad_survival_expression
) <- luad_expression_map$patient_id

luad_survival_expression <-
  log2(
    luad_survival_expression + 1
  )


# ============================================================
# PART E — LUSC SURVIVAL EXPRESSION
# ============================================================


# ------------------------------------------------------------
# 22. LUSC survival metadata
# ------------------------------------------------------------

lusc_survival_metadata <- lusc_survival |>
  dplyr::filter(
    submitter_id %in%
      lusc_survival_overlap
  )


# ------------------------------------------------------------
# 23. Map LUSC patients to selected RNA aliquots
# ------------------------------------------------------------

lusc_expression_map <- lusc_tumor_one_patient |>
  dplyr::filter(
    cases.submitter_id %in%
      lusc_survival_metadata$submitter_id
  ) |>
  dplyr::select(
    patient_id = cases.submitter_id,
    sample_id = sample.submitter_id,
    aliquot_id = cases
  )

lusc_expression_map <-
  lusc_expression_map[
    match(
      lusc_survival_metadata$submitter_id,
      lusc_expression_map$patient_id
    ),
  ]

stopifnot(
  identical(
    lusc_survival_metadata$submitter_id,
    lusc_expression_map$patient_id
  )
)


# ------------------------------------------------------------
# 24. Locate LUSC TPM columns
# ------------------------------------------------------------

lusc_tpm_columns <- match(
  lusc_expression_map$aliquot_id,
  colnames(lusc_tpm)
)

stopifnot(
  !anyNA(
    lusc_tpm_columns
  )
)


# ------------------------------------------------------------
# 25. Locate LUSC candidate genes
# ------------------------------------------------------------

lusc_candidate_rows <- match(
  survival_candidate_genes,
  lusc_gene_annotation$gene_name
)

stopifnot(
  !anyNA(
    lusc_candidate_rows
  )
)


# ------------------------------------------------------------
# 26. Extract and transform LUSC expression
# ------------------------------------------------------------

lusc_survival_expression <-
  lusc_tpm[
    lusc_candidate_rows,
    lusc_tpm_columns,
    drop = FALSE
  ]

rownames(
  lusc_survival_expression
) <- survival_candidate_genes

colnames(
  lusc_survival_expression
) <- lusc_expression_map$patient_id

lusc_survival_expression <-
  log2(
    lusc_survival_expression + 1
  )


# ============================================================
# PART F — FINAL QC
# ============================================================


stopifnot(
  identical(
    dim(luad_survival_expression),
    c(47L, 504L)
  )
)

stopifnot(
  identical(
    dim(lusc_survival_expression),
    c(47L, 493L)
  )
)

stopifnot(
  identical(
    colnames(luad_survival_expression),
    luad_survival_metadata$submitter_id
  )
)

stopifnot(
  identical(
    colnames(lusc_survival_expression),
    lusc_survival_metadata$submitter_id
  )
)

stopifnot(
  !anyNA(
    luad_survival_expression
  )
)

stopifnot(
  !anyNA(
    lusc_survival_expression
  )
)

stopifnot(
  all(
    luad_survival_metadata$OS_time > 0
  )
)

stopifnot(
  all(
    lusc_survival_metadata$OS_time > 0
  )
)

stopifnot(
  all(
    luad_survival_metadata$OS_status %in%
      c(0, 1)
  )
)

stopifnot(
  all(
    lusc_survival_metadata$OS_status %in%
      c(0, 1)
  )
)


# ------------------------------------------------------------
# 26A. Final matrix/metadata alignment QC
# ------------------------------------------------------------

stopifnot(
  identical(
    colnames(luad_survival_expression),
    luad_survival_metadata$submitter_id
  ),
  identical(
    colnames(lusc_survival_expression),
    lusc_survival_metadata$submitter_id
  ),
  identical(
    dim(luad_survival_expression),
    c(47L, 504L)
  ),
  identical(
    dim(lusc_survival_expression),
    c(47L, 493L)
  ),
  !anyNA(luad_survival_expression),
  !anyNA(lusc_survival_expression),
  !anyDuplicated(luad_survival_metadata$submitter_id),
  !anyDuplicated(lusc_survival_metadata$submitter_id)
)


# ============================================================
# PART G — SAVE EXACT FINAL OUTPUTS
# ============================================================


saveRDS(
  luad_survival_expression,
  file.path(
    survival_processed_dir,
    "TCGA_LUAD_PCD_survival_expression.rds"
  )
)

write.csv(
  luad_survival_metadata,
  file.path(
    survival_processed_dir,
    "TCGA_LUAD_PCD_survival_metadata.csv"
  ),
  row.names = FALSE
)

saveRDS(
  lusc_survival_expression,
  file.path(
    survival_processed_dir,
    "TCGA_LUSC_PCD_survival_expression.rds"
  )
)

write.csv(
  lusc_survival_metadata,
  file.path(
    survival_processed_dir,
    "TCGA_LUSC_PCD_survival_metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART H — COMPLETION SUMMARY
# ============================================================

cat(
  "\n====================================================\n"
)

cat(
  "Survival-expression preparation completed.\n"
)

cat(
  "====================================================\n\n"
)

cat(
  "Candidate PCD genes: 47\n\n"
)

cat(
  "LUAD:\n",
  "  Tumor-expression patients before clinical matching: 517\n",
  "  Survival-matched patients: 504\n",
  "  Expression matrix: 47 x 504\n\n",
  sep = ""
)

cat(
  "LUSC:\n",
  "  Tumor-expression patients before clinical matching: 501\n",
  "  Survival-matched patients: 493\n",
  "  Expression matrix: 47 x 493\n\n",
  sep = ""
)

cat(
  "Expression scale: log2(TPM + 1)\n"
)

cat(
  "\n========================================\n",
  "SCRIPT 12 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================