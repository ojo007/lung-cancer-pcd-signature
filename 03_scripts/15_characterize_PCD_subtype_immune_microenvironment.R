# ============================================================
# Script 15: Characterize PCD subtype immune microenvironment
#
# Thesis:
# Exploring a Specialized Programmed Cell Death Pattern to
# Predict Prognosis and Treatment Sensitivity in Lung Cancer
# by Machine Learning and Multi-Omics Analysis
#
# Cohorts:
# TCGA-LUAD
# TCGA-LUSC
#
# Main analyses:
#   1. ESTIMATE
#      - StromalScore
#      - ImmuneScore
#      - ESTIMATEScore
#
#   2. MCP-counter
#      - 8 immune populations
#      - 2 stromal populations
#
# PCD clusters were defined previously in Script 14.
#
# IMPORTANT:
#   The full TCGA TPM matrices are sample-level matrices.
#   Therefore, this script reproduces the same one-primary-
#   tumor-per-patient selection used for PCD clustering.
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "tibble",
  "estimate",
  "MCPcounter"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

expression_dir <- file.path(
  processed_dir,
  "expression"
)

clustering_dir <- file.path(
  results_dir,
  "clustering"
)

immune_dir <- file.path(
  results_dir,
  "immune"
)

estimate_dir <- file.path(
  immune_dir,
  "estimate_temp"
)

ensure_dir(immune_dir)
ensure_dir(estimate_dir)


# ------------------------------------------------------------
# 3. Define exact processed-data files
# ------------------------------------------------------------

luad_tpm_file <- file.path(
  expression_dir,
  "TCGA_LUAD_TPM.rds"
)

lusc_tpm_file <- file.path(
  expression_dir,
  "TCGA_LUSC_TPM.rds"
)

luad_counts_file <- file.path(
  expression_dir,
  "TCGA_LUAD_raw_counts.rds"
)

lusc_counts_file <- file.path(
  expression_dir,
  "TCGA_LUSC_raw_counts.rds"
)

luad_metadata_file <- file.path(
  expression_dir,
  "TCGA_LUAD_sample_metadata.rds"
)

lusc_metadata_file <- file.path(
  expression_dir,
  "TCGA_LUSC_sample_metadata.rds"
)

luad_annotation_file <- file.path(
  expression_dir,
  "TCGA_LUAD_gene_annotation.rds"
)

lusc_annotation_file <- file.path(
  expression_dir,
  "TCGA_LUSC_gene_annotation.rds"
)

luad_cluster_file <- file.path(
  clustering_dir,
  "TCGA_LUAD_PCD_cluster_assignments.csv"
)

lusc_cluster_file <- file.path(
  clustering_dir,
  "TCGA_LUSC_PCD_cluster_assignments.csv"
)

required_files <- c(
  luad_tpm_file,
  lusc_tpm_file,
  luad_counts_file,
  lusc_counts_file,
  luad_metadata_file,
  lusc_metadata_file,
  luad_annotation_file,
  lusc_annotation_file,
  luad_cluster_file,
  lusc_cluster_file
)

if (!all(file.exists(required_files))) {
  
  stop(
    "One or more required processed expression files are missing:\n",
    paste(
      required_files[
        !file.exists(required_files)
      ],
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 4. Load expression, counts, metadata and annotations
# ------------------------------------------------------------

luad_tpm <- readRDS(
  luad_tpm_file
)

lusc_tpm <- readRDS(
  lusc_tpm_file
)

luad_counts <- readRDS(
  luad_counts_file
)

lusc_counts <- readRDS(
  lusc_counts_file
)

luad_sample_metadata <- readRDS(
  luad_metadata_file
)

lusc_sample_metadata <- readRDS(
  lusc_metadata_file
)

luad_gene_annotation <- readRDS(
  luad_annotation_file
)

lusc_gene_annotation <- readRDS(
  lusc_annotation_file
)


# ------------------------------------------------------------
# 5. Validate expression objects
# ------------------------------------------------------------

stopifnot(
  nrow(luad_tpm) == nrow(luad_counts)
)

stopifnot(
  nrow(lusc_tpm) == nrow(lusc_counts)
)

stopifnot(
  identical(
    colnames(luad_tpm),
    colnames(luad_counts)
  )
)

stopifnot(
  identical(
    colnames(lusc_tpm),
    colnames(lusc_counts)
  )
)

stopifnot(
  identical(
    colnames(luad_counts),
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
  "sample_type" %in%
    colnames(
      luad_sample_metadata
    )
)

stopifnot(
  "sample_type" %in%
    colnames(
      lusc_sample_metadata
    )
)


# ------------------------------------------------------------
# 6. Identify gene-symbol column from annotation
# ------------------------------------------------------------

get_gene_symbols <- function(
    annotation,
    expected_n
) {
  
  stopifnot(
    nrow(annotation) == expected_n
  )
  
  candidate_columns <- c(
    "gene_name",
    "gene_symbol",
    "symbol",
    "GeneSymbol",
    "external_gene_name"
  )
  
  symbol_column <- candidate_columns[
    candidate_columns %in%
      colnames(annotation)
  ]
  
  if (length(symbol_column) == 0) {
    
    stop(
      "Could not identify a gene-symbol column in annotation.\n",
      "Available columns: ",
      paste(
        colnames(annotation),
        collapse = ", "
      )
    )
  }
  
  symbols <- as.character(
    annotation[[
      symbol_column[1]
    ]]
  )
  
  symbols
}


luad_gene_symbols <- get_gene_symbols(
  luad_gene_annotation,
  nrow(luad_tpm)
)

lusc_gene_symbols <- get_gene_symbols(
  lusc_gene_annotation,
  nrow(lusc_tpm)
)


# ------------------------------------------------------------
# 7. Load PCD cluster assignments from Script 14
# ------------------------------------------------------------

luad_cluster_assignment <- read.csv(
  luad_cluster_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lusc_cluster_assignment <- read.csv(
  lusc_cluster_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster"
    ) %in%
      colnames(
        luad_cluster_assignment
      )
  )
)

stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster"
    ) %in%
      colnames(
        lusc_cluster_assignment
      )
  )
)

stopifnot(
  nrow(luad_cluster_assignment) == 517
)

stopifnot(
  nrow(lusc_cluster_assignment) == 501
)

stopifnot(
  !anyDuplicated(
    luad_cluster_assignment$patient_id
  )
)

stopifnot(
  !anyDuplicated(
    lusc_cluster_assignment$patient_id
  )
)


# ------------------------------------------------------------
# 8. Select one Primary Tumor sample per patient
#
# Same biological selection rule used for clustering:
#
#   1. Primary Tumor samples only
#   2. Patient ID = first 12 characters of TCGA barcode
#   3. If >1 sample exists for a patient, retain the sample
#      with the largest raw-count library size
#   4. Restrict to patients present in the final PCD clusters
# ------------------------------------------------------------

select_one_tumor_per_patient <- function(
    counts,
    tpm,
    metadata,
    cluster_assignment,
    gene_symbols
) {
  
  stopifnot(
    identical(
      colnames(counts),
      colnames(tpm)
    )
  )
  
  stopifnot(
    identical(
      colnames(counts),
      metadata$cases
    )
  )
  
  stopifnot(
    length(gene_symbols) ==
      nrow(tpm)
  )
  
  sample_info <- metadata |>
    dplyr::mutate(
      
      patient_id = substr(
        cases,
        1,
        12
      ),
      
      library_size =
        colSums(counts)
    ) |>
    dplyr::filter(
      sample_type ==
        "Primary Tumor"
    ) |>
    dplyr::filter(
      patient_id %in%
        cluster_assignment$patient_id
    ) |>
    dplyr::arrange(
      patient_id,
      dplyr::desc(
        library_size
      )
    ) |>
    dplyr::group_by(
      patient_id
    ) |>
    dplyr::slice_head(
      n = 1
    ) |>
    dplyr::ungroup()
  
  stopifnot(
    !anyDuplicated(
      sample_info$patient_id
    )
  )
  
  stopifnot(
    setequal(
      sample_info$patient_id,
      cluster_assignment$patient_id
    )
  )
  
  expr <- tpm[
    ,
    sample_info$cases,
    drop = FALSE
  ]
  
  rownames(expr) <- gene_symbols
  
  colnames(expr) <-
    sample_info$patient_id
  
  expr <- expr[
    ,
    cluster_assignment$patient_id,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      colnames(expr),
      cluster_assignment$patient_id
    )
  )
  
  expr
}


luad_immune_expression <-
  select_one_tumor_per_patient(
    counts = luad_counts,
    tpm = luad_tpm,
    metadata = luad_sample_metadata,
    cluster_assignment =
      luad_cluster_assignment,
    gene_symbols =
      luad_gene_symbols
  )

lusc_immune_expression <-
  select_one_tumor_per_patient(
    counts = lusc_counts,
    tpm = lusc_tpm,
    metadata = lusc_sample_metadata,
    cluster_assignment =
      lusc_cluster_assignment,
    gene_symbols =
      lusc_gene_symbols
  )


# ------------------------------------------------------------
# 9. Validate patient-level tumor expression matrices
# ------------------------------------------------------------

stopifnot(
  ncol(luad_immune_expression) == 517
)

stopifnot(
  ncol(lusc_immune_expression) == 501
)

stopifnot(
  identical(
    colnames(
      luad_immune_expression
    ),
    luad_cluster_assignment$patient_id
  )
)

stopifnot(
  identical(
    colnames(
      lusc_immune_expression
    ),
    lusc_cluster_assignment$patient_id
  )
)

cat(
  "\nPatient-level tumor expression matrices:\n"
)

cat(
  "LUAD:",
  nrow(luad_immune_expression),
  "genes x",
  ncol(luad_immune_expression),
  "patients\n"
)

cat(
  "LUSC:",
  nrow(lusc_immune_expression),
  "genes x",
  ncol(lusc_immune_expression),
  "patients\n"
)


# ------------------------------------------------------------
# 10. Collapse duplicated gene symbols
#
# For duplicated symbols, retain the row with the highest
# mean TPM within each cohort.
# ------------------------------------------------------------

collapse_duplicate_genes <- function(expr) {
  
  stopifnot(
    !is.null(
      rownames(expr)
    )
  )
  
  keep <-
    !is.na(
      rownames(expr)
    ) &
    rownames(expr) != ""
  
  expr <- expr[
    keep,
    ,
    drop = FALSE
  ]
  
  gene_mean <- rowMeans(
    expr,
    na.rm = TRUE
  )
  
  ord <- order(
    gene_mean,
    decreasing = TRUE
  )
  
  expr <- expr[
    ord,
    ,
    drop = FALSE
  ]
  
  expr <- expr[
    !duplicated(
      rownames(expr)
    ),
    ,
    drop = FALSE
  ]
  
  expr
}


luad_immune_unique <-
  collapse_duplicate_genes(
    luad_immune_expression
  )

lusc_immune_unique <-
  collapse_duplicate_genes(
    lusc_immune_expression
  )

stopifnot(
  !anyDuplicated(
    rownames(
      luad_immune_unique
    )
  )
)

stopifnot(
  !anyDuplicated(
    rownames(
      lusc_immune_unique
    )
  )
)

cat(
  "\nUnique gene-symbol matrices:\n"
)

cat(
  "LUAD:",
  nrow(luad_immune_unique),
  "genes\n"
)

cat(
  "LUSC:",
  nrow(lusc_immune_unique),
  "genes\n"
)


# ------------------------------------------------------------
# 11. log2(TPM + 1)
# ------------------------------------------------------------

luad_immune_log2 <- log2(
  luad_immune_unique + 1
)

lusc_immune_log2 <- log2(
  lusc_immune_unique + 1
)

stopifnot(
  !anyNA(
    luad_immune_log2
  )
)

stopifnot(
  !anyNA(
    lusc_immune_log2
  )
)


# ============================================================
# PART A: ESTIMATE
# ============================================================


# ------------------------------------------------------------
# 12. Load ESTIMATE common-gene dataset
# ------------------------------------------------------------

data(
  "common_genes",
  package = "estimate"
)

stopifnot(
  exists(
    "common_genes"
  )
)

stopifnot(
  nrow(common_genes) == 10412
)


# ------------------------------------------------------------
# 13. Write plain expression tables for ESTIMATE
#
# IMPORTANT:
#
# Do NOT call outputGCT() before filterCommonGenes().
#
# filterCommonGenes() itself creates the correctly formatted
# GCT file.
# ------------------------------------------------------------

luad_estimate_plain_file <- file.path(
  estimate_dir,
  "LUAD_estimate_expression.txt"
)

lusc_estimate_plain_file <- file.path(
  estimate_dir,
  "LUSC_estimate_expression.txt"
)

write.table(
  luad_immune_log2,
  file = luad_estimate_plain_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  lusc_immune_log2,
  file = lusc_estimate_plain_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


# ------------------------------------------------------------
# 14. Filter expression matrices to ESTIMATE common genes
# ------------------------------------------------------------

luad_estimate_common_file <- file.path(
  estimate_dir,
  "LUAD_estimate_common_genes.gct"
)

lusc_estimate_common_file <- file.path(
  estimate_dir,
  "LUSC_estimate_common_genes.gct"
)

estimate::filterCommonGenes(
  input.f =
    luad_estimate_plain_file,
  output.f =
    luad_estimate_common_file,
  id = "GeneSymbol"
)

estimate::filterCommonGenes(
  input.f =
    lusc_estimate_plain_file,
  output.f =
    lusc_estimate_common_file,
  id = "GeneSymbol"
)


# ------------------------------------------------------------
# 15. Calculate ESTIMATE scores
# ------------------------------------------------------------

luad_estimate_score_file <- file.path(
  estimate_dir,
  "LUAD_estimate_scores.gct"
)

lusc_estimate_score_file <- file.path(
  estimate_dir,
  "LUSC_estimate_scores.gct"
)

# ------------------------------------------------------------
# ESTIMATE compatibility fix
#
# Some current installations of the legacy `estimate` package
# expose estimateScore() without making its bundled SI_geneset
# object available in the function environment.  Build a local
# copy of estimateScore() and explicitly load the package data
# object into that copy's enclosing environment.
# ------------------------------------------------------------

estimate_score_env <- new.env(
  parent = environment(
    estimate::estimateScore
  )
)

utils::data(
  "SI_geneset",
  package = "estimate",
  envir = estimate_score_env
)

if (!exists(
  "SI_geneset",
  envir = estimate_score_env,
  inherits = FALSE
)) {
  stop(
    "The estimate package is installed, but its bundled ",
    "SI_geneset data object could not be loaded."
  )
}

estimateScore_local <- estimate::estimateScore
environment(estimateScore_local) <- estimate_score_env

estimateScore_local(
  luad_estimate_common_file,
  luad_estimate_score_file,
  platform = "illumina"
)

estimateScore_local(
  lusc_estimate_common_file,
  lusc_estimate_score_file,
  platform = "illumina"
)


# ------------------------------------------------------------
# 16. Read ESTIMATE score files
# ------------------------------------------------------------

read_estimate_scores <- function(
    score_file
) {
  
  x <- read.delim(
    score_file,
    skip = 2,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  stopifnot(
    nrow(x) == 3
  )
  
  stopifnot(
    all(
      c(
        "StromalScore",
        "ImmuneScore",
        "ESTIMATEScore"
      ) %in%
        x$NAME
    )
  )
  
  rownames(x) <- x$NAME
  
  x <- x[
    ,
    -c(1, 2),
    drop = FALSE
  ]
  
  x <- as.data.frame(
    t(x),
    check.names = FALSE
  )
  
  x[] <- lapply(
    x,
    as.numeric
  )
  
  x <- x |>
    tibble::rownames_to_column(
      "patient_id"
    )
  
  # ESTIMATE changes hyphens in TCGA IDs to periods.
  x$patient_id <- gsub(
    "\\.",
    "-",
    x$patient_id
  )
  
  x
}


luad_estimate_scores <-
  read_estimate_scores(
    luad_estimate_score_file
  )

lusc_estimate_scores <-
  read_estimate_scores(
    lusc_estimate_score_file
  )

stopifnot(
  nrow(luad_estimate_scores) == 517
)

stopifnot(
  nrow(lusc_estimate_scores) == 501
)


# ------------------------------------------------------------
# 17. Join ESTIMATE scores to PCD clusters
# ------------------------------------------------------------

luad_estimate_scores <-
  luad_estimate_scores |>
  dplyr::left_join(
    luad_cluster_assignment,
    by = "patient_id"
  )

lusc_estimate_scores <-
  lusc_estimate_scores |>
  dplyr::left_join(
    lusc_cluster_assignment,
    by = "patient_id"
  )

stopifnot(
  sum(
    is.na(
      luad_estimate_scores$PCD_cluster
    )
  ) == 0
)

stopifnot(
  sum(
    is.na(
      lusc_estimate_scores$PCD_cluster
    )
  ) == 0
)


# ------------------------------------------------------------
# 18. ESTIMATE cluster comparisons
# ------------------------------------------------------------

estimate_cluster_test <- function(
    score_table
) {
  
  score_table |>
    tidyr::pivot_longer(
      cols = c(
        StromalScore,
        ImmuneScore,
        ESTIMATEScore
      ),
      names_to =
        "estimate_metric",
      values_to =
        "score"
    ) |>
    dplyr::group_by(
      estimate_metric
    ) |>
    dplyr::summarise(
      
      C1_mean = mean(
        score[
          PCD_cluster ==
            "PCD_C1"
        ],
        na.rm = TRUE
      ),
      
      C2_mean = mean(
        score[
          PCD_cluster ==
            "PCD_C2"
        ],
        na.rm = TRUE
      ),
      
      C2_minus_C1 =
        C2_mean -
        C1_mean,
      
      pvalue =
        wilcox.test(
          score ~
            PCD_cluster
        )$p.value,
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      padj = p.adjust(
        pvalue,
        method = "BH"
      )
    ) |>
    dplyr::arrange(
      padj
    )
}


luad_estimate_tests <-
  estimate_cluster_test(
    luad_estimate_scores
  )

lusc_estimate_tests <-
  estimate_cluster_test(
    lusc_estimate_scores
  )


# ============================================================
# PART B: MCP-COUNTER
# ============================================================


# ------------------------------------------------------------
# 19. Calculate MCP-counter abundance scores
# ------------------------------------------------------------

luad_mcp <-
  MCPcounter::MCPcounter.estimate(
    luad_immune_log2,
    featuresType =
      "HUGO_symbols"
  )

lusc_mcp <-
  MCPcounter::MCPcounter.estimate(
    lusc_immune_log2,
    featuresType =
      "HUGO_symbols"
  )

stopifnot(
  nrow(luad_mcp) == 10,
  ncol(luad_mcp) == 517
)

stopifnot(
  nrow(lusc_mcp) == 10,
  ncol(lusc_mcp) == 501
)


# ------------------------------------------------------------
# 20. Convert MCP-counter output to patient tables
# ------------------------------------------------------------

luad_mcp_scores <-
  as.data.frame(
    t(luad_mcp),
    check.names = FALSE
  ) |>
  tibble::rownames_to_column(
    "patient_id"
  )

lusc_mcp_scores <-
  as.data.frame(
    t(lusc_mcp),
    check.names = FALSE
  ) |>
  tibble::rownames_to_column(
    "patient_id"
  )


# ------------------------------------------------------------
# 21. Join MCP-counter scores to PCD clusters
# ------------------------------------------------------------

luad_mcp_scores <-
  luad_mcp_scores |>
  dplyr::left_join(
    luad_cluster_assignment,
    by = "patient_id"
  )

lusc_mcp_scores <-
  lusc_mcp_scores |>
  dplyr::left_join(
    lusc_cluster_assignment,
    by = "patient_id"
  )

stopifnot(
  sum(
    is.na(
      luad_mcp_scores$PCD_cluster
    )
  ) == 0
)

stopifnot(
  sum(
    is.na(
      lusc_mcp_scores$PCD_cluster
    )
  ) == 0
)


# ------------------------------------------------------------
# 22. MCP-counter cluster comparisons
# ------------------------------------------------------------

mcp_cluster_test <- function(
    score_table
) {
  
  score_table |>
    tidyr::pivot_longer(
      cols = -c(
        patient_id,
        PCD_cluster
      ),
      names_to =
        "cell_population",
      values_to =
        "score"
    ) |>
    dplyr::group_by(
      cell_population
    ) |>
    dplyr::summarise(
      
      C1_mean = mean(
        score[
          PCD_cluster ==
            "PCD_C1"
        ],
        na.rm = TRUE
      ),
      
      C2_mean = mean(
        score[
          PCD_cluster ==
            "PCD_C2"
        ],
        na.rm = TRUE
      ),
      
      C2_minus_C1 =
        C2_mean -
        C1_mean,
      
      pvalue =
        wilcox.test(
          score ~
            PCD_cluster
        )$p.value,
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      padj = p.adjust(
        pvalue,
        method = "BH"
      )
    ) |>
    dplyr::arrange(
      padj
    )
}


luad_mcp_tests <-
  mcp_cluster_test(
    luad_mcp_scores
  )

lusc_mcp_tests <-
  mcp_cluster_test(
    lusc_mcp_scores
  )


# ------------------------------------------------------------
# 23. Cross-cohort MCP-counter comparison
# ------------------------------------------------------------

mcp_cross_cohort <-
  luad_mcp_tests |>
  dplyr::select(
    cell_population,
    LUAD_C1_mean =
      C1_mean,
    LUAD_C2_mean =
      C2_mean,
    LUAD_difference =
      C2_minus_C1,
    LUAD_padj =
      padj
  ) |>
  dplyr::left_join(
    
    lusc_mcp_tests |>
      dplyr::select(
        cell_population,
        LUSC_C1_mean =
          C1_mean,
        LUSC_C2_mean =
          C2_mean,
        LUSC_difference =
          C2_minus_C1,
        LUSC_padj =
          padj
      ),
    
    by = "cell_population"
  ) |>
  dplyr::mutate(
    
    LUAD_significant =
      LUAD_padj < 0.05,
    
    LUSC_significant =
      LUSC_padj < 0.05,
    
    direction_concordant =
      sign(
        LUAD_difference
      ) ==
      sign(
        LUSC_difference
      ),
    
    significant_both =
      LUAD_significant &
      LUSC_significant,
    
    robust_cross_cohort =
      significant_both &
      direction_concordant
  )


# ------------------------------------------------------------
# 24. Save final ESTIMATE results
# ------------------------------------------------------------

write.csv(
  luad_estimate_scores,
  file.path(
    immune_dir,
    "TCGA_LUAD_ESTIMATE_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_estimate_scores,
  file.path(
    immune_dir,
    "TCGA_LUSC_ESTIMATE_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  luad_estimate_tests,
  file.path(
    immune_dir,
    "TCGA_LUAD_ESTIMATE_cluster_tests.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_estimate_tests,
  file.path(
    immune_dir,
    "TCGA_LUSC_ESTIMATE_cluster_tests.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. Save final MCP-counter results
# ------------------------------------------------------------

write.csv(
  luad_mcp_scores,
  file.path(
    immune_dir,
    "TCGA_LUAD_MCPcounter_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_mcp_scores,
  file.path(
    immune_dir,
    "TCGA_LUSC_MCPcounter_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  luad_mcp_tests,
  file.path(
    immune_dir,
    "TCGA_LUAD_MCPcounter_cluster_tests.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_mcp_tests,
  file.path(
    immune_dir,
    "TCGA_LUSC_MCPcounter_cluster_tests.csv"
  ),
  row.names = FALSE
)

write.csv(
  mcp_cross_cohort,
  file.path(
    immune_dir,
    "TCGA_LUAD_LUSC_MCPcounter_cross_cohort.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Save cross-cohort robust MCP-counter populations
# ------------------------------------------------------------

mcp_robust_cross_cohort <-
  mcp_cross_cohort |>
  dplyr::filter(
    robust_cross_cohort
  )

write.csv(
  mcp_robust_cross_cohort,
  file.path(
    immune_dir,
    "TCGA_LUAD_LUSC_MCPcounter_robust_populations.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Final validation summary
# ------------------------------------------------------------

cat(
  "\n==========================================================\n"
)

cat(
  "IMMUNE ANALYSIS COMPLETE\n"
)

cat(
  "==========================================================\n"
)

cat(
  "\nLUAD cluster counts:\n"
)

print(
  table(
    luad_estimate_scores$PCD_cluster
  )
)

cat(
  "\nLUSC cluster counts:\n"
)

print(
  table(
    lusc_estimate_scores$PCD_cluster
  )
)

cat(
  "\nLUAD ESTIMATE tests:\n"
)

print(
  luad_estimate_tests
)

cat(
  "\nLUSC ESTIMATE tests:\n"
)

print(
  lusc_estimate_tests
)

cat(
  "\nLUAD MCP-counter tests:\n"
)

print(
  luad_mcp_tests
)

cat(
  "\nLUSC MCP-counter tests:\n"
)

print(
  lusc_mcp_tests
)

cat(
  "\nCross-cohort MCP-counter comparison:\n"
)

print(
  mcp_cross_cohort
)

cat(
  "\nRobust cross-cohort MCP-counter populations:\n"
)

print(
  mcp_robust_cross_cohort
)

cat(
  "\nFinal files saved in:\n"
)

cat(
  immune_dir,
  "\n"
)

cat(
  "\n==========================================================\n"
)

cat(
  "\n========================================\n",
  "SCRIPT 15 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "Outputs: ",
  immune_dir,
  "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================