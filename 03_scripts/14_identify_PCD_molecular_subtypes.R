# ============================================================
# 14_identify_PCD_molecular_subtypes.R
#
# Purpose:
# Identify and characterize PCD molecular subtypes separately
# in TCGA-LUAD and TCGA-LUSC.
#
# Clustering input:
# 47 concordant LUAD/LUSC differentially expressed PCD genes
#
# Expression:
# Tumor-only log2(TPM + 1), one sample per patient
#
# Consensus clustering:
# ConsensusClusterPlus
# k = 2:6
# reps = 1000
# pItem = 0.8
# pFeature = 1
# hierarchical clustering
# Pearson distance
# average linkage
#
# Final selected solution:
# LUAD: k = 2
# LUSC: k = 2
#
# Final cluster labels remain neutral:
# PCD_C1
# PCD_C2
#
# Biological characterization:
# Full curated 296-gene PCD catalogue
# - Apoptosis
# - Necroptosis
# - Pyroptosis
# - Ferroptosis
# - Cuproptosis
#
# Final CSV outputs:
#
# 04_results/clustering/
# 1. TCGA_LUAD_PCD_cluster_assignments.csv
# 2. TCGA_LUSC_PCD_cluster_assignments.csv
# 3. TCGA_LUAD_PCD_mechanism_scores.csv
# 4. TCGA_LUSC_PCD_mechanism_scores.csv
# 5. TCGA_LUAD_PCD_mechanism_cluster_tests.csv
# 6. TCGA_LUSC_PCD_mechanism_cluster_tests.csv
#
# Consensus diagnostic PNGs are produced automatically under:
# 04_results/clustering/LUAD_consensus/
# 04_results/clustering/LUSC_consensus/
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

expression_dir <- file.path(
  processed_dir,
  "expression"
)

pcd_processed_dir <- file.path(
  processed_dir,
  "PCD_genes"
)

de_results_dir <- file.path(
  results_dir,
  "differential_expression"
)

survival_processed_dir <- file.path(
  processed_dir,
  "survival"
)

clustering_results_dir <- file.path(
  results_dir,
  "clustering"
)

luad_ccp_dir <- file.path(
  clustering_results_dir,
  "LUAD_consensus"
)

lusc_ccp_dir <- file.path(
  clustering_results_dir,
  "LUSC_consensus"
)

clustering_figure_dir <- file.path(
  individual_figure_dir,
  "clustering"
)

ensure_dir(clustering_results_dir)
ensure_dir(luad_ccp_dir)
ensure_dir(lusc_ccp_dir)
ensure_dir(clustering_figure_dir)


# ------------------------------------------------------------
# 2. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "survival",
  "cluster",
  "ConsensusClusterPlus",
  "IntNMF"
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
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(cluster)
  library(ConsensusClusterPlus)
  library(IntNMF)
})

# The molecular subtype definition remains expression-based.
# CPI below is therefore implemented as a prediction-based,
# resampling stability index on the same 47-gene expression
# space. MOVICS::getClustNum() is NOT used because MOVICS
# explicitly requires at least two omics data layers.


# ------------------------------------------------------------
# 2A. Clustering parameters
# ------------------------------------------------------------

K_RANGE <- 2:8
SELECTED_K <- 2L
CONSENSUS_REPS <- 1000L
CPI_RUNS <- 30L
CPI_FOLDS <- 5L
GAP_BOOTSTRAPS <- 100L
RANDOM_SEED <- 12345L


# ============================================================
# PART A — LOAD DATA
# ============================================================


# ------------------------------------------------------------
# 3. RNA-seq matrices
# ------------------------------------------------------------

luad_counts <- readRDS(
  file.path(
    expression_dir,
    "TCGA_LUAD_raw_counts.rds"
  )
)

luad_tpm <- readRDS(
  file.path(
    expression_dir,
    "TCGA_LUAD_TPM.rds"
  )
)

lusc_counts <- readRDS(
  file.path(
    expression_dir,
    "TCGA_LUSC_raw_counts.rds"
  )
)

lusc_tpm <- readRDS(
  file.path(
    expression_dir,
    "TCGA_LUSC_TPM.rds"
  )
)


# ------------------------------------------------------------
# 4. Gene annotations
# ------------------------------------------------------------

luad_gene_annotation <- readRDS(
  file.path(
    expression_dir,
    "TCGA_LUAD_gene_annotation.rds"
  )
)

lusc_gene_annotation <- readRDS(
  file.path(
    expression_dir,
    "TCGA_LUSC_gene_annotation.rds"
  )
)


# ------------------------------------------------------------
# 5. Sample metadata
# ------------------------------------------------------------

luad_sample_metadata <- readRDS(
  file.path(
    expression_dir,
    "TCGA_LUAD_sample_metadata.rds"
  )
)

lusc_sample_metadata <- readRDS(
  file.path(
    expression_dir,
    "TCGA_LUSC_sample_metadata.rds"
  )
)


# ------------------------------------------------------------
# 6. Concordant 47-gene PCD candidate set
# ------------------------------------------------------------

pcd_candidates <- read.csv(
  file.path(
    de_results_dir,
    "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv"
  ),
  stringsAsFactors = FALSE
)

cluster_genes <-
  pcd_candidates$gene_name

stopifnot(
  length(cluster_genes) == 47
)

stopifnot(
  !anyDuplicated(cluster_genes)
)


# ------------------------------------------------------------
# 7. Full PCD gene catalogue
# ------------------------------------------------------------

pcd_master_gene <- read.csv(
  file.path(
    pcd_processed_dir,
    "PCD_master_unique_genes.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(pcd_master_gene) == 296
)


# ------------------------------------------------------------
# 8. Expand multi-PCD memberships
# ------------------------------------------------------------

pcd_full_membership <- pcd_master_gene |>
  dplyr::select(
    gene_symbol,
    PCD_types
  ) |>
  tidyr::separate_rows(
    PCD_types,
    sep = ";\\s*"
  ) |>
  dplyr::rename(
    PCD_type = PCD_types
  ) |>
  dplyr::distinct(
    gene_symbol,
    PCD_type
  )


# ------------------------------------------------------------
# 9. Validate pathway sizes
# ------------------------------------------------------------

pcd_type_counts <- table(
  pcd_full_membership$PCD_type
)

stopifnot(
  pcd_type_counts["Apoptosis"] == 166
)

stopifnot(
  pcd_type_counts["Cuproptosis"] == 12
)

stopifnot(
  pcd_type_counts["Ferroptosis"] == 89
)

stopifnot(
  pcd_type_counts["Necroptosis"] == 32
)

stopifnot(
  pcd_type_counts["Pyroptosis"] == 27
)


# ============================================================
# PART B — SELECT ONE PRIMARY TUMOR PROFILE PER PATIENT
# ============================================================


# ------------------------------------------------------------
# 10. Validate matrix/metadata alignment
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
# 11. Raw-count library sizes
# ------------------------------------------------------------

luad_sample_metadata$library_size <-
  colSums(luad_counts)

lusc_sample_metadata$library_size <-
  colSums(lusc_counts)


# ------------------------------------------------------------
# 12. LUAD one aliquot/sample/patient
# ------------------------------------------------------------

luad_tumor_one_patient <- luad_sample_metadata |>
  dplyr::filter(
    sample_type == "Primary Tumor"
  ) |>
  dplyr::group_by(
    sample.submitter_id
  ) |>
  dplyr::slice_max(
    library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(
    cases.submitter_id
  ) |>
  dplyr::slice_max(
    library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()


# ------------------------------------------------------------
# 13. LUSC one aliquot/sample/patient
# ------------------------------------------------------------

lusc_tumor_one_patient <- lusc_sample_metadata |>
  dplyr::filter(
    sample_type == "Primary Tumor"
  ) |>
  dplyr::group_by(
    sample.submitter_id
  ) |>
  dplyr::slice_max(
    library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(
    cases.submitter_id
  ) |>
  dplyr::slice_max(
    library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()


# ------------------------------------------------------------
# 14. Validate patient counts
# ------------------------------------------------------------

stopifnot(
  nrow(luad_tumor_one_patient) == 517
)

stopifnot(
  nrow(lusc_tumor_one_patient) == 501
)

stopifnot(
  !anyDuplicated(
    luad_tumor_one_patient$cases.submitter_id
  )
)

stopifnot(
  !anyDuplicated(
    lusc_tumor_one_patient$cases.submitter_id
  )
)


# ============================================================
# PART C — PREPARE 47-GENE CLUSTERING MATRICES
# ============================================================


# ------------------------------------------------------------
# 15. Match TPM columns
# ------------------------------------------------------------

luad_cluster_tpm_cols <- match(
  luad_tumor_one_patient$cases,
  colnames(luad_tpm)
)

lusc_cluster_tpm_cols <- match(
  lusc_tumor_one_patient$cases,
  colnames(lusc_tpm)
)

stopifnot(
  !anyNA(luad_cluster_tpm_cols)
)

stopifnot(
  !anyNA(lusc_cluster_tpm_cols)
)


# ------------------------------------------------------------
# 16. Match candidate genes
# ------------------------------------------------------------

luad_cluster_gene_rows <- match(
  cluster_genes,
  luad_gene_annotation$gene_name
)

lusc_cluster_gene_rows <- match(
  cluster_genes,
  lusc_gene_annotation$gene_name
)

stopifnot(
  !anyNA(luad_cluster_gene_rows)
)

stopifnot(
  !anyNA(lusc_cluster_gene_rows)
)


# ------------------------------------------------------------
# 17. LUAD 47-gene expression
# ------------------------------------------------------------

luad_cluster_expression <- luad_tpm[
  luad_cluster_gene_rows,
  luad_cluster_tpm_cols,
  drop = FALSE
]

rownames(
  luad_cluster_expression
) <- cluster_genes

colnames(
  luad_cluster_expression
) <- luad_tumor_one_patient$cases.submitter_id

luad_cluster_expression <-
  log2(
    luad_cluster_expression + 1
  )


# ------------------------------------------------------------
# 18. LUSC 47-gene expression
# ------------------------------------------------------------

lusc_cluster_expression <- lusc_tpm[
  lusc_cluster_gene_rows,
  lusc_cluster_tpm_cols,
  drop = FALSE
]

rownames(
  lusc_cluster_expression
) <- cluster_genes

colnames(
  lusc_cluster_expression
) <- lusc_tumor_one_patient$cases.submitter_id

lusc_cluster_expression <-
  log2(
    lusc_cluster_expression + 1
  )


# ------------------------------------------------------------
# 19. QC
# ------------------------------------------------------------

stopifnot(
  identical(
    dim(luad_cluster_expression),
    c(47L, 517L)
  )
)

stopifnot(
  identical(
    dim(lusc_cluster_expression),
    c(47L, 501L)
  )
)

stopifnot(
  !anyNA(luad_cluster_expression)
)

stopifnot(
  !anyNA(lusc_cluster_expression)
)

stopifnot(
  all(
    apply(
      luad_cluster_expression,
      1,
      sd
    ) > 0
  )
)

stopifnot(
  all(
    apply(
      lusc_cluster_expression,
      1,
      sd
    ) > 0
  )
)


# ============================================================
# PART D — STANDARDIZE EXPRESSION
# ============================================================


luad_cluster_z <- t(
  scale(
    t(luad_cluster_expression)
  )
)

lusc_cluster_z <- t(
  scale(
    t(lusc_cluster_expression)
  )
)

stopifnot(
  !anyNA(luad_cluster_z)
)

stopifnot(
  !anyNA(lusc_cluster_z)
)


# ============================================================
# PART E — CONSENSUS CLUSTERING
# ============================================================


# ------------------------------------------------------------
# 20. LUAD
# ------------------------------------------------------------

set.seed(RANDOM_SEED)

luad_ccp <- ConsensusClusterPlus(
  as.matrix(luad_cluster_z),
  maxK = max(K_RANGE),
  reps = CONSENSUS_REPS,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "pearson",
  innerLinkage = "average",
  finalLinkage = "average",
  seed = RANDOM_SEED,
  plot = "png",
  title = luad_ccp_dir,
  verbose = TRUE
)


# ------------------------------------------------------------
# 21. LUSC
# ------------------------------------------------------------

set.seed(RANDOM_SEED)

lusc_ccp <- ConsensusClusterPlus(
  as.matrix(lusc_cluster_z),
  maxK = max(K_RANGE),
  reps = CONSENSUS_REPS,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "pearson",
  innerLinkage = "average",
  finalLinkage = "average",
  seed = RANDOM_SEED,
  plot = "png",
  title = lusc_ccp_dir,
  verbose = TRUE
)


# ============================================================
# PART F — CLUSTER-NUMBER VALIDATION
#
# Supervisor-requested validation:
# 1. Consensus/PAC across k
# 2. Cluster Prediction Index (CPI) using IntNMF::nmf.opt.k()
# 3. Standard Gap statistic using cluster::clusGap()
# 4. Mean silhouette across k
# 5. Individual-sample silhouette plot for selected k = 2
#
# IMPORTANT:
# The molecular subtypes themselves remain defined by
# ConsensusClusterPlus on the 47-gene expression matrix.
# CPI and Gap are independent cluster-number diagnostics.
# ============================================================


# ------------------------------------------------------------
# 22. PAC from ConsensusClusterPlus
# ------------------------------------------------------------

calculate_pac <- function(
    ccp_result,
    k,
    lower = 0.1,
    upper = 0.9
) {
  
  cm <- ccp_result[[k]]$consensusMatrix
  
  values <- cm[
    upper.tri(cm)
  ]
  
  mean(
    values > lower &
      values < upper
  )
}

luad_pac <- data.frame(
  k = K_RANGE,
  PAC = vapply(
    K_RANGE,
    function(k) {
      calculate_pac(
        luad_ccp,
        k
      )
    },
    numeric(1)
  )
)

lusc_pac <- data.frame(
  k = K_RANGE,
  PAC = vapply(
    K_RANGE,
    function(k) {
      calculate_pac(
        lusc_ccp,
        k
      )
    },
    numeric(1)
  )
)


# ------------------------------------------------------------
# 23. Pearson-distance matrices for silhouette analysis
# ------------------------------------------------------------

make_pearson_distance <- function(
    feature_by_sample_matrix
) {
  
  cor_mat <- stats::cor(
    feature_by_sample_matrix,
    method = "pearson"
  )
  
  cor_mat[is.na(cor_mat)] <- 0
  diag(cor_mat) <- 1
  
  cor_mat[cor_mat > 1] <- 1
  cor_mat[cor_mat < -1] <- -1
  
  stats::as.dist(
    1 - cor_mat
  )
}

luad_distance <- make_pearson_distance(
  luad_cluster_z
)

lusc_distance <- make_pearson_distance(
  lusc_cluster_z
)

stopifnot(
  attr(luad_distance, "Size") == ncol(luad_cluster_z),
  attr(lusc_distance, "Size") == ncol(lusc_cluster_z)
)


# ------------------------------------------------------------
# 24. Mean silhouette for every candidate k = 2:8
# ------------------------------------------------------------

calculate_mean_silhouette <- function(
    ccp_result,
    distance_matrix,
    k
) {
  
  sil <- cluster::silhouette(
    ccp_result[[k]]$consensusClass,
    distance_matrix
  )
  
  mean(
    sil[, "sil_width"]
  )
}

luad_silhouette <- data.frame(
  k = K_RANGE,
  mean_silhouette = vapply(
    K_RANGE,
    function(k) {
      calculate_mean_silhouette(
        luad_ccp,
        luad_distance,
        k
      )
    },
    numeric(1)
  )
)

lusc_silhouette <- data.frame(
  k = K_RANGE,
  mean_silhouette = vapply(
    K_RANGE,
    function(k) {
      calculate_mean_silhouette(
        lusc_ccp,
        lusc_distance,
        k
      )
    },
    numeric(1)
  )
)


# ------------------------------------------------------------
# 25. Prepare non-negative input for IntNMF CPI
#
# IntNMF requires samples on rows and features on columns and
# non-negative values. The 47-gene z-score matrix is therefore
# shifted/rescaled to [0, 1] without changing sample ordering.
# ------------------------------------------------------------

prepare_intnmf_input <- function(
    feature_by_sample_matrix
) {
  
  x <- t(
    feature_by_sample_matrix
  )
  
  x <- x - min(
    x,
    na.rm = TRUE
  )
  
  max_x <- max(
    x,
    na.rm = TRUE
  )
  
  if (!is.finite(max_x) || max_x <= 0) {
    stop(
      "Unable to create non-negative IntNMF input."
    )
  }
  
  x <- x / max_x
  
  x
}

luad_intnmf_input <- prepare_intnmf_input(
  luad_cluster_z
)

lusc_intnmf_input <- prepare_intnmf_input(
  lusc_cluster_z
)

stopifnot(
  nrow(luad_intnmf_input) == ncol(luad_cluster_z),
  nrow(lusc_intnmf_input) == ncol(lusc_cluster_z),
  all(luad_intnmf_input >= 0),
  all(lusc_intnmf_input >= 0)
)


# ------------------------------------------------------------
# 26. Established CPI using IntNMF::nmf.opt.k
#
# IntNMF documentation defines nmf.opt.k() for either a single
# dataset or multiple datasets. We use the single expression
# dataset here because subtype discovery is expression-based.
# ------------------------------------------------------------

set.seed(RANDOM_SEED)

luad_cpi_matrix <- IntNMF::nmf.opt.k(
  dat = luad_intnmf_input,
  n.runs = CPI_RUNS,
  n.fold = CPI_FOLDS,
  k.range = K_RANGE,
  result = TRUE,
  make.plot = FALSE,
  progress = TRUE,
  st.count = 10,
  maxiter = 100
)

set.seed(RANDOM_SEED)

lusc_cpi_matrix <- IntNMF::nmf.opt.k(
  dat = lusc_intnmf_input,
  n.runs = CPI_RUNS,
  n.fold = CPI_FOLDS,
  k.range = K_RANGE,
  result = TRUE,
  make.plot = FALSE,
  progress = TRUE,
  st.count = 10,
  maxiter = 100
)

if (
  nrow(luad_cpi_matrix) != length(K_RANGE) ||
  nrow(lusc_cpi_matrix) != length(K_RANGE)
) {
  stop(
    "Unexpected IntNMF CPI output dimensions."
  )
}

luad_cpi <- data.frame(
  k = K_RANGE,
  CPI = rowMeans(
    luad_cpi_matrix,
    na.rm = TRUE
  ),
  CPI_SD = apply(
    luad_cpi_matrix,
    1,
    stats::sd,
    na.rm = TRUE
  )
)

lusc_cpi <- data.frame(
  k = K_RANGE,
  CPI = rowMeans(
    lusc_cpi_matrix,
    na.rm = TRUE
  ),
  CPI_SD = apply(
    lusc_cpi_matrix,
    1,
    stats::sd,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# 27. Standard Gap statistic using cluster::clusGap
#
# The clustering function uses the same Pearson-distance /
# average-linkage rule used in the molecular subtype analysis.
# clusGap supplies the standard reference-distribution procedure.
# ------------------------------------------------------------

hc_pearson_gap <- function(
    x,
    k
) {
  
  if (k == 1) {
    return(
      list(
        cluster = rep(
          1L,
          nrow(x)
        )
      )
    )
  }
  
  sample_cor <- stats::cor(
    t(x),
    method = "pearson"
  )
  
  sample_cor[is.na(sample_cor)] <- 0
  diag(sample_cor) <- 1
  
  sample_cor[sample_cor > 1] <- 1
  sample_cor[sample_cor < -1] <- -1
  
  hc <- stats::hclust(
    stats::as.dist(
      1 - sample_cor
    ),
    method = "average"
  )
  
  list(
    cluster = stats::cutree(
      hc,
      k = k
    )
  )
}

set.seed(RANDOM_SEED)

luad_gap_object <- cluster::clusGap(
  x = t(luad_cluster_z),
  FUNcluster = hc_pearson_gap,
  K.max = max(K_RANGE),
  B = GAP_BOOTSTRAPS,
  d.power = 2,
  spaceH0 = "scaledPCA",
  verbose = TRUE
)

set.seed(RANDOM_SEED)

lusc_gap_object <- cluster::clusGap(
  x = t(lusc_cluster_z),
  FUNcluster = hc_pearson_gap,
  K.max = max(K_RANGE),
  B = GAP_BOOTSTRAPS,
  d.power = 2,
  spaceH0 = "scaledPCA",
  verbose = TRUE
)

luad_gap <- data.frame(
  k = seq_len(
    nrow(luad_gap_object$Tab)
  ),
  Gap = luad_gap_object$Tab[, "gap"],
  Gap_SE = luad_gap_object$Tab[, "SE.sim"]
)

lusc_gap <- data.frame(
  k = seq_len(
    nrow(lusc_gap_object$Tab)
  ),
  Gap = lusc_gap_object$Tab[, "gap"],
  Gap_SE = lusc_gap_object$Tab[, "SE.sim"]
)

# Keep k = 2:8 in the manuscript diagnostic table.
luad_gap_plot <- luad_gap[
  luad_gap$k %in% K_RANGE,
  ,
  drop = FALSE
]

lusc_gap_plot <- lusc_gap[
  lusc_gap$k %in% K_RANGE,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 28. Gap-statistic selected k
#
# We record both:
# - the absolute global maximum, and
# - Tibshirani's 1-SE rule, which is intended to identify the
#   earliest adequate/inflection solution rather than simply
#   chasing a later small increase in Gap.
# ------------------------------------------------------------

luad_gap_global_k <- which.max(
  luad_gap$Gap
)

lusc_gap_global_k <- which.max(
  lusc_gap$Gap
)

luad_gap_tibshirani_k <- cluster::maxSE(
  f = luad_gap$Gap,
  SE.f = luad_gap$Gap_SE,
  method = "Tibs2001SEmax",
  SE.factor = 1
)

lusc_gap_tibshirani_k <- cluster::maxSE(
  f = lusc_gap$Gap,
  SE.f = lusc_gap$Gap_SE,
  method = "Tibs2001SEmax",
  SE.factor = 1
)


# ------------------------------------------------------------
# 29. Summarize all cluster-number diagnostics
# ------------------------------------------------------------

luad_cluster_number_diagnostics <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "k",
      all = TRUE,
      sort = TRUE
    )
  },
  list(
    luad_pac,
    luad_silhouette,
    luad_cpi,
    luad_gap_plot
  )
)

lusc_cluster_number_diagnostics <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "k",
      all = TRUE,
      sort = TRUE
    )
  },
  list(
    lusc_pac,
    lusc_silhouette,
    lusc_cpi,
    lusc_gap_plot
  )
)

write.csv(
  luad_cluster_number_diagnostics,
  file.path(
    clustering_results_dir,
    "TCGA_LUAD_cluster_number_diagnostics_k2_k8.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_cluster_number_diagnostics,
  file.path(
    clustering_results_dir,
    "TCGA_LUSC_cluster_number_diagnostics_k2_k8.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Diagnostic optima
# ------------------------------------------------------------

luad_best_pac_k <- luad_pac$k[
  which.min(
    luad_pac$PAC
  )
]

lusc_best_pac_k <- lusc_pac$k[
  which.min(
    lusc_pac$PAC
  )
]

luad_best_silhouette_k <- luad_silhouette$k[
  which.max(
    luad_silhouette$mean_silhouette
  )
]

lusc_best_silhouette_k <- lusc_silhouette$k[
  which.max(
    lusc_silhouette$mean_silhouette
  )
]

luad_best_cpi_k <- luad_cpi$k[
  which.max(
    luad_cpi$CPI
  )
]

lusc_best_cpi_k <- lusc_cpi$k[
  which.max(
    lusc_cpi$CPI
  )
]


# ------------------------------------------------------------
# 31. Individual-sample silhouette objects for selected k = 2
# ------------------------------------------------------------

luad_selected_silhouette <- cluster::silhouette(
  luad_ccp[[SELECTED_K]]$consensusClass,
  luad_distance
)

lusc_selected_silhouette <- cluster::silhouette(
  lusc_ccp[[SELECTED_K]]$consensusClass,
  lusc_distance
)


# ------------------------------------------------------------
# 32. Save supervisor-requested diagnostic figures
# ------------------------------------------------------------

save_three_formats <- function(
    filename_stem,
    plot_function,
    width = 7,
    height = 5
) {
  
  png(
    file.path(
      clustering_figure_dir,
      paste0(
        filename_stem,
        ".png"
      )
    ),
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  plot_function()
  dev.off()
  
  pdf(
    file.path(
      clustering_figure_dir,
      paste0(
        filename_stem,
        ".pdf"
      )
    ),
    width = width,
    height = height
  )
  plot_function()
  dev.off()
  
  tiff(
    file.path(
      clustering_figure_dir,
      paste0(
        filename_stem,
        ".tiff"
      )
    ),
    width = width,
    height = height,
    units = "in",
    res = 300,
    compression = "lzw"
  )
  plot_function()
  dev.off()
}

plot_cpi_gap <- function(
    cpi_table,
    gap_table,
    cohort_label
) {
  
  old_par <- par(
    no.readonly = TRUE
  )
  
  on.exit(
    par(old_par)
  )
  
  plot(
    cpi_table$k,
    cpi_table$CPI,
    type = "b",
    pch = 19,
    lty = 2,
    col = "steelblue4",
    xlab = "Number of clusters (k)",
    ylab = "Cluster Prediction Index (CPI)",
    main = paste(
      cohort_label,
      "CPI and Gap statistic"
    ),
    xaxt = "n"
  )
  
  axis(
    1,
    at = K_RANGE
  )
  
  par(
    new = TRUE
  )
  
  plot(
    gap_table$k,
    gap_table$Gap,
    type = "b",
    pch = 17,
    lty = 2,
    col = "firebrick3",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  
  axis(
    4,
    col.axis = "firebrick3"
  )
  
  mtext(
    "Gap statistic",
    side = 4,
    line = 3,
    col = "firebrick3"
  )
  
  abline(
    v = SELECTED_K,
    lty = 3
  )
  
  legend(
    "topright",
    legend = c(
      "CPI",
      "Gap statistic"
    ),
    col = c(
      "steelblue4",
      "firebrick3"
    ),
    pch = c(
      19,
      17
    ),
    lty = 2,
    bty = "n"
  )
}

plot_silhouette_selected <- function(
    sil_object,
    cohort_label
) {
  ordered_rows <- order(
    sil_object[, "cluster"],
    -sil_object[, "sil_width"]
  )
  
  cluster_id <- sil_object[ordered_rows, "cluster"]
  widths <- sil_object[ordered_rows, "sil_width"]
  n_samples <- length(widths)
  
  stopifnot(
    n_samples > 0,
    all(is.finite(widths)),
    length(unique(cluster_id)) == 2L
  )
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(4.5, 9, 3.5, 1))
  
  plot(
    NA,
    xlim = c(min(-0.2, widths), 1),
    ylim = c(n_samples + 0.5, 0.5),
    xlab = "Silhouette width",
    ylab = "",
    yaxt = "n",
    xaxs = "i",
    main = paste(cohort_label, "silhouette plot (k = 2)")
  )
  
  bar_colours <- ifelse(
    cluster_id == 1,
    "#3975B7",
    "#D97850"
  )
  
  rect(
    xleft = pmin(0, widths),
    ybottom = seq_len(n_samples) - 0.5,
    xright = pmax(0, widths),
    ytop = seq_len(n_samples) + 0.5,
    col = bar_colours,
    border = NA
  )
  
  abline(v = 0, col = "grey30", lwd = 0.8)
  
  cluster_positions <- tapply(
    seq_len(n_samples),
    cluster_id,
    mean
  )
  
  cluster_labels <- vapply(
    names(cluster_positions),
    function(id) {
      values <- widths[cluster_id == as.numeric(id)]
      sprintf(
        "PCD_C%s (n = %d; mean = %.2f)",
        id,
        length(values),
        mean(values)
      )
    },
    character(1)
  )
  
  axis(
    side = 2,
    at = as.numeric(cluster_positions),
    labels = cluster_labels,
    las = 1,
    tick = FALSE,
    cex.axis = 0.85
  )
  
  mtext(
    sprintf("Overall mean silhouette width = %.3f", mean(widths)),
    side = 3,
    line = 0.3,
    adj = 0,
    cex = 0.85
  )
}

save_three_formats(
  "TCGA_LUAD_CPI_Gap_cluster_number",
  function() {
    plot_cpi_gap(
      luad_cpi,
      luad_gap_plot,
      "TCGA-LUAD"
    )
  }
)

save_three_formats(
  "TCGA_LUSC_CPI_Gap_cluster_number",
  function() {
    plot_cpi_gap(
      lusc_cpi,
      lusc_gap_plot,
      "TCGA-LUSC"
    )
  }
)

save_three_formats(
  "TCGA_LUAD_silhouette_selected_k2",
  function() {
    plot_silhouette_selected(
      luad_selected_silhouette,
      "TCGA-LUAD"
    )
  },
  width = 7,
  height = 6
)

save_three_formats(
  "TCGA_LUSC_silhouette_selected_k2",
  function() {
    plot_silhouette_selected(
      lusc_selected_silhouette,
      "TCGA-LUSC"
    )
  },
  width = 7,
  height = 6
)


# ------------------------------------------------------------
# 33. Do not force k = 2
#
# The script reports what the established diagnostics say.
# The thesis should describe k = 2 as supported only if the
# resulting CPI / Gap / silhouette patterns genuinely support
# the selected two-cluster solution.
# ------------------------------------------------------------

if (
  luad_best_cpi_k != SELECTED_K ||
  luad_gap_tibshirani_k != SELECTED_K ||
  luad_best_silhouette_k != SELECTED_K
) {
  
  warning(
    "LUAD established CPI/Gap/silhouette diagnostics do not ",
    "all support k = ",
    SELECTED_K,
    ". Review the generated curves before finalizing."
  )
}

if (
  lusc_best_cpi_k != SELECTED_K ||
  lusc_gap_tibshirani_k != SELECTED_K ||
  lusc_best_silhouette_k != SELECTED_K
) {
  
  warning(
    "LUSC established CPI/Gap/silhouette diagnostics do not ",
    "all support k = ",
    SELECTED_K,
    ". Review the generated curves before finalizing."
  )
}


# ============================================================
# PART G — FINAL SELECTED CLUSTER ASSIGNMENTS
# ============================================================


# ------------------------------------------------------------
# 27. LUAD assignments
# ------------------------------------------------------------

luad_cluster_assignment <- data.frame(
  patient_id =
    names(
      luad_ccp[[SELECTED_K]]$consensusClass
    ),
  
  PCD_cluster =
    paste0(
      "PCD_C",
      luad_ccp[[SELECTED_K]]$consensusClass
    ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 28. LUSC assignments
# ------------------------------------------------------------

lusc_cluster_assignment <- data.frame(
  patient_id =
    names(
      lusc_ccp[[SELECTED_K]]$consensusClass
    ),
  
  PCD_cluster =
    paste0(
      "PCD_C",
      lusc_ccp[[SELECTED_K]]$consensusClass
    ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 29. Validate final cluster sizes
# ------------------------------------------------------------

stopifnot(
  sum(
    luad_cluster_assignment$PCD_cluster ==
      "PCD_C1"
  ) == 284
)

stopifnot(
  sum(
    luad_cluster_assignment$PCD_cluster ==
      "PCD_C2"
  ) == 233
)

stopifnot(
  sum(
    lusc_cluster_assignment$PCD_cluster ==
      "PCD_C1"
  ) == 286
)

stopifnot(
  sum(
    lusc_cluster_assignment$PCD_cluster ==
      "PCD_C2"
  ) == 215
)


# ============================================================
# PART H — SURVIVAL CHARACTERIZATION
# ============================================================


# ------------------------------------------------------------
# 30. Load survival metadata
# ------------------------------------------------------------

luad_survival_metadata <- read.csv(
  file.path(
    survival_processed_dir,
    "TCGA_LUAD_PCD_survival_metadata.csv"
  ),
  stringsAsFactors = FALSE
)

lusc_survival_metadata <- read.csv(
  file.path(
    survival_processed_dir,
    "TCGA_LUSC_PCD_survival_metadata.csv"
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 31. LUAD cluster-survival dataset
# ------------------------------------------------------------

luad_cluster_survival <- luad_survival_metadata |>
  dplyr::select(
    patient_id = submitter_id,
    OS_time,
    OS_status
  ) |>
  dplyr::inner_join(
    luad_cluster_assignment,
    by = "patient_id"
  )

stopifnot(
  nrow(luad_cluster_survival) == 504
)


# ------------------------------------------------------------
# 32. LUAD log-rank
# ------------------------------------------------------------

luad_logrank <- survival::survdiff(
  survival::Surv(
    OS_time,
    OS_status
  ) ~ PCD_cluster,
  data = luad_cluster_survival
)

luad_logrank_p <- 1 - pchisq(
  luad_logrank$chisq,
  df = 1
)


# ------------------------------------------------------------
# 33. LUAD cluster Cox
# ------------------------------------------------------------

luad_cluster_cox <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_status
  ) ~ PCD_cluster,
  data = luad_cluster_survival
)


# ------------------------------------------------------------
# 34. LUSC cluster-survival dataset
# ------------------------------------------------------------

lusc_cluster_survival <- lusc_survival_metadata |>
  dplyr::select(
    patient_id = submitter_id,
    OS_time,
    OS_status
  ) |>
  dplyr::inner_join(
    lusc_cluster_assignment,
    by = "patient_id"
  )

stopifnot(
  nrow(lusc_cluster_survival) == 493
)


# ------------------------------------------------------------
# 35. LUSC log-rank
# ------------------------------------------------------------

lusc_logrank <- survival::survdiff(
  survival::Surv(
    OS_time,
    OS_status
  ) ~ PCD_cluster,
  data = lusc_cluster_survival
)

lusc_logrank_p <- 1 - pchisq(
  lusc_logrank$chisq,
  df = 1
)


# ------------------------------------------------------------
# 36. LUSC cluster Cox
# ------------------------------------------------------------

lusc_cluster_cox <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_status
  ) ~ PCD_cluster,
  data = lusc_cluster_survival
)


# ============================================================
# PART I — FULL 296-GENE PCD MECHANISM PROFILES
# ============================================================


# ------------------------------------------------------------
# 37. All curated PCD genes
# ------------------------------------------------------------

full_pcd_genes <- unique(
  pcd_master_gene$gene_symbol
)

stopifnot(
  length(full_pcd_genes) == 296
)


# ------------------------------------------------------------
# 38. Match full PCD catalogue
# ------------------------------------------------------------

luad_full_pcd_rows <- match(
  full_pcd_genes,
  luad_gene_annotation$gene_name
)

lusc_full_pcd_rows <- match(
  full_pcd_genes,
  lusc_gene_annotation$gene_name
)

stopifnot(
  !anyNA(luad_full_pcd_rows)
)

stopifnot(
  !anyNA(lusc_full_pcd_rows)
)


# ------------------------------------------------------------
# 39. LUAD full PCD expression
# ------------------------------------------------------------

luad_full_pcd_expression <- luad_tpm[
  luad_full_pcd_rows,
  luad_cluster_tpm_cols,
  drop = FALSE
]

rownames(
  luad_full_pcd_expression
) <- full_pcd_genes

colnames(
  luad_full_pcd_expression
) <- luad_tumor_one_patient$cases.submitter_id

luad_full_pcd_expression <-
  log2(
    luad_full_pcd_expression + 1
  )


# ------------------------------------------------------------
# 40. LUSC full PCD expression
# ------------------------------------------------------------

lusc_full_pcd_expression <- lusc_tpm[
  lusc_full_pcd_rows,
  lusc_cluster_tpm_cols,
  drop = FALSE
]

rownames(
  lusc_full_pcd_expression
) <- full_pcd_genes

colnames(
  lusc_full_pcd_expression
) <- lusc_tumor_one_patient$cases.submitter_id

lusc_full_pcd_expression <-
  log2(
    lusc_full_pcd_expression + 1
  )


# ------------------------------------------------------------
# 41. Standardize every gene
# ------------------------------------------------------------

luad_full_pcd_z <- t(
  scale(
    t(luad_full_pcd_expression)
  )
)

lusc_full_pcd_z <- t(
  scale(
    t(lusc_full_pcd_expression)
  )
)


# ============================================================
# PART J — MECHANISM SCORE FUNCTION
# ============================================================


calculate_full_pcd_scores <- function(
    z_matrix,
    membership_table,
    cluster_assignment
) {
  
  scores <- lapply(
    sort(
      unique(
        membership_table$PCD_type
      )
    ),
    function(pcd_type) {
      
      genes <- membership_table |>
        dplyr::filter(
          PCD_type == pcd_type
        ) |>
        pull(
          gene_symbol
        ) |>
        unique()
      
      genes <- intersect(
        genes,
        rownames(z_matrix)
      )
      
      data.frame(
        patient_id =
          colnames(z_matrix),
        
        PCD_type =
          pcd_type,
        
        n_genes =
          length(genes),
        
        score =
          colMeans(
            z_matrix[
              genes,
              ,
              drop = FALSE
            ]
          ),
        
        stringsAsFactors = FALSE
      )
    }
  ) |>
    dplyr::bind_rows() |>
    dplyr::left_join(
      cluster_assignment,
      by = "patient_id"
    )
  
  scores
}


# ------------------------------------------------------------
# 42. Calculate mechanism scores
# ------------------------------------------------------------

luad_full_pcd_scores <- calculate_full_pcd_scores(
  luad_full_pcd_z,
  pcd_full_membership,
  luad_cluster_assignment
)

lusc_full_pcd_scores <- calculate_full_pcd_scores(
  lusc_full_pcd_z,
  pcd_full_membership,
  lusc_cluster_assignment
)


# ------------------------------------------------------------
# 43. Validate score mapping
# ------------------------------------------------------------

stopifnot(
  sum(
    is.na(
      luad_full_pcd_scores$PCD_cluster
    )
  ) == 0
)

stopifnot(
  sum(
    is.na(
      lusc_full_pcd_scores$PCD_cluster
    )
  ) == 0
)


# ============================================================
# PART K — MECHANISM CLUSTER TESTS
# ============================================================


# ------------------------------------------------------------
# 44. LUAD
# ------------------------------------------------------------

luad_full_pcd_tests <- luad_full_pcd_scores |>
  dplyr::group_by(
    PCD_type
  ) |>
  dplyr::summarise(
    C1_mean =
      mean(
        score[
          PCD_cluster == "PCD_C1"
        ]
      ),
    
    C2_mean =
      mean(
        score[
          PCD_cluster == "PCD_C2"
        ]
      ),
    
    C2_minus_C1 =
      C2_mean -
      C1_mean,
    
    pvalue =
      wilcox.test(
        score ~ PCD_cluster
      )$p.value,
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    padj =
      p.adjust(
        pvalue,
        method = "BH"
      )
  ) |>
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# 45. LUSC
# ------------------------------------------------------------

lusc_full_pcd_tests <- lusc_full_pcd_scores |>
  dplyr::group_by(
    PCD_type
  ) |>
  dplyr::summarise(
    C1_mean =
      mean(
        score[
          PCD_cluster == "PCD_C1"
        ]
      ),
    
    C2_mean =
      mean(
        score[
          PCD_cluster == "PCD_C2"
        ]
      ),
    
    C2_minus_C1 =
      C2_mean -
      C1_mean,
    
    pvalue =
      wilcox.test(
        score ~ PCD_cluster
      )$p.value,
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    padj =
      p.adjust(
        pvalue,
        method = "BH"
      )
  ) |>
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# 46. Validate mechanism results
# ------------------------------------------------------------

stopifnot(
  sum(
    luad_full_pcd_tests$padj < 0.05
  ) == 2
)

stopifnot(
  all(
    c(
      "Necroptosis",
      "Pyroptosis"
    ) %in%
      luad_full_pcd_tests$PCD_type[
        luad_full_pcd_tests$padj < 0.05
      ]
  )
)

stopifnot(
  sum(
    lusc_full_pcd_tests$padj < 0.05
  ) == 3
)

stopifnot(
  all(
    c(
      "Ferroptosis",
      "Apoptosis",
      "Cuproptosis"
    ) %in%
      lusc_full_pcd_tests$PCD_type[
        lusc_full_pcd_tests$padj < 0.05
      ]
  )
)


# ============================================================
# PART L — SAVE EXACT FINAL OUTPUTS
# ============================================================


# ------------------------------------------------------------
# 47. Cluster assignments
# ------------------------------------------------------------

write.csv(
  luad_cluster_assignment,
  file.path(
    clustering_results_dir,
    "TCGA_LUAD_PCD_cluster_assignments.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_cluster_assignment,
  file.path(
    clustering_results_dir,
    "TCGA_LUSC_PCD_cluster_assignments.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 48. Patient-level mechanism scores
# ------------------------------------------------------------

write.csv(
  luad_full_pcd_scores,
  file.path(
    clustering_results_dir,
    "TCGA_LUAD_PCD_mechanism_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_full_pcd_scores,
  file.path(
    clustering_results_dir,
    "TCGA_LUSC_PCD_mechanism_scores.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 49. Mechanism comparison tests
# ------------------------------------------------------------

write.csv(
  luad_full_pcd_tests,
  file.path(
    clustering_results_dir,
    "TCGA_LUAD_PCD_mechanism_cluster_tests.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_full_pcd_tests,
  file.path(
    clustering_results_dir,
    "TCGA_LUSC_PCD_mechanism_cluster_tests.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART M — FINAL OUTPUT VALIDATION
# ============================================================


expected_files <- c(
  "TCGA_LUAD_PCD_cluster_assignments.csv",
  "TCGA_LUSC_PCD_cluster_assignments.csv",
  "TCGA_LUAD_PCD_mechanism_scores.csv",
  "TCGA_LUSC_PCD_mechanism_scores.csv",
  "TCGA_LUAD_PCD_mechanism_cluster_tests.csv",
  "TCGA_LUSC_PCD_mechanism_cluster_tests.csv"
)

stopifnot(
  all(
    file.exists(
      file.path(
        clustering_results_dir,
        expected_files
      )
    )
  )
)


# ============================================================
# PART N — COMPLETION SUMMARY
# ============================================================


luad_cluster_cox_summary <-
  summary(
    luad_cluster_cox
  )

lusc_cluster_cox_summary <-
  summary(
    lusc_cluster_cox
  )


cat(
  "\n====================================================\n"
)

cat(
  "PCD consensus clustering completed.\n"
)

cat(
  "====================================================\n\n"
)

cat(
  "Clustering gene set: 47 concordant PCD genes\n"
)

cat(
  "Consensus resamplings: 1000\n"
)

cat(
  "Final k: 2 for LUAD and LUSC\n\n"
)


cat(
  "LUAD cluster sizes:\n"
)

print(
  table(
    luad_cluster_assignment$PCD_cluster
  )
)

cat(
  "\nLUSC cluster sizes:\n"
)

print(
  table(
    lusc_cluster_assignment$PCD_cluster
  )
)


cat(
  "\n\nLUAD mean silhouette:\n"
)

print(
  luad_silhouette
)

cat(
  "\nLUSC mean silhouette:\n"
)

print(
  lusc_silhouette
)


cat(
  "\nLUAD overall-survival comparison:\n"
)

cat(
  "  Log-rank p =",
  luad_logrank_p,
  "\n"
)

cat(
  "  C2 vs C1 HR =",
  luad_cluster_cox_summary$coefficients[
    1,
    "exp(coef)"
  ],
  "\n"
)

cat(
  "  Cox p =",
  luad_cluster_cox_summary$coefficients[
    1,
    "Pr(>|z|)"
  ],
  "\n\n"
)


cat(
  "LUSC overall-survival comparison:\n"
)

cat(
  "  Log-rank p =",
  lusc_logrank_p,
  "\n"
)

cat(
  "  C2 vs C1 HR =",
  lusc_cluster_cox_summary$coefficients[
    1,
    "exp(coef)"
  ],
  "\n"
)

cat(
  "  Cox p =",
  lusc_cluster_cox_summary$coefficients[
    1,
    "Pr(>|z|)"
  ],
  "\n\n"
)


cat(
  "LUAD significant PCD mechanisms:\n"
)

print(
  luad_full_pcd_tests |>
    dplyr::filter(
      padj < 0.05
    )
)


cat(
  "\nLUSC significant PCD mechanisms:\n"
)

print(
  lusc_full_pcd_tests |>
    dplyr::filter(
      padj < 0.05
    )
)


cat(
  "\nFinal clustering CSV files: 6\n"
)

cat(
  "\nCluster-number diagnostics:\n",
  "LUAD best PAC k: ", luad_best_pac_k, "\n",
  "LUAD CPI peak k: ", luad_best_cpi_k, "\n",
  "LUAD Gap global maximum k: ", luad_gap_global_k, "\n",
  "LUAD Gap Tibshirani 1-SE k: ", luad_gap_tibshirani_k, "\n",
  "LUAD best silhouette k: ", luad_best_silhouette_k, "\n\n",
  "LUSC best PAC k: ", lusc_best_pac_k, "\n",
  "LUSC CPI peak k: ", lusc_best_cpi_k, "\n",
  "LUSC Gap global maximum k: ", lusc_gap_global_k, "\n",
  "LUSC Gap Tibshirani 1-SE k: ", lusc_gap_tibshirani_k, "\n",
  "LUSC best silhouette k: ", lusc_best_silhouette_k, "\n",
  sep = ""
)

cat(
  "\n========================================\n",
  "SCRIPT 14 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================