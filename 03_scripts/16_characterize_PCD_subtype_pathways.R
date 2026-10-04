# ============================================================
# Script 16: Characterize PCD subtype pathways
#
# Thesis:
# Exploring a Specialized Programmed Cell Death Pattern to
# Predict Prognosis and Treatment Sensitivity in Lung Cancer
# by Machine Learning and Multi-Omics Analysis
#
# Cohorts:
#   TCGA-LUAD
#   TCGA-LUSC
#
# Objective:
#   Characterize biological pathways associated with the
#   previously defined PCD clusters.
#
# Strategy:
#   1. Reconstruct one Primary Tumor sample per patient
#   2. Perform full-transcriptome differential expression
#      between PCD_C2 and PCD_C1 using DESeq2
#   3. Rank genes using the DESeq2 Wald statistic
#   4. Run Hallmark pre-ranked GSEA using fgsea
#   5. Identify significant and directionally concordant
#      pathways across LUAD and LUSC
#
# IMPORTANT DIRECTION:
#
#   Contrast = PCD_C2 vs PCD_C1
#
#   log2FC > 0 / NES > 0  -> enriched in PCD_C2
#   log2FC < 0 / NES < 0  -> enriched in PCD_C1
#
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "DESeq2",
  "dplyr",
  "tibble",
  "msigdbr",
  "fgsea"
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

enrichment_dir <- file.path(
  results_dir,
  "enrichment"
)

ensure_dir(enrichment_dir)


# ------------------------------------------------------------
# 3. Define required input files
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
    "One or more required files are missing:\n",
    paste(
      required_files[
        !file.exists(required_files)
      ],
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 4. Load processed expression data
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

luad_metadata <- readRDS(
  luad_metadata_file
)

lusc_metadata <- readRDS(
  lusc_metadata_file
)

luad_annotation <- readRDS(
  luad_annotation_file
)

lusc_annotation <- readRDS(
  lusc_annotation_file
)


# ------------------------------------------------------------
# 5. Load PCD cluster assignments
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
# 6. Validate sample-level expression objects
# ------------------------------------------------------------

stopifnot(
  identical(
    colnames(luad_counts),
    colnames(luad_tpm)
  )
)

stopifnot(
  identical(
    colnames(lusc_counts),
    colnames(lusc_tpm)
  )
)

stopifnot(
  identical(
    colnames(luad_counts),
    luad_metadata$cases
  )
)

stopifnot(
  identical(
    colnames(lusc_counts),
    lusc_metadata$cases
  )
)


# ------------------------------------------------------------
# 7. Select one Primary Tumor sample per patient
#
# Same biological selection rule used in Scripts 14 and 15:
#
#   - Primary Tumor only
#   - patient ID = first 12 characters of TCGA barcode
#   - if multiple samples exist for a patient, retain
#     the sample with the highest raw-count library size
#   - retain only patients belonging to the final PCD clusters
# ------------------------------------------------------------

select_cluster_tumors <- function(
    counts,
    tpm,
    metadata,
    cluster_assignment
) {
  
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
  
  counts_selected <- counts[
    ,
    sample_info$cases,
    drop = FALSE
  ]
  
  tpm_selected <- tpm[
    ,
    sample_info$cases,
    drop = FALSE
  ]
  
  colnames(counts_selected) <-
    sample_info$patient_id
  
  colnames(tpm_selected) <-
    sample_info$patient_id
  
  counts_selected <-
    counts_selected[
      ,
      cluster_assignment$patient_id,
      drop = FALSE
    ]
  
  tpm_selected <-
    tpm_selected[
      ,
      cluster_assignment$patient_id,
      drop = FALSE
    ]
  
  stopifnot(
    identical(
      colnames(counts_selected),
      cluster_assignment$patient_id
    )
  )
  
  stopifnot(
    identical(
      colnames(tpm_selected),
      cluster_assignment$patient_id
    )
  )
  
  list(
    counts = counts_selected,
    tpm = tpm_selected
  )
}


luad_cluster_data <-
  select_cluster_tumors(
    counts = luad_counts,
    tpm = luad_tpm,
    metadata = luad_metadata,
    cluster_assignment =
      luad_cluster_assignment
  )

lusc_cluster_data <-
  select_cluster_tumors(
    counts = lusc_counts,
    tpm = lusc_tpm,
    metadata = lusc_metadata,
    cluster_assignment =
      lusc_cluster_assignment
  )


# ------------------------------------------------------------
# 8. Validate reconstructed cluster datasets
# ------------------------------------------------------------

stopifnot(
  nrow(luad_cluster_data$counts) == 60660,
  ncol(luad_cluster_data$counts) == 517
)

stopifnot(
  nrow(lusc_cluster_data$counts) == 60660,
  ncol(lusc_cluster_data$counts) == 501
)

stopifnot(
  identical(
    colnames(
      luad_cluster_data$counts
    ),
    luad_cluster_assignment$patient_id
  )
)

stopifnot(
  identical(
    colnames(
      lusc_cluster_data$counts
    ),
    lusc_cluster_assignment$patient_id
  )
)

cat(
  "\nPatient-level expression reconstructed:\n"
)

cat(
  "LUAD:",
  nrow(luad_cluster_data$counts),
  "genes x",
  ncol(luad_cluster_data$counts),
  "patients\n"
)

cat(
  "LUSC:",
  nrow(lusc_cluster_data$counts),
  "genes x",
  ncol(lusc_cluster_data$counts),
  "patients\n"
)


# ============================================================
# PART A: FULL-TRANSCRIPTOME DIFFERENTIAL EXPRESSION
# ============================================================


# ------------------------------------------------------------
# 9. LUAD DESeq2 analysis
#
# Reference = PCD_C1
# Contrast  = PCD_C2 vs PCD_C1
# ------------------------------------------------------------

luad_coldata <- data.frame(
  
  row.names =
    luad_cluster_assignment$patient_id,
  
  PCD_cluster = factor(
    luad_cluster_assignment$PCD_cluster,
    levels = c(
      "PCD_C1",
      "PCD_C2"
    )
  )
)

stopifnot(
  identical(
    colnames(
      luad_cluster_data$counts
    ),
    rownames(
      luad_coldata
    )
  )
)

luad_dds <-
  DESeq2::DESeqDataSetFromMatrix(
    
    countData = round(
      luad_cluster_data$counts
    ),
    
    colData =
      luad_coldata,
    
    design =
      ~ PCD_cluster
  )

luad_keep <- rowSums(
  DESeq2::counts(
    luad_dds
  )
) >= 10

luad_dds <- luad_dds[
  luad_keep,
]

luad_dds <- DESeq2::DESeq(
  luad_dds,
  quiet = TRUE
)

luad_res <- DESeq2::results(
  luad_dds,
  contrast = c(
    "PCD_cluster",
    "PCD_C2",
    "PCD_C1"
  ),
  alpha = 0.05
)


# ------------------------------------------------------------
# 10. LUSC DESeq2 analysis
# ------------------------------------------------------------

lusc_coldata <- data.frame(
  
  row.names =
    lusc_cluster_assignment$patient_id,
  
  PCD_cluster = factor(
    lusc_cluster_assignment$PCD_cluster,
    levels = c(
      "PCD_C1",
      "PCD_C2"
    )
  )
)

stopifnot(
  identical(
    colnames(
      lusc_cluster_data$counts
    ),
    rownames(
      lusc_coldata
    )
  )
)

lusc_dds <-
  DESeq2::DESeqDataSetFromMatrix(
    
    countData = round(
      lusc_cluster_data$counts
    ),
    
    colData =
      lusc_coldata,
    
    design =
      ~ PCD_cluster
  )

lusc_keep <- rowSums(
  DESeq2::counts(
    lusc_dds
  )
) >= 10

lusc_dds <- lusc_dds[
  lusc_keep,
]

lusc_dds <- DESeq2::DESeq(
  lusc_dds,
  quiet = TRUE
)

lusc_res <- DESeq2::results(
  lusc_dds,
  contrast = c(
    "PCD_cluster",
    "PCD_C2",
    "PCD_C1"
  ),
  alpha = 0.05
)


# ------------------------------------------------------------
# 11. Convert DESeq2 results to data frames
# ------------------------------------------------------------

luad_de <- as.data.frame(
  luad_res
) |>
  tibble::rownames_to_column(
    "gene_id"
  )

lusc_de <- as.data.frame(
  lusc_res
) |>
  tibble::rownames_to_column(
    "gene_id"
  )


# ------------------------------------------------------------
# 12. Attach gene annotation
# ------------------------------------------------------------

luad_de <- luad_de |>
  dplyr::left_join(
    
    luad_annotation |>
      dplyr::select(
        gene_id,
        gene_name,
        gene_type,
        ensembl_id
      ),
    
    by = "gene_id"
  )

lusc_de <- lusc_de |>
  dplyr::left_join(
    
    lusc_annotation |>
      dplyr::select(
        gene_id,
        gene_name,
        gene_type,
        ensembl_id
      ),
    
    by = "gene_id"
  )

stopifnot(
  sum(
    is.na(
      luad_de$gene_name
    )
  ) == 0
)

stopifnot(
  sum(
    is.na(
      lusc_de$gene_name
    )
  ) == 0
)


# ------------------------------------------------------------
# 13. Define stringent cluster DE genes
#
# Criteria:
#   BH-adjusted p-value < 0.05
#   absolute log2 fold change >= 1
# ------------------------------------------------------------

luad_sig_cluster_de <- luad_de |>
  dplyr::filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )

lusc_sig_cluster_de <- lusc_de |>
  dplyr::filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )


# ------------------------------------------------------------
# 14. Differential-expression summary
# ------------------------------------------------------------

luad_de_summary <- c(
  
  significant =
    nrow(
      luad_sig_cluster_de
    ),
  
  C2_up = sum(
    luad_sig_cluster_de$log2FoldChange > 0
  ),
  
  C1_up = sum(
    luad_sig_cluster_de$log2FoldChange < 0
  )
)

lusc_de_summary <- c(
  
  significant =
    nrow(
      lusc_sig_cluster_de
    ),
  
  C2_up = sum(
    lusc_sig_cluster_de$log2FoldChange > 0
  ),
  
  C1_up = sum(
    lusc_sig_cluster_de$log2FoldChange < 0
  )
)


# ============================================================
# PART B: PREPARE GSEA RANKINGS
# ============================================================


# ------------------------------------------------------------
# 15. Create one Wald statistic per gene symbol
#
# GSEA ranking metric:
#   DESeq2 Wald statistic
#
# For duplicated gene symbols:
#   retain the row having the largest absolute Wald statistic.
# ------------------------------------------------------------

prepare_gsea_rank <- function(
    de_table
) {
  
  ranked <- de_table |>
    dplyr::filter(
      !is.na(gene_name),
      gene_name != "",
      !is.na(stat)
    ) |>
    dplyr::mutate(
      abs_stat =
        abs(stat)
    ) |>
    dplyr::arrange(
      gene_name,
      dplyr::desc(
        abs_stat
      )
    ) |>
    dplyr::group_by(
      gene_name
    ) |>
    dplyr::slice_head(
      n = 1
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      dplyr::desc(
        stat
      )
    )
  
  gene_list <-
    ranked$stat
  
  names(gene_list) <-
    ranked$gene_name
  
  gene_list
}


luad_ranked_genes <-
  prepare_gsea_rank(
    luad_de
  )

lusc_ranked_genes <-
  prepare_gsea_rank(
    lusc_de
  )


# ------------------------------------------------------------
# 16. Validate ranked gene vectors
# ------------------------------------------------------------

stopifnot(
  length(luad_ranked_genes) > 0,
  length(lusc_ranked_genes) > 0
)

cat(
  "\nGSEA ranked-gene counts:\n",
  "LUAD: ", length(luad_ranked_genes), "\n",
  "LUSC: ", length(lusc_ranked_genes), "\n",
  sep = ""
)

stopifnot(
  sum(
    duplicated(
      names(
        luad_ranked_genes
      )
    )
  ) == 0
)

stopifnot(
  sum(
    duplicated(
      names(
        lusc_ranked_genes
      )
    )
  ) == 0
)

stopifnot(
  !anyNA(
    luad_ranked_genes
  )
)

stopifnot(
  !anyNA(
    lusc_ranked_genes
  )
)


# ============================================================
# PART C: HALLMARK GSEA
# ============================================================


# ------------------------------------------------------------
# 17. Retrieve human Hallmark gene sets
# ------------------------------------------------------------

hallmark_df <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)

hallmark_provenance <- data.frame(
  analysis = "Hallmark GSEA",
  species = "Homo sapiens",
  collection = "H",
  n_gene_sets = dplyr::n_distinct(hallmark_df$gs_name),
  msigdbr_version = as.character(
    utils::packageVersion("msigdbr")
  ),
  fgsea_version = as.character(
    utils::packageVersion("fgsea")
  ),
  R_version = R.version.string,
  analysis_date = as.character(Sys.Date()),
  stringsAsFactors = FALSE
)

write.csv(
  hallmark_provenance,
  file.path(
    enrichment_dir,
    "Hallmark_GSEA_provenance.csv"
  ),
  row.names = FALSE
)

stopifnot(
  length(
    unique(
      hallmark_df$gs_name
    )
  ) == 50
)


# ------------------------------------------------------------
# 18. Convert Hallmark collection to pathway list
# ------------------------------------------------------------

hallmark_pathways <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

stopifnot(
  length(
    hallmark_pathways
  ) == 50
)


# ------------------------------------------------------------
# 19. Run LUAD Hallmark GSEA
#
# NES > 0 -> PCD_C2 enriched
# NES < 0 -> PCD_C1 enriched
# ------------------------------------------------------------

set.seed(
  12345
)

luad_hallmark_gsea <- fgsea::fgsea(
  pathways =
    hallmark_pathways,
  stats =
    luad_ranked_genes,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

luad_hallmark_gsea <-
  as.data.frame(
    luad_hallmark_gsea
  ) |>
  dplyr::mutate(
    enriched_cluster = dplyr::if_else(
      NES > 0,
      "PCD_C2",
      "PCD_C1"
    )
  ) |>
  dplyr::arrange(
    padj,
    dplyr::desc(
      abs(NES)
    )
  )


# ------------------------------------------------------------
# 20. Run LUSC Hallmark GSEA
# ------------------------------------------------------------

set.seed(
  12345
)

lusc_hallmark_gsea <- fgsea::fgsea(
  pathways =
    hallmark_pathways,
  stats =
    lusc_ranked_genes,
  minSize = 15,
  maxSize = 500,
  eps = 0
)

lusc_hallmark_gsea <-
  as.data.frame(
    lusc_hallmark_gsea
  ) |>
  dplyr::mutate(
    enriched_cluster = dplyr::if_else(
      NES > 0,
      "PCD_C2",
      "PCD_C1"
    )
  ) |>
  dplyr::arrange(
    padj,
    dplyr::desc(
      abs(NES)
    )
  )


# ------------------------------------------------------------
# 21. Significant Hallmark pathways
# ------------------------------------------------------------

luad_hallmark_sig <-
  luad_hallmark_gsea |>
  dplyr::filter(
    padj < 0.05
  )

lusc_hallmark_sig <-
  lusc_hallmark_gsea |>
  dplyr::filter(
    padj < 0.05
  )


# ------------------------------------------------------------
# 22. Cross-cohort Hallmark comparison
# ------------------------------------------------------------

hallmark_cross_cohort <-
  luad_hallmark_gsea |>
  dplyr::select(
    pathway,
    LUAD_NES = NES,
    LUAD_padj = padj
  ) |>
  dplyr::left_join(
    
    lusc_hallmark_gsea |>
      dplyr::select(
        pathway,
        LUSC_NES = NES,
        LUSC_padj = padj
      ),
    
    by = "pathway"
  ) |>
  dplyr::mutate(
    
    LUAD_significant =
      LUAD_padj < 0.05,
    
    LUSC_significant =
      LUSC_padj < 0.05,
    
    significant_both =
      LUAD_significant &
      LUSC_significant,
    
    direction_concordant =
      sign(
        LUAD_NES
      ) ==
      sign(
        LUSC_NES
      ),
    
    robust_cross_cohort =
      significant_both &
      direction_concordant
  ) |>
  dplyr::arrange(
    dplyr::desc(
      robust_cross_cohort
    ),
    LUAD_padj
  )


# ------------------------------------------------------------
# 23. Robust cross-cohort pathways
# ------------------------------------------------------------

hallmark_robust_cross_cohort <-
  hallmark_cross_cohort |>
  dplyr::filter(
    robust_cross_cohort
  )


# ------------------------------------------------------------
# 24. Add biological direction label
# ------------------------------------------------------------

hallmark_robust_cross_cohort <-
  hallmark_robust_cross_cohort |>
  dplyr::mutate(
    
    enriched_cluster =
      dplyr::if_else(
        LUAD_NES > 0,
        "PCD_C2",
        "PCD_C1"
      )
  )


# ============================================================
# PART D: SAVE RESULTS
# ============================================================


# ------------------------------------------------------------
# 25. Helper to convert fgsea leadingEdge list-column to text
# ------------------------------------------------------------

prepare_fgsea_for_csv <- function(
    x
) {
  
  x <- as.data.frame(
    x
  )
  
  if (
    "leadingEdge" %in%
    colnames(x)
  ) {
    
    x$leadingEdge <- vapply(
      x$leadingEdge,
      function(genes) {
        
        paste(
          genes,
          collapse = ";"
        )
      },
      character(1)
    )
  }
  
  x
}


luad_hallmark_out <-
  prepare_fgsea_for_csv(
    luad_hallmark_gsea
  )

lusc_hallmark_out <-
  prepare_fgsea_for_csv(
    lusc_hallmark_gsea
  )


# ------------------------------------------------------------
# 26. Save full DESeq2 results
# ------------------------------------------------------------

write.csv(
  luad_de,
  file.path(
    enrichment_dir,
    "TCGA_LUAD_PCD_cluster_full_DE.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_de,
  file.path(
    enrichment_dir,
    "TCGA_LUSC_PCD_cluster_full_DE.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Save stringent DE genes
# ------------------------------------------------------------

write.csv(
  luad_sig_cluster_de,
  file.path(
    enrichment_dir,
    "TCGA_LUAD_PCD_cluster_significant_DE.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_sig_cluster_de,
  file.path(
    enrichment_dir,
    "TCGA_LUSC_PCD_cluster_significant_DE.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Save complete Hallmark GSEA results
# ------------------------------------------------------------

write.csv(
  luad_hallmark_out,
  file.path(
    enrichment_dir,
    "TCGA_LUAD_Hallmark_GSEA.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_hallmark_out,
  file.path(
    enrichment_dir,
    "TCGA_LUSC_Hallmark_GSEA.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Save significant Hallmark pathways
# ------------------------------------------------------------

write.csv(
  luad_hallmark_out |>
    dplyr::filter(
      padj < 0.05
    ),
  file.path(
    enrichment_dir,
    "TCGA_LUAD_Hallmark_GSEA_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_hallmark_out |>
    dplyr::filter(
      padj < 0.05
    ),
  file.path(
    enrichment_dir,
    "TCGA_LUSC_Hallmark_GSEA_significant.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Save cross-cohort pathway comparison
# ------------------------------------------------------------

write.csv(
  hallmark_cross_cohort,
  file.path(
    enrichment_dir,
    "TCGA_LUAD_LUSC_Hallmark_cross_cohort.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Save robust cross-cohort pathways
# ------------------------------------------------------------

write.csv(
  hallmark_robust_cross_cohort,
  file.path(
    enrichment_dir,
    "TCGA_LUAD_LUSC_Hallmark_robust_pathways.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART E: FINAL VALIDATION
# ============================================================


# ------------------------------------------------------------
# 32. Final summary
# ------------------------------------------------------------

cat(
  "\n==========================================================\n"
)

cat(
  "FUNCTIONAL ENRICHMENT ANALYSIS COMPLETE\n"
)

cat(
  "==========================================================\n"
)

cat(
  "\nCluster counts:\n"
)

cat(
  "\nLUAD:\n"
)

print(
  table(
    luad_cluster_assignment$PCD_cluster
  )
)

cat(
  "\nLUSC:\n"
)

print(
  table(
    lusc_cluster_assignment$PCD_cluster
  )
)


cat(
  "\nFull-transcriptome DE summary:\n"
)

cat(
  "\nLUAD:\n"
)

print(
  luad_de_summary
)

cat(
  "\nLUSC:\n"
)

print(
  lusc_de_summary
)


cat(
  "\nRanked genes:\n"
)

cat(
  "LUAD:",
  length(
    luad_ranked_genes
  ),
  "\n"
)

cat(
  "LUSC:",
  length(
    lusc_ranked_genes
  ),
  "\n"
)


cat(
  "\nSignificant Hallmark pathways:\n"
)

cat(
  "LUAD:",
  nrow(
    luad_hallmark_sig
  ),
  "\n"
)

cat(
  "LUSC:",
  nrow(
    lusc_hallmark_sig
  ),
  "\n"
)


cat(
  "\nRobust cross-cohort Hallmark pathways:",
  nrow(
    hallmark_robust_cross_cohort
  ),
  "\n"
)


cat(
  "\nRobust pathways enriched in PCD_C1:",
  sum(
    hallmark_robust_cross_cohort$enriched_cluster ==
      "PCD_C1"
  ),
  "\n"
)

cat(
  "Robust pathways enriched in PCD_C2:",
  sum(
    hallmark_robust_cross_cohort$enriched_cluster ==
      "PCD_C2"
  ),
  "\n"
)


cat(
  "\nRobust cross-cohort pathways:\n"
)

print(
  hallmark_robust_cross_cohort |>
    dplyr::select(
      pathway,
      LUAD_NES,
      LUAD_padj,
      LUSC_NES,
      LUSC_padj,
      enriched_cluster
    )
)


cat(
  "\nResults saved to:\n"
)

cat(
  enrichment_dir,
  "\n"
)

cat(
  "\n==========================================================\n"
)

cat(
  "\n========================================\n",
  "SCRIPT 16 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "Outputs: ",
  enrichment_dir,
  "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================