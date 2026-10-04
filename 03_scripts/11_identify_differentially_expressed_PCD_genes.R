# ============================================================
# 11_identify_differentially_expressed_PCD_genes.R
#
# Purpose:
# Differential-expression analysis of curated PCD genes
# in TCGA-LUAD and TCGA-LUSC.
#
# Workflow:
# 1. Load full RNA-seq raw-count matrices and metadata
# 2. Collapse technical RNA-seq aliquots
# 3. Retain one LUAD Primary Tumor sample per patient
# 4. Run transcriptome-wide DESeq2
# 5. Extract curated PCD-gene results
# 6. Apply:
#       padj < 0.05
#       |log2FoldChange| >= 1
# 7. Compare LUAD and LUSC
# 8. Save final differential-expression outputs
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

pcd_processed_dir <- file.path(
  processed_dir,
  "PCD_genes"
)

pcd_raw_dir <- file.path(
  raw_dir,
  "PCD_genes"
)

de_results_dir <- file.path(
  results_dir,
  "differential_expression"
)

ensure_dir(de_results_dir)


# ------------------------------------------------------------
# 2. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "DESeq2",
  "dplyr",
  "tidyr"
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
  library(DESeq2)
  library(dplyr)
  library(tidyr)
})

# NOTE:
# Bioconductor packages loaded by DESeq2 can export functions with names
# that overlap tidyverse verbs. Downstream calls therefore use explicit
# dplyr:: / tidyr:: namespaces to prevent function masking.


# ------------------------------------------------------------
# 2A. Required input files
# ------------------------------------------------------------

required_input_files <- c(
  file.path(pcd_processed_dir, "PCD_master_unique_genes.csv"),
  file.path(pcd_raw_dir, "PCD_gene_sources_master_long.csv"),
  file.path(expression_processed_dir, "TCGA_LUAD_gene_annotation.rds"),
  file.path(expression_processed_dir, "TCGA_LUSC_gene_annotation.rds"),
  file.path(expression_processed_dir, "TCGA_LUAD_PCD_gene_mapping.csv"),
  file.path(expression_processed_dir, "TCGA_LUSC_PCD_gene_mapping.csv"),
  file.path(expression_processed_dir, "TCGA_LUAD_raw_counts.rds"),
  file.path(expression_processed_dir, "TCGA_LUSC_raw_counts.rds"),
  file.path(expression_processed_dir, "TCGA_LUAD_sample_metadata.rds"),
  file.path(expression_processed_dir, "TCGA_LUSC_sample_metadata.rds")
)

check_files_exist(
  required_input_files,
  label = "Script 11 input file(s)"
)


# ------------------------------------------------------------
# 3. Load PCD reference tables
# ------------------------------------------------------------

pcd_master_gene <- read.csv(
  file.path(
    pcd_processed_dir,
    "PCD_master_unique_genes.csv"
  ),
  stringsAsFactors = FALSE
)

pcd_master_long_final <- read.csv(
  file.path(
    pcd_raw_dir,
    "PCD_gene_sources_master_long.csv"
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(pcd_master_gene) == 296
)

stopifnot(
  !anyDuplicated(
    pcd_master_gene$gene_symbol
  )
)


# ------------------------------------------------------------
# 4. Load RNA-seq gene annotations
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
# 5. Load PCD gene mappings
# ------------------------------------------------------------

luad_pcd_match <- read.csv(
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_PCD_gene_mapping.csv"
  ),
  stringsAsFactors = FALSE
)

lusc_pcd_match <- read.csv(
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_PCD_gene_mapping.csv"
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(luad_pcd_match) == 296,
  nrow(lusc_pcd_match) == 296
)


# ------------------------------------------------------------
# 6. Load full raw-count matrices
# ------------------------------------------------------------

luad_counts <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_raw_counts.rds"
  )
)

lusc_counts <- readRDS(
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_raw_counts.rds"
  )
)


# ------------------------------------------------------------
# 7. Load sample metadata
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
# 8. Validate original count/metadata alignment
# ------------------------------------------------------------

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


# ============================================================
# PART A — PREPARE LUAD
# ============================================================


# ------------------------------------------------------------
# 9. Collapse LUAD technical aliquots
# ------------------------------------------------------------

luad_sample_ids <-
  luad_sample_metadata$sample.submitter_id

luad_counts_collapsed <- t(
  rowsum(
    t(luad_counts),
    group = luad_sample_ids,
    reorder = FALSE
  )
)

luad_metadata_collapsed <-
  luad_sample_metadata |>
  dplyr::distinct(
    sample.submitter_id,
    .keep_all = TRUE
  )

luad_metadata_collapsed <-
  luad_metadata_collapsed[
    match(
      colnames(luad_counts_collapsed),
      luad_metadata_collapsed$sample.submitter_id
    ),
  ]

stopifnot(
  identical(
    colnames(luad_counts_collapsed),
    luad_metadata_collapsed$sample.submitter_id
  )
)


# ------------------------------------------------------------
# 10. Calculate LUAD library sizes
# ------------------------------------------------------------

luad_metadata_collapsed$library_size <-
  colSums(
    luad_counts_collapsed
  )


# ------------------------------------------------------------
# 11. Select one LUAD Primary Tumor per patient
#
# Rule:
# If more than one distinct Primary Tumor sample exists
# for a patient, retain the sample with the largest
# raw-count library size.
# ------------------------------------------------------------

luad_tumor_keep <-
  luad_metadata_collapsed |>
  dplyr::filter(
    sample_type == "Primary Tumor"
  ) |>
  group_by(
    cases.submitter_id
  ) |>
  slice_max(
    order_by = library_size,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup()

luad_normal_keep <-
  luad_metadata_collapsed |>
  dplyr::filter(
    sample_type == "Solid Tissue Normal"
  )

luad_metadata_deseq <- bind_rows(
  luad_tumor_keep,
  luad_normal_keep
)


# ------------------------------------------------------------
# 12. Record removed LUAD tumor samples
#
# This is a preprocessing audit file, not one of the
# 16 differential-expression result files.
# ------------------------------------------------------------

luad_removed_biological_samples <-
  luad_metadata_collapsed |>
  dplyr::filter(
    sample_type == "Primary Tumor",
    !sample.submitter_id %in%
      luad_tumor_keep$sample.submitter_id
  )

write.csv(
  luad_removed_biological_samples,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_removed_duplicate_tumor_samples.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Reorder LUAD samples
# ------------------------------------------------------------

luad_keep_samples <- colnames(
  luad_counts_collapsed
)[
  colnames(luad_counts_collapsed) %in%
    luad_metadata_deseq$sample.submitter_id
]

luad_metadata_deseq <-
  luad_metadata_deseq[
    match(
      luad_keep_samples,
      luad_metadata_deseq$sample.submitter_id
    ),
  ]

luad_counts_deseq <-
  luad_counts_collapsed[
    ,
    luad_keep_samples,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 14. Validate final LUAD sample set
# ------------------------------------------------------------

stopifnot(
  ncol(luad_counts_deseq) == 576
)

stopifnot(
  sum(
    luad_metadata_deseq$sample_type ==
      "Primary Tumor"
  ) == 517
)

stopifnot(
  sum(
    luad_metadata_deseq$sample_type ==
      "Solid Tissue Normal"
  ) == 59
)

stopifnot(
  identical(
    colnames(luad_counts_deseq),
    luad_metadata_deseq$sample.submitter_id
  )
)

stopifnot(
  !anyDuplicated(
    luad_metadata_deseq$cases.submitter_id[
      luad_metadata_deseq$sample_type ==
        "Primary Tumor"
    ]
  )
)


# ============================================================
# PART B — LUAD DESEQ2
# ============================================================


# ------------------------------------------------------------
# 15. LUAD DESeq2 metadata
# ------------------------------------------------------------

luad_coldata <- data.frame(
  sample_id =
    luad_metadata_deseq$sample.submitter_id,
  
  condition =
    ifelse(
      luad_metadata_deseq$sample_type ==
        "Primary Tumor",
      "Tumor",
      "Normal"
    ),
  
  row.names =
    luad_metadata_deseq$sample.submitter_id,
  
  stringsAsFactors = FALSE
)

luad_coldata$condition <- factor(
  luad_coldata$condition,
  levels = c(
    "Normal",
    "Tumor"
  )
)

stopifnot(
  identical(
    colnames(luad_counts_deseq),
    rownames(luad_coldata)
  )
)


# ------------------------------------------------------------
# 16. Create LUAD DESeq2 object
# ------------------------------------------------------------

dds_luad <- DESeqDataSetFromMatrix(
  countData = round(
    luad_counts_deseq
  ),
  colData = luad_coldata,
  design = ~ condition
)


# ------------------------------------------------------------
# 17. LUAD low-expression filtering
#
# Keep genes with >=10 counts in >=10 samples.
# ------------------------------------------------------------

luad_keep_genes <- rowSums(
  counts(dds_luad) >= 10
) >= 10

dds_luad <- dds_luad[
  luad_keep_genes,
]

stopifnot(
  nrow(dds_luad) == 33311
)


# ------------------------------------------------------------
# 18. Run LUAD DESeq2
# ------------------------------------------------------------

dds_luad <- DESeq(
  dds_luad
)

res_luad <- results(
  dds_luad,
  contrast = c(
    "condition",
    "Tumor",
    "Normal"
  ),
  alpha = 0.05
)


# ------------------------------------------------------------
# 19. Annotate LUAD results
# ------------------------------------------------------------

luad_de_all <- as.data.frame(
  res_luad
)

luad_de_all$gene_id <- rownames(
  luad_de_all
)

luad_de_all <- luad_de_all |>
  dplyr::left_join(
    luad_gene_annotation,
    by = "gene_id"
  )

stopifnot(
  sum(
    is.na(
      luad_de_all$gene_name
    )
  ) == 0
)


# ------------------------------------------------------------
# 20. Identify LUAD PCD genes filtered out
# ------------------------------------------------------------

luad_filtered_gene_ids <- rownames(
  dds_luad
)

luad_pcd_filtered_out <- luad_pcd_match[
  !luad_pcd_match$gene_id %in%
    luad_filtered_gene_ids,
]

stopifnot(
  nrow(luad_pcd_filtered_out) == 2
)

stopifnot(
  setequal(
    luad_pcd_filtered_out$gene_symbol,
    c(
      "FTMT",
      "LGALS13"
    )
  )
)


# ------------------------------------------------------------
# 21. Extract LUAD PCD DE results
# ------------------------------------------------------------

luad_pcd_de <- luad_de_all |>
  dplyr::filter(
    gene_name %in%
      pcd_master_gene$gene_symbol
  ) |>
  dplyr::left_join(
    pcd_master_gene,
    by = c(
      "gene_name" = "gene_symbol"
    )
  ) |>
  dplyr::arrange(
    padj
  )

stopifnot(
  nrow(luad_pcd_de) == 294
)


# ------------------------------------------------------------
# 22. Classify LUAD differential expression
#
# Significant:
# padj < 0.05 AND |log2FC| >= 1
# ------------------------------------------------------------

luad_pcd_de <- luad_pcd_de |>
  dplyr::mutate(
    DE_status = case_when(
      
      padj < 0.05 &
        log2FoldChange >= 1 ~
        "Up",
      
      padj < 0.05 &
        log2FoldChange <= -1 ~
        "Down",
      
      TRUE ~
        "Not significant"
    )
  )

stopifnot(
  sum(
    luad_pcd_de$DE_status == "Up"
  ) == 36
)

stopifnot(
  sum(
    luad_pcd_de$DE_status == "Down"
  ) == 21
)


# ------------------------------------------------------------
# 23. LUAD significant PCD genes
# ------------------------------------------------------------

luad_pcd_significant <- luad_pcd_de |>
  dplyr::filter(
    DE_status != "Not significant"
  ) |>
  dplyr::arrange(
    padj
  )

stopifnot(
  nrow(luad_pcd_significant) == 57
)


# ------------------------------------------------------------
# 24. LUAD summary by PCD type
# ------------------------------------------------------------

luad_pcd_de_by_type <- luad_pcd_de |>
  dplyr::select(
    gene_name,
    log2FoldChange,
    pvalue,
    padj,
    DE_status
  ) |>
  dplyr::left_join(
    pcd_master_long_final |>
      dplyr::distinct(
        gene_symbol,
        PCD_type
      ),
    by = c(
      "gene_name" = "gene_symbol"
    )
  )

luad_pcd_de_summary <- luad_pcd_de_by_type |>
  dplyr::count(
    PCD_type,
    DE_status
  ) |>
  tidyr::pivot_wider(
    names_from = DE_status,
    values_from = n,
    values_fill = 0
  )


# ============================================================
# PART C — PREPARE LUSC
# ============================================================


# ------------------------------------------------------------
# 25. Collapse LUSC technical aliquots
# ------------------------------------------------------------

lusc_sample_ids <-
  lusc_sample_metadata$sample.submitter_id

lusc_counts_collapsed <- t(
  rowsum(
    t(lusc_counts),
    group = lusc_sample_ids,
    reorder = FALSE
  )
)

lusc_metadata_collapsed <-
  lusc_sample_metadata |>
  dplyr::distinct(
    sample.submitter_id,
    .keep_all = TRUE
  )

lusc_metadata_collapsed <-
  lusc_metadata_collapsed[
    match(
      colnames(lusc_counts_collapsed),
      lusc_metadata_collapsed$sample.submitter_id
    ),
  ]

stopifnot(
  identical(
    colnames(lusc_counts_collapsed),
    lusc_metadata_collapsed$sample.submitter_id
  )
)


# ------------------------------------------------------------
# 26. LUSC final sample set
# ------------------------------------------------------------

lusc_counts_deseq <-
  lusc_counts_collapsed

lusc_metadata_deseq <-
  lusc_metadata_collapsed

stopifnot(
  ncol(lusc_counts_deseq) == 552
)

stopifnot(
  sum(
    lusc_metadata_deseq$sample_type ==
      "Primary Tumor"
  ) == 501
)

stopifnot(
  sum(
    lusc_metadata_deseq$sample_type ==
      "Solid Tissue Normal"
  ) == 51
)


# ============================================================
# PART D — LUSC DESEQ2
# ============================================================


# ------------------------------------------------------------
# 27. LUSC DESeq2 metadata
# ------------------------------------------------------------

lusc_coldata <- data.frame(
  sample_id =
    lusc_metadata_deseq$sample.submitter_id,
  
  condition =
    ifelse(
      lusc_metadata_deseq$sample_type ==
        "Primary Tumor",
      "Tumor",
      "Normal"
    ),
  
  row.names =
    lusc_metadata_deseq$sample.submitter_id,
  
  stringsAsFactors = FALSE
)

lusc_coldata$condition <- factor(
  lusc_coldata$condition,
  levels = c(
    "Normal",
    "Tumor"
  )
)

stopifnot(
  identical(
    colnames(lusc_counts_deseq),
    rownames(lusc_coldata)
  )
)


# ------------------------------------------------------------
# 28. Create LUSC DESeq2 object
# ------------------------------------------------------------

dds_lusc <- DESeqDataSetFromMatrix(
  countData = round(
    lusc_counts_deseq
  ),
  colData = lusc_coldata,
  design = ~ condition
)


# ------------------------------------------------------------
# 29. LUSC low-expression filtering
# ------------------------------------------------------------

lusc_keep_genes <- rowSums(
  counts(dds_lusc) >= 10
) >= 10

dds_lusc <- dds_lusc[
  lusc_keep_genes,
]

stopifnot(
  nrow(dds_lusc) == 34495
)


# ------------------------------------------------------------
# 30. Run LUSC DESeq2
# ------------------------------------------------------------

dds_lusc <- DESeq(
  dds_lusc
)

res_lusc <- results(
  dds_lusc,
  contrast = c(
    "condition",
    "Tumor",
    "Normal"
  ),
  alpha = 0.05
)


# ------------------------------------------------------------
# 31. Annotate LUSC results
# ------------------------------------------------------------

lusc_de_all <- as.data.frame(
  res_lusc
)

lusc_de_all$gene_id <- rownames(
  lusc_de_all
)

lusc_de_all <- lusc_de_all |>
  dplyr::left_join(
    lusc_gene_annotation,
    by = "gene_id"
  )

stopifnot(
  sum(
    is.na(
      lusc_de_all$gene_name
    )
  ) == 0
)


# ------------------------------------------------------------
# 32. Identify LUSC PCD genes filtered out
# ------------------------------------------------------------

lusc_filtered_gene_ids <- rownames(
  dds_lusc
)

lusc_pcd_filtered_out <- lusc_pcd_match[
  !lusc_pcd_match$gene_id %in%
    lusc_filtered_gene_ids,
]

stopifnot(
  nrow(lusc_pcd_filtered_out) == 2
)

stopifnot(
  setequal(
    lusc_pcd_filtered_out$gene_symbol,
    c(
      "FTMT",
      "LGALS13"
    )
  )
)


# ------------------------------------------------------------
# 33. Extract LUSC PCD DE results
# ------------------------------------------------------------

lusc_pcd_de <- lusc_de_all |>
  dplyr::filter(
    gene_name %in%
      pcd_master_gene$gene_symbol
  ) |>
  dplyr::left_join(
    pcd_master_gene,
    by = c(
      "gene_name" = "gene_symbol"
    )
  ) |>
  dplyr::arrange(
    padj
  )

stopifnot(
  nrow(lusc_pcd_de) == 294
)


# ------------------------------------------------------------
# 34. Classify LUSC differential expression
# ------------------------------------------------------------

lusc_pcd_de <- lusc_pcd_de |>
  dplyr::mutate(
    DE_status = case_when(
      
      padj < 0.05 &
        log2FoldChange >= 1 ~
        "Up",
      
      padj < 0.05 &
        log2FoldChange <= -1 ~
        "Down",
      
      TRUE ~
        "Not significant"
    )
  )

stopifnot(
  sum(
    lusc_pcd_de$DE_status == "Up"
  ) == 61
)

stopifnot(
  sum(
    lusc_pcd_de$DE_status == "Down"
  ) == 34
)


# ------------------------------------------------------------
# 35. LUSC significant PCD genes
# ------------------------------------------------------------

lusc_pcd_significant <- lusc_pcd_de |>
  dplyr::filter(
    DE_status != "Not significant"
  ) |>
  dplyr::arrange(
    padj
  )

stopifnot(
  nrow(lusc_pcd_significant) == 95
)


# ------------------------------------------------------------
# 36. LUSC summary by PCD type
# ------------------------------------------------------------

lusc_pcd_de_by_type <- lusc_pcd_de |>
  dplyr::select(
    gene_name,
    log2FoldChange,
    pvalue,
    padj,
    DE_status
  ) |>
  dplyr::left_join(
    pcd_master_long_final |>
      dplyr::distinct(
        gene_symbol,
        PCD_type
      ),
    by = c(
      "gene_name" = "gene_symbol"
    )
  )

lusc_pcd_de_summary <- lusc_pcd_de_by_type |>
  dplyr::count(
    PCD_type,
    DE_status
  ) |>
  tidyr::pivot_wider(
    names_from = DE_status,
    values_from = n,
    values_fill = 0
  )


# ============================================================
# PART E — LUAD vs LUSC COMPARISON
# ============================================================


# ------------------------------------------------------------
# 37. Build direct PCD comparison table
# ------------------------------------------------------------

pcd_de_comparison <- luad_pcd_de |>
  dplyr::select(
    gene_name,
    PCD_types,
    LUAD_log2FC = log2FoldChange,
    LUAD_padj = padj,
    LUAD_status = DE_status
  ) |>
  full_join(
    lusc_pcd_de |>
      dplyr::select(
        gene_name,
        LUSC_log2FC = log2FoldChange,
        LUSC_padj = padj,
        LUSC_status = DE_status
      ),
    by = "gene_name"
  )


# ------------------------------------------------------------
# 38. Classify LUAD/LUSC relationship
# ------------------------------------------------------------

pcd_de_comparison <- pcd_de_comparison |>
  dplyr::mutate(
    comparison_class = case_when(
      
      LUAD_status == "Up" &
        LUSC_status == "Up" ~
        "Shared Up",
      
      LUAD_status == "Down" &
        LUSC_status == "Down" ~
        "Shared Down",
      
      LUAD_status == "Up" &
        LUSC_status == "Down" ~
        "Opposite: LUAD Up / LUSC Down",
      
      LUAD_status == "Down" &
        LUSC_status == "Up" ~
        "Opposite: LUAD Down / LUSC Up",
      
      LUAD_status != "Not significant" &
        LUSC_status == "Not significant" ~
        "LUAD specific",
      
      LUAD_status == "Not significant" &
        LUSC_status != "Not significant" ~
        "LUSC specific",
      
      TRUE ~
        "Not significant in either"
    )
  )


# ------------------------------------------------------------
# 39. Concordant PCD candidate genes
# ------------------------------------------------------------

concordant_pcd_candidates <- pcd_de_comparison |>
  dplyr::filter(
    comparison_class %in% c(
      "Shared Up",
      "Shared Down"
    )
  ) |>
  dplyr::mutate(
    mean_abs_log2FC =
      (
        abs(LUAD_log2FC) +
          abs(LUSC_log2FC)
      ) / 2
  ) |>
  dplyr::arrange(
    comparison_class,
    desc(mean_abs_log2FC)
  )

stopifnot(
  nrow(concordant_pcd_candidates) == 47
)

stopifnot(
  sum(
    concordant_pcd_candidates$comparison_class ==
      "Shared Up"
  ) == 31
)

stopifnot(
  sum(
    concordant_pcd_candidates$comparison_class ==
      "Shared Down"
  ) == 16
)


# ------------------------------------------------------------
# 40. Concordant candidates by PCD type
# ------------------------------------------------------------

concordant_pcd_by_type <-
  concordant_pcd_candidates |>
  dplyr::select(
    gene_name,
    comparison_class
  ) |>
  dplyr::left_join(
    pcd_master_long_final |>
      dplyr::distinct(
        gene_symbol,
        PCD_type
      ),
    by = c(
      "gene_name" = "gene_symbol"
    )
  )

concordant_pcd_type_summary <-
  concordant_pcd_by_type |>
  dplyr::count(
    PCD_type,
    comparison_class
  ) |>
  tidyr::pivot_wider(
    names_from = comparison_class,
    values_from = n,
    values_fill = 0
  )


# ------------------------------------------------------------
# 41. Opposite-direction candidates
# ------------------------------------------------------------

opposite_pcd_candidates <-
  pcd_de_comparison |>
  dplyr::filter(
    grepl(
      "^Opposite",
      comparison_class
    )
  )

stopifnot(
  nrow(opposite_pcd_candidates) == 3
)

stopifnot(
  setequal(
    opposite_pcd_candidates$gene_name,
    c(
      "TFRC",
      "DPP4",
      "IL1A"
    )
  )
)


# ------------------------------------------------------------
# 42. Subtype-specific PCD genes
# ------------------------------------------------------------

luad_specific_pcd <-
  pcd_de_comparison |>
  dplyr::filter(
    comparison_class == "LUAD specific"
  )

lusc_specific_pcd <-
  pcd_de_comparison |>
  dplyr::filter(
    comparison_class == "LUSC specific"
  )

stopifnot(
  nrow(luad_specific_pcd) == 7
)

stopifnot(
  nrow(lusc_specific_pcd) == 45
)


# ============================================================
# PART F — SAVE FINAL RESULTS
# ============================================================


# ------------------------------------------------------------
# LUAD — 6 files
# ------------------------------------------------------------

write.csv(
  luad_de_all,
  file.path(
    de_results_dir,
    "TCGA_LUAD_DESeq2_all_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  luad_pcd_de_summary,
  file.path(
    de_results_dir,
    "TCGA_LUAD_PCD_DE_summary_by_type.csv"
  ),
  row.names = FALSE
)

write.csv(
  luad_pcd_de,
  file.path(
    de_results_dir,
    "TCGA_LUAD_PCD_DESeq2_all_tested.csv"
  ),
  row.names = FALSE
)

write.csv(
  luad_pcd_significant,
  file.path(
    de_results_dir,
    "TCGA_LUAD_PCD_DESeq2_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  luad_pcd_filtered_out,
  file.path(
    de_results_dir,
    "TCGA_LUAD_PCD_genes_filtered_low_expression.csv"
  ),
  row.names = FALSE
)

write.csv(
  luad_specific_pcd,
  file.path(
    de_results_dir,
    "TCGA_LUAD_specific_PCD_genes.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# LUSC — 6 files
# ------------------------------------------------------------

write.csv(
  lusc_de_all,
  file.path(
    de_results_dir,
    "TCGA_LUSC_DESeq2_all_genes.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_pcd_de_summary,
  file.path(
    de_results_dir,
    "TCGA_LUSC_PCD_DE_summary_by_type.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_pcd_de,
  file.path(
    de_results_dir,
    "TCGA_LUSC_PCD_DESeq2_all_tested.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_pcd_significant,
  file.path(
    de_results_dir,
    "TCGA_LUSC_PCD_DESeq2_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_pcd_filtered_out,
  file.path(
    de_results_dir,
    "TCGA_LUSC_PCD_genes_filtered_low_expression.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_specific_pcd,
  file.path(
    de_results_dir,
    "TCGA_LUSC_specific_PCD_genes.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Cross-cohort — 4 files
# ------------------------------------------------------------

write.csv(
  pcd_de_comparison,
  file.path(
    de_results_dir,
    "TCGA_LUAD_vs_LUSC_PCD_DE_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  concordant_pcd_candidates,
  file.path(
    de_results_dir,
    "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv"
  ),
  row.names = FALSE
)

write.csv(
  concordant_pcd_type_summary,
  file.path(
    de_results_dir,
    "TCGA_concordant_PCD_candidates_summary_by_type.csv"
  ),
  row.names = FALSE
)

write.csv(
  opposite_pcd_candidates,
  file.path(
    de_results_dir,
    "TCGA_opposite_direction_PCD_candidates.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART G — FINAL VALIDATION
# ============================================================

expected_de_files <- c(
  "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv",
  "TCGA_concordant_PCD_candidates_summary_by_type.csv",
  "TCGA_LUAD_DESeq2_all_genes.csv",
  "TCGA_LUAD_PCD_DE_summary_by_type.csv",
  "TCGA_LUAD_PCD_DESeq2_all_tested.csv",
  "TCGA_LUAD_PCD_DESeq2_significant.csv",
  "TCGA_LUAD_PCD_genes_filtered_low_expression.csv",
  "TCGA_LUAD_specific_PCD_genes.csv",
  "TCGA_LUAD_vs_LUSC_PCD_DE_comparison.csv",
  "TCGA_LUSC_DESeq2_all_genes.csv",
  "TCGA_LUSC_PCD_DE_summary_by_type.csv",
  "TCGA_LUSC_PCD_DESeq2_all_tested.csv",
  "TCGA_LUSC_PCD_DESeq2_significant.csv",
  "TCGA_LUSC_PCD_genes_filtered_low_expression.csv",
  "TCGA_LUSC_specific_PCD_genes.csv",
  "TCGA_opposite_direction_PCD_candidates.csv"
)

stopifnot(
  all(
    file.exists(
      file.path(
        de_results_dir,
        expected_de_files
      )
    )
  )
)


# ------------------------------------------------------------
# 43. Completion summary
# ------------------------------------------------------------

cat(
  "\n====================================================\n"
)

cat(
  "PCD differential-expression analysis completed.\n"
)

cat(
  "====================================================\n\n"
)

cat(
  "LUAD:\n",
  "  DESeq2 samples: 576\n",
  "  Tumor: 517\n",
  "  Normal: 59\n",
  "  Tested PCD genes: 294\n",
  "  Significant PCD genes: 57\n",
  "  Up: 36\n",
  "  Down: 21\n\n",
  sep = ""
)

cat(
  "LUSC:\n",
  "  DESeq2 samples: 552\n",
  "  Tumor: 501\n",
  "  Normal: 51\n",
  "  Tested PCD genes: 294\n",
  "  Significant PCD genes: 95\n",
  "  Up: 61\n",
  "  Down: 34\n\n",
  sep = ""
)

cat(
  "Cross-cohort:\n",
  "  Concordant genes: 47\n",
  "  Shared Up: 31\n",
  "  Shared Down: 16\n",
  "  Opposite-direction genes: 3\n",
  "  LUAD-specific genes: 7\n",
  "  LUSC-specific genes: 45\n\n",
  sep = ""
)

cat(
  "Final differential-expression files:",
  length(expected_de_files),
  "\n"
)

cat(
  "Results directory:\n",
  de_results_dir,
  "\n"
)

cat(
  "\n========================================\n",
  "SCRIPT 11 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================