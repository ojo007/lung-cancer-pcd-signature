# ============================================================
# 10_extract_TCGA_PCD_expression.R
#
# Purpose:
# Match curated PCD genes to TCGA LUAD/LUSC RNA-seq annotation
# and extract PCD-specific raw-count and TPM matrices.
# ============================================================


# -----------------------------
# 0. Shared project configuration
# -----------------------------

source("03_scripts/00_project_config.R")


# -----------------------------
# 1. Project paths
# -----------------------------

expression_processed_dir <- file.path(
  processed_dir,
  "expression"
)

pcd_processed_dir <- file.path(
  processed_dir,
  "PCD_genes"
)


# -----------------------------
# 2. Required package
# -----------------------------

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

suppressPackageStartupMessages(
  library(dplyr)
)


# -----------------------------
# 3. Load curated PCD genes
# -----------------------------

pcd_master_gene <- read.csv(
  file.path(
    pcd_processed_dir,
    "PCD_master_unique_genes.csv"
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


# -----------------------------
# 4. Load gene annotations
# -----------------------------

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

stopifnot(
  nrow(luad_gene_annotation) == 60660,
  nrow(lusc_gene_annotation) == 60660
)

stopifnot(
  identical(
    colnames(luad_gene_annotation),
    colnames(lusc_gene_annotation)
  )
)


# -----------------------------
# 5. Match PCD genes
# -----------------------------

luad_pcd_match <- pcd_master_gene |>
  left_join(
    luad_gene_annotation |>
      select(
        gene_name,
        gene_id,
        gene_type,
        ensembl_id
      ),
    by = c(
      "gene_symbol" = "gene_name"
    )
  )

lusc_pcd_match <- pcd_master_gene |>
  left_join(
    lusc_gene_annotation |>
      select(
        gene_name,
        gene_id,
        gene_type,
        ensembl_id
      ),
    by = c(
      "gene_symbol" = "gene_name"
    )
  )


# -----------------------------
# 6. Validate matching
# -----------------------------

stopifnot(
  sum(
    is.na(
      luad_pcd_match$gene_id
    )
  ) == 0
)

stopifnot(
  sum(
    is.na(
      lusc_pcd_match$gene_id
    )
  ) == 0
)

stopifnot(
  !anyDuplicated(
    luad_pcd_match$gene_symbol
  )
)

stopifnot(
  !anyDuplicated(
    lusc_pcd_match$gene_symbol
  )
)

stopifnot(
  identical(
    luad_pcd_match$ensembl_id,
    lusc_pcd_match$ensembl_id
  )
)


# -----------------------------
# 7. Load RNA-seq matrices
# -----------------------------

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


# -----------------------------
# 8. Validate annotation order
# -----------------------------

stopifnot(
  identical(
    rownames(luad_counts),
    luad_gene_annotation$gene_id
  )
)

stopifnot(
  identical(
    rownames(luad_tpm),
    luad_gene_annotation$gene_id
  )
)

stopifnot(
  identical(
    rownames(lusc_counts),
    lusc_gene_annotation$gene_id
  )
)

stopifnot(
  identical(
    rownames(lusc_tpm),
    lusc_gene_annotation$gene_id
  )
)


# -----------------------------
# 9. Match PCD row indices
# -----------------------------

luad_pcd_index <- match(
  luad_pcd_match$gene_id,
  luad_gene_annotation$gene_id
)

lusc_pcd_index <- match(
  lusc_pcd_match$gene_id,
  lusc_gene_annotation$gene_id
)

stopifnot(
  !anyNA(
    luad_pcd_index
  )
)

stopifnot(
  !anyNA(
    lusc_pcd_index
  )
)


# -----------------------------
# 10. Extract PCD matrices
# -----------------------------

luad_pcd_counts <- luad_counts[
  luad_pcd_index,
  ,
  drop = FALSE
]

luad_pcd_tpm <- luad_tpm[
  luad_pcd_index,
  ,
  drop = FALSE
]

lusc_pcd_counts <- lusc_counts[
  lusc_pcd_index,
  ,
  drop = FALSE
]

lusc_pcd_tpm <- lusc_tpm[
  lusc_pcd_index,
  ,
  drop = FALSE
]


# -----------------------------
# 10A. Validate extracted matrix dimensions and sample alignment
# -----------------------------

stopifnot(
  nrow(luad_pcd_counts) == 296,
  nrow(luad_pcd_tpm) == 296,
  nrow(lusc_pcd_counts) == 296,
  nrow(lusc_pcd_tpm) == 296,
  identical(
    colnames(luad_pcd_counts),
    colnames(luad_pcd_tpm)
  ),
  identical(
    colnames(lusc_pcd_counts),
    colnames(lusc_pcd_tpm)
  )
)


# -----------------------------
# 11. Replace row names
# with gene symbols
# -----------------------------

rownames(
  luad_pcd_counts
) <- luad_pcd_match$gene_symbol

rownames(
  luad_pcd_tpm
) <- luad_pcd_match$gene_symbol

rownames(
  lusc_pcd_counts
) <- lusc_pcd_match$gene_symbol

rownames(
  lusc_pcd_tpm
) <- lusc_pcd_match$gene_symbol


# -----------------------------
# 12. Final QC
# -----------------------------

stopifnot(
  nrow(luad_pcd_counts) == 296,
  nrow(lusc_pcd_counts) == 296
)

stopifnot(
  !anyNA(
    luad_pcd_counts
  ),
  !anyNA(
    luad_pcd_tpm
  ),
  !anyNA(
    lusc_pcd_counts
  ),
  !anyNA(
    lusc_pcd_tpm
  )
)

stopifnot(
  identical(
    rownames(luad_pcd_counts),
    rownames(lusc_pcd_counts)
  )
)


# -----------------------------
# 13. Save final outputs
# -----------------------------

saveRDS(
  luad_pcd_counts,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_PCD_raw_counts.rds"
  )
)

saveRDS(
  luad_pcd_tpm,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_PCD_TPM.rds"
  )
)

saveRDS(
  lusc_pcd_counts,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_PCD_raw_counts.rds"
  )
)

saveRDS(
  lusc_pcd_tpm,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_PCD_TPM.rds"
  )
)

write.csv(
  luad_pcd_match,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_PCD_gene_mapping.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_pcd_match,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_PCD_gene_mapping.csv"
  ),
  row.names = FALSE
)


cat(
  "
========================================
",
  "SCRIPT 10 COMPLETED SUCCESSFULLY
",
  "========================================
",
  sep = ""
)

cat(
  "LUAD PCD matrix: ",
  nrow(luad_pcd_counts),
  " genes x ",
  ncol(luad_pcd_counts),
  " samples
",
  sep = ""
)

cat(
  "LUSC PCD matrix: ",
  nrow(lusc_pcd_counts),
  " genes x ",
  ncol(lusc_pcd_counts),
  " samples
",
  sep = ""
)

cat(
  "All 296 curated PCD genes mapped in both histologies.
"
)

# ============================================================
# End of script
# ============================================================