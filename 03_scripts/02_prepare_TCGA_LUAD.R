# ============================================================
# MSc Lung Cancer PCD Project
# Script: 02_prepare_TCGA_LUAD.R
#
# Purpose:
#   Read canonical TCGA-LUAD STAR-count files, construct raw-count
#   and TPM matrices, create gene annotation and sample metadata,
#   perform basic integrity checks, and save processed datasets.
#
# Inputs:
#   01_raw_data/TCGA_LUAD/RNAseq/
#
# Outputs:
#   02_processed_data/expression/
# ============================================================

# ------------------------------------------------------------
# 0. Load shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")

# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "TCGAbiolinks",
  "data.table"
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

# ------------------------------------------------------------
# 2. Canonical directories
# ------------------------------------------------------------

luad_raw_dir <- file.path(
  raw_dir,
  "TCGA_LUAD",
  "RNAseq"
)

expression_processed_dir <- file.path(
  processed_dir,
  "expression"
)

check_dirs_exist(
  luad_raw_dir,
  label = "TCGA-LUAD RNA-seq directory"
)

ensure_dir(expression_processed_dir)

# ------------------------------------------------------------
# 3. Locate canonical STAR-count files
# ------------------------------------------------------------

luad_files <- list.files(
  luad_raw_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat(
  "\nNumber of canonical LUAD STAR-count files:",
  length(luad_files),
  "\n"
)

if (length(luad_files) == 0) {
  stop(
    "No TCGA-LUAD STAR-count TSV files were found under:\n",
    luad_raw_dir
  )
}

if (anyDuplicated(basename(luad_files)) > 0) {
  stop(
    "Duplicate STAR-count filenames were found in the canonical LUAD raw-data directory."
  )
}

# ------------------------------------------------------------
# 4. Recreate GDC query to recover authoritative metadata
# ------------------------------------------------------------

query_luad_exp <- TCGAbiolinks::GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c(
    "Primary Tumor",
    "Solid Tissue Normal"
  )
)

luad_query_results <- TCGAbiolinks::getResults(
  query_luad_exp
)

required_metadata_columns <- c(
  "file_id",
  "file_name",
  "cases",
  "sample_type",
  "cases.submitter_id",
  "sample.submitter_id"
)

stopifnot(
  all(
    required_metadata_columns %in%
      colnames(luad_query_results)
  )
)

# ------------------------------------------------------------
# 5. Match canonical files to GDC metadata without reordering
# ------------------------------------------------------------

luad_downloaded <- data.frame(
  path = luad_files,
  file_name = basename(luad_files),
  stringsAsFactors = FALSE
)

metadata_match <- match(
  luad_downloaded$file_name,
  luad_query_results$file_name
)

if (anyNA(metadata_match)) {
  stop(
    "At least one canonical LUAD STAR-count file could not be matched ",
    "to the current GDC query metadata."
  )
}

luad_metadata_matched <- luad_query_results[
  metadata_match,
  required_metadata_columns,
  drop = FALSE
]

# Validate filename alignment before attaching the remaining metadata.
stopifnot(
  identical(
    luad_downloaded$file_name,
    luad_metadata_matched$file_name
  )
)

# Attach metadata columns without duplicating file_name.
luad_downloaded <- cbind(
  luad_downloaded,
  luad_metadata_matched[
    ,
    setdiff(
      colnames(luad_metadata_matched),
      "file_name"
    ),
    drop = FALSE
  ]
)

# Confirm complete query coverage in canonical storage
query_files_missing_from_disk <- setdiff(
  luad_query_results$file_name,
  luad_downloaded$file_name
)

disk_files_not_in_query <- setdiff(
  luad_downloaded$file_name,
  luad_query_results$file_name
)

if (length(query_files_missing_from_disk) > 0) {
  stop(
    "Canonical LUAD RNA-seq storage is missing ",
    length(query_files_missing_from_disk),
    " file(s) returned by the current GDC query."
  )
}

if (length(disk_files_not_in_query) > 0) {
  stop(
    "Canonical LUAD RNA-seq storage contains ",
    length(disk_files_not_in_query),
    " unexpected STAR-count file(s)."
  )
}

# ------------------------------------------------------------
# 6. Metadata QC
# ------------------------------------------------------------

cat("\nMetadata dimensions:\n")
print(dim(luad_downloaded))

cat("\nMissing sample types:\n")
print(sum(is.na(luad_downloaded$sample_type)))

cat("\nSample-type distribution:\n")
print(table(luad_downloaded$sample_type))

cat("\nDuplicated case identifiers:\n")
print(sum(duplicated(luad_downloaded$cases)))

cat("\nDuplicated biological sample IDs:\n")
print(
  sum(
    duplicated(
      luad_downloaded$sample.submitter_id
    )
  )
)

stopifnot(
  !anyNA(luad_downloaded$sample_type),
  !anyNA(luad_downloaded$cases),
  !anyNA(luad_downloaded$sample.submitter_id)
)

# ------------------------------------------------------------
# 7. Read first STAR file and establish reference gene order
# ------------------------------------------------------------

luad_test <- data.table::fread(
  luad_files[1],
  data.table = FALSE,
  skip = 1
)

required_star_columns <- c(
  "gene_id",
  "gene_name",
  "gene_type",
  "unstranded",
  "tpm_unstranded"
)

if (!all(required_star_columns %in% colnames(luad_test))) {
  stop(
    "STAR-count file is missing required column(s): ",
    paste(
      setdiff(required_star_columns, colnames(luad_test)),
      collapse = ", "
    )
  )
}

luad_genes <- luad_test[
  !grepl("^N_", luad_test$gene_id),
  ,
  drop = FALSE
]

reference_luad_genes <- luad_genes$gene_id

if (anyDuplicated(reference_luad_genes) > 0) {
  stop("Duplicate gene_id values were found in the reference LUAD STAR-count file.")
}

cat("\nNumber of gene rows after removing STAR summary rows:\n")
print(nrow(luad_genes))

# ------------------------------------------------------------
# 8. Gene annotation
# ------------------------------------------------------------

luad_gene_annotation <- luad_genes[
  ,
  c(
    "gene_id",
    "gene_name",
    "gene_type"
  ),
  drop = FALSE
]

luad_gene_annotation$ensembl_id <- sub(
  "\\..*$",
  "",
  luad_gene_annotation$gene_id
)

# ------------------------------------------------------------
# 9. Allocate matrices
#
# Counts and TPM are populated in a single pass so every raw
# STAR-count file is read only once.
# ------------------------------------------------------------

n_genes <- nrow(luad_genes)
n_samples <- nrow(luad_downloaded)

luad_counts <- matrix(
  NA_real_,
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(
    reference_luad_genes,
    luad_downloaded$cases
  )
)

luad_tpm <- matrix(
  NA_real_,
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(
    reference_luad_genes,
    luad_downloaded$cases
  )
)

# ------------------------------------------------------------
# 10. Read all samples once and validate gene order
# ------------------------------------------------------------

all_luad_gene_order_ok <- TRUE

for (i in seq_len(n_samples)) {
  
  x <- data.table::fread(
    luad_downloaded$path[i],
    data.table = FALSE,
    skip = 1
  )
  
  if (!all(required_star_columns %in% colnames(x))) {
    stop(
      "Required STAR-count columns are missing in file:\n",
      luad_downloaded$path[i]
    )
  }
  
  x <- x[
    !grepl("^N_", x$gene_id),
    ,
    drop = FALSE
  ]
  
  if (!identical(reference_luad_genes, x$gene_id)) {
    all_luad_gene_order_ok <- FALSE
    stop(
      "Gene order mismatch in LUAD STAR-count file:\n",
      luad_downloaded$path[i]
    )
  }
  
  luad_counts[, i] <- x$unstranded
  luad_tpm[, i] <- x$tpm_unstranded
  
  if (i %% 50 == 0 || i == n_samples) {
    message(
      "Processed LUAD expression files: ",
      i,
      " of ",
      n_samples
    )
  }
}

# ------------------------------------------------------------
# 11. Final QC
# ------------------------------------------------------------

cat("\nRaw-count matrix dimensions:\n")
print(dim(luad_counts))

cat("\nTPM matrix dimensions:\n")
print(dim(luad_tpm))

cat("\nMissing raw-count values:\n")
print(sum(is.na(luad_counts)))

cat("\nMissing TPM values:\n")
print(sum(is.na(luad_tpm)))

cat("\nNegative raw counts:\n")
print(sum(luad_counts < 0))

cat("\nNegative TPM values:\n")
print(sum(luad_tpm < 0))

cat("\nGene order consistent across all files:\n")
print(all_luad_gene_order_ok)

stopifnot(
  all_luad_gene_order_ok,
  sum(is.na(luad_counts)) == 0,
  sum(is.na(luad_tpm)) == 0,
  sum(luad_counts < 0) == 0,
  sum(luad_tpm < 0) == 0,
  identical(
    rownames(luad_counts),
    luad_gene_annotation$gene_id
  ),
  identical(
    rownames(luad_tpm),
    luad_gene_annotation$gene_id
  ),
  identical(
    colnames(luad_counts),
    luad_downloaded$cases
  ),
  identical(
    colnames(luad_tpm),
    luad_downloaded$cases
  )
)

# ------------------------------------------------------------
# 12. Save processed datasets
# ------------------------------------------------------------

saveRDS(
  luad_counts,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_raw_counts.rds"
  )
)

saveRDS(
  luad_tpm,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_TPM.rds"
  )
)

saveRDS(
  luad_gene_annotation,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_gene_annotation.rds"
  )
)

saveRDS(
  luad_downloaded,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_sample_metadata.rds"
  )
)

write.csv(
  luad_gene_annotation,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_gene_annotation.csv"
  ),
  row.names = FALSE
)

write.csv(
  luad_downloaded,
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_sample_metadata.csv"
  ),
  row.names = FALSE
)

cat(
  "\n========================================\n",
  "SCRIPT 02 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "Processed samples: ", ncol(luad_counts), "\n",
  "Genes: ", nrow(luad_counts), "\n",
  "Outputs: ", expression_processed_dir, "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================