# ============================================================
# MSc Lung Cancer PCD Project
# Script: 04_prepare_TCGA_LUSC.R
#
# Purpose:
#   Read canonical TCGA-LUSC STAR-count files, construct raw-count
#   and TPM matrices, create gene annotation and sample metadata,
#   perform basic integrity checks, and save processed datasets.
#
# Inputs:
#   01_raw_data/TCGA_LUSC/RNAseq/
#
# Outputs:
#   02_processed_data/expression/
# ============================================================

source("03_scripts/00_project_config.R")

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

lusc_raw_dir <- file.path(
  raw_dir,
  "TCGA_LUSC",
  "RNAseq"
)

expression_processed_dir <- file.path(
  processed_dir,
  "expression"
)

check_dirs_exist(
  lusc_raw_dir,
  label = "TCGA-LUSC RNA-seq directory"
)

ensure_dir(expression_processed_dir)

lusc_files <- list.files(
  lusc_raw_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat(
  "\nNumber of canonical LUSC STAR-count files:",
  length(lusc_files),
  "\n"
)

if (length(lusc_files) == 0) {
  stop(
    "No TCGA-LUSC STAR-count TSV files were found under:\n",
    lusc_raw_dir
  )
}

if (anyDuplicated(basename(lusc_files)) > 0) {
  stop(
    "Duplicate STAR-count filenames were found in the canonical LUSC raw-data directory."
  )
}

query_lusc_exp <- TCGAbiolinks::GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c(
    "Primary Tumor",
    "Solid Tissue Normal"
  )
)

lusc_query_results <- TCGAbiolinks::getResults(
  query_lusc_exp
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
      colnames(lusc_query_results)
  )
)

lusc_downloaded <- data.frame(
  path = lusc_files,
  file_name = basename(lusc_files),
  stringsAsFactors = FALSE
)

metadata_match <- match(
  lusc_downloaded$file_name,
  lusc_query_results$file_name
)

if (anyNA(metadata_match)) {
  stop(
    "At least one canonical LUSC STAR-count file could not be matched ",
    "to the current GDC query metadata."
  )
}

lusc_metadata_matched <- lusc_query_results[
  metadata_match,
  required_metadata_columns,
  drop = FALSE
]

stopifnot(
  identical(
    lusc_downloaded$file_name,
    lusc_metadata_matched$file_name
  )
)

lusc_downloaded <- cbind(
  lusc_downloaded,
  lusc_metadata_matched[
    ,
    setdiff(
      colnames(lusc_metadata_matched),
      "file_name"
    ),
    drop = FALSE
  ]
)

query_files_missing_from_disk <- setdiff(
  lusc_query_results$file_name,
  lusc_downloaded$file_name
)

disk_files_not_in_query <- setdiff(
  lusc_downloaded$file_name,
  lusc_query_results$file_name
)

if (length(query_files_missing_from_disk) > 0) {
  stop(
    "Canonical LUSC RNA-seq storage is missing ",
    length(query_files_missing_from_disk),
    " file(s) returned by the current GDC query."
  )
}

if (length(disk_files_not_in_query) > 0) {
  stop(
    "Canonical LUSC RNA-seq storage contains ",
    length(disk_files_not_in_query),
    " unexpected STAR-count file(s)."
  )
}

cat("\nMetadata dimensions:\n")
print(dim(lusc_downloaded))

cat("\nMissing sample types:\n")
print(sum(is.na(lusc_downloaded$sample_type)))

cat("\nSample-type distribution:\n")
print(table(lusc_downloaded$sample_type))

cat("\nDuplicated case identifiers:\n")
print(sum(duplicated(lusc_downloaded$cases)))

cat("\nDuplicated biological sample IDs:\n")
print(
  sum(
    duplicated(
      lusc_downloaded$sample.submitter_id
    )
  )
)

stopifnot(
  !anyNA(lusc_downloaded$sample_type),
  !anyNA(lusc_downloaded$cases),
  !anyNA(lusc_downloaded$sample.submitter_id)
)

lusc_test <- data.table::fread(
  lusc_files[1],
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

if (!all(required_star_columns %in% colnames(lusc_test))) {
  stop(
    "STAR-count file is missing required column(s): ",
    paste(
      setdiff(required_star_columns, colnames(lusc_test)),
      collapse = ", "
    )
  )
}

lusc_genes <- lusc_test[
  !grepl("^N_", lusc_test$gene_id),
  ,
  drop = FALSE
]

reference_lusc_genes <- lusc_genes$gene_id

if (anyDuplicated(reference_lusc_genes) > 0) {
  stop(
    "Duplicate gene_id values were found in the reference LUSC STAR-count file."
  )
}

cat(
  "\nNumber of gene rows after removing STAR summary rows:\n"
)
print(nrow(lusc_genes))

lusc_gene_annotation <- lusc_genes[
  ,
  c(
    "gene_id",
    "gene_name",
    "gene_type"
  ),
  drop = FALSE
]

lusc_gene_annotation$ensembl_id <- sub(
  "\\..*$",
  "",
  lusc_gene_annotation$gene_id
)

n_genes <- nrow(lusc_genes)
n_samples <- nrow(lusc_downloaded)

lusc_counts <- matrix(
  NA_real_,
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(
    reference_lusc_genes,
    lusc_downloaded$cases
  )
)

lusc_tpm <- matrix(
  NA_real_,
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(
    reference_lusc_genes,
    lusc_downloaded$cases
  )
)

all_lusc_gene_order_ok <- TRUE

for (i in seq_len(n_samples)) {
  
  x <- data.table::fread(
    lusc_downloaded$path[i],
    data.table = FALSE,
    skip = 1
  )
  
  if (!all(required_star_columns %in% colnames(x))) {
    stop(
      "Required STAR-count columns are missing in file:\n",
      lusc_downloaded$path[i]
    )
  }
  
  x <- x[
    !grepl("^N_", x$gene_id),
    ,
    drop = FALSE
  ]
  
  if (!identical(reference_lusc_genes, x$gene_id)) {
    all_lusc_gene_order_ok <- FALSE
    stop(
      "Gene order mismatch in LUSC STAR-count file:\n",
      lusc_downloaded$path[i]
    )
  }
  
  lusc_counts[, i] <- x$unstranded
  lusc_tpm[, i] <- x$tpm_unstranded
  
  if (i %% 50 == 0 || i == n_samples) {
    message(
      "Processed LUSC expression files: ",
      i,
      " of ",
      n_samples
    )
  }
}

cat("\nRaw-count matrix dimensions:\n")
print(dim(lusc_counts))

cat("\nTPM matrix dimensions:\n")
print(dim(lusc_tpm))

cat("\nMissing raw-count values:\n")
print(sum(is.na(lusc_counts)))

cat("\nMissing TPM values:\n")
print(sum(is.na(lusc_tpm)))

cat("\nNegative raw counts:\n")
print(sum(lusc_counts < 0))

cat("\nNegative TPM values:\n")
print(sum(lusc_tpm < 0))

cat("\nGene order consistent across all files:\n")
print(all_lusc_gene_order_ok)

stopifnot(
  all_lusc_gene_order_ok,
  sum(is.na(lusc_counts)) == 0,
  sum(is.na(lusc_tpm)) == 0,
  sum(lusc_counts < 0) == 0,
  sum(lusc_tpm < 0) == 0,
  identical(
    rownames(lusc_counts),
    lusc_gene_annotation$gene_id
  ),
  identical(
    rownames(lusc_tpm),
    lusc_gene_annotation$gene_id
  ),
  identical(
    colnames(lusc_counts),
    lusc_downloaded$cases
  ),
  identical(
    colnames(lusc_tpm),
    lusc_downloaded$cases
  )
)

saveRDS(
  lusc_counts,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_raw_counts.rds"
  )
)

saveRDS(
  lusc_tpm,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_TPM.rds"
  )
)

saveRDS(
  lusc_gene_annotation,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_gene_annotation.rds"
  )
)

saveRDS(
  lusc_downloaded,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_sample_metadata.rds"
  )
)

write.csv(
  lusc_gene_annotation,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_gene_annotation.csv"
  ),
  row.names = FALSE
)

write.csv(
  lusc_downloaded,
  file.path(
    expression_processed_dir,
    "TCGA_LUSC_sample_metadata.csv"
  ),
  row.names = FALSE
)

cat(
  "\n========================================\n",
  "SCRIPT 04 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "Processed samples: ", ncol(lusc_counts), "\n",
  "Genes: ", nrow(lusc_counts), "\n",
  "Outputs: ", expression_processed_dir, "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================