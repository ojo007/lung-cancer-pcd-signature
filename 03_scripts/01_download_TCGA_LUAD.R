# ============================================================
# MSc Lung Cancer PCD Project
# Script: 01_download_TCGA_LUAD.R
#
# Purpose:
#   Query TCGA-LUAD STAR-count RNA-seq data, create the GDC
#   manifest, optionally download through a short Windows staging
#   path, and synchronize validated files into the canonical
#   project raw-data directory.
# ============================================================

source("03_scripts/00_project_config.R")

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  stop("Package 'TCGAbiolinks' is required.")
}

luad_raw_base <- file.path(raw_dir, "TCGA_LUAD")
luad_rnaseq_dir <- file.path(luad_raw_base, "RNAseq")

ensure_dir(luad_raw_base)
ensure_dir(luad_rnaseq_dir)

query_luad_exp <- TCGAbiolinks::GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

luad_query_results <- TCGAbiolinks::getResults(query_luad_exp)

cat("\nNumber of files returned:\n")
print(dim(luad_query_results))
cat("\nSample types:\n")
print(table(luad_query_results$sample_type))

stopifnot(
  nrow(luad_query_results) > 0,
  all(c("file_id", "file_name", "md5sum", "file_size", "sample_type") %in%
        colnames(luad_query_results))
)

luad_manifest <- data.frame(
  id = luad_query_results$file_id,
  filename = luad_query_results$file_name,
  md5 = luad_query_results$md5sum,
  size = luad_query_results$file_size,
  state = "validated",
  stringsAsFactors = FALSE
)

luad_manifest_file <- file.path(
  luad_raw_base,
  "gdc_manifest_TCGA_LUAD_RNAseq.txt"
)

write.table(
  luad_manifest,
  file = luad_manifest_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nManifest saved to:\n", luad_manifest_file, "\n")

RUN_DOWNLOAD <- FALSE

if (RUN_DOWNLOAD) {
  check_files_exist(gdc_client, label = "GDC client")
  ensure_dir(gdc_luad_rnaseq_stage)
  
  download_status <- system2(
    command = gdc_client,
    args = c(
      "download",
      "-m", shQuote(luad_manifest_file),
      "-d", shQuote(gdc_luad_rnaseq_stage)
    )
  )
  
  if (!identical(download_status, 0L)) {
    stop("gdc-client returned a non-zero exit status: ", download_status)
  }
}

SYNC_STAGING_TO_PROJECT <- TRUE

if (SYNC_STAGING_TO_PROJECT && dir.exists(gdc_luad_rnaseq_stage)) {
  
  staged_tsv <- list.files(
    gdc_luad_rnaseq_stage,
    pattern = "\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  staged_partial <- list.files(
    gdc_luad_rnaseq_stage,
    pattern = "\\.partial$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  cat("\nStaged TSV files:", length(staged_tsv), "\n")
  cat("Staged partial files:", length(staged_partial), "\n")
  
  if (length(staged_partial) > 0) {
    warning("Partial GDC downloads were detected in the staging directory.")
  }
  
  stage_map <- data.frame(
    path = staged_tsv,
    filename = basename(staged_tsv),
    stringsAsFactors = FALSE
  )
  
  stage_map <- merge(
    luad_manifest[, c("id", "filename")],
    stage_map,
    by = "filename",
    all.x = TRUE,
    sort = FALSE
  )
  
  missing_from_stage <- stage_map[is.na(stage_map$path), , drop = FALSE]
  
  if (RUN_DOWNLOAD && nrow(missing_from_stage) > 0) {
    stop(
      "The download completed, but ",
      nrow(missing_from_stage),
      " manifest file(s) were not found in staging."
    )
  }
  
  available_stage <- stage_map[!is.na(stage_map$path), , drop = FALSE]
  
  if (nrow(available_stage) > 0) {
    for (i in seq_len(nrow(available_stage))) {
      
      destination_dir <- file.path(
        luad_rnaseq_dir,
        available_stage$id[i]
      )
      ensure_dir(destination_dir)
      
      destination_file <- file.path(
        destination_dir,
        available_stage$filename[i]
      )
      
      if (!file.exists(destination_file)) {
        copied <- file.copy(
          from = available_stage$path[i],
          to = destination_file,
          overwrite = FALSE,
          copy.mode = TRUE,
          copy.date = TRUE
        )
        
        if (!copied) {
          stop(
            "Failed to copy staged file:\n",
            available_stage$path[i],
            "\nto:\n",
            destination_file
          )
        }
      }
    }
  }
}

canonical_tsv <- list.files(
  luad_rnaseq_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

canonical_partial <- list.files(
  luad_rnaseq_dir,
  pattern = "\\.partial$",
  recursive = TRUE,
  full.names = TRUE
)

canonical_names <- basename(canonical_tsv)
manifest_present <- luad_manifest$filename %in% canonical_names

cat("\nCanonical LUAD RNA-seq TSV files:", length(canonical_tsv), "\n")
cat("Canonical partial files:", length(canonical_partial), "\n")
cat(
  "Manifest files present in canonical storage:",
  sum(manifest_present), "/", nrow(luad_manifest), "\n"
)

cat(
  "Total canonical RNA-seq size (GiB):",
  sum(file.info(canonical_tsv)$size, na.rm = TRUE) / 1024^3,
  "\n"
)

if (length(canonical_tsv) > 0) {
  if (length(canonical_partial) > 0) {
    stop("Partial files exist in the canonical LUAD RNA-seq directory.")
  }
  
  missing_manifest_files <- luad_manifest$filename[!manifest_present]
  
  if (length(missing_manifest_files) > 0) {
    stop(
      "Canonical LUAD RNA-seq storage is incomplete. Missing ",
      length(missing_manifest_files),
      " manifest file(s)."
    )
  }
}

cat(
  "\n========================================\n",
  "SCRIPT 01 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)

# ============================================================
# End of script
# ============================================================