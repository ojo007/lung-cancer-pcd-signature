# ============================================================
# MSc Lung Cancer PCD Project
# Script: 03_download_TCGA_LUSC.R
#
# Purpose:
#   Query TCGA-LUSC STAR-count RNA-seq data, create the GDC
#   manifest, optionally download through a short Windows staging
#   path, and synchronize validated files into the canonical
#   project raw-data directory.
# ============================================================

source("03_scripts/00_project_config.R")

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  stop("Package 'TCGAbiolinks' is required.")
}

lusc_raw_base <- file.path(raw_dir, "TCGA_LUSC")
lusc_rnaseq_dir <- file.path(lusc_raw_base, "RNAseq")

ensure_dir(lusc_raw_base)
ensure_dir(lusc_rnaseq_dir)

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

cat("\nNumber of files returned:\n")
print(dim(lusc_query_results))

cat("\nSample types:\n")
print(table(lusc_query_results$sample_type))

stopifnot(
  nrow(lusc_query_results) > 0,
  all(
    c(
      "file_id",
      "file_name",
      "md5sum",
      "file_size",
      "sample_type"
    ) %in% colnames(lusc_query_results)
  )
)

lusc_manifest <- data.frame(
  id = lusc_query_results$file_id,
  filename = lusc_query_results$file_name,
  md5 = lusc_query_results$md5sum,
  size = lusc_query_results$file_size,
  state = "validated",
  stringsAsFactors = FALSE
)

lusc_manifest_file <- file.path(
  lusc_raw_base,
  "gdc_manifest_TCGA_LUSC_RNAseq.txt"
)

write.table(
  lusc_manifest,
  file = lusc_manifest_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nManifest saved to:\n",
  lusc_manifest_file,
  "\n"
)

RUN_DOWNLOAD <- FALSE

if (RUN_DOWNLOAD) {
  check_files_exist(
    gdc_client,
    label = "GDC client"
  )
  
  ensure_dir(gdc_lusc_rnaseq_stage)
  
  download_status <- system2(
    command = gdc_client,
    args = c(
      "download",
      "-m",
      shQuote(lusc_manifest_file),
      "-d",
      shQuote(gdc_lusc_rnaseq_stage)
    )
  )
  
  if (!identical(download_status, 0L)) {
    stop(
      "gdc-client returned a non-zero exit status: ",
      download_status
    )
  }
}

SYNC_STAGING_TO_PROJECT <- TRUE

if (
  SYNC_STAGING_TO_PROJECT &&
  dir.exists(gdc_lusc_rnaseq_stage)
) {
  
  staged_tsv <- list.files(
    gdc_lusc_rnaseq_stage,
    pattern = "\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  staged_partial <- list.files(
    gdc_lusc_rnaseq_stage,
    pattern = "\\.partial$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  cat(
    "\nStaged TSV files:",
    length(staged_tsv),
    "\n"
  )
  
  cat(
    "Staged partial files:",
    length(staged_partial),
    "\n"
  )
  
  if (length(staged_partial) > 0) {
    warning(
      "Partial GDC downloads were detected in the staging directory."
    )
  }
  
  stage_map <- data.frame(
    path = staged_tsv,
    filename = basename(staged_tsv),
    stringsAsFactors = FALSE
  )
  
  stage_map <- merge(
    lusc_manifest[, c("id", "filename")],
    stage_map,
    by = "filename",
    all.x = TRUE,
    sort = FALSE
  )
  
  missing_from_stage <- stage_map[
    is.na(stage_map$path),
    ,
    drop = FALSE
  ]
  
  if (
    RUN_DOWNLOAD &&
    nrow(missing_from_stage) > 0
  ) {
    stop(
      "The download completed, but ",
      nrow(missing_from_stage),
      " manifest file(s) were not found in staging."
    )
  }
  
  available_stage <- stage_map[
    !is.na(stage_map$path),
    ,
    drop = FALSE
  ]
  
  if (nrow(available_stage) > 0) {
    
    for (i in seq_len(nrow(available_stage))) {
      
      destination_dir <- file.path(
        lusc_rnaseq_dir,
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
  lusc_rnaseq_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

canonical_partial <- list.files(
  lusc_rnaseq_dir,
  pattern = "\\.partial$",
  recursive = TRUE,
  full.names = TRUE
)

canonical_names <- basename(canonical_tsv)

manifest_present <- lusc_manifest$filename %in% canonical_names

cat(
  "\nCanonical LUSC RNA-seq TSV files:",
  length(canonical_tsv),
  "\n"
)

cat(
  "Canonical partial files:",
  length(canonical_partial),
  "\n"
)

cat(
  "Manifest files present in canonical storage:",
  sum(manifest_present),
  "/",
  nrow(lusc_manifest),
  "\n"
)

cat(
  "Total canonical RNA-seq size (GiB):",
  sum(
    file.info(canonical_tsv)$size,
    na.rm = TRUE
  ) / 1024^3,
  "\n"
)

if (length(canonical_tsv) > 0) {
  
  if (length(canonical_partial) > 0) {
    stop(
      "Partial files exist in the canonical LUSC RNA-seq directory."
    )
  }
  
  missing_manifest_files <- lusc_manifest$filename[
    !manifest_present
  ]
  
  if (length(missing_manifest_files) > 0) {
    stop(
      "Canonical LUSC RNA-seq storage is incomplete. Missing ",
      length(missing_manifest_files),
      " manifest file(s)."
    )
  }
}

cat(
  "\n========================================\n",
  "SCRIPT 03 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)

# ============================================================
# End of script
# ============================================================