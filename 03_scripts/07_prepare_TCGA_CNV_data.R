# ============================================================
# MSc Lung Cancer PCD Project
# Script: 07_prepare_TCGA_CNV_data.R
#
# Purpose:
#   Query, organize, combine, and prepare TCGA-LUAD and
#   TCGA-LUSC Masked Copy Number Segment data.
#
# Important:
#   - Only Primary Tumor samples are retained.
#   - Duplicate patients/samples are retained here.
#   - Final sample resolution occurs later during multi-omics
#     harmonization.
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

cnv_processed_dir <- file.path(
  processed_dir,
  "cnv"
)

ensure_dir(cnv_processed_dir)

sync_gdc_txt_stage <- function(
    manifest,
    staging_dir,
    canonical_dir
) {
  
  if (!dir.exists(staging_dir)) {
    return(invisible(NULL))
  }
  
  staged_files <- list.files(
    staging_dir,
    pattern = "\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  staged_partial <- list.files(
    staging_dir,
    pattern = "\\.partial$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  cat("\nStaged CNV files:", length(staged_files), "\n")
  cat("Staged partial files:", length(staged_partial), "\n")
  
  if (length(staged_partial) > 0) {
    warning("Partial CNV downloads were detected in staging.")
  }
  
  if (length(staged_files) == 0) {
    return(invisible(NULL))
  }
  
  stage_map <- data.frame(
    path = staged_files,
    filename = basename(staged_files),
    stringsAsFactors = FALSE
  )
  
  stage_map <- merge(
    manifest[, c("id", "filename")],
    stage_map,
    by = "filename",
    all.x = TRUE,
    sort = FALSE
  )
  
  available_stage <- stage_map[
    !is.na(stage_map$path),
    ,
    drop = FALSE
  ]
  
  for (i in seq_len(nrow(available_stage))) {
    
    destination_dir <- file.path(
      canonical_dir,
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
          "Failed to copy staged CNV file:\n",
          available_stage$path[i]
        )
      }
    }
  }
  
  invisible(NULL)
}

prepare_tcga_cnv <- function(
    project,
    cancer_name,
    raw_dir_cancer,
    staging_dir,
    run_download = FALSE,
    sync_staging = TRUE
) {
  
  cat(
    "\n========================================\n",
    "Preparing CNV data for ",
    project,
    "\n",
    "========================================\n",
    sep = ""
  )
  
  # ----------------------------------------------------------
  # 1. Query GDC and retain Primary Tumor
  # ----------------------------------------------------------
  
  query_cnv <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Copy Number Variation",
    data.type = "Masked Copy Number Segment",
    access = "open"
  )
  
  cnv_results <- TCGAbiolinks::getResults(
    query_cnv
  )
  
  stopifnot(
    nrow(cnv_results) > 0,
    all(
      c(
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "sample_type",
        "cases.submitter_id",
        "sample.submitter_id",
        "cases"
      ) %in% colnames(cnv_results)
    )
  )
  
  cnv_tumor_results <- cnv_results[
    cnv_results$sample_type == "Primary Tumor",
    ,
    drop = FALSE
  ]
  
  cat("\nPrimary Tumor CNV query dimensions:\n")
  print(dim(cnv_tumor_results))
  
  cat(
    "\nUnique patients: ",
    length(
      unique(
        cnv_tumor_results$cases.submitter_id
      )
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Unique biological samples: ",
    length(
      unique(
        cnv_tumor_results$sample.submitter_id
      )
    ),
    "\n",
    sep = ""
  )
  
  # ----------------------------------------------------------
  # 2. Raw directory and manifest
  # ----------------------------------------------------------
  
  cnv_raw_dir <- file.path(
    raw_dir,
    raw_dir_cancer,
    "cnv"
  )
  
  ensure_dir(cnv_raw_dir)
  
  cnv_manifest <- data.frame(
    id = cnv_tumor_results$file_id,
    filename = cnv_tumor_results$file_name,
    md5 = cnv_tumor_results$md5sum,
    size = cnv_tumor_results$file_size,
    state = "validated",
    stringsAsFactors = FALSE
  )
  
  manifest_file <- file.path(
    cnv_raw_dir,
    paste0(
      "gdc_manifest_",
      cancer_name,
      "_CNV.txt"
    )
  )
  
  write.table(
    cnv_manifest,
    file = manifest_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  cat(
    "\nExpected CNV download size (MiB): ",
    sum(
      cnv_tumor_results$file_size,
      na.rm = TRUE
    ) / 1024^2,
    "\n",
    sep = ""
  )
  
  # ----------------------------------------------------------
  # 3. Optional GDC download
  # ----------------------------------------------------------
  
  if (run_download) {
    
    check_files_exist(
      gdc_client,
      label = "GDC client"
    )
    
    ensure_dir(staging_dir)
    
    download_status <- system2(
      command = gdc_client,
      args = c(
        "download",
        "-m",
        shQuote(manifest_file),
        "-d",
        shQuote(staging_dir)
      )
    )
    
    if (!identical(download_status, 0L)) {
      stop(
        "gdc-client returned a non-zero exit status for ",
        project,
        ": ",
        download_status
      )
    }
  }
  
  # ----------------------------------------------------------
  # 4. Synchronize staging -> canonical raw storage
  # ----------------------------------------------------------
  
  if (sync_staging) {
    sync_gdc_txt_stage(
      manifest = cnv_manifest,
      staging_dir = staging_dir,
      canonical_dir = cnv_raw_dir
    )
  }
  
  # ----------------------------------------------------------
  # 5. Locate exact canonical manifest files
  # ----------------------------------------------------------
  
  cnv_raw_txt <- list.files(
    cnv_raw_dir,
    pattern = "\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  cnv_files <- cnv_raw_txt[
    basename(cnv_raw_txt) %in%
      cnv_manifest$filename
  ]
  
  unexpected_files <- cnv_raw_txt[
    !basename(cnv_raw_txt) %in%
      cnv_manifest$filename
  ]
  
  missing_manifest_files <- setdiff(
    cnv_manifest$filename,
    basename(cnv_files)
  )
  
  cat(
    "\nCanonical manifest-matching CNV files: ",
    length(cnv_files),
    "\n",
    sep = ""
  )
  
  cat(
    "Unexpected canonical CNV files: ",
    length(unexpected_files),
    "\n",
    sep = ""
  )
  
  cat(
    "Manifest files missing from canonical storage: ",
    length(missing_manifest_files),
    "\n",
    sep = ""
  )
  
  if (length(unexpected_files) > 0) {
    warning(
      length(unexpected_files),
      " unexpected CNV file(s) exist in ",
      cnv_raw_dir,
      ". They will not be used."
    )
  }
  
  if (length(missing_manifest_files) > 0) {
    stop(
      "Canonical CNV storage for ",
      project,
      " is incomplete. Missing ",
      length(missing_manifest_files),
      " manifest file(s)."
    )
  }
  
  cnv_match <- match(
    cnv_manifest$filename,
    basename(cnv_files)
  )
  
  stopifnot(!anyNA(cnv_match))
  
  cnv_files <- cnv_files[
    cnv_match
  ]
  
  # ----------------------------------------------------------
  # 6. Read and combine segment files
  # ----------------------------------------------------------
  
  cnv_list <- vector(
    "list",
    length(cnv_files)
  )
  
  for (i in seq_along(cnv_files)) {
    
    x <- data.table::fread(
      cnv_files[i],
      data.table = FALSE
    )
    
    if (nrow(x) == 0) {
      stop(
        "Unexpected empty CNV file:\n",
        cnv_files[i]
      )
    }
    
    required_segment_columns <- c(
      "Chromosome",
      "Start",
      "End",
      "Num_Probes",
      "Segment_Mean",
      "GDC_Aliquot"
    )
    
    if (!all(required_segment_columns %in% colnames(x))) {
      stop(
        "Unexpected CNV file structure in:\n",
        cnv_files[i]
      )
    }
    
    x$source_file <- basename(
      cnv_files[i]
    )
    
    cnv_list[[i]] <- x
    
    if (i %% 50 == 0 || i == length(cnv_files)) {
      message(
        "Processed ",
        cancer_name,
        " CNV files: ",
        i,
        " of ",
        length(cnv_files)
      )
    }
  }
  
  cnv_combined <- data.table::rbindlist(
    cnv_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  cnv_combined <- as.data.frame(
    cnv_combined
  )
  
  # ----------------------------------------------------------
  # 7. Attach authoritative GDC metadata
  # ----------------------------------------------------------
  
  cnv_metadata <- cnv_tumor_results[
    ,
    c(
      "file_name",
      "cases.submitter_id",
      "sample.submitter_id",
      "cases",
      "sample_type"
    ),
    drop = FALSE
  ]
  
  colnames(cnv_metadata)[
    colnames(cnv_metadata) ==
      "cases.submitter_id"
  ] <- "patient_id"
  
  if (anyDuplicated(cnv_metadata$file_name) > 0) {
    stop(
      project,
      " CNV query contains duplicated file_name values."
    )
  }
  
  metadata_index <- match(
    cnv_combined$source_file,
    cnv_metadata$file_name
  )
  
  if (anyNA(metadata_index)) {
    stop(
      "At least one combined CNV row could not be matched ",
      "to GDC metadata."
    )
  }
  
  metadata_for_rows <- cnv_metadata[
    metadata_index,
    setdiff(
      colnames(cnv_metadata),
      "file_name"
    ),
    drop = FALSE
  ]
  
  cnv_combined <- cbind(
    cnv_combined,
    metadata_for_rows
  )
  
  stopifnot(
    !anyNA(cnv_combined$patient_id),
    all(cnv_combined$sample_type == "Primary Tumor")
  )
  
  # ----------------------------------------------------------
  # 8. Sample summary
  # ----------------------------------------------------------
  
  cnv_sample_summary <- unique(
    cnv_combined[
      ,
      c(
        "source_file",
        "patient_id",
        "sample.submitter_id",
        "cases",
        "sample_type",
        "GDC_Aliquot"
      ),
      drop = FALSE
    ]
  )
  
  segment_counts <- as.data.frame(
    table(
      cnv_combined$source_file
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(segment_counts) <- c(
    "source_file",
    "Segment_Count"
  )
  
  cnv_sample_summary <- merge(
    cnv_sample_summary,
    segment_counts,
    by = "source_file",
    all.x = TRUE,
    sort = FALSE
  )
  
  # ----------------------------------------------------------
  # 9. QC
  # ----------------------------------------------------------
  
  cat("\nCombined CNV dimensions:\n")
  print(dim(cnv_combined))
  
  cat(
    "\nUnique patients in combined CNV data: ",
    length(
      unique(
        cnv_combined$patient_id
      )
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Unique biological samples in combined CNV data: ",
    length(
      unique(
        cnv_combined$sample.submitter_id
      )
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Patients with multiple CNV samples: ",
    sum(
      table(
        cnv_sample_summary$patient_id
      ) > 1
    ),
    "\n",
    sep = ""
  )
  
  stopifnot(
    nrow(cnv_combined) > 0,
    !anyNA(cnv_combined$GDC_Aliquot),
    !anyNA(cnv_combined$Segment_Mean),
    all(cnv_combined$End >= cnv_combined$Start)
  )
  
  # ----------------------------------------------------------
  # 10. Save outputs
  # ----------------------------------------------------------
  
  saveRDS(
    cnv_combined,
    file.path(
      cnv_processed_dir,
      paste0(
        cancer_name,
        "_CNV_raw_combined.rds"
      )
    )
  )
  
  write.csv(
    cnv_metadata,
    file.path(
      cnv_processed_dir,
      paste0(
        cancer_name,
        "_CNV_sample_metadata.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    cnv_sample_summary,
    file.path(
      cnv_processed_dir,
      paste0(
        cancer_name,
        "_CNV_sample_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  invisible(
    list(
      manifest = cnv_manifest,
      combined = cnv_combined,
      metadata = cnv_metadata,
      sample_summary = cnv_sample_summary
    )
  )
}

luad_cnv_objects <- prepare_tcga_cnv(
  project = "TCGA-LUAD",
  cancer_name = "TCGA_LUAD",
  raw_dir_cancer = "TCGA_LUAD",
  staging_dir = gdc_luad_cnv_stage,
  run_download = FALSE,
  sync_staging = TRUE
)

lusc_cnv_objects <- prepare_tcga_cnv(
  project = "TCGA-LUSC",
  cancer_name = "TCGA_LUSC",
  raw_dir_cancer = "TCGA_LUSC",
  staging_dir = gdc_lusc_cnv_stage,
  run_download = FALSE,
  sync_staging = TRUE
)

cat(
  "\n========================================\n",
  "SCRIPT 07 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "LUAD CNV rows: ",
  nrow(luad_cnv_objects$combined),
  "\n",
  "LUSC CNV rows: ",
  nrow(lusc_cnv_objects$combined),
  "\n",
  "Outputs: ",
  cnv_processed_dir,
  "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================