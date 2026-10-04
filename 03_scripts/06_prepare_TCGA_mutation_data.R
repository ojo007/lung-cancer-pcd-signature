# ============================================================
# MSc Lung Cancer PCD Project
# Script: 06_prepare_TCGA_mutation_data.R
#
# Purpose:
#   Query, organize, combine, and prepare TCGA-LUAD and
#   TCGA-LUSC Masked Somatic Mutation (MAF) data.
#
# Important:
#   - Raw MAF files are retained in canonical raw storage.
#   - Empty MAF files are skipped when building combined tables.
#   - Multiple tumor aliquots per patient are retained here.
#   - Sample resolution is performed later during multi-omics
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

mutation_processed_dir <- file.path(
  processed_dir,
  "mutation"
)

ensure_dir(mutation_processed_dir)

# ------------------------------------------------------------
# Shared helper: synchronize staged GDC files to canonical raw
# storage using the GDC file ID as the subdirectory name.
# ------------------------------------------------------------

sync_gdc_stage <- function(
    manifest,
    staging_dir,
    canonical_dir,
    pattern = "\\.maf\\.gz$"
) {
  
  if (!dir.exists(staging_dir)) {
    return(invisible(NULL))
  }
  
  staged_files <- list.files(
    staging_dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  staged_partial <- list.files(
    staging_dir,
    pattern = "\\.partial$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  cat(
    "\nStaged MAF files:",
    length(staged_files),
    "\n"
  )
  
  cat(
    "Staged partial files:",
    length(staged_partial),
    "\n"
  )
  
  if (length(staged_partial) > 0) {
    warning(
      "Partial mutation downloads were detected in staging."
    )
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
          "Failed to copy staged mutation file:\n",
          available_stage$path[i]
        )
      }
    }
  }
  
  invisible(NULL)
}

# ------------------------------------------------------------
# Shared helper: prepare mutation data for one TCGA project
# ------------------------------------------------------------

prepare_tcga_mutation <- function(
    project,
    cancer_name,
    raw_dir_cancer,
    staging_dir,
    run_download = FALSE,
    sync_staging = TRUE
) {
  
  cat(
    "\n========================================\n",
    "Preparing mutation data for ",
    project,
    "\n",
    "========================================\n",
    sep = ""
  )
  
  # ----------------------------------------------------------
  # 1. Query GDC
  # ----------------------------------------------------------
  
  query_mut <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    workflow.type =
      "Aliquot Ensemble Somatic Variant Merging and Masking",
    access = "open"
  )
  
  mut_results <- TCGAbiolinks::getResults(
    query_mut
  )
  
  cat("\nMutation query dimensions:\n")
  print(dim(mut_results))
  
  stopifnot(
    nrow(mut_results) > 0,
    all(
      c(
        "file_id",
        "file_name",
        "md5sum",
        "file_size"
      ) %in% colnames(mut_results)
    )
  )
  
  # ----------------------------------------------------------
  # 2. Canonical raw directory + manifest
  # ----------------------------------------------------------
  
  mutation_raw_dir <- file.path(
    raw_dir,
    raw_dir_cancer,
    "mutation"
  )
  
  ensure_dir(mutation_raw_dir)
  
  mutation_manifest <- data.frame(
    id = mut_results$file_id,
    filename = mut_results$file_name,
    md5 = mut_results$md5sum,
    size = mut_results$file_size,
    state = "validated",
    stringsAsFactors = FALSE
  )
  
  manifest_file <- file.path(
    mutation_raw_dir,
    paste0(
      "gdc_manifest_",
      cancer_name,
      "_mutation.txt"
    )
  )
  
  write.table(
    mutation_manifest,
    file = manifest_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  cat("\nManifest dimensions:\n")
  print(dim(mutation_manifest))
  
  cat(
    "\nExpected mutation download size (MiB): ",
    sum(mut_results$file_size, na.rm = TRUE) / 1024^2,
    "\n",
    sep = ""
  )
  
  # ----------------------------------------------------------
  # 3. Optional download to short staging path
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
    sync_gdc_stage(
      manifest = mutation_manifest,
      staging_dir = staging_dir,
      canonical_dir = mutation_raw_dir
    )
  }
  
  # ----------------------------------------------------------
  # 5. Locate canonical MAF files
  # ----------------------------------------------------------
  
  maf_files_all <- list.files(
    mutation_raw_dir,
    pattern = "\\.maf\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(maf_files_all) == 0) {
    stop(
      "No canonical MAF files found for ",
      project,
      " under:\n",
      mutation_raw_dir
    )
  }
  
  # Restrict to files belonging to the current query manifest.
  maf_files <- maf_files_all[
    basename(maf_files_all) %in%
      mutation_manifest$filename
  ]
  
  unexpected_maf <- maf_files_all[
    !basename(maf_files_all) %in%
      mutation_manifest$filename
  ]
  
  missing_manifest_files <- setdiff(
    mutation_manifest$filename,
    basename(maf_files)
  )
  
  cat(
    "\nCanonical manifest-matching MAF files: ",
    length(maf_files),
    "\n",
    sep = ""
  )
  
  cat(
    "Unexpected canonical MAF files: ",
    length(unexpected_maf),
    "\n",
    sep = ""
  )
  
  cat(
    "Manifest files missing from canonical storage: ",
    length(missing_manifest_files),
    "\n",
    sep = ""
  )
  
  if (length(unexpected_maf) > 0) {
    warning(
      length(unexpected_maf),
      " unexpected MAF file(s) exist in ",
      mutation_raw_dir,
      ". They will not be used."
    )
  }
  
  if (length(missing_manifest_files) > 0) {
    stop(
      "Canonical mutation storage for ",
      project,
      " is incomplete. Missing ",
      length(missing_manifest_files),
      " manifest file(s)."
    )
  }
  
  # Order files according to manifest order for reproducibility.
  maf_match <- match(
    mutation_manifest$filename,
    basename(maf_files)
  )
  
  stopifnot(
    !anyNA(maf_match)
  )
  
  maf_files <- maf_files[
    maf_match
  ]
  
  # ----------------------------------------------------------
  # 6. Read and combine MAF files
  # ----------------------------------------------------------
  
  maf_list <- vector(
    "list",
    length(maf_files)
  )
  
  empty_maf_files <- character()
  
  for (i in seq_along(maf_files)) {
    
    x <- tryCatch(
      data.table::fread(
        maf_files[i],
        skip = "Hugo_Symbol",
        data.table = FALSE
      ),
      error = function(e) {
        data.frame()
      }
    )
    
    if (nrow(x) == 0) {
      
      empty_maf_files <- c(
        empty_maf_files,
        maf_files[i]
      )
      
      next
    }
    
    required_maf_columns <- c(
      "Hugo_Symbol",
      "Variant_Classification",
      "Tumor_Sample_Barcode"
    )
    
    if (!all(required_maf_columns %in% colnames(x))) {
      stop(
        "Unexpected MAF structure in file:\n",
        maf_files[i]
      )
    }
    
    x$source_file <- basename(
      maf_files[i]
    )
    
    maf_list[[i]] <- x
    
    if (i %% 50 == 0 || i == length(maf_files)) {
      message(
        "Processed ",
        cancer_name,
        " MAF files: ",
        i,
        " of ",
        length(maf_files)
      )
    }
  }
  
  maf_list <- maf_list[
    !vapply(
      maf_list,
      is.null,
      logical(1)
    )
  ]
  
  if (length(maf_list) == 0) {
    stop(
      "All MAF files were empty for ",
      project,
      "."
    )
  }
  
  mutation_combined <- data.table::rbindlist(
    maf_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  mutation_combined <- as.data.frame(
    mutation_combined
  )
  
  # ----------------------------------------------------------
  # 7. Patient/sample identifiers
  # ----------------------------------------------------------
  
  mutation_combined$patient_id <- substr(
    mutation_combined$Tumor_Sample_Barcode,
    1,
    12
  )
  
  patient_sample_map <- unique(
    mutation_combined[
      ,
      c(
        "patient_id",
        "Tumor_Sample_Barcode"
      ),
      drop = FALSE
    ]
  )
  
  mutation_sample_summary <- as.data.frame(
    table(
      mutation_combined$Tumor_Sample_Barcode
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(
    mutation_sample_summary
  ) <- c(
    "Tumor_Sample_Barcode",
    "Mutation_Count"
  )
  
  mutation_sample_summary$patient_id <- substr(
    mutation_sample_summary$Tumor_Sample_Barcode,
    1,
    12
  )
  
  mutation_sample_summary$sample_type_code <- substr(
    mutation_sample_summary$Tumor_Sample_Barcode,
    14,
    15
  )
  
  multi_patients <- names(
    which(
      table(
        mutation_sample_summary$patient_id
      ) > 1
    )
  )
  
  # ----------------------------------------------------------
  # 8. QC
  # ----------------------------------------------------------
  
  cat("\nCombined mutation dimensions:\n")
  print(dim(mutation_combined))
  
  cat(
    "\nEmpty MAF files: ",
    length(empty_maf_files),
    "\n",
    sep = ""
  )
  
  cat(
    "Mutation-bearing samples: ",
    length(
      unique(
        mutation_combined$Tumor_Sample_Barcode
      )
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Unique patients: ",
    length(
      unique(
        mutation_combined$patient_id
      )
    ),
    "\n",
    sep = ""
  )
  
  cat("\nSample type codes:\n")
  print(
    table(
      mutation_sample_summary$sample_type_code
    )
  )
  
  cat(
    "\nPatients with multiple mutation samples: ",
    length(multi_patients),
    "\n",
    sep = ""
  )
  
  cat("\nMutation classifications:\n")
  print(
    sort(
      table(
        mutation_combined$Variant_Classification
      ),
      decreasing = TRUE
    )
  )
  
  stopifnot(
    nrow(mutation_combined) > 0,
    !anyNA(mutation_combined$Tumor_Sample_Barcode),
    !anyNA(mutation_combined$patient_id),
    all(nchar(mutation_combined$patient_id) == 12)
  )
  
  # ----------------------------------------------------------
  # 9. Save outputs
  # ----------------------------------------------------------
  
  saveRDS(
    mutation_combined,
    file.path(
      mutation_processed_dir,
      paste0(
        cancer_name,
        "_mutation_raw_combined.rds"
      )
    )
  )
  
  write.csv(
    patient_sample_map,
    file.path(
      mutation_processed_dir,
      paste0(
        cancer_name,
        "_mutation_sample_map.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    mutation_sample_summary,
    file.path(
      mutation_processed_dir,
      paste0(
        cancer_name,
        "_mutation_sample_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  invisible(
    list(
      manifest = mutation_manifest,
      combined = mutation_combined,
      sample_map = patient_sample_map,
      sample_summary = mutation_sample_summary,
      empty_files = empty_maf_files
    )
  )
}

# ------------------------------------------------------------
# TCGA-LUAD
# ------------------------------------------------------------

luad_mutation_objects <- prepare_tcga_mutation(
  project = "TCGA-LUAD",
  cancer_name = "TCGA_LUAD",
  raw_dir_cancer = "TCGA_LUAD",
  staging_dir = gdc_luad_mut_stage,
  run_download = FALSE,
  sync_staging = TRUE
)

# ------------------------------------------------------------
# TCGA-LUSC
# ------------------------------------------------------------

lusc_mutation_objects <- prepare_tcga_mutation(
  project = "TCGA-LUSC",
  cancer_name = "TCGA_LUSC",
  raw_dir_cancer = "TCGA_LUSC",
  staging_dir = gdc_lusc_mut_stage,
  run_download = FALSE,
  sync_staging = TRUE
)

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "SCRIPT 06 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "LUAD mutation rows: ",
  nrow(luad_mutation_objects$combined),
  "\n",
  "LUSC mutation rows: ",
  nrow(lusc_mutation_objects$combined),
  "\n",
  "Outputs: ",
  mutation_processed_dir,
  "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================