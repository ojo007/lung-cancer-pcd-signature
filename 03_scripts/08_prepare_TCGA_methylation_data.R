# ============================================================
# MSc Lung Cancer PCD Project
# Script: 08_prepare_TCGA_methylation_data.R
#
# Purpose:
#   Prepare TCGA-LUAD and TCGA-LUSC Illumina Human
#   Methylation 450 (HM450) beta-value data.
#
# Main processing:
#   - Restrict to Primary Tumor HM450 samples
#   - Verify raw files against GDC manifests
#   - Check CpG probe structure/order
#   - Calculate sample- and probe-level missingness
#   - Retain CpGs with <=10% missing values
#   - Construct beta-value matrices
#   - Align sample metadata to matrix columns
#
# LUSC-specific QC:
#   - One extreme methylation sample was excluded:
#       TCGA-LA-A7SW
#   - Raw beta missingness = 31.99%
#   - Missingness after initial probe filtering = 17.08%
#   - Final LUSC cohort = 369 samples
#
# Notes:
#   - Remaining missing beta values are NOT imputed here.
#   - Cross-omics patient harmonization is performed later.
# ============================================================


# ============================================================
# 0. LOAD SHARED PROJECT CONFIGURATION
# ============================================================

source("03_scripts/00_project_config.R")


# ============================================================
# 1. REQUIRED PACKAGES AND OUTPUT DIRECTORY
# ============================================================

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

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
})

methylation_processed_dir <- file.path(
  processed_dir,
  "methylation"
)

ensure_dir(
  methylation_processed_dir
)


# ============================================================
# PART A — TCGA-LUAD HM450
# ============================================================


# ------------------------------------------------------------
# A1. Query methylation data
# ------------------------------------------------------------

query_luad_meth <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  access = "open"
)


luad_meth_results <- getResults(
  query_luad_meth
)


# ------------------------------------------------------------
# A2. Keep Primary Tumor HM450 samples
# ------------------------------------------------------------

luad_meth_450_results <- luad_meth_results[
  luad_meth_results$sample_type == "Primary Tumor" &
    luad_meth_results$platform ==
    "Illumina Human Methylation 450",
]


cat("\nLUAD HM450 Primary Tumor files:\n")

print(
  dim(
    luad_meth_450_results
  )
)


cat("\nLUAD HM450 unique patients:\n")

print(
  length(
    unique(
      luad_meth_450_results$cases.submitter_id
    )
  )
)


cat("\nLUAD HM450 unique samples:\n")

print(
  length(
    unique(
      luad_meth_450_results$sample.submitter_id
    )
  )
)


# ------------------------------------------------------------
# A3. Raw directory
# ------------------------------------------------------------

luad_meth_raw_dir <- file.path(
  raw_dir,
  "TCGA_LUAD",
  "methylation"
)


ensure_dir(
  luad_meth_raw_dir
)


# ------------------------------------------------------------
# A4. Manifest
# ------------------------------------------------------------

luad_meth_manifest <- data.frame(
  id = luad_meth_450_results$file_id,
  filename = luad_meth_450_results$file_name,
  md5 = luad_meth_450_results$md5sum,
  size = luad_meth_450_results$file_size,
  state = "validated",
  stringsAsFactors = FALSE
)


luad_manifest_file <- file.path(
  luad_meth_raw_dir,
  "gdc_manifest_TCGA_LUAD_HM450_methylation.txt"
)


write.table(
  luad_meth_manifest,
  file = luad_manifest_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# A5. External download
# ------------------------------------------------------------
#
# Download performed with gdc-client:
#
# gdc-client.exe download -m
# "...\gdc_manifest_TCGA_LUAD_HM450_methylation.txt"
# -d "C:\Users\ibrah\GDC\GDC_LUAD_METH"
#
# Verified:
#   473 expected files
#   0 missing
#   0 duplicate filenames
#   0 partial downloads
#   Total size = 5.783575 GiB
#
# Verified files were copied into:
# 01_raw_data/TCGA_LUAD/methylation/
# ------------------------------------------------------------


# ------------------------------------------------------------
# A6. Locate expected raw files
# ------------------------------------------------------------

luad_meth_expected_names <-
  luad_meth_manifest$filename


luad_meth_raw_files <- list.files(
  luad_meth_raw_dir,
  recursive = TRUE,
  full.names = TRUE
)


luad_meth_raw_correct <-
  luad_meth_raw_files[
    basename(
      luad_meth_raw_files
    ) %in%
      luad_meth_expected_names
  ]


stopifnot(
  length(
    luad_meth_raw_correct
  ) == 473
)


# ------------------------------------------------------------
# A7. Read reference HM450 file
# ------------------------------------------------------------

luad_meth_test <- fread(
  luad_meth_raw_correct[1],
  data.table = FALSE
)


stopifnot(
  nrow(
    luad_meth_test
  ) == 486427
)


luad_reference_probes <-
  luad_meth_test$V1


# ------------------------------------------------------------
# A8. Sample metadata
# ------------------------------------------------------------

luad_meth_metadata <-
  luad_meth_450_results[, c(
    "file_name",
    "file_id",
    "cases.submitter_id",
    "sample.submitter_id",
    "cases",
    "sample_type",
    "platform"
  )]


colnames(
  luad_meth_metadata
)[
  colnames(
    luad_meth_metadata
  ) == "cases.submitter_id"
] <- "patient_id"


write.csv(
  luad_meth_metadata,
  file.path(
    methylation_processed_dir,
    "TCGA_LUAD_HM450_sample_metadata.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# A9. Probe and sample missingness scan
# ------------------------------------------------------------

n_luad_meth_samples <-
  length(
    luad_meth_raw_correct
  )


n_luad_meth_probes <-
  length(
    luad_reference_probes
  )


luad_probe_missing_count <-
  integer(
    n_luad_meth_probes
  )


luad_meth_sample_qc <- data.frame(
  file = basename(
    luad_meth_raw_correct
  ),
  probes = NA_integer_,
  missing_beta = NA_integer_,
  missing_percent = NA_real_,
  stringsAsFactors = FALSE
)


for (
  i in seq_along(
    luad_meth_raw_correct
  )
) {
  
  x <- fread(
    luad_meth_raw_correct[i],
    select = 2,
    data.table = FALSE
  )
  
  
  beta <- x[[1]]
  
  
  if (
    length(beta) !=
    n_luad_meth_probes
  ) {
    
    stop(
      "Unexpected probe count in LUAD file: ",
      basename(
        luad_meth_raw_correct[i]
      )
    )
  }
  
  
  missing_here <-
    is.na(
      beta
    )
  
  
  luad_probe_missing_count <-
    luad_probe_missing_count +
    missing_here
  
  
  luad_meth_sample_qc$probes[i] <-
    length(
      beta
    )
  
  
  luad_meth_sample_qc$missing_beta[i] <-
    sum(
      missing_here
    )
  
  
  luad_meth_sample_qc$missing_percent[i] <-
    100 *
    mean(
      missing_here
    )
  
  
  if (
    i %% 25 == 0
  ) {
    
    message(
      "Scanned LUAD methylation files: ",
      i,
      " of ",
      length(
        luad_meth_raw_correct
      )
    )
  }
  
  
  rm(
    x,
    beta,
    missing_here
  )
}


# ------------------------------------------------------------
# A10. Probe QC
# ------------------------------------------------------------

luad_meth_probe_qc <- data.frame(
  CpG = luad_reference_probes,
  
  missing_count =
    luad_probe_missing_count,
  
  missing_percent =
    100 *
    luad_probe_missing_count /
    n_luad_meth_samples,
  
  stringsAsFactors = FALSE
)


write.csv(
  luad_meth_probe_qc,
  file.path(
    methylation_processed_dir,
    "TCGA_LUAD_HM450_probe_missingness_QC.csv"
  ),
  row.names = FALSE
)


write.csv(
  luad_meth_sample_qc,
  file.path(
    methylation_processed_dir,
    "TCGA_LUAD_HM450_sample_QC.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# A11. Retain probes with <=10% missingness
# ------------------------------------------------------------

luad_meth_keep_probes <-
  luad_meth_probe_qc$CpG[
    luad_meth_probe_qc$missing_percent <= 10
  ]


stopifnot(
  length(
    luad_meth_keep_probes
  ) == 406026
)


luad_meth_probe_filter <-
  luad_meth_probe_qc[
    luad_meth_probe_qc$missing_percent <= 10,
  ]


write.csv(
  luad_meth_probe_filter,
  file.path(
    methylation_processed_dir,
    "TCGA_LUAD_HM450_retained_probes_le10pct_missing.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# A12. Retained probe index
# ------------------------------------------------------------

luad_keep_index <- match(
  luad_meth_keep_probes,
  luad_reference_probes
)


stopifnot(
  sum(
    is.na(
      luad_keep_index
    )
  ) == 0
)


saveRDS(
  luad_keep_index,
  file.path(
    methylation_processed_dir,
    "TCGA_LUAD_HM450_retained_probe_index.rds"
  )
)


# ------------------------------------------------------------
# A13. Build filtered beta matrix
# ------------------------------------------------------------

luad_meth_matrix <- matrix(
  NA_real_,
  nrow =
    length(
      luad_meth_keep_probes
    ),
  ncol =
    length(
      luad_meth_raw_correct
    )
)


rownames(
  luad_meth_matrix
) <- luad_meth_keep_probes


colnames(
  luad_meth_matrix
) <- basename(
  luad_meth_raw_correct
)


for (
  i in seq_along(
    luad_meth_raw_correct
  )
) {
  
  x <- fread(
    luad_meth_raw_correct[i],
    select = 2,
    data.table = FALSE
  )
  
  
  beta <- x[[1]]
  
  
  luad_meth_matrix[, i] <-
    beta[
      luad_keep_index
    ]
  
  
  rm(
    x,
    beta
  )
  
  
  if (
    i %% 25 == 0
  ) {
    
    gc()
    
    
    message(
      "Loaded LUAD methylation samples: ",
      i,
      " of ",
      length(
        luad_meth_raw_correct
      )
    )
  }
}


# ------------------------------------------------------------
# A14. Matrix metadata alignment
# ------------------------------------------------------------

luad_meth_matrix_metadata <-
  luad_meth_metadata[
    match(
      colnames(
        luad_meth_matrix
      ),
      luad_meth_metadata$file_name
    ),
  ]


stopifnot(
  sum(
    is.na(
      luad_meth_matrix_metadata$patient_id
    )
  ) == 0
)


stopifnot(
  identical(
    colnames(
      luad_meth_matrix
    ),
    luad_meth_matrix_metadata$file_name
  )
)


# ------------------------------------------------------------
# A15. Save LUAD final outputs
# ------------------------------------------------------------

saveRDS(
  luad_meth_matrix,
  file.path(
    methylation_processed_dir,
    "TCGA_LUAD_HM450_beta_matrix.rds"
  ),
  compress = FALSE
)


write.csv(
  luad_meth_matrix_metadata,
  file.path(
    methylation_processed_dir,
    "TCGA_LUAD_HM450_matrix_column_metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART B — TCGA-LUSC HM450
# ============================================================


# ------------------------------------------------------------
# B1. Query methylation data
# ------------------------------------------------------------

query_lusc_meth <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  access = "open"
)


lusc_meth_results <- getResults(
  query_lusc_meth
)


# ------------------------------------------------------------
# B2. Keep Primary Tumor HM450
# ------------------------------------------------------------

lusc_meth_450_results <-
  lusc_meth_results[
    lusc_meth_results$sample_type ==
      "Primary Tumor" &
      lusc_meth_results$platform ==
      "Illumina Human Methylation 450",
  ]


stopifnot(
  nrow(
    lusc_meth_450_results
  ) == 370
)


# ------------------------------------------------------------
# B3. Raw directory
# ------------------------------------------------------------

lusc_meth_raw_dir <- file.path(
  raw_dir,
  "TCGA_LUSC",
  "methylation"
)


ensure_dir(
  lusc_meth_raw_dir
)


# ------------------------------------------------------------
# B4. Manifest
# ------------------------------------------------------------

lusc_meth_manifest <- data.frame(
  id = lusc_meth_450_results$file_id,
  filename = lusc_meth_450_results$file_name,
  md5 = lusc_meth_450_results$md5sum,
  size = lusc_meth_450_results$file_size,
  state = "validated",
  stringsAsFactors = FALSE
)


lusc_manifest_file <- file.path(
  lusc_meth_raw_dir,
  "gdc_manifest_TCGA_LUSC_HM450_methylation.txt"
)


write.table(
  lusc_meth_manifest,
  file = lusc_manifest_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# B5. External download
# ------------------------------------------------------------
#
# gdc-client download:
#
# 370 expected files
# 0 missing
# 0 duplicate expected filenames
# 0 partial files
# Total = 4.512992 GiB
#
# One gdc-client checksum warning occurred during download,
# but manual MD5 verification confirmed the downloaded file
# matched the GDC manifest exactly.
# ------------------------------------------------------------


# ------------------------------------------------------------
# B6. Locate raw files
# ------------------------------------------------------------

lusc_meth_expected_names <-
  lusc_meth_manifest$filename


lusc_meth_raw_files <- list.files(
  lusc_meth_raw_dir,
  recursive = TRUE,
  full.names = TRUE
)


lusc_meth_raw_correct <-
  lusc_meth_raw_files[
    basename(
      lusc_meth_raw_files
    ) %in%
      lusc_meth_expected_names
  ]


stopifnot(
  length(
    lusc_meth_raw_correct
  ) == 370
)


# ------------------------------------------------------------
# B7. Reference CpG structure
# ------------------------------------------------------------

lusc_meth_test <- fread(
  lusc_meth_raw_correct[1],
  data.table = FALSE
)


stopifnot(
  nrow(
    lusc_meth_test
  ) == 486427
)


lusc_reference_probes <-
  lusc_meth_test$V1


# ------------------------------------------------------------
# B8. Metadata
# ------------------------------------------------------------

lusc_meth_metadata <-
  lusc_meth_450_results[, c(
    "file_name",
    "file_id",
    "cases.submitter_id",
    "sample.submitter_id",
    "cases",
    "sample_type",
    "platform"
  )]


colnames(
  lusc_meth_metadata
)[
  colnames(
    lusc_meth_metadata
  ) == "cases.submitter_id"
] <- "patient_id"


write.csv(
  lusc_meth_metadata,
  file.path(
    methylation_processed_dir,
    "TCGA_LUSC_HM450_sample_metadata.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# B9. Missingness scan across 370 samples
# ------------------------------------------------------------

n_lusc_meth_samples <-
  length(
    lusc_meth_raw_correct
  )


n_lusc_meth_probes <-
  length(
    lusc_reference_probes
  )


lusc_probe_missing_count <-
  integer(
    n_lusc_meth_probes
  )


lusc_meth_sample_qc <- data.frame(
  file = basename(
    lusc_meth_raw_correct
  ),
  probes = NA_integer_,
  missing_beta = NA_integer_,
  missing_percent = NA_real_,
  stringsAsFactors = FALSE
)


for (
  i in seq_along(
    lusc_meth_raw_correct
  )
) {
  
  x <- fread(
    lusc_meth_raw_correct[i],
    select = 2,
    data.table = FALSE
  )
  
  
  beta <- x[[1]]
  
  
  if (
    length(beta) !=
    n_lusc_meth_probes
  ) {
    
    stop(
      "Unexpected probe count in LUSC file: ",
      basename(
        lusc_meth_raw_correct[i]
      )
    )
  }
  
  
  missing_here <-
    is.na(
      beta
    )
  
  
  lusc_probe_missing_count <-
    lusc_probe_missing_count +
    missing_here
  
  
  lusc_meth_sample_qc$probes[i] <-
    length(
      beta
    )
  
  
  lusc_meth_sample_qc$missing_beta[i] <-
    sum(
      missing_here
    )
  
  
  lusc_meth_sample_qc$missing_percent[i] <-
    100 *
    mean(
      missing_here
    )
  
  
  if (
    i %% 25 == 0
  ) {
    
    message(
      "Scanned LUSC methylation files: ",
      i,
      " of ",
      length(
        lusc_meth_raw_correct
      )
    )
  }
  
  
  rm(
    x,
    beta,
    missing_here
  )
}


# ------------------------------------------------------------
# B10. Identify sample-level QC outlier
# ------------------------------------------------------------

lusc_outlier_file <-
  paste0(
    "8e8ef4bc-c87f-45c3-9915-",
    "fa2364dc254f.methylation_array.",
    "sesame.level3betas.txt"
  )


lusc_outlier_index <- match(
  lusc_outlier_file,
  basename(
    lusc_meth_raw_correct
  )
)


stopifnot(
  !is.na(
    lusc_outlier_index
  )
)


lusc_outlier_beta <- fread(
  lusc_meth_raw_correct[
    lusc_outlier_index
  ],
  select = 2,
  data.table = FALSE
)[[1]]


lusc_outlier_missing <-
  is.na(
    lusc_outlier_beta
  )


# Raw missingness should be approximately 31.99%
cat(
  "\nLUSC excluded sample raw missingness (%):\n"
)

print(
  mean(
    lusc_outlier_missing
  ) * 100
)


# ------------------------------------------------------------
# B11. Recalculate probe QC without outlier
# ------------------------------------------------------------

lusc_probe_missing_count_369 <-
  lusc_probe_missing_count -
  as.integer(
    lusc_outlier_missing
  )


lusc_meth_probe_qc_369 <- data.frame(
  CpG = lusc_reference_probes,
  
  missing_count =
    lusc_probe_missing_count_369,
  
  missing_percent =
    100 *
    lusc_probe_missing_count_369 /
    369,
  
  stringsAsFactors = FALSE
)


lusc_meth_keep_probes_369 <-
  lusc_meth_probe_qc_369$CpG[
    lusc_meth_probe_qc_369$missing_percent <= 10
  ]


stopifnot(
  length(
    lusc_meth_keep_probes_369
  ) == 398357
)


# ------------------------------------------------------------
# B12. Final retained probe files
# ------------------------------------------------------------

lusc_meth_probe_filter_final <-
  lusc_meth_probe_qc_369[
    lusc_meth_probe_qc_369$missing_percent <= 10,
  ]


write.csv(
  lusc_meth_probe_qc_369,
  file.path(
    methylation_processed_dir,
    "TCGA_LUSC_HM450_probe_missingness_QC.csv"
  ),
  row.names = FALSE
)


write.csv(
  lusc_meth_probe_filter_final,
  file.path(
    methylation_processed_dir,
    "TCGA_LUSC_HM450_retained_probes_le10pct_missing.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# B13. Final retained probe index
# ------------------------------------------------------------

lusc_final_retained_probe_index <- match(
  lusc_meth_keep_probes_369,
  lusc_reference_probes
)


stopifnot(
  sum(
    is.na(
      lusc_final_retained_probe_index
    )
  ) == 0
)


saveRDS(
  lusc_final_retained_probe_index,
  file.path(
    methylation_processed_dir,
    "TCGA_LUSC_HM450_retained_probe_index.rds"
  )
)


# ------------------------------------------------------------
# B14. Remove excluded raw sample
# ------------------------------------------------------------

lusc_meth_good_files <-
  lusc_meth_raw_correct[
    basename(
      lusc_meth_raw_correct
    ) !=
      lusc_outlier_file
  ]


stopifnot(
  length(
    lusc_meth_good_files
  ) == 369
)


# ------------------------------------------------------------
# B15. Build final 369-sample beta matrix directly
# ------------------------------------------------------------

lusc_meth_matrix <- matrix(
  NA_real_,
  nrow =
    length(
      lusc_meth_keep_probes_369
    ),
  ncol =
    length(
      lusc_meth_good_files
    )
)


rownames(
  lusc_meth_matrix
) <- lusc_meth_keep_probes_369


colnames(
  lusc_meth_matrix
) <- basename(
  lusc_meth_good_files
)


for (
  i in seq_along(
    lusc_meth_good_files
  )
) {
  
  x <- fread(
    lusc_meth_good_files[i],
    select = 2,
    data.table = FALSE
  )
  
  
  beta <- x[[1]]
  
  
  lusc_meth_matrix[, i] <-
    beta[
      lusc_final_retained_probe_index
    ]
  
  
  rm(
    x,
    beta
  )
  
  
  if (
    i %% 25 == 0
  ) {
    
    gc()
    
    
    message(
      "Loaded final LUSC methylation samples: ",
      i,
      " of ",
      length(
        lusc_meth_good_files
      )
    )
  }
}


stopifnot(
  identical(
    dim(
      lusc_meth_matrix
    ),
    c(
      398357L,
      369L
    )
  )
)


# ------------------------------------------------------------
# B16. Final sample QC
# ------------------------------------------------------------

lusc_final_sample_missing <- data.frame(
  file =
    colnames(
      lusc_meth_matrix
    ),
  
  missing_beta =
    colSums(
      is.na(
        lusc_meth_matrix
      )
    ),
  
  missing_percent =
    100 *
    colMeans(
      is.na(
        lusc_meth_matrix
      )
    ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# B17. Final matrix metadata
# ------------------------------------------------------------

lusc_meth_matrix_metadata <-
  lusc_meth_metadata[
    match(
      colnames(
        lusc_meth_matrix
      ),
      lusc_meth_metadata$file_name
    ),
  ]


stopifnot(
  sum(
    is.na(
      lusc_meth_matrix_metadata$patient_id
    )
  ) == 0
)


stopifnot(
  identical(
    colnames(
      lusc_meth_matrix
    ),
    lusc_meth_matrix_metadata$file_name
  )
)


stopifnot(
  length(
    unique(
      lusc_meth_matrix_metadata$patient_id
    )
  ) == 369
)


# ------------------------------------------------------------
# B18. Record excluded sample
# ------------------------------------------------------------

lusc_meth_excluded_samples <-
  lusc_meth_metadata[
    lusc_meth_metadata$file_name ==
      lusc_outlier_file,
  ]


lusc_meth_excluded_samples$raw_missing_percent <-
  mean(
    lusc_outlier_missing
  ) * 100


lusc_meth_excluded_samples$reason <-
  paste(
    "Excluded as methylation sample-level QC outlier;",
    "31.99% missing beta values before probe filtering",
    "and 17.08% missing after initial <=10%",
    "probe-missingness filtering."
  )


# ------------------------------------------------------------
# B19. Save final LUSC outputs
# ------------------------------------------------------------

saveRDS(
  lusc_meth_matrix,
  file.path(
    methylation_processed_dir,
    "TCGA_LUSC_HM450_beta_matrix.rds"
  ),
  compress = FALSE
)


write.csv(
  lusc_meth_matrix_metadata,
  file.path(
    methylation_processed_dir,
    "TCGA_LUSC_HM450_matrix_column_metadata.csv"
  ),
  row.names = FALSE
)


write.csv(
  lusc_final_sample_missing,
  file.path(
    methylation_processed_dir,
    "TCGA_LUSC_HM450_sample_QC.csv"
  ),
  row.names = FALSE
)


write.csv(
  lusc_meth_excluded_samples,
  file.path(
    methylation_processed_dir,
    "TCGA_LUSC_HM450_excluded_samples_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# FINAL QC SUMMARY
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "Methylation preparation complete.\n"
)

cat(
  "========================================\n"
)


cat(
  "\nLUAD HM450 matrix dimensions:\n"
)

print(
  dim(
    luad_meth_matrix
  )
)


cat(
  "\nLUAD overall remaining missingness (%):\n"
)

print(
  mean(
    is.na(
      luad_meth_matrix
    )
  ) * 100
)


cat(
  "\nLUAD unique patients represented:\n"
)

print(
  length(
    unique(
      luad_meth_matrix_metadata$patient_id
    )
  )
)


cat(
  "\nLUSC HM450 final matrix dimensions:\n"
)

print(
  dim(
    lusc_meth_matrix
  )
)


cat(
  "\nLUSC overall remaining missingness (%):\n"
)

print(
  mean(
    is.na(
      lusc_meth_matrix
    )
  ) * 100
)


cat(
  "\nLUSC unique patients represented:\n"
)

print(
  length(
    unique(
      lusc_meth_matrix_metadata$patient_id
    )
  )
)


cat(
  "\nLUSC maximum sample missingness after QC (%):\n"
)

print(
  max(
    lusc_final_sample_missing$missing_percent
  )
)


cat(
  "\nProcessed methylation files:\n"
)

print(
  list.files(
    methylation_processed_dir,
    pattern = "TCGA_(LUAD|LUSC)_HM450"
  )
)


cat(
  "\n========================================\n",
  "SCRIPT 08 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  sep = ""
)


# ============================================================
# End of script
# ============================================================