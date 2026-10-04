# ============================================================
# MSc Lung Cancer PCD Project
# Script: 05_prepare_TCGA_clinical_survival_data.R
#
# Purpose:
#   Download, clean, validate, and prepare TCGA-LUAD and
#   TCGA-LUSC clinical data for downstream survival analyses.
#
# Outputs:
#   02_processed_data/clinical/
#     - TCGA_LUAD_clinical_raw.rds
#     - TCGA_LUAD_clinical_clean.rds
#     - TCGA_LUAD_clinical_clean.csv
#     - TCGA_LUAD_survival_final.rds
#     - TCGA_LUAD_survival_final.csv
#     - TCGA_LUSC_clinical_raw.rds
#     - TCGA_LUSC_clinical_clean.rds
#     - TCGA_LUSC_clinical_clean.csv
#     - TCGA_LUSC_survival_final.rds
#     - TCGA_LUSC_survival_final.csv
# ============================================================

# ------------------------------------------------------------
# 0. Load shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")

# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  stop("Package 'TCGAbiolinks' is required.")
}

# ------------------------------------------------------------
# 2. Output directory
# ------------------------------------------------------------

clinical_dir <- file.path(
  processed_dir,
  "clinical"
)

ensure_dir(clinical_dir)

# ------------------------------------------------------------
# 3. Shared clinical-preparation function
# ------------------------------------------------------------

prepare_tcga_clinical <- function(
    project,
    cancer_name
) {
  
  cat(
    "\n========================================\n",
    "Preparing clinical data for ",
    project,
    "\n",
    "========================================\n",
    sep = ""
  )
  
  # ----------------------------------------------------------
  # 3.1 Download clinical data
  # ----------------------------------------------------------
  
  clinical_raw <- TCGAbiolinks::GDCquery_clinic(
    project = project,
    type = "clinical"
  )
  
  if (!is.data.frame(clinical_raw)) {
    clinical_raw <- as.data.frame(
      clinical_raw,
      stringsAsFactors = FALSE
    )
  }
  
  cat("\nRaw clinical dimensions:\n")
  print(dim(clinical_raw))
  
  if (nrow(clinical_raw) == 0) {
    stop("No clinical records returned for ", project)
  }
  
  # Preserve the full object for provenance.
  saveRDS(
    clinical_raw,
    file.path(
      clinical_dir,
      paste0(
        cancer_name,
        "_clinical_raw.rds"
      )
    )
  )
  
  # ----------------------------------------------------------
  # 3.2 Validate key identifiers / survival variables
  # ----------------------------------------------------------
  
  required_core <- c(
    "submitter_id",
    "vital_status",
    "days_to_death",
    "days_to_last_follow_up"
  )
  
  missing_core <- setdiff(
    required_core,
    colnames(clinical_raw)
  )
  
  if (length(missing_core) > 0) {
    stop(
      project,
      " clinical data are missing required field(s): ",
      paste(missing_core, collapse = ", ")
    )
  }
  
  if (anyNA(clinical_raw$submitter_id)) {
    stop(project, " contains missing submitter_id values.")
  }
  
  # ----------------------------------------------------------
  # 3.3 Select variables used downstream
  # ----------------------------------------------------------
  
  variables_to_keep <- c(
    "submitter_id",
    "bcr_patient_barcode",
    "vital_status",
    "days_to_death",
    "days_to_last_follow_up",
    "age_at_diagnosis",
    "age_at_index",
    "sex_at_birth",
    "gender",
    "race",
    "ethnicity",
    "ajcc_pathologic_stage",
    "ajcc_pathologic_t",
    "ajcc_pathologic_n",
    "ajcc_pathologic_m",
    "primary_diagnosis",
    "year_of_diagnosis",
    "tissue_or_organ_of_origin",
    "laterality",
    "residual_disease",
    "prior_malignancy",
    "prior_treatment",
    "pack_years_smoked",
    "cigarettes_per_day",
    "tobacco_smoking_status"
  )
  
  variables_to_keep <- intersect(
    variables_to_keep,
    colnames(clinical_raw)
  )
  
  clinical_clean <- clinical_raw[
    ,
    variables_to_keep,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # 3.4 Ensure one row per patient
  #
  # GDCquery_clinic normally returns case-level rows. If an
  # unexpected duplicate patient is returned, retain the first
  # row only after reporting it. No aggregation is performed.
  # ----------------------------------------------------------
  
  duplicate_patients <- clinical_clean$submitter_id[
    duplicated(clinical_clean$submitter_id)
  ]
  
  cat(
    "\nDuplicated patient IDs in raw clinical data:",
    length(unique(duplicate_patients)),
    "\n"
  )
  
  if (length(duplicate_patients) > 0) {
    clinical_clean <- clinical_clean[
      !duplicated(clinical_clean$submitter_id),
      ,
      drop = FALSE
    ]
  }
  
  # ----------------------------------------------------------
  # 3.5 Standardize sex field for downstream scripts
  #
  # Keep the original field(s), but ensure a consistent `sex`
  # variable exists when possible.
  # ----------------------------------------------------------
  
  if ("sex_at_birth" %in% colnames(clinical_clean)) {
    
    clinical_clean$sex <- clinical_clean$sex_at_birth
    
  } else if ("gender" %in% colnames(clinical_clean)) {
    
    clinical_clean$sex <- clinical_clean$gender
    
  } else {
    
    clinical_clean$sex <- NA_character_
  }
  
  # ----------------------------------------------------------
  # 3.6 Overall-survival event
  #
  # 1 = Dead
  # 0 = Alive / censored
  # NA = unknown
  # ----------------------------------------------------------
  
  vital_status_normalized <- trimws(
    as.character(
      clinical_clean$vital_status
    )
  )
  
  clinical_clean$OS_status <- ifelse(
    vital_status_normalized == "Dead",
    1L,
    ifelse(
      vital_status_normalized == "Alive",
      0L,
      NA_integer_
    )
  )
  
  # ----------------------------------------------------------
  # 3.7 Overall-survival time in days
  #
  # Dead  -> days_to_death
  # Alive -> days_to_last_follow_up
  # ----------------------------------------------------------
  
  clinical_clean$days_to_death <- suppressWarnings(
    as.numeric(
      clinical_clean$days_to_death
    )
  )
  
  clinical_clean$days_to_last_follow_up <- suppressWarnings(
    as.numeric(
      clinical_clean$days_to_last_follow_up
    )
  )
  
  clinical_clean$OS_time <- ifelse(
    clinical_clean$OS_status == 1L,
    clinical_clean$days_to_death,
    ifelse(
      clinical_clean$OS_status == 0L,
      clinical_clean$days_to_last_follow_up,
      NA_real_
    )
  )
  
  # ----------------------------------------------------------
  # 3.8 Age harmonization
  #
  # GDC age_at_diagnosis is stored in days. age_at_index,
  # where available, is already in years.
  # ----------------------------------------------------------
  
  if ("age_at_diagnosis" %in% colnames(clinical_clean)) {
    
    clinical_clean$age_at_diagnosis <- suppressWarnings(
      as.numeric(
        clinical_clean$age_at_diagnosis
      )
    )
    
    clinical_clean$age_at_diagnosis_years <-
      clinical_clean$age_at_diagnosis / 365.25
    
  } else {
    
    clinical_clean$age_at_diagnosis_years <- NA_real_
  }
  
  if ("age_at_index" %in% colnames(clinical_clean)) {
    
    clinical_clean$age_at_index <- suppressWarnings(
      as.numeric(
        clinical_clean$age_at_index
      )
    )
  }
  
  # ----------------------------------------------------------
  # 3.9 Save cleaned clinical data
  # ----------------------------------------------------------
  
  saveRDS(
    clinical_clean,
    file.path(
      clinical_dir,
      paste0(
        cancer_name,
        "_clinical_clean.rds"
      )
    )
  )
  
  write.csv(
    clinical_clean,
    file.path(
      clinical_dir,
      paste0(
        cancer_name,
        "_clinical_clean.csv"
      )
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 3.10 Create survival-ready dataset
  # ----------------------------------------------------------
  
  survival_final <- clinical_clean[
    !is.na(clinical_clean$OS_status) &
      !is.na(clinical_clean$OS_time) &
      clinical_clean$OS_time > 0,
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # 3.11 Survival QC
  # ----------------------------------------------------------
  
  cat("\nClean clinical dimensions:\n")
  print(dim(clinical_clean))
  
  cat("\nVital-status distribution:\n")
  print(
    table(
      clinical_clean$vital_status,
      useNA = "ifany"
    )
  )
  
  cat(
    "\nMissing OS status:",
    sum(is.na(clinical_clean$OS_status)),
    "\n"
  )
  
  cat(
    "Missing OS time:",
    sum(is.na(clinical_clean$OS_time)),
    "\n"
  )
  
  cat(
    "OS time <= 0 before survival filtering:",
    sum(
      !is.na(clinical_clean$OS_time) &
        clinical_clean$OS_time <= 0
    ),
    "\n"
  )
  
  cat("\nFinal survival dimensions:\n")
  print(dim(survival_final))
  
  cat("\nFinal OS-status distribution:\n")
  print(
    table(
      survival_final$OS_status
    )
  )
  
  cat(
    "\nFinal unique patients:",
    length(
      unique(
        survival_final$submitter_id
      )
    ),
    "\n"
  )
  
  stopifnot(
    nrow(survival_final) > 0,
    !anyNA(survival_final$submitter_id),
    !anyDuplicated(survival_final$submitter_id),
    all(survival_final$OS_time > 0),
    all(survival_final$OS_status %in% c(0L, 1L))
  )
  
  # ----------------------------------------------------------
  # 3.12 Save survival-ready data
  # ----------------------------------------------------------
  
  saveRDS(
    survival_final,
    file.path(
      clinical_dir,
      paste0(
        cancer_name,
        "_survival_final.rds"
      )
    )
  )
  
  write.csv(
    survival_final,
    file.path(
      clinical_dir,
      paste0(
        cancer_name,
        "_survival_final.csv"
      )
    ),
    row.names = FALSE
  )
  
  cat(
    "\n",
    cancer_name,
    " clinical/survival preparation complete.\n",
    sep = ""
  )
  
  invisible(
    list(
      raw = clinical_raw,
      clean = clinical_clean,
      survival = survival_final
    )
  )
}

# ------------------------------------------------------------
# 4. Prepare TCGA-LUAD
# ------------------------------------------------------------

luad_clinical_objects <- prepare_tcga_clinical(
  project = "TCGA-LUAD",
  cancer_name = "TCGA_LUAD"
)

# ------------------------------------------------------------
# 5. Prepare TCGA-LUSC
# ------------------------------------------------------------

lusc_clinical_objects <- prepare_tcga_clinical(
  project = "TCGA-LUSC",
  cancer_name = "TCGA_LUSC"
)

# ------------------------------------------------------------
# 6. Final summary
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "SCRIPT 05 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "LUAD survival patients: ",
  nrow(luad_clinical_objects$survival),
  "\n",
  "LUSC survival patients: ",
  nrow(lusc_clinical_objects$survival),
  "\n",
  "Outputs: ",
  clinical_dir,
  "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================