# ============================================================
# SCRIPT 27
# TIDE / IPS immunotherapy-surrogate feasibility
#
# TCGA-LUAD
#
# FINAL CONSOLIDATED CLEAN-RUN VERSION
#
# IMPORTANT:
# - TIDE is an in silico immunotherapy-response surrogate.
# - Predicted responder status is NOT observed clinical response.
# - TIDE components are model-derived expression scores.
# - Nonsynonymous mutation burden is NOT formal TMB because
#   no callable-megabase denominator was available.
# - IPS was assessed for feasibility but was not available
#   locally for direct analysis.
# ============================================================


source("03_scripts/00_project_config.R")


cat(
  "\n========================================\n",
  "SCRIPT 27: TIDE / IPS IMMUNOTHERAPY SURROGATE\n",
  "========================================\n"
)


# ============================================================
# SECTION 1
# Define inputs and output directories
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 1: INPUT AND OUTPUT SETUP\n",
  "========================================\n"
)


immunotherapy_processed_dir_s31 <-
  "02_processed_data/immunotherapy"


immunotherapy_results_dir_s31 <-
  "04_results/immunotherapy"


script26_indicator_dir_s31 <-
  "04_results/mpcds/immunotherapy_indicators"



dir.create(
  immunotherapy_processed_dir_s31,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  immunotherapy_results_dir_s31,
  recursive = TRUE,
  showWarnings = FALSE
)



# ------------------------------------------------------------
# Exact validated Script 27 processed files
# ------------------------------------------------------------


patient_tpm_file_s31 <-
  file.path(
    immunotherapy_processed_dir_s31,
    "TCGA_LUAD_patient_protein_coding_TPM.rds"
  )


tide_log2_file_s31 <-
  file.path(
    immunotherapy_processed_dir_s31,
    "TCGA_LUAD_TIDE_log2TPM.rds"
  )


tide_centered_file_s31 <-
  file.path(
    immunotherapy_processed_dir_s31,
    "TCGA_LUAD_TIDE_centered_expression.rds"
  )


tide_input_file_s31 <-
  file.path(
    immunotherapy_processed_dir_s31,
    "TCGA_LUAD_TIDE_input.tsv"
  )


tide_results_file_s31 <-
  file.path(
    immunotherapy_processed_dir_s31,
    "TCGA_LUAD_TIDE_results.tsv"
  )


tide_clean_results_file_s31 <-
  file.path(
    immunotherapy_processed_dir_s31,
    "TCGA_LUAD_TIDE_patient_results.csv"
  )


subtype_file_s31 <-
  "04_results/clustering/TCGA_LUAD_PCD_cluster_assignments.csv"


mpcds_file_s31 <-
  "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv"


required_files_initial_s31 <-
  c(
    patient_tpm_file_s31,
    tide_results_file_s31,
    subtype_file_s31,
    mpcds_file_s31
  )


cat(
  "\nRequired starting files:\n"
)


print(
  data.frame(
    file = required_files_initial_s31,
    exists = file.exists(required_files_initial_s31),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


stopifnot(
  all(
    file.exists(
      required_files_initial_s31
    )
  )
)


# ============================================================
# SECTION 2
# Validate patient-level protein-coding expression
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 2: PATIENT-LEVEL EXPRESSION\n",
  "========================================\n"
)


luad_pc_tpm_by_symbol_s31 <-
  readRDS(
    patient_tpm_file_s31
  )


cat(
  "\nExpression class:\n"
)


print(
  class(
    luad_pc_tpm_by_symbol_s31
  )
)


cat(
  "\nExpression dimensions:\n"
)


print(
  dim(
    luad_pc_tpm_by_symbol_s31
  )
)


stopifnot(
  is.matrix(
    luad_pc_tpm_by_symbol_s31
  )
)


stopifnot(
  nrow(
    luad_pc_tpm_by_symbol_s31
  ) ==
    19938
)


stopifnot(
  ncol(
    luad_pc_tpm_by_symbol_s31
  ) ==
    517
)


stopifnot(
  !any(
    duplicated(
      rownames(
        luad_pc_tpm_by_symbol_s31
      )
    )
  )
)


stopifnot(
  !any(
    duplicated(
      colnames(
        luad_pc_tpm_by_symbol_s31
      )
    )
  )
)


stopifnot(
  !anyNA(
    luad_pc_tpm_by_symbol_s31
  )
)


stopifnot(
  all(
    is.finite(
      luad_pc_tpm_by_symbol_s31
    )
  )
)


stopifnot(
  all(
    luad_pc_tpm_by_symbol_s31 >=
      0
  )
)


# ============================================================
# SECTION 3
# Reconstruct validated TIDE normalization
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 3: TIDE NORMALIZATION\n",
  "========================================\n"
)


luad_tide_log2_s31 <-
  log2(
    luad_pc_tpm_by_symbol_s31 +
      1
  )


luad_tide_centered_all_s31 <-
  sweep(
    luad_tide_log2_s31,
    1,
    rowMeans(
      luad_tide_log2_s31
    ),
    "-"
  )


gene_sd_s31 <-
  apply(
    luad_tide_centered_all_s31,
    1,
    sd
  )


zero_sd_genes_s31 <-
  names(
    gene_sd_s31
  )[
    is.finite(
      gene_sd_s31
    ) &
      gene_sd_s31 ==
      0
  ]


cat(
  "\nZero-SD genes removed:\n"
)


print(
  length(
    zero_sd_genes_s31
  )
)


luad_tide_input_s31 <-
  luad_tide_centered_all_s31[
    is.finite(
      gene_sd_s31
    ) &
      gene_sd_s31 >
      0,
    ,
    drop = FALSE
  ]


cat(
  "\nFinal TIDE input dimensions:\n"
)


print(
  dim(
    luad_tide_input_s31
  )
)


stopifnot(
  nrow(
    luad_tide_input_s31
  ) ==
    19496
)


stopifnot(
  ncol(
    luad_tide_input_s31
  ) ==
    517
)


stopifnot(
  !anyNA(
    luad_tide_input_s31
  )
)


stopifnot(
  all(
    is.finite(
      luad_tide_input_s31
    )
  )
)


stopifnot(
  !any(
    duplicated(
      rownames(
        luad_tide_input_s31
      )
    )
  )
)


stopifnot(
  !any(
    duplicated(
      colnames(
        luad_tide_input_s31
      )
    )
  )
)


max_abs_centered_mean_s31 <-
  max(
    abs(
      rowMeans(
        luad_tide_input_s31
      )
    )
  )


cat(
  "\nMaximum absolute post-centering gene mean:\n"
)


print(
  max_abs_centered_mean_s31
)


saveRDS(
  luad_tide_log2_s31,
  tide_log2_file_s31
)


saveRDS(
  luad_tide_input_s31,
  tide_centered_file_s31
)


tide_input_dataframe_s31 <-
  data.frame(
    Gene =
      rownames(
        luad_tide_input_s31
      ),
    luad_tide_input_s31,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


write.table(
  tide_input_dataframe_s31,
  tide_input_file_s31,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "\nSaved normalized TIDE input:\n"
)


print(
  tide_input_file_s31
)


# ============================================================
# SECTION 4
# External TIDE execution note
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 4: TIDE EXECUTION STATUS\n",
  "========================================\n"
)


cat(
  paste0(
    "\nThe full TIDE analysis was executed externally with ",
    "TIDEpy using the NSCLC model and --ignore_norm.\n"
  )
)


cat(
  "\nExpected TIDE result file:\n"
)


print(
  tide_results_file_s31
)


stopifnot(
  file.exists(
    tide_results_file_s31
  )
)


# ============================================================
# SECTION 5
# Read and validate full TIDE result
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 5: VALIDATE FULL TIDE RESULTS\n",
  "========================================\n"
)


tide_results_raw_s31 <-
  read.delim(
    tide_results_file_s31,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


cat(
  "\nRaw TIDE result dimensions:\n"
)


print(
  dim(
    tide_results_raw_s31
  )
)


required_tide_fields_s31 <-
  c(
    "No benefits",
    "Responder",
    "TIDE",
    "IFNG",
    "MSI Score",
    "CD274",
    "CD8",
    "CTL.flag",
    "Dysfunction",
    "Exclusion",
    "MDSC",
    "CAF",
    "TAM M2",
    "CTL"
  )


missing_tide_fields_s31 <-
  setdiff(
    required_tide_fields_s31,
    colnames(
      tide_results_raw_s31
    )
  )


cat(
  "\nMissing TIDE fields:\n"
)


print(
  missing_tide_fields_s31
)


stopifnot(
  length(
    missing_tide_fields_s31
  ) ==
    0
)


expected_patients_s31 <-
  colnames(
    luad_tide_input_s31
  )


cat(
  "\nExpected TIDE patients:\n"
)


print(
  length(
    expected_patients_s31
  )
)


cat(
  "\nReturned TIDE patients:\n"
)


print(
  nrow(
    tide_results_raw_s31
  )
)


stopifnot(
  nrow(
    tide_results_raw_s31
  ) ==
    517
)


stopifnot(
  sum(
    duplicated(
      rownames(
        tide_results_raw_s31
      )
    )
  ) ==
    0
)


missing_tide_patients_s31 <-
  setdiff(
    expected_patients_s31,
    rownames(
      tide_results_raw_s31
    )
  )


unexpected_tide_patients_s31 <-
  setdiff(
    rownames(
      tide_results_raw_s31
    ),
    expected_patients_s31
  )


cat(
  "\nMissing expected patients:\n"
)


print(
  length(
    missing_tide_patients_s31
  )
)


cat(
  "\nUnexpected patients:\n"
)


print(
  length(
    unexpected_tide_patients_s31
  )
)


stopifnot(
  length(
    missing_tide_patients_s31
  ) ==
    0
)


stopifnot(
  length(
    unexpected_tide_patients_s31
  ) ==
    0
)


tide_results_ordered_s31 <-
  tide_results_raw_s31[
    expected_patients_s31,
    required_tide_fields_s31,
    drop = FALSE
  ]


stopifnot(
  identical(
    rownames(
      tide_results_ordered_s31
    ),
    expected_patients_s31
  )
)


# ------------------------------------------------------------
# FIX:
# Convert TIDEpy text boolean fields to R logical values.
#
# TIDEpy output contains:
# "True" / "False"
#
# Downstream R analyses require:
# TRUE / FALSE
# ------------------------------------------------------------


convert_tide_boolean_s31 <-
  function(x) {
    
    x_character <-
      trimws(
        as.character(
          x
        )
      )
    
    
    x_lower <-
      tolower(
        x_character
      )
    
    
    valid_values <-
      c(
        "true",
        "false"
      )
    
    
    invalid_values_s31 <-
      !x_lower %in%
      valid_values
    
    
    if (
      any(
        invalid_values_s31
      )
    ) {
      
      stop(
        paste0(
          "Unexpected TIDE boolean value(s): ",
          paste(
            unique(
              x_character[
                invalid_values_s31
              ]
            ),
            collapse = ", "
          )
        )
      )
    }
    
    
    x_lower ==
      "true"
  }


tide_results_ordered_s31$Responder <-
  convert_tide_boolean_s31(
    tide_results_ordered_s31$Responder
  )


tide_results_ordered_s31[["No benefits"]] <-
  convert_tide_boolean_s31(
    tide_results_ordered_s31[["No benefits"]]
  )


tide_results_ordered_s31[["CTL.flag"]] <-
  convert_tide_boolean_s31(
    tide_results_ordered_s31[["CTL.flag"]]
  )


cat(
  "\nTIDE boolean field classes after conversion:\n"
)


print(
  c(
    Responder =
      class(
        tide_results_ordered_s31$Responder
      ),
    
    No_benefits =
      class(
        tide_results_ordered_s31[["No benefits"]]
      ),
    
    CTL_flag =
      class(
        tide_results_ordered_s31[["CTL.flag"]]
      )
  )
)


stopifnot(
  is.logical(
    tide_results_ordered_s31$Responder
  )
)


stopifnot(
  is.logical(
    tide_results_ordered_s31[["No benefits"]]
  )
)


stopifnot(
  is.logical(
    tide_results_ordered_s31[["CTL.flag"]]
  )
)


# ------------------------------------------------------------
# Validate numeric TIDE fields
# ------------------------------------------------------------


numeric_tide_fields_s31 <-
  c(
    "TIDE",
    "IFNG",
    "MSI Score",
    "CD274",
    "CD8",
    "Dysfunction",
    "Exclusion",
    "MDSC",
    "CAF",
    "TAM M2",
    "CTL"
  )


for (
  field_s31 in
  numeric_tide_fields_s31
) {
  
  stopifnot(
    is.numeric(
      tide_results_ordered_s31[[field_s31]]
    )
  )
  
  
  stopifnot(
    !anyNA(
      tide_results_ordered_s31[[field_s31]]
    )
  )
  
  
  stopifnot(
    all(
      is.finite(
        tide_results_ordered_s31[[field_s31]]
      )
    )
  )
}


cat(
  "\nTIDE summary:\n"
)


print(
  summary(
    tide_results_ordered_s31$TIDE
  )
)


cat(
  "\nPredicted responder distribution:\n"
)


print(
  table(
    tide_results_ordered_s31$Responder
  )
)


cat(
  "\nNo-benefits distribution:\n"
)


print(
  table(
    tide_results_ordered_s31[["No benefits"]]
  )
)


cat(
  "\nCTL.flag distribution:\n"
)


print(
  table(
    tide_results_ordered_s31[["CTL.flag"]]
  )
)


cat(
  "\nTIDE < 0 versus predicted responder:\n"
)


print(
  table(
    TIDE_below_zero =
      tide_results_ordered_s31$TIDE <
      0,
    Responder =
      tide_results_ordered_s31$Responder
  )
)


stopifnot(
  sum(
    tide_results_ordered_s31$Responder
  ) ==
    193
)


stopifnot(
  sum(
    !tide_results_ordered_s31$Responder
  ) ==
    324
)


stopifnot(
  sum(
    tide_results_ordered_s31[["No benefits"]]
  ) ==
    47
)


stopifnot(
  sum(
    tide_results_ordered_s31[["CTL.flag"]]
  ) ==
    141
)


tide_clean_s31 <-
  data.frame(
    patient_id =
      rownames(
        tide_results_ordered_s31
      ),
    tide_results_ordered_s31,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


write.csv(
  tide_clean_s31,
  tide_clean_results_file_s31,
  row.names = FALSE
)


# ============================================================
# SECTION 6
# Read PCD subtype and MPCDS data
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 6: LOAD SUBTYPE AND MPCDS DATA\n",
  "========================================\n"
)


luad_subtypes_s31 <-
  read.csv(
    subtype_file_s31,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


luad_mpcds_s31 <-
  read.csv(
    mpcds_file_s31,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat(
  "\nSubtype columns:\n"
)


print(
  colnames(
    luad_subtypes_s31
  )
)


cat(
  "\nMPCDS columns:\n"
)


print(
  colnames(
    luad_mpcds_s31
  )
)


stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster"
    ) %in%
      colnames(
        luad_subtypes_s31
      )
  )
)


stopifnot(
  "patient_id" %in%
    colnames(
      luad_mpcds_s31
    )
)


stopifnot(
  nrow(
    luad_subtypes_s31
  ) ==
    517
)


stopifnot(
  nrow(
    luad_mpcds_s31
  ) ==
    504
)


cat(
  "\nPCD subtype distribution:\n"
)


print(
  table(
    luad_subtypes_s31$PCD_cluster
  )
)


# ============================================================
# SECTION 7
# Build clean TIDE patient table
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 7: BUILD TIDE PATIENT TABLE\n",
  "========================================\n"
)


tide_patient_s31 <-
  data.frame(
    patient_id =
      rownames(
        tide_results_ordered_s31
      ),
    
    No_benefits =
      tide_results_ordered_s31[["No benefits"]],
    
    Responder =
      tide_results_ordered_s31$Responder,
    
    TIDE =
      tide_results_ordered_s31$TIDE,
    
    IFNG =
      tide_results_ordered_s31$IFNG,
    
    MSI_Score =
      tide_results_ordered_s31[["MSI Score"]],
    
    CD274 =
      tide_results_ordered_s31$CD274,
    
    CD8 =
      tide_results_ordered_s31$CD8,
    
    CTL_flag =
      tide_results_ordered_s31[["CTL.flag"]],
    
    Dysfunction =
      tide_results_ordered_s31$Dysfunction,
    
    Exclusion =
      tide_results_ordered_s31$Exclusion,
    
    MDSC =
      tide_results_ordered_s31$MDSC,
    
    CAF =
      tide_results_ordered_s31$CAF,
    
    TAM_M2 =
      tide_results_ordered_s31[["TAM M2"]],
    
    CTL =
      tide_results_ordered_s31$CTL,
    
    stringsAsFactors = FALSE
  )


stopifnot(
  nrow(
    tide_patient_s31
  ) ==
    517
)


stopifnot(
  is.logical(
    tide_patient_s31$Responder
  )
)


stopifnot(
  is.logical(
    tide_patient_s31$No_benefits
  )
)


stopifnot(
  is.logical(
    tide_patient_s31$CTL_flag
  )
)


# ============================================================
# SECTION 8
# Overall TIDE summary
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 8: OVERALL TIDE SUMMARY\n",
  "========================================\n"
)


overall_tide_summary_s31 <-
  data.frame(
    metric =
      c(
        "Patients",
        "Predicted responders",
        "Predicted nonresponders",
        "Predicted responder percentage",
        "Median TIDE",
        "Mean TIDE",
        "Minimum TIDE",
        "Maximum TIDE"
      ),
    
    value =
      c(
        nrow(
          tide_patient_s31
        ),
        
        sum(
          tide_patient_s31$Responder
        ),
        
        sum(
          !tide_patient_s31$Responder
        ),
        
        100 *
          mean(
            tide_patient_s31$Responder
          ),
        
        median(
          tide_patient_s31$TIDE
        ),
        
        mean(
          tide_patient_s31$TIDE
        ),
        
        min(
          tide_patient_s31$TIDE
        ),
        
        max(
          tide_patient_s31$TIDE
        )
      ),
    
    stringsAsFactors = FALSE
  )


print(
  overall_tide_summary_s31,
  row.names = FALSE
)


# ============================================================
# SECTION 9
# TIDE associations with PCD subtype and MPCDS
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 9: TIDE VS PCD SUBTYPE AND MPCDS\n",
  "========================================\n"
)


# ------------------------------------------------------------
# Identify MPCDS score safely
# ------------------------------------------------------------


candidate_mpcds_columns_s31 <-
  setdiff(
    colnames(
      luad_mpcds_s31
    ),
    "patient_id"
  )


numeric_candidate_mpcds_s31 <-
  candidate_mpcds_columns_s31[
    vapply(
      luad_mpcds_s31[
        ,
        candidate_mpcds_columns_s31,
        drop = FALSE
      ],
      is.numeric,
      logical(1)
    )
  ]


preferred_mpcds_names_s31 <-
  c(
    "MPCDS",
    "MPCDS_score",
    "mpcds",
    "mpcds_score",
    "risk_score",
    "RiskScore",
    "score"
  )


matched_mpcds_name_s31 <-
  intersect(
    preferred_mpcds_names_s31,
    numeric_candidate_mpcds_s31
  )


if (
  length(
    matched_mpcds_name_s31
  ) ==
  1
) {
  
  mpcds_score_col_s31 <-
    matched_mpcds_name_s31
  
} else {
  
  stop(
    paste0(
      "Unable to uniquely identify MPCDS score column. ",
      "Preferred numeric matches: ",
      paste(
        matched_mpcds_name_s31,
        collapse = ", "
      )
    )
  )
}


cat(
  "\nSelected MPCDS score column:\n"
)


print(
  mpcds_score_col_s31
)


# ------------------------------------------------------------
# PCD subtype merge
# ------------------------------------------------------------


tide_subtype_s31 <-
  merge(
    tide_patient_s31,
    luad_subtypes_s31[
      ,
      c(
        "patient_id",
        "PCD_cluster"
      ),
      drop = FALSE
    ],
    by = "patient_id",
    all = FALSE
  )


stopifnot(
  nrow(
    tide_subtype_s31
  ) ==
    517
)


cat(
  "\nSubtype distribution:\n"
)


print(
  table(
    tide_subtype_s31$PCD_cluster
  )
)


# ------------------------------------------------------------
# Continuous TIDE score by subtype
# ------------------------------------------------------------


tide_subtype_wilcox_s31 <-
  wilcox.test(
    TIDE ~
      PCD_cluster,
    data = tide_subtype_s31,
    exact = FALSE
  )


cat(
  "\nWilcoxon test: TIDE by PCD subtype\n"
)


print(
  tide_subtype_wilcox_s31
)


# ------------------------------------------------------------
# Predicted responder by subtype
# ------------------------------------------------------------


tide_responder_table_s31 <-
  table(
    Subtype =
      tide_subtype_s31$PCD_cluster,
    Predicted_responder =
      tide_subtype_s31$Responder
  )


cat(
  "\nPredicted responder counts:\n"
)


print(
  tide_responder_table_s31
)


cat(
  "\nPredicted responder percentages:\n"
)


print(
  round(
    prop.table(
      tide_responder_table_s31,
      margin = 1
    ) *
      100,
    2
  )
)


tide_responder_fisher_s31 <-
  fisher.test(
    tide_responder_table_s31
  )


cat(
  "\nFisher test: responder by subtype\n"
)


print(
  tide_responder_fisher_s31
)


# ------------------------------------------------------------
# MPCDS merge
# ------------------------------------------------------------


mpcds_for_tide_s31 <-
  data.frame(
    patient_id =
      luad_mpcds_s31$patient_id,
    
    MPCDS =
      luad_mpcds_s31[[mpcds_score_col_s31]],
    
    stringsAsFactors = FALSE
  )


stopifnot(
  !anyNA(
    mpcds_for_tide_s31$MPCDS
  )
)


stopifnot(
  all(
    is.finite(
      mpcds_for_tide_s31$MPCDS
    )
  )
)


tide_mpcds_s31 <-
  merge(
    tide_patient_s31,
    mpcds_for_tide_s31,
    by = "patient_id",
    all = FALSE
  )


stopifnot(
  nrow(
    tide_mpcds_s31
  ) ==
    504
)


tide_mpcds_spearman_s31 <-
  cor.test(
    tide_mpcds_s31$MPCDS,
    tide_mpcds_s31$TIDE,
    method = "spearman",
    exact = FALSE
  )


cat(
  "\nSpearman correlation: MPCDS vs TIDE\n"
)


print(
  tide_mpcds_spearman_s31
)


# ------------------------------------------------------------
# MPCDS correlations with TIDE components
# ------------------------------------------------------------


tide_component_names_s31 <-
  c(
    "TIDE",
    "IFNG",
    "MSI_Score",
    "CD274",
    "CD8",
    "Dysfunction",
    "Exclusion",
    "MDSC",
    "CAF",
    "TAM_M2",
    "CTL"
  )


tide_component_correlations_s31 <-
  do.call(
    rbind,
    lapply(
      tide_component_names_s31,
      function(component_name) {
        
        component_values <-
          tide_mpcds_s31[[component_name]]
        
        
        complete_rows <-
          complete.cases(
            tide_mpcds_s31$MPCDS,
            component_values
          )
        
        
        test_result <-
          cor.test(
            tide_mpcds_s31$MPCDS[
              complete_rows
            ],
            component_values[
              complete_rows
            ],
            method = "spearman",
            exact = FALSE
          )
        
        
        data.frame(
          component =
            component_name,
          
          n =
            sum(
              complete_rows
            ),
          
          rho =
            unname(
              test_result$estimate
            ),
          
          p_value =
            test_result$p.value,
          
          stringsAsFactors = FALSE
        )
      }
    )
  )


tide_component_correlations_s31$FDR <-
  p.adjust(
    tide_component_correlations_s31$p_value,
    method = "BH"
  )


cat(
  "\nMPCDS correlations with TIDE components:\n"
)


print(
  tide_component_correlations_s31[
    order(
      tide_component_correlations_s31$FDR
    ),
    ,
    drop = FALSE
  ],
  row.names = FALSE
)


write.csv(
  tide_subtype_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_with_PCD_subtype.csv"
  ),
  row.names = FALSE
)


write.csv(
  tide_mpcds_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_with_MPCDS.csv"
  ),
  row.names = FALSE
)


write.csv(
  tide_component_correlations_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_MPCDS_TIDE_component_correlations.csv"
  ),
  row.names = FALSE
)


# ============================================================
# SECTION 10
# TIDE-component profiles across PCD subtypes
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 10: TIDE COMPONENTS BY PCD SUBTYPE\n",
  "========================================\n"
)


subtype_tide_components_s31 <-
  c(
    "TIDE",
    "IFNG",
    "MSI_Score",
    "CD274",
    "CD8",
    "Dysfunction",
    "Exclusion",
    "MDSC",
    "CAF",
    "TAM_M2",
    "CTL"
  )


summarize_subtype_component_s31 <-
  function(component_name) {
    
    x_c1 <-
      tide_subtype_s31[
        tide_subtype_s31$PCD_cluster ==
          "PCD_C1",
        component_name
      ]
    
    
    x_c2 <-
      tide_subtype_s31[
        tide_subtype_s31$PCD_cluster ==
          "PCD_C2",
        component_name
      ]
    
    
    x_c1 <-
      x_c1[
        is.finite(
          x_c1
        )
      ]
    
    
    x_c2 <-
      x_c2[
        is.finite(
          x_c2
        )
      ]
    
    
    test_result <-
      wilcox.test(
        x_c1,
        x_c2,
        exact = FALSE
      )
    
    
    data.frame(
      component =
        component_name,
      
      n_C1 =
        length(
          x_c1
        ),
      
      median_C1 =
        median(
          x_c1
        ),
      
      q1_C1 =
        unname(
          quantile(
            x_c1,
            0.25
          )
        ),
      
      q3_C1 =
        unname(
          quantile(
            x_c1,
            0.75
          )
        ),
      
      n_C2 =
        length(
          x_c2
        ),
      
      median_C2 =
        median(
          x_c2
        ),
      
      q1_C2 =
        unname(
          quantile(
            x_c2,
            0.25
          )
        ),
      
      q3_C2 =
        unname(
          quantile(
            x_c2,
            0.75
          )
        ),
      
      median_difference_C2_minus_C1 =
        median(
          x_c2
        ) -
        median(
          x_c1
        ),
      
      W =
        unname(
          test_result$statistic
        ),
      
      p_value =
        test_result$p.value,
      
      stringsAsFactors = FALSE
    )
  }


tide_subtype_component_tests_s31 <-
  do.call(
    rbind,
    lapply(
      subtype_tide_components_s31,
      summarize_subtype_component_s31
    )
  )


tide_subtype_component_tests_s31$FDR <-
  p.adjust(
    tide_subtype_component_tests_s31$p_value,
    method = "BH"
  )


tide_subtype_component_tests_s31$significant_FDR_0.05 <-
  tide_subtype_component_tests_s31$FDR <
  0.05


cat(
  "\nTIDE component comparisons:\n"
)


print(
  tide_subtype_component_tests_s31[
    order(
      tide_subtype_component_tests_s31$FDR
    ),
    ,
    drop = FALSE
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# CTL flag by subtype
# ------------------------------------------------------------


ctl_flag_table_s31 <-
  table(
    Subtype =
      tide_subtype_s31$PCD_cluster,
    CTL_flag =
      tide_subtype_s31$CTL_flag
  )


ctl_flag_fisher_s31 <-
  fisher.test(
    ctl_flag_table_s31
  )


cat(
  "\nCTL.flag by subtype:\n"
)


print(
  ctl_flag_table_s31
)


cat(
  "\nCTL.flag row percentages:\n"
)


print(
  round(
    prop.table(
      ctl_flag_table_s31,
      margin = 1
    ) *
      100,
    2
  )
)


print(
  ctl_flag_fisher_s31
)


# ------------------------------------------------------------
# No-benefits flag by subtype
# ------------------------------------------------------------


no_benefits_table_s31 <-
  table(
    Subtype =
      tide_subtype_s31$PCD_cluster,
    No_benefits =
      tide_subtype_s31$No_benefits
  )


no_benefits_fisher_s31 <-
  fisher.test(
    no_benefits_table_s31
  )


cat(
  "\nNo-benefits flag by subtype:\n"
)


print(
  no_benefits_table_s31
)


cat(
  "\nNo-benefits row percentages:\n"
)


print(
  round(
    prop.table(
      no_benefits_table_s31,
      margin = 1
    ) *
      100,
    2
  )
)


print(
  no_benefits_fisher_s31
)


write.csv(
  tide_subtype_component_tests_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_PCD_subtype_TIDE_component_comparisons.csv"
  ),
  row.names = FALSE
)


write.csv(
  as.data.frame(
    ctl_flag_table_s31
  ),
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_PCD_subtype_CTL_flag_counts.csv"
  ),
  row.names = FALSE
)


write.csv(
  as.data.frame(
    no_benefits_table_s31
  ),
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_PCD_subtype_no_benefits_counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# SECTION 11
# Load independent Script 30 immunotherapy indicators
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 11: LOAD INDEPENDENT INDICATORS\n",
  "========================================\n"
)


cd274_file_s31 <-
  file.path(
    script26_indicator_dir_s31,
    "TCGA_LUAD_patient_CD274_with_PCD_subtype.csv"
  )


mantis_file_s31 <-
  file.path(
    script26_indicator_dir_s31,
    "TCGA_LUAD_patient_MANTIS_with_PCD_subtype.csv"
  )


mutation_file_s31 <-
  file.path(
    script26_indicator_dir_s31,
    "TCGA_LUAD_patient_mutation_burden_with_PCD_subtype.csv"
  )


independent_indicator_files_s31 <-
  c(
    cd274_file_s31,
    mantis_file_s31,
    mutation_file_s31
  )


cat(
  "\nIndependent indicator files:\n"
)


print(
  data.frame(
    file =
      independent_indicator_files_s31,
    
    exists =
      file.exists(
        independent_indicator_files_s31
      ),
    
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


stopifnot(
  all(
    file.exists(
      independent_indicator_files_s31
    )
  )
)


independent_cd274_s31 <-
  read.csv(
    cd274_file_s31,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


independent_mantis_s31 <-
  read.csv(
    mantis_file_s31,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


independent_mutation_s31 <-
  read.csv(
    mutation_file_s31,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


cat(
  "\nIndependent indicator dimensions:\n"
)


print(
  c(
    CD274 =
      nrow(
        independent_cd274_s31
      ),
    
    MANTIS =
      nrow(
        independent_mantis_s31
      ),
    
    Mutation =
      nrow(
        independent_mutation_s31
      )
  )
)


stopifnot(
  all(
    c(
      "patient_id",
      "CD274_TPM",
      "CD274_log2TPM",
      "PCD_cluster"
    ) %in%
      colnames(
        independent_cd274_s31
      )
  )
)


stopifnot(
  all(
    c(
      "patient_id",
      "MANTIS_score",
      "MSI_class",
      "PCD_cluster"
    ) %in%
      colnames(
        independent_mantis_s31
      )
  )
)


stopifnot(
  all(
    c(
      "patient_id",
      "n_nonsyn_mutations",
      "log2_nonsyn_mutations",
      "PCD_cluster"
    ) %in%
      colnames(
        independent_mutation_s31
      )
  )
)


stopifnot(
  !any(
    duplicated(
      independent_cd274_s31$patient_id
    )
  )
)


stopifnot(
  !any(
    duplicated(
      independent_mantis_s31$patient_id
    )
  )
)


stopifnot(
  !any(
    duplicated(
      independent_mutation_s31$patient_id
    )
  )
)


# ============================================================
# SECTION 12
# Integrate TIDE with independent indicators
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 12: TIDE + INDEPENDENT BIOMARKERS\n",
  "========================================\n"
)


# ------------------------------------------------------------
# Build minimal tables
# ------------------------------------------------------------


cd274_minimal_s31 <-
  independent_cd274_s31[
    ,
    c(
      "patient_id",
      "CD274_TPM",
      "CD274_log2TPM"
    ),
    drop = FALSE
  ]


names(
  cd274_minimal_s31
)[
  names(
    cd274_minimal_s31
  ) ==
    "CD274_log2TPM"
] <-
  "Independent_CD274_log2TPM"


mantis_minimal_s31 <-
  independent_mantis_s31[
    ,
    c(
      "patient_id",
      "MANTIS_score",
      "MSI_class"
    ),
    drop = FALSE
  ]


mutation_minimal_s31 <-
  independent_mutation_s31[
    ,
    c(
      "patient_id",
      "n_nonsyn_mutations",
      "log2_nonsyn_mutations"
    ),
    drop = FALSE
  ]


names(
  mutation_minimal_s31
)[
  names(
    mutation_minimal_s31
  ) ==
    "log2_nonsyn_mutations"
] <-
  "Independent_log2_nonsyn_mutations"


# ------------------------------------------------------------
# TIDE + independent CD274
# ------------------------------------------------------------


tide_cd274_independent_s31 <-
  merge(
    tide_patient_s31,
    cd274_minimal_s31,
    by = "patient_id",
    all = FALSE
  )


stopifnot(
  nrow(
    tide_cd274_independent_s31
  ) ==
    517
)


tide_vs_cd274_s31 <-
  cor.test(
    tide_cd274_independent_s31$TIDE,
    tide_cd274_independent_s31$
      Independent_CD274_log2TPM,
    method = "spearman",
    exact = FALSE
  )


tide_internal_vs_independent_cd274_s31 <-
  cor.test(
    tide_cd274_independent_s31$CD274,
    tide_cd274_independent_s31$
      Independent_CD274_log2TPM,
    method = "spearman",
    exact = FALSE
  )


cat(
  "\nTIDE vs independent CD274:\n"
)


print(
  tide_vs_cd274_s31
)


cat(
  "\nTIDE-derived CD274 vs independent CD274:\n"
)


print(
  tide_internal_vs_independent_cd274_s31
)


# ------------------------------------------------------------
# TIDE + mutation burden
# ------------------------------------------------------------


tide_mutation_independent_s31 <-
  merge(
    tide_patient_s31,
    mutation_minimal_s31,
    by = "patient_id",
    all = FALSE
  )


stopifnot(
  nrow(
    tide_mutation_independent_s31
  ) ==
    505
)


tide_vs_mutation_s31 <-
  cor.test(
    tide_mutation_independent_s31$TIDE,
    tide_mutation_independent_s31$
      Independent_log2_nonsyn_mutations,
    method = "spearman",
    exact = FALSE
  )


cat(
  "\nTIDE vs nonsynonymous mutation burden:\n"
)


print(
  tide_vs_mutation_s31
)


# ------------------------------------------------------------
# TIDE + MANTIS
# ------------------------------------------------------------


tide_mantis_independent_s31 <-
  merge(
    tide_patient_s31,
    mantis_minimal_s31,
    by = "patient_id",
    all = FALSE
  )


stopifnot(
  nrow(
    tide_mantis_independent_s31
  ) ==
    514
)


tide_vs_mantis_s31 <-
  cor.test(
    tide_mantis_independent_s31$TIDE,
    tide_mantis_independent_s31$MANTIS_score,
    method = "spearman",
    exact = FALSE
  )


tide_msi_vs_mantis_s31 <-
  cor.test(
    tide_mantis_independent_s31$MSI_Score,
    tide_mantis_independent_s31$MANTIS_score,
    method = "spearman",
    exact = FALSE
  )


cat(
  "\nTIDE vs independent MANTIS:\n"
)


print(
  tide_vs_mantis_s31
)


cat(
  "\nTIDE MSI Score vs independent MANTIS:\n"
)


print(
  tide_msi_vs_mantis_s31
)


# ------------------------------------------------------------
# Predicted responder comparisons
# ------------------------------------------------------------


responder_cd274_wilcox_s31 <-
  wilcox.test(
    Independent_CD274_log2TPM ~
      Responder,
    data =
      tide_cd274_independent_s31,
    exact = FALSE
  )


responder_mutation_wilcox_s31 <-
  wilcox.test(
    Independent_log2_nonsyn_mutations ~
      Responder,
    data =
      tide_mutation_independent_s31,
    exact = FALSE
  )


responder_mantis_wilcox_s31 <-
  wilcox.test(
    MANTIS_score ~
      Responder,
    data =
      tide_mantis_independent_s31,
    exact = FALSE
  )


cat(
  "\nCD274 median by predicted responder:\n"
)


print(
  aggregate(
    Independent_CD274_log2TPM ~
      Responder,
    data =
      tide_cd274_independent_s31,
    FUN = median
  )
)


print(
  responder_cd274_wilcox_s31
)


cat(
  "\nMutation burden median by predicted responder:\n"
)


print(
  aggregate(
    Independent_log2_nonsyn_mutations ~
      Responder,
    data =
      tide_mutation_independent_s31,
    FUN = median
  )
)


print(
  responder_mutation_wilcox_s31
)


cat(
  "\nMANTIS median by predicted responder:\n"
)


print(
  aggregate(
    MANTIS_score ~
      Responder,
    data =
      tide_mantis_independent_s31,
    FUN = median
  )
)


print(
  responder_mantis_wilcox_s31
)


# ------------------------------------------------------------
# Compact correlation summary
# ------------------------------------------------------------


tide_independent_correlations_s31 <-
  data.frame(
    comparison =
      c(
        "TIDE_vs_independent_CD274",
        "TIDE_CD274_vs_independent_CD274",
        "TIDE_vs_nonsynonymous_mutation_burden",
        "TIDE_vs_MANTIS",
        "TIDE_MSI_Score_vs_MANTIS"
      ),
    
    n =
      c(
        nrow(
          tide_cd274_independent_s31
        ),
        
        nrow(
          tide_cd274_independent_s31
        ),
        
        nrow(
          tide_mutation_independent_s31
        ),
        
        nrow(
          tide_mantis_independent_s31
        ),
        
        nrow(
          tide_mantis_independent_s31
        )
      ),
    
    rho =
      c(
        unname(
          tide_vs_cd274_s31$estimate
        ),
        
        unname(
          tide_internal_vs_independent_cd274_s31$estimate
        ),
        
        unname(
          tide_vs_mutation_s31$estimate
        ),
        
        unname(
          tide_vs_mantis_s31$estimate
        ),
        
        unname(
          tide_msi_vs_mantis_s31$estimate
        )
      ),
    
    p_value =
      c(
        tide_vs_cd274_s31$p.value,
        tide_internal_vs_independent_cd274_s31$p.value,
        tide_vs_mutation_s31$p.value,
        tide_vs_mantis_s31$p.value,
        tide_msi_vs_mantis_s31$p.value
      ),
    
    stringsAsFactors = FALSE
  )


tide_independent_correlations_s31$FDR <-
  p.adjust(
    tide_independent_correlations_s31$p_value,
    method = "BH"
  )


cat(
  "\nTIDE-independent correlation summary:\n"
)


print(
  tide_independent_correlations_s31,
  row.names = FALSE
)


write.csv(
  tide_independent_correlations_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_independent_indicator_correlations.csv"
  ),
  row.names = FALSE
)


write.csv(
  tide_cd274_independent_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_with_independent_CD274.csv"
  ),
  row.names = FALSE
)


write.csv(
  tide_mutation_independent_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_with_nonsynonymous_mutation_burden.csv"
  ),
  row.names = FALSE
)


write.csv(
  tide_mantis_independent_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_with_MANTIS.csv"
  ),
  row.names = FALSE
)


# ============================================================
# SECTION 13
# IPS feasibility
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 13: IPS FEASIBILITY\n",
  "========================================\n"
)


all_project_files_ips_s31 <-
  list.files(
    ".",
    recursive = TRUE,
    full.names = TRUE
  )


ips_pattern_s31 <-
  paste(
    c(
      "IPS",
      "Immunophenoscore",
      "Immunophenotype",
      "TCIA",
      "Cancer_Immunome",
      "CancerImmunome",
      "CTLA4",
      "PD1",
      "PD-1"
    ),
    collapse = "|"
  )


ips_candidate_files_s31 <-
  all_project_files_ips_s31[
    grepl(
      ips_pattern_s31,
      all_project_files_ips_s31,
      ignore.case = TRUE
    )
  ]


ips_tabular_files_s31 <-
  ips_candidate_files_s31[
    grepl(
      "\\.(csv|tsv|txt|xlsx|xls|rds)$",
      ips_candidate_files_s31,
      ignore.case = TRUE
    )
  ]


# ------------------------------------------------------------
# Exclude Script 27's own generated IPS-named files.
#
# Otherwise a second run could incorrectly interpret its own
# feasibility output as an actual IPS dataset.
# ------------------------------------------------------------


ips_output_exclusion_pattern_s31 <-
  paste(
    c(
      "TIDE_IPS_feasibility",
      "TIDE_IPS_interpretation",
      "Script27_analysis_audit"
    ),
    collapse = "|"
  )


ips_tabular_files_s31 <-
  ips_tabular_files_s31[
    !grepl(
      ips_output_exclusion_pattern_s31,
      ips_tabular_files_s31,
      ignore.case = TRUE
    )
  ]


project_csv_files_ips_s31 <-
  all_project_files_ips_s31[
    grepl(
      "\\.csv$",
      all_project_files_ips_s31,
      ignore.case = TRUE
    )
  ]


project_csv_files_ips_s31 <-
  project_csv_files_ips_s31[
    !grepl(
      ips_output_exclusion_pattern_s31,
      project_csv_files_ips_s31,
      ignore.case = TRUE
    )
  ]


ips_header_hits_s31 <-
  list()


if (
  length(
    project_csv_files_ips_s31
  ) >
  0
) {
  
  for (
    current_file_ips_s31 in
    project_csv_files_ips_s31
  ) {
    
    current_header_ips_s31 <-
      tryCatch(
        read.csv(
          current_file_ips_s31,
          nrows = 1,
          check.names = FALSE,
          stringsAsFactors = FALSE
        ),
        error =
          function(e) {
            
            NULL
          }
      )
    
    
    if (
      !is.null(
        current_header_ips_s31
      )
    ) {
      
      current_columns_ips_s31 <-
        colnames(
          current_header_ips_s31
        )
      
      
      header_match_ips_s31 <-
        grepl(
          paste(
            c(
              "^IPS$",
              "IPS_",
              "_IPS",
              "Immunophenoscore",
              "ips_ctla4",
              "ips_pd1",
              "ctla4_neg_pd1_neg",
              "ctla4_pos_pd1_neg",
              "ctla4_neg_pd1_pos",
              "ctla4_pos_pd1_pos"
            ),
            collapse = "|"
          ),
          current_columns_ips_s31,
          ignore.case = TRUE
        )
      
      
      if (
        any(
          header_match_ips_s31
        )
      ) {
        
        ips_header_hits_s31[[current_file_ips_s31]] <-
          current_columns_ips_s31[
            header_match_ips_s31
          ]
      }
    }
  }
}


tide_feasible_s31 <-
  nrow(
    tide_results_ordered_s31
  ) ==
  517


ips_local_available_s31 <-
  length(
    ips_tabular_files_s31
  ) >
  0 ||
  length(
    ips_header_hits_s31
  ) >
  0


cat(
  "\nPossible independent/local IPS files:\n"
)


if (
  length(
    ips_tabular_files_s31
  ) ==
  0
) {
  
  cat(
    "None\n"
  )
  
} else {
  
  print(
    ips_tabular_files_s31
  )
}


cat(
  "\nIPS-like columns found:\n"
)


if (
  length(
    ips_header_hits_s31
  ) ==
  0
) {
  
  cat(
    "None\n"
  )
  
} else {
  
  print(
    ips_header_hits_s31
  )
}


cat(
  "\nTIDE completed:\n"
)


print(
  tide_feasible_s31
)


cat(
  "\nIndependent/local IPS candidate available:\n"
)


print(
  ips_local_available_s31
)


ips_feasibility_s31 <-
  data.frame(
    framework =
      c(
        "TIDE",
        "IPS"
      ),
    
    status =
      c(
        ifelse(
          tide_feasible_s31,
          "Completed",
          "Not completed"
        ),
        
        ifelse(
          ips_local_available_s31,
          "Local candidate data found",
          "No local candidate data found"
        )
      ),
    
    interpretation =
      c(
        paste0(
          "Expression-based in silico immunotherapy surrogate ",
          "completed for TCGA-LUAD."
        ),
        
        ifelse(
          ips_local_available_s31,
          paste0(
            "Candidate local IPS-related data require ",
            "validation before analysis."
          ),
          paste0(
            "No existing IPS dataset detected in the project; ",
            "external Cancer Immunome Atlas access would be ",
            "required for direct IPS analysis."
          )
        )
      ),
    
    stringsAsFactors = FALSE
  )


cat(
  "\nTIDE / IPS feasibility summary:\n"
)


print(
  ips_feasibility_s31,
  row.names = FALSE
)


write.csv(
  ips_feasibility_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_IPS_feasibility.csv"
  ),
  row.names = FALSE
)


# ============================================================
# SECTION 14
# Final summaries and colored Figure S14
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 14: FINAL SUMMARY AND FIGURE S14\n",
  "========================================\n"
)


# ------------------------------------------------------------
# Final PCD subtype summary
# ------------------------------------------------------------


pcd_c1_tide_s31 <-
  tide_subtype_s31$TIDE[
    tide_subtype_s31$PCD_cluster ==
      "PCD_C1"
  ]


pcd_c2_tide_s31 <-
  tide_subtype_s31$TIDE[
    tide_subtype_s31$PCD_cluster ==
      "PCD_C2"
  ]


subtype_final_summary_s31 <-
  data.frame(
    metric =
      c(
        "n_PCD_C1",
        "n_PCD_C2",
        "Median_TIDE_PCD_C1",
        "Median_TIDE_PCD_C2",
        "TIDE_Wilcoxon_p",
        "Predicted_responder_percent_PCD_C1",
        "Predicted_responder_percent_PCD_C2",
        "Responder_Fisher_p",
        "No_benefits_percent_PCD_C1",
        "No_benefits_percent_PCD_C2",
        "No_benefits_Fisher_p"
      ),
    
    value =
      c(
        sum(
          tide_subtype_s31$PCD_cluster ==
            "PCD_C1"
        ),
        
        sum(
          tide_subtype_s31$PCD_cluster ==
            "PCD_C2"
        ),
        
        median(
          pcd_c1_tide_s31
        ),
        
        median(
          pcd_c2_tide_s31
        ),
        
        tide_subtype_wilcox_s31$p.value,
        
        100 *
          tide_responder_table_s31[
            "PCD_C1",
            "TRUE"
          ] /
          sum(
            tide_responder_table_s31[
              "PCD_C1",
            ]
          ),
        
        100 *
          tide_responder_table_s31[
            "PCD_C2",
            "TRUE"
          ] /
          sum(
            tide_responder_table_s31[
              "PCD_C2",
            ]
          ),
        
        tide_responder_fisher_s31$p.value,
        
        100 *
          no_benefits_table_s31[
            "PCD_C1",
            "TRUE"
          ] /
          sum(
            no_benefits_table_s31[
              "PCD_C1",
            ]
          ),
        
        100 *
          no_benefits_table_s31[
            "PCD_C2",
            "TRUE"
          ] /
          sum(
            no_benefits_table_s31[
              "PCD_C2",
            ]
          ),
        
        no_benefits_fisher_s31$p.value
      ),
    
    stringsAsFactors = FALSE
  )


cat(
  "\nPCD subtype final summary:\n"
)


print(
  subtype_final_summary_s31,
  row.names = FALSE
)


# ------------------------------------------------------------
# MPCDS final summary
# ------------------------------------------------------------


mpcds_tide_final_summary_s31 <-
  data.frame(
    metric =
      c(
        "MPCDS_TIDE_n",
        "MPCDS_TIDE_Spearman_rho",
        "MPCDS_TIDE_p"
      ),
    
    value =
      c(
        nrow(
          tide_mpcds_s31
        ),
        
        unname(
          tide_mpcds_spearman_s31$estimate
        ),
        
        tide_mpcds_spearman_s31$p.value
      ),
    
    stringsAsFactors = FALSE
  )


cat(
  "\nMPCDS-TIDE final summary:\n"
)


print(
  mpcds_tide_final_summary_s31,
  row.names = FALSE
)


# ------------------------------------------------------------
# Significant MPCDS-associated TIDE components
# ------------------------------------------------------------


significant_mpcds_components_s31 <-
  tide_component_correlations_s31[
    tide_component_correlations_s31$FDR <
      0.05,
    ,
    drop = FALSE
  ]


significant_mpcds_components_s31 <-
  significant_mpcds_components_s31[
    order(
      significant_mpcds_components_s31$FDR
    ),
    ,
    drop = FALSE
  ]


cat(
  "\nSignificant MPCDS-associated TIDE components:\n"
)


print(
  significant_mpcds_components_s31,
  row.names = FALSE
)


# ------------------------------------------------------------
# Significant subtype-differential components
# ------------------------------------------------------------


significant_subtype_components_final_s31 <-
  tide_subtype_component_tests_s31[
    tide_subtype_component_tests_s31$FDR <
      0.05,
    ,
    drop = FALSE
  ]


significant_subtype_components_final_s31 <-
  significant_subtype_components_final_s31[
    order(
      significant_subtype_components_final_s31$FDR
    ),
    ,
    drop = FALSE
  ]


cat(
  "\nSignificant subtype-differential TIDE components:\n"
)


print(
  significant_subtype_components_final_s31,
  row.names = FALSE
)


# ------------------------------------------------------------
# Interpretation summary
# ------------------------------------------------------------


interpretation_summary_s31 <-
  data.frame(
    topic =
      c(
        "Overall TIDE feasibility",
        "PCD subtype overall TIDE",
        "PCD subtype predicted responder",
        "PCD subtype immune contexture",
        "MPCDS and TIDE",
        "MPCDS and immune exclusion",
        "Independent CD274",
        "Independent mutation burden",
        "Independent MANTIS",
        "IPS feasibility",
        "Clinical interpretation"
      ),
    
    conclusion =
      c(
        paste0(
          "TIDE analysis was technically feasible and completed ",
          "for 517 TCGA-LUAD patients."
        ),
        
        paste0(
          "PCD_C1 and PCD_C2 did not differ significantly in ",
          "overall continuous TIDE score."
        ),
        
        paste0(
          "PCD_C1 and PCD_C2 did not differ significantly in ",
          "TIDE-predicted responder frequency."
        ),
        
        paste0(
          "The subtypes showed distinct TIDE-derived immune ",
          "contexture profiles despite similar overall TIDE scores."
        ),
        
        paste0(
          "Higher MPCDS was associated with higher TIDE score, ",
          "indicating a less favorable in silico TIDE profile."
        ),
        
        paste0(
          "The MPCDS association was strongest for TIDE-derived ",
          "MDSC and immune-exclusion components."
        ),
        
        paste0(
          "Overall TIDE was not significantly associated with ",
          "independently extracted CD274 expression."
        ),
        
        paste0(
          "Overall TIDE showed a weak positive association with ",
          "nonsynonymous mutation burden."
        ),
        
        paste0(
          "Overall TIDE was not significantly associated with ",
          "continuous MANTIS score."
        ),
        
        paste0(
          "No local IPS dataset was available; direct IPS analysis ",
          "would require external Cancer Immunome Atlas access."
        ),
        
        paste0(
          "All TIDE findings are exploratory in silico ",
          "immunotherapy-response surrogates and do not constitute ",
          "clinical treatment-response validation."
        )
      ),
    
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# Save final tables
# ------------------------------------------------------------


write.csv(
  overall_tide_summary_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_overall_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  subtype_final_summary_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_PCD_subtype_final_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  mpcds_tide_final_summary_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_MPCDS_TIDE_final_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  significant_mpcds_components_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_MPCDS_significant_TIDE_components.csv"
  ),
  row.names = FALSE
)


write.csv(
  significant_subtype_components_final_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_PCD_subtype_significant_TIDE_components.csv"
  ),
  row.names = FALSE
)


write.csv(
  interpretation_summary_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_TIDE_IPS_interpretation_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Figure generation
# ------------------------------------------------------------

# Final manuscript figures are intentionally deferred to Script 32.
# Script 27 produces validated TIDE analysis tables only.


# ------------------------------------------------------------
# Regression checks against established TIDE analysis
# ------------------------------------------------------------

stopifnot(
  nrow(
    tide_patient_s31
  ) ==
    517,
  sum(
    tide_patient_s31$Responder
  ) ==
    193,
  sum(
    !tide_patient_s31$Responder
  ) ==
    324,
  sum(
    tide_patient_s31$No_benefits
  ) ==
    47,
  abs(
    unname(
      tide_mpcds_spearman_s31$estimate
    ) -
      0.3540614
  ) <
    1e-5,
  abs(
    tide_mpcds_spearman_s31$p.value -
      2.487166e-16
  ) <
    1e-18,
  abs(
    100 *
      no_benefits_table_s31[
        "PCD_C1",
        "TRUE"
      ] /
      sum(
        no_benefits_table_s31[
          "PCD_C1",
        ]
      ) -
      4.58
  ) <
    0.02,
  abs(
    100 *
      no_benefits_table_s31[
        "PCD_C2",
        "TRUE"
      ] /
      sum(
        no_benefits_table_s31[
          "PCD_C2",
        ]
      ) -
      14.59
  ) <
    0.02,
  abs(
    no_benefits_fisher_s31$p.value -
      9.392e-05
  ) <
    1e-06,
  abs(
    unname(
      no_benefits_fisher_s31$estimate
    ) -
      3.553
  ) <
    0.02
)

cat(
  "\nEstablished TIDE regression checks: PASSED\n"
)


# ============================================================
# SECTION 15
# Final scientific audit
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 15: FINAL SCIENTIFIC AUDIT\n",
  "========================================\n"
)


script31_audit_s31 <-
  data.frame(
    item =
      c(
        "TIDE technical execution",
        "Full TCGA-LUAD cohort",
        "PCD subtype comparison",
        "MPCDS association",
        "TIDE component analysis",
        "Independent CD274 integration",
        "Independent mutation burden integration",
        "Independent MANTIS integration",
        "IPS feasibility assessment",
        "Clinical-response validation"
      ),
    
    status =
      c(
        "Completed",
        "Completed",
        "Completed",
        "Completed",
        "Completed",
        "Completed",
        "Completed",
        "Completed",
        "Completed - external IPS data unavailable locally",
        "Not performed - no ICI-treated clinical cohort"
      ),
    
    stringsAsFactors = FALSE
  )


write.csv(
  script31_audit_s31,
  file.path(
    immunotherapy_results_dir_s31,
    "TCGA_LUAD_Script27_analysis_audit.csv"
  ),
  row.names = FALSE
)


cat(
  "\nFinal Script 27 audit:\n"
)


print(
  script31_audit_s31,
  row.names = FALSE
)


# ============================================================
# FINAL SUCCESS MARKER
# ============================================================


cat(
  "\n========================================\n",
  "SCRIPT 27 TIDE IMMUNOTHERAPY-SURROGATE ANALYSIS COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)