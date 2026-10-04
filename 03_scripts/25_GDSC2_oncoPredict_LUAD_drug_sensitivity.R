# ============================================================
# SCRIPT 25
# GDSC2 / oncoPredict drug-sensitivity analysis
#
# TCGA-LUAD programmed-cell-death project
#
# Main objectives:
# 1. Validate official GDSC2 training matrices
# 2. Predict TCGA-LUAD drug sensitivity using oncoPredict
# 3. Test predicted IC50 associations with MPCDS
# 4. Test predicted IC50 differences by PCD subtype
# 5. Compare MPCDS-drug associations with DepMap/PRISM
# 6. Generate final summaries and Figure S15
#
# Interpretation:
# - Lower predicted GDSC2 IC50 = greater predicted sensitivity
# - Lower PRISM AUC = greater sensitivity
# - All GDSC2/oncoPredict results are in silico predictions
# - PRISM is an independent cell-line pharmacogenomic analysis
# - Cross-platform agreement does not constitute clinical
#   treatment-response validation
# ============================================================


# ============================================================
# SECTION 1
# Environment, paths, packages, and input validation
# ============================================================


source("03_scripts/00_project_config.R")


cat(
  "\n========================================\n",
  "SCRIPT 25: GDSC2 / oncoPredict\n",
  "SECTION 1: ENVIRONMENT AND INPUTS\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 1.1 Stable download/repository options
# ------------------------------------------------------------


options(
  timeout = 1200
)


options(
  download.file.method = "libcurl"
)


# ------------------------------------------------------------
# 1.2 Required packages
# ------------------------------------------------------------


required_packages_s32 <-
  c(
    "oncoPredict",
    "sva",
    "limma",
    "glmnet",
    "ridge",
    "pls"
  )


package_status_s32 <-
  vapply(
    required_packages_s32,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )


cat(
  "\nRequired package status:\n"
)


print(
  package_status_s32
)


if (
  !all(
    package_status_s32
  )
) {
  
  stop(
    paste0(
      "Missing required package(s): ",
      paste(
        names(
          package_status_s32
        )[
          !package_status_s32
        ],
        collapse = ", "
      )
    )
  )
}


cat(
  "\noncoPredict version:\n"
)


print(
  as.character(
    packageVersion(
      "oncoPredict"
    )
  )
)


calcPhenotype_available_s32 <-
  exists(
    "calcPhenotype",
    where = asNamespace(
      "oncoPredict"
    ),
    inherits = FALSE
  )


cat(
  "\ncalcPhenotype available:\n"
)


print(
  calcPhenotype_available_s32
)


stopifnot(
  calcPhenotype_available_s32
)


# ------------------------------------------------------------
# 1.3 Define project inputs
# ------------------------------------------------------------


gdsc2_expr_file_s32 <-
  "01_raw_data/GDSC/GDSC2_Expr (RMA Normalized and Log Transformed).rds"


gdsc2_res_file_s32 <-
  "01_raw_data/GDSC/GDSC2_Res.rds"


patient_tpm_file_s32 <-
  "02_processed_data/immunotherapy/TCGA_LUAD_patient_protein_coding_TPM.rds"


mpcds_file_s32 <-
  "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv"


subtype_file_s32 <-
  "04_results/clustering/TCGA_LUAD_PCD_cluster_assignments.csv"


prism_file_s32 <-
  "04_results/drug_sensitivity/DepMap_PRISM_MPCDS_all_drug_associations.csv"


input_files_s32 <-
  c(
    gdsc2_expr_file_s32,
    gdsc2_res_file_s32,
    patient_tpm_file_s32,
    mpcds_file_s32,
    subtype_file_s32,
    prism_file_s32
  )


cat(
  "\nInput file check:\n"
)


print(
  data.frame(
    file = input_files_s32,
    exists = file.exists(
      input_files_s32
    ),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


stopifnot(
  all(
    file.exists(
      input_files_s32
    )
  )
)


# ------------------------------------------------------------
# 1.4 Define output directories
# ------------------------------------------------------------


gdsc_processed_dir_s32 <-
  "02_processed_data/drug_sensitivity/GDSC2"


prediction_dir_s32 <-
  file.path(
    gdsc_processed_dir_s32,
    "calcPhenotype_Output"
  )


gdsc_results_dir_s32 <-
  "04_results/drug_sensitivity/GDSC2"




dir.create(
  gdsc_processed_dir_s32,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  prediction_dir_s32,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  gdsc_results_dir_s32,
  recursive = TRUE,
  showWarnings = FALSE
)



cat(
  "\n========================================\n",
  "SECTION 1 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 2
# Load and validate official GDSC2 training matrices
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 2: VALIDATE GDSC2 MATRICES\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 2.1 Load matrices
# ------------------------------------------------------------


gdsc2_expr_s32 <-
  readRDS(
    gdsc2_expr_file_s32
  )


gdsc2_res_log_s32 <-
  readRDS(
    gdsc2_res_file_s32
  )


stopifnot(
  is.matrix(
    gdsc2_expr_s32
  ) ||
    is.data.frame(
      gdsc2_expr_s32
    )
)


stopifnot(
  is.matrix(
    gdsc2_res_log_s32
  ) ||
    is.data.frame(
      gdsc2_res_log_s32
    )
)


gdsc2_expr_s32 <-
  as.matrix(
    gdsc2_expr_s32
  )


gdsc2_res_log_s32 <-
  as.matrix(
    gdsc2_res_log_s32
  )


# ------------------------------------------------------------
# 2.2 Expression matrix QC
# ------------------------------------------------------------


cat(
  "\nGDSC2 expression dimensions:\n"
)


print(
  dim(
    gdsc2_expr_s32
  )
)


cat(
  "\nGDSC2 expression first genes:\n"
)


print(
  head(
    rownames(
      gdsc2_expr_s32
    ),
    10
  )
)


cat(
  "\nGDSC2 expression first cell lines:\n"
)


print(
  head(
    colnames(
      gdsc2_expr_s32
    ),
    10
  )
)


cat(
  "\nExpression NA count:\n"
)


print(
  sum(
    is.na(
      gdsc2_expr_s32
    )
  )
)


cat(
  "\nExpression non-finite count:\n"
)


print(
  sum(
    !is.finite(
      gdsc2_expr_s32
    )
  )
)


cat(
  "\nDuplicate expression genes:\n"
)


print(
  sum(
    duplicated(
      rownames(
        gdsc2_expr_s32
      )
    )
  )
)


cat(
  "\nDuplicate expression cell lines:\n"
)


print(
  sum(
    duplicated(
      colnames(
        gdsc2_expr_s32
      )
    )
  )
)


cat(
  "\nExpression summary:\n"
)


print(
  summary(
    as.vector(
      gdsc2_expr_s32
    )
  )
)


stopifnot(
  dim(
    gdsc2_expr_s32
  )[
    1
  ] ==
    17419
)


stopifnot(
  dim(
    gdsc2_expr_s32
  )[
    2
  ] ==
    805
)


stopifnot(
  !anyNA(
    gdsc2_expr_s32
  )
)


stopifnot(
  all(
    is.finite(
      gdsc2_expr_s32
    )
  )
)


stopifnot(
  !any(
    duplicated(
      rownames(
        gdsc2_expr_s32
      )
    )
  )
)


stopifnot(
  !any(
    duplicated(
      colnames(
        gdsc2_expr_s32
      )
    )
  )
)


# ------------------------------------------------------------
# 2.3 Response matrix QC
# ------------------------------------------------------------


cat(
  "\nGDSC2 response dimensions:\n"
)


print(
  dim(
    gdsc2_res_log_s32
  )
)


cat(
  "\nFirst response cell lines:\n"
)


print(
  head(
    rownames(
      gdsc2_res_log_s32
    ),
    10
  )
)


cat(
  "\nFirst GDSC2 drugs:\n"
)


print(
  head(
    colnames(
      gdsc2_res_log_s32
    ),
    10
  )
)


cat(
  "\nResponse NA count:\n"
)


print(
  sum(
    is.na(
      gdsc2_res_log_s32
    )
  )
)


cat(
  "\nResponse log-scale summary:\n"
)


print(
  summary(
    as.vector(
      gdsc2_res_log_s32
    )
  )
)


stopifnot(
  dim(
    gdsc2_res_log_s32
  )[
    1
  ] ==
    805
)


stopifnot(
  dim(
    gdsc2_res_log_s32
  )[
    2
  ] ==
    198
)


stopifnot(
  !any(
    duplicated(
      rownames(
        gdsc2_res_log_s32
      )
    )
  )
)


stopifnot(
  !any(
    duplicated(
      colnames(
        gdsc2_res_log_s32
      )
    )
  )
)


# ------------------------------------------------------------
# 2.4 Confirm cell-line alignment
# ------------------------------------------------------------


gdsc2_common_cell_lines_s32 <-
  intersect(
    colnames(
      gdsc2_expr_s32
    ),
    rownames(
      gdsc2_res_log_s32
    )
  )


cat(
  "\nGDSC2 expression/response common cell lines:\n"
)


print(
  length(
    gdsc2_common_cell_lines_s32
  )
)


stopifnot(
  length(
    gdsc2_common_cell_lines_s32
  ) ==
    805
)


stopifnot(
  identical(
    colnames(
      gdsc2_expr_s32
    ),
    rownames(
      gdsc2_res_log_s32
    )
  )
)


# ------------------------------------------------------------
# 2.5 Exponentiate stored log-scale response
# ------------------------------------------------------------


gdsc2_res_s32 <-
  exp(
    gdsc2_res_log_s32
  )


cat(
  "\nExponentiated GDSC2 response summary:\n"
)


print(
  summary(
    as.vector(
      gdsc2_res_s32
    )
  )
)


# ------------------------------------------------------------
# 2.6 Save validated GDSC2 matrices
# ------------------------------------------------------------


saveRDS(
  gdsc2_expr_s32,
  file.path(
    gdsc_processed_dir_s32,
    "GDSC2_expression_validated.rds"
  )
)


saveRDS(
  gdsc2_res_log_s32,
  file.path(
    gdsc_processed_dir_s32,
    "GDSC2_response_logIC50_validated.rds"
  )
)


saveRDS(
  gdsc2_res_s32,
  file.path(
    gdsc_processed_dir_s32,
    "GDSC2_response_IC50_exp_transformed.rds"
  )
)


cat(
  "\n========================================\n",
  "SECTION 2 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 3
# Prepare TCGA-LUAD expression for prediction
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 3: PREPARE TCGA-LUAD EXPRESSION\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 3.1 Load TCGA-LUAD protein-coding TPM
# ------------------------------------------------------------


luad_patient_tpm_s32 <-
  readRDS(
    patient_tpm_file_s32
  )


stopifnot(
  is.matrix(
    luad_patient_tpm_s32
  ) ||
    is.data.frame(
      luad_patient_tpm_s32
    )
)


luad_patient_tpm_s32 <-
  as.matrix(
    luad_patient_tpm_s32
  )


stopifnot(
  !anyNA(
    luad_patient_tpm_s32
  )
)


stopifnot(
  all(
    is.finite(
      luad_patient_tpm_s32
    )
  )
)


stopifnot(
  all(
    luad_patient_tpm_s32 >=
      0
  )
)


# ------------------------------------------------------------
# 3.2 log2(TPM + 1)
# ------------------------------------------------------------


luad_log2_tpm_s32 <-
  log2(
    luad_patient_tpm_s32 +
      1
  )


cat(
  "\nTCGA-LUAD log2(TPM+1) dimensions:\n"
)


print(
  dim(
    luad_log2_tpm_s32
  )
)


# ------------------------------------------------------------
# 3.3 Identify overlapping genes
# ------------------------------------------------------------


common_genes_s32 <-
  intersect(
    rownames(
      gdsc2_expr_s32
    ),
    rownames(
      luad_log2_tpm_s32
    )
  )


cat(
  "\nGDSC2 genes:\n"
)


print(
  nrow(
    gdsc2_expr_s32
  )
)


cat(
  "\nTCGA-LUAD genes:\n"
)


print(
  nrow(
    luad_log2_tpm_s32
  )
)


cat(
  "\nCommon genes:\n"
)


print(
  length(
    common_genes_s32
  )
)


cat(
  "\nPercentage of GDSC2 genes available in TCGA-LUAD:\n"
)


print(
  100 *
    length(
      common_genes_s32
    ) /
    nrow(
      gdsc2_expr_s32
    )
)


stopifnot(
  length(
    common_genes_s32
  ) ==
    15979
)


# ------------------------------------------------------------
# 3.4 Restrict both matrices to identical genes/order
# ------------------------------------------------------------


training_expr_s32 <-
  gdsc2_expr_s32[
    common_genes_s32,
    ,
    drop = FALSE
  ]


test_expr_s32 <-
  luad_log2_tpm_s32[
    common_genes_s32,
    ,
    drop = FALSE
  ]


stopifnot(
  identical(
    rownames(
      training_expr_s32
    ),
    rownames(
      test_expr_s32
    )
  )
)


stopifnot(
  identical(
    colnames(
      training_expr_s32
    ),
    rownames(
      gdsc2_res_s32
    )
  )
)


cat(
  "\nTraining expression dimensions:\n"
)


print(
  dim(
    training_expr_s32
  )
)


cat(
  "\nTCGA-LUAD test expression dimensions:\n"
)


print(
  dim(
    test_expr_s32
  )
)


cat(
  "\n========================================\n",
  "SECTION 3 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 4
# Predict TCGA-LUAD GDSC2 drug sensitivity
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 4: RUN calcPhenotype()\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 4.1 Preserve project working directory
# ------------------------------------------------------------


original_wd_s32 <-
  getwd()


prediction_dir_absolute_s32 <-
  normalizePath(
    prediction_dir_s32,
    winslash = "/",
    mustWork = TRUE
  )


# ------------------------------------------------------------
# 4.2 Run oncoPredict
#
# trainingExprData:
#   GDSC2 RMA normalized/log-transformed expression
#
# trainingPtype:
#   exponentiated GDSC2 response matrix
#
# testExprData:
#   TCGA-LUAD log2(TPM+1)
#
# batchCorrect = "eb":
#   empirical-Bayes/ComBat correction
# ------------------------------------------------------------


gdsc2_predictions_s32 <-
  tryCatch(
    
    {
      
      setwd(
        prediction_dir_absolute_s32
      )
      
      
      oncoPredict::calcPhenotype(
        
        trainingExprData =
          training_expr_s32,
        
        trainingPtype =
          gdsc2_res_s32,
        
        testExprData =
          test_expr_s32,
        
        batchCorrect =
          "eb",
        
        powerTransformPhenotype =
          TRUE,
        
        removeLowVaryingGenes =
          0.2,
        
        minNumSamples =
          10,
        
        selection =
          1,
        
        printOutput =
          TRUE,
        
        pcr =
          FALSE,
        
        removeLowVaringGenesFrom =
          "rawData",
        
        report_pc =
          FALSE,
        
        cc =
          FALSE,
        
        percent =
          80,
        
        rsq =
          FALSE,
        
        folder =
          FALSE,
        
        parallel =
          FALSE,
        
        cores =
          1
      )
    },
    
    finally = {
      
      setwd(
        original_wd_s32
      )
    }
  )


# ------------------------------------------------------------
# 4.3 Prediction QC
# ------------------------------------------------------------


cat(
  "\nPrediction object class:\n"
)


print(
  class(
    gdsc2_predictions_s32
  )
)


cat(
  "\nPrediction dimensions:\n"
)


print(
  dim(
    gdsc2_predictions_s32
  )
)


cat(
  "\nPrediction NA count:\n"
)


print(
  sum(
    is.na(
      gdsc2_predictions_s32
    )
  )
)


cat(
  "\nPrediction non-finite count:\n"
)


print(
  sum(
    !is.finite(
      gdsc2_predictions_s32
    )
  )
)


cat(
  "\nPrediction summary:\n"
)


print(
  summary(
    as.vector(
      gdsc2_predictions_s32
    )
  )
)


stopifnot(
  is.matrix(
    gdsc2_predictions_s32
  )
)


stopifnot(
  identical(
    dim(
      gdsc2_predictions_s32
    ),
    c(
      517L,
      198L
    )
  )
)


stopifnot(
  !anyNA(
    gdsc2_predictions_s32
  )
)


stopifnot(
  all(
    is.finite(
      gdsc2_predictions_s32
    )
  )
)


missing_prediction_patients_s32 <-
  setdiff(
    colnames(
      test_expr_s32
    ),
    rownames(
      gdsc2_predictions_s32
    )
  )


unexpected_prediction_patients_s32 <-
  setdiff(
    rownames(
      gdsc2_predictions_s32
    ),
    colnames(
      test_expr_s32
    )
  )


cat(
  "\nMissing prediction patients:\n"
)


print(
  length(
    missing_prediction_patients_s32
  )
)


cat(
  "\nUnexpected prediction patients:\n"
)


print(
  length(
    unexpected_prediction_patients_s32
  )
)


stopifnot(
  length(
    missing_prediction_patients_s32
  ) ==
    0
)


stopifnot(
  length(
    unexpected_prediction_patients_s32
  ) ==
    0
)


# ------------------------------------------------------------
# 4.4 Save predictions
# ------------------------------------------------------------


prediction_rds_s32 <-
  file.path(
    gdsc_processed_dir_s32,
    "TCGA_LUAD_GDSC2_predicted_IC50.rds"
  )


prediction_csv_s32 <-
  file.path(
    gdsc_processed_dir_s32,
    "TCGA_LUAD_GDSC2_predicted_IC50.csv"
  )


saveRDS(
  gdsc2_predictions_s32,
  prediction_rds_s32
)


write.csv(
  gdsc2_predictions_s32,
  prediction_csv_s32,
  row.names = TRUE
)


stopifnot(
  file.exists(
    prediction_rds_s32
  )
)


stopifnot(
  file.exists(
    prediction_csv_s32
  )
)


cat(
  "\n========================================\n",
  "SECTION 4 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 5
# MPCDS and PCD-subtype drug associations
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 5: MPCDS AND SUBTYPE ASSOCIATIONS\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 5.1 Load MPCDS and subtype assignments
# ------------------------------------------------------------


mpcds_s32 <-
  read.csv(
    mpcds_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


subtype_s32 <-
  read.csv(
    subtype_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


stopifnot(
  all(
    c(
      "patient_id",
      "MPCDS"
    ) %in%
      colnames(
        mpcds_s32
      )
  )
)


stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster"
    ) %in%
      colnames(
        subtype_s32
      )
  )
)


# ------------------------------------------------------------
# 5.2 Parse GDSC2 drug names and IDs
# ------------------------------------------------------------


drug_full_names_s32 <-
  colnames(
    gdsc2_predictions_s32
  )


drug_ids_s32 <-
  sub(
    "^.*_([0-9]+)$",
    "\\1",
    drug_full_names_s32
  )


drug_names_s32 <-
  sub(
    "_[0-9]+$",
    "",
    drug_full_names_s32
  )


drug_metadata_s32 <-
  data.frame(
    
    drug_column =
      drug_full_names_s32,
    
    drug_name =
      drug_names_s32,
    
    gdsc_drug_id =
      drug_ids_s32,
    
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# 5.3 MPCDS analysis cohort
# ------------------------------------------------------------


mpcds_match_s32 <-
  match(
    rownames(
      gdsc2_predictions_s32
    ),
    mpcds_s32$patient_id
  )


mpcds_available_s32 <-
  !is.na(
    mpcds_match_s32
  )


cat(
  "\nPatients with both predictions and MPCDS:\n"
)


print(
  sum(
    mpcds_available_s32
  )
)


stopifnot(
  sum(
    mpcds_available_s32
  ) ==
    504
)


mpcds_values_s32 <-
  mpcds_s32$MPCDS[
    mpcds_match_s32[
      mpcds_available_s32
    ]
  ]


mpcds_prediction_matrix_s32 <-
  gdsc2_predictions_s32[
    mpcds_available_s32,
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 5.4 Spearman MPCDS-drug associations
# ------------------------------------------------------------


mpcds_drug_results_s32 <-
  do.call(
    rbind,
    lapply(
      seq_len(
        ncol(
          mpcds_prediction_matrix_s32
        )
      ),
      function(j_s32) {
        
        drug_values_s32 <-
          mpcds_prediction_matrix_s32[
            ,
            j_s32
          ]
        
        
        complete_s32 <-
          complete.cases(
            mpcds_values_s32,
            drug_values_s32
          )
        
        
        cor_test_s32 <-
          cor.test(
            mpcds_values_s32[
              complete_s32
            ],
            drug_values_s32[
              complete_s32
            ],
            method = "spearman",
            exact = FALSE
          )
        
        
        rho_s32 <-
          unname(
            cor_test_s32$estimate
          )
        
        
        direction_s32 <-
          if (
            rho_s32 <
            0
          ) {
            
            "Higher MPCDS associated with lower predicted IC50 (greater predicted sensitivity)"
            
          } else if (
            rho_s32 >
            0
          ) {
            
            "Higher MPCDS associated with higher predicted IC50 (greater predicted resistance)"
            
          } else {
            
            "No directional association"
          }
        
        
        data.frame(
          
          drug_column =
            colnames(
              mpcds_prediction_matrix_s32
            )[
              j_s32
            ],
          
          n =
            sum(
              complete_s32
            ),
          
          rho =
            rho_s32,
          
          p_value =
            cor_test_s32$p.value,
          
          direction =
            direction_s32,
          
          stringsAsFactors = FALSE
        )
      }
    )
  )


mpcds_drug_results_s32$FDR <-
  p.adjust(
    mpcds_drug_results_s32$p_value,
    method = "BH"
  )


mpcds_drug_results_s32 <-
  merge(
    drug_metadata_s32,
    mpcds_drug_results_s32,
    by = "drug_column",
    all.y = TRUE,
    sort = FALSE
  )


mpcds_drug_results_s32 <-
  mpcds_drug_results_s32[
    match(
      drug_full_names_s32,
      mpcds_drug_results_s32$drug_column
    ),
    ,
    drop = FALSE
  ]


mpcds_drug_results_s32$significant_FDR_0.05 <-
  mpcds_drug_results_s32$FDR <
  0.05


mpcds_drug_results_s32$nominal_p_0.05 <-
  mpcds_drug_results_s32$p_value <
  0.05


cat(
  "\n========================================\n",
  "MPCDS / GDSC2 ASSOCIATION SUMMARY\n",
  "========================================\n"
)


cat(
  "\nDrugs tested:\n"
)


print(
  nrow(
    mpcds_drug_results_s32
  )
)


cat(
  "\nFDR-significant drugs:\n"
)


print(
  sum(
    mpcds_drug_results_s32$significant_FDR_0.05
  )
)


cat(
  "\nNominal p < 0.05 drugs:\n"
)


print(
  sum(
    mpcds_drug_results_s32$nominal_p_0.05
  )
)


cat(
  "\nTop 20 MPCDS associations:\n"
)


print(
  head(
    mpcds_drug_results_s32[
      order(
        mpcds_drug_results_s32$FDR,
        mpcds_drug_results_s32$p_value
      ),
      ,
      drop = FALSE
    ],
    20
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5.5 PCD-subtype analysis cohort
# ------------------------------------------------------------


subtype_match_s32 <-
  match(
    rownames(
      gdsc2_predictions_s32
    ),
    subtype_s32$patient_id
  )


subtype_available_s32 <-
  !is.na(
    subtype_match_s32
  )


stopifnot(
  sum(
    subtype_available_s32
  ) ==
    517
)


subtype_values_s32 <-
  subtype_s32$PCD_cluster[
    subtype_match_s32[
      subtype_available_s32
    ]
  ]


subtype_prediction_matrix_s32 <-
  gdsc2_predictions_s32[
    subtype_available_s32,
    ,
    drop = FALSE
  ]


cat(
  "\nSubtype distribution:\n"
)


print(
  table(
    subtype_values_s32
  )
)


stopifnot(
  all(
    c(
      "PCD_C1",
      "PCD_C2"
    ) %in%
      subtype_values_s32
  )
)


# ------------------------------------------------------------
# 5.6 Wilcoxon subtype comparisons
# ------------------------------------------------------------


subtype_drug_results_s32 <-
  do.call(
    rbind,
    lapply(
      seq_len(
        ncol(
          subtype_prediction_matrix_s32
        )
      ),
      function(j_s32) {
        
        drug_values_s32 <-
          subtype_prediction_matrix_s32[
            ,
            j_s32
          ]
        
        
        c1_s32 <-
          drug_values_s32[
            subtype_values_s32 ==
              "PCD_C1"
          ]
        
        
        c2_s32 <-
          drug_values_s32[
            subtype_values_s32 ==
              "PCD_C2"
          ]
        
        
        c1_s32 <-
          c1_s32[
            is.finite(
              c1_s32
            )
          ]
        
        
        c2_s32 <-
          c2_s32[
            is.finite(
              c2_s32
            )
          ]
        
        
        wilcox_s32 <-
          wilcox.test(
            c1_s32,
            c2_s32,
            exact = FALSE
          )
        
        
        median_c1_s32 <-
          median(
            c1_s32
          )
        
        
        median_c2_s32 <-
          median(
            c2_s32
          )
        
        
        median_difference_s32 <-
          median_c2_s32 -
          median_c1_s32
        
        
        direction_s32 <-
          if (
            median_difference_s32 <
            0
          ) {
            
            "PCD_C2 lower predicted IC50 than PCD_C1 (greater predicted sensitivity in PCD_C2)"
            
          } else if (
            median_difference_s32 >
            0
          ) {
            
            "PCD_C2 higher predicted IC50 than PCD_C1 (greater predicted resistance in PCD_C2)"
            
          } else {
            
            "No median difference"
          }
        
        
        data.frame(
          
          drug_column =
            colnames(
              subtype_prediction_matrix_s32
            )[
              j_s32
            ],
          
          n_C1 =
            length(
              c1_s32
            ),
          
          n_C2 =
            length(
              c2_s32
            ),
          
          median_C1 =
            median_c1_s32,
          
          median_C2 =
            median_c2_s32,
          
          median_difference_C2_minus_C1 =
            median_difference_s32,
          
          W =
            unname(
              wilcox_s32$statistic
            ),
          
          p_value =
            wilcox_s32$p.value,
          
          direction =
            direction_s32,
          
          stringsAsFactors = FALSE
        )
      }
    )
  )


subtype_drug_results_s32$FDR <-
  p.adjust(
    subtype_drug_results_s32$p_value,
    method = "BH"
  )


subtype_drug_results_s32 <-
  merge(
    drug_metadata_s32,
    subtype_drug_results_s32,
    by = "drug_column",
    all.y = TRUE,
    sort = FALSE
  )


subtype_drug_results_s32 <-
  subtype_drug_results_s32[
    match(
      drug_full_names_s32,
      subtype_drug_results_s32$drug_column
    ),
    ,
    drop = FALSE
  ]


subtype_drug_results_s32$significant_FDR_0.05 <-
  subtype_drug_results_s32$FDR <
  0.05


subtype_drug_results_s32$nominal_p_0.05 <-
  subtype_drug_results_s32$p_value <
  0.05


cat(
  "\n========================================\n",
  "PCD SUBTYPE / GDSC2 ASSOCIATION SUMMARY\n",
  "========================================\n"
)


cat(
  "\nDrugs tested:\n"
)


print(
  nrow(
    subtype_drug_results_s32
  )
)


cat(
  "\nFDR-significant drugs:\n"
)


print(
  sum(
    subtype_drug_results_s32$significant_FDR_0.05
  )
)


cat(
  "\nNominal p < 0.05 drugs:\n"
)


print(
  sum(
    subtype_drug_results_s32$nominal_p_0.05
  )
)


# ------------------------------------------------------------
# 5.7 Direction summaries
# ------------------------------------------------------------


mpcds_direction_summary_s32 <-
  data.frame(
    
    category =
      c(
        "FDR significant: higher MPCDS -> sensitivity",
        "FDR significant: higher MPCDS -> resistance",
        "Nominal p<0.05: higher MPCDS -> sensitivity",
        "Nominal p<0.05: higher MPCDS -> resistance"
      ),
    
    n =
      c(
        
        sum(
          mpcds_drug_results_s32$FDR <
            0.05 &
            mpcds_drug_results_s32$rho <
            0
        ),
        
        sum(
          mpcds_drug_results_s32$FDR <
            0.05 &
            mpcds_drug_results_s32$rho >
            0
        ),
        
        sum(
          mpcds_drug_results_s32$p_value <
            0.05 &
            mpcds_drug_results_s32$rho <
            0
        ),
        
        sum(
          mpcds_drug_results_s32$p_value <
            0.05 &
            mpcds_drug_results_s32$rho >
            0
        )
      ),
    
    stringsAsFactors = FALSE
  )


subtype_direction_summary_s32 <-
  data.frame(
    
    category =
      c(
        "FDR significant: PCD_C2 more sensitive",
        "FDR significant: PCD_C2 more resistant",
        "Nominal p<0.05: PCD_C2 more sensitive",
        "Nominal p<0.05: PCD_C2 more resistant"
      ),
    
    n =
      c(
        
        sum(
          subtype_drug_results_s32$FDR <
            0.05 &
            subtype_drug_results_s32[[
              "median_difference_C2_minus_C1"
            ]] <
            0
        ),
        
        sum(
          subtype_drug_results_s32$FDR <
            0.05 &
            subtype_drug_results_s32[[
              "median_difference_C2_minus_C1"
            ]] >
            0
        ),
        
        sum(
          subtype_drug_results_s32$p_value <
            0.05 &
            subtype_drug_results_s32[[
              "median_difference_C2_minus_C1"
            ]] <
            0
        ),
        
        sum(
          subtype_drug_results_s32$p_value <
            0.05 &
            subtype_drug_results_s32[[
              "median_difference_C2_minus_C1"
            ]] >
            0
        )
      ),
    
    stringsAsFactors = FALSE
  )


cat(
  "\nMPCDS direction summary:\n"
)


print(
  mpcds_direction_summary_s32,
  row.names = FALSE
)


cat(
  "\nSubtype direction summary:\n"
)


print(
  subtype_direction_summary_s32,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5.8 Save full association results
# ------------------------------------------------------------


write.csv(
  mpcds_drug_results_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_MPCDS_all_drug_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  subtype_drug_results_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PCD_subtype_all_drug_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  mpcds_drug_results_s32[
    mpcds_drug_results_s32$FDR <
      0.05,
    ,
    drop = FALSE
  ],
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_MPCDS_FDR_significant.csv"
  ),
  row.names = FALSE
)


write.csv(
  subtype_drug_results_s32[
    subtype_drug_results_s32$FDR <
      0.05,
    ,
    drop = FALSE
  ],
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PCD_subtype_FDR_significant.csv"
  ),
  row.names = FALSE
)


write.csv(
  mpcds_drug_results_s32[
    mpcds_drug_results_s32$p_value <
      0.05,
    ,
    drop = FALSE
  ],
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_MPCDS_nominal_p_lt_0.05.csv"
  ),
  row.names = FALSE
)


write.csv(
  subtype_drug_results_s32[
    subtype_drug_results_s32$p_value <
      0.05,
    ,
    drop = FALSE
  ],
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PCD_subtype_nominal_p_lt_0.05.csv"
  ),
  row.names = FALSE
)


write.csv(
  mpcds_direction_summary_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_MPCDS_direction_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  subtype_direction_summary_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PCD_subtype_direction_summary.csv"
  ),
  row.names = FALSE
)


cat(
  "\n========================================\n",
  "SECTION 5 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 6
# GDSC2 vs DepMap/PRISM cross-platform comparison
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 6: GDSC2 / PRISM CROSS-PLATFORM ANALYSIS\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 6.1 Load PRISM results
# ------------------------------------------------------------


prism_s32 <-
  read.csv(
    prism_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


required_prism_columns_s32 <-
  c(
    "drug_full_name",
    "drug_name",
    "BRD_ID",
    "n",
    "spearman_rho",
    "p_value",
    "FDR"
  )


stopifnot(
  all(
    required_prism_columns_s32 %in%
      colnames(
        prism_s32
      )
  )
)


# ------------------------------------------------------------
# 6.2 Conservative drug-name normalization
# ------------------------------------------------------------


normalize_drug_name_s32 <-
  function(x) {
    
    x <-
      toupper(
        trimws(
          as.character(
            x
          )
        )
      )
    
    
    x <-
      gsub(
        "^\\(\\+\\)-",
        "",
        x
      )
    
    
    x <-
      gsub(
        "^\\(-\\)-",
        "",
        x
      )
    
    
    x <-
      gsub(
        "[^A-Z0-9]",
        "",
        x
      )
    
    
    x
  }


gdsc_s32 <-
  mpcds_drug_results_s32


gdsc_s32$drug_key <-
  normalize_drug_name_s32(
    gdsc_s32$drug_name
  )


prism_s32$drug_key <-
  normalize_drug_name_s32(
    prism_s32$drug_name
  )


# ------------------------------------------------------------
# 6.3 Identify ambiguous normalized names
# ------------------------------------------------------------


gdsc_duplicate_keys_s32 <-
  unique(
    gdsc_s32$drug_key[
      duplicated(
        gdsc_s32$drug_key
      ) |
        duplicated(
          gdsc_s32$drug_key,
          fromLast = TRUE
        )
    ]
  )


prism_duplicate_keys_s32 <-
  unique(
    prism_s32$drug_key[
      duplicated(
        prism_s32$drug_key
      ) |
        duplicated(
          prism_s32$drug_key,
          fromLast = TRUE
        )
    ]
  )


cat(
  "\nDuplicated GDSC2 normalized drug keys:\n"
)


print(
  length(
    gdsc_duplicate_keys_s32
  )
)


cat(
  "\nDuplicated PRISM normalized drug keys:\n"
)


print(
  length(
    prism_duplicate_keys_s32
  )
)


# ------------------------------------------------------------
# 6.4 Strict unambiguous matching
# ------------------------------------------------------------


gdsc_unique_s32 <-
  gdsc_s32[
    !gdsc_s32$drug_key %in%
      gdsc_duplicate_keys_s32,
    ,
    drop = FALSE
  ]


prism_unique_s32 <-
  prism_s32[
    !prism_s32$drug_key %in%
      prism_duplicate_keys_s32,
    ,
    drop = FALSE
  ]


common_drug_keys_s32 <-
  intersect(
    gdsc_unique_s32$drug_key,
    prism_unique_s32$drug_key
  )


cat(
  "\n========================================\n",
  "STRICT CROSS-PLATFORM DRUG OVERLAP\n",
  "========================================\n"
)


cat(
  "\nStrict normalized-name overlap:\n"
)


print(
  length(
    common_drug_keys_s32
  )
)


cat(
  "\nPercentage of GDSC2 entries represented in PRISM:\n"
)


print(
  100 *
    length(
      common_drug_keys_s32
    ) /
    nrow(
      gdsc_s32
    )
)


stopifnot(
  length(
    common_drug_keys_s32
  ) ==
    73
)


# ------------------------------------------------------------
# 6.5 Build matched table
# ------------------------------------------------------------


gdsc_match_s32 <-
  gdsc_unique_s32[
    match(
      common_drug_keys_s32,
      gdsc_unique_s32$drug_key
    ),
    ,
    drop = FALSE
  ]


prism_match_s32 <-
  prism_unique_s32[
    match(
      common_drug_keys_s32,
      prism_unique_s32$drug_key
    ),
    ,
    drop = FALSE
  ]


stopifnot(
  identical(
    gdsc_match_s32$drug_key,
    prism_match_s32$drug_key
  )
)


cross_platform_all_s32 <-
  data.frame(
    
    drug_key =
      common_drug_keys_s32,
    
    GDSC_drug_name =
      gdsc_match_s32$drug_name,
    
    GDSC_drug_id =
      gdsc_match_s32$gdsc_drug_id,
    
    GDSC_rho =
      gdsc_match_s32$rho,
    
    GDSC_p =
      gdsc_match_s32$p_value,
    
    GDSC_FDR =
      gdsc_match_s32$FDR,
    
    PRISM_drug_name =
      prism_match_s32$drug_name,
    
    PRISM_BRD_ID =
      prism_match_s32$BRD_ID,
    
    PRISM_n =
      prism_match_s32$n,
    
    PRISM_rho =
      prism_match_s32$spearman_rho,
    
    PRISM_p =
      prism_match_s32$p_value,
    
    PRISM_FDR =
      prism_match_s32$FDR,
    
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# 6.6 Testability filtering
# ------------------------------------------------------------


cross_platform_all_s32$testable_both_platforms <-
  is.finite(
    cross_platform_all_s32$GDSC_rho
  ) &
  is.finite(
    cross_platform_all_s32$PRISM_rho
  )


non_testable_matched_s32 <-
  cross_platform_all_s32[
    !cross_platform_all_s32$testable_both_platforms,
    ,
    drop = FALSE
  ]


cross_platform_s32 <-
  cross_platform_all_s32[
    cross_platform_all_s32$testable_both_platforms,
    ,
    drop = FALSE
  ]


cat(
  "\n========================================\n",
  "CROSS-PLATFORM TESTABILITY CHECK\n",
  "========================================\n"
)


cat(
  "\nMatched before filtering:\n"
)


print(
  nrow(
    cross_platform_all_s32
  )
)


cat(
  "\nTestable on both platforms:\n"
)


print(
  nrow(
    cross_platform_s32
  )
)


cat(
  "\nExcluded as non-testable:\n"
)


print(
  nrow(
    non_testable_matched_s32
  )
)


if (
  nrow(
    non_testable_matched_s32
  ) >
  0
) {
  
  print(
    non_testable_matched_s32[
      ,
      c(
        "GDSC_drug_name",
        "GDSC_rho",
        "PRISM_drug_name",
        "PRISM_rho",
        "PRISM_n"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}


stopifnot(
  nrow(
    cross_platform_s32
  ) ==
    72
)


# ------------------------------------------------------------
# 6.7 Biological direction
#
# Lower GDSC2 IC50 = greater sensitivity
# Lower PRISM AUC   = greater sensitivity
#
# Therefore:
# rho < 0 = higher MPCDS associated with greater sensitivity
# on both platforms.
# ------------------------------------------------------------


cross_platform_s32$GDSC_direction <-
  ifelse(
    cross_platform_s32$GDSC_rho <
      0,
    "Sensitivity",
    ifelse(
      cross_platform_s32$GDSC_rho >
        0,
      "Resistance",
      "Neutral"
    )
  )


cross_platform_s32$PRISM_direction <-
  ifelse(
    cross_platform_s32$PRISM_rho <
      0,
    "Sensitivity",
    ifelse(
      cross_platform_s32$PRISM_rho >
        0,
      "Resistance",
      "Neutral"
    )
  )


cross_platform_s32$direction_concordant <-
  sign(
    cross_platform_s32$GDSC_rho
  ) ==
  sign(
    cross_platform_s32$PRISM_rho
  )


stopifnot(
  !anyNA(
    cross_platform_s32$direction_concordant
  )
)


# ------------------------------------------------------------
# 6.8 Significance classifications
# ------------------------------------------------------------


cross_platform_s32$GDSC_FDR_significant <-
  cross_platform_s32$GDSC_FDR <
  0.05


cross_platform_s32$GDSC_nominal <-
  cross_platform_s32$GDSC_p <
  0.05


cross_platform_s32$PRISM_FDR_significant <-
  cross_platform_s32$PRISM_FDR <
  0.05


cross_platform_s32$PRISM_nominal <-
  cross_platform_s32$PRISM_p <
  0.05


cross_platform_s32$cross_platform_support <-
  ifelse(
    
    cross_platform_s32$GDSC_FDR_significant &
      cross_platform_s32$PRISM_nominal &
      cross_platform_s32$direction_concordant,
    
    "GDSC2 FDR-significant + PRISM nominal directionally concordant",
    
    ifelse(
      
      cross_platform_s32$GDSC_FDR_significant &
        cross_platform_s32$PRISM_nominal &
        !cross_platform_s32$direction_concordant,
      
      "GDSC2 FDR-significant + PRISM nominal but directionally discordant",
      
      ifelse(
        
        cross_platform_s32$GDSC_FDR_significant &
          cross_platform_s32$direction_concordant,
        
        "GDSC2 FDR-significant; PRISM direction concordant but not nominal",
        
        ifelse(
          
          cross_platform_s32$GDSC_FDR_significant,
          
          "GDSC2 FDR-significant only",
          
          ifelse(
            
            cross_platform_s32$PRISM_nominal,
            
            "PRISM nominal only",
            
            "Neither cross-platform significance criterion"
          )
        )
      )
    )
  )


# ------------------------------------------------------------
# 6.9 Directional concordance
# ------------------------------------------------------------


n_cross_platform_s32 <-
  nrow(
    cross_platform_s32
  )


n_cross_concordant_s32 <-
  sum(
    cross_platform_s32$direction_concordant
  )


n_cross_discordant_s32 <-
  sum(
    !cross_platform_s32$direction_concordant
  )


cross_concordance_s32 <-
  n_cross_concordant_s32 /
  n_cross_platform_s32


cross_binom_s32 <-
  binom.test(
    
    x =
      n_cross_concordant_s32,
    
    n =
      n_cross_platform_s32,
    
    p =
      0.5,
    
    alternative =
      "greater"
  )


cat(
  "\n========================================\n",
  "CROSS-PLATFORM DIRECTIONAL CONCORDANCE\n",
  "========================================\n"
)


cat(
  "\nTestable drugs:\n"
)


print(
  n_cross_platform_s32
)


cat(
  "\nConcordant:\n"
)


print(
  n_cross_concordant_s32
)


cat(
  "\nDiscordant:\n"
)


print(
  n_cross_discordant_s32
)


cat(
  "\nConcordance proportion:\n"
)


print(
  cross_concordance_s32
)


cat(
  "\nBinomial p-value:\n"
)


print(
  cross_binom_s32$p.value
)


# ------------------------------------------------------------
# 6.10 Cross-platform effect-size correlation
# ------------------------------------------------------------


cross_effect_test_s32 <-
  cor.test(
    cross_platform_s32$GDSC_rho,
    cross_platform_s32$PRISM_rho,
    method = "spearman",
    exact = FALSE
  )


cat(
  "\n========================================\n",
  "CROSS-PLATFORM EFFECT-SIZE CORRELATION\n",
  "========================================\n"
)


print(
  cross_effect_test_s32
)


# ------------------------------------------------------------
# 6.11 GDSC2-FDR matched drugs
# ------------------------------------------------------------


gdsc_fdr_matched_s32 <-
  cross_platform_s32[
    cross_platform_s32$GDSC_FDR_significant,
    ,
    drop = FALSE
  ]


supported_s32 <-
  cross_platform_s32[
    cross_platform_s32$GDSC_FDR_significant &
      cross_platform_s32$PRISM_nominal &
      cross_platform_s32$direction_concordant,
    ,
    drop = FALSE
  ]


discordant_supported_s32 <-
  cross_platform_s32[
    cross_platform_s32$GDSC_FDR_significant &
      cross_platform_s32$PRISM_nominal &
      !cross_platform_s32$direction_concordant,
    ,
    drop = FALSE
  ]


cat(
  "\n========================================\n",
  "GDSC2 FDR-SIGNIFICANT MATCHED DRUGS\n",
  "========================================\n"
)


cat(
  "\nGDSC2 FDR-significant matched drugs:\n"
)


print(
  nrow(
    gdsc_fdr_matched_s32
  )
)


cat(
  "\nDirectionally concordant with PRISM:\n"
)


print(
  sum(
    gdsc_fdr_matched_s32$direction_concordant
  )
)


cat(
  "\nDirectionally discordant with PRISM:\n"
)


print(
  sum(
    !gdsc_fdr_matched_s32$direction_concordant
  )
)


cat(
  "\nPRISM nominal among GDSC2 FDR hits:\n"
)


print(
  sum(
    gdsc_fdr_matched_s32$PRISM_nominal
  )
)


cat(
  "\nDirectionally supported GDSC2-FDR + PRISM-nominal drugs:\n"
)


print(
  nrow(
    supported_s32
  )
)


# ------------------------------------------------------------
# 6.12 Cross-platform support summary
# ------------------------------------------------------------


support_summary_s32 <-
  as.data.frame(
    table(
      cross_platform_s32$cross_platform_support
    ),
    stringsAsFactors = FALSE
  )


colnames(
  support_summary_s32
) <-
  c(
    "classification",
    "n_drugs"
  )


cat(
  "\n========================================\n",
  "CROSS-PLATFORM SUPPORT SUMMARY\n",
  "========================================\n"
)


print(
  support_summary_s32,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6.13 Save cross-platform outputs
# ------------------------------------------------------------


write.csv(
  cross_platform_all_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PRISM_MPCDS_all_strict_name_matches.csv"
  ),
  row.names = FALSE
)


write.csv(
  cross_platform_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PRISM_MPCDS_cross_platform_testable_drugs.csv"
  ),
  row.names = FALSE
)


write.csv(
  non_testable_matched_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PRISM_MPCDS_non_testable_matched_drugs.csv"
  ),
  row.names = FALSE
)


write.csv(
  supported_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PRISM_MPCDS_directionally_supported_drugs.csv"
  ),
  row.names = FALSE
)


write.csv(
  discordant_supported_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PRISM_MPCDS_discordant_drugs.csv"
  ),
  row.names = FALSE
)


write.csv(
  support_summary_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_PRISM_MPCDS_cross_platform_summary.csv"
  ),
  row.names = FALSE
)


gdsc_ambiguous_s32 <-
  gdsc_s32[
    gdsc_s32$drug_key %in%
      gdsc_duplicate_keys_s32,
    ,
    drop = FALSE
  ]


prism_ambiguous_s32 <-
  prism_s32[
    prism_s32$drug_key %in%
      prism_duplicate_keys_s32,
    ,
    drop = FALSE
  ]


write.csv(
  gdsc_ambiguous_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_ambiguous_normalized_drug_names.csv"
  ),
  row.names = FALSE
)


write.csv(
  prism_ambiguous_s32,
  file.path(
    gdsc_results_dir_s32,
    "DepMap_PRISM_ambiguous_normalized_drug_names.csv"
  ),
  row.names = FALSE
)


cat(
  "\n========================================\n",
  "SECTION 6 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 7
# Final summaries, Figure S15, interpretation, and audit
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 7: FINAL SUMMARY AND FIGURE\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 7.1 Reconstruct final statistics
# ------------------------------------------------------------


n_gdsc_drugs_s32 <-
  nrow(
    mpcds_drug_results_s32
  )


n_mpcds_fdr_s32 <-
  sum(
    mpcds_drug_results_s32$FDR <
      0.05
  )


n_mpcds_sensitivity_s32 <-
  sum(
    mpcds_drug_results_s32$FDR <
      0.05 &
      mpcds_drug_results_s32$rho <
      0
  )


n_mpcds_resistance_s32 <-
  sum(
    mpcds_drug_results_s32$FDR <
      0.05 &
      mpcds_drug_results_s32$rho >
      0
  )


n_subtype_fdr_s32 <-
  sum(
    subtype_drug_results_s32$FDR <
      0.05
  )


n_c2_sensitive_s32 <-
  sum(
    subtype_drug_results_s32$FDR <
      0.05 &
      subtype_drug_results_s32[[
        "median_difference_C2_minus_C1"
      ]] <
      0
  )


n_c2_resistant_s32 <-
  sum(
    subtype_drug_results_s32$FDR <
      0.05 &
      subtype_drug_results_s32[[
        "median_difference_C2_minus_C1"
      ]] >
      0
  )


n_gdsc_fdr_matched_s32 <-
  sum(
    cross_platform_s32$GDSC_FDR_significant
  )


n_gdsc_fdr_prism_nominal_s32 <-
  sum(
    cross_platform_s32$GDSC_FDR_significant &
      cross_platform_s32$PRISM_nominal
  )


# ------------------------------------------------------------
# 7.2 Final summary output
# ------------------------------------------------------------


cat(
  "\n========================================\n",
  "FINAL GDSC2 SUMMARY\n",
  "========================================\n"
)


cat(
  "\nTotal GDSC2 drugs tested:\n"
)


print(
  n_gdsc_drugs_s32
)


cat(
  "\nMPCDS FDR-significant drugs:\n"
)


print(
  n_mpcds_fdr_s32
)


cat(
  "\nHigher MPCDS -> greater predicted sensitivity:\n"
)


print(
  n_mpcds_sensitivity_s32
)


cat(
  "\nHigher MPCDS -> greater predicted resistance:\n"
)


print(
  n_mpcds_resistance_s32
)


cat(
  "\nSubtype FDR-significant drugs:\n"
)


print(
  n_subtype_fdr_s32
)


cat(
  "\nPCD_C2 more sensitive:\n"
)


print(
  n_c2_sensitive_s32
)


cat(
  "\nPCD_C2 more resistant:\n"
)


print(
  n_c2_resistant_s32
)


cat(
  "\nCross-platform testable matched drugs:\n"
)


print(
  n_cross_platform_s32
)


cat(
  "\nCross-platform concordant drugs:\n"
)


print(
  n_cross_concordant_s32
)


cat(
  "\nCross-platform discordant drugs:\n"
)


print(
  n_cross_discordant_s32
)


cat(
  "\nCross-platform concordance proportion:\n"
)


print(
  cross_concordance_s32
)


cat(
  "\nCross-platform concordance binomial p-value:\n"
)


print(
  cross_binom_s32$p.value
)


cat(
  "\nCross-platform effect-size Spearman rho:\n"
)


print(
  unname(
    cross_effect_test_s32$estimate
  )
)


cat(
  "\nCross-platform effect-size p-value:\n"
)


print(
  cross_effect_test_s32$p.value
)


cat(
  "\nGDSC2-FDR matched drugs:\n"
)


print(
  n_gdsc_fdr_matched_s32
)


cat(
  "\nGDSC2-FDR drugs also PRISM nominal:\n"
)


print(
  n_gdsc_fdr_prism_nominal_s32
)


# ------------------------------------------------------------
# 7.3 Compact summary table
# ------------------------------------------------------------


final_summary_s32 <-
  data.frame(
    
    metric =
      c(
        "GDSC2 drugs tested",
        "MPCDS FDR-significant",
        "MPCDS FDR: higher MPCDS -> predicted sensitivity",
        "MPCDS FDR: higher MPCDS -> predicted resistance",
        "Subtype FDR-significant",
        "Subtype FDR: PCD_C2 more sensitive",
        "Subtype FDR: PCD_C2 more resistant",
        "Cross-platform testable matched drugs",
        "Cross-platform directionally concordant",
        "Cross-platform directionally discordant",
        "Cross-platform concordance proportion",
        "Cross-platform concordance binomial p",
        "Cross-platform Spearman rho",
        "Cross-platform Spearman p",
        "GDSC2-FDR matched drugs",
        "GDSC2-FDR + PRISM nominal drugs"
      ),
    
    value =
      c(
        n_gdsc_drugs_s32,
        n_mpcds_fdr_s32,
        n_mpcds_sensitivity_s32,
        n_mpcds_resistance_s32,
        n_subtype_fdr_s32,
        n_c2_sensitive_s32,
        n_c2_resistant_s32,
        n_cross_platform_s32,
        n_cross_concordant_s32,
        n_cross_discordant_s32,
        cross_concordance_s32,
        cross_binom_s32$p.value,
        unname(
          cross_effect_test_s32$estimate
        ),
        cross_effect_test_s32$p.value,
        n_gdsc_fdr_matched_s32,
        n_gdsc_fdr_prism_nominal_s32
      ),
    
    stringsAsFactors = FALSE
  )


write.csv(
  final_summary_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_final_analysis_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7.4 Top association tables
# ------------------------------------------------------------


top_mpcds_s32 <-
  mpcds_drug_results_s32[
    order(
      mpcds_drug_results_s32$FDR,
      -abs(
        mpcds_drug_results_s32$rho
      )
    ),
    ,
    drop = FALSE
  ]


top_mpcds_s32 <-
  head(
    top_mpcds_s32,
    20
  )


write.csv(
  top_mpcds_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_top20_MPCDS_associations.csv"
  ),
  row.names = FALSE
)


top_subtype_s32 <-
  subtype_drug_results_s32[
    order(
      subtype_drug_results_s32$FDR,
      -abs(
        subtype_drug_results_s32[[
          "median_difference_C2_minus_C1"
        ]]
      )
    ),
    ,
    drop = FALSE
  ]


top_subtype_s32 <-
  head(
    top_subtype_s32,
    20
  )


write.csv(
  top_subtype_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_top20_PCD_subtype_associations.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7.5 Plotting packages
# ------------------------------------------------------------


library(
  ggplot2
)


library(
  patchwork
)


# ------------------------------------------------------------
# 7.6 Project colors
# ------------------------------------------------------------


PCD_C1_COLOR_s32 <-
  "#325074"


PCD_C2_COLOR_s32 <-
  "#FEC12C"


NEUTRAL_COLOR_s32 <-
  "grey55"


# ------------------------------------------------------------
# 7.7 Panel A
# MPCDS-GDSC2 landscape
# ------------------------------------------------------------


panel_a_data_s32 <-
  mpcds_drug_results_s32


panel_a_data_s32$significance <-
  ifelse(
    panel_a_data_s32$FDR <
      0.05,
    "FDR < 0.05",
    "Not FDR significant"
  )


panel_a_s32 <-
  ggplot(
    panel_a_data_s32,
    aes(
      x = rho,
      y =
        -log10(
          pmax(
            FDR,
            .Machine$double.xmin
          )
        )
    )
  ) +
  geom_point(
    aes(
      shape = significance
    ),
    size = 2.1,
    alpha = 0.75
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_hline(
    yintercept =
      -log10(
        0.05
      ),
    linetype = 3
  ) +
  labs(
    title =
      "A. MPCDS–GDSC2 association landscape",
    x =
      "Spearman rho",
    y =
      expression(
        -log[10](
          "BH FDR"
        )
      ),
    shape = NULL
  ) +
  theme_classic(
    base_size = 11
  )


# ------------------------------------------------------------
# 7.8 Panel B
# PCD subtype drug differences
# ------------------------------------------------------------


panel_b_data_s32 <-
  subtype_drug_results_s32


panel_b_data_s32$direction_group <-
  ifelse(
    
    panel_b_data_s32$FDR <
      0.05 &
      panel_b_data_s32[[
        "median_difference_C2_minus_C1"
      ]] <
      0,
    
    "PCD_C2 more sensitive",
    
    ifelse(
      
      panel_b_data_s32$FDR <
        0.05 &
        panel_b_data_s32[[
          "median_difference_C2_minus_C1"
        ]] >
        0,
      
      "PCD_C2 more resistant",
      
      "Not FDR significant"
    )
  )


panel_b_s32 <-
  ggplot(
    panel_b_data_s32,
    aes(
      x =
        median_difference_C2_minus_C1,
      y =
        -log10(
          pmax(
            FDR,
            .Machine$double.xmin
          )
        ),
      fill =
        direction_group
    )
  ) +
  geom_point(
    shape = 21,
    size = 2.3,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_hline(
    yintercept =
      -log10(
        0.05
      ),
    linetype = 3
  ) +
  scale_fill_manual(
    values =
      c(
        "PCD_C2 more sensitive" =
          PCD_C2_COLOR_s32,
        "PCD_C2 more resistant" =
          PCD_C1_COLOR_s32,
        "Not FDR significant" =
          NEUTRAL_COLOR_s32
      )
  ) +
  labs(
    title =
      "B. PCD subtype predicted IC50 differences",
    x =
      "Median predicted IC50: PCD_C2 − PCD_C1",
    y =
      expression(
        -log[10](
          "BH FDR"
        )
      ),
    fill = NULL
  ) +
  theme_classic(
    base_size = 11
  )


# ------------------------------------------------------------
# 7.9 Panel C
# GDSC2 vs PRISM effect sizes
# ------------------------------------------------------------


panel_c_s32 <-
  ggplot(
    cross_platform_s32,
    aes(
      x = GDSC_rho,
      y = PRISM_rho
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_point(
    size = 2.3,
    alpha = 0.75
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7
  ) +
  labs(
    title =
      "C. GDSC2 vs PRISM MPCDS effects",
    subtitle =
      paste0(
        "Spearman rho = ",
        sprintf(
          "%.2f",
          unname(
            cross_effect_test_s32$estimate
          )
        ),
        "; p = ",
        format.pval(
          cross_effect_test_s32$p.value,
          digits = 2
        )
      ),
    x =
      "GDSC2 MPCDS–predicted IC50 rho",
    y =
      "PRISM MPCDS–AUC rho"
  ) +
  theme_classic(
    base_size = 11
  )


# ------------------------------------------------------------
# 7.10 Panel D
# Directional concordance
# ------------------------------------------------------------


panel_d_data_s32 <-
  data.frame(
    
    direction =
      c(
        "Concordant",
        "Discordant"
      ),
    
    n =
      c(
        n_cross_concordant_s32,
        n_cross_discordant_s32
      ),
    
    stringsAsFactors = FALSE
  )


panel_d_s32 <-
  ggplot(
    panel_d_data_s32,
    aes(
      x = direction,
      y = n,
      fill = direction
    )
  ) +
  geom_col(
    width = 0.65
  ) +
  geom_text(
    aes(
      label = n
    ),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(
    values =
      c(
        "Concordant" =
          PCD_C1_COLOR_s32,
        "Discordant" =
          PCD_C2_COLOR_s32
      )
  ) +
  labs(
    title =
      "D. Cross-platform directional agreement",
    subtitle =
      paste0(
        "Concordance = ",
        sprintf(
          "%.1f%%",
          100 *
            cross_concordance_s32
        ),
        "; binomial p = ",
        format.pval(
          cross_binom_s32$p.value,
          digits = 2
        )
      ),
    x = NULL,
    y =
      "Number of matched drugs",
    fill = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none"
  )


# ------------------------------------------------------------
# 7.11 Combine Figure S15
# ------------------------------------------------------------
# Figure generation is intentionally deferred to Script 32.
# Script 25 produces analysis/result tables only.
# ------------------------------------------------------------


# ------------------------------------------------------------
# 7.12 Interpretation summary
# ------------------------------------------------------------


interpretation_s32 <-
  data.frame(
    
    topic =
      c(
        "GDSC2 MPCDS",
        "GDSC2 PCD subtype",
        "Cross-platform overlap",
        "Cross-platform concordance",
        "Cross-platform effect correlation",
        "Cross-platform replication",
        "Scientific conclusion"
      ),
    
    interpretation =
      c(
        
        paste0(
          n_mpcds_fdr_s32,
          " of ",
          n_gdsc_drugs_s32,
          " predicted GDSC2 drug responses were associated with MPCDS at BH FDR < 0.05; ",
          n_mpcds_sensitivity_s32,
          " indicated greater predicted sensitivity with higher MPCDS and ",
          n_mpcds_resistance_s32,
          " indicated greater predicted resistance."
        ),
        
        paste0(
          n_subtype_fdr_s32,
          " predicted drug responses differed between PCD_C1 and PCD_C2 at BH FDR < 0.05; ",
          n_c2_sensitive_s32,
          " had lower predicted IC50 in PCD_C2 and ",
          n_c2_resistant_s32,
          " had higher predicted IC50 in PCD_C2."
        ),
        
        paste0(
          n_cross_platform_s32,
          " unambiguous drugs were testable in both GDSC2 and PRISM."
        ),
        
        paste0(
          n_cross_concordant_s32,
          " of ",
          n_cross_platform_s32,
          " testable matched drugs were directionally concordant (",
          sprintf(
            "%.1f%%",
            100 *
              cross_concordance_s32
          ),
          "); concordance did not exceed chance expectation."
        ),
        
        paste0(
          "GDSC2 and PRISM MPCDS effect sizes showed negligible correlation (Spearman rho = ",
          sprintf(
            "%.3f",
            unname(
              cross_effect_test_s32$estimate
            )
          ),
          ", p = ",
          format.pval(
            cross_effect_test_s32$p.value,
            digits = 3
          ),
          ")."
        ),
        
        "No GDSC2 FDR-significant matched drug also reached nominal PRISM p < 0.05.",
        
        paste0(
          "GDSC2/oncoPredict findings therefore represent in silico drug-sensitivity hypotheses rather than independently validated therapeutic-response biomarkers."
        )
      ),
    
    stringsAsFactors = FALSE
  )


write.csv(
  interpretation_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_GDSC2_interpretation_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7.13 Analysis audit
# ------------------------------------------------------------


audit_s32 <-
  data.frame(
    
    analysis_component =
      c(
        "oncoPredict installation",
        "GDSC2 expression validation",
        "GDSC2 response validation",
        "GDSC2 TCGA-LUAD prediction",
        "MPCDS association analysis",
        "PCD subtype association analysis",
        "PRISM cross-platform comparison",
        "Cross-platform drug replication",
        "Clinical treatment-response validation"
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
        "No replicated GDSC2-FDR plus PRISM-nominal signal",
        "Not performed - no treated clinical response cohort"
      ),
    
    stringsAsFactors = FALSE
  )


write.csv(
  audit_s32,
  file.path(
    gdsc_results_dir_s32,
    "TCGA_LUAD_Script25_analysis_audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7.14 Final output-file verification
# ------------------------------------------------------------


final_files_s32 <-
  c(
    
    prediction_rds_s32,
    
    prediction_csv_s32,
    
    file.path(
      gdsc_results_dir_s32,
      "TCGA_LUAD_GDSC2_MPCDS_all_drug_associations.csv"
    ),
    
    file.path(
      gdsc_results_dir_s32,
      "TCGA_LUAD_GDSC2_PCD_subtype_all_drug_associations.csv"
    ),
    
    file.path(
      gdsc_results_dir_s32,
      "TCGA_LUAD_GDSC2_PRISM_MPCDS_cross_platform_testable_drugs.csv"
    ),
    
    file.path(
      gdsc_results_dir_s32,
      "TCGA_LUAD_GDSC2_final_analysis_summary.csv"
    ),
    
    file.path(
      gdsc_results_dir_s32,
      "TCGA_LUAD_GDSC2_top20_MPCDS_associations.csv"
    ),
    
    file.path(
      gdsc_results_dir_s32,
      "TCGA_LUAD_GDSC2_top20_PCD_subtype_associations.csv"
    ),
    
    file.path(
      gdsc_results_dir_s32,
      "TCGA_LUAD_GDSC2_interpretation_summary.csv"
    ),
    
    file.path(
      gdsc_results_dir_s32,
      "TCGA_LUAD_Script25_analysis_audit.csv"
    )
  )


cat(
  "\n========================================\n",
  "FINAL OUTPUT FILE CHECK\n",
  "========================================\n"
)


print(
  data.frame(
    
    file =
      final_files_s32,
    
    exists =
      file.exists(
        final_files_s32
      ),
    
    size_bytes =
      file.info(
        final_files_s32
      )$size,
    
    stringsAsFactors = FALSE
  ),
  
  row.names = FALSE
)


stopifnot(
  all(
    file.exists(
      final_files_s32
    )
  )
)


# ------------------------------------------------------------
# 7.15 Expected-result audit
# ------------------------------------------------------------


stopifnot(
  n_gdsc_drugs_s32 ==
    198
)


stopifnot(
  n_mpcds_fdr_s32 ==
    110
)


stopifnot(
  n_mpcds_sensitivity_s32 ==
    59
)


stopifnot(
  n_mpcds_resistance_s32 ==
    51
)


stopifnot(
  n_subtype_fdr_s32 ==
    104
)


stopifnot(
  n_c2_sensitive_s32 ==
    38
)


stopifnot(
  n_c2_resistant_s32 ==
    66
)


stopifnot(
  n_cross_platform_s32 ==
    72
)


stopifnot(
  n_cross_concordant_s32 ==
    37
)


stopifnot(
  n_cross_discordant_s32 ==
    35
)


stopifnot(
  n_gdsc_fdr_matched_s32 ==
    39
)


stopifnot(
  n_gdsc_fdr_prism_nominal_s32 ==
    0
)


# ------------------------------------------------------------
# 7.16 Final scientific interpretation
# ------------------------------------------------------------


# ------------------------------------------------------------
# Regression checks against established GDSC2/oncoPredict results
# ------------------------------------------------------------

stopifnot(
  n_gdsc_drugs_s32 == 198,
  n_mpcds_fdr_s32 == 110,
  n_mpcds_sensitivity_s32 == 59,
  n_mpcds_resistance_s32 == 51,
  n_subtype_fdr_s32 == 104,
  n_c2_sensitive_s32 == 38,
  n_c2_resistant_s32 == 66
)

cat(
  "\nEstablished GDSC2 regression checks: PASSED\n"
)


cat(
  "\n========================================\n",
  "SCRIPT 25 SCIENTIFIC CONCLUSION\n",
  "========================================\n"
)


cat(
  paste0(
    "\nGDSC2/oncoPredict identified multiple transcriptome-derived ",
    "predicted drug-sensitivity associations with MPCDS and PCD ",
    "subtype. However, these associations did not show meaningful ",
    "cross-platform replication in the independent DepMap/PRISM ",
    "cell-line screen. Therefore, the GDSC2 findings should be ",
    "interpreted as in silico therapeutic hypotheses rather than ",
    "validated treatment-response biomarkers.\n"
  )
)


# ------------------------------------------------------------
# 7.17 Success marker
# ------------------------------------------------------------


cat(
  "\n========================================\n",
  "SCRIPT 25 GDSC2 / ONCOPREDICT ANALYSIS COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)