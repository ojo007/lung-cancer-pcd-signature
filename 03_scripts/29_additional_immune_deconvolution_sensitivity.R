# ============================================================
# SCRIPT 29
# Additional immune-deconvolution sensitivity analysis
#
# TCGA-LUAD
#
# Methods:
#   xCell
#   quanTIseq
#   TIMER
#   EPIC
#   CIBERSORT / LM22
#
# Analyses:
#   1. Full-cohort immune deconvolution
#   2. PCD_C1 versus PCD_C2 associations
#   3. Continuous MPCDS associations
#   4. Cross-method biological mapping
#   5. Cross-method replication analysis
#   6. MCP-counter supportive comparison
#   7. Final summary figure and tables
# ============================================================


source("03_scripts/00_project_config.R")


cat(
  "\n========================================\n",
  "SCRIPT 29\n",
  "ADDITIONAL IMMUNE-DECONVOLUTION SENSITIVITY ANALYSIS\n",
  "========================================\n"
)


# ============================================================
# SECTION 1
# Required packages
# ============================================================


required_packages_s34 <-
  c(
    "IOBR"
  )


package_available_s34 <-
  vapply(
    required_packages_s34,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )


cat(
  "\nPackage availability:\n"
)


print(
  data.frame(
    package =
      required_packages_s34,
    available =
      unname(
        package_available_s34
      ),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


stopifnot(
  all(
    package_available_s34
  )
)


cat(
  "\nIOBR version:\n"
)


print(
  as.character(
    packageVersion(
      "IOBR"
    )
  )
)


# ============================================================
# SECTION 2
# Paths
# ============================================================


patient_tpm_file_s34 <-
  "02_processed_data/immunotherapy/TCGA_LUAD_patient_protein_coding_TPM.rds"


subtype_file_s34 <-
  "04_results/clustering/TCGA_LUAD_PCD_cluster_assignments.csv"


mpcds_file_s34 <-
  "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv"


mcp_file_s34 <-
  "04_results/immune/TCGA_LUAD_LUSC_MCPcounter_robust_populations.csv"


results_dir_s34 <-
  "04_results/immune/additional_deconvolution_sensitivity"


processed_dir_s34 <-
  "02_processed_data/immune/additional_deconvolution_sensitivity"



dir.create(
  results_dir_s34,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  processed_dir_s34,
  recursive = TRUE,
  showWarnings = FALSE
)



required_input_files_s34 <-
  c(
    patient_tpm_file_s34,
    subtype_file_s34,
    mpcds_file_s34,
    mcp_file_s34
  )


cat(
  "\nInput-file check:\n"
)


print(
  data.frame(
    file =
      required_input_files_s34,
    exists =
      file.exists(
        required_input_files_s34
      ),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


stopifnot(
  all(
    file.exists(
      required_input_files_s34
    )
  )
)


# ============================================================
# SECTION 3
# Load and validate input data
# ============================================================


patient_tpm_s34 <-
  as.matrix(
    readRDS(
      patient_tpm_file_s34
    )
  )


subtype_s34 <-
  read.csv(
    subtype_file_s34,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


mpcds_s34 <-
  read.csv(
    mpcds_file_s34,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


mcp_existing_s34 <-
  read.csv(
    mcp_file_s34,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat(
  "\n========================================\n",
  "INPUT QC\n",
  "========================================\n"
)


cat(
  "\nTPM dimensions:\n"
)


print(
  dim(
    patient_tpm_s34
  )
)


cat(
  "\nTPM range:\n"
)


print(
  range(
    patient_tpm_s34,
    na.rm = TRUE
  )
)


cat(
  "\nSubtype dimensions:\n"
)


print(
  dim(
    subtype_s34
  )
)


cat(
  "\nMPCDS dimensions:\n"
)


print(
  dim(
    mpcds_s34
  )
)


cat(
  "\nMCP-counter dimensions:\n"
)


print(
  dim(
    mcp_existing_s34
  )
)


stopifnot(
  identical(
    dim(
      patient_tpm_s34
    ),
    c(
      19938L,
      517L
    )
  )
)


stopifnot(
  !anyNA(
    patient_tpm_s34
  )
)


stopifnot(
  all(
    is.finite(
      patient_tpm_s34
    )
  )
)


stopifnot(
  !any(
    duplicated(
      rownames(
        patient_tpm_s34
      )
    )
  )
)


stopifnot(
  !any(
    duplicated(
      colnames(
        patient_tpm_s34
      )
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
        subtype_s34
      )
  )
)


stopifnot(
  all(
    c(
      "patient_id",
      "MPCDS",
      "MPCDS_z"
    ) %in%
      colnames(
        mpcds_s34
      )
  )
)


stopifnot(
  nrow(
    subtype_s34
  ) ==
    517L
)


stopifnot(
  nrow(
    mpcds_s34
  ) ==
    504L
)


stopifnot(
  nrow(
    mcp_existing_s34
  ) ==
    7L
)


stopifnot(
  !any(
    duplicated(
      subtype_s34[[
        "patient_id"
      ]]
    )
  )
)


stopifnot(
  !any(
    duplicated(
      mpcds_s34[[
        "patient_id"
      ]]
    )
  )
)


expected_patients_s34 <-
  colnames(
    patient_tpm_s34
  )


stopifnot(
  setequal(
    expected_patients_s34,
    subtype_s34[[
      "patient_id"
    ]]
  )
)


stopifnot(
  length(
    intersect(
      expected_patients_s34,
      mpcds_s34[[
        "patient_id"
      ]]
    )
  ) ==
    504L
)


cat(
  "\nSubtype counts:\n"
)


print(
  table(
    subtype_s34[[
      "PCD_cluster"
    ]]
  )
)


stopifnot(
  sum(
    subtype_s34[[
      "PCD_cluster"
    ]] ==
      "PCD_C1"
  ) ==
    284L
)


stopifnot(
  sum(
    subtype_s34[[
      "PCD_cluster"
    ]] ==
      "PCD_C2"
  ) ==
    233L
)


# ============================================================
# SECTION 4
# Helper for full-cohort deconvolution methods
# ============================================================


run_deconvolution_s34 <-
  function(
    method_name_s34,
    method_function_s34
  ) {
    
    warning_messages_s34 <-
      character(0)
    
    
    start_time_s34 <-
      Sys.time()
    
    
    result_s34 <-
      tryCatch(
        
        withCallingHandlers(
          
          method_function_s34(),
          
          warning =
            function(w_s34) {
              
              warning_messages_s34 <<-
                c(
                  warning_messages_s34,
                  conditionMessage(
                    w_s34
                  )
                )
              
              
              invokeRestart(
                "muffleWarning"
              )
            }
        ),
        
        error =
          function(e_s34) {
            
            structure(
              list(
                message =
                  conditionMessage(
                    e_s34
                  )
              ),
              class =
                "s34_method_error"
            )
          }
      )
    
    
    elapsed_s34 <-
      as.numeric(
        difftime(
          Sys.time(),
          start_time_s34,
          units =
            "secs"
        )
      )
    
    
    success_s34 <-
      !inherits(
        result_s34,
        "s34_method_error"
      )
    
    
    if (
      !success_s34
    ) {
      
      stop(
        paste0(
          method_name_s34,
          " failed: ",
          result_s34$message
        )
      )
    }
    
    
    result_s34 <-
      as.data.frame(
        result_s34,
        check.names = FALSE
      )
    
    
    stopifnot(
      nrow(
        result_s34
      ) ==
        517L
    )
    
    
    stopifnot(
      "ID" %in%
        colnames(
          result_s34
        )
    )
    
    
    stopifnot(
      !any(
        duplicated(
          result_s34[[
            "ID"
          ]]
        )
      )
    )
    
    
    stopifnot(
      setequal(
        result_s34[[
          "ID"
        ]],
        expected_patients_s34
      )
    )
    
    
    result_s34 <-
      result_s34[
        match(
          expected_patients_s34,
          result_s34[[
            "ID"
          ]]
        ),
        ,
        drop = FALSE
      ]
    
    
    stopifnot(
      identical(
        result_s34[[
          "ID"
        ]],
        expected_patients_s34
      )
    )
    
    
    list(
      result =
        result_s34,
      elapsed_seconds =
        elapsed_s34,
      warnings =
        warning_messages_s34
    )
  }


# ============================================================
# SECTION 5
# Full xCell analysis
# ============================================================


cat(
  "\n========================================\n",
  "RUNNING xCell\n",
  "========================================\n"
)


xcell_run_s34 <-
  run_deconvolution_s34(
    method_name_s34 =
      "xCell",
    method_function_s34 =
      function() {
        
        IOBR::deconvo_xcell(
          eset =
            patient_tpm_s34,
          project =
            "TCGA-LUAD",
          arrays =
            FALSE
        )
      }
  )


xcell_s34 <-
  xcell_run_s34$result


write.csv(
  xcell_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_xCell_full_scores.csv"
  ),
  row.names = FALSE
)


saveRDS(
  xcell_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_xCell_full_scores.rds"
  )
)


# ============================================================
# SECTION 6
# Full quanTIseq analysis
# ============================================================


cat(
  "\n========================================\n",
  "RUNNING quanTIseq\n",
  "========================================\n"
)


quantiseq_run_s34 <-
  run_deconvolution_s34(
    method_name_s34 =
      "quanTIseq",
    method_function_s34 =
      function() {
        
        IOBR::deconvo_quantiseq(
          eset =
            patient_tpm_s34,
          project =
            "TCGA-LUAD",
          tumor =
            TRUE,
          arrays =
            FALSE,
          scale_mrna =
            TRUE
        )
      }
  )


quantiseq_s34 <-
  quantiseq_run_s34$result


write.csv(
  quantiseq_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_quanTIseq_full_scores.csv"
  ),
  row.names = FALSE
)


saveRDS(
  quantiseq_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_quanTIseq_full_scores.rds"
  )
)


# ============================================================
# SECTION 7
# Full TIMER analysis
# ============================================================


cat(
  "\n========================================\n",
  "RUNNING TIMER\n",
  "========================================\n"
)


timer_indications_s34 <-
  rep(
    "luad",
    ncol(
      patient_tpm_s34
    )
  )


timer_run_s34 <-
  run_deconvolution_s34(
    method_name_s34 =
      "TIMER",
    method_function_s34 =
      function() {
        
        IOBR::deconvo_timer(
          eset =
            patient_tpm_s34,
          project =
            "TCGA-LUAD",
          indications =
            timer_indications_s34
        )
      }
  )


timer_s34 <-
  timer_run_s34$result


write.csv(
  timer_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_TIMER_full_scores.csv"
  ),
  row.names = FALSE
)


saveRDS(
  timer_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_TIMER_full_scores.rds"
  )
)


# ============================================================
# SECTION 8
# Full EPIC analysis
# ============================================================


cat(
  "\n========================================\n",
  "RUNNING EPIC\n",
  "========================================\n"
)


epic_run_s34 <-
  run_deconvolution_s34(
    method_name_s34 =
      "EPIC",
    method_function_s34 =
      function() {
        
        IOBR::deconvo_epic(
          eset =
            patient_tpm_s34,
          project =
            "TCGA-LUAD",
          tumor =
            TRUE
        )
      }
  )


epic_s34 <-
  epic_run_s34$result


write.csv(
  epic_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_EPIC_full_scores.csv"
  ),
  row.names = FALSE
)


saveRDS(
  epic_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_EPIC_full_scores.rds"
  )
)


# ============================================================
# SECTION 9
# Full CIBERSORT / LM22 analysis
#
# Setting fixed from the preceding technical sensitivity test:
# relative mode, 100 permutations.
#
# The earlier 10-sample comparison showed identical fractions
# for 100 and 1000 permutations. That technical benchmark is
# not repeated in this production script.
# ============================================================


cat(
  "\n========================================\n",
  "RUNNING CIBERSORT / LM22\n",
  "========================================\n"
)


cibersort_permutations_s34 <-
  100L


cibersort_seed_s34 <-
  34002L


cibersort_warning_s34 <-
  character(0)


cibersort_start_s34 <-
  Sys.time()


cibersort_s34 <-
  tryCatch(
    
    withCallingHandlers(
      
      IOBR::deconvo_cibersort(
        eset =
          patient_tpm_s34,
        project =
          "TCGA-LUAD",
        arrays =
          FALSE,
        perm =
          cibersort_permutations_s34,
        absolute =
          FALSE,
        parallel =
          FALSE,
        seed =
          cibersort_seed_s34
      ),
      
      warning =
        function(w_s34) {
          
          cibersort_warning_s34 <<-
            c(
              cibersort_warning_s34,
              conditionMessage(
                w_s34
              )
            )
          
          
          invokeRestart(
            "muffleWarning"
          )
        }
    ),
    
    error =
      function(e_s34) {
        
        stop(
          paste0(
            "CIBERSORT failed: ",
            conditionMessage(
              e_s34
            )
          )
        )
      }
  )


cibersort_elapsed_s34 <-
  as.numeric(
    difftime(
      Sys.time(),
      cibersort_start_s34,
      units =
        "secs"
    )
  )


cibersort_s34 <-
  as.data.frame(
    cibersort_s34,
    check.names = FALSE
  )


stopifnot(
  nrow(
    cibersort_s34
  ) ==
    517L
)


stopifnot(
  "ID" %in%
    colnames(
      cibersort_s34
    )
)


stopifnot(
  setequal(
    cibersort_s34[[
      "ID"
    ]],
    expected_patients_s34
  )
)


cibersort_s34 <-
  cibersort_s34[
    match(
      expected_patients_s34,
      cibersort_s34[[
        "ID"
      ]]
    ),
    ,
    drop = FALSE
  ]


stopifnot(
  identical(
    cibersort_s34[[
      "ID"
    ]],
    expected_patients_s34
  )
)


cibersort_diagnostic_columns_s34 <-
  intersect(
    c(
      "P-value_CIBERSORT",
      "Correlation_CIBERSORT",
      "RMSE_CIBERSORT"
    ),
    colnames(
      cibersort_s34
    )
  )


cibersort_feature_columns_s34 <-
  setdiff(
    colnames(
      cibersort_s34
    ),
    c(
      "ID",
      "ProjectID",
      cibersort_diagnostic_columns_s34
    )
  )


stopifnot(
  length(
    cibersort_feature_columns_s34
  ) ==
    22L
)


cibersort_fraction_matrix_s34 <-
  as.matrix(
    cibersort_s34[
      ,
      cibersort_feature_columns_s34,
      drop = FALSE
    ]
  )


stopifnot(
  !anyNA(
    cibersort_fraction_matrix_s34
  )
)


stopifnot(
  all(
    is.finite(
      cibersort_fraction_matrix_s34
    )
  )
)


stopifnot(
  min(
    cibersort_fraction_matrix_s34
  ) >=
    0
)


stopifnot(
  max(
    cibersort_fraction_matrix_s34
  ) <=
    1
)


cibersort_fraction_sums_s34 <-
  rowSums(
    cibersort_fraction_matrix_s34
  )


cibersort_max_sum_deviation_s34 <-
  max(
    abs(
      cibersort_fraction_sums_s34 -
        1
    )
  )


stopifnot(
  cibersort_max_sum_deviation_s34 <
    1e-6
)


write.csv(
  cibersort_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_CIBERSORT_LM22_full_scores.csv"
  ),
  row.names = FALSE
)


saveRDS(
  cibersort_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_CIBERSORT_LM22_full_scores.rds"
  )
)


if (
  length(
    cibersort_diagnostic_columns_s34
  ) >
  0L
) {
  
  write.csv(
    cibersort_s34[
      ,
      c(
        "ID",
        intersect(
          "ProjectID",
          colnames(
            cibersort_s34
          )
        ),
        cibersort_diagnostic_columns_s34
      ),
      drop = FALSE
    ],
    file.path(
      results_dir_s34,
      "TCGA_LUAD_CIBERSORT_LM22_patient_diagnostics.csv"
    ),
    row.names = FALSE
  )
}


# ============================================================
# SECTION 10
# Full-method audit
# ============================================================


method_dimensions_s34 <-
  data.frame(
    method =
      c(
        "xCell",
        "quanTIseq",
        "TIMER",
        "EPIC",
        "CIBERSORT"
      ),
    patients =
      c(
        nrow(
          xcell_s34
        ),
        nrow(
          quantiseq_s34
        ),
        nrow(
          timer_s34
        ),
        nrow(
          epic_s34
        ),
        nrow(
          cibersort_s34
        )
      ),
    output_columns =
      c(
        ncol(
          xcell_s34
        ),
        ncol(
          quantiseq_s34
        ),
        ncol(
          timer_s34
        ),
        ncol(
          epic_s34
        ),
        ncol(
          cibersort_s34
        )
      ),
    stringsAsFactors = FALSE
  )


cat(
  "\n========================================\n",
  "FULL DECONVOLUTION DIMENSIONS\n",
  "========================================\n"
)


print(
  method_dimensions_s34,
  row.names = FALSE
)


stopifnot(
  all(
    method_dimensions_s34[[
      "patients"
    ]] ==
      517L
  )
)


method_runtime_audit_s34 <-
  data.frame(
    method =
      c(
        "xCell",
        "quanTIseq",
        "TIMER",
        "EPIC",
        "CIBERSORT"
      ),
    elapsed_seconds =
      c(
        xcell_run_s34$elapsed_seconds,
        quantiseq_run_s34$elapsed_seconds,
        timer_run_s34$elapsed_seconds,
        epic_run_s34$elapsed_seconds,
        cibersort_elapsed_s34
      ),
    warning_count =
      c(
        length(
          xcell_run_s34$warnings
        ),
        length(
          quantiseq_run_s34$warnings
        ),
        length(
          timer_run_s34$warnings
        ),
        length(
          epic_run_s34$warnings
        ),
        length(
          cibersort_warning_s34
        )
      ),
    stringsAsFactors = FALSE
  )


write.csv(
  method_dimensions_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_full_deconvolution_method_dimensions.csv"
  ),
  row.names = FALSE
)


write.csv(
  method_runtime_audit_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_full_deconvolution_method_audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# SECTION 11
# Define immune feature columns
# ============================================================


xcell_feature_columns_s34 <-
  setdiff(
    colnames(
      xcell_s34
    ),
    c(
      "ID",
      "ProjectID"
    )
  )


quantiseq_feature_columns_s34 <-
  setdiff(
    colnames(
      quantiseq_s34
    ),
    c(
      "ID",
      "ProjectID"
    )
  )


timer_feature_columns_s34 <-
  setdiff(
    colnames(
      timer_s34
    ),
    c(
      "ID",
      "ProjectID"
    )
  )


epic_feature_columns_s34 <-
  setdiff(
    colnames(
      epic_s34
    ),
    c(
      "ID",
      "ProjectID"
    )
  )


stopifnot(
  length(
    xcell_feature_columns_s34
  ) ==
    67L
)


stopifnot(
  length(
    quantiseq_feature_columns_s34
  ) ==
    11L
)


stopifnot(
  length(
    timer_feature_columns_s34
  ) ==
    6L
)


stopifnot(
  length(
    epic_feature_columns_s34
  ) ==
    8L
)


stopifnot(
  length(
    cibersort_feature_columns_s34
  ) ==
    22L
)


feature_inventory_s34 <-
  rbind(
    data.frame(
      method =
        "xCell",
      feature =
        xcell_feature_columns_s34,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method =
        "quanTIseq",
      feature =
        quantiseq_feature_columns_s34,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method =
        "TIMER",
      feature =
        timer_feature_columns_s34,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method =
        "EPIC",
      feature =
        epic_feature_columns_s34,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method =
        "CIBERSORT",
      feature =
        cibersort_feature_columns_s34,
      stringsAsFactors = FALSE
    )
  )


stopifnot(
  nrow(
    feature_inventory_s34
  ) ==
    114L
)


# ============================================================
# SECTION 12
# Feature-level QC
# ============================================================


method_tables_s34 <-
  list(
    xCell =
      xcell_s34,
    quanTIseq =
      quantiseq_s34,
    TIMER =
      timer_s34,
    EPIC =
      epic_s34,
    CIBERSORT =
      cibersort_s34
  )


feature_qc_s34 <-
  do.call(
    rbind,
    lapply(
      names(
        method_tables_s34
      ),
      function(method_s34) {
        
        x_s34 <-
          method_tables_s34[[
            method_s34
          ]]
        
        
        features_s34 <-
          feature_inventory_s34[[
            "feature"
          ]][
            feature_inventory_s34[[
              "method"
            ]] ==
              method_s34
          ]
        
        
        do.call(
          rbind,
          lapply(
            features_s34,
            function(feature_s34) {
              
              values_s34 <-
                x_s34[[
                  feature_s34
                ]]
              
              
              data.frame(
                method =
                  method_s34,
                feature =
                  feature_s34,
                numeric =
                  is.numeric(
                    values_s34
                  ),
                n_NA =
                  sum(
                    is.na(
                      values_s34
                    )
                  ),
                n_nonfinite =
                  if (
                    is.numeric(
                      values_s34
                    )
                  ) {
                    sum(
                      !is.finite(
                        values_s34
                      )
                    )
                  } else {
                    NA_integer_
                  },
                n_unique =
                  length(
                    unique(
                      values_s34
                    )
                  ),
                minimum =
                  min(
                    values_s34,
                    na.rm = TRUE
                  ),
                maximum =
                  max(
                    values_s34,
                    na.rm = TRUE
                  ),
                stringsAsFactors = FALSE
              )
            }
          )
        )
      }
    )
  )


rownames(
  feature_qc_s34
) <-
  NULL


stopifnot(
  all(
    feature_qc_s34[[
      "numeric"
    ]]
  )
)


stopifnot(
  all(
    feature_qc_s34[[
      "n_NA"
    ]] ==
      0L
  )
)


stopifnot(
  all(
    feature_qc_s34[[
      "n_nonfinite"
    ]] ==
      0L
  )
)


feature_qc_s34$constant <-
  feature_qc_s34[[
    "n_unique"
  ]] <=
  1L


feature_qc_s34$very_low_variation <-
  feature_qc_s34[[
    "n_unique"
  ]] <=
  5L


feature_qc_s34$test_eligible <-
  !feature_qc_s34[[
    "constant"
  ]] &
  !feature_qc_s34[[
    "very_low_variation"
  ]]


feature_qc_s34$exclusion_reason <-
  ifelse(
    feature_qc_s34[[
      "constant"
    ]],
    "constant_feature",
    ifelse(
      feature_qc_s34[[
        "very_low_variation"
      ]],
      "very_low_variation_le_5_unique_values",
      NA_character_
    )
  )


cat(
  "\n========================================\n",
  "FEATURE QC\n",
  "========================================\n"
)


cat(
  "\nEligible features:\n"
)


print(
  sum(
    feature_qc_s34[[
      "test_eligible"
    ]]
  )
)


cat(
  "\nExcluded features:\n"
)


print(
  feature_qc_s34[
    !feature_qc_s34[[
      "test_eligible"
    ]],
    c(
      "method",
      "feature",
      "n_unique",
      "exclusion_reason"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


stopifnot(
  sum(
    feature_qc_s34[[
      "test_eligible"
    ]]
  ) ==
    113L
)


excluded_feature_s34 <-
  feature_qc_s34[[
    "feature"
  ]][
    !feature_qc_s34[[
      "test_eligible"
    ]]
  ]


stopifnot(
  length(
    excluded_feature_s34
  ) ==
    1L
)


stopifnot(
  excluded_feature_s34 ==
    "T_cells_CD4_naive_CIBERSORT"
)


write.csv(
  feature_inventory_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_multimethod_immune_feature_inventory.csv"
  ),
  row.names = FALSE
)


write.csv(
  feature_qc_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_multimethod_immune_feature_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# SECTION 13
# Harmonized long-format immune data
# ============================================================


method_to_long_s34 <-
  function(
    x_s34,
    method_s34,
    feature_columns_s34
  ) {
    
    do.call(
      rbind,
      lapply(
        feature_columns_s34,
        function(feature_s34) {
          
          data.frame(
            patient_id =
              x_s34[[
                "ID"
              ]],
            method =
              method_s34,
            feature =
              feature_s34,
            score =
              x_s34[[
                feature_s34
              ]],
            stringsAsFactors = FALSE
          )
        }
      )
    )
  }


immune_long_s34 <-
  rbind(
    method_to_long_s34(
      xcell_s34,
      "xCell",
      xcell_feature_columns_s34
    ),
    method_to_long_s34(
      quantiseq_s34,
      "quanTIseq",
      quantiseq_feature_columns_s34
    ),
    method_to_long_s34(
      timer_s34,
      "TIMER",
      timer_feature_columns_s34
    ),
    method_to_long_s34(
      epic_s34,
      "EPIC",
      epic_feature_columns_s34
    ),
    method_to_long_s34(
      cibersort_s34,
      "CIBERSORT",
      cibersort_feature_columns_s34
    )
  )


stopifnot(
  nrow(
    immune_long_s34
  ) ==
    517L *
    114L
)


immune_subtype_s34 <-
  merge(
    immune_long_s34,
    subtype_s34[
      ,
      c(
        "patient_id",
        "PCD_cluster"
      ),
      drop = FALSE
    ],
    by =
      "patient_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )


stopifnot(
  !anyNA(
    immune_subtype_s34[[
      "PCD_cluster"
    ]]
  )
)


immune_mpcds_s34 <-
  merge(
    immune_long_s34,
    mpcds_s34[
      ,
      c(
        "patient_id",
        "MPCDS",
        "MPCDS_z"
      ),
      drop = FALSE
    ],
    by =
      "patient_id",
    all =
      FALSE,
    sort =
      FALSE
  )


stopifnot(
  nrow(
    immune_mpcds_s34
  ) ==
    504L *
    114L
)


saveRDS(
  immune_long_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_multimethod_immune_long_517.rds"
  )
)


saveRDS(
  immune_subtype_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_multimethod_immune_subtype_517.rds"
  )
)


saveRDS(
  immune_mpcds_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_multimethod_immune_MPCDS_504.rds"
  )
)


# ============================================================
# SECTION 14
# PCD_C1 versus PCD_C2 associations
#
# Positive rank-biserial:
#   higher in PCD_C2
#
# Negative rank-biserial:
#   higher in PCD_C1
# ============================================================


run_subtype_test_s34 <-
  function(
    x_s34,
    eligible_s34
  ) {
    
    c1_s34 <-
      x_s34[[
        "score"
      ]][
        x_s34[[
          "PCD_cluster"
        ]] ==
          "PCD_C1"
      ]
    
    
    c2_s34 <-
      x_s34[[
        "score"
      ]][
        x_s34[[
          "PCD_cluster"
        ]] ==
          "PCD_C2"
      ]
    
    
    stopifnot(
      length(
        c1_s34
      ) ==
        284L
    )
    
    
    stopifnot(
      length(
        c2_s34
      ) ==
        233L
    )
    
    
    c1_mean_s34 <-
      mean(
        c1_s34
      )
    
    
    c2_mean_s34 <-
      mean(
        c2_s34
      )
    
    
    c1_median_s34 <-
      median(
        c1_s34
      )
    
    
    c2_median_s34 <-
      median(
        c2_s34
      )
    
    
    if (
      !eligible_s34
    ) {
      
      return(
        data.frame(
          n_C1 =
            length(
              c1_s34
            ),
          n_C2 =
            length(
              c2_s34
            ),
          C1_mean =
            c1_mean_s34,
          C2_mean =
            c2_mean_s34,
          mean_difference_C2_minus_C1 =
            c2_mean_s34 -
            c1_mean_s34,
          C1_median =
            c1_median_s34,
          C2_median =
            c2_median_s34,
          median_difference_C2_minus_C1 =
            c2_median_s34 -
            c1_median_s34,
          wilcoxon_W =
            NA_real_,
          rank_biserial_C2_vs_C1 =
            NA_real_,
          p_value =
            NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
    
    
    test_s34 <-
      suppressWarnings(
        wilcox.test(
          c2_s34,
          c1_s34,
          alternative =
            "two.sided",
          exact =
            FALSE,
          correct =
            FALSE
        )
      )
    
    
    U_s34 <-
      unname(
        test_s34$statistic
      )
    
    
    rank_biserial_s34 <-
      2 *
      U_s34 /
      (
        length(
          c2_s34
        ) *
          length(
            c1_s34
          )
      ) -
      1
    
    
    data.frame(
      n_C1 =
        length(
          c1_s34
        ),
      n_C2 =
        length(
          c2_s34
        ),
      C1_mean =
        c1_mean_s34,
      C2_mean =
        c2_mean_s34,
      mean_difference_C2_minus_C1 =
        c2_mean_s34 -
        c1_mean_s34,
      C1_median =
        c1_median_s34,
      C2_median =
        c2_median_s34,
      median_difference_C2_minus_C1 =
        c2_median_s34 -
        c1_median_s34,
      wilcoxon_W =
        U_s34,
      rank_biserial_C2_vs_C1 =
        rank_biserial_s34,
      p_value =
        test_s34$p.value,
      stringsAsFactors = FALSE
    )
  }


subtype_split_s34 <-
  split(
    immune_subtype_s34,
    interaction(
      immune_subtype_s34[[
        "method"
      ]],
      immune_subtype_s34[[
        "feature"
      ]],
      drop = TRUE,
      lex.order = TRUE
    )
  )


stopifnot(
  length(
    subtype_split_s34
  ) ==
    114L
)


subtype_results_s34 <-
  do.call(
    rbind,
    lapply(
      subtype_split_s34,
      function(x_s34) {
        
        method_s34 <-
          unique(
            x_s34[[
              "method"
            ]]
          )
        
        
        feature_s34 <-
          unique(
            x_s34[[
              "feature"
            ]]
          )
        
        
        qc_s34 <-
          feature_qc_s34[
            feature_qc_s34[[
              "method"
            ]] ==
              method_s34 &
              feature_qc_s34[[
                "feature"
              ]] ==
              feature_s34,
            ,
            drop = FALSE
          ]
        
        
        stopifnot(
          nrow(
            qc_s34
          ) ==
            1L
        )
        
        
        test_s34 <-
          run_subtype_test_s34(
            x_s34,
            qc_s34[[
              "test_eligible"
            ]]
          )
        
        
        data.frame(
          method =
            method_s34,
          feature =
            feature_s34,
          n_unique =
            qc_s34[[
              "n_unique"
            ]],
          test_eligible =
            qc_s34[[
              "test_eligible"
            ]],
          exclusion_reason =
            qc_s34[[
              "exclusion_reason"
            ]],
          test_s34,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )
  )


rownames(
  subtype_results_s34
) <-
  NULL


subtype_results_s34$FDR_BH <-
  NA_real_


for (
  method_s34 in unique(
    subtype_results_s34[[
      "method"
    ]]
  )
) {
  
  index_s34 <-
    which(
      subtype_results_s34[[
        "method"
      ]] ==
        method_s34 &
        subtype_results_s34[[
          "test_eligible"
        ]] &
        !is.na(
          subtype_results_s34[[
            "p_value"
          ]]
        )
    )
  
  
  subtype_results_s34[[
    "FDR_BH"
  ]][
    index_s34
  ] <-
    p.adjust(
      subtype_results_s34[[
        "p_value"
      ]][
        index_s34
      ],
      method =
        "BH"
    )
}


subtype_results_s34$significant_FDR_0_05 <-
  !is.na(
    subtype_results_s34[[
      "FDR_BH"
    ]]
  ) &
  subtype_results_s34[[
    "FDR_BH"
  ]] <
  0.05


subtype_results_s34$significant_FDR_0_10 <-
  !is.na(
    subtype_results_s34[[
      "FDR_BH"
    ]]
  ) &
  subtype_results_s34[[
    "FDR_BH"
  ]] <
  0.10


subtype_results_s34$direction <-
  ifelse(
    !subtype_results_s34[[
      "test_eligible"
    ]],
    "Not tested",
    ifelse(
      subtype_results_s34[[
        "rank_biserial_C2_vs_C1"
      ]] >
        0,
      "Higher in PCD_C2",
      ifelse(
        subtype_results_s34[[
          "rank_biserial_C2_vs_C1"
        ]] <
          0,
        "Higher in PCD_C1",
        "No difference"
      )
    )
  )


cat(
  "\n========================================\n",
  "PCD SUBTYPE ASSOCIATION SUMMARY\n",
  "========================================\n"
)


subtype_method_summary_s34 <-
  do.call(
    rbind,
    lapply(
      c(
        "xCell",
        "quanTIseq",
        "TIMER",
        "EPIC",
        "CIBERSORT"
      ),
      function(method_s34) {
        
        x_s34 <-
          subtype_results_s34[
            subtype_results_s34[[
              "method"
            ]] ==
              method_s34,
            ,
            drop = FALSE
          ]
        
        
        data.frame(
          method =
            method_s34,
          n_features =
            nrow(
              x_s34
            ),
          n_tested =
            sum(
              x_s34[[
                "test_eligible"
              ]]
            ),
          FDR_lt_0_05 =
            sum(
              x_s34[[
                "significant_FDR_0_05"
              ]]
            ),
          FDR_lt_0_10 =
            sum(
              x_s34[[
                "significant_FDR_0_10"
              ]]
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )


print(
  subtype_method_summary_s34,
  row.names = FALSE
)


cat(
  "\nTotal subtype FDR < 0.05:\n"
)


print(
  sum(
    subtype_results_s34[[
      "significant_FDR_0_05"
    ]]
  )
)


write.csv(
  subtype_results_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_multimethod_PCD_subtype_immune_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  subtype_results_s34[
    subtype_results_s34[[
      "significant_FDR_0_05"
    ]],
    ,
    drop = FALSE
  ],
  file.path(
    results_dir_s34,
    "TCGA_LUAD_multimethod_PCD_subtype_immune_FDR05.csv"
  ),
  row.names = FALSE
)


write.csv(
  subtype_method_summary_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_multimethod_PCD_subtype_method_summary.csv"
  ),
  row.names = FALSE
)


saveRDS(
  subtype_results_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_multimethod_PCD_subtype_immune_associations.rds"
  )
)


# ============================================================
# SECTION 15
# Continuous MPCDS associations
# ============================================================


run_mpcds_test_s34 <-
  function(
    x_s34,
    eligible_s34
  ) {
    
    immune_s34 <-
      x_s34[[
        "score"
      ]]
    
    
    mpcds_z_s34 <-
      x_s34[[
        "MPCDS_z"
      ]]
    
    
    stopifnot(
      length(
        immune_s34
      ) ==
        504L
    )
    
    
    if (
      !eligible_s34
    ) {
      
      return(
        data.frame(
          n =
            504L,
          n_unique_immune =
            length(
              unique(
                immune_s34
              )
            ),
          spearman_rho =
            NA_real_,
          p_value =
            NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
    
    
    test_s34 <-
      suppressWarnings(
        cor.test(
          immune_s34,
          mpcds_z_s34,
          method =
            "spearman",
          exact =
            FALSE
        )
      )
    
    
    data.frame(
      n =
        504L,
      n_unique_immune =
        length(
          unique(
            immune_s34
          )
        ),
      spearman_rho =
        unname(
          test_s34$estimate
        ),
      p_value =
        test_s34$p.value,
      stringsAsFactors = FALSE
    )
  }


mpcds_split_s34 <-
  split(
    immune_mpcds_s34,
    interaction(
      immune_mpcds_s34[[
        "method"
      ]],
      immune_mpcds_s34[[
        "feature"
      ]],
      drop = TRUE,
      lex.order = TRUE
    )
  )


stopifnot(
  length(
    mpcds_split_s34
  ) ==
    114L
)


mpcds_results_s34 <-
  do.call(
    rbind,
    lapply(
      mpcds_split_s34,
      function(x_s34) {
        
        method_s34 <-
          unique(
            x_s34[[
              "method"
            ]]
          )
        
        
        feature_s34 <-
          unique(
            x_s34[[
              "feature"
            ]]
          )
        
        
        qc_s34 <-
          feature_qc_s34[
            feature_qc_s34[[
              "method"
            ]] ==
              method_s34 &
              feature_qc_s34[[
                "feature"
              ]] ==
              feature_s34,
            ,
            drop = FALSE
          ]
        
        
        test_s34 <-
          run_mpcds_test_s34(
            x_s34,
            qc_s34[[
              "test_eligible"
            ]]
          )
        
        
        data.frame(
          method =
            method_s34,
          feature =
            feature_s34,
          test_eligible =
            qc_s34[[
              "test_eligible"
            ]],
          exclusion_reason =
            qc_s34[[
              "exclusion_reason"
            ]],
          test_s34,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )
  )


rownames(
  mpcds_results_s34
) <-
  NULL


mpcds_results_s34$FDR_BH <-
  NA_real_


for (
  method_s34 in unique(
    mpcds_results_s34[[
      "method"
    ]]
  )
) {
  
  index_s34 <-
    which(
      mpcds_results_s34[[
        "method"
      ]] ==
        method_s34 &
        mpcds_results_s34[[
          "test_eligible"
        ]] &
        !is.na(
          mpcds_results_s34[[
            "p_value"
          ]]
        )
    )
  
  
  mpcds_results_s34[[
    "FDR_BH"
  ]][
    index_s34
  ] <-
    p.adjust(
      mpcds_results_s34[[
        "p_value"
      ]][
        index_s34
      ],
      method =
        "BH"
    )
}


mpcds_results_s34$significant_FDR_0_05 <-
  !is.na(
    mpcds_results_s34[[
      "FDR_BH"
    ]]
  ) &
  mpcds_results_s34[[
    "FDR_BH"
  ]] <
  0.05


mpcds_results_s34$significant_FDR_0_10 <-
  !is.na(
    mpcds_results_s34[[
      "FDR_BH"
    ]]
  ) &
  mpcds_results_s34[[
    "FDR_BH"
  ]] <
  0.10


mpcds_results_s34$direction <-
  ifelse(
    !mpcds_results_s34[[
      "test_eligible"
    ]],
    "Not tested",
    ifelse(
      mpcds_results_s34[[
        "spearman_rho"
      ]] >
        0,
      "Positive with MPCDS",
      ifelse(
        mpcds_results_s34[[
          "spearman_rho"
        ]] <
          0,
        "Negative with MPCDS",
        "No monotonic association"
      )
    )
  )


mpcds_results_s34$absolute_rho <-
  abs(
    mpcds_results_s34[[
      "spearman_rho"
    ]]
  )


mpcds_method_summary_s34 <-
  do.call(
    rbind,
    lapply(
      c(
        "xCell",
        "quanTIseq",
        "TIMER",
        "EPIC",
        "CIBERSORT"
      ),
      function(method_s34) {
        
        x_s34 <-
          mpcds_results_s34[
            mpcds_results_s34[[
              "method"
            ]] ==
              method_s34,
            ,
            drop = FALSE
          ]
        
        
        data.frame(
          method =
            method_s34,
          n_features =
            nrow(
              x_s34
            ),
          n_tested =
            sum(
              x_s34[[
                "test_eligible"
              ]]
            ),
          FDR_lt_0_05 =
            sum(
              x_s34[[
                "significant_FDR_0_05"
              ]]
            ),
          FDR_lt_0_10 =
            sum(
              x_s34[[
                "significant_FDR_0_10"
              ]]
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )


cat(
  "\n========================================\n",
  "MPCDS ASSOCIATION SUMMARY\n",
  "========================================\n"
)


print(
  mpcds_method_summary_s34,
  row.names = FALSE
)


cat(
  "\nTotal MPCDS FDR < 0.05:\n"
)


print(
  sum(
    mpcds_results_s34[[
      "significant_FDR_0_05"
    ]]
  )
)


write.csv(
  mpcds_results_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_multimethod_MPCDS_immune_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  mpcds_results_s34[
    mpcds_results_s34[[
      "significant_FDR_0_05"
    ]],
    ,
    drop = FALSE
  ],
  file.path(
    results_dir_s34,
    "TCGA_LUAD_multimethod_MPCDS_immune_FDR05.csv"
  ),
  row.names = FALSE
)


write.csv(
  mpcds_method_summary_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_multimethod_MPCDS_method_summary.csv"
  ),
  row.names = FALSE
)


saveRDS(
  mpcds_results_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_multimethod_MPCDS_immune_associations.rds"
  )
)


# ============================================================
# SECTION 16
# Manual strict biological mapping
# ============================================================


map_row_s34 <-
  function(
    biological_population,
    method,
    feature,
    match_tier
  ) {
    
    data.frame(
      biological_population =
        biological_population,
      method =
        method,
      feature =
        feature,
      match_tier =
        match_tier,
      stringsAsFactors = FALSE
    )
  }


biological_map_s34 <-
  do.call(
    rbind,
    list(
      
      map_row_s34(
        "B cells",
        "xCell",
        "B-cells_xCell",
        "strict"
      ),
      map_row_s34(
        "B cells",
        "quanTIseq",
        "B_cells_quantiseq",
        "strict"
      ),
      map_row_s34(
        "B cells",
        "TIMER",
        "B_cell_TIMER",
        "strict"
      ),
      map_row_s34(
        "B cells",
        "EPIC",
        "Bcells_EPIC",
        "strict"
      ),
      
      map_row_s34(
        "Naive B cells",
        "xCell",
        "naive_B-cells_xCell",
        "strict"
      ),
      map_row_s34(
        "Naive B cells",
        "CIBERSORT",
        "B_cells_naive_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "Memory B cells",
        "xCell",
        "Memory_B-cells_xCell",
        "strict"
      ),
      map_row_s34(
        "Memory B cells",
        "CIBERSORT",
        "B_cells_memory_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "Plasma cells",
        "xCell",
        "Plasma_cells_xCell",
        "strict"
      ),
      map_row_s34(
        "Plasma cells",
        "CIBERSORT",
        "Plasma_cells_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "CD4 T cells",
        "xCell",
        "CD4+_T-cells_xCell",
        "strict"
      ),
      map_row_s34(
        "CD4 T cells",
        "quanTIseq",
        "T_cells_CD4_quantiseq",
        "strict"
      ),
      map_row_s34(
        "CD4 T cells",
        "TIMER",
        "T_cell_CD4_TIMER",
        "strict"
      ),
      map_row_s34(
        "CD4 T cells",
        "EPIC",
        "CD4_Tcells_EPIC",
        "strict"
      ),
      
      map_row_s34(
        "CD8 T cells",
        "xCell",
        "CD8+_T-cells_xCell",
        "strict"
      ),
      map_row_s34(
        "CD8 T cells",
        "quanTIseq",
        "T_cells_CD8_quantiseq",
        "strict"
      ),
      map_row_s34(
        "CD8 T cells",
        "TIMER",
        "T_cell_CD8_TIMER",
        "strict"
      ),
      map_row_s34(
        "CD8 T cells",
        "EPIC",
        "CD8_Tcells_EPIC",
        "strict"
      ),
      map_row_s34(
        "CD8 T cells",
        "CIBERSORT",
        "T_cells_CD8_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "Regulatory T cells",
        "xCell",
        "Tregs_xCell",
        "strict"
      ),
      map_row_s34(
        "Regulatory T cells",
        "quanTIseq",
        "Tregs_quantiseq",
        "strict"
      ),
      map_row_s34(
        "Regulatory T cells",
        "CIBERSORT",
        "T_cells_regulatory_(Tregs)_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "Gamma-delta T cells",
        "xCell",
        "Tgd_cells_xCell",
        "strict"
      ),
      map_row_s34(
        "Gamma-delta T cells",
        "CIBERSORT",
        "T_cells_gamma_delta_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "NK cells",
        "xCell",
        "NK_cells_xCell",
        "strict"
      ),
      map_row_s34(
        "NK cells",
        "quanTIseq",
        "NK_cells_quantiseq",
        "strict"
      ),
      map_row_s34(
        "NK cells",
        "EPIC",
        "NKcells_EPIC",
        "strict"
      ),
      
      map_row_s34(
        "Monocytes",
        "xCell",
        "Monocytes_xCell",
        "strict"
      ),
      map_row_s34(
        "Monocytes",
        "quanTIseq",
        "Monocytes_quantiseq",
        "strict"
      ),
      map_row_s34(
        "Monocytes",
        "CIBERSORT",
        "Monocytes_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "Macrophages",
        "xCell",
        "Macrophages_xCell",
        "strict"
      ),
      map_row_s34(
        "Macrophages",
        "TIMER",
        "Macrophage_TIMER",
        "strict"
      ),
      map_row_s34(
        "Macrophages",
        "EPIC",
        "Macrophages_EPIC",
        "strict"
      ),
      
      map_row_s34(
        "M1 macrophages",
        "xCell",
        "Macrophages_M1_xCell",
        "strict"
      ),
      map_row_s34(
        "M1 macrophages",
        "quanTIseq",
        "Macrophages_M1_quantiseq",
        "strict"
      ),
      map_row_s34(
        "M1 macrophages",
        "CIBERSORT",
        "Macrophages_M1_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "M2 macrophages",
        "xCell",
        "Macrophages_M2_xCell",
        "strict"
      ),
      map_row_s34(
        "M2 macrophages",
        "quanTIseq",
        "Macrophages_M2_quantiseq",
        "strict"
      ),
      map_row_s34(
        "M2 macrophages",
        "CIBERSORT",
        "Macrophages_M2_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "Dendritic cells",
        "xCell",
        "DC_xCell",
        "strict"
      ),
      map_row_s34(
        "Dendritic cells",
        "quanTIseq",
        "Dendritic_cells_quantiseq",
        "strict"
      ),
      map_row_s34(
        "Dendritic cells",
        "TIMER",
        "DC_TIMER",
        "strict"
      ),
      
      map_row_s34(
        "Neutrophils",
        "xCell",
        "Neutrophils_xCell",
        "strict"
      ),
      map_row_s34(
        "Neutrophils",
        "quanTIseq",
        "Neutrophils_quantiseq",
        "strict"
      ),
      map_row_s34(
        "Neutrophils",
        "TIMER",
        "Neutrophil_TIMER",
        "strict"
      ),
      map_row_s34(
        "Neutrophils",
        "CIBERSORT",
        "Neutrophils_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "Eosinophils",
        "xCell",
        "Eosinophils_xCell",
        "strict"
      ),
      map_row_s34(
        "Eosinophils",
        "CIBERSORT",
        "Eosinophils_CIBERSORT",
        "strict"
      ),
      
      map_row_s34(
        "Endothelial cells",
        "xCell",
        "Endothelial_cells_xCell",
        "strict"
      ),
      map_row_s34(
        "Endothelial cells",
        "EPIC",
        "Endothelial_EPIC",
        "strict"
      ),
      
      map_row_s34(
        "Fibroblast/CAF axis",
        "xCell",
        "Fibroblasts_xCell",
        "approximate"
      ),
      map_row_s34(
        "Fibroblast/CAF axis",
        "EPIC",
        "CAFs_EPIC",
        "approximate"
      )
    )
  )


rownames(
  biological_map_s34
) <-
  NULL


inventory_key_s34 <-
  paste(
    feature_inventory_s34[[
      "method"
    ]],
    feature_inventory_s34[[
      "feature"
    ]],
    sep =
      "|||"
  )


map_key_s34 <-
  paste(
    biological_map_s34[[
      "method"
    ]],
    biological_map_s34[[
      "feature"
    ]],
    sep =
      "|||"
  )


stopifnot(
  all(
    map_key_s34 %in%
      inventory_key_s34
  )
)


stopifnot(
  !any(
    duplicated(
      map_key_s34
    )
  )
)


stopifnot(
  nrow(
    biological_map_s34
  ) ==
    52L
)


stopifnot(
  sum(
    biological_map_s34[[
      "match_tier"
    ]] ==
      "strict"
  ) ==
    50L
)


stopifnot(
  sum(
    biological_map_s34[[
      "match_tier"
    ]] ==
      "approximate"
  ) ==
    2L
)


# ============================================================
# SECTION 17
# Attach subtype and MPCDS statistics to biological map
# ============================================================


crossmethod_mapping_s34 <-
  merge(
    biological_map_s34,
    subtype_results_s34[
      ,
      c(
        "method",
        "feature",
        "test_eligible",
        "rank_biserial_C2_vs_C1",
        "FDR_BH",
        "significant_FDR_0_05",
        "direction"
      ),
      drop = FALSE
    ],
    by =
      c(
        "method",
        "feature"
      ),
    all.x =
      TRUE,
    sort =
      FALSE
  )


names(
  crossmethod_mapping_s34
)[
  names(
    crossmethod_mapping_s34
  ) ==
    "FDR_BH"
] <-
  "subtype_FDR_BH"


names(
  crossmethod_mapping_s34
)[
  names(
    crossmethod_mapping_s34
  ) ==
    "significant_FDR_0_05"
] <-
  "subtype_FDR05"


names(
  crossmethod_mapping_s34
)[
  names(
    crossmethod_mapping_s34
  ) ==
    "direction"
] <-
  "subtype_direction"


crossmethod_mapping_s34 <-
  merge(
    crossmethod_mapping_s34,
    mpcds_results_s34[
      ,
      c(
        "method",
        "feature",
        "spearman_rho",
        "FDR_BH",
        "significant_FDR_0_05",
        "direction"
      ),
      drop = FALSE
    ],
    by =
      c(
        "method",
        "feature"
      ),
    all.x =
      TRUE,
    sort =
      FALSE
  )


names(
  crossmethod_mapping_s34
)[
  names(
    crossmethod_mapping_s34
  ) ==
    "FDR_BH"
] <-
  "MPCDS_FDR_BH"


names(
  crossmethod_mapping_s34
)[
  names(
    crossmethod_mapping_s34
  ) ==
    "significant_FDR_0_05"
] <-
  "MPCDS_FDR05"


names(
  crossmethod_mapping_s34
)[
  names(
    crossmethod_mapping_s34
  ) ==
    "direction"
] <-
  "MPCDS_direction"


stopifnot(
  nrow(
    crossmethod_mapping_s34
  ) ==
    52L
)


write.csv(
  biological_map_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_manual_biological_mapping.csv"
  ),
  row.names = FALSE
)


write.csv(
  crossmethod_mapping_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_mapping_with_associations.csv"
  ),
  row.names = FALSE
)


saveRDS(
  crossmethod_mapping_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_crossmethod_mapping_with_associations.rds"
  )
)


# ============================================================
# SECTION 18
# Cross-method replication analysis
#
# Replication rule:
#
# - at least two FDR < 0.05 method-level results;
# - at least two significant methods in the same direction;
# - no FDR-significant method in the opposite direction.
#
# Nonsignificant methods are not treated as contradictions.
# ============================================================


strict_mapping_s34 <-
  crossmethod_mapping_s34[
    crossmethod_mapping_s34[[
      "match_tier"
    ]] ==
      "strict",
    ,
    drop = FALSE
  ]


approximate_mapping_s34 <-
  crossmethod_mapping_s34[
    crossmethod_mapping_s34[[
      "match_tier"
    ]] ==
      "approximate",
    ,
    drop = FALSE
  ]


strict_mapping_s34$subtype_sign <-
  sign(
    strict_mapping_s34[[
      "rank_biserial_C2_vs_C1"
    ]]
  )


strict_mapping_s34$MPCDS_sign <-
  sign(
    strict_mapping_s34[[
      "spearman_rho"
    ]]
  )


approximate_mapping_s34$subtype_sign <-
  sign(
    approximate_mapping_s34[[
      "rank_biserial_C2_vs_C1"
    ]]
  )


approximate_mapping_s34$MPCDS_sign <-
  sign(
    approximate_mapping_s34[[
      "spearman_rho"
    ]]
  )


summarize_replication_s34 <-
  function(
    x_s34
  ) {
    
    subtype_sig_s34 <-
      x_s34[[
        "subtype_FDR05"
      ]]
    
    
    subtype_C1_s34 <-
      sum(
        subtype_sig_s34 &
          x_s34[[
            "subtype_sign"
          ]] <
          0,
        na.rm = TRUE
      )
    
    
    subtype_C2_s34 <-
      sum(
        subtype_sig_s34 &
          x_s34[[
            "subtype_sign"
          ]] >
          0,
        na.rm = TRUE
      )
    
    
    subtype_conflict_s34 <-
      subtype_C1_s34 >
      0L &
      subtype_C2_s34 >
      0L
    
    
    subtype_direction_s34 <-
      if (
        subtype_C1_s34 >=
        2L &
        subtype_C2_s34 ==
        0L
      ) {
        
        "Higher in PCD_C1"
        
      } else if (
        subtype_C2_s34 >=
        2L &
        subtype_C1_s34 ==
        0L
      ) {
        
        "Higher in PCD_C2"
        
      } else {
        
        NA_character_
      }
    
    
    mpcds_sig_s34 <-
      x_s34[[
        "MPCDS_FDR05"
      ]]
    
    
    mpcds_negative_s34 <-
      sum(
        mpcds_sig_s34 &
          x_s34[[
            "MPCDS_sign"
          ]] <
          0,
        na.rm = TRUE
      )
    
    
    mpcds_positive_s34 <-
      sum(
        mpcds_sig_s34 &
          x_s34[[
            "MPCDS_sign"
          ]] >
          0,
        na.rm = TRUE
      )
    
    
    mpcds_conflict_s34 <-
      mpcds_negative_s34 >
      0L &
      mpcds_positive_s34 >
      0L
    
    
    mpcds_direction_s34 <-
      if (
        mpcds_negative_s34 >=
        2L &
        mpcds_positive_s34 ==
        0L
      ) {
        
        "Negative with MPCDS"
        
      } else if (
        mpcds_positive_s34 >=
        2L &
        mpcds_negative_s34 ==
        0L
      ) {
        
        "Positive with MPCDS"
        
      } else {
        
        NA_character_
      }
    
    
    subtype_all_signs_s34 <-
      unique(
        x_s34[[
          "subtype_sign"
        ]][
          !is.na(
            x_s34[[
              "subtype_sign"
            ]]
          ) &
            x_s34[[
              "subtype_sign"
            ]] !=
            0
        ]
      )
    
    
    mpcds_all_signs_s34 <-
      unique(
        x_s34[[
          "MPCDS_sign"
        ]][
          !is.na(
            x_s34[[
              "MPCDS_sign"
            ]]
          ) &
            x_s34[[
              "MPCDS_sign"
            ]] !=
            0
        ]
      )
    
    
    data.frame(
      n_methods_mapped =
        length(
          unique(
            x_s34[[
              "method"
            ]]
          )
        ),
      
      methods =
        paste(
          x_s34[[
            "method"
          ]],
          collapse =
            "; "
        ),
      
      subtype_n_FDR05 =
        sum(
          subtype_sig_s34,
          na.rm = TRUE
        ),
      
      subtype_n_FDR05_higher_C1 =
        subtype_C1_s34,
      
      subtype_n_FDR05_higher_C2 =
        subtype_C2_s34,
      
      subtype_significant_conflict =
        subtype_conflict_s34,
      
      subtype_replicated_FDR05 =
        !is.na(
          subtype_direction_s34
        ),
      
      subtype_replication_direction =
        subtype_direction_s34,
      
      subtype_all_methods_direction_concordant =
        length(
          subtype_all_signs_s34
        ) <=
        1L,
      
      MPCDS_n_FDR05 =
        sum(
          mpcds_sig_s34,
          na.rm = TRUE
        ),
      
      MPCDS_n_FDR05_negative =
        mpcds_negative_s34,
      
      MPCDS_n_FDR05_positive =
        mpcds_positive_s34,
      
      MPCDS_significant_conflict =
        mpcds_conflict_s34,
      
      MPCDS_replicated_FDR05 =
        !is.na(
          mpcds_direction_s34
        ),
      
      MPCDS_replication_direction =
        mpcds_direction_s34,
      
      MPCDS_all_methods_direction_concordant =
        length(
          mpcds_all_signs_s34
        ) <=
        1L,
      
      stringsAsFactors = FALSE
    )
  }


strict_population_order_s34 <-
  unique(
    biological_map_s34[[
      "biological_population"
    ]][
      biological_map_s34[[
        "match_tier"
      ]] ==
        "strict"
    ]
  )


strict_replication_s34 <-
  do.call(
    rbind,
    lapply(
      strict_population_order_s34,
      function(population_s34) {
        
        x_s34 <-
          strict_mapping_s34[
            strict_mapping_s34[[
              "biological_population"
            ]] ==
              population_s34,
            ,
            drop = FALSE
          ]
        
        
        cbind(
          data.frame(
            biological_population =
              population_s34,
            match_tier =
              "strict",
            stringsAsFactors = FALSE
          ),
          summarize_replication_s34(
            x_s34
          )
        )
      }
    )
  )


rownames(
  strict_replication_s34
) <-
  NULL


approximate_replication_s34 <-
  do.call(
    rbind,
    lapply(
      unique(
        approximate_mapping_s34[[
          "biological_population"
        ]]
      ),
      function(population_s34) {
        
        x_s34 <-
          approximate_mapping_s34[
            approximate_mapping_s34[[
              "biological_population"
            ]] ==
              population_s34,
            ,
            drop = FALSE
          ]
        
        
        cbind(
          data.frame(
            biological_population =
              population_s34,
            match_tier =
              "approximate",
            stringsAsFactors = FALSE
          ),
          summarize_replication_s34(
            x_s34
          )
        )
      }
    )
  )


stopifnot(
  nrow(
    strict_replication_s34
  ) ==
    17L
)


stopifnot(
  nrow(
    approximate_replication_s34
  ) ==
    1L
)


subtype_replicated_s34 <-
  strict_replication_s34[
    strict_replication_s34[[
      "subtype_replicated_FDR05"
    ]],
    ,
    drop = FALSE
  ]


mpcds_replicated_s34 <-
  strict_replication_s34[
    strict_replication_s34[[
      "MPCDS_replicated_FDR05"
    ]],
    ,
    drop = FALSE
  ]


cat(
  "\n========================================\n",
  "STRICT CROSS-METHOD REPLICATION\n",
  "========================================\n"
)


cat(
  "\nSubtype-replicated populations:\n"
)


print(
  subtype_replicated_s34[
    ,
    c(
      "biological_population",
      "subtype_n_FDR05",
      "subtype_replication_direction"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


cat(
  "\nMPCDS-replicated populations:\n"
)


print(
  mpcds_replicated_s34[
    ,
    c(
      "biological_population",
      "MPCDS_n_FDR05",
      "MPCDS_replication_direction"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


cat(
  "\nSubtype conflicts:\n"
)


print(
  strict_replication_s34[
    strict_replication_s34[[
      "subtype_significant_conflict"
    ]],
    c(
      "biological_population",
      "subtype_n_FDR05_higher_C1",
      "subtype_n_FDR05_higher_C2"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


cat(
  "\nMPCDS conflicts:\n"
)


print(
  strict_replication_s34[
    strict_replication_s34[[
      "MPCDS_significant_conflict"
    ]],
    c(
      "biological_population",
      "MPCDS_n_FDR05_negative",
      "MPCDS_n_FDR05_positive"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


write.csv(
  strict_replication_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_strict_replication_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  subtype_replicated_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_subtype_replicated_FDR05.csv"
  ),
  row.names = FALSE
)


write.csv(
  mpcds_replicated_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_MPCDS_replicated_FDR05.csv"
  ),
  row.names = FALSE
)


write.csv(
  approximate_replication_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_approximate_mapping_sensitivity.csv"
  ),
  row.names = FALSE
)


saveRDS(
  strict_replication_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_crossmethod_strict_replication_summary.rds"
  )
)


# ============================================================
# SECTION 19
# MCP-counter supportive mapping
#
# MCP-counter results are an existing LUAD/LUSC cross-cohort
# robustness layer and are not counted as one of the five new
# LUAD deconvolution methods.
# ============================================================


mcp_biological_map_s34 <-
  data.frame(
    MCPcounter_population =
      c(
        "CD8 T cells",
        "Endothelial cells",
        "B lineage",
        "Myeloid dendritic cells",
        "Monocytic lineage",
        "Fibroblasts",
        "T cells"
      ),
    biological_population =
      c(
        "CD8 T cells",
        "Endothelial cells",
        "B cells",
        "Dendritic cells",
        "Monocytes",
        "Fibroblast/CAF axis",
        NA_character_
      ),
    MCP_match_tier =
      c(
        "strict",
        "strict",
        "approximate",
        "approximate",
        "approximate",
        "approximate",
        "unmapped_broad"
      ),
    stringsAsFactors = FALSE
  )


stopifnot(
  setequal(
    mcp_biological_map_s34[[
      "MCPcounter_population"
    ]],
    mcp_existing_s34[[
      "cell_population"
    ]]
  )
)


mcp_mapping_audit_s34 <-
  data.frame(
    MCPcounter_population =
      mcp_existing_s34[[
        "cell_population"
      ]],
    LUAD_C1_mean =
      mcp_existing_s34[[
        "LUAD_C1_mean"
      ]],
    LUAD_C2_mean =
      mcp_existing_s34[[
        "LUAD_C2_mean"
      ]],
    LUAD_difference_C2_minus_C1 =
      mcp_existing_s34[[
        "LUAD_difference"
      ]],
    LUAD_FDR =
      mcp_existing_s34[[
        "LUAD_padj"
      ]],
    LUAD_significant =
      mcp_existing_s34[[
        "LUAD_significant"
      ]],
    robust_cross_cohort =
      mcp_existing_s34[[
        "robust_cross_cohort"
      ]],
    LUAD_direction =
      ifelse(
        mcp_existing_s34[[
          "LUAD_difference"
        ]] >
          0,
        "Higher in PCD_C2",
        ifelse(
          mcp_existing_s34[[
            "LUAD_difference"
          ]] <
            0,
          "Higher in PCD_C1",
          "No difference"
        )
      ),
    stringsAsFactors = FALSE
  )


mcp_support_s34 <-
  merge(
    mcp_biological_map_s34,
    mcp_mapping_audit_s34,
    by =
      "MCPcounter_population",
    all.x =
      TRUE,
    sort =
      FALSE
  )


mcp_mappable_s34 <-
  mcp_support_s34[
    !is.na(
      mcp_support_s34[[
        "biological_population"
      ]]
    ),
    ,
    drop = FALSE
  ]


strict_with_mcp_s34 <-
  merge(
    strict_replication_s34,
    mcp_mappable_s34[
      ,
      c(
        "biological_population",
        "MCPcounter_population",
        "MCP_match_tier",
        "LUAD_difference_C2_minus_C1",
        "LUAD_FDR",
        "LUAD_significant",
        "robust_cross_cohort",
        "LUAD_direction"
      ),
      drop = FALSE
    ],
    by =
      "biological_population",
    all.x =
      TRUE,
    sort =
      FALSE
  )


write.csv(
  mcp_mapping_audit_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_MCPcounter_crossmethod_mapping_audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  mcp_support_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_MCPcounter_support.csv"
  ),
  row.names = FALSE
)


write.csv(
  strict_with_mcp_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_strict_replication_with_MCPcounter.csv"
  ),
  row.names = FALSE
)


# ============================================================
# SECTION 20
# Final replication summary table
# ============================================================


final_summary_s34 <-
  strict_replication_s34[
    ,
    c(
      "biological_population",
      "n_methods_mapped",
      "subtype_n_FDR05",
      "subtype_replicated_FDR05",
      "subtype_replication_direction",
      "subtype_significant_conflict",
      "MPCDS_n_FDR05",
      "MPCDS_replicated_FDR05",
      "MPCDS_replication_direction",
      "MPCDS_significant_conflict"
    ),
    drop = FALSE
  ]


final_summary_s34$subtype_status <-
  ifelse(
    final_summary_s34[[
      "subtype_significant_conflict"
    ]],
    "Conflicting significant directions",
    ifelse(
      final_summary_s34[[
        "subtype_replicated_FDR05"
      ]],
      final_summary_s34[[
        "subtype_replication_direction"
      ]],
      "Not replicated"
    )
  )


final_summary_s34$MPCDS_status <-
  ifelse(
    final_summary_s34[[
      "MPCDS_significant_conflict"
    ]],
    "Conflicting significant directions",
    ifelse(
      final_summary_s34[[
        "MPCDS_replicated_FDR05"
      ]],
      final_summary_s34[[
        "MPCDS_replication_direction"
      ]],
      "Not replicated"
    )
  )


mcp_summary_s34 <-
  strict_with_mcp_s34[
    ,
    c(
      "biological_population",
      "MCPcounter_population",
      "MCP_match_tier",
      "LUAD_FDR",
      "LUAD_direction",
      "robust_cross_cohort"
    ),
    drop = FALSE
  ]


final_summary_s34 <-
  merge(
    final_summary_s34,
    mcp_summary_s34,
    by =
      "biological_population",
    all.x =
      TRUE,
    sort =
      FALSE
  )


# ============================================================
# SECTION 21
# Figure generation
# ============================================================

# Final manuscript figures are intentionally deferred to Script 32.
# Script 29 produces the validated multi-method immune-deconvolution
# tables and replication summaries only.


# ============================================================
# SECTION 22
# Save final summary outputs
# ============================================================


write.csv(
  final_summary_s34,
  file.path(
    results_dir_s34,
    "TCGA_LUAD_crossmethod_final_immune_summary.csv"
  ),
  row.names = FALSE
)


saveRDS(
  final_summary_s34,
  file.path(
    processed_dir_s34,
    "TCGA_LUAD_crossmethod_final_immune_summary.rds"
  )
)


# ============================================================
# SECTION 23
# Final expected-result integrity audit
#
# These counts come from the completed section-by-section run.
# If a clean source produces different values, the script
# should stop so the difference can be investigated.
# ============================================================


final_audit_s34 <-
  data.frame(
    metric =
      c(
        "patients",
        "immune_features_total",
        "immune_features_tested",
        "subtype_FDR05_features",
        "MPCDS_FDR05_features",
        "strict_populations",
        "subtype_replicated_populations",
        "subtype_conflicted_populations",
        "MPCDS_replicated_populations",
        "MPCDS_conflicted_populations"
      ),
    observed =
      c(
        517L,
        nrow(
          feature_qc_s34
        ),
        sum(
          feature_qc_s34[[
            "test_eligible"
          ]]
        ),
        sum(
          subtype_results_s34[[
            "significant_FDR_0_05"
          ]]
        ),
        sum(
          mpcds_results_s34[[
            "significant_FDR_0_05"
          ]]
        ),
        nrow(
          strict_replication_s34
        ),
        sum(
          strict_replication_s34[[
            "subtype_replicated_FDR05"
          ]]
        ),
        sum(
          strict_replication_s34[[
            "subtype_significant_conflict"
          ]]
        ),
        sum(
          strict_replication_s34[[
            "MPCDS_replicated_FDR05"
          ]]
        ),
        sum(
          strict_replication_s34[[
            "MPCDS_significant_conflict"
          ]]
        )
      ),
    expected =
      c(
        517L,
        114L,
        113L,
        45L,
        38L,
        17L,
        6L,
        3L,
        4L,
        1L
      ),
    stringsAsFactors = FALSE
  )


final_audit_s34$matches_expected <-
  final_audit_s34[[
    "observed"
  ]] ==
  final_audit_s34[[
    "expected"
  ]]


cat(
  "\n========================================\n",
  "FINAL SCRIPT 29 AUDIT\n",
  "========================================\n"
)


print(
  final_audit_s34,
  row.names = FALSE
)


stopifnot(
  all(
    final_audit_s34[[
      "matches_expected"
    ]]
  )
)

cat(
  "\nEstablished immune-deconvolution regression checks: PASSED\n"
)


# ============================================================
# SECTION 24
# Final output-file verification
# ============================================================


final_files_s34 <-
  c(
    file.path(
      results_dir_s34,
      "TCGA_LUAD_xCell_full_scores.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_quanTIseq_full_scores.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_TIMER_full_scores.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_EPIC_full_scores.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_CIBERSORT_LM22_full_scores.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_multimethod_immune_feature_inventory.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_multimethod_immune_feature_QC.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_multimethod_PCD_subtype_immune_associations.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_multimethod_MPCDS_immune_associations.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_crossmethod_mapping_with_associations.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_crossmethod_strict_replication_summary.csv"
    ),
    file.path(
      results_dir_s34,
      "TCGA_LUAD_crossmethod_final_immune_summary.csv"
    )
  )


final_file_check_s34 <-
  data.frame(
    file =
      final_files_s34,
    exists =
      file.exists(
        final_files_s34
      ),
    size_bytes =
      file.info(
        final_files_s34
      )$size,
    stringsAsFactors = FALSE
  )


cat(
  "\n========================================\n",
  "FINAL OUTPUT FILE CHECK\n",
  "========================================\n"
)


print(
  final_file_check_s34,
  row.names = FALSE
)


stopifnot(
  all(
    final_file_check_s34[[
      "exists"
    ]]
  )
)


stopifnot(
  all(
    final_file_check_s34[[
      "size_bytes"
    ]] >
      0
  )
)


# ============================================================
# SECTION 25
# Completion marker
# ============================================================


cat(
  "\n========================================\n",
  "SCRIPT 29 IMMUNE-DECONVOLUTION SENSITIVITY ANALYSIS COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)