# ============================================================
# SCRIPT 28
# TCGA-LUAD 70/30 train-test + LASSO-Cox sensitivity analysis
#
# PURPOSE
# Proposal-aligned methodological sensitivity analysis.
#
# IMPORTANT:
# - This does NOT replace the primary repeated nested-CV ridge MPCDS.
# - Candidate genes are fixed before modelling.
# - The train/test split is fixed before fitting.
# - Scaling is derived from training data only.
# - Lambda selection uses training-only cross-validation.
# - The held-out test cohort is used only after the sparse model
#   has been locked.
# ============================================================


source("03_scripts/00_project_config.R")


cat(
  "\n========================================\n",
  "SCRIPT 28\n",
  "TCGA-LUAD 70/30 LASSO-COX SENSITIVITY\n",
  "========================================\n"
)


# ============================================================
# SECTION 1
# Environment, inputs, packages, and data inventory
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 1: INPUT INVENTORY\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 1.1 Required packages
# ------------------------------------------------------------


required_packages_s33 <-
  c(
    "survival",
    "glmnet",
    "timeROC"
  )


package_status_s33 <-
  vapply(
    required_packages_s33,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )


cat(
  "\nPackage status:\n"
)


print(
  package_status_s33
)


if (
  !all(
    package_status_s33
  )
) {
  
  stop(
    paste0(
      "Missing required package(s): ",
      paste(
        names(
          package_status_s33
        )[
          !package_status_s33
        ],
        collapse = ", "
      )
    )
  )
}


cat(
  "\nPackage versions:\n"
)


print(
  data.frame(
    package =
      required_packages_s33,
    version =
      vapply(
        required_packages_s33,
        function(pkg_s33) {
          as.character(
            packageVersion(
              pkg_s33
            )
          )
        },
        character(1)
      ),
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 1.2 Input files
# ------------------------------------------------------------


mpcds_file_s33 <-
  "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv"


weight_file_s33 <-
  file.path(
    "04_results/mpcds/external_validation",
    "TCGA_LUAD_MPCDS_standardized_weights_47.csv"
  )


patient_tpm_file_s33 <-
  "02_processed_data/immunotherapy/TCGA_LUAD_patient_protein_coding_TPM.rds"


input_files_s33 <-
  c(
    mpcds_file_s33,
    weight_file_s33,
    patient_tpm_file_s33
  )


cat(
  "\nInput file check:\n"
)


print(
  data.frame(
    file =
      input_files_s33,
    exists =
      file.exists(
        input_files_s33
      ),
    size_bytes =
      file.info(
        input_files_s33
      )$size,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


stopifnot(
  all(
    file.exists(
      input_files_s33
    )
  )
)


# ------------------------------------------------------------
# 1.3 Output directories
# ------------------------------------------------------------


results_dir_s33 <-
  "04_results/mpcds/lasso_70_30_sensitivity"


processed_dir_s33 <-
  "02_processed_data/survival/lasso_70_30_sensitivity"



dir.create(
  results_dir_s33,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  processed_dir_s33,
  recursive = TRUE,
  showWarnings = FALSE
)



# ------------------------------------------------------------
# 1.4 Load inputs
# ------------------------------------------------------------


mpcds_s33 <-
  read.csv(
    mpcds_file_s33,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


weights_s33 <-
  read.csv(
    weight_file_s33,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


patient_tpm_s33 <-
  as.matrix(
    readRDS(
      patient_tpm_file_s33
    )
  )


# ------------------------------------------------------------
# 1.5 Survival cohort validation
# ------------------------------------------------------------


stopifnot(
  all(
    c(
      "patient_id",
      "OS_time",
      "OS_status"
    ) %in%
      colnames(
        mpcds_s33
      )
  )
)


stopifnot(
  nrow(
    mpcds_s33
  ) ==
    504
)


stopifnot(
  sum(
    mpcds_s33[[
      "OS_status"
    ]] ==
      1
  ) ==
    182
)


stopifnot(
  !any(
    duplicated(
      mpcds_s33[[
        "patient_id"
      ]]
    )
  )
)


stopifnot(
  all(
    mpcds_s33[[
      "OS_status"
    ]] %in%
      c(
        0,
        1
      )
  )
)


stopifnot(
  all(
    is.finite(
      mpcds_s33[[
        "OS_time"
      ]]
    )
  )
)


stopifnot(
  all(
    mpcds_s33[[
      "OS_time"
    ]] >
      0
  )
)


# ------------------------------------------------------------
# 1.6 Candidate gene validation
# ------------------------------------------------------------


stopifnot(
  "gene" %in%
    colnames(
      weights_s33
    )
)


candidate_genes_s33 <-
  unique(
    as.character(
      weights_s33[[
        "gene"
      ]]
    )
  )


stopifnot(
  length(
    candidate_genes_s33
  ) ==
    47
)


missing_genes_s33 <-
  setdiff(
    candidate_genes_s33,
    rownames(
      patient_tpm_s33
    )
  )


missing_patients_s33 <-
  setdiff(
    mpcds_s33[[
      "patient_id"
    ]],
    colnames(
      patient_tpm_s33
    )
  )


cat(
  "\nCandidate genes:\n"
)


print(
  length(
    candidate_genes_s33
  )
)


cat(
  "\nMissing candidate genes:\n"
)


print(
  length(
    missing_genes_s33
  )
)


cat(
  "\nSurvival patients missing from expression:\n"
)


print(
  length(
    missing_patients_s33
  )
)


stopifnot(
  length(
    missing_genes_s33
  ) ==
    0
)


stopifnot(
  length(
    missing_patients_s33
  ) ==
    0
)


cat(
  "\n========================================\n",
  "SECTION 1 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 2
# Build locked modelling matrix and 70/30 split
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 2: MODELLING MATRIX + 70/30 SPLIT\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 2.1 Construct patient x gene matrix
# ------------------------------------------------------------


expression_47_tpm_s33 <-
  t(
    patient_tpm_s33[
      candidate_genes_s33,
      mpcds_s33[[
        "patient_id"
      ]],
      drop = FALSE
    ]
  )


stopifnot(
  identical(
    dim(
      expression_47_tpm_s33
    ),
    c(
      504L,
      47L
    )
  )
)


stopifnot(
  identical(
    rownames(
      expression_47_tpm_s33
    ),
    mpcds_s33[[
      "patient_id"
    ]]
  )
)


stopifnot(
  identical(
    colnames(
      expression_47_tpm_s33
    ),
    candidate_genes_s33
  )
)


stopifnot(
  !anyNA(
    expression_47_tpm_s33
  )
)


stopifnot(
  all(
    is.finite(
      expression_47_tpm_s33
    )
  )
)


# ------------------------------------------------------------
# 2.2 log2(TPM + 1)
# ------------------------------------------------------------


expression_47_log_s33 <-
  log2(
    expression_47_tpm_s33 +
      1
  )


full_gene_sd_s33 <-
  apply(
    expression_47_log_s33,
    2,
    sd
  )


stopifnot(
  all(
    is.finite(
      full_gene_sd_s33
    )
  )
)


stopifnot(
  all(
    full_gene_sd_s33 >
      0
  )
)


# ------------------------------------------------------------
# 2.3 Reproducible event-stratified 70/30 split
# ------------------------------------------------------------


split_seed_s33 <-
  33001L


set.seed(
  split_seed_s33
)


event_indices_s33 <-
  which(
    mpcds_s33[[
      "OS_status"
    ]] ==
      1
  )


censored_indices_s33 <-
  which(
    mpcds_s33[[
      "OS_status"
    ]] ==
      0
  )


n_train_events_s33 <-
  floor(
    0.70 *
      length(
        event_indices_s33
      )
  )


n_train_censored_s33 <-
  floor(
    0.70 *
      length(
        censored_indices_s33
      )
  )


train_event_indices_s33 <-
  sample(
    event_indices_s33,
    n_train_events_s33,
    replace = FALSE
  )


train_censored_indices_s33 <-
  sample(
    censored_indices_s33,
    n_train_censored_s33,
    replace = FALSE
  )


train_indices_s33 <-
  sort(
    c(
      train_event_indices_s33,
      train_censored_indices_s33
    )
  )


test_indices_s33 <-
  setdiff(
    seq_len(
      nrow(
        mpcds_s33
      )
    ),
    train_indices_s33
  )


stopifnot(
  length(
    intersect(
      train_indices_s33,
      test_indices_s33
    )
  ) ==
    0
)


stopifnot(
  length(
    union(
      train_indices_s33,
      test_indices_s33
    )
  ) ==
    504
)


train_metadata_s33 <-
  mpcds_s33[
    train_indices_s33,
    c(
      "patient_id",
      "OS_time",
      "OS_status"
    ),
    drop = FALSE
  ]


test_metadata_s33 <-
  mpcds_s33[
    test_indices_s33,
    c(
      "patient_id",
      "OS_time",
      "OS_status"
    ),
    drop = FALSE
  ]


split_summary_s33 <-
  data.frame(
    cohort =
      c(
        "Full",
        "Training",
        "Test"
      ),
    n =
      c(
        504,
        nrow(
          train_metadata_s33
        ),
        nrow(
          test_metadata_s33
        )
      ),
    events =
      c(
        182,
        sum(
          train_metadata_s33[[
            "OS_status"
          ]]
        ),
        sum(
          test_metadata_s33[[
            "OS_status"
          ]]
        )
      ),
    censored =
      c(
        322,
        sum(
          train_metadata_s33[[
            "OS_status"
          ]] ==
            0
        ),
        sum(
          test_metadata_s33[[
            "OS_status"
          ]] ==
            0
        )
      ),
    event_fraction =
      c(
        mean(
          mpcds_s33[[
            "OS_status"
          ]]
        ),
        mean(
          train_metadata_s33[[
            "OS_status"
          ]]
        ),
        mean(
          test_metadata_s33[[
            "OS_status"
          ]]
        )
      ),
    stringsAsFactors = FALSE
  )


cat(
  "\n70/30 split summary:\n"
)


print(
  split_summary_s33,
  row.names = FALSE
)


stopifnot(
  nrow(
    train_metadata_s33
  ) ==
    352
)


stopifnot(
  nrow(
    test_metadata_s33
  ) ==
    152
)


stopifnot(
  sum(
    train_metadata_s33[[
      "OS_status"
    ]]
  ) ==
    127
)


stopifnot(
  sum(
    test_metadata_s33[[
      "OS_status"
    ]]
  ) ==
    55
)


# ------------------------------------------------------------
# 2.4 Split expression
# ------------------------------------------------------------


x_train_raw_s33 <-
  expression_47_log_s33[
    train_indices_s33,
    ,
    drop = FALSE
  ]


x_test_raw_s33 <-
  expression_47_log_s33[
    test_indices_s33,
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 2.5 Training-only scaling
# ------------------------------------------------------------


train_gene_mean_s33 <-
  colMeans(
    x_train_raw_s33
  )


train_gene_sd_s33 <-
  apply(
    x_train_raw_s33,
    2,
    sd
  )


stopifnot(
  all(
    is.finite(
      train_gene_sd_s33
    )
  )
)


stopifnot(
  all(
    train_gene_sd_s33 >
      0
  )
)


x_train_z_s33 <-
  sweep(
    x_train_raw_s33,
    2,
    train_gene_mean_s33,
    FUN = "-"
  )


x_train_z_s33 <-
  sweep(
    x_train_z_s33,
    2,
    train_gene_sd_s33,
    FUN = "/"
  )


x_test_z_s33 <-
  sweep(
    x_test_raw_s33,
    2,
    train_gene_mean_s33,
    FUN = "-"
  )


x_test_z_s33 <-
  sweep(
    x_test_z_s33,
    2,
    train_gene_sd_s33,
    FUN = "/"
  )


stopifnot(
  max(
    abs(
      colMeans(
        x_train_z_s33
      )
    )
  ) <
    1e-10
)


stopifnot(
  max(
    abs(
      apply(
        x_train_z_s33,
        2,
        sd
      ) -
        1
    )
  ) <
    1e-10
)


# ------------------------------------------------------------
# 2.6 Save split and scaling information
# ------------------------------------------------------------


split_assignment_s33 <-
  data.frame(
    patient_id =
      mpcds_s33[[
        "patient_id"
      ]],
    OS_time =
      mpcds_s33[[
        "OS_time"
      ]],
    OS_status =
      mpcds_s33[[
        "OS_status"
      ]],
    split =
      ifelse(
        seq_len(
          nrow(
            mpcds_s33
          )
        ) %in%
          train_indices_s33,
        "Training",
        "Test"
      ),
    split_seed =
      split_seed_s33,
    stringsAsFactors = FALSE
  )


feature_table_s33 <-
  data.frame(
    gene =
      candidate_genes_s33,
    training_mean_log2TPM1 =
      train_gene_mean_s33[
        candidate_genes_s33
      ],
    training_SD_log2TPM1 =
      train_gene_sd_s33[
        candidate_genes_s33
      ],
    stringsAsFactors = FALSE
  )


write.csv(
  split_assignment_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_70_30_patient_split.csv"
  ),
  row.names = FALSE
)


write.csv(
  feature_table_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_47_gene_training_scaling_parameters.csv"
  ),
  row.names = FALSE
)


saveRDS(
  x_train_z_s33,
  file.path(
    processed_dir_s33,
    "TCGA_LUAD_LASSO_x_train_z_47.rds"
  )
)


saveRDS(
  x_test_z_s33,
  file.path(
    processed_dir_s33,
    "TCGA_LUAD_LASSO_x_test_z_47.rds"
  )
)


saveRDS(
  train_metadata_s33,
  file.path(
    processed_dir_s33,
    "TCGA_LUAD_LASSO_train_survival_metadata.rds"
  )
)


saveRDS(
  test_metadata_s33,
  file.path(
    processed_dir_s33,
    "TCGA_LUAD_LASSO_test_survival_metadata.rds"
  )
)


cat(
  "\n========================================\n",
  "SECTION 2 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 3
# Training-only LASSO-Cox tuning and pre-test model lock
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 3: TRAINING-ONLY LASSO-COX\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 3.1 Training survival object
# ------------------------------------------------------------


y_train_s33 <-
  survival::Surv(
    time =
      train_metadata_s33[[
        "OS_time"
      ]],
    event =
      train_metadata_s33[[
        "OS_status"
      ]]
  )


# ------------------------------------------------------------
# 3.2 Initial LASSO path
#
# Explicitly lock Breslow tie handling because glmnet 5.1
# plans to change the default to Efron.
# ------------------------------------------------------------


lasso_path_s33 <-
  glmnet::glmnet(
    x =
      x_train_z_s33,
    y =
      y_train_s33,
    family =
      "cox",
    alpha =
      1,
    standardize =
      FALSE,
    cox.ties =
      "breslow"
  )


# ------------------------------------------------------------
# 3.3 Reproducible 10-fold event-stratified CV
# ------------------------------------------------------------


cv_seed_s33 <-
  33002L


nfolds_s33 <-
  10L


set.seed(
  cv_seed_s33
)


event_train_indices_s33 <-
  which(
    train_metadata_s33[[
      "OS_status"
    ]] ==
      1
  )


censored_train_indices_s33 <-
  which(
    train_metadata_s33[[
      "OS_status"
    ]] ==
      0
  )


event_fold_labels_s33 <-
  sample(
    rep(
      seq_len(
        nfolds_s33
      ),
      length.out =
        length(
          event_train_indices_s33
        )
    )
  )


censored_fold_labels_s33 <-
  sample(
    rep(
      seq_len(
        nfolds_s33
      ),
      length.out =
        length(
          censored_train_indices_s33
        )
    )
  )


foldid_s33 <-
  integer(
    nrow(
      train_metadata_s33
    )
  )


foldid_s33[
  event_train_indices_s33
] <-
  event_fold_labels_s33


foldid_s33[
  censored_train_indices_s33
] <-
  censored_fold_labels_s33


fold_balance_s33 <-
  do.call(
    rbind,
    lapply(
      seq_len(
        nfolds_s33
      ),
      function(k_s33) {
        
        idx_s33 <-
          which(
            foldid_s33 ==
              k_s33
          )
        
        
        data.frame(
          fold =
            k_s33,
          n =
            length(
              idx_s33
            ),
          events =
            sum(
              train_metadata_s33[[
                "OS_status"
              ]][
                idx_s33
              ] ==
                1
            ),
          censored =
            sum(
              train_metadata_s33[[
                "OS_status"
              ]][
                idx_s33
              ] ==
                0
            ),
          event_fraction =
            mean(
              train_metadata_s33[[
                "OS_status"
              ]][
                idx_s33
              ] ==
                1
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )


cat(
  "\nCross-validation fold balance:\n"
)


print(
  fold_balance_s33,
  row.names = FALSE
)


# ------------------------------------------------------------
# 3.4 Cross-validated LASSO-Cox
# ------------------------------------------------------------


cv_lasso_s33 <-
  glmnet::cv.glmnet(
    x =
      x_train_z_s33,
    y =
      y_train_s33,
    family =
      "cox",
    alpha =
      1,
    foldid =
      foldid_s33,
    type.measure =
      "deviance",
    standardize =
      FALSE,
    cox.ties =
      "breslow"
  )


cat(
  "\nlambda.min:\n"
)


print(
  cv_lasso_s33$lambda.min
)


cat(
  "\nlambda.1se:\n"
)


print(
  cv_lasso_s33$lambda.1se
)


# ------------------------------------------------------------
# 3.5 Extract coefficients
# ------------------------------------------------------------


coef_min_matrix_s33 <-
  as.matrix(
    coef(
      cv_lasso_s33,
      s = "lambda.min"
    )
  )


coef_1se_matrix_s33 <-
  as.matrix(
    coef(
      cv_lasso_s33,
      s = "lambda.1se"
    )
  )


coef_min_s33 <-
  data.frame(
    gene =
      rownames(
        coef_min_matrix_s33
      ),
    coefficient =
      as.numeric(
        coef_min_matrix_s33[
          ,
          1
        ]
      ),
    stringsAsFactors = FALSE
  )


coef_min_s33$selected <-
  coef_min_s33[[
    "coefficient"
  ]] !=
  0


coef_1se_s33 <-
  data.frame(
    gene =
      rownames(
        coef_1se_matrix_s33
      ),
    coefficient =
      as.numeric(
        coef_1se_matrix_s33[
          ,
          1
        ]
      ),
    stringsAsFactors = FALSE
  )


coef_1se_s33$selected <-
  coef_1se_s33[[
    "coefficient"
  ]] !=
  0


selected_min_s33 <-
  coef_min_s33[
    coef_min_s33[[
      "selected"
    ]],
    ,
    drop = FALSE
  ]


selected_1se_s33 <-
  coef_1se_s33[
    coef_1se_s33[[
      "selected"
    ]],
    ,
    drop = FALSE
  ]


selected_min_s33 <-
  selected_min_s33[
    order(
      -abs(
        selected_min_s33[[
          "coefficient"
        ]]
      )
    ),
    ,
    drop = FALSE
  ]


selected_1se_s33 <-
  selected_1se_s33[
    order(
      -abs(
        selected_1se_s33[[
          "coefficient"
        ]]
      )
    ),
    ,
    drop = FALSE
  ]


cat(
  "\nGenes selected at lambda.min:\n"
)


print(
  selected_min_s33,
  row.names = FALSE
)


cat(
  "\nGenes selected at lambda.1se:\n"
)


print(
  selected_1se_s33,
  row.names = FALSE
)


stopifnot(
  nrow(
    selected_min_s33
  ) ==
    4
)


stopifnot(
  nrow(
    selected_1se_s33
  ) ==
    0
)


stopifnot(
  all(
    selected_min_s33[[
      "gene"
    ]] %in%
      c(
        "GCLC",
        "EGLN3",
        "DSG2",
        "CLSPN"
      )
  )
)


# ------------------------------------------------------------
# 3.6 Apparent training performance
# ------------------------------------------------------------


train_lp_min_s33 <-
  as.numeric(
    predict(
      cv_lasso_s33,
      newx =
        x_train_z_s33,
      s =
        "lambda.min",
      type =
        "link"
    )
  )


train_lp_1se_s33 <-
  as.numeric(
    predict(
      cv_lasso_s33,
      newx =
        x_train_z_s33,
      s =
        "lambda.1se",
      type =
        "link"
    )
  )


train_eval_s33 <-
  train_metadata_s33


train_eval_s33$LP_min <-
  train_lp_min_s33


train_cox_min_s33 <-
  survival::coxph(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      LP_min,
    data =
      train_eval_s33,
    ties =
      "breslow",
    x =
      TRUE
  )


train_cox_min_summary_s33 <-
  summary(
    train_cox_min_s33
  )


apparent_training_summary_s33 <-
  data.frame(
    model =
      c(
        "lambda.min",
        "lambda.1se"
      ),
    n_genes =
      c(
        nrow(
          selected_min_s33
        ),
        nrow(
          selected_1se_s33
        )
      ),
    lambda =
      c(
        cv_lasso_s33$lambda.min,
        cv_lasso_s33$lambda.1se
      ),
    HR_per_1_unit_LP =
      c(
        exp(
          coef(
            train_cox_min_s33
          )[
            1
          ]
        ),
        NA_real_
      ),
    Cox_p =
      c(
        train_cox_min_summary_s33$coefficients[
          1,
          "Pr(>|z|)"
        ],
        NA_real_
      ),
    apparent_C_index =
      c(
        train_cox_min_summary_s33$concordance[
          1
        ],
        0.5
      ),
    stringsAsFactors = FALSE
  )


cat(
  "\nApparent training performance:\n"
)


print(
  apparent_training_summary_s33,
  row.names = FALSE
)


# ------------------------------------------------------------
# 3.7 Pre-test model lock
#
# lambda.1se is a null model.
# lambda.min is therefore locked BEFORE held-out evaluation.
# ------------------------------------------------------------


locked_lambda_rule_s33 <-
  "lambda.min"


locked_lambda_value_s33 <-
  cv_lasso_s33$lambda.min


locked_coef_s33 <-
  selected_min_s33[
    ,
    c(
      "gene",
      "coefficient"
    ),
    drop = FALSE
  ]


locked_model_reason_s33 <-
  paste0(
    "lambda.min was locked before held-out test evaluation because ",
    "lambda.1se selected zero predictors and therefore represented ",
    "the null model. lambda.min selected four predictors and was ",
    "retained as the non-null sparse LASSO-Cox sensitivity model. ",
    "The held-out test cohort was not used for this decision."
  )


model_lock_summary_s33 <-
  data.frame(
    item =
      c(
        "locked_lambda_rule",
        "locked_lambda_value",
        "n_selected_locked_model",
        "lambda_1se",
        "n_selected_lambda_1se",
        "model_selection_timing",
        "model_selection_reason"
      ),
    value =
      c(
        locked_lambda_rule_s33,
        as.character(
          locked_lambda_value_s33
        ),
        as.character(
          nrow(
            locked_coef_s33
          )
        ),
        as.character(
          cv_lasso_s33$lambda.1se
        ),
        as.character(
          nrow(
            selected_1se_s33
          )
        ),
        "Before held-out test evaluation",
        locked_model_reason_s33
      ),
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# 3.8 Save training/model-lock outputs
# ------------------------------------------------------------


cv_fold_assignment_s33 <-
  data.frame(
    patient_id =
      train_metadata_s33[[
        "patient_id"
      ]],
    OS_status =
      train_metadata_s33[[
        "OS_status"
      ]],
    fold =
      foldid_s33,
    cv_seed =
      cv_seed_s33,
    stringsAsFactors = FALSE
  )


write.csv(
  cv_fold_assignment_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_training_CV_fold_assignment.csv"
  ),
  row.names = FALSE
)


write.csv(
  coef_min_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_lambda_min_all_coefficients.csv"
  ),
  row.names = FALSE
)


write.csv(
  selected_min_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_lambda_min_selected_genes.csv"
  ),
  row.names = FALSE
)


write.csv(
  coef_1se_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_lambda_1se_all_coefficients.csv"
  ),
  row.names = FALSE
)


write.csv(
  selected_1se_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_lambda_1se_selected_genes.csv"
  ),
  row.names = FALSE
)


write.csv(
  apparent_training_summary_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_apparent_training_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  locked_coef_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_locked_lambda_min_coefficients.csv"
  ),
  row.names = FALSE
)


write.csv(
  model_lock_summary_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_pretest_model_lock.csv"
  ),
  row.names = FALSE
)


lambda_summary_s33 <-
  data.frame(
    metric =
      c(
        "lambda.min",
        "lambda.1se",
        "n_selected_lambda.min",
        "n_selected_lambda.1se",
        "CV_seed",
        "CV_folds"
      ),
    value =
      c(
        cv_lasso_s33$lambda.min,
        cv_lasso_s33$lambda.1se,
        nrow(
          selected_min_s33
        ),
        nrow(
          selected_1se_s33
        ),
        cv_seed_s33,
        nfolds_s33
      ),
    stringsAsFactors = FALSE
  )


write.csv(
  lambda_summary_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_lambda_summary.csv"
  ),
  row.names = FALSE
)


saveRDS(
  cv_lasso_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_training_cv_model.rds"
  )
)


cat(
  "\n========================================\n",
  "SECTION 3 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 4
# Held-out evaluation of locked four-gene model
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 4: HELD-OUT TEST EVALUATION\n",
  "========================================\n"
)


# ------------------------------------------------------------
# 4.1 Locked coefficients
# ------------------------------------------------------------


locked_genes_s33 <-
  locked_coef_s33[[
    "gene"
  ]]


locked_beta_s33 <-
  locked_coef_s33[[
    "coefficient"
  ]]


names(
  locked_beta_s33
) <-
  locked_genes_s33


# ------------------------------------------------------------
# 4.2 Locked linear predictors
# ------------------------------------------------------------


train_lp_s33 <-
  as.numeric(
    x_train_z_s33[
      ,
      locked_genes_s33,
      drop = FALSE
    ] %*%
      locked_beta_s33[
        locked_genes_s33
      ]
  )


test_lp_s33 <-
  as.numeric(
    x_test_z_s33[
      ,
      locked_genes_s33,
      drop = FALSE
    ] %*%
      locked_beta_s33[
        locked_genes_s33
      ]
  )


stopifnot(
  length(
    train_lp_s33
  ) ==
    352
)


stopifnot(
  length(
    test_lp_s33
  ) ==
    152
)


# ------------------------------------------------------------
# 4.3 Continuous held-out Cox association
# ------------------------------------------------------------


test_eval_s33 <-
  test_metadata_s33


test_eval_s33$LASSO_LP <-
  test_lp_s33


test_cox_continuous_s33 <-
  survival::coxph(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      LASSO_LP,
    data =
      test_eval_s33,
    ties =
      "breslow",
    x =
      TRUE
  )


test_cox_summary_s33 <-
  summary(
    test_cox_continuous_s33
  )


test_hr_s33 <-
  unname(
    exp(
      coef(
        test_cox_continuous_s33
      )[
        1
      ]
    )
  )


test_ci_s33 <-
  unname(
    exp(
      confint(
        test_cox_continuous_s33
      )[
        1,
      ]
    )
  )


test_p_s33 <-
  test_cox_summary_s33$coefficients[
    1,
    "Pr(>|z|)"
  ]


test_c_index_s33 <-
  unname(
    test_cox_summary_s33$concordance[
      1
    ]
  )


test_c_index_se_s33 <-
  unname(
    test_cox_summary_s33$concordance[
      2
    ]
  )


cat(
  "\nHeld-out continuous performance:\n"
)


print(
  data.frame(
    HR =
      test_hr_s33,
    CI_lower =
      test_ci_s33[
        1
      ],
    CI_upper =
      test_ci_s33[
        2
      ],
    Cox_p =
      test_p_s33,
    C_index =
      test_c_index_s33,
    C_index_SE =
      test_c_index_se_s33
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 4.4 PH diagnostic
# ------------------------------------------------------------


test_ph_s33 <-
  survival::cox.zph(
    test_cox_continuous_s33,
    transform =
      "km"
  )


cat(
  "\nHeld-out PH diagnostic:\n"
)


print(
  test_ph_s33
)


test_ph_table_s33 <-
  as.data.frame(
    test_ph_s33$table
  )


test_ph_table_s33$term <-
  rownames(
    test_ph_table_s33
  )


rownames(
  test_ph_table_s33
) <-
  NULL


test_ph_table_s33 <-
  test_ph_table_s33[
    ,
    c(
      "term",
      setdiff(
        colnames(
          test_ph_table_s33
        ),
        "term"
      )
    ),
    drop = FALSE
  ]


# ------------------------------------------------------------
# 4.5 Training-derived median cutoff
# ------------------------------------------------------------


training_cutoff_s33 <-
  median(
    train_lp_s33
  )


train_risk_group_s33 <-
  factor(
    ifelse(
      train_lp_s33 >
        training_cutoff_s33,
      "High LASSO",
      "Low LASSO"
    ),
    levels =
      c(
        "Low LASSO",
        "High LASSO"
      )
  )


test_risk_group_s33 <-
  factor(
    ifelse(
      test_lp_s33 >
        training_cutoff_s33,
      "High LASSO",
      "Low LASSO"
    ),
    levels =
      c(
        "Low LASSO",
        "High LASSO"
      )
  )


cat(
  "\nTraining-derived cutoff:\n"
)


print(
  training_cutoff_s33
)


cat(
  "\nHeld-out risk groups:\n"
)


print(
  table(
    test_risk_group_s33
  )
)


# ------------------------------------------------------------
# 4.6 Held-out KM / log-rank
# ------------------------------------------------------------


test_eval_s33$risk_group <-
  test_risk_group_s33


test_km_s33 <-
  survival::survfit(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    data =
      test_eval_s33
  )


test_logrank_s33 <-
  survival::survdiff(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    data =
      test_eval_s33
  )


test_logrank_p_s33 <-
  stats::pchisq(
    test_logrank_s33$chisq,
    df = 1,
    lower.tail = FALSE
  )


test_cox_group_s33 <-
  survival::coxph(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    data =
      test_eval_s33,
    ties =
      "breslow",
    x =
      TRUE
  )


test_cox_group_summary_s33 <-
  summary(
    test_cox_group_s33
  )


group_hr_s33 <-
  unname(
    exp(
      coef(
        test_cox_group_s33
      )[
        1
      ]
    )
  )


group_ci_s33 <-
  unname(
    exp(
      confint(
        test_cox_group_s33
      )[
        1,
      ]
    )
  )


group_p_s33 <-
  test_cox_group_summary_s33$coefficients[
    1,
    "Pr(>|z|)"
  ]


# ------------------------------------------------------------
# 4.7 Risk-group event summary
# ------------------------------------------------------------


test_group_summary_s33 <-
  do.call(
    rbind,
    lapply(
      levels(
        test_risk_group_s33
      ),
      function(group_s33) {
        
        idx_s33 <-
          test_risk_group_s33 ==
          group_s33
        
        
        data.frame(
          risk_group =
            group_s33,
          n =
            sum(
              idx_s33
            ),
          events =
            sum(
              test_metadata_s33[[
                "OS_status"
              ]][
                idx_s33
              ] ==
                1
            ),
          censored =
            sum(
              test_metadata_s33[[
                "OS_status"
              ]][
                idx_s33
              ] ==
                0
            ),
          median_LP =
            median(
              test_lp_s33[
                idx_s33
              ]
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )


# ------------------------------------------------------------
# 4.8 Time-dependent ROC
# ------------------------------------------------------------


roc_times_s33 <-
  c(
    365,
    1095,
    1825
  )


roc_labels_s33 <-
  c(
    "1-year",
    "3-year",
    "5-year"
  )


roc_support_s33 <-
  data.frame(
    horizon =
      roc_labels_s33,
    time_days =
      roc_times_s33,
    observed_events_by_horizon =
      vapply(
        roc_times_s33,
        function(t_s33) {
          sum(
            test_metadata_s33[[
              "OS_status"
            ]] ==
              1 &
              test_metadata_s33[[
                "OS_time"
              ]] <=
              t_s33
          )
        },
        numeric(1)
      ),
    known_event_free_beyond_horizon =
      vapply(
        roc_times_s33,
        function(t_s33) {
          sum(
            test_metadata_s33[[
              "OS_time"
            ]] >
              t_s33
          )
        },
        numeric(1)
      ),
    stringsAsFactors = FALSE
  )


test_timeROC_s33 <-
  timeROC::timeROC(
    T =
      test_metadata_s33[[
        "OS_time"
      ]],
    delta =
      test_metadata_s33[[
        "OS_status"
      ]],
    marker =
      test_lp_s33,
    cause =
      1,
    weighting =
      "marginal",
    times =
      roc_times_s33,
    iid =
      TRUE
  )


test_timeROC_ci_s33 <-
  confint(
    test_timeROC_s33,
    level = 0.95
  )


auc_ci_raw_s33 <-
  test_timeROC_ci_s33$CI_AUC


# ------------------------------------------------------------
# 4.9 Correct timeROC CI scale when returned as percentages
# ------------------------------------------------------------


if (
  max(
    auc_ci_raw_s33,
    na.rm = TRUE
  ) >
  1
) {
  
  auc_ci_proportion_s33 <-
    auc_ci_raw_s33 /
    100
  
  ci_scale_detected_s33 <-
    "Percentage scale converted to proportion scale"
  
} else {
  
  auc_ci_proportion_s33 <-
    auc_ci_raw_s33
  
  ci_scale_detected_s33 <-
    "Already on proportion scale"
}


stopifnot(
  all(
    auc_ci_proportion_s33 >=
      0,
    na.rm = TRUE
  )
)


stopifnot(
  all(
    auc_ci_proportion_s33 <=
      1,
    na.rm = TRUE
  )
)


time_auc_summary_s33 <-
  data.frame(
    horizon =
      roc_labels_s33,
    time_days =
      roc_times_s33,
    AUC =
      as.numeric(
        test_timeROC_s33$AUC
      ),
    CI_lower =
      as.numeric(
        auc_ci_proportion_s33[
          ,
          1
        ]
      ),
    CI_upper =
      as.numeric(
        auc_ci_proportion_s33[
          ,
          2
        ]
      ),
    observed_events_by_horizon =
      roc_support_s33[[
        "observed_events_by_horizon"
      ]],
    known_event_free_beyond_horizon =
      roc_support_s33[[
        "known_event_free_beyond_horizon"
      ]],
    stringsAsFactors = FALSE
  )


cat(
  "\nHeld-out time-dependent AUC:\n"
)


print(
  time_auc_summary_s33,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4.10 Performance tables
# ------------------------------------------------------------


continuous_performance_s33 <-
  data.frame(
    cohort =
      "Held-out test",
    n =
      nrow(
        test_metadata_s33
      ),
    events =
      sum(
        test_metadata_s33[[
          "OS_status"
        ]]
      ),
    locked_model =
      "lambda.min",
    n_genes =
      nrow(
        locked_coef_s33
      ),
    HR_per_1_unit_LP =
      test_hr_s33,
    HR_CI_lower =
      test_ci_s33[
        1
      ],
    HR_CI_upper =
      test_ci_s33[
        2
      ],
    Cox_p =
      test_p_s33,
    Harrell_C_index =
      test_c_index_s33,
    C_index_SE =
      test_c_index_se_s33,
    stringsAsFactors = FALSE
  )


grouped_performance_s33 <-
  data.frame(
    cutoff_source =
      "Training-set median LP",
    cutoff =
      training_cutoff_s33,
    n_low =
      sum(
        test_risk_group_s33 ==
          "Low LASSO"
      ),
    n_high =
      sum(
        test_risk_group_s33 ==
          "High LASSO"
      ),
    HR_high_vs_low =
      group_hr_s33,
    HR_CI_lower =
      group_ci_s33[
        1
      ],
    HR_CI_upper =
      group_ci_s33[
        2
      ],
    Cox_p =
      group_p_s33,
    logrank_p =
      test_logrank_p_s33,
    stringsAsFactors = FALSE
  )


# ------------------------------------------------------------
# 4.11 Patient-level scores
# ------------------------------------------------------------


train_scores_s33 <-
  data.frame(
    patient_id =
      train_metadata_s33[[
        "patient_id"
      ]],
    OS_time =
      train_metadata_s33[[
        "OS_time"
      ]],
    OS_status =
      train_metadata_s33[[
        "OS_status"
      ]],
    LASSO_LP =
      train_lp_s33,
    risk_group =
      as.character(
        train_risk_group_s33
      ),
    cohort =
      "Training",
    cutoff_source =
      "Training median",
    cutoff =
      training_cutoff_s33,
    stringsAsFactors = FALSE
  )


test_scores_s33 <-
  data.frame(
    patient_id =
      test_metadata_s33[[
        "patient_id"
      ]],
    OS_time =
      test_metadata_s33[[
        "OS_time"
      ]],
    OS_status =
      test_metadata_s33[[
        "OS_status"
      ]],
    LASSO_LP =
      test_lp_s33,
    risk_group =
      as.character(
        test_risk_group_s33
      ),
    cohort =
      "Held-out test",
    cutoff_source =
      "Training median",
    cutoff =
      training_cutoff_s33,
    stringsAsFactors = FALSE
  )


all_scores_s33 <-
  rbind(
    train_scores_s33,
    test_scores_s33
  )


# ------------------------------------------------------------
# 4.12 Save held-out outputs
# ------------------------------------------------------------


write.csv(
  train_scores_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_training_patient_scores.csv"
  ),
  row.names = FALSE
)


write.csv(
  test_scores_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_test_patient_scores.csv"
  ),
  row.names = FALSE
)


write.csv(
  all_scores_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_all_patient_scores_70_30.csv"
  ),
  row.names = FALSE
)


write.csv(
  continuous_performance_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_continuous_performance.csv"
  ),
  row.names = FALSE
)


write.csv(
  grouped_performance_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_grouped_performance.csv"
  ),
  row.names = FALSE
)


write.csv(
  time_auc_summary_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_time_dependent_AUC.csv"
  ),
  row.names = FALSE
)


write.csv(
  roc_support_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_ROC_horizon_support.csv"
  ),
  row.names = FALSE
)


write.csv(
  test_group_summary_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_risk_group_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  test_ph_table_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_PH_diagnostics.csv"
  ),
  row.names = FALSE
)


ci_audit_s33 <-
  data.frame(
    item =
      c(
        "timeROC_AUC_scale",
        "timeROC_confint_raw_scale",
        "saved_CI_scale",
        "correction"
      ),
    value =
      c(
        "0-1 proportion",
        ifelse(
          max(
            auc_ci_raw_s33,
            na.rm = TRUE
          ) >
            1,
          "0-100 percentage",
          "0-1 proportion"
        ),
        "0-1 proportion",
        ci_scale_detected_s33
      ),
    stringsAsFactors = FALSE
  )


write.csv(
  ci_audit_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_timeROC_CI_scale_audit.csv"
  ),
  row.names = FALSE
)


saveRDS(
  test_km_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_KM_fit.rds"
  )
)


saveRDS(
  test_timeROC_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_heldout_timeROC_object.rds"
  )
)


cat(
  "\n========================================\n",
  "SECTION 4 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# SECTION 5
# Final summaries, interpretation, and audit
# ============================================================


cat(
  "\n========================================\n",
  "SECTION 5: FINAL SUMMARY AND AUDIT\n",
  "========================================\n"
)


# ------------------------------------------------------------
# Figure generation
# ------------------------------------------------------------

# Final manuscript figures are intentionally deferred to Script 32.
# Script 28 produces the validated LASSO sensitivity-analysis
# tables and numerical outputs only.


# ------------------------------------------------------------
# 5.1 Final numerical summary
# ------------------------------------------------------------




final_summary_s33 <-
  data.frame(
    metric =
      c(
        "Full TCGA-LUAD cohort",
        "Training cohort",
        "Held-out test cohort",
        "Training events",
        "Held-out test events",
        "Initial candidate genes",
        "lambda.min",
        "lambda.1se",
        "Genes selected at lambda.min",
        "Genes selected at lambda.1se",
        "Apparent training C-index",
        "Held-out test C-index",
        "Held-out continuous Cox p",
        "Held-out High-vs-Low HR",
        "Held-out High-vs-Low Cox p",
        "Held-out log-rank p",
        "Held-out 1-year AUC",
        "Held-out 3-year AUC",
        "Held-out 5-year AUC"
      ),
    value =
      c(
        504,
        352,
        152,
        127,
        55,
        47,
        cv_lasso_s33$lambda.min,
        cv_lasso_s33$lambda.1se,
        nrow(
          selected_min_s33
        ),
        nrow(
          selected_1se_s33
        ),
        apparent_training_summary_s33$apparent_C_index[
          apparent_training_summary_s33$model == "lambda.min"
        ],
        test_c_index_s33,
        test_p_s33,
        group_hr_s33,
        group_p_s33,
        test_logrank_p_s33,
        time_auc_summary_s33[[
          "AUC"
        ]][
          1
        ],
        time_auc_summary_s33[[
          "AUC"
        ]][
          2
        ],
        time_auc_summary_s33[[
          "AUC"
        ]][
          3
        ]
      ),
    stringsAsFactors = FALSE
  )


cat(
  "\nFINAL SCRIPT 28 NUMERICAL SUMMARY:\n"
)


print(
  final_summary_s33,
  row.names = FALSE
)


write.csv(
  final_summary_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_70_30_final_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5.2 Interpretation summary
# ------------------------------------------------------------


final_interpretation_s33 <-
  data.frame(
    topic =
      c(
        "Purpose",
        "Feature selection",
        "Training performance",
        "Held-out continuous performance",
        "Held-out grouped performance",
        "Time-dependent discrimination",
        "Methodological interpretation",
        "Relationship to primary MPCDS"
      ),
    interpretation =
      c(
        "This analysis was performed as a proposal-aligned sensitivity analysis using a reproducible 70/30 train-test split and LASSO-Cox modelling.",
        "Training-only cross-validation selected four genes at lambda.min: GCLC, EGLN3, DSG2, and CLSPN. lambda.1se selected no predictors.",
        paste0(
          "The four-gene lambda.min model showed apparent training discrimination with C-index ",
          sprintf(
            "%.3f",
            apparent_training_summary_s33$apparent_C_index[
              apparent_training_summary_s33$model == "lambda.min"
            ]
          ),
          ". This estimate is apparent and is not a validation result."
        ),
        paste0(
          "In the untouched held-out test cohort, the continuous LASSO score showed weak discrimination (C-index ",
          sprintf(
            "%.3f",
            test_c_index_s33
          ),
          ") and no significant association with overall survival (p = ",
          format.pval(
            test_p_s33,
            digits = 3
          ),
          ")."
        ),
        paste0(
          "Using the training-derived median cutoff, held-out high- and low-risk groups showed essentially no survival separation (HR ",
          sprintf(
            "%.3f",
            group_hr_s33
          ),
          "; log-rank p = ",
          format.pval(
            test_logrank_p_s33,
            digits = 3
          ),
          ")."
        ),
        paste0(
          "Held-out time-dependent AUCs were ",
          sprintf(
            "%.3f",
            time_auc_summary_s33[[
              "AUC"
            ]][
              1
            ]
          ),
          " at 1 year, ",
          sprintf(
            "%.3f",
            time_auc_summary_s33[[
              "AUC"
            ]][
              2
            ]
          ),
          " at 3 years, and ",
          sprintf(
            "%.3f",
            time_auc_summary_s33[[
              "AUC"
            ]][
              3
            ]
          ),
          " at 5 years, indicating limited held-out discrimination."
        ),
        "The simple single-split LASSO workflow did not demonstrate robust prognostic generalization despite selecting a sparse four-gene model in the training cohort.",
        "These findings support retaining the repeated nested-cross-validation ridge MPCDS as the primary prognostic model. Script 33 is a methodological sensitivity analysis and does not replace or invalidate the primary MPCDS framework."
      ),
    stringsAsFactors = FALSE
  )


write.csv(
  final_interpretation_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_LASSO_70_30_interpretation_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5.3 Analysis audit
# ------------------------------------------------------------


analysis_audit_s33 <-
  data.frame(
    analysis_component =
      c(
        "47-gene candidate set",
        "70/30 split",
        "Split stratification",
        "Training-only standardization",
        "Training-only LASSO tuning",
        "Cox tie handling",
        "lambda.min model",
        "lambda.1se model",
        "Pre-test model lock",
        "Held-out continuous evaluation",
        "Training-derived risk cutoff",
        "Held-out KM evaluation",
        "Held-out time-dependent AUC",
        "timeROC CI scale",
        "Test-set leakage",
        "Primary MPCDS replacement",
        "Figure S16"
      ),
    status =
      c(
        "Completed - fixed before modelling",
        "Completed - 352 training / 152 test",
        "Completed - OS event status",
        "Completed",
        "Completed - 10-fold CV",
        "Breslow explicitly locked",
        "4 genes selected",
        "0 genes selected - null model",
        "Completed before test evaluation",
        "Completed",
        "Completed",
        "Completed",
        "Completed",
        "Converted to 0-1 scale when required",
        "None detected",
        "No - sensitivity analysis only",
        "Completed"
      ),
    stringsAsFactors = FALSE
  )


write.csv(
  analysis_audit_s33,
  file.path(
    results_dir_s33,
    "TCGA_LUAD_Script28_analysis_audit.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5.4 Expected-result audit
# ------------------------------------------------------------


stopifnot(
  abs(
    cv_lasso_s33$lambda.min -
      0.06064625
  ) <
    1e-6
)


stopifnot(
  abs(
    cv_lasso_s33$lambda.1se -
      0.1276545
  ) <
    1e-6
)


stopifnot(
  nrow(
    locked_coef_s33
  ) ==
    4
)


stopifnot(
  abs(
    test_c_index_s33 -
      0.5242081
  ) <
    1e-6
)


stopifnot(
  abs(
    test_logrank_p_s33 -
      0.9887106
  ) <
    1e-6
)


stopifnot(
  abs(
    time_auc_summary_s33[[
      "AUC"
    ]][
      1
    ] -
      0.5506775
  ) <
    1e-6
)


stopifnot(
  abs(
    time_auc_summary_s33[[
      "AUC"
    ]][
      2
    ] -
      0.4622812
  ) <
    1e-6
)


stopifnot(
  abs(
    time_auc_summary_s33[[
      "AUC"
    ]][
      3
    ] -
      0.4782228
  ) <
    1e-6
)



stopifnot(
  nrow(
    train_metadata_s33
  ) ==
    352,
  sum(
    train_metadata_s33$OS_status
  ) ==
    127,
  nrow(
    test_metadata_s33
  ) ==
    152,
  sum(
    test_metadata_s33$OS_status
  ) ==
    55,
  identical(
    sort(
      locked_coef_s33$gene
    ),
    sort(
      c(
        "GCLC",
        "EGLN3",
        "DSG2",
        "CLSPN"
      )
    )
  ),
  abs(
    apparent_training_summary_s33$apparent_C_index[
      apparent_training_summary_s33$model == "lambda.min"
    ] -
      0.65746773
  ) <
    1e-6,
  abs(
    test_p_s33 -
      0.69979044
  ) <
    1e-6,
  abs(
    group_hr_s33 -
      0.99612687
  ) <
    1e-6,
  abs(
    group_p_s33 -
      0.98871850
  ) <
    1e-6
)

cat(
  "\\nEstablished LASSO sensitivity regression checks: PASSED\\n"
)


# ------------------------------------------------------------
# 5.5 Final output verification
# ------------------------------------------------------------


final_files_s33 <-
  c(
    file.path(
      results_dir_s33,
      "TCGA_LUAD_LASSO_70_30_final_summary.csv"
    ),
    file.path(
      results_dir_s33,
      "TCGA_LUAD_LASSO_70_30_interpretation_summary.csv"
    ),
    file.path(
      results_dir_s33,
      "TCGA_LUAD_Script28_analysis_audit.csv"
    ),
    file.path(
      results_dir_s33,
      "TCGA_LUAD_LASSO_heldout_time_dependent_AUC.csv"
    ),
    file.path(
      results_dir_s33,
      "TCGA_LUAD_LASSO_timeROC_CI_scale_audit.csv"
    )
  )


cat(
  "\n========================================\n",
  "FINAL SCRIPT 28 OUTPUT FILE CHECK\n",
  "========================================\n"
)


print(
  data.frame(
    file =
      final_files_s33,
    exists =
      file.exists(
        final_files_s33
      ),
    size_bytes =
      file.info(
        final_files_s33
      )$size,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


stopifnot(
  all(
    file.exists(
      final_files_s33
    )
  )
)


# ------------------------------------------------------------
# 5.6 Scientific conclusion
# ------------------------------------------------------------


cat(
  "\n========================================\n",
  "SCRIPT 28 SCIENTIFIC CONCLUSION\n",
  "========================================\n"
)


cat(
  paste0(
    "\nA proposal-aligned 70/30 train-test LASSO-Cox sensitivity ",
    "analysis selected a four-gene training model comprising GCLC, ",
    "EGLN3, DSG2, and CLSPN. Although the model showed apparent ",
    "training discrimination, it demonstrated weak held-out ",
    "prognostic performance, with C-index approximately 0.52, ",
    "no significant continuous survival association, no meaningful ",
    "training-cutoff risk-group separation, and time-dependent AUCs ",
    "near chance levels. These results indicate limited generalization ",
    "of the single-split sparse LASSO model and support retaining the ",
    "repeated nested-cross-validation ridge MPCDS as the primary ",
    "prognostic framework.\n"
  )
)


# ------------------------------------------------------------
# FINAL SUCCESS MARKER
# ------------------------------------------------------------


cat(
  "\n========================================\n",
  "SCRIPT 28 LASSO 70/30 SENSITIVITY ANALYSIS COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)