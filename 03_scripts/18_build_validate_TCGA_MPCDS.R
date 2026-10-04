# ============================================================
# Script 18
# Build and internally validate the TCGA MPCDS prognostic model
#
# Thesis:
# Exploring a Specialized Programmed Cell Death Pattern to
# Predict Prognosis and Treatment Sensitivity in Lung Cancer
# by Machine Learning and Multi-Omics Analysis
#
# Purpose:
# 1. Load the 47 concordant PCD candidate genes.
# 2. Construct LUAD and LUSC survival-expression datasets.
# 3. Internally validate penalized Cox models using repeated
#    nested cross-validation.
# 4. Assess histology-specific prognostic performance.
# 5. Construct the final LUAD ridge-based MPCDS.
# 6. Assess clinical association and PH assumptions.
# 7. Quantify leakage-safe incremental discrimination beyond
#    age, sex, and pathological stage.
# 8. Generate development-cohort risk groups.
# 9. Save all final Script 18 outputs.
#
# Important:
# - LUAD and LUSC are modelled separately.
# - No LUSC prognostic MPCDS is forced if performance is poor.
# - Repeated nested CV provides the main internal validation.
# - Development-cohort KM/Cox results are descriptive and are
#   NOT external validation.
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


# ------------------------------------------------------------
# 1. Required packages and output directory
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "survival",
  "glmnet"
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
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

mpcds_dir <- file.path(
  results_dir,
  "mpcds"
)

survival_dir <- file.path(
  processed_dir,
  "survival"
)

de_dir <- file.path(
  results_dir,
  "differential_expression"
)

ensure_dir(mpcds_dir)


# ------------------------------------------------------------
# 2. Required input files
# ------------------------------------------------------------

luad_surv_expr_file <- file.path(
  survival_dir,
  "TCGA_LUAD_PCD_survival_expression.rds"
)

lusc_surv_expr_file <- file.path(
  survival_dir,
  "TCGA_LUSC_PCD_survival_expression.rds"
)

luad_surv_meta_file <- file.path(
  survival_dir,
  "TCGA_LUAD_PCD_survival_metadata.csv"
)

lusc_surv_meta_file <- file.path(
  survival_dir,
  "TCGA_LUSC_PCD_survival_metadata.csv"
)

candidate_gene_file <- file.path(
  de_dir,
  "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv"
)

required_files <- c(
  luad_surv_expr_file,
  lusc_surv_expr_file,
  luad_surv_meta_file,
  lusc_surv_meta_file,
  candidate_gene_file
)

check_files_exist(required_files)


# ------------------------------------------------------------
# 2A. Load survival-expression datasets
# ------------------------------------------------------------

luad_surv_expr <- readRDS(
  luad_surv_expr_file
)

lusc_surv_expr <- readRDS(
  lusc_surv_expr_file
)

luad_surv_meta <- read.csv(
  luad_surv_meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lusc_surv_meta <- read.csv(
  lusc_surv_meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Load the 47 concordant PCD candidate genes
# ------------------------------------------------------------

candidate_genes <- read.csv(
  candidate_gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

model_genes <- candidate_genes$gene_name


stopifnot(
  length(model_genes) == 47,
  !anyDuplicated(model_genes),
  !anyNA(model_genes),
  all(nzchar(model_genes))
)

stopifnot(
  all(
    model_genes %in%
      rownames(luad_surv_expr)
  )
)

stopifnot(
  all(
    model_genes %in%
      rownames(lusc_surv_expr)
  )
)


# ------------------------------------------------------------
# 4. Validate expression-survival alignment
# ------------------------------------------------------------

stopifnot(
  identical(
    colnames(luad_surv_expr),
    luad_surv_meta$submitter_id
  )
)

stopifnot(
  identical(
    colnames(lusc_surv_expr),
    lusc_surv_meta$submitter_id
  )
)


stopifnot(
  !any(
    is.na(
      luad_surv_meta$OS_time
    )
  )
)

stopifnot(
  !any(
    is.na(
      lusc_surv_meta$OS_time
    )
  )
)

stopifnot(
  !any(
    is.na(
      luad_surv_meta$OS_status
    )
  )
)

stopifnot(
  !any(
    is.na(
      lusc_surv_meta$OS_status
    )
  )
)

stopifnot(
  all(
    luad_surv_meta$OS_time > 0
  )
)

stopifnot(
  all(
    lusc_surv_meta$OS_time > 0
  )
)


# ------------------------------------------------------------
# 5. Construct modelling matrices
# ------------------------------------------------------------

luad_model_expr <- luad_surv_expr[
  model_genes,
  ,
  drop = FALSE
]

lusc_model_expr <- lusc_surv_expr[
  model_genes,
  ,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(luad_model_expr),
    model_genes
  )
)

stopifnot(
  identical(
    rownames(lusc_model_expr),
    model_genes
  )
)


luad_x <- t(
  luad_model_expr
)

lusc_x <- t(
  lusc_model_expr
)


stopifnot(
  identical(
    rownames(luad_x),
    luad_surv_meta$submitter_id
  )
)

stopifnot(
  identical(
    rownames(lusc_x),
    lusc_surv_meta$submitter_id
  )
)

stopifnot(
  ncol(luad_x) == 47
)

stopifnot(
  ncol(lusc_x) == 47
)

stopifnot(
  !any(
    is.na(luad_x)
  )
)

stopifnot(
  !any(
    is.na(lusc_x)
  )
)

stopifnot(
  !any(
    is.infinite(luad_x)
  )
)

stopifnot(
  !any(
    is.infinite(lusc_x)
  )
)


cat(
  "\nMPCDS modelling cohorts:\n",
  "LUAD: ", nrow(luad_x), " patients; events = ",
  sum(luad_surv_meta$OS_status), "\n",
  "LUSC: ", nrow(lusc_x), " patients; events = ",
  sum(lusc_surv_meta$OS_status), "\n",
  "Candidate genes: ", ncol(luad_x), "\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. Survival outcomes
# ------------------------------------------------------------

luad_y <- survival::Surv(
  time =
    luad_surv_meta$OS_time,
  event =
    luad_surv_meta$OS_status
)

lusc_y <- survival::Surv(
  time =
    lusc_surv_meta$OS_time,
  event =
    lusc_surv_meta$OS_status
)


# ------------------------------------------------------------
# 7. Modelling constants
# ------------------------------------------------------------

RANDOM_SEED <- 12345L
NESTED_REPEATS <- 20L
OUTER_FOLDS <- 5L
INNER_FOLDS <- 5L
FINAL_CV_FOLDS <- 10L

alpha_grid <- seq(
  0,
  1,
  by = 0.1
)


# ------------------------------------------------------------
# 8. Stratified fold generator
# ------------------------------------------------------------

make_stratified_folds <- function(
    status,
    k = 5,
    seed = 12345
) {
  
  set.seed(seed)
  
  foldid <- integer(
    length(status)
  )
  
  for (
    event_value in
    sort(
      unique(status)
    )
  ) {
    
    idx <- which(
      status ==
        event_value
    )
    
    idx <- sample(
      idx
    )
    
    foldid[idx] <- rep(
      1:k,
      length.out =
        length(idx)
    )
  }
  
  foldid
}


# ------------------------------------------------------------
# 9. One stratified nested-CV run
# ------------------------------------------------------------

nested_cox_cv_stratified <- function(
    x,
    y,
    status,
    alpha_grid,
    outer_folds = 5,
    inner_folds = 5,
    seed = 12345
) {
  
  outer_foldid <-
    make_stratified_folds(
      status =
        status,
      k =
        outer_folds,
      seed =
        seed
    )
  
  
  outer_results <- vector(
    "list",
    outer_folds
  )
  
  
  for (
    outer_k in
    seq_len(
      outer_folds
    )
  ) {
    
    test_index <- which(
      outer_foldid ==
        outer_k
    )
    
    train_index <- which(
      outer_foldid !=
        outer_k
    )
    
    
    x_train <- x[
      train_index,
      ,
      drop = FALSE
    ]
    
    x_test <- x[
      test_index,
      ,
      drop = FALSE
    ]
    
    y_train <- y[
      train_index
    ]
    
    y_test <- y[
      test_index
    ]
    
    status_train <- status[
      train_index
    ]
    
    
    inner_foldid <-
      make_stratified_folds(
        status =
          status_train,
        k =
          inner_folds,
        seed =
          seed +
          outer_k
      )
    
    
    inner_summary <- vector(
      "list",
      length(alpha_grid)
    )
    
    
    for (
      i in
      seq_along(
        alpha_grid
      )
    ) {
      
      alpha_value <-
        alpha_grid[i]
      
      
      cv_fit <-
        glmnet::cv.glmnet(
          
          x =
            x_train,
          
          y =
            y_train,
          
          family =
            "cox",
          
          alpha =
            alpha_value,
          
          foldid =
            inner_foldid,
          
          type.measure =
            "C",
          
          standardize =
            TRUE,
          
          nlambda =
            100,
          
          cox.ties =
            "breslow"
          
        )
      
      
      lambda_index <-
        which.min(
          abs(
            cv_fit$lambda -
              cv_fit$lambda.min
          )
        )
      
      
      inner_summary[[i]] <-
        data.frame(
          
          alpha =
            alpha_value,
          
          lambda =
            cv_fit$lambda.min,
          
          inner_C =
            cv_fit$cvm[
              lambda_index
            ],
          
          nzero =
            cv_fit$nzero[
              lambda_index
            ]
        )
    }
    
    
    inner_summary <-
      dplyr::bind_rows(
        inner_summary
      )
    
    
    best_index <-
      which.max(
        inner_summary$inner_C
      )
    
    
    best_alpha <-
      inner_summary$alpha[
        best_index
      ]
    
    best_lambda <-
      inner_summary$lambda[
        best_index
      ]
    
    
    outer_fit <-
      glmnet::glmnet(
        
        x =
          x_train,
        
        y =
          y_train,
        
        family =
          "cox",
        
        alpha =
          best_alpha,
        
        lambda =
          best_lambda,
        
        standardize =
          TRUE,
        
        cox.ties =
          "breslow"
      )
    
    
    risk_test <- as.numeric(
      predict(
        outer_fit,
        newx =
          x_test,
        type =
          "link"
      )
    )
    
    
    outer_concordance <-
      survival::concordance(
        
        y_test ~
          risk_test,
        
        reverse =
          TRUE
      )
    
    
    outer_results[[outer_k]] <-
      data.frame(
        
        outer_fold =
          outer_k,
        
        n_train =
          length(
            train_index
          ),
        
        n_test =
          length(
            test_index
          ),
        
        events_test =
          sum(
            status[
              test_index
            ]
          ),
        
        best_alpha =
          best_alpha,
        
        best_lambda =
          best_lambda,
        
        inner_C =
          inner_summary$inner_C[
            best_index
          ],
        
        nzero =
          inner_summary$nzero[
            best_index
          ],
        
        outer_C =
          outer_concordance$concordance
      )
  }
  
  
  dplyr::bind_rows(
    outer_results
  )
}


# ------------------------------------------------------------
# 10. Repeated nested cross-validation
# ------------------------------------------------------------

repeated_nested_cox_cv <- function(
    x,
    y,
    status,
    alpha_grid,
    repeats = 20,
    outer_folds = 5,
    inner_folds = 5,
    seed = 12345
) {
  
  repeat_results <- vector(
    "list",
    repeats
  )
  
  
  for (
    r in
    seq_len(repeats)
  ) {
    
    cat(
      "\nRepeated nested CV:",
      r,
      "of",
      repeats,
      "\n"
    )
    
    
    repeat_seed <-
      seed +
      (r - 1) *
      1000
    
    
    result <-
      nested_cox_cv_stratified(
        
        x =
          x,
        
        y =
          y,
        
        status =
          status,
        
        alpha_grid =
          alpha_grid,
        
        outer_folds =
          outer_folds,
        
        inner_folds =
          inner_folds,
        
        seed =
          repeat_seed
      )
    
    
    result$repeat_id <-
      r
    
    
    repeat_results[[r]] <-
      result
  }
  
  
  dplyr::bind_rows(
    repeat_results
  )
}


# ------------------------------------------------------------
# 11. LUAD repeated nested CV
# ------------------------------------------------------------

set.seed(
  12345
)

luad_repeated_nested <-
  repeated_nested_cox_cv(
    
    x =
      luad_x,
    
    y =
      luad_y,
    
    status =
      luad_surv_meta$OS_status,
    
    alpha_grid =
      alpha_grid,
    
    repeats =
      NESTED_REPEATS,
    
    outer_folds =
      OUTER_FOLDS,
    
    inner_folds =
      INNER_FOLDS,
    
    seed =
      12345
  )


# ------------------------------------------------------------
# 12. LUSC repeated nested CV
# ------------------------------------------------------------

set.seed(
  54321
)

lusc_repeated_nested <-
  repeated_nested_cox_cv(
    
    x =
      lusc_x,
    
    y =
      lusc_y,
    
    status =
      lusc_surv_meta$OS_status,
    
    alpha_grid =
      alpha_grid,
    
    repeats =
      NESTED_REPEATS,
    
    outer_folds =
      OUTER_FOLDS,
    
    inner_folds =
      INNER_FOLDS,
    
    seed =
      54321
  )


stopifnot(
  nrow(
    luad_repeated_nested
  ) == 100
)

stopifnot(
  nrow(
    lusc_repeated_nested
  ) == 100
)


# ------------------------------------------------------------
# 13. Repeat-level nested-CV performance
# ------------------------------------------------------------

luad_repeat_performance <-
  luad_repeated_nested |>
  dplyr::group_by(
    repeat_id
  ) |>
  dplyr::summarise(
    
    mean_C =
      mean(
        outer_C
      ),
    
    median_C =
      median(
        outer_C
      ),
    
    SD_C =
      sd(
        outer_C
      ),
    
    .groups =
      "drop"
  )


lusc_repeat_performance <-
  lusc_repeated_nested |>
  dplyr::group_by(
    repeat_id
  ) |>
  dplyr::summarise(
    
    mean_C =
      mean(
        outer_C
      ),
    
    median_C =
      median(
        outer_C
      ),
    
    SD_C =
      sd(
        outer_C
      ),
    
    .groups =
      "drop"
  )


# ------------------------------------------------------------
# 14. Overall internal-validation summary
# ------------------------------------------------------------

luad_repeated_summary <-
  luad_repeat_performance |>
  dplyr::summarise(
    
    overall_mean_C =
      mean(
        mean_C
      ),
    
    overall_median_C =
      median(
        mean_C
      ),
    
    SD_across_repeats =
      sd(
        mean_C
      ),
    
    Q025 =
      quantile(
        mean_C,
        0.025
      ),
    
    Q975 =
      quantile(
        mean_C,
        0.975
      ),
    
    min_C =
      min(
        mean_C
      ),
    
    max_C =
      max(
        mean_C
      )
  )


lusc_repeated_summary <-
  lusc_repeat_performance |>
  dplyr::summarise(
    
    overall_mean_C =
      mean(
        mean_C
      ),
    
    overall_median_C =
      median(
        mean_C
      ),
    
    SD_across_repeats =
      sd(
        mean_C
      ),
    
    Q025 =
      quantile(
        mean_C,
        0.025
      ),
    
    Q975 =
      quantile(
        mean_C,
        0.975
      ),
    
    min_C =
      min(
        mean_C
      ),
    
    max_C =
      max(
        mean_C
      )
  )


mpcds_internal_validation_summary <-
  data.frame(
    
    histology =
      c(
        "LUAD",
        "LUSC"
      ),
    
    mean_C =
      c(
        luad_repeated_summary$overall_mean_C,
        lusc_repeated_summary$overall_mean_C
      ),
    
    median_C =
      c(
        luad_repeated_summary$overall_median_C,
        lusc_repeated_summary$overall_median_C
      ),
    
    SD_across_repeats =
      c(
        luad_repeated_summary$SD_across_repeats,
        lusc_repeated_summary$SD_across_repeats
      ),
    
    Q025 =
      c(
        luad_repeated_summary$Q025,
        lusc_repeated_summary$Q025
      ),
    
    Q975 =
      c(
        luad_repeated_summary$Q975,
        lusc_repeated_summary$Q975
      ),
    
    min_C =
      c(
        luad_repeated_summary$min_C,
        lusc_repeated_summary$min_C
      ),
    
    max_C =
      c(
        luad_repeated_summary$max_C,
        lusc_repeated_summary$max_C
      )
  )


# ------------------------------------------------------------
# 15. Final LUAD ridge model
#
# Repeated nested CV predominantly selected ridge-like models
# in LUAD. LUSC showed near-random discrimination and is not
# used to construct a prognostic MPCDS.
# ------------------------------------------------------------

set.seed(
  12345
)

luad_final_foldid <-
  make_stratified_folds(
    
    status =
      luad_surv_meta$OS_status,
    
    k =
      FINAL_CV_FOLDS,
    
    seed =
      12345
  )


set.seed(
  12345
)

luad_final_cv <-
  glmnet::cv.glmnet(
    
    x =
      luad_x,
    
    y =
      luad_y,
    
    family =
      "cox",
    
    alpha =
      0,
    
    foldid =
      luad_final_foldid,
    
    type.measure =
      "C",
    
    standardize =
      TRUE,
    
    nlambda =
      100,
    
    cox.ties =
      "breslow"
  )


# ------------------------------------------------------------
# 16. Final LUAD MPCDS coefficients
# ------------------------------------------------------------

luad_final_coef <-
  as.matrix(
    coef(
      luad_final_cv,
      s =
        "lambda.min"
    )
  )


luad_mpcds_coefficients <-
  data.frame(
    
    gene =
      rownames(
        luad_final_coef
      ),
    
    coefficient =
      as.numeric(
        luad_final_coef[
          ,
          1
        ]
      ),
    
    stringsAsFactors =
      FALSE
  ) |>
  dplyr::arrange(
    dplyr::desc(
      abs(
        coefficient
      )
    )
  )


stopifnot(
  nrow(
    luad_mpcds_coefficients
  ) == 47,
  sum(luad_mpcds_coefficients$coefficient != 0) == 47
)


# ------------------------------------------------------------
# 17. Calculate full-development-cohort LUAD MPCDS
# ------------------------------------------------------------

luad_mpcds <- as.numeric(
  predict(
    luad_final_cv,
    newx =
      luad_x,
    s =
      "lambda.min",
    type =
      "link"
  )
)


luad_mpcds_data <-
  data.frame(
    
    patient_id =
      rownames(
        luad_x
      ),
    
    OS_time =
      luad_surv_meta$OS_time,
    
    OS_status =
      luad_surv_meta$OS_status,
    
    MPCDS =
      luad_mpcds,
    
    stringsAsFactors =
      FALSE
  )


stopifnot(
  identical(
    luad_mpcds_data$patient_id,
    luad_surv_meta$submitter_id
  )
)

stopifnot(
  !any(
    is.na(
      luad_mpcds_data$MPCDS
    )
  )
)

stopifnot(
  !any(
    is.infinite(
      luad_mpcds_data$MPCDS
    )
  )
)


# ------------------------------------------------------------
# 18. Standardized MPCDS
# ------------------------------------------------------------

luad_mpcds_data$MPCDS_z <-
  as.numeric(
    scale(
      luad_mpcds_data$MPCDS
    )
  )


# ------------------------------------------------------------
# 19. Clinical annotation
# ------------------------------------------------------------

luad_mpcds_clinical <-
  luad_mpcds_data |>
  dplyr::left_join(
    
    luad_surv_meta |>
      dplyr::select(
        
        submitter_id,
        
        age_at_diagnosis_years,
        
        sex_at_birth,
        
        ajcc_pathologic_stage
      ),
    
    by =
      c(
        "patient_id" =
          "submitter_id"
      )
  )


# ------------------------------------------------------------
# 20. Harmonize pathological stage
# ------------------------------------------------------------

luad_mpcds_clinical <-
  luad_mpcds_clinical |>
  dplyr::mutate(
    
    stage_group =
      dplyr::case_when(
        
        ajcc_pathologic_stage %in%
          c(
            "Stage I",
            "Stage IA",
            "Stage IB"
          ) ~
          "Stage I",
        
        ajcc_pathologic_stage %in%
          c(
            "Stage II",
            "Stage IIA",
            "Stage IIB"
          ) ~
          "Stage II",
        
        ajcc_pathologic_stage %in%
          c(
            "Stage IIIA",
            "Stage IIIB"
          ) ~
          "Stage III",
        
        ajcc_pathologic_stage ==
          "Stage IV" ~
          "Stage IV",
        
        TRUE ~
          NA_character_
      )
  )


luad_mpcds_clinical$stage_group <-
  factor(
    
    luad_mpcds_clinical$stage_group,
    
    levels =
      c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      )
  )


luad_mpcds_clinical$sex_at_birth <-
  factor(
    
    luad_mpcds_clinical$sex_at_birth,
    
    levels =
      c(
        "female",
        "male"
      )
  )


# ------------------------------------------------------------
# 21. Complete-case clinical cohort
# ------------------------------------------------------------

luad_mpcds_complete <-
  luad_mpcds_clinical |>
  dplyr::filter(
    
    !is.na(
      MPCDS_z
    ),
    
    !is.na(
      age_at_diagnosis_years
    ),
    
    !is.na(
      sex_at_birth
    ),
    
    !is.na(
      stage_group
    )
  ) |>
  dplyr::mutate(
    
    age_10yr =
      age_at_diagnosis_years /
      10
  )


# ------------------------------------------------------------
# 22. Conventional clinical + MPCDS Cox model
# ------------------------------------------------------------

luad_clinical_mpcds_cox <-
  survival::coxph(
    
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      MPCDS_z +
      age_at_diagnosis_years +
      sex_at_birth +
      stage_group,
    
    data =
      luad_mpcds_complete,
    
    x =
      TRUE
  )


# ------------------------------------------------------------
# 23. Stage-stratified clinical-adjustment sensitivity model
# ------------------------------------------------------------

luad_stage_stratified_cox <-
  survival::coxph(
    
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      MPCDS_z +
      age_10yr +
      sex_at_birth +
      strata(
        stage_group
      ),
    
    data =
      luad_mpcds_complete,
    
    x =
      TRUE
  )


luad_stage_stratified_result <-
  data.frame(
    
    term =
      "MPCDS_z",
    
    HR =
      exp(
        coef(
          luad_stage_stratified_cox
        )[
          "MPCDS_z"
        ]
      ),
    
    CI_lower =
      exp(
        confint(
          luad_stage_stratified_cox
        )[
          "MPCDS_z",
          1
        ]
      ),
    
    CI_upper =
      exp(
        confint(
          luad_stage_stratified_cox
        )[
          "MPCDS_z",
          2
        ]
      ),
    
    p_value =
      summary(
        luad_stage_stratified_cox
      )$coefficients[
        "MPCDS_z",
        "Pr(>|z|)"
      ]
  )


# ------------------------------------------------------------
# 24. PH diagnostics
# ------------------------------------------------------------

luad_stage_stratified_ph <-
  survival::cox.zph(
    luad_stage_stratified_cox
  )


luad_stage_stratified_ph_table <-
  as.data.frame(
    luad_stage_stratified_ph$table
  )


luad_stage_stratified_ph_table$term <-
  rownames(
    luad_stage_stratified_ph_table
  )


rownames(
  luad_stage_stratified_ph_table
) <- NULL


# ------------------------------------------------------------
# 25. Time-varying MPCDS sensitivity model
# ------------------------------------------------------------

luad_mpcds_timevarying_cox <-
  survival::coxph(
    
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      MPCDS_z +
      tt(MPCDS_z) +
      age_10yr +
      sex_at_birth +
      strata(
        stage_group
      ),
    
    data =
      luad_mpcds_complete,
    
    tt =
      function(
    x,
    t,
    ...
      ) {
        
        x *
          log(
            t /
              365
          )
      },
    
    x =
      TRUE
  )


tv_coef <-
  coef(
    luad_mpcds_timevarying_cox
  )


tv_vcov <-
  vcov(
    luad_mpcds_timevarying_cox
  )


beta_main <-
  tv_coef[
    "MPCDS_z"
  ]


beta_time <-
  tv_coef[
    "tt(MPCDS_z)"
  ]


var_main <-
  tv_vcov[
    "MPCDS_z",
    "MPCDS_z"
  ]


var_time <-
  tv_vcov[
    "tt(MPCDS_z)",
    "tt(MPCDS_z)"
  ]


cov_main_time <-
  tv_vcov[
    "MPCDS_z",
    "tt(MPCDS_z)"
  ]


mpcds_hr_at_time <- function(
    days
) {
  
  log_time <-
    log(
      days /
        365
    )
  
  
  log_hr <-
    beta_main +
    beta_time *
    log_time
  
  
  variance <-
    var_main +
    (
      log_time^2
    ) *
    var_time +
    2 *
    log_time *
    cov_main_time
  
  
  se_log_hr <-
    sqrt(
      variance
    )
  
  
  data.frame(
    
    time_days =
      days,
    
    time_years =
      days /
      365,
    
    HR =
      exp(
        log_hr
      ),
    
    CI_lower =
      exp(
        log_hr -
          1.96 *
          se_log_hr
      ),
    
    CI_upper =
      exp(
        log_hr +
          1.96 *
          se_log_hr
      )
  )
}


luad_mpcds_time_specific_hr <-
  dplyr::bind_rows(
    
    mpcds_hr_at_time(
      365
    ),
    
    mpcds_hr_at_time(
      3 *
        365
    ),
    
    mpcds_hr_at_time(
      5 *
        365
    )
  )


# ------------------------------------------------------------
# 26. Prepare clinical CV dataset
# ------------------------------------------------------------

luad_clinical_cv_data <-
  luad_mpcds_complete |>
  dplyr::select(
    
    patient_id,
    
    OS_time,
    
    OS_status,
    
    age_10yr,
    
    sex_at_birth,
    
    stage_group
  )


rownames(
  luad_clinical_cv_data
) <-
  luad_clinical_cv_data$patient_id


luad_x_clinical <-
  luad_x[
    luad_clinical_cv_data$patient_id,
    ,
    drop = FALSE
  ]


stopifnot(
  identical(
    rownames(
      luad_x_clinical
    ),
    luad_clinical_cv_data$patient_id
  )
)


# ------------------------------------------------------------
# 27. Correct leakage-safe clinical-vs-MPCDS CV function
#
# Stage is entered as an explicit predictor so that stage
# contributes to the held-out linear predictor.
# ------------------------------------------------------------

compare_clinical_mpcds_cv_corrected <- function(
    x,
    clinical_data,
    outer_folds = 5,
    seed = 12345
) {
  
  foldid <-
    make_stratified_folds(
      
      status =
        clinical_data$OS_status,
      
      k =
        outer_folds,
      
      seed =
        seed
    )
  
  
  results <- vector(
    "list",
    outer_folds
  )
  
  
  for (
    k in
    seq_len(
      outer_folds
    )
  ) {
    
    test_idx <- which(
      foldid ==
        k
    )
    
    train_idx <- which(
      foldid !=
        k
    )
    
    
    x_train <- x[
      train_idx,
      ,
      drop = FALSE
    ]
    
    x_test <- x[
      test_idx,
      ,
      drop = FALSE
    ]
    
    
    train_data <- clinical_data[
      train_idx,
      ,
      drop = FALSE
    ]
    
    test_data <- clinical_data[
      test_idx,
      ,
      drop = FALSE
    ]
    
    
    train_data$stage_group <-
      factor(
        
        train_data$stage_group,
        
        levels =
          c(
            "Stage I",
            "Stage II",
            "Stage III",
            "Stage IV"
          )
      )
    
    
    test_data$stage_group <-
      factor(
        
        test_data$stage_group,
        
        levels =
          c(
            "Stage I",
            "Stage II",
            "Stage III",
            "Stage IV"
          )
      )
    
    
    train_data$sex_at_birth <-
      factor(
        
        train_data$sex_at_birth,
        
        levels =
          c(
            "female",
            "male"
          )
      )
    
    
    test_data$sex_at_birth <-
      factor(
        
        test_data$sex_at_birth,
        
        levels =
          c(
            "female",
            "male"
          )
      )
    
    
    y_train <-
      survival::Surv(
        
        train_data$OS_time,
        
        train_data$OS_status
      )
    
    
    inner_foldid <-
      make_stratified_folds(
        
        status =
          train_data$OS_status,
        
        k =
          5,
        
        seed =
          seed +
          k
      )
    
    
    ridge_cv <-
      glmnet::cv.glmnet(
        
        x =
          x_train,
        
        y =
          y_train,
        
        family =
          "cox",
        
        alpha =
          0,
        
        foldid =
          inner_foldid,
        
        type.measure =
          "C",
        
        standardize =
          TRUE,
        
        nlambda =
          100,
        
        cox.ties =
          "breslow"
      )
    
    
    ridge_fit <-
      glmnet::glmnet(
        
        x =
          x_train,
        
        y =
          y_train,
        
        family =
          "cox",
        
        alpha =
          0,
        
        lambda =
          ridge_cv$lambda.min,
        
        standardize =
          TRUE,
        
        cox.ties =
          "breslow"
      )
    
    
    train_data$MPCDS_cv <-
      as.numeric(
        predict(
          ridge_fit,
          newx =
            x_train,
          type =
            "link"
        )
      )
    
    
    test_data$MPCDS_cv <-
      as.numeric(
        predict(
          ridge_fit,
          newx =
            x_test,
          type =
            "link"
        )
      )
    
    
    train_mean <-
      mean(
        train_data$MPCDS_cv
      )
    
    
    train_sd <-
      sd(
        train_data$MPCDS_cv
      )
    
    
    train_data$MPCDS_cv_z <-
      (
        train_data$MPCDS_cv -
          train_mean
      ) /
      train_sd
    
    
    test_data$MPCDS_cv_z <-
      (
        test_data$MPCDS_cv -
          train_mean
      ) /
      train_sd
    
    
    clinical_fit <-
      survival::coxph(
        
        survival::Surv(
          OS_time,
          OS_status
        ) ~
          age_10yr +
          sex_at_birth +
          stage_group,
        
        data =
          train_data,
        
        ties =
          "efron"
      )
    
    
    combined_fit <-
      survival::coxph(
        
        survival::Surv(
          OS_time,
          OS_status
        ) ~
          MPCDS_cv_z +
          age_10yr +
          sex_at_birth +
          stage_group,
        
        data =
          train_data,
        
        ties =
          "efron"
      )
    
    
    clinical_lp <-
      as.numeric(
        predict(
          clinical_fit,
          newdata =
            test_data,
          type =
            "lp"
        )
      )
    
    
    combined_lp <-
      as.numeric(
        predict(
          combined_fit,
          newdata =
            test_data,
          type =
            "lp"
        )
      )
    
    
    test_surv <-
      survival::Surv(
        
        test_data$OS_time,
        
        test_data$OS_status
      )
    
    
    clinical_c <-
      survival::concordance(
        
        test_surv ~
          clinical_lp,
        
        reverse =
          TRUE
        
      )$concordance
    
    
    combined_c <-
      survival::concordance(
        
        test_surv ~
          combined_lp,
        
        reverse =
          TRUE
        
      )$concordance
    
    
    results[[k]] <-
      data.frame(
        
        outer_fold =
          k,
        
        n_test =
          length(
            test_idx
          ),
        
        events =
          sum(
            test_data$OS_status
          ),
        
        lambda =
          ridge_cv$lambda.min,
        
        clinical_C =
          clinical_c,
        
        combined_C =
          combined_c,
        
        delta_C =
          combined_c -
          clinical_c
      )
  }
  
  
  dplyr::bind_rows(
    results
  )
}


# ------------------------------------------------------------
# 28. Repeated clinical-vs-MPCDS CV
# ------------------------------------------------------------

repeated_clinical_mpcds_cv <- function(
    x,
    clinical_data,
    repeats = 20,
    outer_folds = 5,
    seed = 12345
) {
  
  repeat_results <- vector(
    "list",
    repeats
  )
  
  
  for (
    r in
    seq_len(
      repeats
    )
  ) {
    
    cat(
      "\nClinical incremental CV:",
      r,
      "of",
      repeats,
      "\n"
    )
    
    
    repeat_seed <-
      seed +
      (
        r - 1
      ) *
      1000
    
    
    result <-
      compare_clinical_mpcds_cv_corrected(
        
        x =
          x,
        
        clinical_data =
          clinical_data,
        
        outer_folds =
          outer_folds,
        
        seed =
          repeat_seed
      )
    
    
    result$repeat_id <-
      r
    
    
    repeat_results[[r]] <-
      result
  }
  
  
  dplyr::bind_rows(
    repeat_results
  )
}


# ------------------------------------------------------------
# 29. Run repeated clinical incremental-value analysis
# ------------------------------------------------------------

set.seed(
  12345
)

luad_clinical_mpcds_repeated <-
  repeated_clinical_mpcds_cv(
    
    x =
      luad_x_clinical,
    
    clinical_data =
      luad_clinical_cv_data,
    
    repeats =
      20,
    
    outer_folds =
      OUTER_FOLDS,
    
    seed =
      12345
  )


stopifnot(
  nrow(
    luad_clinical_mpcds_repeated
  ) == 100
)


# ------------------------------------------------------------
# 30. Repeat-level clinical incremental performance
# ------------------------------------------------------------

luad_clinical_mpcds_repeat_summary <-
  luad_clinical_mpcds_repeated |>
  dplyr::group_by(
    repeat_id
  ) |>
  dplyr::summarise(
    
    mean_clinical_C =
      mean(
        clinical_C
      ),
    
    mean_combined_C =
      mean(
        combined_C
      ),
    
    mean_delta_C =
      mean(
        delta_C
      ),
    
    median_delta_C =
      median(
        delta_C
      ),
    
    folds_improved =
      sum(
        delta_C >
          0
      ),
    
    .groups =
      "drop"
  )


luad_incremental_summary <-
  luad_clinical_mpcds_repeat_summary |>
  dplyr::summarise(
    
    overall_clinical_C =
      mean(
        mean_clinical_C
      ),
    
    overall_combined_C =
      mean(
        mean_combined_C
      ),
    
    overall_delta_C =
      mean(
        mean_delta_C
      ),
    
    median_delta_C =
      median(
        mean_delta_C
      ),
    
    SD_delta_C =
      sd(
        mean_delta_C
      ),
    
    Q025_delta_C =
      quantile(
        mean_delta_C,
        0.025
      ),
    
    Q975_delta_C =
      quantile(
        mean_delta_C,
        0.975
      ),
    
    min_delta_C =
      min(
        mean_delta_C
      ),
    
    max_delta_C =
      max(
        mean_delta_C
      ),
    
    repeats_positive =
      sum(
        mean_delta_C >
          0
      ),
    
    repeats_total =
      dplyr::n()
  )


# ------------------------------------------------------------
# 31. Development-cohort MPCDS risk groups
# ------------------------------------------------------------

luad_mpcds_cutoff <-
  median(
    luad_mpcds_data$MPCDS
  )


luad_mpcds_data$risk_group <-
  ifelse(
    
    luad_mpcds_data$MPCDS >
      luad_mpcds_cutoff,
    
    "High MPCDS",
    
    "Low MPCDS"
  )


luad_mpcds_data$risk_group <-
  factor(
    
    luad_mpcds_data$risk_group,
    
    levels =
      c(
        "Low MPCDS",
        "High MPCDS"
      )
  )


# ------------------------------------------------------------
# 32. Development-cohort KM analysis
# ------------------------------------------------------------

luad_mpcds_km <-
  survival::survfit(
    
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    
    data =
      luad_mpcds_data
  )


luad_mpcds_logrank <-
  survival::survdiff(
    
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    
    data =
      luad_mpcds_data
  )


luad_mpcds_logrank_p <-
  1 -
  pchisq(
    
    luad_mpcds_logrank$chisq,
    
    df =
      length(
        luad_mpcds_logrank$n
      ) -
      1
  )


luad_mpcds_group_cox <-
  survival::coxph(
    
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    
    data =
      luad_mpcds_data
  )


luad_mpcds_group_result <-
  data.frame(
    
    comparison =
      "High MPCDS vs Low MPCDS",
    
    HR =
      exp(
        coef(
          luad_mpcds_group_cox
        )[1]
      ),
    
    CI_lower =
      exp(
        confint(
          luad_mpcds_group_cox
        )[
          1,
          1
        ]
      ),
    
    CI_upper =
      exp(
        confint(
          luad_mpcds_group_cox
        )[
          1,
          2
        ]
      ),
    
    p_value =
      summary(
        luad_mpcds_group_cox
      )$coefficients[
        1,
        "Pr(>|z|)"
      ]
  )


# ------------------------------------------------------------
# 33. KM summary table
# ------------------------------------------------------------

luad_km_table <-
  summary(
    luad_mpcds_km
  )$table


luad_mpcds_km_summary <-
  data.frame(
    
    risk_group =
      c(
        "Low MPCDS",
        "High MPCDS"
      ),
    
    n =
      as.numeric(
        luad_km_table[
          ,
          "records"
        ]
      ),
    
    events =
      as.numeric(
        luad_km_table[
          ,
          "events"
        ]
      ),
    
    median_OS_days =
      as.numeric(
        luad_km_table[
          ,
          "median"
        ]
      ),
    
    median_OS_lower95 =
      as.numeric(
        luad_km_table[
          ,
          "0.95LCL"
        ]
      ),
    
    median_OS_upper95 =
      as.numeric(
        luad_km_table[
          ,
          "0.95UCL"
        ]
      ),
    
    median_cutoff =
      luad_mpcds_cutoff,
    
    logrank_p =
      luad_mpcds_logrank_p
  )


# ------------------------------------------------------------
# 34. Save final results
# ------------------------------------------------------------

write.csv(
  luad_mpcds_coefficients,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_ridge_coefficients.csv"),
  row.names = FALSE
)


write.csv(
  luad_mpcds_data,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_patient_scores.csv"),
  row.names = FALSE
)


write.csv(
  luad_repeated_nested,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_repeated_nested_CV.csv"),
  row.names = FALSE
)


write.csv(
  lusc_repeated_nested,
  file.path(mpcds_dir, "TCGA_LUSC_MPCDS_repeated_nested_CV.csv"),
  row.names = FALSE
)


write.csv(
  luad_repeat_performance,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_repeat_performance.csv"),
  row.names = FALSE
)


write.csv(
  lusc_repeat_performance,
  file.path(mpcds_dir, "TCGA_LUSC_MPCDS_repeat_performance.csv"),
  row.names = FALSE
)


write.csv(
  mpcds_internal_validation_summary,
  file.path(mpcds_dir, "TCGA_MPCDS_internal_validation_summary.csv"),
  row.names = FALSE
)


write.csv(
  luad_stage_stratified_result,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_stage_stratified_Cox.csv"),
  row.names = FALSE
)


write.csv(
  luad_stage_stratified_ph_table,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_PH_diagnostics.csv"),
  row.names = FALSE
)


write.csv(
  luad_mpcds_time_specific_hr,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_time_specific_HR_sensitivity.csv"),
  row.names = FALSE
)


write.csv(
  luad_clinical_mpcds_repeated,
  file.path(mpcds_dir, "TCGA_LUAD_clinical_MPCDS_repeated_CV.csv"),
  row.names = FALSE
)


write.csv(
  luad_clinical_mpcds_repeat_summary,
  file.path(mpcds_dir, "TCGA_LUAD_clinical_MPCDS_repeat_summary.csv"),
  row.names = FALSE
)


write.csv(
  luad_incremental_summary,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_incremental_prediction_summary.csv"),
  row.names = FALSE
)


write.csv(
  luad_mpcds_km_summary,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_KM_summary.csv"),
  row.names = FALSE
)


write.csv(
  luad_mpcds_group_result,
  file.path(mpcds_dir, "TCGA_LUAD_MPCDS_high_vs_low_Cox.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 34A. Save model provenance
# ------------------------------------------------------------

mpcds_model_provenance <- data.frame(
  model = "LUAD ridge-Cox MPCDS",
  candidate_genes = length(model_genes),
  nested_repeats = NESTED_REPEATS,
  outer_folds = OUTER_FOLDS,
  inner_folds = INNER_FOLDS,
  alpha_grid = paste(alpha_grid, collapse = ";"),
  final_alpha = 0,
  final_lambda = luad_final_cv$lambda.min,
  final_cv_folds = FINAL_CV_FOLDS,
  standardize = TRUE,
  glmnet_version = as.character(
    utils::packageVersion("glmnet")
  ),
  survival_version = as.character(
    utils::packageVersion("survival")
  ),
  R_version = R.version.string,
  analysis_date = as.character(Sys.Date()),
  stringsAsFactors = FALSE
)

write.csv(
  mpcds_model_provenance,
  file.path(
    mpcds_dir,
    "TCGA_LUAD_MPCDS_model_provenance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 35. Final verification
# ------------------------------------------------------------

expected_files <-
  c(
    
    "TCGA_LUAD_clinical_MPCDS_repeat_summary.csv",
    
    "TCGA_LUAD_clinical_MPCDS_repeated_CV.csv",
    
    "TCGA_LUAD_MPCDS_high_vs_low_Cox.csv",
    
    "TCGA_LUAD_MPCDS_incremental_prediction_summary.csv",
    
    "TCGA_LUAD_MPCDS_KM_summary.csv",
    
    "TCGA_LUAD_MPCDS_patient_scores.csv",
    
    "TCGA_LUAD_MPCDS_PH_diagnostics.csv",
    
    "TCGA_LUAD_MPCDS_repeat_performance.csv",
    
    "TCGA_LUAD_MPCDS_repeated_nested_CV.csv",
    
    "TCGA_LUAD_MPCDS_ridge_coefficients.csv",
    
    "TCGA_LUAD_MPCDS_stage_stratified_Cox.csv",
    
    "TCGA_LUAD_MPCDS_time_specific_HR_sensitivity.csv",
    
    "TCGA_LUSC_MPCDS_repeat_performance.csv",
    
    "TCGA_LUSC_MPCDS_repeated_nested_CV.csv",
    
    "TCGA_MPCDS_internal_validation_summary.csv"
  )


saved_files <-
  list.files(
    mpcds_dir,
    full.names = FALSE
  )


stopifnot(
  all(
    expected_files %in%
      saved_files
  )
)


# ------------------------------------------------------------
# 36. Final console summary
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "SCRIPT 18 COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "LUAD modelling patients:",
  nrow(
    luad_x
  ),
  "\n"
)


cat(
  "LUSC modelling patients:",
  nrow(
    lusc_x
  ),
  "\n"
)


cat(
  "Candidate PCD genes:",
  ncol(
    luad_x
  ),
  "\n"
)


cat(
  "Final LUAD model:",
  "ridge Cox, alpha = 0",
  "\n"
)


cat(
  "Final LUAD lambda.min:",
  luad_final_cv$lambda.min,
  "\n"
)


cat(
  "LUAD repeated nested-CV mean C-index:",
  round(
    luad_repeated_summary$overall_mean_C,
    4
  ),
  "\n"
)


cat(
  "LUSC repeated nested-CV mean C-index:",
  round(
    lusc_repeated_summary$overall_mean_C,
    4
  ),
  "\n"
)


cat(
  "LUAD clinical-only repeated CV C-index:",
  round(
    luad_incremental_summary$overall_clinical_C,
    4
  ),
  "\n"
)


cat(
  "LUAD clinical + MPCDS repeated CV C-index:",
  round(
    luad_incremental_summary$overall_combined_C,
    4
  ),
  "\n"
)


cat(
  "LUAD incremental Delta C:",
  round(
    luad_incremental_summary$overall_delta_C,
    4
  ),
  "\n"
)


cat(
  "LUAD stage-stratified MPCDS HR:",
  round(
    luad_stage_stratified_result$HR,
    4
  ),
  "\n"
)


cat(
  "LUAD high-vs-low MPCDS HR:",
  round(
    luad_mpcds_group_result$HR,
    4
  ),
  "\n"
)


cat(
  "LUAD KM log-rank p:",
  format(
    luad_mpcds_logrank_p,
    scientific = TRUE
  ),
  "\n"
)


cat(
  "Number of saved MPCDS result files:",
  length(
    expected_files
  ),
  "\n"
)


cat(
  "All expected output files verified:",
  all(
    expected_files %in%
      saved_files
  ),
  "\n"
)


cat(
  "============================================================\n"
)

cat(
  "\n========================================\n",
  "SCRIPT 18 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "Outputs: ", mpcds_dir, "\n",
  "Interpretation rule: nested CV is the primary internal-validation estimate;\n",
  "full-cohort KM/Cox analyses are development-cohort descriptive analyses.\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================