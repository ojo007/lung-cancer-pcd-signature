# ============================================================
# 23_decision_curve_analysis_LUAD_MPCDS.R
#
# Decision curve analysis of the LUAD MPCDS
#
# Primary comparison:
#   1. Clinical model
#   2. Clinical + continuous MPCDS
#
# Evaluation horizons:
#   - 1 year
#   - 3 years
#   - 5 years
#
# Important:
#   - Reuse the established TCGA-LUAD prognostic framework.
#   - Do NOT optimize a new MPCDS cutoff.
#   - Do NOT refit or redefine the MPCDS itself.
#   - DCA is survival-aware and accounts for censoring.
#   - Results are apparent development-cohort estimates.
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


cat(
  "\n========================================\n",
  "SCRIPT 23: DECISION CURVE ANALYSIS OF LUAD MPCDS\n",
  "========================================\n"
)


# ============================================================
# 1. Inspect existing clinical + MPCDS inputs
# ============================================================


candidate_files <-
  c(
    "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv",
    "04_results/mpcds/TCGA_LUAD_MPCDS_incremental_prediction_summary.csv",
    "04_results/mpcds/TCGA_LUAD_clinical_MPCDS_repeat_summary.csv",
    "02_processed_data/survival/TCGA_LUAD_PCD_survival_metadata.csv"
  )


cat(
  "\nCandidate existing inputs:\n"
)


print(
  data.frame(
    file =
      candidate_files,
    exists =
      file.exists(
        candidate_files
      ),
    stringsAsFactors =
      FALSE
  )
)


# ------------------------------------------------------------
# Load TCGA-LUAD survival metadata
# ------------------------------------------------------------

luad_survival_metadata_file <-
  "02_processed_data/survival/TCGA_LUAD_PCD_survival_metadata.csv"


stopifnot(
  file.exists(
    luad_survival_metadata_file
  )
)


luad_survival_metadata <-
  read.csv(
    luad_survival_metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat(
  "\n========================================\n",
  "TCGA-LUAD SURVIVAL METADATA\n",
  "========================================\n"
)


cat(
  "\nDimensions:\n"
)


print(
  dim(
    luad_survival_metadata
  )
)


cat(
  "\nColumn names:\n"
)


print(
  colnames(
    luad_survival_metadata
  )
)


cat(
  "\nFirst five rows:\n"
)


print(
  utils::head(
    luad_survival_metadata,
    5
  )
)


# ------------------------------------------------------------
# Load existing MPCDS patient-level scores
# ------------------------------------------------------------

luad_mpcds_file <-
  "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv"


stopifnot(
  file.exists(
    luad_mpcds_file
  )
)


luad_mpcds_scores <-
  read.csv(
    luad_mpcds_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


cat(
  "\n========================================\n",
  "TCGA-LUAD MPCDS PATIENT SCORES\n",
  "========================================\n"
)


cat(
  "\nDimensions:\n"
)


print(
  dim(
    luad_mpcds_scores
  )
)


cat(
  "\nColumn names:\n"
)


print(
  colnames(
    luad_mpcds_scores
  )
)


# ------------------------------------------------------------
# Script 18 provenance note
# ------------------------------------------------------------
# The canonical Script 18 outputs are already verified above:
#   - TCGA_LUAD_MPCDS_patient_scores.csv
#   - TCGA_LUAD_MPCDS_incremental_prediction_summary.csv
#   - TCGA_LUAD_clinical_MPCDS_repeat_summary.csv
#
# Script 23 deliberately reconstructs the exact clinical and
# clinical + MPCDS Cox models from these verified outputs and
# the survival metadata, so it does not depend on a historical
# Script 18 filename.
# ------------------------------------------------------------


# ------------------------------------------------------------
# Identify likely clinical variables
# ------------------------------------------------------------

clinical_candidate_columns <-
  grep(
    paste0(
      "age|sex|gender|stage|",
      "OS|surv|status|time|",
      "submitter|patient"
    ),
    colnames(
      luad_survival_metadata
    ),
    value = TRUE,
    ignore.case = TRUE
  )


cat(
  "\n========================================\n",
  "LIKELY CLINICAL / SURVIVAL COLUMNS\n",
  "========================================\n"
)


print(
  clinical_candidate_columns
)


# ------------------------------------------------------------
# Patient alignment
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "PATIENT ALIGNMENT CHECK\n",
  "========================================\n"
)


cat(
  "\nMPCDS patient count:\n"
)


print(
  nrow(
    luad_mpcds_scores
  )
)


cat(
  "\nSurvival metadata row count:\n"
)


print(
  nrow(
    luad_survival_metadata
  )
)


cat(
  "\nMPCDS patients found in survival metadata:\n"
)


print(
  sum(
    luad_mpcds_scores$patient_id %in%
      luad_survival_metadata$submitter_id
  )
)


cat(
  "\nMPCDS patients NOT found in survival metadata:\n"
)


print(
  sum(
    !luad_mpcds_scores$patient_id %in%
      luad_survival_metadata$submitter_id
  )
)


stopifnot(
  all(
    luad_mpcds_scores$patient_id %in%
      luad_survival_metadata$submitter_id
  )
)


# ============================================================
# 2. Reconstruct the exact Script 18 clinical MPCDS cohort
# ============================================================


luad_dca_data <-
  merge(
    luad_mpcds_scores,
    luad_survival_metadata[
      ,
      c(
        "submitter_id",
        "age_at_diagnosis_years",
        "sex_at_birth",
        "ajcc_pathologic_stage"
      )
    ],
    by.x =
      "patient_id",
    by.y =
      "submitter_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )


cat(
  "\n========================================\n",
  "JOINED LUAD CLINICAL + MPCDS DATA\n",
  "========================================\n"
)


print(
  dim(
    luad_dca_data
  )
)


# ------------------------------------------------------------
# Survival alignment
# ------------------------------------------------------------

survival_alignment_check <-
  match(
    luad_dca_data$patient_id,
    luad_survival_metadata$submitter_id
  )


stopifnot(
  !anyNA(
    survival_alignment_check
  )
)


stopifnot(
  identical(
    as.numeric(
      luad_dca_data$OS_time
    ),
    as.numeric(
      luad_survival_metadata$OS_time[
        survival_alignment_check
      ]
    )
  )
)


stopifnot(
  identical(
    as.numeric(
      luad_dca_data$OS_status
    ),
    as.numeric(
      luad_survival_metadata$OS_status[
        survival_alignment_check
      ]
    )
  )
)


cat(
  "\nSurvival alignment with original metadata: PASSED\n"
)


# ------------------------------------------------------------
# Reproduce Script 18 stage harmonization exactly
# ------------------------------------------------------------

luad_dca_data$stage_group <-
  ifelse(
    luad_dca_data$ajcc_pathologic_stage %in%
      c(
        "Stage I",
        "Stage IA",
        "Stage IB"
      ),
    "Stage I",
    ifelse(
      luad_dca_data$ajcc_pathologic_stage %in%
        c(
          "Stage II",
          "Stage IIA",
          "Stage IIB"
        ),
      "Stage II",
      ifelse(
        luad_dca_data$ajcc_pathologic_stage %in%
          c(
            "Stage IIIA",
            "Stage IIIB"
          ),
        "Stage III",
        ifelse(
          luad_dca_data$ajcc_pathologic_stage ==
            "Stage IV",
          "Stage IV",
          NA_character_
        )
      )
    )
  )


luad_dca_data$stage_group <-
  factor(
    luad_dca_data$stage_group,
    levels =
      c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      )
  )


luad_dca_data$sex_at_birth <-
  factor(
    luad_dca_data$sex_at_birth,
    levels =
      c(
        "female",
        "male"
      )
  )


cat(
  "\n========================================\n",
  "STAGE HARMONIZATION BEFORE FILTERING\n",
  "========================================\n"
)


print(
  table(
    luad_dca_data$stage_group,
    useNA =
      "ifany"
  )
)


# ------------------------------------------------------------
# Exact Script 18 complete-case cohort
# ------------------------------------------------------------

luad_dca_complete <-
  luad_dca_data[
    !is.na(
      luad_dca_data$MPCDS_z
    ) &
      !is.na(
        luad_dca_data$age_at_diagnosis_years
      ) &
      !is.na(
        luad_dca_data$sex_at_birth
      ) &
      !is.na(
        luad_dca_data$stage_group
      ),
    ,
    drop = FALSE
  ]


luad_dca_complete$age_10yr <-
  luad_dca_complete$age_at_diagnosis_years /
  10


cat(
  "\n========================================\n",
  "RECONSTRUCTED SCRIPT 18 COMPLETE-CASE COHORT\n",
  "========================================\n"
)


cat(
  "\nN:\n"
)


print(
  nrow(
    luad_dca_complete
  )
)


cat(
  "\nOverall survival status:\n"
)


print(
  table(
    luad_dca_complete$OS_status
  )
)


cat(
  "\nStage distribution:\n"
)


print(
  table(
    luad_dca_complete$stage_group
  )
)


cat(
  "\nSex distribution:\n"
)


print(
  table(
    luad_dca_complete$sex_at_birth
  )
)


cat(
  "\nAge summary:\n"
)


print(
  summary(
    luad_dca_complete$age_at_diagnosis_years
  )
)


cat(
  "\nMPCDS_z summary:\n"
)


print(
  summary(
    luad_dca_complete$MPCDS_z
  )
)


cat(
  "\nOS follow-up summary in days:\n"
)


print(
  summary(
    luad_dca_complete$OS_time
  )
)


# ------------------------------------------------------------
# Exact integrity checks
# ------------------------------------------------------------

stopifnot(
  nrow(
    luad_dca_complete
  ) ==
    475
)


stopifnot(
  sum(
    luad_dca_complete$OS_status ==
      1
  ) ==
    171
)


stopifnot(
  sum(
    luad_dca_complete$OS_status ==
      0
  ) ==
    304
)


expected_stage_counts <-
  c(
    "Stage I" =
      258,
    "Stage II" =
      112,
    "Stage III" =
      80,
    "Stage IV" =
      25
  )


observed_stage_counts <-
  table(
    luad_dca_complete$stage_group
  )


stopifnot(
  identical(
    as.integer(
      observed_stage_counts[
        names(
          expected_stage_counts
        )
      ]
    ),
    as.integer(
      expected_stage_counts
    )
  )
)


cat(
  "\nScript 18 cohort reconstruction: PASSED\n"
)


# ------------------------------------------------------------
# Reproduce exact clinical Cox model
# ------------------------------------------------------------

luad_dca_clinical_cox <-
  survival::coxph(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      age_at_diagnosis_years +
      sex_at_birth +
      stage_group,
    data =
      luad_dca_complete,
    x =
      TRUE,
    y =
      TRUE
  )


# ------------------------------------------------------------
# Reproduce exact clinical + MPCDS Cox model
# ------------------------------------------------------------

luad_dca_clinical_mpcds_cox <-
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
      luad_dca_complete,
    x =
      TRUE,
    y =
      TRUE
  )


cat(
  "\n========================================\n",
  "REPRODUCED CLINICAL COX MODEL\n",
  "========================================\n"
)


print(
  summary(
    luad_dca_clinical_cox
  )
)


cat(
  "\n========================================\n",
  "REPRODUCED CLINICAL + MPCDS COX MODEL\n",
  "========================================\n"
)


print(
  summary(
    luad_dca_clinical_mpcds_cox
  )
)


cat(
  "\n========================================\n",
  "LIKELIHOOD-RATIO MODEL COMPARISON\n",
  "========================================\n"
)


print(
  anova(
    luad_dca_clinical_cox,
    luad_dca_clinical_mpcds_cox,
    test =
      "LRT"
  )
)


dca_model_reproduction_summary <-
  data.frame(
    model =
      c(
        "Clinical",
        "Clinical + MPCDS"
      ),
    C_index =
      c(
        as.numeric(
          summary(
            luad_dca_clinical_cox
          )$concordance[
            1
          ]
        ),
        as.numeric(
          summary(
            luad_dca_clinical_mpcds_cox
          )$concordance[
            1
          ]
        )
      ),
    stringsAsFactors =
      FALSE
  )


cat(
  "\n========================================\n",
  "MODEL REPRODUCTION C-INDEX CHECK\n",
  "========================================\n"
)


print(
  dca_model_reproduction_summary,
  row.names = FALSE
)


stopifnot(
  abs(
    dca_model_reproduction_summary$C_index[
      dca_model_reproduction_summary$model ==
        "Clinical"
    ] -
      0.6740567
  ) <
    0.001
)


stopifnot(
  abs(
    dca_model_reproduction_summary$C_index[
      dca_model_reproduction_summary$model ==
        "Clinical + MPCDS"
    ] -
      0.7243107
  ) <
    0.001
)


cat(
  "\n========================================\n",
  "SECTION 2 MODEL REPRODUCTION SUCCESSFUL\n",
  "========================================\n"
)


# ============================================================
# 3. Assess 1-, 3-, and 5-year outcome incidence for DCA
# ============================================================


dca_times_days <-
  c(
    365,
    3 * 365,
    5 * 365
  )


dca_time_labels <-
  c(
    "1 year",
    "3 years",
    "5 years"
  )


luad_dca_km <-
  survival::survfit(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      1,
    data =
      luad_dca_complete
  )


luad_dca_km_summary <-
  summary(
    luad_dca_km,
    times =
      dca_times_days,
    extend =
      TRUE
  )


dca_horizon_summary <-
  data.frame(
    
    timepoint =
      dca_time_labels,
    
    time_days =
      dca_times_days,
    
    KM_survival =
      as.numeric(
        luad_dca_km_summary$surv
      ),
    
    KM_event_probability =
      1 -
      as.numeric(
        luad_dca_km_summary$surv
      ),
    
    events_observed_by_time =
      vapply(
        dca_times_days,
        function(
    horizon
        ) {
          
          sum(
            luad_dca_complete$OS_status ==
              1 &
              luad_dca_complete$OS_time <=
              horizon
          )
        },
    numeric(1)
      ),
    
    observed_beyond_time =
      vapply(
        dca_times_days,
        function(
    horizon
        ) {
          
          sum(
            luad_dca_complete$OS_time >
              horizon
          )
        },
    numeric(1)
      ),
    
    censored_by_time =
      vapply(
        dca_times_days,
        function(
    horizon
        ) {
          
          sum(
            luad_dca_complete$OS_status ==
              0 &
              luad_dca_complete$OS_time <=
              horizon
          )
        },
    numeric(1)
      ),
    
    total =
      nrow(
        luad_dca_complete
      ),
    
    stringsAsFactors =
      FALSE
  )


cat(
  "\n========================================\n",
  "DCA HORIZON OUTCOME SUMMARY\n",
  "========================================\n"
)


print(
  dca_horizon_summary,
  row.names = FALSE
)


dca_horizon_summary_display <-
  dca_horizon_summary


dca_horizon_summary_display$KM_survival_percent <-
  100 *
  dca_horizon_summary_display$KM_survival


dca_horizon_summary_display$KM_event_percent <-
  100 *
  dca_horizon_summary_display$KM_event_probability


cat(
  "\n========================================\n",
  "DCA EVENT PROBABILITY BY HORIZON\n",
  "========================================\n"
)


print(
  dca_horizon_summary_display[
    ,
    c(
      "timepoint",
      "KM_survival_percent",
      "KM_event_percent",
      "events_observed_by_time",
      "observed_beyond_time",
      "censored_by_time",
      "total"
    )
  ],
  row.names =
    FALSE
)


cat(
  "\nMaximum follow-up in years:\n"
)


print(
  max(
    luad_dca_complete$OS_time
  ) /
    365
)


stopifnot(
  length(
    luad_dca_km_summary$surv
  ) ==
    3
)


stopifnot(
  all(
    is.finite(
      dca_horizon_summary$KM_event_probability
    )
  )
)


stopifnot(
  all(
    dca_horizon_summary$KM_event_probability >
      0 &
      dca_horizon_summary$KM_event_probability <
      1
  )
)


cat(
  "\n========================================\n",
  "SECTION 3 DCA HORIZON ASSESSMENT SUCCESSFUL\n",
  "========================================\n"
)


# ============================================================
# 4. Generate model-predicted mortality probabilities for DCA
# ============================================================


dca_thresholds <-
  seq(
    from =
      0.05,
    to =
      0.70,
    by =
      0.01
  )


cat(
  "\n========================================\n",
  "DCA THRESHOLD PROBABILITIES\n",
  "========================================\n"
)


print(
  range(
    dca_thresholds
  )
)


cat(
  "\nNumber of thresholds:\n"
)


print(
  length(
    dca_thresholds
  )
)


# ------------------------------------------------------------
# Cox mortality-risk prediction helper
# ------------------------------------------------------------

predict_cox_mortality_risk <-
  function(
    model,
    newdata,
    horizon_days
  ) {
    
    baseline_hazard <-
      survival::basehaz(
        model,
        centered =
          FALSE
      )
    
    
    valid_hazard_rows <-
      which(
        baseline_hazard$time <=
          horizon_days
      )
    
    
    if (
      length(
        valid_hazard_rows
      ) ==
      0
    ) {
      
      stop(
        paste0(
          "No baseline hazard estimate available at or before ",
          horizon_days,
          " days."
        )
      )
    }
    
    
    horizon_row <-
      max(
        valid_hazard_rows
      )
    
    
    baseline_cumulative_hazard <-
      baseline_hazard$hazard[
        horizon_row
      ]
    
    
    linear_predictor <-
      as.numeric(
        stats::predict(
          model,
          newdata =
            newdata,
          type =
            "lp",
          reference =
            "zero"
        )
      )
    
    
    subject_cumulative_hazard <-
      baseline_cumulative_hazard *
      exp(
        linear_predictor
      )
    
    
    predicted_risk <-
      1 -
      exp(
        -subject_cumulative_hazard
      )
    
    
    predicted_risk[
      predicted_risk <
        0
    ] <-
      0
    
    
    predicted_risk[
      predicted_risk >
        1
    ] <-
      1
    
    
    predicted_risk
  }


# ------------------------------------------------------------
# Clinical model predicted risks
# ------------------------------------------------------------

luad_dca_complete$clinical_risk_1y <-
  predict_cox_mortality_risk(
    model =
      luad_dca_clinical_cox,
    newdata =
      luad_dca_complete,
    horizon_days =
      365
  )


luad_dca_complete$clinical_risk_3y <-
  predict_cox_mortality_risk(
    model =
      luad_dca_clinical_cox,
    newdata =
      luad_dca_complete,
    horizon_days =
      3 * 365
  )


luad_dca_complete$clinical_risk_5y <-
  predict_cox_mortality_risk(
    model =
      luad_dca_clinical_cox,
    newdata =
      luad_dca_complete,
    horizon_days =
      5 * 365
  )


# ------------------------------------------------------------
# Clinical + MPCDS predicted risks
# ------------------------------------------------------------

luad_dca_complete$clinical_mpcds_risk_1y <-
  predict_cox_mortality_risk(
    model =
      luad_dca_clinical_mpcds_cox,
    newdata =
      luad_dca_complete,
    horizon_days =
      365
  )


luad_dca_complete$clinical_mpcds_risk_3y <-
  predict_cox_mortality_risk(
    model =
      luad_dca_clinical_mpcds_cox,
    newdata =
      luad_dca_complete,
    horizon_days =
      3 * 365
  )


luad_dca_complete$clinical_mpcds_risk_5y <-
  predict_cox_mortality_risk(
    model =
      luad_dca_clinical_mpcds_cox,
    newdata =
      luad_dca_complete,
    horizon_days =
      5 * 365
  )


dca_risk_columns <-
  c(
    "clinical_risk_1y",
    "clinical_mpcds_risk_1y",
    "clinical_risk_3y",
    "clinical_mpcds_risk_3y",
    "clinical_risk_5y",
    "clinical_mpcds_risk_5y"
  )


stopifnot(
  all(
    dca_risk_columns %in%
      colnames(
        luad_dca_complete
      )
  )
)


stopifnot(
  all(
    vapply(
      luad_dca_complete[
        ,
        dca_risk_columns,
        drop = FALSE
      ],
      function(
    x
      ) {
        all(
          is.finite(
            x
          )
        )
      },
    logical(1)
    )
  )
)


stopifnot(
  all(
    vapply(
      luad_dca_complete[
        ,
        dca_risk_columns,
        drop = FALSE
      ],
      function(
    x
      ) {
        all(
          x >=
            0 &
            x <=
            1
        )
      },
    logical(1)
    )
  )
)


cat(
  "\n========================================\n",
  "PREDICTED MORTALITY-RISK SUMMARIES\n",
  "========================================\n"
)


print(
  summary(
    luad_dca_complete[
      ,
      dca_risk_columns,
      drop = FALSE
    ]
  )
)


dca_mean_predicted_risk <-
  data.frame(
    
    timepoint =
      c(
        "1 year",
        "1 year",
        "3 years",
        "3 years",
        "5 years",
        "5 years"
      ),
    
    model =
      c(
        "Clinical",
        "Clinical + MPCDS",
        "Clinical",
        "Clinical + MPCDS",
        "Clinical",
        "Clinical + MPCDS"
      ),
    
    mean_predicted_risk =
      c(
        mean(
          luad_dca_complete$clinical_risk_1y
        ),
        mean(
          luad_dca_complete$clinical_mpcds_risk_1y
        ),
        mean(
          luad_dca_complete$clinical_risk_3y
        ),
        mean(
          luad_dca_complete$clinical_mpcds_risk_3y
        ),
        mean(
          luad_dca_complete$clinical_risk_5y
        ),
        mean(
          luad_dca_complete$clinical_mpcds_risk_5y
        )
      ),
    
    stringsAsFactors =
      FALSE
  )


cat(
  "\n========================================\n",
  "MEAN PREDICTED MORTALITY RISKS\n",
  "========================================\n"
)


print(
  dca_mean_predicted_risk,
  row.names = FALSE
)


dca_calibration_sanity <-
  data.frame(
    
    timepoint =
      dca_time_labels,
    
    KM_event_probability =
      dca_horizon_summary$KM_event_probability,
    
    clinical_mean_risk =
      c(
        mean(
          luad_dca_complete$clinical_risk_1y
        ),
        mean(
          luad_dca_complete$clinical_risk_3y
        ),
        mean(
          luad_dca_complete$clinical_risk_5y
        )
      ),
    
    clinical_mpcds_mean_risk =
      c(
        mean(
          luad_dca_complete$clinical_mpcds_risk_1y
        ),
        mean(
          luad_dca_complete$clinical_mpcds_risk_3y
        ),
        mean(
          luad_dca_complete$clinical_mpcds_risk_5y
        )
      ),
    
    stringsAsFactors =
      FALSE
  )


cat(
  "\n========================================\n",
  "DCA CALIBRATION SANITY CHECK\n",
  "========================================\n"
)


print(
  dca_calibration_sanity,
  row.names = FALSE
)


dca_risk_correlations <-
  data.frame(
    
    timepoint =
      dca_time_labels,
    
    correlation =
      c(
        stats::cor(
          luad_dca_complete$clinical_risk_1y,
          luad_dca_complete$clinical_mpcds_risk_1y,
          method =
            "spearman"
        ),
        stats::cor(
          luad_dca_complete$clinical_risk_3y,
          luad_dca_complete$clinical_mpcds_risk_3y,
          method =
            "spearman"
        ),
        stats::cor(
          luad_dca_complete$clinical_risk_5y,
          luad_dca_complete$clinical_mpcds_risk_5y,
          method =
            "spearman"
        )
      ),
    
    stringsAsFactors =
      FALSE
  )


cat(
  "\n========================================\n",
  "CLINICAL VS CLINICAL + MPCDS RISK CORRELATION\n",
  "========================================\n"
)


print(
  dca_risk_correlations,
  row.names = FALSE
)


cat(
  "\n========================================\n",
  "SECTION 4 PREDICTED RISKS GENERATED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# 5. Run survival decision curve analysis
# ============================================================


if (
  !requireNamespace(
    "dcurves",
    quietly = TRUE
  )
) {
  stop(
    "Package 'dcurves' is required. Install it with ",
    "install.packages('dcurves') and rerun Script 23."
  )
}

if (
  !requireNamespace(
    "tibble",
    quietly = TRUE
  )
) {
  stop(
    "Package 'tibble' is required. Install it with ",
    "install.packages('tibble') and rerun Script 23."
  )
}


cat(
  "\n========================================\n",
  "DECISION CURVE ANALYSIS PACKAGE\n",
  "========================================\n"
)


cat(
  "\ndcurves version:\n"
)


print(
  as.character(
    utils::packageVersion(
      "dcurves"
    )
  )
)


# ------------------------------------------------------------
# Compact DCA input dataset
# ------------------------------------------------------------

luad_dca_analysis <-
  luad_dca_complete[
    ,
    c(
      "patient_id",
      "OS_time",
      "OS_status",
      "clinical_risk_1y",
      "clinical_mpcds_risk_1y",
      "clinical_risk_3y",
      "clinical_mpcds_risk_3y",
      "clinical_risk_5y",
      "clinical_mpcds_risk_5y"
    ),
    drop = FALSE
  ]


stopifnot(
  nrow(
    luad_dca_analysis
  ) ==
    475
)


stopifnot(
  sum(
    luad_dca_analysis$OS_status ==
      1
  ) ==
    171
)


# ------------------------------------------------------------
# 1-year survival DCA
# ------------------------------------------------------------

luad_dca_1y <-
  dcurves::dca(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      clinical_risk_1y +
      clinical_mpcds_risk_1y,
    data =
      luad_dca_analysis,
    thresholds =
      dca_thresholds,
    time =
      365,
    label =
      list(
        clinical_risk_1y =
          "Clinical",
        clinical_mpcds_risk_1y =
          "Clinical + MPCDS"
      )
  )


# ------------------------------------------------------------
# 3-year survival DCA
# ------------------------------------------------------------

luad_dca_3y <-
  dcurves::dca(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      clinical_risk_3y +
      clinical_mpcds_risk_3y,
    data =
      luad_dca_analysis,
    thresholds =
      dca_thresholds,
    time =
      3 * 365,
    label =
      list(
        clinical_risk_3y =
          "Clinical",
        clinical_mpcds_risk_3y =
          "Clinical + MPCDS"
      )
  )


# ------------------------------------------------------------
# 5-year survival DCA
# ------------------------------------------------------------

luad_dca_5y <-
  dcurves::dca(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      clinical_risk_5y +
      clinical_mpcds_risk_5y,
    data =
      luad_dca_analysis,
    thresholds =
      dca_thresholds,
    time =
      5 * 365,
    label =
      list(
        clinical_risk_5y =
          "Clinical",
        clinical_mpcds_risk_5y =
          "Clinical + MPCDS"
      )
  )


# ------------------------------------------------------------
# Correct exported tibble conversion
# ------------------------------------------------------------

luad_dca_1y_table <-
  tibble::as_tibble(
    luad_dca_1y
  )


luad_dca_3y_table <-
  tibble::as_tibble(
    luad_dca_3y
  )


luad_dca_5y_table <-
  tibble::as_tibble(
    luad_dca_5y
  )


cat(
  "\n========================================\n",
  "DCA TABLE COLUMN NAMES\n",
  "========================================\n"
)


print(
  colnames(
    luad_dca_1y_table
  )
)


cat(
  "\n1-year model labels:\n"
)


print(
  unique(
    luad_dca_1y_table$label
  )
)


stopifnot(
  inherits(
    luad_dca_1y,
    "dca"
  ),
  inherits(
    luad_dca_3y,
    "dca"
  ),
  inherits(
    luad_dca_5y,
    "dca"
  )
)


stopifnot(
  nrow(
    luad_dca_1y_table
  ) >
    0,
  nrow(
    luad_dca_3y_table
  ) >
    0,
  nrow(
    luad_dca_5y_table
  ) >
    0
)


cat(
  "\n========================================\n",
  "SECTION 5 SURVIVAL DCA COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)


# ============================================================
# 6. Quantify incremental net benefit of adding MPCDS
# ============================================================


build_dca_comparison <-
  function(
    dca_table,
    timepoint
  ) {
    
    dca_df <-
      as.data.frame(
        dca_table
      )
    
    
    dca_df <-
      dca_df[
        ,
        c(
          "label",
          "threshold",
          "net_benefit"
        ),
        drop = FALSE
      ]
    
    
    dca_df$label <-
      as.character(
        dca_df$label
      )
    
    
    treat_all <-
      dca_df[
        dca_df$label ==
          "Treat All",
        c(
          "threshold",
          "net_benefit"
        ),
        drop = FALSE
      ]
    
    
    names(
      treat_all
    )[
      names(
        treat_all
      ) ==
        "net_benefit"
    ] <-
      "NB_treat_all"
    
    
    treat_none <-
      dca_df[
        dca_df$label ==
          "Treat None",
        c(
          "threshold",
          "net_benefit"
        ),
        drop = FALSE
      ]
    
    
    names(
      treat_none
    )[
      names(
        treat_none
      ) ==
        "net_benefit"
    ] <-
      "NB_treat_none"
    
    
    clinical <-
      dca_df[
        dca_df$label ==
          "Clinical",
        c(
          "threshold",
          "net_benefit"
        ),
        drop = FALSE
      ]
    
    
    names(
      clinical
    )[
      names(
        clinical
      ) ==
        "net_benefit"
    ] <-
      "NB_clinical"
    
    
    clinical_mpcds <-
      dca_df[
        dca_df$label ==
          "Clinical + MPCDS",
        c(
          "threshold",
          "net_benefit"
        ),
        drop = FALSE
      ]
    
    
    names(
      clinical_mpcds
    )[
      names(
        clinical_mpcds
      ) ==
        "net_benefit"
    ] <-
      "NB_clinical_mpcds"
    
    
    comparison <-
      Reduce(
        function(
    x,
    y
        ) {
          
          merge(
            x,
            y,
            by =
              "threshold",
            all =
              TRUE,
            sort =
              TRUE
          )
        },
    list(
      treat_all,
      treat_none,
      clinical,
      clinical_mpcds
    )
      )
    
    
    comparison$timepoint <-
      timepoint
    
    
    comparison$NB_best_default <-
      pmax(
        comparison$NB_treat_all,
        comparison$NB_treat_none,
        na.rm =
          TRUE
      )
    
    
    comparison$delta_NB_mpcds_vs_clinical <-
      comparison$NB_clinical_mpcds -
      comparison$NB_clinical
    
    
    comparison$delta_NB_mpcds_vs_default <-
      comparison$NB_clinical_mpcds -
      comparison$NB_best_default
    
    
    comparison$clinical_finite <-
      is.finite(
        comparison$NB_clinical
      )
    
    
    comparison$clinical_mpcds_finite <-
      is.finite(
        comparison$NB_clinical_mpcds
      )
    
    
    comparison$both_models_finite <-
      comparison$clinical_finite &
      comparison$clinical_mpcds_finite
    
    
    comparison$mpcds_better_than_clinical <-
      comparison$both_models_finite &
      comparison$NB_clinical_mpcds >
      comparison$NB_clinical
    
    
    comparison$mpcds_better_than_default <-
      comparison$clinical_mpcds_finite &
      comparison$NB_clinical_mpcds >
      comparison$NB_best_default
    
    
    comparison$incremental_clinical_utility <-
      comparison$mpcds_better_than_clinical &
      comparison$mpcds_better_than_default
    
    
    comparison <-
      comparison[
        ,
        c(
          "timepoint",
          "threshold",
          "NB_treat_all",
          "NB_treat_none",
          "NB_best_default",
          "NB_clinical",
          "NB_clinical_mpcds",
          "delta_NB_mpcds_vs_clinical",
          "delta_NB_mpcds_vs_default",
          "clinical_finite",
          "clinical_mpcds_finite",
          "both_models_finite",
          "mpcds_better_than_clinical",
          "mpcds_better_than_default",
          "incremental_clinical_utility"
        )
      ]
    
    
    comparison
  }


luad_dca_comparison_1y <-
  build_dca_comparison(
    luad_dca_1y_table,
    "1 year"
  )


luad_dca_comparison_3y <-
  build_dca_comparison(
    luad_dca_3y_table,
    "3 years"
  )


luad_dca_comparison_5y <-
  build_dca_comparison(
    luad_dca_5y_table,
    "5 years"
  )


luad_dca_comparison_all <-
  rbind(
    luad_dca_comparison_1y,
    luad_dca_comparison_3y,
    luad_dca_comparison_5y
  )


rownames(
  luad_dca_comparison_all
) <-
  NULL


stopifnot(
  nrow(
    luad_dca_comparison_1y
  ) ==
    length(
      dca_thresholds
    ),
  nrow(
    luad_dca_comparison_3y
  ) ==
    length(
      dca_thresholds
    ),
  nrow(
    luad_dca_comparison_5y
  ) ==
    length(
      dca_thresholds
    )
)


# ------------------------------------------------------------
# Non-finite estimate summary
# ------------------------------------------------------------

dca_nonfinite_summary <-
  data.frame(
    
    timepoint =
      dca_time_labels,
    
    n_thresholds =
      c(
        nrow(
          luad_dca_comparison_1y
        ),
        nrow(
          luad_dca_comparison_3y
        ),
        nrow(
          luad_dca_comparison_5y
        )
      ),
    
    n_clinical_finite =
      c(
        sum(
          luad_dca_comparison_1y$clinical_finite
        ),
        sum(
          luad_dca_comparison_3y$clinical_finite
        ),
        sum(
          luad_dca_comparison_5y$clinical_finite
        )
      ),
    
    n_clinical_mpcds_finite =
      c(
        sum(
          luad_dca_comparison_1y$clinical_mpcds_finite
        ),
        sum(
          luad_dca_comparison_3y$clinical_mpcds_finite
        ),
        sum(
          luad_dca_comparison_5y$clinical_mpcds_finite
        )
      ),
    
    n_both_models_finite =
      c(
        sum(
          luad_dca_comparison_1y$both_models_finite
        ),
        sum(
          luad_dca_comparison_3y$both_models_finite
        ),
        sum(
          luad_dca_comparison_5y$both_models_finite
        )
      ),
    
    stringsAsFactors =
      FALSE
  )


cat(
  "\n========================================\n",
  "NON-FINITE DCA ESTIMATES\n",
  "========================================\n"
)


print(
  dca_nonfinite_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# Horizon summary helper
# ------------------------------------------------------------

summarize_dca_horizon <-
  function(
    comparison
  ) {
    
    valid <-
      comparison[
        comparison$both_models_finite,
        ,
        drop = FALSE
      ]
    
    
    better_clinical <-
      valid[
        valid$mpcds_better_than_clinical,
        ,
        drop = FALSE
      ]
    
    
    useful <-
      valid[
        valid$incremental_clinical_utility,
        ,
        drop = FALSE
      ]
    
    
    if (
      nrow(
        valid
      ) >
      0
    ) {
      
      max_delta_index <-
        which.max(
          valid$delta_NB_mpcds_vs_clinical
        )
      
      
      maximum_delta <-
        valid$delta_NB_mpcds_vs_clinical[
          max_delta_index
        ]
      
      
      maximum_delta_threshold <-
        valid$threshold[
          max_delta_index
        ]
      
    } else {
      
      maximum_delta <-
        NA_real_
      
      maximum_delta_threshold <-
        NA_real_
    }
    
    
    if (
      nrow(
        better_clinical
      ) >
      0
    ) {
      
      min_better_clinical <-
        min(
          better_clinical$threshold
        )
      
      max_better_clinical <-
        max(
          better_clinical$threshold
        )
      
    } else {
      
      min_better_clinical <-
        NA_real_
      
      max_better_clinical <-
        NA_real_
    }
    
    
    if (
      nrow(
        useful
      ) >
      0
    ) {
      
      min_useful <-
        min(
          useful$threshold
        )
      
      max_useful <-
        max(
          useful$threshold
        )
      
    } else {
      
      min_useful <-
        NA_real_
      
      max_useful <-
        NA_real_
    }
    
    
    data.frame(
      
      timepoint =
        unique(
          comparison$timepoint
        )[1],
      
      thresholds_total =
        nrow(
          comparison
        ),
      
      thresholds_both_models_finite =
        nrow(
          valid
        ),
      
      thresholds_mpcds_better_than_clinical =
        nrow(
          better_clinical
        ),
      
      proportion_mpcds_better_than_clinical =
        ifelse(
          nrow(
            valid
          ) >
            0,
          nrow(
            better_clinical
          ) /
            nrow(
              valid
            ),
          NA_real_
        ),
      
      min_threshold_mpcds_better_than_clinical =
        min_better_clinical,
      
      max_threshold_mpcds_better_than_clinical =
        max_better_clinical,
      
      thresholds_incremental_clinical_utility =
        nrow(
          useful
        ),
      
      proportion_incremental_clinical_utility =
        ifelse(
          nrow(
            valid
          ) >
            0,
          nrow(
            useful
          ) /
            nrow(
              valid
            ),
          NA_real_
        ),
      
      min_threshold_incremental_clinical_utility =
        min_useful,
      
      max_threshold_incremental_clinical_utility =
        max_useful,
      
      maximum_delta_NB_vs_clinical =
        maximum_delta,
      
      threshold_at_maximum_delta =
        maximum_delta_threshold,
      
      stringsAsFactors =
        FALSE
    )
  }


luad_dca_incremental_summary <-
  rbind(
    
    summarize_dca_horizon(
      luad_dca_comparison_1y
    ),
    
    summarize_dca_horizon(
      luad_dca_comparison_3y
    ),
    
    summarize_dca_horizon(
      luad_dca_comparison_5y
    )
  )


rownames(
  luad_dca_incremental_summary
) <-
  NULL


cat(
  "\n========================================\n",
  "INCREMENTAL DCA SUMMARY\n",
  "========================================\n"
)


print(
  luad_dca_incremental_summary,
  row.names = FALSE
)


cat(
  "\n========================================\n",
  "SECTION 6 INCREMENTAL NET BENEFIT ANALYSIS SUCCESSFUL\n",
  "========================================\n"
)


# ============================================================
# 7. Save DCA analysis results
# ============================================================


dca_result_dir <-
  "04_results/mpcds/decision_curve_analysis"



dir.create(
  dca_result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Output paths
# ------------------------------------------------------------

dca_patient_risk_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_patient_predicted_risks.csv"
  )


dca_horizon_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_horizon_summary.csv"
  )


dca_calibration_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_calibration_sanity_check.csv"
  )


dca_raw_1y_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_1year_raw.csv"
  )


dca_raw_3y_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_3year_raw.csv"
  )


dca_raw_5y_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_5year_raw.csv"
  )


dca_comparison_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_threshold_level_comparison.csv"
  )


dca_incremental_summary_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_incremental_net_benefit_summary.csv"
  )


dca_nonfinite_file <-
  file.path(
    dca_result_dir,
    "TCGA_LUAD_DCA_nonfinite_estimate_summary.csv"
  )


# ------------------------------------------------------------
# Save tables
# ------------------------------------------------------------

write.csv(
  luad_dca_analysis,
  dca_patient_risk_file,
  row.names = FALSE
)


write.csv(
  dca_horizon_summary,
  dca_horizon_file,
  row.names = FALSE
)


write.csv(
  dca_calibration_sanity,
  dca_calibration_file,
  row.names = FALSE
)


write.csv(
  as.data.frame(
    luad_dca_1y_table
  ),
  dca_raw_1y_file,
  row.names = FALSE
)


write.csv(
  as.data.frame(
    luad_dca_3y_table
  ),
  dca_raw_3y_file,
  row.names = FALSE
)


write.csv(
  as.data.frame(
    luad_dca_5y_table
  ),
  dca_raw_5y_file,
  row.names = FALSE
)


write.csv(
  luad_dca_comparison_all,
  dca_comparison_file,
  row.names = FALSE
)


write.csv(
  luad_dca_incremental_summary,
  dca_incremental_summary_file,
  row.names = FALSE
)


write.csv(
  dca_nonfinite_summary,
  dca_nonfinite_file,
  row.names = FALSE
)


saved_dca_tables <-
  c(
    dca_patient_risk_file,
    dca_horizon_file,
    dca_calibration_file,
    dca_raw_1y_file,
    dca_raw_3y_file,
    dca_raw_5y_file,
    dca_comparison_file,
    dca_incremental_summary_file,
    dca_nonfinite_file
  )


cat(
  "\n========================================\n",
  "SAVED DCA RESULT TABLES\n",
  "========================================\n"
)


print(
  data.frame(
    file =
      saved_dca_tables,
    exists =
      file.exists(
        saved_dca_tables
      ),
    stringsAsFactors =
      FALSE
  )
)


stopifnot(
  all(
    file.exists(
      saved_dca_tables
    )
  )
)


# ------------------------------------------------------------
# Figure generation is intentionally deferred to Script 32.
# Script 23 produces analysis tables only.
# ------------------------------------------------------------


# ------------------------------------------------------------
# Final summaries
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "FINAL DCA HORIZON SUMMARY\n",
  "========================================\n"
)


print(
  dca_horizon_summary,
  row.names = FALSE
)


cat(
  "\n========================================\n",
  "FINAL INCREMENTAL NET-BENEFIT SUMMARY\n",
  "========================================\n"
)


print(
  luad_dca_incremental_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# Scientific interpretation
#
# - Continuous MPCDS was used; no cutoff was optimized.
#
# - Net benefit was evaluated from 5% to 70% threshold risk.
#
# - Maximum apparent incremental NB versus Clinical:
#
#     1 year:
#       threshold ~0.23
#       delta NB ~0.0175
#
#     3 years:
#       threshold ~0.64
#       delta NB ~0.0891
#
#     5 years:
#       threshold ~0.52
#       delta NB ~0.0747
#
# - Minimum/maximum useful thresholds do NOT imply benefit
#   at every intervening threshold. Use the threshold-level
#   comparison table for exact interpretation.
#
# - These are apparent TCGA-LUAD development-cohort estimates,
#   not externally validated evidence of clinical utility.
#
# - The 5-year analysis requires particular caution because
#   censoring before 5 years is substantial.
# ------------------------------------------------------------


cat(
  "\n============================================================\n",
  "SCRIPT 23 LUAD MPCDS DECISION CURVE ANALYSIS SUMMARY\n",
  "============================================================\n",
  sep = ""
)

cat(
  "Analysis cohort: ",
  nrow(luad_dca_complete),
  " patients; ",
  sum(luad_dca_complete$OS_status),
  " events\n",
  "Clinical C-index: ",
  round(
    dca_model_reproduction_summary$C_index[
      dca_model_reproduction_summary$model == "Clinical"
    ],
    4
  ),
  "\n",
  "Clinical + MPCDS C-index: ",
  round(
    dca_model_reproduction_summary$C_index[
      dca_model_reproduction_summary$model == "Clinical + MPCDS"
    ],
    4
  ),
  "\n",
  "Threshold range: 0.05 to 0.70\n",
  "Outputs: ",
  dca_result_dir,
  "\n",
  sep = ""
)


cat(
  "\n========================================\n",
  "SCRIPT 23 COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)