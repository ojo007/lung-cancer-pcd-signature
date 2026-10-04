# ============================================================
# Script 21: Validate LUAD MPCDS across external GEO cohorts
#
# Purpose:
#   1. Apply the fixed TCGA-LUAD 47-gene ridge-Cox MPCDS to
#      independent LUAD GEO cohorts prepared by Script 19.
#   2. Preserve the TCGA-derived coefficients; no outcome-based
#      refitting or feature selection is performed externally.
#   3. Adapt to cross-platform expression by applying the
#      TCGA coefficient weights to within-cohort gene z-scores,
#      following the established GSE68465 validation strategy.
#   4. Evaluate continuous MPCDS using Harrell's C-index and
#      Cox regression.
#   5. Use median high/low groups only for descriptive KM
#      visualization; no outcome-optimized cutoff is derived.
#   6. Retain the established clinically adjusted GSE68465
#      sensitivity analysis.
#
# External LUAD cohorts:
#   GSE68465
#   GSE72094
#   GSE31210
#   GSE50081
#
# Important:
#   - This is external validation of the LUAD prognostic MPCDS.
#   - No LUSC MPCDS is advanced because Script 18 showed
#     near-chance internal discrimination in LUSC.
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "survival"
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


# ------------------------------------------------------------
# 2. Canonical directories
# ------------------------------------------------------------

geo_prepared_dir <- file.path(
  processed_dir,
  "GEO",
  "external_validation"
)

mpcds_dir <- file.path(
  results_dir,
  "mpcds"
)

external_mpcds_dir <- file.path(
  mpcds_dir,
  "external_validation"
)

survival_dir <- file.path(
  processed_dir,
  "survival"
)

ensure_dir(
  external_mpcds_dir
)


# ------------------------------------------------------------
# 3. Required TCGA model inputs
# ------------------------------------------------------------

tcga_expression_file <- file.path(
  survival_dir,
  "TCGA_LUAD_PCD_survival_expression.rds"
)

tcga_scores_file <- file.path(
  mpcds_dir,
  "TCGA_LUAD_MPCDS_patient_scores.csv"
)

tcga_coefficients_file <- file.path(
  mpcds_dir,
  "TCGA_LUAD_MPCDS_ridge_coefficients.csv"
)

manifest_file <- file.path(
  geo_prepared_dir,
  "external_GEO_cohort_manifest.csv"
)

check_files_exist(
  c(
    tcga_expression_file,
    tcga_scores_file,
    tcga_coefficients_file,
    manifest_file
  )
)


# ------------------------------------------------------------
# 4. Load fixed TCGA-LUAD MPCDS
# ------------------------------------------------------------

tcga_expression <- readRDS(
  tcga_expression_file
)

tcga_scores <- read.csv(
  tcga_scores_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tcga_coefficients <- read.csv(
  tcga_coefficients_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(tcga_coefficients) == 47,
  all(
    c(
      "gene",
      "coefficient"
    ) %in%
      colnames(
        tcga_coefficients
      )
  ),
  all(
    tcga_scores$patient_id %in%
      colnames(
        tcga_expression
      )
  )
)

tcga_x <- t(
  tcga_expression[
    ,
    tcga_scores$patient_id,
    drop = FALSE
  ]
)

stopifnot(
  nrow(tcga_x) == 504,
  ncol(tcga_x) == 47,
  identical(
    rownames(tcga_x),
    tcga_scores$patient_id
  )
)

beta_47 <- tcga_coefficients$coefficient[
  match(
    colnames(tcga_x),
    tcga_coefficients$gene
  )
]

names(beta_47) <- colnames(
  tcga_x
)

stopifnot(
  !anyNA(
    beta_47
  )
)


# ------------------------------------------------------------
# 5. Verify exact reconstruction of TCGA MPCDS
# ------------------------------------------------------------

tcga_manual_score <- as.numeric(
  tcga_x %*%
    beta_47
)

tcga_score_correlation <- stats::cor(
  tcga_manual_score,
  tcga_scores$MPCDS,
  method = "pearson"
)

stopifnot(
  is.finite(
    tcga_score_correlation
  ),
  tcga_score_correlation >
    0.999999
)


# ------------------------------------------------------------
# 6. TCGA standardized-expression weights
#
# glmnet returns coefficients on the original predictor scale.
# For transfer to independently standardized external data:
#
#   standardized weight = raw beta * TCGA gene SD
#
# This reproduces the established GSE68465 validation method.
# ------------------------------------------------------------

tcga_gene_sd <- apply(
  tcga_x,
  2,
  stats::sd
)

stopifnot(
  !anyNA(
    tcga_gene_sd
  ),
  all(
    tcga_gene_sd >
      0
  )
)

standardized_weights_47 <- beta_47 *
  tcga_gene_sd[
    names(
      beta_47
    )
  ]

tcga_x_z <- scale(
  tcga_x
)

tcga_score_zweight <- as.numeric(
  tcga_x_z %*%
    standardized_weights_47
)

stopifnot(
  stats::cor(
    tcga_score_zweight,
    tcga_manual_score,
    method = "pearson"
  ) >
    0.999999
)

write.csv(
  data.frame(
    gene =
      names(
        beta_47
      ),
    TCGA_raw_coefficient =
      as.numeric(
        beta_47
      ),
    TCGA_gene_SD =
      as.numeric(
        tcga_gene_sd[
          names(
            beta_47
          )
        ]
      ),
    standardized_weight =
      as.numeric(
        standardized_weights_47
      ),
    stringsAsFactors =
      FALSE
  ),
  file.path(
    external_mpcds_dir,
    "TCGA_LUAD_MPCDS_standardized_weights_47.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. LUAD external cohort manifest
# ------------------------------------------------------------

all_manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

luad_manifest <- all_manifest |>
  dplyr::filter(
    histology ==
      "LUAD"
  )

stopifnot(
  nrow(
    luad_manifest
  ) == 4,
  setequal(
    luad_manifest$accession,
    c(
      "GSE68465",
      "GSE72094",
      "GSE31210",
      "GSE50081"
    )
  )
)


# ------------------------------------------------------------
# 8. General helpers for phenotype columns
# ------------------------------------------------------------

normalize_text <- function(
    x
) {
  
  x <- trimws(
    as.character(
      x
    )
  )
  
  x[
    x %in%
      c(
        "",
        "NA",
        "N/A",
        "na",
        "n/a",
        "--",
        "unknown",
        "Unknown"
      )
  ] <- NA_character_
  
  x
}


find_first_column <- function(
    data,
    exact_names = character(),
    regex_patterns = character()
) {
  
  nms <- colnames(
    data
  )
  
  for (
    nm in exact_names
  ) {
    hit <- which(
      tolower(
        nms
      ) ==
        tolower(
          nm
        )
    )
    
    if (
      length(
        hit
      ) >
      0
    ) {
      return(
        nms[
          hit[
            1
          ]
        ]
      )
    }
  }
  
  for (
    pattern in regex_patterns
  ) {
    
    hit <- grep(
      pattern,
      nms,
      ignore.case = TRUE,
      perl = TRUE
    )
    
    if (
      length(
        hit
      ) >
      0
    ) {
      return(
        nms[
          hit[
            1
          ]
        ]
      )
    }
  }
  
  NA_character_
}


numeric_from_text <- function(
    x
) {
  
  x <- normalize_text(
    x
  )
  
  suppressWarnings(
    as.numeric(
      gsub(
        "[^0-9eE+.-]",
        "",
        x
      )
    )
  )
}


parse_event_status <- function(
    x
) {
  
  raw <- tolower(
    normalize_text(
      x
    )
  )
  
  out <- rep(
    NA_real_,
    length(
      raw
    )
  )
  
  out[
    raw %in%
      c(
        "dead",
        "deceased",
        "death",
        "died",
        "1",
        "yes",
        "event"
      )
  ] <- 1
  
  out[
    raw %in%
      c(
        "alive",
        "living",
        "0",
        "no",
        "censored",
        "censor"
      )
  ] <- 0
  
  # Handle phrases such as "dead: 1", "alive (0)", etc.
  out[
    is.na(
      out
    ) &
      grepl(
        "dead|deceased|died",
        raw
      )
  ] <- 1
  
  out[
    is.na(
      out
    ) &
      grepl(
        "alive|living|censor",
        raw
      )
  ] <- 0
  
  out
}


# ------------------------------------------------------------
# 9. Survival-column audit
#
# Cohort-specific exact names are tried first. Generic patterns
# are fallback only. The chosen fields are written to disk so
# that survival-variable provenance is explicit.
# ------------------------------------------------------------

extract_external_os <- function(
    pheno,
    accession
) {
  
  if (
    !"geo_accession" %in%
    colnames(
      pheno
    )
  ) {
    stop(
      accession,
      ": phenotype table lacks geo_accession."
    )
  }
  
  # Known GSE68465 fields from the prior validated workflow.
  if (
    accession ==
    "GSE68465"
  ) {
    
    time_col <-
      "months_to_last_contact_or_death:ch1"
    
    event_col <-
      "vital_status:ch1"
    
    if (
      !all(
        c(
          time_col,
          event_col
        ) %in%
        colnames(
          pheno
        )
      )
    ) {
      stop(
        accession,
        ": expected validated OS columns were not found."
      )
    }
    
    os_time <- numeric_from_text(
      pheno[[time_col]]
    )
    
    os_status <- parse_event_status(
      pheno[[event_col]]
    )
    
    time_unit <- "months"
    time_strategy <-
      "single_OS_time_column"
    
  } else if (
    accession ==
    "GSE31210"
  ) {
    
    time_col <-
      "days before death/censor:ch1"
    
    event_col <-
      "death:ch1"
    
    if (
      !all(
        c(
          time_col,
          event_col
        ) %in%
        colnames(
          pheno
        )
      )
    ) {
      stop(
        accession,
        ": expected OS columns were not found."
      )
    }
    
    os_time <- numeric_from_text(
      pheno[[time_col]]
    )
    
    os_status <- parse_event_status(
      pheno[[event_col]]
    )
    
    time_unit <- "days"
    time_strategy <-
      "GSE31210_days_before_death_or_censor"
    
  } else if (
    accession ==
    "GSE50081"
  ) {
    
    time_col <-
      "survival time:ch1"
    
    event_col <-
      "status:ch1"
    
    if (
      !all(
        c(
          time_col,
          event_col
        ) %in%
        colnames(
          pheno
        )
      )
    ) {
      stop(
        accession,
        ": expected OS columns were not found."
      )
    }
    
    os_time <- numeric_from_text(
      pheno[[time_col]]
    )
    
    os_status <- parse_event_status(
      pheno[[event_col]]
    )
    
    time_unit <- "years"
    time_strategy <-
      "GSE50081_survival_time_years"
    
  } else {
    
    event_col <- find_first_column(
      pheno,
      exact_names =
        c(
          "vital_status:ch1",
          "vital status:ch1",
          "os_status:ch1",
          "os status:ch1",
          "survival_status:ch1",
          "survival status:ch1",
          "death_event:ch1",
          "death event:ch1"
        ),
      regex_patterns =
        c(
          "vital.*status",
          "(^|[^a-z])os.*status",
          "survival.*status",
          "death.*event"
        )
    )
    
    # First preference: one explicit OS time/follow-up column.
    months_col <- find_first_column(
      pheno,
      exact_names =
        c(
          "overall_survival_months:ch1",
          "overall survival months:ch1",
          "os_months:ch1",
          "os months:ch1",
          "survival_months:ch1",
          "survival months:ch1",
          "months_to_last_contact_or_death:ch1"
        ),
      regex_patterns =
        c(
          "(overall.*survival|\\bos\\b).*month",
          "survival.*month",
          "month.*(survival|follow|death)"
        )
    )
    
    days_col <- find_first_column(
      pheno,
      exact_names =
        c(
          "overall_survival_days:ch1",
          "overall survival days:ch1",
          "os_days:ch1",
          "os days:ch1",
          "survival_days:ch1",
          "survival days:ch1"
        ),
      regex_patterns =
        c(
          "(overall.*survival|\\bos\\b).*day",
          "survival.*day"
        )
    )
    
    # Alternative: separate day-to-death and day-to-follow-up.
    death_days_col <- find_first_column(
      pheno,
      exact_names =
        c(
          "days_to_death:ch1",
          "days to death:ch1"
        ),
      regex_patterns =
        c(
          "day.*to.*death",
          "days.*death"
        )
    )
    
    follow_days_col <- find_first_column(
      pheno,
      exact_names =
        c(
          "days_to_last_followup:ch1",
          "days_to_last_follow_up:ch1",
          "days to last followup:ch1",
          "days to last follow up:ch1"
        ),
      regex_patterns =
        c(
          "day.*last.*follow",
          "days.*follow"
        )
    )
    
    # Generic fallback for OS/follow-up columns.
    generic_time_col <- find_first_column(
      pheno,
      regex_patterns =
        c(
          "(^|[^a-z])os([^a-z]|$)",
          "overall.*survival",
          "follow.*up",
          "survival.*time"
        )
    )
    
    if (
      is.na(
        event_col
      )
    ) {
      stop(
        accession,
        ": could not identify an OS event/vital-status column. ",
        "Inspect phenotype fields containing status/death/alive."
      )
    }
    
    os_status <- parse_event_status(
      pheno[[event_col]]
    )
    
    if (
      !is.na(
        months_col
      )
    ) {
      
      os_time <- numeric_from_text(
        pheno[[months_col]]
      )
      
      time_col <- months_col
      time_unit <- "months"
      time_strategy <-
        "single_OS_months_column"
      
    } else if (
      !is.na(
        days_col
      )
    ) {
      
      os_time <- numeric_from_text(
        pheno[[days_col]]
      )
      
      time_col <- days_col
      time_unit <- "days"
      time_strategy <-
        "single_OS_days_column"
      
    } else if (
      !is.na(
        death_days_col
      ) &&
      !is.na(
        follow_days_col
      )
    ) {
      
      death_days <- numeric_from_text(
        pheno[[death_days_col]]
      )
      
      follow_days <- numeric_from_text(
        pheno[[follow_days_col]]
      )
      
      os_time <- ifelse(
        os_status ==
          1,
        death_days,
        follow_days
      )
      
      # If the status-specific field is absent but the other
      # field is available, retain the available follow-up.
      missing_time <- is.na(
        os_time
      )
      
      os_time[
        missing_time
      ] <- dplyr::coalesce(
        death_days[
          missing_time
        ],
        follow_days[
          missing_time
        ]
      )
      
      time_col <- paste(
        death_days_col,
        follow_days_col,
        sep = " + "
      )
      
      time_unit <- "days"
      time_strategy <-
        "status_specific_death_followup_days"
      
    } else if (
      !is.na(
        generic_time_col
      )
    ) {
      
      os_time <- numeric_from_text(
        pheno[[generic_time_col]]
      )
      
      time_col <-
        generic_time_col
      
      # Determine likely units conservatively from the name.
      if (
        grepl(
          "month",
          generic_time_col,
          ignore.case = TRUE
        )
      ) {
        time_unit <- "months"
      } else if (
        grepl(
          "day",
          generic_time_col,
          ignore.case = TRUE
        )
      ) {
        time_unit <- "days"
      } else {
        
        # GEO LUAD survival datasets frequently encode OS in
        # months when the unit is omitted. Do not silently
        # guess: retain "unspecified" and use the numeric value
        # only for rank/Cox analyses, which are invariant to a
        # constant time-unit conversion.
        time_unit <- "unspecified"
      }
      
      time_strategy <-
        "generic_survival_time_column"
      
    } else {
      
      stop(
        accession,
        ": could not identify a usable OS/follow-up time column."
      )
    }
  }
  
  # Standardize survival time to months for downstream
  # cross-cohort comparability and later time-dependent ROC.
  if (
    identical(
      time_unit,
      "days"
    )
  ) {
    
    os_time_months <-
      os_time /
      30.4375
    
  } else if (
    identical(
      time_unit,
      "years"
    )
  ) {
    
    os_time_months <-
      os_time *
      12
    
  } else {
    
    os_time_months <-
      os_time
  }
  
  clinical <- data.frame(
    sample_id =
      pheno$geo_accession,
    OS_time =
      os_time_months,
    OS_status =
      os_status,
    stringsAsFactors =
      FALSE
  )
  
  clinical <- clinical |>
    dplyr::filter(
      !is.na(
        OS_time
      ),
      !is.na(
        OS_status
      ),
      OS_time >
        0,
      OS_status %in%
        c(
          0,
          1
        )
    )
  
  if (
    nrow(
      clinical
    ) <
    30
  ) {
    stop(
      accession,
      ": only ",
      nrow(
        clinical
      ),
      " samples had usable OS after parsing. ",
      "Review the selected survival fields."
    )
  }
  
  audit <- data.frame(
    accession =
      accession,
    time_column =
      time_col,
    event_column =
      event_col,
    source_time_unit =
      time_unit,
    analysis_time_unit =
      "months",
    time_strategy =
      time_strategy,
    n_phenotype_samples =
      nrow(
        pheno
      ),
    n_OS_usable =
      nrow(
        clinical
      ),
    n_events =
      sum(
        clinical$OS_status
      ),
    stringsAsFactors =
      FALSE
  )
  
  list(
    clinical =
      clinical,
    audit =
      audit
  )
}


# ------------------------------------------------------------
# 10. Generic MPCDS validation helper
# ------------------------------------------------------------

validate_mpcds_cohort <- function(
    accession
) {
  
  cat(
    "\n============================================================\n",
    "VALIDATING LUAD MPCDS: ",
    accession,
    "\n",
    "============================================================\n",
    sep = ""
  )
  
  cohort_dir <- file.path(
    geo_prepared_dir,
    accession
  )
  
  cohort_result_dir <- file.path(
    external_mpcds_dir,
    accession
  )
  
  ensure_dir(
    cohort_result_dir
  )
  
  z_file <- file.path(
    cohort_dir,
    paste0(
      accession,
      "_candidate47_expression_z.rds"
    )
  )
  
  pheno_file <- file.path(
    cohort_dir,
    paste0(
      accession,
      "_phenotype_selected.csv"
    )
  )
  
  check_files_exist(
    c(
      z_file,
      pheno_file
    )
  )
  
  external_z <- readRDS(
    z_file
  )
  
  pheno <- read.csv(
    pheno_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  available_genes <- intersect(
    names(
      standardized_weights_47
    ),
    rownames(
      external_z
    )
  )
  
  missing_genes <- setdiff(
    names(
      standardized_weights_47
    ),
    available_genes
  )
  
  if (
    length(
      available_genes
    ) <
    38
  ) {
    stop(
      accession,
      ": insufficient MPCDS gene coverage: ",
      length(
        available_genes
      ),
      "/47."
    )
  }
  
  external_z <- external_z[
    available_genes,
    ,
    drop = FALSE
  ]
  
  weights <- standardized_weights_47[
    available_genes
  ]
  
  stopifnot(
    identical(
      rownames(
        external_z
      ),
      names(
        weights
      )
    )
  )
  
  score_all <- as.numeric(
    t(
      external_z
    ) %*%
      weights
  )
  
  names(
    score_all
  ) <- colnames(
    external_z
  )
  
  survival_result <- extract_external_os(
    pheno =
      pheno,
    accession =
      accession
  )
  
  clinical <- survival_result$clinical
  
  write.csv(
    survival_result$audit,
    file.path(
      cohort_result_dir,
      paste0(
        accession,
        "_survival_column_audit.csv"
      )
    ),
    row.names = FALSE
  )
  
  if (
    !all(
      clinical$sample_id %in%
      names(
        score_all
      )
    )
  ) {
    missing_samples <- setdiff(
      clinical$sample_id,
      names(
        score_all
      )
    )
    
    stop(
      accession,
      ": ",
      length(
        missing_samples
      ),
      " survival samples are absent from the prepared expression matrix."
    )
  }
  
  validation <- clinical |>
    dplyr::mutate(
      MPCDS =
        score_all[
          match(
            sample_id,
            names(
              score_all
            )
          )
        ]
    )
  
  stopifnot(
    !anyNA(
      validation$MPCDS
    ),
    !any(
      is.infinite(
        validation$MPCDS
      )
    )
  )
  
  validation$MPCDS_z <- as.numeric(
    scale(
      validation$MPCDS
    )
  )
  
  # Primary C-index.
  cindex_object <- survival::concordance(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      MPCDS,
    data =
      validation,
    reverse =
      TRUE
  )
  
  cindex <- cindex_object$concordance
  cindex_se <- sqrt(
    cindex_object$var
  )
  
  cindex_lower <- max(
    0,
    cindex -
      1.96 *
      cindex_se
  )
  
  cindex_upper <- min(
    1,
    cindex +
      1.96 *
      cindex_se
  )
  
  # Continuous Cox HR per external-cohort SD.
  continuous_cox <- survival::coxph(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      MPCDS_z,
    data =
      validation,
    x =
      TRUE
  )
  
  continuous_summary <- summary(
    continuous_cox
  )
  
  continuous_hr <- continuous_summary$conf.int[
    "MPCDS_z",
    "exp(coef)"
  ]
  
  continuous_lower <- continuous_summary$conf.int[
    "MPCDS_z",
    "lower .95"
  ]
  
  continuous_upper <- continuous_summary$conf.int[
    "MPCDS_z",
    "upper .95"
  ]
  
  continuous_p <- continuous_summary$coefficients[
    "MPCDS_z",
    "Pr(>|z|)"
  ]
  
  ph_test <- survival::cox.zph(
    continuous_cox
  )
  
  ph_p <- ph_test$table[
    "MPCDS_z",
    "p"
  ]
  
  ph_table <- as.data.frame(
    ph_test$table
  )
  
  ph_table$term <- rownames(
    ph_table
  )
  
  rownames(
    ph_table
  ) <- NULL
  
  # Descriptive median split only.
  cutoff <- stats::median(
    validation$MPCDS
  )
  
  validation$risk_group <- ifelse(
    validation$MPCDS >
      cutoff,
    "High MPCDS",
    "Low MPCDS"
  )
  
  validation$risk_group <- factor(
    validation$risk_group,
    levels =
      c(
        "Low MPCDS",
        "High MPCDS"
      )
  )
  
  km_fit <- survival::survfit(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    data =
      validation
  )
  
  logrank <- survival::survdiff(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    data =
      validation
  )
  
  logrank_p <- 1 -
    stats::pchisq(
      logrank$chisq,
      df =
        length(
          logrank$n
        ) -
        1
    )
  
  group_cox <- survival::coxph(
    survival::Surv(
      OS_time,
      OS_status
    ) ~
      risk_group,
    data =
      validation
  )
  
  group_summary <- summary(
    group_cox
  )
  
  group_hr <- group_summary$conf.int[
    "risk_groupHigh MPCDS",
    "exp(coef)"
  ]
  
  group_lower <- group_summary$conf.int[
    "risk_groupHigh MPCDS",
    "lower .95"
  ]
  
  group_upper <- group_summary$conf.int[
    "risk_groupHigh MPCDS",
    "upper .95"
  ]
  
  group_p <- group_summary$coefficients[
    "risk_groupHigh MPCDS",
    "Pr(>|z|)"
  ]
  
  raw_weight_coverage <- 100 *
    sum(
      abs(
        beta_47[
          available_genes
        ]
      )
    ) /
    sum(
      abs(
        beta_47
      )
    )
  
  standardized_weight_coverage <- 100 *
    sum(
      abs(
        standardized_weights_47[
          available_genes
        ]
      )
    ) /
    sum(
      abs(
        standardized_weights_47
      )
    )
  
  summary_row <- data.frame(
    accession =
      accession,
    n_expression_samples =
      ncol(
        external_z
      ),
    n_OS_samples =
      nrow(
        validation
      ),
    n_events =
      sum(
        validation$OS_status
      ),
    n_genes_available =
      length(
        available_genes
      ),
    missing_genes =
      ifelse(
        length(
          missing_genes
        ) >
          0,
        paste(
          missing_genes,
          collapse = ";"
        ),
        "None"
      ),
    raw_coefficient_weight_coverage_percent =
      raw_weight_coverage,
    standardized_weight_coverage_percent =
      standardized_weight_coverage,
    C_index =
      cindex,
    C_index_lower95 =
      cindex_lower,
    C_index_upper95 =
      cindex_upper,
    continuous_HR_per_SD =
      continuous_hr,
    continuous_HR_lower95 =
      continuous_lower,
    continuous_HR_upper95 =
      continuous_upper,
    continuous_p =
      continuous_p,
    PH_p =
      ph_p,
    median_cutoff =
      cutoff,
    high_vs_low_HR =
      group_hr,
    high_vs_low_lower95 =
      group_lower,
    high_vs_low_upper95 =
      group_upper,
    high_vs_low_Cox_p =
      group_p,
    logrank_p =
      logrank_p,
    stringsAsFactors =
      FALSE
  )
  
  write.csv(
    validation,
    file.path(
      cohort_result_dir,
      paste0(
        accession,
        "_MPCDS_patient_scores.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      gene =
        available_genes,
      standardized_weight =
        as.numeric(
          weights
        ),
      stringsAsFactors =
        FALSE
    ),
    file.path(
      cohort_result_dir,
      paste0(
        accession,
        "_MPCDS_platform_weights.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    ph_table,
    file.path(
      cohort_result_dir,
      paste0(
        accession,
        "_MPCDS_PH_diagnostics.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    summary_row,
    file.path(
      cohort_result_dir,
      paste0(
        accession,
        "_MPCDS_validation_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  km_table <- as.data.frame(
    summary(
      km_fit
    )$table
  )
  
  km_table$risk_group <- gsub(
    "risk_group=",
    "",
    rownames(
      km_table
    )
  )
  
  rownames(
    km_table
  ) <- NULL
  
  write.csv(
    km_table,
    file.path(
      cohort_result_dir,
      paste0(
        accession,
        "_MPCDS_KM_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  cat(
    "Expression samples: ",
    summary_row$n_expression_samples,
    "\n",
    "Usable OS samples: ",
    summary_row$n_OS_samples,
    "\n",
    "Events: ",
    summary_row$n_events,
    "\n",
    "MPCDS genes: ",
    summary_row$n_genes_available,
    "/47\n",
    "C-index: ",
    round(
      summary_row$C_index,
      4
    ),
    "\n",
    "Continuous HR/SD: ",
    round(
      summary_row$continuous_HR_per_SD,
      4
    ),
    ", p = ",
    signif(
      summary_row$continuous_p,
      4
    ),
    "\n",
    "PH p: ",
    signif(
      summary_row$PH_p,
      4
    ),
    "\n",
    "Median-split log-rank p: ",
    signif(
      summary_row$logrank_p,
      4
    ),
    "\n",
    sep = ""
  )
  
  list(
    summary =
      summary_row,
    validation =
      validation,
    phenotype =
      pheno
  )
}


# ------------------------------------------------------------
# 11. Run all four LUAD external validations
# ------------------------------------------------------------

cohort_results <- vector(
  "list",
  nrow(
    luad_manifest
  )
)

for (
  i in seq_len(
    nrow(
      luad_manifest
    )
  )
) {
  
  accession <- luad_manifest$accession[
    i
  ]
  
  cohort_results[[i]] <- tryCatch(
    {
      validate_mpcds_cohort(
        accession
      )
    },
    error = function(
    e
    ) {
      
      warning(
        accession,
        " MPCDS validation failed: ",
        conditionMessage(
          e
        )
      )
      
      list(
        summary =
          data.frame(
            accession =
              accession,
            n_expression_samples =
              NA_integer_,
            n_OS_samples =
              NA_integer_,
            n_events =
              NA_integer_,
            n_genes_available =
              NA_integer_,
            missing_genes =
              NA_character_,
            raw_coefficient_weight_coverage_percent =
              NA_real_,
            standardized_weight_coverage_percent =
              NA_real_,
            C_index =
              NA_real_,
            C_index_lower95 =
              NA_real_,
            C_index_upper95 =
              NA_real_,
            continuous_HR_per_SD =
              NA_real_,
            continuous_HR_lower95 =
              NA_real_,
            continuous_HR_upper95 =
              NA_real_,
            continuous_p =
              NA_real_,
            PH_p =
              NA_real_,
            median_cutoff =
              NA_real_,
            high_vs_low_HR =
              NA_real_,
            high_vs_low_lower95 =
              NA_real_,
            high_vs_low_upper95 =
              NA_real_,
            high_vs_low_Cox_p =
              NA_real_,
            logrank_p =
              NA_real_,
            validation_status =
              paste0(
                "FAILED: ",
                conditionMessage(
                  e
                )
              ),
            stringsAsFactors =
              FALSE
          ),
        validation =
          NULL,
        phenotype =
          NULL
      )
    }
  )
}


# ------------------------------------------------------------
# 12. Combine primary external-validation results
# ------------------------------------------------------------

external_summary <- dplyr::bind_rows(
  lapply(
    cohort_results,
    function(
    x
    ) {
      x$summary
    }
  )
)

if (
  !"validation_status" %in%
  colnames(
    external_summary
  )
) {
  external_summary$validation_status <-
    "SUCCESS"
} else {
  external_summary$validation_status[
    is.na(
      external_summary$validation_status
    )
  ] <- "SUCCESS"
}

write.csv(
  external_summary,
  file.path(
    external_mpcds_dir,
    "LUAD_MPCDS_external_GEO_validation_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. GSE68465 clinically adjusted sensitivity analysis
#
# Preserves the previously validated clinical model:
#   age + sex + pT + pN
#
# This block is deliberately GSE68465-specific rather than
# forcing non-equivalent clinical covariates across cohorts.
# ------------------------------------------------------------

gse68465_index <- match(
  "GSE68465",
  luad_manifest$accession
)

gse68465_result <- cohort_results[[gse68465_index]]

if (
  !is.null(
    gse68465_result$validation
  ) &&
  !is.null(
    gse68465_result$phenotype
  )
) {
  
  validation <- gse68465_result$validation
  pheno <- gse68465_result$phenotype
  
  required_gse68465_columns <- c(
    "geo_accession",
    "age:ch1",
    "Sex:ch1",
    "disease_stage:ch1"
  )
  
  if (
    all(
      required_gse68465_columns %in%
      colnames(
        pheno
      )
    )
  ) {
    
    clinical_extra <- data.frame(
      sample_id =
        pheno$geo_accession,
      age =
        numeric_from_text(
          pheno[["age:ch1"]]
        ),
      sex =
        as.character(
          pheno[["Sex:ch1"]]
        ),
      stage_original =
        as.character(
          pheno[["disease_stage:ch1"]]
        ),
      stringsAsFactors =
        FALSE
    )
    
    adjusted_data <- validation |>
      dplyr::left_join(
        clinical_extra,
        by =
          "sample_id"
      ) |>
      dplyr::mutate(
        pN =
          dplyr::case_when(
            grepl(
              "^pN0",
              stage_original
            ) ~
              "N0",
            grepl(
              "^pN1",
              stage_original
            ) ~
              "N1",
            grepl(
              "^pN2",
              stage_original
            ) ~
              "N2",
            TRUE ~
              NA_character_
          ),
        pT =
          dplyr::case_when(
            grepl(
              "pT1$",
              stage_original
            ) ~
              "T1",
            grepl(
              "pT2$",
              stage_original
            ) ~
              "T2",
            grepl(
              "pT3$",
              stage_original
            ) ~
              "T3",
            grepl(
              "pT4$",
              stage_original
            ) ~
              "T4",
            TRUE ~
              NA_character_
          ),
        sex =
          factor(
            sex,
            levels =
              c(
                "Female",
                "Male"
              )
          ),
        pN =
          factor(
            pN,
            levels =
              c(
                "N0",
                "N1",
                "N2"
              )
          ),
        pT =
          factor(
            pT,
            levels =
              c(
                "T1",
                "T2",
                "T3",
                "T4"
              )
          ),
        age_10yr =
          age /
          10
      ) |>
      dplyr::filter(
        !is.na(
          MPCDS_z
        ),
        !is.na(
          age_10yr
        ),
        !is.na(
          sex
        ),
        !is.na(
          pN
        ),
        !is.na(
          pT
        )
      )
    
    if (
      nrow(
        adjusted_data
      ) >=
      100
    ) {
      
      clinical_fit <- survival::coxph(
        survival::Surv(
          OS_time,
          OS_status
        ) ~
          age_10yr +
          sex +
          pT +
          pN,
        data =
          adjusted_data,
        x =
          TRUE
      )
      
      combined_fit <- survival::coxph(
        survival::Surv(
          OS_time,
          OS_status
        ) ~
          MPCDS_z +
          age_10yr +
          sex +
          pT +
          pN,
        data =
          adjusted_data,
        x =
          TRUE
      )
      
      lrt <- stats::anova(
        clinical_fit,
        combined_fit,
        test =
          "LRT"
      )
      
      combined_summary <- summary(
        combined_fit
      )
      
      combined_ph <- survival::cox.zph(
        combined_fit
      )
      
      clinical_c <- summary(
        clinical_fit
      )$concordance[
        1
      ]
      
      combined_c <- combined_summary$concordance[
        1
      ]
      
      adjusted_summary <- data.frame(
        accession =
          "GSE68465",
        n =
          nrow(
            adjusted_data
          ),
        events =
          sum(
            adjusted_data$OS_status
          ),
        adjusted_HR_per_SD =
          combined_summary$conf.int[
            "MPCDS_z",
            "exp(coef)"
          ],
        CI_lower =
          combined_summary$conf.int[
            "MPCDS_z",
            "lower .95"
          ],
        CI_upper =
          combined_summary$conf.int[
            "MPCDS_z",
            "upper .95"
          ],
        p_value =
          combined_summary$coefficients[
            "MPCDS_z",
            "Pr(>|z|)"
          ],
        clinical_C =
          clinical_c,
        clinical_plus_MPCDS_C =
          combined_c,
        delta_C =
          combined_c -
          clinical_c,
        LRT_p_value =
          lrt[
            2,
            "Pr(>|Chi|)"
          ],
        MPCDS_PH_p_value =
          combined_ph$table[
            "MPCDS_z",
            "p"
          ],
        stringsAsFactors =
          FALSE
      )
      
      write.csv(
        adjusted_summary,
        file.path(
          external_mpcds_dir,
          "GSE68465_clinical_adjusted_MPCDS_summary.csv"
        ),
        row.names = FALSE
      )
      
      cat(
        "\nGSE68465 clinical-adjusted sensitivity:\n"
      )
      
      print(
        adjusted_summary
      )
    }
  }
}


# ------------------------------------------------------------
# 14. Cross-cohort descriptive consistency summary
# ------------------------------------------------------------

successful_summary <- external_summary |>
  dplyr::filter(
    validation_status ==
      "SUCCESS"
  )

if (
  nrow(
    successful_summary
  ) >
  0
) {
  
  cross_cohort_summary <- data.frame(
    n_LUAD_external_cohorts =
      nrow(
        successful_summary
      ),
    total_OS_samples =
      sum(
        successful_summary$n_OS_samples
      ),
    total_events =
      sum(
        successful_summary$n_events
      ),
    cohorts_HR_above_1 =
      sum(
        successful_summary$continuous_HR_per_SD >
          1
      ),
    cohorts_continuous_p_below_0_05 =
      sum(
        successful_summary$continuous_p <
          0.05
      ),
    median_C_index =
      stats::median(
        successful_summary$C_index
      ),
    min_C_index =
      min(
        successful_summary$C_index
      ),
    max_C_index =
      max(
        successful_summary$C_index
      ),
    stringsAsFactors =
      FALSE
  )
  
} else {
  
  cross_cohort_summary <- data.frame(
    n_LUAD_external_cohorts =
      0L,
    total_OS_samples =
      0L,
    total_events =
      0L,
    cohorts_HR_above_1 =
      0L,
    cohorts_continuous_p_below_0_05 =
      0L,
    median_C_index =
      NA_real_,
    min_C_index =
      NA_real_,
    max_C_index =
      NA_real_,
    stringsAsFactors =
      FALSE
  )
}

write.csv(
  cross_cohort_summary,
  file.path(
    external_mpcds_dir,
    "LUAD_MPCDS_external_GEO_cross_cohort_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Final verification and console summary
# ------------------------------------------------------------

n_failed <- sum(
  external_summary$validation_status !=
    "SUCCESS"
)

cat(
  "\n============================================================\n",
  "SCRIPT 21 LUAD MPCDS EXTERNAL VALIDATION SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  external_summary
)

cat(
  "\nCross-cohort summary:\n"
)

print(
  cross_cohort_summary
)

cat(
  "\nTCGA score reconstruction correlation: ",
  format(
    tcga_score_correlation,
    digits = 8
  ),
  "\n",
  "LUAD GEO cohorts attempted: ",
  nrow(
    luad_manifest
  ),
  "\n",
  "Successfully validated: ",
  nrow(
    luad_manifest
  ) -
    n_failed,
  "\n",
  "Failed: ",
  n_failed,
  "\n",
  "Outputs: ",
  external_mpcds_dir,
  "\n",
  sep = ""
)

if (
  n_failed >
  0
) {
  
  cat(
    "\nIMPORTANT: One or more cohorts failed survival parsing or ",
    "MPCDS validation. Inspect the warning and the phenotype ",
    "columns before interpreting cross-cohort results.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "\nAll four independent LUAD GEO cohorts were evaluated using ",
    "fixed TCGA-derived MPCDS weights without outcome refitting.\n",
    sep = ""
  )
}

cat(
  "============================================================\n",
  "\n========================================\n",
  "SCRIPT 21 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================