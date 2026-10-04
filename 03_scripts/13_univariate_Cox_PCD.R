# ============================================================
# 13_univariate_Cox_PCD.R
#
# Purpose:
# Perform univariate overall-survival Cox regression for the
# 47 concordant PCD candidates in TCGA-LUAD and TCGA-LUSC,
# followed by cross-cohort fixed-effect meta-analysis and
# heterogeneity assessment.
#
# Inputs:
# 02_processed_data/survival/
# - TCGA_LUAD_PCD_survival_expression.rds
# - TCGA_LUAD_PCD_survival_metadata.csv
# - TCGA_LUSC_PCD_survival_expression.rds
# - TCGA_LUSC_PCD_survival_metadata.csv
#
# Outputs:
# 04_results/survival/
# - TCGA_LUAD_univariate_Cox_47_PCD.csv
# - TCGA_LUSC_univariate_Cox_47_PCD.csv
# - TCGA_LUAD_FDR_significant_PCD_genes.csv
# - TCGA_LUAD_LUSC_Cox_meta_analysis.csv
#
# Survival endpoint:
# Overall Survival (OS)
#
# Expression:
# log2(TPM + 1)
#
# Multiple testing:
# Benjamini-Hochberg FDR
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

survival_processed_dir <- file.path(
  processed_dir,
  "survival"
)

survival_results_dir <- file.path(
  results_dir,
  "survival"
)

ensure_dir(survival_results_dir)


# ------------------------------------------------------------
# 2. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "survival",
  "dplyr"
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
  library(survival)
  library(dplyr)
})

# Explicit namespaces are used for dplyr verbs below to avoid
# function masking in interactive sessions.


# ------------------------------------------------------------
# 2A. Required input files
# ------------------------------------------------------------

required_input_files <- c(
  file.path(
    survival_processed_dir,
    "TCGA_LUAD_PCD_survival_expression.rds"
  ),
  file.path(
    survival_processed_dir,
    "TCGA_LUAD_PCD_survival_metadata.csv"
  ),
  file.path(
    survival_processed_dir,
    "TCGA_LUSC_PCD_survival_expression.rds"
  ),
  file.path(
    survival_processed_dir,
    "TCGA_LUSC_PCD_survival_metadata.csv"
  )
)

check_files_exist(
  required_input_files,
  label = "Script 13 input file(s)"
)


# ------------------------------------------------------------
# 3. Load survival-expression datasets
# ------------------------------------------------------------

luad_survival_expression <- readRDS(
  file.path(
    survival_processed_dir,
    "TCGA_LUAD_PCD_survival_expression.rds"
  )
)

lusc_survival_expression <- readRDS(
  file.path(
    survival_processed_dir,
    "TCGA_LUSC_PCD_survival_expression.rds"
  )
)

luad_survival_metadata <- read.csv(
  file.path(
    survival_processed_dir,
    "TCGA_LUAD_PCD_survival_metadata.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lusc_survival_metadata <- read.csv(
  file.path(
    survival_processed_dir,
    "TCGA_LUSC_PCD_survival_metadata.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 4. Validate input dimensions/alignment
# ------------------------------------------------------------

stopifnot(
  identical(
    dim(luad_survival_expression),
    c(47L, 504L)
  )
)

stopifnot(
  identical(
    dim(lusc_survival_expression),
    c(47L, 493L)
  )
)

stopifnot(
  identical(
    colnames(luad_survival_expression),
    luad_survival_metadata$submitter_id
  )
)

stopifnot(
  identical(
    colnames(lusc_survival_expression),
    lusc_survival_metadata$submitter_id
  )
)

stopifnot(
  !anyNA(
    luad_survival_expression
  ),
  !anyNA(
    lusc_survival_expression
  )
)

stopifnot(
  all(
    luad_survival_metadata$OS_time > 0
  ),
  all(
    lusc_survival_metadata$OS_time > 0
  )
)

stopifnot(
  all(
    luad_survival_metadata$OS_status %in% c(0, 1)
  ),
  all(
    lusc_survival_metadata$OS_status %in% c(0, 1)
  )
)


# ------------------------------------------------------------
# 5. Candidate genes
# ------------------------------------------------------------

survival_candidate_genes <-
  rownames(
    luad_survival_expression
  )

stopifnot(
  length(survival_candidate_genes) == 47
)

stopifnot(
  identical(
    survival_candidate_genes,
    rownames(lusc_survival_expression)
  )
)


# ------------------------------------------------------------
# 6. Build Cox datasets
# ------------------------------------------------------------

luad_cox_data <- data.frame(
  patient_id =
    luad_survival_metadata$submitter_id,
  
  OS_time =
    luad_survival_metadata$OS_time,
  
  OS_status =
    luad_survival_metadata$OS_status,
  
  t(luad_survival_expression),
  
  check.names = FALSE
)

lusc_cox_data <- data.frame(
  patient_id =
    lusc_survival_metadata$submitter_id,
  
  OS_time =
    lusc_survival_metadata$OS_time,
  
  OS_status =
    lusc_survival_metadata$OS_status,
  
  t(lusc_survival_expression),
  
  check.names = FALSE
)

stopifnot(
  identical(
    dim(luad_cox_data),
    c(504L, 50L)
  )
)

stopifnot(
  identical(
    dim(lusc_cox_data),
    c(493L, 50L)
  )
)


# ============================================================
# PART A — LUAD UNIVARIATE COX
# ============================================================


# ------------------------------------------------------------
# 7. Run LUAD Cox models
# ------------------------------------------------------------

luad_univ_cox_full <- lapply(
  survival_candidate_genes,
  function(gene) {
    
    fit <- coxph(
      Surv(
        OS_time,
        OS_status
      ) ~ luad_cox_data[[gene]],
      data = luad_cox_data
    )
    
    fit_summary <- summary(
      fit
    )
    
    data.frame(
      gene = gene,
      
      beta =
        fit_summary$coefficients[
          1,
          "coef"
        ],
      
      SE =
        fit_summary$coefficients[
          1,
          "se(coef)"
        ],
      
      HR =
        fit_summary$coefficients[
          1,
          "exp(coef)"
        ],
      
      lower_95CI =
        fit_summary$conf.int[
          1,
          "lower .95"
        ],
      
      upper_95CI =
        fit_summary$conf.int[
          1,
          "upper .95"
        ],
      
      pvalue =
        fit_summary$coefficients[
          1,
          "Pr(>|z|)"
        ],
      
      stringsAsFactors = FALSE
    )
  }
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    padj = p.adjust(
      pvalue,
      method = "BH"
    ),
    
    direction =
      dplyr::if_else(
        HR > 1,
        "Risk",
        "Protective"
      )
  ) |>
  dplyr::arrange(
    pvalue
  )


# ------------------------------------------------------------
# 8. LUAD validation
# ------------------------------------------------------------

stopifnot(
  nrow(luad_univ_cox_full) == 47
)

stopifnot(
  sum(
    luad_univ_cox_full$pvalue < 0.05
  ) == 15
)

stopifnot(
  sum(
    luad_univ_cox_full$padj < 0.05
  ) == 9
)

stopifnot(
  sum(
    is.na(luad_univ_cox_full$SE) |
      luad_univ_cox_full$SE <= 0
  ) == 0
)


# ------------------------------------------------------------
# 9. LUAD FDR-significant genes
# ------------------------------------------------------------

luad_fdr_significant_pcd <-
  luad_univ_cox_full |>
  dplyr::filter(
    padj < 0.05
  )

stopifnot(
  nrow(luad_fdr_significant_pcd) == 9
)

stopifnot(
  setequal(
    luad_fdr_significant_pcd$gene,
    c(
      "DSG2",
      "DSG3",
      "LMNB1",
      "CLSPN",
      "DAPK2",
      "EGLN3",
      "GCLC",
      "TXNRD1",
      "SLC7A11"
    )
  )
)


# ============================================================
# PART B — LUSC UNIVARIATE COX
# ============================================================


# ------------------------------------------------------------
# 10. Run LUSC Cox models
# ------------------------------------------------------------

lusc_univ_cox_full <- lapply(
  survival_candidate_genes,
  function(gene) {
    
    fit <- coxph(
      Surv(
        OS_time,
        OS_status
      ) ~ lusc_cox_data[[gene]],
      data = lusc_cox_data
    )
    
    fit_summary <- summary(
      fit
    )
    
    data.frame(
      gene = gene,
      
      beta =
        fit_summary$coefficients[
          1,
          "coef"
        ],
      
      SE =
        fit_summary$coefficients[
          1,
          "se(coef)"
        ],
      
      HR =
        fit_summary$coefficients[
          1,
          "exp(coef)"
        ],
      
      lower_95CI =
        fit_summary$conf.int[
          1,
          "lower .95"
        ],
      
      upper_95CI =
        fit_summary$conf.int[
          1,
          "upper .95"
        ],
      
      pvalue =
        fit_summary$coefficients[
          1,
          "Pr(>|z|)"
        ],
      
      stringsAsFactors = FALSE
    )
  }
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    padj = p.adjust(
      pvalue,
      method = "BH"
    ),
    
    direction =
      dplyr::if_else(
        HR > 1,
        "Risk",
        "Protective"
      )
  ) |>
  dplyr::arrange(
    pvalue
  )


# ------------------------------------------------------------
# 11. LUSC validation
# ------------------------------------------------------------

stopifnot(
  nrow(lusc_univ_cox_full) == 47
)

stopifnot(
  sum(
    lusc_univ_cox_full$pvalue < 0.05
  ) == 2
)

stopifnot(
  sum(
    lusc_univ_cox_full$padj < 0.05
  ) == 0
)

stopifnot(
  sum(
    is.na(lusc_univ_cox_full$SE) |
      lusc_univ_cox_full$SE <= 0
  ) == 0
)


# ============================================================
# PART C — CROSS-COHORT META-ANALYSIS
# ============================================================


# ------------------------------------------------------------
# 12. Join LUAD and LUSC Cox estimates
# ------------------------------------------------------------

cox_meta_input <- luad_univ_cox_full |>
  dplyr::select(
    gene,
    LUAD_beta = beta,
    LUAD_SE = SE,
    LUAD_HR = HR,
    LUAD_pvalue = pvalue,
    LUAD_padj = padj
  ) |>
  inner_join(
    lusc_univ_cox_full |>
      dplyr::select(
        gene,
        LUSC_beta = beta,
        LUSC_SE = SE,
        LUSC_HR = HR,
        LUSC_pvalue = pvalue,
        LUSC_padj = padj
      ),
    by = "gene"
  )


# ------------------------------------------------------------
# 13. Cox direction concordance
# ------------------------------------------------------------

cox_meta_input <- cox_meta_input |>
  dplyr::mutate(
    direction_concordance =
      case_when(
        
        LUAD_beta > 0 &
          LUSC_beta > 0 ~
          "Risk in both",
        
        LUAD_beta < 0 &
          LUSC_beta < 0 ~
          "Protective in both",
        
        TRUE ~
          "Opposite direction"
      )
  )

stopifnot(
  sum(
    cox_meta_input$direction_concordance ==
      "Risk in both"
  ) == 14
)

stopifnot(
  sum(
    cox_meta_input$direction_concordance ==
      "Protective in both"
  ) == 3
)

stopifnot(
  sum(
    cox_meta_input$direction_concordance ==
      "Opposite direction"
  ) == 30
)


# ------------------------------------------------------------
# 14. Fixed-effect inverse-variance meta-analysis
# ------------------------------------------------------------

cox_meta_results <- cox_meta_input |>
  rowwise() |>
  dplyr::mutate(
    
    w_LUAD =
      1 / (LUAD_SE^2),
    
    w_LUSC =
      1 / (LUSC_SE^2),
    
    meta_beta =
      (
        w_LUAD * LUAD_beta +
          w_LUSC * LUSC_beta
      ) /
      (
        w_LUAD +
          w_LUSC
      ),
    
    meta_SE =
      sqrt(
        1 /
          (
            w_LUAD +
              w_LUSC
          )
      ),
    
    meta_HR =
      exp(meta_beta),
    
    meta_lower_95CI =
      exp(
        meta_beta -
          1.96 * meta_SE
      ),
    
    meta_upper_95CI =
      exp(
        meta_beta +
          1.96 * meta_SE
      ),
    
    meta_z =
      meta_beta /
      meta_SE,
    
    meta_pvalue =
      2 *
      pnorm(
        abs(meta_z),
        lower.tail = FALSE
      ),
    
    Q =
      w_LUAD *
      (
        LUAD_beta -
          meta_beta
      )^2 +
      w_LUSC *
      (
        LUSC_beta -
          meta_beta
      )^2,
    
    heterogeneity_p =
      pchisq(
        Q,
        df = 1,
        lower.tail = FALSE
      ),
    
    I2 =
      ifelse(
        Q > 1,
        100 *
          (
            Q - 1
          ) /
          Q,
        0
      )
  ) |>
  ungroup() |>
  dplyr::mutate(
    
    meta_padj =
      p.adjust(
        meta_pvalue,
        method = "BH"
      ),
    
    meta_direction =
      dplyr::if_else(
        meta_HR > 1,
        "Risk",
        "Protective"
      )
  ) |>
  dplyr::arrange(
    meta_pvalue
  )


# ------------------------------------------------------------
# 15. Meta-analysis validation
# ------------------------------------------------------------

stopifnot(
  nrow(cox_meta_results) == 47
)

stopifnot(
  sum(
    cox_meta_results$meta_pvalue < 0.05
  ) == 8
)

stopifnot(
  sum(
    cox_meta_results$meta_padj < 0.05
  ) == 0
)


# ============================================================
# PART D — SAVE EXACT FINAL OUTPUTS
# ============================================================


# ------------------------------------------------------------
# 16. LUAD univariate Cox
# ------------------------------------------------------------

write.csv(
  luad_univ_cox_full,
  file.path(
    survival_results_dir,
    "TCGA_LUAD_univariate_Cox_47_PCD.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. LUSC univariate Cox
# ------------------------------------------------------------

write.csv(
  lusc_univ_cox_full,
  file.path(
    survival_results_dir,
    "TCGA_LUSC_univariate_Cox_47_PCD.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. LUAD FDR-significant PCD genes
# ------------------------------------------------------------

write.csv(
  luad_fdr_significant_pcd,
  file.path(
    survival_results_dir,
    "TCGA_LUAD_FDR_significant_PCD_genes.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. LUAD/LUSC Cox meta-analysis
# ------------------------------------------------------------

write.csv(
  cox_meta_results,
  file.path(
    survival_results_dir,
    "TCGA_LUAD_LUSC_Cox_meta_analysis.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART E — FINAL OUTPUT VALIDATION
# ============================================================

expected_survival_files <- c(
  "TCGA_LUAD_univariate_Cox_47_PCD.csv",
  "TCGA_LUSC_univariate_Cox_47_PCD.csv",
  "TCGA_LUAD_FDR_significant_PCD_genes.csv",
  "TCGA_LUAD_LUSC_Cox_meta_analysis.csv"
)

stopifnot(
  all(
    file.exists(
      file.path(
        survival_results_dir,
        expected_survival_files
      )
    )
  )
)


# ============================================================
# PART F — COMPLETION SUMMARY
# ============================================================

cat(
  "\n====================================================\n"
)

cat(
  "PCD survival analysis completed.\n"
)

cat(
  "====================================================\n\n"
)

cat(
  "LUAD:\n",
  "  Patients: 504\n",
  "  Candidate genes tested: 47\n",
  "  Nominal p < 0.05: 15\n",
  "  BH-FDR < 0.05: 9\n\n",
  sep = ""
)

cat(
  "LUSC:\n",
  "  Patients: 493\n",
  "  Candidate genes tested: 47\n",
  "  Nominal p < 0.05: 2\n",
  "  BH-FDR < 0.05: 0\n\n",
  sep = ""
)

cat(
  "Cross-cohort Cox directions:\n",
  "  Risk in both: 14\n",
  "  Protective in both: 3\n",
  "  Opposite direction: 30\n\n",
  sep = ""
)

cat(
  "Fixed-effect meta-analysis:\n",
  "  Nominal p < 0.05: 8\n",
  "  BH-FDR < 0.05: 0\n\n",
  sep = ""
)

cat(
  "Final survival result files: 4\n"
)

cat(
  "\n========================================\n",
  "SCRIPT 13 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================