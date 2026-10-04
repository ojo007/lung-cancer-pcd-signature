# ============================================================
# SCRIPT 24
# DepMap / PRISM pharmacogenomic analysis of LUAD MPCDS
#
# Thesis:
# Exploring a Specialized Programmed Cell Death Pattern to
# Predict Prognosis and Treatment Sensitivity in Lung Cancer
# by Machine Learning and Multi-Omics Analysis
#
# Purpose:
# Transfer the fixed TCGA-LUAD 47-gene MPCDS to DepMap LUAD
# cell lines and test its association with PRISM Repurposing
# Secondary drug-response AUC.
#
# Primary analysis:
#   Spearman correlation between continuous MPCDS and AUC
#
# Secondary analysis:
#   Linear regression of AUC on standardized MPCDS
#
# Important interpretation:
#   Lower PRISM AUC = greater relative drug sensitivity.
#
#   Negative association:
#       higher MPCDS -> lower AUC -> greater sensitivity
#
#   Positive association:
#       higher MPCDS -> higher AUC -> greater resistance
#
# This analysis is pharmacogenomic and hypothesis-generating.
# It is NOT clinical treatment-response validation.
# ============================================================


# ------------------------------------------------------------
# 1. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")

required_directories <-
  c(
    "01_raw_data/DepMap",
    "02_processed_data/drug_sensitivity",
    "04_results/drug_sensitivity"
  )

for (
  directory_path in
  required_directories
) {
  ensure_dir(
    directory_path
  )
}

stopifnot(
  all(
    dir.exists(
      required_directories
    )
  )
)


# ------------------------------------------------------------
# 2. Required packages
# ------------------------------------------------------------

required_packages <-
  c(
    "data.table",
    "dplyr"
  )


missing_packages <-
  required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      FUN.VALUE = logical(1),
      quietly = TRUE
    )
  ]


if (
  length(
    missing_packages
  ) > 0
) {
  
  stop(
    paste(
      "Missing required package(s):",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

depmap_model_file <-
  "01_raw_data/DepMap/Model.csv"


depmap_expression_file <-
  paste0(
    "01_raw_data/DepMap/",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
  )


depmap_prism_file <-
  paste0(
    "01_raw_data/DepMap/",
    "Drug_sensitivity_AUC_",
    "(PRISM_Repurposing_Secondary_Screen).csv"
  )


depmap_weight_file <-
  file.path(
    "04_results/mpcds/external_validation",
    "TCGA_LUAD_MPCDS_standardized_weights_47.csv"
  )


required_files <-
  c(
    depmap_model_file,
    depmap_expression_file,
    depmap_prism_file,
    depmap_weight_file
  )


if (
  any(
    !file.exists(
      required_files
    )
  )
) {
  
  stop(
    paste(
      "Missing required file(s):\n",
      paste(
        required_files[
          !file.exists(
            required_files
          )
        ],
        collapse = "\n"
      )
    )
  )
}


cat(
  "All required input files found.\n\n"
)


# ------------------------------------------------------------
# 4. Record input-file metadata
# ------------------------------------------------------------

depmap_input_file_metadata <-
  data.frame(
    
    file =
      basename(
        required_files
      ),
    
    size_bytes =
      file.info(
        required_files
      )$size,
    
    modified_time =
      as.character(
        file.info(
          required_files
        )$mtime
      ),
    
    stringsAsFactors =
      FALSE
  )


write.csv(
  depmap_input_file_metadata,
  paste0(
    "04_results/drug_sensitivity/",
    "DepMap_PRISM_input_file_metadata.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Load fixed TCGA-LUAD MPCDS standardized weights
# ------------------------------------------------------------

depmap_platform_weights_47 <-
  data.table::fread(
    depmap_weight_file,
    data.table = FALSE
  )


required_weight_columns <-
  c(
    "gene",
    "TCGA_raw_coefficient",
    "TCGA_gene_SD",
    "standardized_weight"
  )


if (
  !all(
    required_weight_columns %in%
    colnames(
      depmap_platform_weights_47
    )
  )
) {
  
  stop(
    "Standardized MPCDS weight file has unexpected columns."
  )
}


stopifnot(
  nrow(
    depmap_platform_weights_47
  ) == 47
)


stopifnot(
  anyDuplicated(
    depmap_platform_weights_47$gene
  ) == 0
)


stopifnot(
  !any(
    is.na(
      depmap_platform_weights_47$standardized_weight
    )
  )
)


depmap_mpcds_genes <-
  depmap_platform_weights_47$gene


depmap_beta_z_47 <-
  depmap_platform_weights_47$standardized_weight


names(
  depmap_beta_z_47
) <-
  depmap_platform_weights_47$gene


cat(
  "Fixed MPCDS genes loaded:",
  length(
    depmap_mpcds_genes
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 6. Load DepMap model metadata
# ------------------------------------------------------------

depmap_models <-
  data.table::fread(
    depmap_model_file,
    data.table = FALSE
  )


required_model_columns <-
  c(
    "ModelID",
    "CellLineName",
    "OncotreeLineage",
    "OncotreePrimaryDisease",
    "OncotreeSubtype",
    "OncotreeCode",
    "PrimaryOrMetastasis"
  )


if (
  !all(
    required_model_columns %in%
    colnames(
      depmap_models
    )
  )
) {
  
  stop(
    "Model.csv does not contain all required metadata columns."
  )
}


# ------------------------------------------------------------
# 7. Define LUAD models using exact OncoTree annotation
# ------------------------------------------------------------

depmap_luad_models <-
  depmap_models |>
  dplyr::filter(
    OncotreeLineage == "Lung",
    OncotreeCode == "LUAD"
  )


stopifnot(
  anyDuplicated(
    depmap_luad_models$ModelID
  ) == 0
)


cat(
  "DepMap LUAD models in metadata:",
  nrow(
    depmap_luad_models
  ),
  "\n"
)


# ------------------------------------------------------------
# 8. Inspect expression-file header
# ------------------------------------------------------------

depmap_expression_header <-
  data.table::fread(
    depmap_expression_file,
    nrows = 0,
    data.table = FALSE
  )


depmap_expression_columns <-
  colnames(
    depmap_expression_header
  )


required_expression_metadata <-
  c(
    "ModelID",
    "IsDefaultEntryForModel"
  )


if (
  !all(
    required_expression_metadata %in%
    depmap_expression_columns
  )
) {
  
  stop(
    "Expression file does not contain required model metadata."
  )
}


# ------------------------------------------------------------
# 9. Parse DepMap expression gene symbols
# ------------------------------------------------------------

expression_metadata_columns <-
  c(
    "V1",
    "SequencingID",
    "ModelConditionID",
    "ModelID",
    "IsDefaultEntryForMC",
    "IsDefaultEntryForModel"
  )


depmap_expression_gene_columns <-
  setdiff(
    depmap_expression_columns,
    expression_metadata_columns
  )


depmap_expression_gene_symbols <-
  sub(
    "\\s+\\([0-9]+\\)$",
    "",
    depmap_expression_gene_columns
  )


if (
  anyDuplicated(
    depmap_expression_gene_symbols
  ) > 0
) {
  
  stop(
    "Duplicate gene symbols detected after parsing expression header."
  )
}


depmap_gene_column_map <-
  data.frame(
    
    depmap_column =
      depmap_expression_gene_columns,
    
    gene =
      depmap_expression_gene_symbols,
    
    stringsAsFactors =
      FALSE
  )


depmap_mpcds_missing_47 <-
  setdiff(
    depmap_mpcds_genes,
    depmap_expression_gene_symbols
  )


if (
  length(
    depmap_mpcds_missing_47
  ) > 0
) {
  
  stop(
    paste(
      "Missing MPCDS gene(s) in DepMap expression:",
      paste(
        depmap_mpcds_missing_47,
        collapse = ", "
      )
    )
  )
}


depmap_mpcds_column_map <-
  depmap_gene_column_map |>
  dplyr::filter(
    gene %in%
      depmap_mpcds_genes
  ) |>
  dplyr::arrange(
    match(
      gene,
      depmap_mpcds_genes
    )
  )


stopifnot(
  nrow(
    depmap_mpcds_column_map
  ) == 47
)


stopifnot(
  anyDuplicated(
    depmap_mpcds_column_map$gene
  ) == 0
)


cat(
  "MPCDS genes available in DepMap expression: 47 / 47\n"
)


# ------------------------------------------------------------
# 10. Read only required expression columns
# ------------------------------------------------------------

depmap_expression_columns_needed <-
  c(
    "ModelID",
    "IsDefaultEntryForModel",
    depmap_mpcds_column_map$depmap_column
  )


depmap_expression_mpcds <-
  data.table::fread(
    depmap_expression_file,
    select =
      depmap_expression_columns_needed,
    data.table = FALSE
  )


# ------------------------------------------------------------
# 11. Rename expression columns to gene symbols
# ------------------------------------------------------------

for (
  i in
  seq_len(
    nrow(
      depmap_mpcds_column_map
    )
  )
) {
  
  old_name <-
    depmap_mpcds_column_map$depmap_column[
      i
    ]
  
  new_name <-
    depmap_mpcds_column_map$gene[
      i
    ]
  
  colnames(
    depmap_expression_mpcds
  )[
    colnames(
      depmap_expression_mpcds
    ) ==
      old_name
  ] <-
    new_name
}


stopifnot(
  all(
    depmap_mpcds_genes %in%
      colnames(
        depmap_expression_mpcds
      )
  )
)


# ------------------------------------------------------------
# 12. Retain default expression entry per model
# ------------------------------------------------------------

depmap_expression_default <-
  depmap_expression_mpcds |>
  dplyr::filter(
    IsDefaultEntryForModel == "Yes"
  )


stopifnot(
  anyDuplicated(
    depmap_expression_default$ModelID
  ) == 0
)


cat(
  "Default expression profiles:",
  nrow(
    depmap_expression_default
  ),
  "\n"
)


# ------------------------------------------------------------
# 13. Load PRISM Secondary AUC
# ------------------------------------------------------------

depmap_prism_auc <-
  data.table::fread(
    depmap_prism_file,
    data.table = FALSE
  )


colnames(
  depmap_prism_auc
)[1] <-
  "ModelID"


stopifnot(
  anyDuplicated(
    depmap_prism_auc$ModelID
  ) == 0
)


depmap_prism_model_ids <-
  depmap_prism_auc$ModelID


cat(
  "PRISM models:",
  nrow(
    depmap_prism_auc
  ),
  "\n"
)


cat(
  "PRISM compounds:",
  ncol(
    depmap_prism_auc
  ) - 1,
  "\n\n"
)


# ------------------------------------------------------------
# 14. Determine LUAD overlap
# ------------------------------------------------------------

depmap_luad_ids <-
  depmap_luad_models$ModelID


depmap_expression_default_ids <-
  depmap_expression_default$ModelID


n_luad_expression <-
  sum(
    depmap_luad_ids %in%
      depmap_expression_default_ids
  )


n_luad_prism <-
  sum(
    depmap_luad_ids %in%
      depmap_prism_model_ids
  )


depmap_luad_final_ids <-
  Reduce(
    intersect,
    list(
      depmap_luad_ids,
      depmap_expression_default_ids,
      depmap_prism_model_ids
    )
  )


cat(
  "LUAD with default expression:",
  n_luad_expression,
  "\n"
)


cat(
  "LUAD with PRISM:",
  n_luad_prism,
  "\n"
)


cat(
  "Complete LUAD expression + PRISM:",
  length(
    depmap_luad_final_ids
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 15. Final LUAD metadata
# ------------------------------------------------------------

depmap_luad_final_metadata <-
  depmap_luad_models |>
  dplyr::filter(
    ModelID %in%
      depmap_luad_final_ids
  ) |>
  dplyr::select(
    ModelID,
    CellLineName,
    OncotreeLineage,
    OncotreePrimaryDisease,
    OncotreeSubtype,
    OncotreeCode,
    PrimaryOrMetastasis
  )


stopifnot(
  anyDuplicated(
    depmap_luad_final_metadata$ModelID
  ) == 0
)


# ------------------------------------------------------------
# 16. Final 47-gene LUAD expression matrix
# ------------------------------------------------------------

depmap_luad_expression <-
  depmap_expression_default |>
  dplyr::filter(
    ModelID %in%
      depmap_luad_final_ids
  )


depmap_luad_expression <-
  depmap_luad_expression[
    match(
      depmap_luad_final_metadata$ModelID,
      depmap_luad_expression$ModelID
    ),
    ,
    drop = FALSE
  ]


stopifnot(
  identical(
    depmap_luad_expression$ModelID,
    depmap_luad_final_metadata$ModelID
  )
)


depmap_luad_x <-
  as.matrix(
    depmap_luad_expression[
      ,
      depmap_mpcds_genes,
      drop = FALSE
    ]
  )


storage.mode(
  depmap_luad_x
) <-
  "numeric"


stopifnot(
  ncol(
    depmap_luad_x
  ) == 47
)


stopifnot(
  sum(
    is.na(
      depmap_luad_x
    )
  ) == 0
)


# ------------------------------------------------------------
# 17. Check expression variance
# ------------------------------------------------------------

depmap_luad_gene_sd <-
  apply(
    depmap_luad_x,
    2,
    sd
  )


if (
  any(
    is.na(
      depmap_luad_gene_sd
    ) |
    depmap_luad_gene_sd <= 0
  )
) {
  
  stop(
    "At least one MPCDS gene has zero or undefined variance."
  )
}


# ------------------------------------------------------------
# 18. Standardize DepMap LUAD expression
# ------------------------------------------------------------

depmap_luad_x_z <-
  scale(
    depmap_luad_x
  )


stopifnot(
  !any(
    is.na(
      depmap_luad_x_z
    )
  )
)


# ------------------------------------------------------------
# 19. Align fixed TCGA standardized weights
# ------------------------------------------------------------

depmap_beta_z_final <-
  depmap_beta_z_47[
    colnames(
      depmap_luad_x_z
    )
  ]


stopifnot(
  length(
    depmap_beta_z_final
  ) == 47
)


stopifnot(
  !any(
    is.na(
      depmap_beta_z_final
    )
  )
)


stopifnot(
  identical(
    names(
      depmap_beta_z_final
    ),
    colnames(
      depmap_luad_x_z
    )
  )
)


# ------------------------------------------------------------
# 20. Calculate DepMap LUAD MPCDS
# ------------------------------------------------------------

depmap_luad_mpcds <-
  as.numeric(
    depmap_luad_x_z %*%
      depmap_beta_z_final
  )


depmap_luad_scores <-
  depmap_luad_final_metadata |>
  dplyr::mutate(
    
    MPCDS =
      depmap_luad_mpcds,
    
    MPCDS_z =
      as.numeric(
        scale(
          depmap_luad_mpcds
        )
      )
  )


stopifnot(
  !any(
    is.na(
      depmap_luad_scores$MPCDS
    )
  )
)


# ------------------------------------------------------------
# 21. Align PRISM to final LUAD cohort
# ------------------------------------------------------------

depmap_prism_luad <-
  depmap_prism_auc[
    match(
      depmap_luad_scores$ModelID,
      depmap_prism_auc$ModelID
    ),
    ,
    drop = FALSE
  ]


stopifnot(
  identical(
    depmap_prism_luad$ModelID,
    depmap_luad_scores$ModelID
  )
)


# ------------------------------------------------------------
# 22. Prespecified drug sample-size threshold
# ------------------------------------------------------------

depmap_min_n_per_drug <-
  30


depmap_drug_names <-
  colnames(
    depmap_prism_luad
  )[-1]


depmap_drug_n_final <-
  colSums(
    !is.na(
      depmap_prism_luad[
        ,
        -1,
        drop = FALSE
      ]
    )
  )


n_drugs_n30 <-
  sum(
    depmap_drug_n_final >=
      depmap_min_n_per_drug
  )


cat(
  "Compounds with n >= 30:",
  n_drugs_n30,
  "\n\n"
)


# ------------------------------------------------------------
# 23. Primary Spearman analysis
# ------------------------------------------------------------

depmap_drug_results_list <-
  lapply(
    depmap_drug_names,
    function(drug_name) {
      
      auc_values <-
        depmap_prism_luad[[drug_name]]
      
      keep <-
        !is.na(
          auc_values
        ) &
        !is.na(
          depmap_luad_scores$MPCDS_z
        )
      
      n_available <-
        sum(
          keep
        )
      
      if (
        n_available <
        depmap_min_n_per_drug
      ) {
        
        return(
          data.frame(
            
            drug_full_name =
              drug_name,
            
            n =
              n_available,
            
            auc_unique_values =
              length(
                unique(
                  auc_values[keep]
                )
              ),
            
            auc_SD =
              if (
                n_available > 1
              ) {
                sd(
                  auc_values[keep]
                )
              } else {
                NA_real_
              },
            
            spearman_rho =
              NA_real_,
            
            p_value =
              NA_real_,
            
            stringsAsFactors =
              FALSE
          )
        )
      }
      
      auc_unique_values <-
        length(
          unique(
            auc_values[keep]
          )
        )
      
      auc_sd <-
        sd(
          auc_values[keep]
        )
      
      if (
        auc_unique_values < 2 ||
        is.na(
          auc_sd
        ) ||
        auc_sd == 0
      ) {
        
        return(
          data.frame(
            
            drug_full_name =
              drug_name,
            
            n =
              n_available,
            
            auc_unique_values =
              auc_unique_values,
            
            auc_SD =
              auc_sd,
            
            spearman_rho =
              NA_real_,
            
            p_value =
              NA_real_,
            
            stringsAsFactors =
              FALSE
          )
        )
      }
      
      test_result <-
        suppressWarnings(
          cor.test(
            depmap_luad_scores$MPCDS_z[
              keep
            ],
            auc_values[
              keep
            ],
            method =
              "spearman",
            exact =
              FALSE
          )
        )
      
      data.frame(
        
        drug_full_name =
          drug_name,
        
        n =
          n_available,
        
        auc_unique_values =
          auc_unique_values,
        
        auc_SD =
          auc_sd,
        
        spearman_rho =
          unname(
            test_result$estimate
          ),
        
        p_value =
          test_result$p.value,
        
        stringsAsFactors =
          FALSE
      )
    }
  )


depmap_drug_results <-
  dplyr::bind_rows(
    depmap_drug_results_list
  )


# ------------------------------------------------------------
# 24. BH correction over actually tested compounds only
# ------------------------------------------------------------

depmap_drug_results$FDR <-
  NA_real_


depmap_tested_index <-
  !is.na(
    depmap_drug_results$p_value
  )


depmap_drug_results$FDR[
  depmap_tested_index
] <-
  p.adjust(
    depmap_drug_results$p_value[
      depmap_tested_index
    ],
    method =
      "BH"
  )


# ------------------------------------------------------------
# 25. Parse compound names and association direction
# ------------------------------------------------------------

depmap_drug_results <-
  depmap_drug_results |>
  dplyr::mutate(
    
    drug_name =
      sub(
        "\\s*\\(BRD:.*$",
        "",
        drug_full_name
      ),
    
    BRD_ID =
      ifelse(
        grepl(
          "\\(BRD:",
          drug_full_name
        ),
        sub(
          "^.*\\(BRD:([^\\)]+)\\)$",
          "\\1",
          drug_full_name
        ),
        NA_character_
      ),
    
    association_direction =
      dplyr::case_when(
        
        is.na(
          spearman_rho
        ) ~
          "Not tested",
        
        spearman_rho < 0 ~
          "Higher MPCDS = greater sensitivity",
        
        spearman_rho > 0 ~
          "Higher MPCDS = greater resistance",
        
        TRUE ~
          "No direction"
      )
  )


depmap_drug_results_tested <-
  depmap_drug_results |>
  dplyr::filter(
    !is.na(
      spearman_rho
    )
  ) |>
  dplyr::arrange(
    FDR,
    p_value
  )


# ------------------------------------------------------------
# 26. Identify constant-AUC eligible compounds
# ------------------------------------------------------------

depmap_eligible_but_untestable <-
  depmap_drug_results |>
  dplyr::filter(
    n >=
      depmap_min_n_per_drug,
    is.na(
      spearman_rho
    )
  )


# ------------------------------------------------------------
# 27. Secondary linear regression analysis
# ------------------------------------------------------------

depmap_linear_results_list <-
  lapply(
    depmap_drug_names,
    function(drug_name) {
      
      auc_values <-
        depmap_prism_luad[[drug_name]]
      
      keep <-
        !is.na(
          auc_values
        ) &
        !is.na(
          depmap_luad_scores$MPCDS_z
        )
      
      n_available <-
        sum(
          keep
        )
      
      if (
        n_available <
        depmap_min_n_per_drug
      ) {
        
        return(
          data.frame(
            
            drug_full_name =
              drug_name,
            
            linear_n =
              n_available,
            
            linear_beta =
              NA_real_,
            
            linear_SE =
              NA_real_,
            
            linear_p =
              NA_real_,
            
            stringsAsFactors =
              FALSE
          )
        )
      }
      
      if (
        length(
          unique(
            auc_values[
              keep
            ]
          )
        ) < 2
      ) {
        
        return(
          data.frame(
            
            drug_full_name =
              drug_name,
            
            linear_n =
              n_available,
            
            linear_beta =
              NA_real_,
            
            linear_SE =
              NA_real_,
            
            linear_p =
              NA_real_,
            
            stringsAsFactors =
              FALSE
          )
        )
      }
      
      linear_data <-
        data.frame(
          
          AUC =
            auc_values[
              keep
            ],
          
          MPCDS_z =
            depmap_luad_scores$MPCDS_z[
              keep
            ]
        )
      
      fit <-
        lm(
          AUC ~ MPCDS_z,
          data =
            linear_data
        )
      
      fit_summary <-
        summary(
          fit
        )$coefficients
      
      data.frame(
        
        drug_full_name =
          drug_name,
        
        linear_n =
          n_available,
        
        linear_beta =
          unname(
            fit_summary[
              "MPCDS_z",
              "Estimate"
            ]
          ),
        
        linear_SE =
          unname(
            fit_summary[
              "MPCDS_z",
              "Std. Error"
            ]
          ),
        
        linear_p =
          unname(
            fit_summary[
              "MPCDS_z",
              "Pr(>|t|)"
            ]
          ),
        
        stringsAsFactors =
          FALSE
      )
    }
  )


depmap_linear_results <-
  dplyr::bind_rows(
    depmap_linear_results_list
  )


# ------------------------------------------------------------
# 28. Linear-model BH correction
# ------------------------------------------------------------

depmap_linear_results$linear_FDR <-
  NA_real_


depmap_linear_tested <-
  !is.na(
    depmap_linear_results$linear_p
  )


depmap_linear_results$linear_FDR[
  depmap_linear_tested
] <-
  p.adjust(
    depmap_linear_results$linear_p[
      depmap_linear_tested
    ],
    method =
      "BH"
  )


# ------------------------------------------------------------
# 29. Combine primary and secondary results
# ------------------------------------------------------------

depmap_drug_results_combined <-
  depmap_drug_results |>
  dplyr::left_join(
    depmap_linear_results,
    by =
      "drug_full_name"
  )


# ------------------------------------------------------------
# 30. Spearman-linear agreement
# ------------------------------------------------------------

depmap_direction_agreement <-
  depmap_drug_results_combined |>
  dplyr::filter(
    !is.na(
      spearman_rho
    ),
    !is.na(
      linear_beta
    )
  ) |>
  dplyr::mutate(
    same_direction =
      sign(
        spearman_rho
      ) ==
      sign(
        linear_beta
      )
  )


depmap_direction_agreement_rate <-
  mean(
    depmap_direction_agreement$same_direction
  )


depmap_effect_concordance <-
  cor(
    depmap_direction_agreement$spearman_rho,
    depmap_direction_agreement$linear_beta,
    method =
      "spearman"
  )


# ------------------------------------------------------------
# 31. Exploratory nominal results
# ------------------------------------------------------------

depmap_nominal_p05 <-
  depmap_drug_results_combined |>
  dplyr::filter(
    !is.na(
      p_value
    ),
    p_value < 0.05
  ) |>
  dplyr::arrange(
    p_value
  )


depmap_nominal_p001 <-
  depmap_drug_results_combined |>
  dplyr::filter(
    !is.na(
      p_value
    ),
    p_value < 0.001
  ) |>
  dplyr::arrange(
    p_value
  )


depmap_top_sensitivity <-
  depmap_drug_results_tested |>
  dplyr::arrange(
    spearman_rho
  ) |>
  dplyr::slice_head(
    n = 20
  )


depmap_top_resistance <-
  depmap_drug_results_tested |>
  dplyr::arrange(
    dplyr::desc(
      spearman_rho
    )
  ) |>
  dplyr::slice_head(
    n = 20
  )


# ------------------------------------------------------------
# 32. Analysis summary
# ------------------------------------------------------------

n_total_drugs <-
  length(
    depmap_drug_names
  )


n_tested_spearman <-
  sum(
    !is.na(
      depmap_drug_results$spearman_rho
    )
  )


n_tested_linear <-
  sum(
    !is.na(
      depmap_linear_results$linear_beta
    )
  )


n_constant_auc <-
  nrow(
    depmap_eligible_but_untestable
  )


n_spearman_fdr05 <-
  sum(
    depmap_drug_results$FDR <
      0.05,
    na.rm = TRUE
  )


n_spearman_fdr10 <-
  sum(
    depmap_drug_results$FDR <
      0.10,
    na.rm = TRUE
  )


n_linear_fdr05 <-
  sum(
    depmap_linear_results$linear_FDR <
      0.05,
    na.rm = TRUE
  )


n_linear_fdr10 <-
  sum(
    depmap_linear_results$linear_FDR <
      0.10,
    na.rm = TRUE
  )


n_nominal_p05 <-
  sum(
    depmap_drug_results$p_value <
      0.05,
    na.rm = TRUE
  )


n_nominal_p01 <-
  sum(
    depmap_drug_results$p_value <
      0.01,
    na.rm = TRUE
  )


n_nominal_p001 <-
  sum(
    depmap_drug_results$p_value <
      0.001,
    na.rm = TRUE
  )


min_spearman_fdr <-
  min(
    depmap_drug_results$FDR,
    na.rm = TRUE
  )


min_linear_fdr <-
  min(
    depmap_linear_results$linear_FDR,
    na.rm = TRUE
  )


depmap_analysis_summary <-
  data.frame(
    
    metric =
      c(
        "DepMap LUAD models in metadata",
        "LUAD with default expression",
        "LUAD with PRISM",
        "Complete LUAD expression plus PRISM",
        "MPCDS genes available",
        "Total PRISM compounds",
        "Compounds with n >= 30",
        "Eligible constant-AUC compounds",
        "Spearman-tested compounds",
        "Spearman nominal p < 0.05",
        "Spearman nominal p < 0.01",
        "Spearman nominal p < 0.001",
        "Spearman FDR < 0.05",
        "Spearman FDR < 0.10",
        "Minimum Spearman FDR",
        "Linear-tested compounds",
        "Linear FDR < 0.05",
        "Linear FDR < 0.10",
        "Minimum linear FDR",
        "Spearman-linear direction agreement",
        "Spearman-linear effect concordance"
      ),
    
    value =
      c(
        nrow(
          depmap_luad_models
        ),
        n_luad_expression,
        n_luad_prism,
        length(
          depmap_luad_final_ids
        ),
        length(
          depmap_mpcds_genes
        ),
        n_total_drugs,
        n_drugs_n30,
        n_constant_auc,
        n_tested_spearman,
        n_nominal_p05,
        n_nominal_p01,
        n_nominal_p001,
        n_spearman_fdr05,
        n_spearman_fdr10,
        min_spearman_fdr,
        n_tested_linear,
        n_linear_fdr05,
        n_linear_fdr10,
        min_linear_fdr,
        depmap_direction_agreement_rate,
        depmap_effect_concordance
      ),
    
    stringsAsFactors =
      FALSE
  )


# ------------------------------------------------------------
# 33. Save processed data
# ------------------------------------------------------------

write.csv(
  depmap_luad_scores,
  paste0(
    "02_processed_data/drug_sensitivity/",
    "DepMap_LUAD_MPCDS_scores.csv"
  ),
  row.names = FALSE
)


write.csv(
  depmap_luad_final_metadata,
  paste0(
    "02_processed_data/drug_sensitivity/",
    "DepMap_LUAD_final_model_metadata.csv"
  ),
  row.names = FALSE
)


write.csv(
  depmap_mpcds_column_map,
  paste0(
    "02_processed_data/drug_sensitivity/",
    "DepMap_MPCDS_47_gene_column_mapping.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 34. Save result tables
# ------------------------------------------------------------

write.csv(
  depmap_drug_results_combined,
  paste0(
    "04_results/drug_sensitivity/",
    "DepMap_PRISM_MPCDS_all_drug_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  depmap_nominal_p05,
  paste0(
    "04_results/drug_sensitivity/",
    "DepMap_PRISM_MPCDS_nominal_p_lt_0.05.csv"
  ),
  row.names = FALSE
)


write.csv(
  depmap_nominal_p001,
  paste0(
    "04_results/drug_sensitivity/",
    "DepMap_PRISM_MPCDS_nominal_p_lt_0.001.csv"
  ),
  row.names = FALSE
)


write.csv(
  depmap_top_sensitivity,
  paste0(
    "04_results/drug_sensitivity/",
    "DepMap_PRISM_top_20_sensitivity_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  depmap_top_resistance,
  paste0(
    "04_results/drug_sensitivity/",
    "DepMap_PRISM_top_20_resistance_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  depmap_eligible_but_untestable,
  paste0(
    "04_results/drug_sensitivity/",
    "DepMap_PRISM_untestable_constant_AUC_drugs.csv"
  ),
  row.names = FALSE
)


write.csv(
  depmap_analysis_summary,
  paste0(
    "04_results/drug_sensitivity/",
    "DepMap_PRISM_MPCDS_analysis_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 35. Verify saved outputs
# ------------------------------------------------------------

depmap_expected_outputs <-
  c(
    paste0(
      "02_processed_data/drug_sensitivity/",
      "DepMap_LUAD_MPCDS_scores.csv"
    ),
    paste0(
      "02_processed_data/drug_sensitivity/",
      "DepMap_LUAD_final_model_metadata.csv"
    ),
    paste0(
      "02_processed_data/drug_sensitivity/",
      "DepMap_MPCDS_47_gene_column_mapping.csv"
    ),
    paste0(
      "04_results/drug_sensitivity/",
      "DepMap_PRISM_MPCDS_all_drug_associations.csv"
    ),
    paste0(
      "04_results/drug_sensitivity/",
      "DepMap_PRISM_MPCDS_nominal_p_lt_0.05.csv"
    ),
    paste0(
      "04_results/drug_sensitivity/",
      "DepMap_PRISM_MPCDS_nominal_p_lt_0.001.csv"
    ),
    paste0(
      "04_results/drug_sensitivity/",
      "DepMap_PRISM_top_20_sensitivity_associations.csv"
    ),
    paste0(
      "04_results/drug_sensitivity/",
      "DepMap_PRISM_top_20_resistance_associations.csv"
    ),
    paste0(
      "04_results/drug_sensitivity/",
      "DepMap_PRISM_untestable_constant_AUC_drugs.csv"
    ),
    paste0(
      "04_results/drug_sensitivity/",
      "DepMap_PRISM_MPCDS_analysis_summary.csv"
    ),
    paste0(
      "04_results/drug_sensitivity/",
      "DepMap_PRISM_input_file_metadata.csv"
    )
  )


depmap_output_verified <-
  file.exists(
    depmap_expected_outputs
  )


# ------------------------------------------------------------
# 36. Regression checks against established PRISM analysis
# ------------------------------------------------------------

stopifnot(
  length(
    depmap_luad_final_ids
  ) == 50,
  length(
    depmap_mpcds_genes
  ) == 47,
  n_total_drugs == 1482,
  n_drugs_n30 == 1442,
  n_constant_auc == 2,
  n_tested_spearman == 1440,
  n_nominal_p05 == 121,
  n_nominal_p01 == 29,
  n_nominal_p001 == 5,
  n_spearman_fdr05 == 0,
  n_spearman_fdr10 == 0,
  abs(
    min_spearman_fdr -
      0.1893165
  ) <
    1e-5
)

cat(
  "\nEstablished PRISM regression checks: PASSED\n"
)


# ------------------------------------------------------------
# 37. Final reproducibility summary
# ------------------------------------------------------------

# ------------------------------------------------------------

cat(
  "\n",
  paste(
    rep(
      "=",
      60
    ),
    collapse = ""
  ),
  "\n"
)


cat(
  "SCRIPT 24 PRISM PHARMACOGENOMIC ANALYSIS COMPLETE\n"
)


cat(
  paste(
    rep(
      "=",
      60
    ),
    collapse = ""
  ),
  "\n"
)


cat(
  "DepMap LUAD models in metadata:",
  nrow(
    depmap_luad_models
  ),
  "\n"
)


cat(
  "LUAD with default expression:",
  n_luad_expression,
  "\n"
)


cat(
  "LUAD with PRISM:",
  n_luad_prism,
  "\n"
)


cat(
  "Complete LUAD expression + PRISM:",
  length(
    depmap_luad_final_ids
  ),
  "\n"
)


cat(
  "MPCDS genes available:",
  length(
    depmap_mpcds_genes
  ),
  "/ 47\n"
)


cat(
  "Total PRISM compounds:",
  n_total_drugs,
  "\n"
)


cat(
  "Compounds with n >= 30:",
  n_drugs_n30,
  "\n"
)


cat(
  "Constant-AUC eligible compounds:",
  n_constant_auc,
  "\n"
)


cat(
  "Spearman-tested compounds:",
  n_tested_spearman,
  "\n"
)


cat(
  "Spearman nominal p < 0.05:",
  n_nominal_p05,
  "\n"
)


cat(
  "Spearman nominal p < 0.01:",
  n_nominal_p01,
  "\n"
)


cat(
  "Spearman nominal p < 0.001:",
  n_nominal_p001,
  "\n"
)


cat(
  "Spearman FDR < 0.05:",
  n_spearman_fdr05,
  "\n"
)


cat(
  "Spearman FDR < 0.10:",
  n_spearman_fdr10,
  "\n"
)


cat(
  "Minimum Spearman FDR:",
  round(
    min_spearman_fdr,
    7
  ),
  "\n"
)


cat(
  "Linear-tested compounds:",
  n_tested_linear,
  "\n"
)


cat(
  "Linear FDR < 0.05:",
  n_linear_fdr05,
  "\n"
)


cat(
  "Linear FDR < 0.10:",
  n_linear_fdr10,
  "\n"
)


cat(
  "Minimum linear FDR:",
  round(
    min_linear_fdr,
    7
  ),
  "\n"
)


cat(
  "Spearman-linear direction agreement:",
  round(
    100 *
      depmap_direction_agreement_rate,
    2
  ),
  "%\n"
)


cat(
  "Spearman-linear effect concordance:",
  round(
    depmap_effect_concordance,
    4
  ),
  "\n"
)


cat(
  "Saved output files:",
  length(
    depmap_expected_outputs
  ),
  "\n"
)


cat(
  "All expected outputs verified:",
  all(
    depmap_output_verified
  ),
  "\n"
)


cat(
  paste(
    rep(
      "=",
      60
    ),
    collapse = ""
  ),
  "\n"
)