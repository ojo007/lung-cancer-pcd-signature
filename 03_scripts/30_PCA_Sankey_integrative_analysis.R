# ============================================================
# SCRIPT 30
# Sankey and final integrative visualizations
#
# TCGA-LUAD
#
# Analyses:
#   1. 47-gene discovery-cohort PCA
#   2. Clinical-stage harmonization
#   3. PCD subtype -> MPCDS -> AJCC stage Sankey
#   4. Combined integrative supplementary figure
# ============================================================


source("03_scripts/00_project_config.R")


cat(
  "\n========================================\n",
  "SCRIPT 30\n",
  "PCA AND SANKEY INTEGRATIVE ANALYSIS\n",
  "========================================\n"
)


# ============================================================
# 1. Package requirements
# ============================================================

# PCA uses base R (stats::prcomp).
# Final plotting is intentionally deferred to Script 32.


# ============================================================
# 2. Input files
# ============================================================


subtype_file_s35 <-
  "04_results/clustering/TCGA_LUAD_PCD_cluster_assignments.csv"


mpcds_file_s35 <-
  "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv"


expression_file_s35 <-
  "02_processed_data/immunotherapy/TCGA_LUAD_patient_protein_coding_TPM.rds"


clinical_file_s35 <-
  "02_processed_data/clinical/TCGA_LUAD_clinical_clean.csv"


ridge_file_s35 <-
  "04_results/mpcds/TCGA_LUAD_MPCDS_ridge_coefficients.csv"


input_files_s35 <-
  c(
    subtype_file_s35,
    mpcds_file_s35,
    expression_file_s35,
    clinical_file_s35,
    ridge_file_s35
  )


cat(
  "\nInput-file check:\n"
)


print(
  data.frame(
    file =
      input_files_s35,
    exists =
      file.exists(
        input_files_s35
      ),
    size_bytes =
      file.info(
        input_files_s35
      )$size,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)


stopifnot(
  all(
    file.exists(
      input_files_s35
    )
  )
)


# ============================================================
# 3. Output directories
# ============================================================


results_dir_s35 <-
  "04_results/final_summary"


processed_dir_s35 <-
  "02_processed_data/final_summary"



dir.create(
  results_dir_s35,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  processed_dir_s35,
  recursive = TRUE,
  showWarnings = FALSE
)



# ============================================================
# 4. Load inputs
# ============================================================


subtype_s35 <-
  read.csv(
    subtype_file_s35,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


mpcds_s35 <-
  read.csv(
    mpcds_file_s35,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


expression_s35 <-
  as.matrix(
    readRDS(
      expression_file_s35
    )
  )


clinical_s35 <-
  read.csv(
    clinical_file_s35,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


ridge_s35 <-
  read.csv(
    ridge_file_s35,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


# ============================================================
# 5. Input QC
# ============================================================


cat(
  "\n========================================\n",
  "INPUT QC\n",
  "========================================\n"
)


cat(
  "\nSubtype dimensions:\n"
)


print(
  dim(
    subtype_s35
  )
)


cat(
  "\nMPCDS dimensions:\n"
)


print(
  dim(
    mpcds_s35
  )
)


cat(
  "\nExpression dimensions:\n"
)


print(
  dim(
    expression_s35
  )
)


cat(
  "\nClinical dimensions:\n"
)


print(
  dim(
    clinical_s35
  )
)


cat(
  "\nRidge coefficient dimensions:\n"
)


print(
  dim(
    ridge_s35
  )
)


stopifnot(
  identical(
    dim(
      subtype_s35
    ),
    c(
      517L,
      2L
    )
  )
)


stopifnot(
  nrow(
    mpcds_s35
  ) ==
    504L
)


stopifnot(
  identical(
    dim(
      expression_s35
    ),
    c(
      19938L,
      517L
    )
  )
)


stopifnot(
  nrow(
    ridge_s35
  ) ==
    47L
)


stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster"
    ) %in%
      colnames(
        subtype_s35
      )
  )
)


stopifnot(
  all(
    c(
      "patient_id",
      "OS_time",
      "OS_status",
      "MPCDS",
      "MPCDS_z",
      "risk_group"
    ) %in%
      colnames(
        mpcds_s35
      )
  )
)


stopifnot(
  all(
    c(
      "submitter_id",
      "ajcc_pathologic_stage"
    ) %in%
      colnames(
        clinical_s35
      )
  )
)


stopifnot(
  all(
    c(
      "gene",
      "coefficient"
    ) %in%
      colnames(
        ridge_s35
      )
  )
)


stopifnot(
  !any(
    duplicated(
      subtype_s35[[
        "patient_id"
      ]]
    )
  )
)


stopifnot(
  !any(
    duplicated(
      mpcds_s35[[
        "patient_id"
      ]]
    )
  )
)


stopifnot(
  !any(
    duplicated(
      rownames(
        expression_s35
      )
    )
  )
)


stopifnot(
  !any(
    duplicated(
      colnames(
        expression_s35
      )
    )
  )
)


cat(
  "\nSubtype counts:\n"
)


print(
  table(
    subtype_s35[[
      "PCD_cluster"
    ]]
  )
)


stopifnot(
  sum(
    subtype_s35[[
      "PCD_cluster"
    ]] ==
      "PCD_C1"
  ) ==
    284L
)


stopifnot(
  sum(
    subtype_s35[[
      "PCD_cluster"
    ]] ==
      "PCD_C2"
  ) ==
    233L
)


# ============================================================
# 6. Build core 504-patient integrative table
# ============================================================


integrative_core_s35 <-
  merge(
    subtype_s35[
      ,
      c(
        "patient_id",
        "PCD_cluster"
      ),
      drop = FALSE
    ],
    mpcds_s35[
      ,
      c(
        "patient_id",
        "OS_time",
        "OS_status",
        "MPCDS",
        "MPCDS_z",
        "risk_group"
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
    integrative_core_s35
  ) ==
    504L
)


write.csv(
  integrative_core_s35,
  file.path(
    results_dir_s35,
    "TCGA_LUAD_integrative_visualization_core_patients.csv"
  ),
  row.names = FALSE
)


saveRDS(
  integrative_core_s35,
  file.path(
    processed_dir_s35,
    "TCGA_LUAD_integrative_visualization_core_patients.rds"
  )
)


# ============================================================
# 7. Define 47 PCD candidate genes
#
# Gene identities come from the locked LUAD ridge coefficient
# table. Coefficient magnitudes are not used in the PCA.
# ============================================================


pcd47_genes_s35 <-
  unique(
    trimws(
      as.character(
        ridge_s35[[
          "gene"
        ]]
      )
    )
  )


pcd47_genes_s35 <-
  pcd47_genes_s35[
    !is.na(
      pcd47_genes_s35
    ) &
      nzchar(
        pcd47_genes_s35
      )
  ]


stopifnot(
  length(
    pcd47_genes_s35
  ) ==
    47L
)


pcd47_missing_s35 <-
  setdiff(
    pcd47_genes_s35,
    rownames(
      expression_s35
    )
  )


cat(
  "\n47-gene expression match:\n"
)


print(
  data.frame(
    requested =
      47L,
    present =
      47L -
      length(
        pcd47_missing_s35
      ),
    missing =
      length(
        pcd47_missing_s35
      )
  ),
  row.names = FALSE
)


stopifnot(
  length(
    pcd47_missing_s35
  ) ==
    0L
)


# ============================================================
# 8. Extract 47-gene expression
# ============================================================


pca_expression_tpm_s35 <-
  expression_s35[
    pcd47_genes_s35,
    subtype_s35[[
      "patient_id"
    ]],
    drop = FALSE
  ]


stopifnot(
  identical(
    dim(
      pca_expression_tpm_s35
    ),
    c(
      47L,
      517L
    )
  )
)


stopifnot(
  !anyNA(
    pca_expression_tpm_s35
  )
)


stopifnot(
  all(
    is.finite(
      pca_expression_tpm_s35
    )
  )
)


# ============================================================
# 9. Log2(TPM + 1) and gene-wise standardization
# ============================================================


pca_log_expression_s35 <-
  log2(
    pca_expression_tpm_s35 +
      1
  )


gene_sd_s35 <-
  apply(
    pca_log_expression_s35,
    1,
    sd
  )


zero_variance_genes_s35 <-
  names(
    gene_sd_s35
  )[
    !is.finite(
      gene_sd_s35
    ) |
      gene_sd_s35 ==
      0
  ]


stopifnot(
  length(
    zero_variance_genes_s35
  ) ==
    0L
)


pca_z_expression_s35 <-
  t(
    scale(
      t(
        pca_log_expression_s35
      ),
      center =
        TRUE,
      scale =
        TRUE
    )
  )


stopifnot(
  !anyNA(
    pca_z_expression_s35
  )
)


stopifnot(
  all(
    is.finite(
      pca_z_expression_s35
    )
  )
)


# ============================================================
# 10. PCA
# ============================================================


pca_fit_s35 <-
  prcomp(
    t(
      pca_z_expression_s35
    ),
    center =
      FALSE,
    scale. =
      FALSE
  )


pca_variance_s35 <-
  pca_fit_s35$sdev^2


pca_variance_percent_s35 <-
  100 *
  pca_variance_s35 /
  sum(
    pca_variance_s35
  )


pca_variance_table_s35 <-
  data.frame(
    PC =
      paste0(
        "PC",
        seq_along(
          pca_variance_percent_s35
        )
      ),
    variance_percent =
      pca_variance_percent_s35,
    cumulative_variance_percent =
      cumsum(
        pca_variance_percent_s35
      ),
    stringsAsFactors = FALSE
  )


cat(
  "\n========================================\n",
  "PCA VARIANCE EXPLAINED\n",
  "========================================\n"
)


print(
  utils::head(
    pca_variance_table_s35,
    10
  ),
  row.names = FALSE
)


# ============================================================
# 11. PCA scores
# ============================================================


pca_scores_s35 <-
  data.frame(
    patient_id =
      rownames(
        pca_fit_s35$x
      ),
    PC1 =
      pca_fit_s35$x[
        ,
        1
      ],
    PC2 =
      pca_fit_s35$x[
        ,
        2
      ],
    PC3 =
      pca_fit_s35$x[
        ,
        3
      ],
    stringsAsFactors = FALSE
  )


pca_scores_s35 <-
  merge(
    pca_scores_s35,
    subtype_s35[
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
  nrow(
    pca_scores_s35
  ) ==
    517L
)


stopifnot(
  !anyNA(
    pca_scores_s35[[
      "PCD_cluster"
    ]]
  )
)


# ============================================================
# 12. PCA loadings
# ============================================================


pca_loadings_s35 <-
  data.frame(
    gene =
      rownames(
        pca_fit_s35$rotation
      ),
    PC1_loading =
      pca_fit_s35$rotation[
        ,
        1
      ],
    PC2_loading =
      pca_fit_s35$rotation[
        ,
        2
      ],
    PC3_loading =
      pca_fit_s35$rotation[
        ,
        3
      ],
    stringsAsFactors = FALSE
  )


pca_loadings_s35$abs_PC1_loading <-
  abs(
    pca_loadings_s35[[
      "PC1_loading"
    ]]
  )


pca_loadings_s35$abs_PC2_loading <-
  abs(
    pca_loadings_s35[[
      "PC2_loading"
    ]]
  )


# ============================================================
# 13. PCA figure generation
# ============================================================

# Final PCA figure generation is deferred to Script 32.


# ============================================================
# 14. Save PCA tables
# ============================================================


pca_scores_file_s35 <-
  file.path(
    results_dir_s35,
    "TCGA_LUAD_47gene_PCD_subtype_PCA_scores.csv"
  )


pca_loadings_file_s35 <-
  file.path(
    results_dir_s35,
    "TCGA_LUAD_47gene_PCD_subtype_PCA_loadings.csv"
  )


pca_variance_file_s35 <-
  file.path(
    results_dir_s35,
    "TCGA_LUAD_47gene_PCD_subtype_PCA_variance.csv"
  )


write.csv(
  pca_scores_s35,
  pca_scores_file_s35,
  row.names = FALSE
)


write.csv(
  pca_loadings_s35,
  pca_loadings_file_s35,
  row.names = FALSE
)


write.csv(
  pca_variance_table_s35,
  pca_variance_file_s35,
  row.names = FALSE
)


saveRDS(
  pca_fit_s35,
  file.path(
    processed_dir_s35,
    "TCGA_LUAD_47gene_PCD_subtype_PCA_fit.rds"
  )
)


# ============================================================
# 15. Harmonize AJCC stage
# ============================================================


clinical_stage_s35 <-
  data.frame(
    patient_id =
      as.character(
        clinical_s35[[
          "submitter_id"
        ]]
      ),
    stage_raw =
      as.character(
        clinical_s35[[
          "ajcc_pathologic_stage"
        ]]
      ),
    stringsAsFactors = FALSE
  )


stopifnot(
  !any(
    duplicated(
      clinical_stage_s35[[
        "patient_id"
      ]]
    )
  )
)


stage_normalized_s35 <-
  toupper(
    trimws(
      clinical_stage_s35[[
        "stage_raw"
      ]]
    )
  )


clinical_stage_s35$stage_group <-
  NA_character_


clinical_stage_s35$stage_group[
  grepl(
    "^STAGE I($|A|B|C)",
    stage_normalized_s35
  )
] <-
  "Stage I"


clinical_stage_s35$stage_group[
  grepl(
    "^STAGE II($|A|B|C)",
    stage_normalized_s35
  )
] <-
  "Stage II"


clinical_stage_s35$stage_group[
  grepl(
    "^STAGE III($|A|B|C)",
    stage_normalized_s35
  )
] <-
  "Stage III"


clinical_stage_s35$stage_group[
  grepl(
    "^STAGE IV($|A|B|C)",
    stage_normalized_s35
  )
] <-
  "Stage IV"


# ============================================================
# 16. Stage harmonization audit
# ============================================================


stage_harmonization_audit_s35 <-
  unique(
    clinical_stage_s35[
      ,
      c(
        "stage_raw",
        "stage_group"
      ),
      drop = FALSE
    ]
  )


stage_harmonization_audit_s35 <-
  stage_harmonization_audit_s35[
    order(
      stage_harmonization_audit_s35[[
        "stage_raw"
      ]]
    ),
    ,
    drop = FALSE
  ]


rownames(
  stage_harmonization_audit_s35
) <-
  NULL


unharmonized_stage_s35 <-
  unique(
    clinical_stage_s35[[
      "stage_raw"
    ]][
      !is.na(
        clinical_stage_s35[[
          "stage_raw"
        ]]
      ) &
        nzchar(
          trimws(
            clinical_stage_s35[[
              "stage_raw"
            ]]
          )
        ) &
        is.na(
          clinical_stage_s35[[
            "stage_group"
          ]]
        )
    ]
  )


stopifnot(
  length(
    unharmonized_stage_s35
  ) ==
    0L
)


write.csv(
  stage_harmonization_audit_s35,
  file.path(
    results_dir_s35,
    "TCGA_LUAD_AJCC_stage_harmonization_audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. Merge stage into integrative cohort
# ============================================================


integrative_stage_s35 <-
  merge(
    integrative_core_s35,
    clinical_stage_s35,
    by =
      "patient_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )


stopifnot(
  nrow(
    integrative_stage_s35
  ) ==
    504L
)


write.csv(
  integrative_stage_s35,
  file.path(
    results_dir_s35,
    "TCGA_LUAD_integrative_visualization_with_stage.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. Complete-case Sankey cohort
# ============================================================


sankey_complete_s35 <-
  integrative_stage_s35[
    !is.na(
      integrative_stage_s35[[
        "PCD_cluster"
      ]]
    ) &
      !is.na(
        integrative_stage_s35[[
          "risk_group"
        ]]
      ) &
      !is.na(
        integrative_stage_s35[[
          "stage_group"
        ]]
      ),
    ,
    drop = FALSE
  ]


stopifnot(
  nrow(
    sankey_complete_s35
  ) ==
    496L
)


write.csv(
  sankey_complete_s35,
  file.path(
    results_dir_s35,
    "TCGA_LUAD_Sankey_complete_case_patients.csv"
  ),
  row.names = FALSE
)


saveRDS(
  sankey_complete_s35,
  file.path(
    processed_dir_s35,
    "TCGA_LUAD_Sankey_complete_case_patients.rds"
  )
)


# ============================================================
# 19. Three-way flow counts
# ============================================================


sankey_flow_counts_s35 <-
  as.data.frame(
    table(
      PCD_cluster =
        sankey_complete_s35[[
          "PCD_cluster"
        ]],
      risk_group =
        sankey_complete_s35[[
          "risk_group"
        ]],
      stage_group =
        sankey_complete_s35[[
          "stage_group"
        ]]
    ),
    stringsAsFactors = FALSE
  )


sankey_flow_counts_s35 <-
  sankey_flow_counts_s35[
    sankey_flow_counts_s35[[
      "Freq"
    ]] >
      0,
    ,
    drop = FALSE
  ]


names(
  sankey_flow_counts_s35
)[
  names(
    sankey_flow_counts_s35
  ) ==
    "Freq"
] <-
  "n_patients"


stopifnot(
  sum(
    sankey_flow_counts_s35[[
      "n_patients"
    ]]
  ) ==
    496L
)


write.csv(
  sankey_flow_counts_s35,
  file.path(
    results_dir_s35,
    "TCGA_LUAD_Sankey_three_way_flow_counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. Prepare Sankey plot data
# ============================================================


sankey_plot_data_s35 <-
  sankey_complete_s35


sankey_plot_data_s35$PCD_cluster <-
  factor(
    sankey_plot_data_s35[[
      "PCD_cluster"
    ]],
    levels =
      c(
        "PCD_C1",
        "PCD_C2"
      )
  )


sankey_plot_data_s35$risk_group <-
  factor(
    sankey_plot_data_s35[[
      "risk_group"
    ]],
    levels =
      c(
        "Low MPCDS",
        "High MPCDS"
      )
  )


sankey_plot_data_s35$stage_group <-
  factor(
    sankey_plot_data_s35[[
      "stage_group"
    ]],
    levels =
      c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      )
  )


stopifnot(
  ggalluvial::is_alluvia_form(
    sankey_plot_data_s35,
    axes =
      2:4,
    silent =
      TRUE
  )
)


# ============================================================
# 21. Sankey axis summary
# ============================================================


sankey_axis_summary_s35 <-
  rbind(
    
    data.frame(
      axis =
        "PCD subtype",
      category =
        names(
          table(
            sankey_plot_data_s35[[
              "PCD_cluster"
            ]]
          )
        ),
      n =
        as.integer(
          table(
            sankey_plot_data_s35[[
              "PCD_cluster"
            ]]
          )
        ),
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      axis =
        "MPCDS group",
      category =
        names(
          table(
            sankey_plot_data_s35[[
              "risk_group"
            ]]
          )
        ),
      n =
        as.integer(
          table(
            sankey_plot_data_s35[[
              "risk_group"
            ]]
          )
        ),
      stringsAsFactors = FALSE
    ),
    
    data.frame(
      axis =
        "AJCC stage",
      category =
        names(
          table(
            sankey_plot_data_s35[[
              "stage_group"
            ]]
          )
        ),
      n =
        as.integer(
          table(
            sankey_plot_data_s35[[
              "stage_group"
            ]]
          )
        ),
      stringsAsFactors = FALSE
    )
  )


sankey_axis_summary_s35$percent <-
  100 *
  sankey_axis_summary_s35[[
    "n"
  ]] /
  nrow(
    sankey_plot_data_s35
  )


node_count_lookup_s35 <-
  setNames(
    sankey_axis_summary_s35[[
      "n"
    ]],
    sankey_axis_summary_s35[[
      "category"
    ]]
  )


node_label_s35 <-
  function(x_s35) {
    
    x_chr_s35 <-
      as.character(
        x_s35
      )
    
    
    paste0(
      x_chr_s35,
      "\n(n=",
      node_count_lookup_s35[
        x_chr_s35
      ],
      ")"
    )
  }


# ============================================================
# 22. Sankey figure generation
# ============================================================

# Final Sankey/alluvial figure generation is deferred to Script 32.


# ============================================================
# 23. Sankey flow summary
# ============================================================


sankey_flow_summary_s35 <-
  sankey_flow_counts_s35


sankey_flow_summary_s35$percent_of_total <-
  100 *
  sankey_flow_summary_s35[[
    "n_patients"
  ]] /
  sum(
    sankey_flow_summary_s35[[
      "n_patients"
    ]]
  )


sankey_flow_summary_s35$percent_within_subtype <-
  NA_real_


for (
  subtype_i_s35 in unique(
    sankey_flow_summary_s35[[
      "PCD_cluster"
    ]]
  )
) {
  
  idx_s35 <-
    sankey_flow_summary_s35[[
      "PCD_cluster"
    ]] ==
    subtype_i_s35
  
  
  subtype_total_s35 <-
    sum(
      sankey_flow_summary_s35[[
        "n_patients"
      ]][
        idx_s35
      ]
    )
  
  
  sankey_flow_summary_s35[[
    "percent_within_subtype"
  ]][
    idx_s35
  ] <-
    100 *
    sankey_flow_summary_s35[[
      "n_patients"
    ]][
      idx_s35
    ] /
    subtype_total_s35
}


sankey_axis_file_s35 <-
  file.path(
    results_dir_s35,
    "TCGA_LUAD_Sankey_axis_summary.csv"
  )


sankey_flow_summary_file_s35 <-
  file.path(
    results_dir_s35,
    "TCGA_LUAD_Sankey_flow_summary.csv"
  )


write.csv(
  sankey_axis_summary_s35,
  sankey_axis_file_s35,
  row.names = FALSE
)


write.csv(
  sankey_flow_summary_s35,
  sankey_flow_summary_file_s35,
  row.names = FALSE
)



# ============================================================
# 24. Combined figure generation
# ============================================================

# Combined PCA + Sankey visualization is deferred to Script 32.


# ============================================================
# 26. Final expected-result audit
# ============================================================


script35_audit_s35 <-
  data.frame(
    metric =
      c(
        "LUAD subtype patients",
        "PCA genes",
        "PCA PC1 variance percent",
        "PCA PC2 variance percent",
        "Sankey complete-case patients",
        "Sankey PCD_C1 patients",
        "Sankey PCD_C2 patients",
        "Sankey Low MPCDS patients",
        "Sankey High MPCDS patients",
        "Sankey Stage I patients",
        "Sankey Stage II patients",
        "Sankey Stage III patients",
        "Sankey Stage IV patients"
      ),
    
    observed =
      c(
        nrow(
          pca_scores_s35
        ),
        
        length(
          pcd47_genes_s35
        ),
        
        pca_variance_percent_s35[
          1
        ],
        
        pca_variance_percent_s35[
          2
        ],
        
        nrow(
          sankey_plot_data_s35
        ),
        
        sum(
          sankey_plot_data_s35[[
            "PCD_cluster"
          ]] ==
            "PCD_C1"
        ),
        
        sum(
          sankey_plot_data_s35[[
            "PCD_cluster"
          ]] ==
            "PCD_C2"
        ),
        
        sum(
          sankey_plot_data_s35[[
            "risk_group"
          ]] ==
            "Low MPCDS"
        ),
        
        sum(
          sankey_plot_data_s35[[
            "risk_group"
          ]] ==
            "High MPCDS"
        ),
        
        sum(
          sankey_plot_data_s35[[
            "stage_group"
          ]] ==
            "Stage I"
        ),
        
        sum(
          sankey_plot_data_s35[[
            "stage_group"
          ]] ==
            "Stage II"
        ),
        
        sum(
          sankey_plot_data_s35[[
            "stage_group"
          ]] ==
            "Stage III"
        ),
        
        sum(
          sankey_plot_data_s35[[
            "stage_group"
          ]] ==
            "Stage IV"
        )
      ),
    
    expected =
      c(
        517,
        47,
        15.697042,
        11.770018,
        496,
        270,
        226,
        247,
        249,
        270,
        120,
        81,
        25
      ),
    
    stringsAsFactors = FALSE
  )


script35_audit_s35$tolerance <-
  c(
    0,
    0,
    0.0001,
    0.0001,
    rep(
      0,
      9
    )
  )


script35_audit_s35$matches_expected <-
  abs(
    script35_audit_s35[[
      "observed"
    ]] -
      script35_audit_s35[[
        "expected"
      ]]
  ) <=
  script35_audit_s35[[
    "tolerance"
  ]]


cat(
  "\n========================================\n",
  "FINAL SCRIPT 30 AUDIT\n",
  "========================================\n"
)


print(
  script35_audit_s35,
  row.names = FALSE
)


stopifnot(
  all(
    script35_audit_s35[[
      "matches_expected"
    ]]
  )
)

cat(
  "\nEstablished PCA/Sankey regression checks: PASSED\n"
)


# ============================================================
# 27. Final output verification
# ============================================================


script35_final_files_s35 <-
  c(
    pca_scores_file_s35,
    pca_loadings_file_s35,
    pca_variance_file_s35,
    
    file.path(
      results_dir_s35,
      "TCGA_LUAD_AJCC_stage_harmonization_audit.csv"
    ),
    
    file.path(
      results_dir_s35,
      "TCGA_LUAD_Sankey_complete_case_patients.csv"
    ),
    
    file.path(
      results_dir_s35,
      "TCGA_LUAD_Sankey_three_way_flow_counts.csv"
    ),
    
    sankey_axis_file_s35,
    sankey_flow_summary_file_s35,
    
    file.path(
      processed_dir_s35,
      "TCGA_LUAD_47gene_PCD_subtype_PCA_fit.rds"
    )
  )


script35_final_file_check_s35 <-
  data.frame(
    file =
      script35_final_files_s35,
    exists =
      file.exists(
        script35_final_files_s35
      ),
    size_bytes =
      file.info(
        script35_final_files_s35
      )$size,
    stringsAsFactors = FALSE
  )


cat(
  "\n========================================\n",
  "FINAL OUTPUT FILE CHECK\n",
  "========================================\n"
)


print(
  script35_final_file_check_s35,
  row.names = FALSE
)


stopifnot(
  all(
    script35_final_file_check_s35[[
      "exists"
    ]]
  )
)


stopifnot(
  all(
    script35_final_file_check_s35[[
      "size_bytes"
    ]] >
      0
  )
)


# ============================================================
# 28. Completion marker
# ============================================================


cat(
  "\n========================================\n",
  "SCRIPT 30 PCA AND SANKEY ANALYSIS COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)