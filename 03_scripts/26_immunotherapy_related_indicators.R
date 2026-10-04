# ============================================================
# 26_immunotherapy_related_indicators.R
#
# MSc Lung Cancer Programmed Cell Death Project
#
# TCGA-LUAD:
# Immunotherapy-related indicator characterization of:
#
#   1. Continuous MPCDS
#   2. Expression-derived PCD subtypes
#
# Indicators:
#
#   - CD274 / PD-L1 expression
#   - Nonsynonymous mutation burden
#   - MANTIS microsatellite-instability score
#
# Important interpretation constraints:
#
#   - These are biomarker associations/context only.
#   - They do NOT demonstrate immunotherapy response.
#   - Mutation burden is NOT labelled formal TMB because
#     no callable-megabase denominator is available.
#   - Continuous MANTIS score is the primary MSI-related
#     analysis because MSI-H is extremely rare in LUAD.
# ============================================================


source("03_scripts/00_project_config.R")


# ============================================================
# 1. Setup
# ============================================================


options(
  stringsAsFactors = FALSE
)


required_packages <-
  c(
    "readxl"
  )


missing_packages <-
  required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]


if (
  length(
    missing_packages
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing required packages: ",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


library(readxl)


cat(
  "\n========================================\n",
  "SCRIPT 26: IMMUNOTHERAPY-RELATED INDICATORS\n",
  "========================================\n"
)


# ============================================================
# 2. Define input and output paths
# ============================================================


expression_processed_dir <-
  "02_processed_data/expression"


mutation_processed_dir <-
  "02_processed_data/mutation"


clustering_results_dir <-
  "04_results/clustering"


mutation_results_dir <-
  "04_results/mutation"


mpcds_results_dir <-
  "04_results/mpcds"


msi_raw_dir <-
  "01_raw_data/TCGA_LUAD/MSI"


immunotherapy_results_dir <-
  "04_results/mpcds/immunotherapy_indicators"



dir.create(
  immunotherapy_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# ============================================================
# 3. Define exact input files
# ============================================================


luad_tpm_file <-
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_TPM.rds"
  )


luad_annotation_file <-
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_gene_annotation.csv"
  )


luad_metadata_file <-
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_sample_metadata.csv"
  )


luad_removed_duplicates_file <-
  file.path(
    expression_processed_dir,
    "TCGA_LUAD_removed_duplicate_tumor_samples.csv"
  )


luad_subtype_file <-
  file.path(
    clustering_results_dir,
    "TCGA_LUAD_PCD_cluster_assignments.csv"
  )


luad_mpcds_file <-
  file.path(
    mpcds_results_dir,
    "TCGA_LUAD_MPCDS_patient_scores.csv"
  )


luad_mutation_burden_file <-
  file.path(
    mutation_results_dir,
    "TCGA_LUAD_nonsynonymous_mutation_burden.csv"
  )


msi_source_file <-
  file.path(
    msi_raw_dir,
    "ds_po.17.00073-1.xlsx"
  )


required_input_files <-
  c(
    luad_tpm_file,
    luad_annotation_file,
    luad_metadata_file,
    luad_removed_duplicates_file,
    luad_subtype_file,
    luad_mpcds_file,
    luad_mutation_burden_file,
    msi_source_file
  )


missing_input_files <-
  required_input_files[
    !file.exists(
      required_input_files
    )
  ]


if (
  length(
    missing_input_files
  ) >
  0
) {
  
  stop(
    paste0(
      "Missing required input files:\n",
      paste(
        missing_input_files,
        collapse = "\n"
      )
    )
  )
}


cat(
  "\nAll required input files are present.\n"
)


# ============================================================
# 4. Load canonical subtype and MPCDS datasets
# ============================================================


luad_subtypes <-
  read.csv(
    luad_subtype_file,
    check.names = FALSE
  )


luad_mpcds_section2 <-
  read.csv(
    luad_mpcds_file,
    check.names = FALSE
  )


stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster"
    ) %in%
      colnames(
        luad_subtypes
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
        luad_mpcds_section2
      )
  )
)


stopifnot(
  !any(
    duplicated(
      luad_subtypes$patient_id
    )
  )
)


stopifnot(
  !any(
    duplicated(
      luad_mpcds_section2$patient_id
    )
  )
)


cat(
  "\n========================================\n",
  "CANONICAL COHORTS\n",
  "========================================\n"
)


cat(
  "\nPCD subtype patients:\n"
)


print(
  nrow(
    luad_subtypes
  )
)


cat(
  "\nMPCDS patients:\n"
)


print(
  nrow(
    luad_mpcds_section2
  )
)


cat(
  "\nPCD subtype distribution:\n"
)


print(
  table(
    luad_subtypes$PCD_cluster
  )
)


# ============================================================
# 5. Prepare patient-level CD274 / PD-L1 expression
# ============================================================


# ------------------------------------------------------------
# 5.1 Load full LUAD TPM matrix, annotation, and metadata
# ------------------------------------------------------------


luad_tpm <-
  readRDS(
    luad_tpm_file
  )


luad_gene_annotation <-
  read.csv(
    luad_annotation_file,
    check.names = FALSE
  )


luad_expression_metadata <-
  read.csv(
    luad_metadata_file,
    check.names = FALSE
  )


luad_removed_duplicates <-
  read.csv(
    luad_removed_duplicates_file,
    check.names = FALSE
  )


stopifnot(
  is.matrix(
    luad_tpm
  ) ||
    is.array(
      luad_tpm
    )
)


stopifnot(
  identical(
    colnames(
      luad_tpm
    ),
    luad_expression_metadata$cases
  )
)


stopifnot(
  all(
    c(
      "gene_id",
      "gene_name",
      "gene_type",
      "ensembl_id"
    ) %in%
      colnames(
        luad_gene_annotation
      )
  )
)


# ------------------------------------------------------------
# 5.2 Identify CD274
# ------------------------------------------------------------


cd274_annotation <-
  luad_gene_annotation[
    luad_gene_annotation$gene_name ==
      "CD274" &
      luad_gene_annotation$gene_type ==
      "protein_coding",
    ,
    drop = FALSE
  ]


cat(
  "\n========================================\n",
  "CD274 ANNOTATION\n",
  "========================================\n"
)


print(
  cd274_annotation
)


stopifnot(
  nrow(
    cd274_annotation
  ) ==
    1
)


cd274_gene_id <-
  cd274_annotation$gene_id[1]


stopifnot(
  cd274_gene_id %in%
    rownames(
      luad_tpm
    )
)


# ------------------------------------------------------------
# 5.3 Attach CD274 TPM to sample metadata
# ------------------------------------------------------------


luad_cd274_metadata <-
  luad_expression_metadata


luad_cd274_metadata$CD274_TPM <-
  as.numeric(
    luad_tpm[
      cd274_gene_id,
      ,
      drop = TRUE
    ]
  )


stopifnot(
  nrow(
    luad_cd274_metadata
  ) ==
    ncol(
      luad_tpm
    )
)


stopifnot(
  !any(
    is.na(
      luad_cd274_metadata$CD274_TPM
    )
  )
)


# ------------------------------------------------------------
# 5.4 Restrict to Primary Tumor RNA files
# ------------------------------------------------------------


luad_cd274_primary <-
  luad_cd274_metadata[
    luad_cd274_metadata$sample_type ==
      "Primary Tumor",
    ,
    drop = FALSE
  ]


cat(
  "\n========================================\n",
  "PRIMARY-TUMOR CD274 DATA\n",
  "========================================\n"
)


cat(
  "\nPrimary Tumor RNA-file rows:\n"
)


print(
  nrow(
    luad_cd274_primary
  )
)


cat(
  "\nUnique Primary Tumor patients:\n"
)


print(
  length(
    unique(
      luad_cd274_primary$cases.submitter_id
    )
  )
)


cat(
  "\nUnique biological tumor samples:\n"
)


print(
  length(
    unique(
      luad_cd274_primary$sample.submitter_id
    )
  )
)


# ------------------------------------------------------------
# 5.5 Reproduce Script 11 biological-sample selection
#
# Script 11 first collapsed raw-count technical aliquots
# sharing sample.submitter_id.
#
# For patients with >1 distinct Primary Tumor biological
# sample, Script 11 retained the tumor sample with the
# largest collapsed raw-count library size.
#
# The removed biological samples are already recorded in:
# TCGA_LUAD_removed_duplicate_tumor_samples.csv
# ------------------------------------------------------------


removed_biological_sample_ids <-
  unique(
    luad_removed_duplicates$sample.submitter_id
  )


all_primary_biological_sample_ids <-
  unique(
    luad_cd274_primary$sample.submitter_id
  )


retained_primary_biological_sample_ids <-
  setdiff(
    all_primary_biological_sample_ids,
    removed_biological_sample_ids
  )


cat(
  "\n========================================\n",
  "RETAINED LUAD BIOLOGICAL TUMOR SAMPLES\n",
  "========================================\n"
)


cat(
  "\nDistinct Primary Tumor biological samples before removal:\n"
)


print(
  length(
    all_primary_biological_sample_ids
  )
)


cat(
  "\nBiological tumor samples removed by Script 11:\n"
)


print(
  length(
    removed_biological_sample_ids
  )
)


cat(
  "\nBiological tumor samples retained:\n"
)


print(
  length(
    retained_primary_biological_sample_ids
  )
)


stopifnot(
  length(
    retained_primary_biological_sample_ids
  ) ==
    517
)


# ------------------------------------------------------------
# 5.6 Restrict TPM-level data to retained biological samples
# ------------------------------------------------------------


luad_cd274_retained_aliquots <-
  luad_cd274_primary[
    luad_cd274_primary$sample.submitter_id %in%
      retained_primary_biological_sample_ids,
    ,
    drop = FALSE
  ]


retained_aliquot_counts <-
  table(
    luad_cd274_retained_aliquots$sample.submitter_id
  )


retained_samples_with_multiple_aliquots <-
  retained_aliquot_counts[
    retained_aliquot_counts >
      1
  ]


cat(
  "\nRNA aliquot rows belonging to retained biological tumors:\n"
)


print(
  nrow(
    luad_cd274_retained_aliquots
  )
)


cat(
  "\nRetained biological samples with >1 RNA aliquot:\n"
)


print(
  length(
    retained_samples_with_multiple_aliquots
  )
)


cat(
  "\nNumber of excess technical aliquot rows:\n"
)


print(
  nrow(
    luad_cd274_retained_aliquots
  ) -
    length(
      retained_primary_biological_sample_ids
    )
)


# ------------------------------------------------------------
# 5.7 Collapse technical TPM aliquots
#
# TPM values must NOT be summed.
#
# For the same biological tumor sample represented by
# multiple technical RNA aliquots, use arithmetic mean TPM.
# ------------------------------------------------------------


luad_cd274_biological <-
  aggregate(
    CD274_TPM ~
      cases.submitter_id +
      sample.submitter_id,
    data =
      luad_cd274_retained_aliquots,
    FUN =
      mean
  )


colnames(
  luad_cd274_biological
)[
  colnames(
    luad_cd274_biological
  ) ==
    "cases.submitter_id"
] <-
  "patient_id"


luad_cd274_biological$CD274_log2TPM <-
  log2(
    luad_cd274_biological$CD274_TPM +
      1
  )


stopifnot(
  nrow(
    luad_cd274_biological
  ) ==
    517
)


stopifnot(
  length(
    unique(
      luad_cd274_biological$patient_id
    )
  ) ==
    517
)


stopifnot(
  length(
    unique(
      luad_cd274_biological$sample.submitter_id
    )
  ) ==
    517
)


stopifnot(
  !any(
    is.na(
      luad_cd274_biological$CD274_TPM
    )
  )
)


stopifnot(
  all(
    is.finite(
      luad_cd274_biological$CD274_TPM
    )
  )
)


cat(
  "\n========================================\n",
  "PATIENT-LEVEL CD274 DATASET\n",
  "========================================\n"
)


cat(
  "\nRows / unique patients:\n"
)


print(
  c(
    rows =
      nrow(
        luad_cd274_biological
      ),
    patients =
      length(
        unique(
          luad_cd274_biological$patient_id
        )
      )
  )
)


cat(
  "\nCD274 TPM summary:\n"
)


print(
  summary(
    luad_cd274_biological$CD274_TPM
  )
)


cat(
  "\nCD274 log2(TPM + 1) summary:\n"
)


print(
  summary(
    luad_cd274_biological$CD274_log2TPM
  )
)


# ------------------------------------------------------------
# 5.8 Join CD274 to PCD subtype
# ------------------------------------------------------------


luad_cd274_subtype <-
  merge(
    luad_cd274_biological,
    luad_subtypes,
    by =
      "patient_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )


stopifnot(
  nrow(
    luad_cd274_subtype
  ) ==
    517
)


stopifnot(
  !any(
    is.na(
      luad_cd274_subtype$PCD_cluster
    )
  )
)


luad_cd274_subtype$PCD_cluster <-
  factor(
    luad_cd274_subtype$PCD_cluster,
    levels =
      c(
        "PCD_C1",
        "PCD_C2"
      )
  )


# ------------------------------------------------------------
# 5.9 Join CD274 to MPCDS
# ------------------------------------------------------------


luad_cd274_mpcds <-
  merge(
    luad_cd274_biological,
    luad_mpcds_section2[
      ,
      c(
        "patient_id",
        "MPCDS",
        "MPCDS_z"
      )
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
    luad_cd274_mpcds
  ) ==
    504
)


stopifnot(
  !any(
    is.na(
      luad_cd274_mpcds$MPCDS_z
    )
  )
)


# ============================================================
# 6. CD274 association analyses
# ============================================================


# ------------------------------------------------------------
# 6.1 Continuous MPCDS vs CD274
# ------------------------------------------------------------


cd274_mpcds_cor <-
  cor.test(
    luad_cd274_mpcds$MPCDS_z,
    luad_cd274_mpcds$CD274_log2TPM,
    method =
      "spearman",
    exact =
      FALSE
  )


cat(
  "\n========================================\n",
  "MPCDS VS CD274 EXPRESSION\n",
  "========================================\n"
)


print(
  cd274_mpcds_cor
)


cd274_mpcds_summary <-
  data.frame(
    analysis =
      "MPCDS_z_vs_CD274_log2TPM",
    n =
      nrow(
        luad_cd274_mpcds
      ),
    method =
      "Spearman correlation",
    rho =
      unname(
        cd274_mpcds_cor$estimate
      ),
    statistic =
      unname(
        cd274_mpcds_cor$statistic
      ),
    p_value =
      cd274_mpcds_cor$p.value,
    stringsAsFactors =
      FALSE
  )


# ------------------------------------------------------------
# 6.2 PCD subtype vs CD274
# ------------------------------------------------------------


cd274_subtype_summary <-
  do.call(
    rbind,
    lapply(
      levels(
        luad_cd274_subtype$PCD_cluster
      ),
      function(cluster_name) {
        
        x <-
          luad_cd274_subtype$CD274_log2TPM[
            luad_cd274_subtype$PCD_cluster ==
              cluster_name
          ]
        
        
        data.frame(
          PCD_cluster =
            cluster_name,
          n =
            length(
              x
            ),
          mean =
            mean(
              x
            ),
          sd =
            sd(
              x
            ),
          median =
            median(
              x
            ),
          Q1 =
            unname(
              quantile(
                x,
                probs =
                  0.25
              )
            ),
          Q3 =
            unname(
              quantile(
                x,
                probs =
                  0.75
              )
            ),
          min =
            min(
              x
            ),
          max =
            max(
              x
            ),
          stringsAsFactors =
            FALSE
        )
      }
    )
  )


cd274_subtype_wilcox <-
  wilcox.test(
    CD274_log2TPM ~ PCD_cluster,
    data =
      luad_cd274_subtype,
    exact =
      FALSE,
    conf.int =
      TRUE
  )


cat(
  "\n========================================\n",
  "PCD SUBTYPE VS CD274 EXPRESSION\n",
  "========================================\n"
)


print(
  cd274_subtype_summary,
  row.names =
    FALSE
)


print(
  cd274_subtype_wilcox
)


cd274_subtype_test_summary <-
  data.frame(
    analysis =
      "PCD_C1_vs_PCD_C2_CD274_log2TPM",
    n_PCD_C1 =
      sum(
        luad_cd274_subtype$PCD_cluster ==
          "PCD_C1"
      ),
    n_PCD_C2 =
      sum(
        luad_cd274_subtype$PCD_cluster ==
          "PCD_C2"
      ),
    median_PCD_C1 =
      median(
        luad_cd274_subtype$CD274_log2TPM[
          luad_cd274_subtype$PCD_cluster ==
            "PCD_C1"
        ]
      ),
    median_PCD_C2 =
      median(
        luad_cd274_subtype$CD274_log2TPM[
          luad_cd274_subtype$PCD_cluster ==
            "PCD_C2"
        ]
      ),
    W =
      unname(
        cd274_subtype_wilcox$statistic
      ),
    p_value =
      cd274_subtype_wilcox$p.value,
    stringsAsFactors =
      FALSE
  )


# ============================================================
# 7. Prepare nonsynonymous mutation burden
# ============================================================


luad_mutation_burden <-
  read.csv(
    luad_mutation_burden_file,
    check.names =
      FALSE
  )


stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster",
      "n_nonsyn_mutations"
    ) %in%
      colnames(
        luad_mutation_burden
      )
  )
)


stopifnot(
  !any(
    duplicated(
      luad_mutation_burden$patient_id
    )
  )
)


stopifnot(
  !any(
    is.na(
      luad_mutation_burden$n_nonsyn_mutations
    )
  )
)


stopifnot(
  all(
    is.finite(
      luad_mutation_burden$n_nonsyn_mutations
    )
  )
)


luad_mutation_burden$log2_nonsyn_mutations <-
  log2(
    luad_mutation_burden$n_nonsyn_mutations +
      1
  )


cat(
  "\n========================================\n",
  "NONSYNONYMOUS MUTATION BURDEN COHORT\n",
  "========================================\n"
)


cat(
  "\nRows / unique patients:\n"
)


print(
  c(
    rows =
      nrow(
        luad_mutation_burden
      ),
    patients =
      length(
        unique(
          luad_mutation_burden$patient_id
        )
      )
  )
)


cat(
  "\nMutation burden summary:\n"
)


print(
  summary(
    luad_mutation_burden$n_nonsyn_mutations
  )
)


# ------------------------------------------------------------
# 7.1 Mutation burden + MPCDS
# ------------------------------------------------------------


luad_mutation_mpcds <-
  merge(
    luad_mutation_burden[
      ,
      c(
        "patient_id",
        "n_nonsyn_mutations",
        "log2_nonsyn_mutations"
      )
    ],
    luad_mpcds_section2[
      ,
      c(
        "patient_id",
        "MPCDS",
        "MPCDS_z"
      )
    ],
    by =
      "patient_id",
    all =
      FALSE,
    sort =
      FALSE
  )


stopifnot(
  !any(
    is.na(
      luad_mutation_mpcds$MPCDS_z
    )
  )
)


# ------------------------------------------------------------
# 7.2 Mutation burden + canonical PCD subtype
# ------------------------------------------------------------


luad_mutation_subtype <-
  merge(
    luad_mutation_burden[
      ,
      c(
        "patient_id",
        "n_nonsyn_mutations",
        "log2_nonsyn_mutations"
      )
    ],
    luad_subtypes,
    by =
      "patient_id",
    all =
      FALSE,
    sort =
      FALSE
  )


luad_mutation_subtype$PCD_cluster <-
  factor(
    luad_mutation_subtype$PCD_cluster,
    levels =
      c(
        "PCD_C1",
        "PCD_C2"
      )
  )


# ------------------------------------------------------------
# 7.3 Confirm mutation-file cluster labels agree
# ------------------------------------------------------------


mutation_cluster_check <-
  merge(
    luad_mutation_burden[
      ,
      c(
        "patient_id",
        "PCD_cluster"
      )
    ],
    luad_subtypes,
    by =
      "patient_id",
    suffixes =
      c(
        "_mutation_file",
        "_canonical"
      ),
    all =
      FALSE
  )


stopifnot(
  all(
    mutation_cluster_check$PCD_cluster_mutation_file ==
      mutation_cluster_check$PCD_cluster_canonical
  )
)


cat(
  "\n========================================\n",
  "MUTATION BURDEN COHORT MATCHING\n",
  "========================================\n"
)


cat(
  "\nMutation burden + MPCDS:\n"
)


print(
  nrow(
    luad_mutation_mpcds
  )
)


cat(
  "\nMutation burden + subtype:\n"
)


print(
  nrow(
    luad_mutation_subtype
  )
)


cat(
  "\nSubtype distribution:\n"
)


print(
  table(
    luad_mutation_subtype$PCD_cluster
  )
)


# ============================================================
# 8. Mutation-burden association analyses
# ============================================================


# ------------------------------------------------------------
# 8.1 Continuous MPCDS vs mutation burden
# ------------------------------------------------------------


mutation_mpcds_cor <-
  cor.test(
    luad_mutation_mpcds$MPCDS_z,
    luad_mutation_mpcds$log2_nonsyn_mutations,
    method =
      "spearman",
    exact =
      FALSE
  )


cat(
  "\n========================================\n",
  "MPCDS VS NONSYNONYMOUS MUTATION BURDEN\n",
  "========================================\n"
)


print(
  mutation_mpcds_cor
)


mutation_mpcds_summary <-
  data.frame(
    analysis =
      "MPCDS_z_vs_log2_nonsynonymous_mutation_burden",
    n =
      nrow(
        luad_mutation_mpcds
      ),
    method =
      "Spearman correlation",
    rho =
      unname(
        mutation_mpcds_cor$estimate
      ),
    statistic =
      unname(
        mutation_mpcds_cor$statistic
      ),
    p_value =
      mutation_mpcds_cor$p.value,
    stringsAsFactors =
      FALSE
  )


# ------------------------------------------------------------
# 8.2 PCD subtype vs mutation burden
# ------------------------------------------------------------


mutation_subtype_summary <-
  do.call(
    rbind,
    lapply(
      levels(
        luad_mutation_subtype$PCD_cluster
      ),
      function(cluster_name) {
        
        x <-
          luad_mutation_subtype$log2_nonsyn_mutations[
            luad_mutation_subtype$PCD_cluster ==
              cluster_name
          ]
        
        
        data.frame(
          PCD_cluster =
            cluster_name,
          n =
            length(
              x
            ),
          mean =
            mean(
              x
            ),
          sd =
            sd(
              x
            ),
          median =
            median(
              x
            ),
          Q1 =
            unname(
              quantile(
                x,
                probs =
                  0.25
              )
            ),
          Q3 =
            unname(
              quantile(
                x,
                probs =
                  0.75
              )
            ),
          min =
            min(
              x
            ),
          max =
            max(
              x
            ),
          stringsAsFactors =
            FALSE
        )
      }
    )
  )


mutation_subtype_wilcox <-
  wilcox.test(
    log2_nonsyn_mutations ~ PCD_cluster,
    data =
      luad_mutation_subtype,
    exact =
      FALSE,
    conf.int =
      TRUE
  )


cat(
  "\n========================================\n",
  "PCD SUBTYPE VS NONSYNONYMOUS MUTATION BURDEN\n",
  "========================================\n"
)


print(
  mutation_subtype_summary,
  row.names =
    FALSE
)


print(
  mutation_subtype_wilcox
)


mutation_subtype_test_summary <-
  data.frame(
    analysis =
      "PCD_C1_vs_PCD_C2_log2_nonsynonymous_mutation_burden",
    n_PCD_C1 =
      sum(
        luad_mutation_subtype$PCD_cluster ==
          "PCD_C1"
      ),
    n_PCD_C2 =
      sum(
        luad_mutation_subtype$PCD_cluster ==
          "PCD_C2"
      ),
    median_PCD_C1 =
      median(
        luad_mutation_subtype$log2_nonsyn_mutations[
          luad_mutation_subtype$PCD_cluster ==
            "PCD_C1"
        ]
      ),
    median_PCD_C2 =
      median(
        luad_mutation_subtype$log2_nonsyn_mutations[
          luad_mutation_subtype$PCD_cluster ==
            "PCD_C2"
        ]
      ),
    W =
      unname(
        mutation_subtype_wilcox$statistic
      ),
    p_value =
      mutation_subtype_wilcox$p.value,
    stringsAsFactors =
      FALSE
  )


# ============================================================
# 9. Prepare published MANTIS MSI scores
# ============================================================


# ------------------------------------------------------------
# Source:
# Bonneville et al.
# Landscape of Microsatellite Instability Across 39 Cancer
# Types.
#
# Supplemental File S1 provides sample-level Case ID,
# Cancer Type, and MANTIS Score.
# ------------------------------------------------------------


msi_s1 <-
  readxl::read_excel(
    msi_source_file,
    sheet =
      "Supplemental File S1"
  )


stopifnot(
  all(
    c(
      "Case ID",
      "Cancer Type",
      "MANTIS Score"
    ) %in%
      colnames(
        msi_s1
      )
  )
)


luad_msi_s1 <-
  msi_s1[
    msi_s1$`Cancer Type` ==
      "TCGA-LUAD",
    ,
    drop =
      FALSE
  ]


stopifnot(
  nrow(
    luad_msi_s1
  ) ==
    569
)


stopifnot(
  !any(
    duplicated(
      luad_msi_s1$`Case ID`
    )
  )
)


stopifnot(
  !any(
    is.na(
      luad_msi_s1$`MANTIS Score`
    )
  )
)


luad_msi_core <-
  data.frame(
    patient_id =
      as.character(
        luad_msi_s1$`Case ID`
      ),
    MANTIS_score =
      as.numeric(
        luad_msi_s1$`MANTIS Score`
      ),
    stringsAsFactors =
      FALSE
  )


# ------------------------------------------------------------
# Published MANTIS threshold:
#
# >= 0.4 = MSI-H
# <  0.4 = MSS
#
# Classification is descriptive only in LUAD because the
# MSI-H group is extremely small.
# ------------------------------------------------------------


luad_msi_core$MSI_class <-
  ifelse(
    luad_msi_core$MANTIS_score >=
      0.4,
    "MSI-H",
    "MSS"
  )


cat(
  "\n========================================\n",
  "TCGA-LUAD MANTIS DATA\n",
  "========================================\n"
)


cat(
  "\nRows:\n"
)


print(
  nrow(
    luad_msi_core
  )
)


cat(
  "\nMANTIS score summary:\n"
)


print(
  summary(
    luad_msi_core$MANTIS_score
  )
)


cat(
  "\nMSI classification:\n"
)


print(
  table(
    luad_msi_core$MSI_class
  )
)


# ------------------------------------------------------------
# 9.1 MANTIS + subtype
# ------------------------------------------------------------


luad_msi_subtype <-
  merge(
    luad_msi_core,
    luad_subtypes,
    by =
      "patient_id",
    all =
      FALSE,
    sort =
      FALSE
  )


luad_msi_subtype$PCD_cluster <-
  factor(
    luad_msi_subtype$PCD_cluster,
    levels =
      c(
        "PCD_C1",
        "PCD_C2"
      )
  )


# ------------------------------------------------------------
# 9.2 MANTIS + MPCDS
# ------------------------------------------------------------


luad_msi_mpcds <-
  merge(
    luad_msi_core,
    luad_mpcds_section2[
      ,
      c(
        "patient_id",
        "MPCDS",
        "MPCDS_z"
      )
    ],
    by =
      "patient_id",
    all =
      FALSE,
    sort =
      FALSE
  )


stopifnot(
  !any(
    is.na(
      luad_msi_mpcds$MANTIS_score
    )
  )
)


stopifnot(
  !any(
    is.na(
      luad_msi_mpcds$MPCDS_z
    )
  )
)


cat(
  "\n========================================\n",
  "MANTIS COHORT MATCHING\n",
  "========================================\n"
)


cat(
  "\nMANTIS + subtype:\n"
)


print(
  nrow(
    luad_msi_subtype
  )
)


cat(
  "\nMANTIS + MPCDS:\n"
)


print(
  nrow(
    luad_msi_mpcds
  )
)


cat(
  "\nMSI class by PCD subtype:\n"
)


print(
  table(
    luad_msi_subtype$PCD_cluster,
    luad_msi_subtype$MSI_class
  )
)


# ============================================================
# 10. MANTIS association analyses
# ============================================================


# ------------------------------------------------------------
# 10.1 Continuous MPCDS vs continuous MANTIS score
# ------------------------------------------------------------


mantis_mpcds_cor <-
  cor.test(
    luad_msi_mpcds$MPCDS_z,
    luad_msi_mpcds$MANTIS_score,
    method =
      "spearman",
    exact =
      FALSE
  )


cat(
  "\n========================================\n",
  "MPCDS VS MANTIS SCORE\n",
  "========================================\n"
)


print(
  mantis_mpcds_cor
)


mantis_mpcds_summary <-
  data.frame(
    analysis =
      "MPCDS_z_vs_MANTIS_score",
    n =
      nrow(
        luad_msi_mpcds
      ),
    method =
      "Spearman correlation",
    rho =
      unname(
        mantis_mpcds_cor$estimate
      ),
    statistic =
      unname(
        mantis_mpcds_cor$statistic
      ),
    p_value =
      mantis_mpcds_cor$p.value,
    stringsAsFactors =
      FALSE
  )


# ------------------------------------------------------------
# 10.2 PCD subtype vs continuous MANTIS score
# ------------------------------------------------------------


mantis_subtype_summary <-
  do.call(
    rbind,
    lapply(
      levels(
        luad_msi_subtype$PCD_cluster
      ),
      function(cluster_name) {
        
        x <-
          luad_msi_subtype$MANTIS_score[
            luad_msi_subtype$PCD_cluster ==
              cluster_name
          ]
        
        
        data.frame(
          PCD_cluster =
            cluster_name,
          n =
            length(
              x
            ),
          mean =
            mean(
              x
            ),
          sd =
            sd(
              x
            ),
          median =
            median(
              x
            ),
          Q1 =
            unname(
              quantile(
                x,
                probs =
                  0.25
              )
            ),
          Q3 =
            unname(
              quantile(
                x,
                probs =
                  0.75
              )
            ),
          min =
            min(
              x
            ),
          max =
            max(
              x
            ),
          stringsAsFactors =
            FALSE
        )
      }
    )
  )


mantis_subtype_wilcox <-
  wilcox.test(
    MANTIS_score ~ PCD_cluster,
    data =
      luad_msi_subtype,
    exact =
      FALSE,
    conf.int =
      TRUE
  )


cat(
  "\n========================================\n",
  "PCD SUBTYPE VS MANTIS SCORE\n",
  "========================================\n"
)


print(
  mantis_subtype_summary,
  row.names =
    FALSE
)


print(
  mantis_subtype_wilcox
)


mantis_subtype_test_summary <-
  data.frame(
    analysis =
      "PCD_C1_vs_PCD_C2_MANTIS_score",
    n_PCD_C1 =
      sum(
        luad_msi_subtype$PCD_cluster ==
          "PCD_C1"
      ),
    n_PCD_C2 =
      sum(
        luad_msi_subtype$PCD_cluster ==
          "PCD_C2"
      ),
    median_PCD_C1 =
      median(
        luad_msi_subtype$MANTIS_score[
          luad_msi_subtype$PCD_cluster ==
            "PCD_C1"
        ]
      ),
    median_PCD_C2 =
      median(
        luad_msi_subtype$MANTIS_score[
          luad_msi_subtype$PCD_cluster ==
            "PCD_C2"
        ]
      ),
    W =
      unname(
        mantis_subtype_wilcox$statistic
      ),
    p_value =
      mantis_subtype_wilcox$p.value,
    stringsAsFactors =
      FALSE
  )


cat(
  "\nDescriptive MSI-H counts only:\n"
)


print(
  table(
    luad_msi_subtype$PCD_cluster,
    luad_msi_subtype$MSI_class
  )
)


cat(
  paste0(
    "\nCategorical MSI-H versus MSS testing is not used ",
    "as a primary inferential analysis because only 3 ",
    "MSI-H cases are present in the subtype-matched ",
    "LUAD cohort.\n"
  )
)


# ============================================================
# 11. Consolidate results
# ============================================================


# ------------------------------------------------------------
# 11.1 Continuous MPCDS associations
# ------------------------------------------------------------


continuous_mpcds_associations <-
  data.frame(
    indicator =
      c(
        "CD274_log2TPM",
        "Nonsynonymous_mutation_burden",
        "MANTIS_score"
      ),
    n =
      c(
        nrow(
          luad_cd274_mpcds
        ),
        nrow(
          luad_mutation_mpcds
        ),
        nrow(
          luad_msi_mpcds
        )
      ),
    method =
      rep(
        "Spearman correlation",
        3
      ),
    rho =
      c(
        unname(
          cd274_mpcds_cor$estimate
        ),
        unname(
          mutation_mpcds_cor$estimate
        ),
        unname(
          mantis_mpcds_cor$estimate
        )
      ),
    p_value =
      c(
        cd274_mpcds_cor$p.value,
        mutation_mpcds_cor$p.value,
        mantis_mpcds_cor$p.value
      ),
    stringsAsFactors =
      FALSE
  )


cat(
  "\n========================================\n",
  "CONTINUOUS MPCDS ASSOCIATIONS\n",
  "========================================\n"
)


print(
  continuous_mpcds_associations,
  row.names =
    FALSE
)


# ------------------------------------------------------------
# 11.2 PCD subtype comparisons
# ------------------------------------------------------------


subtype_indicator_comparisons <-
  data.frame(
    indicator =
      c(
        "CD274_log2TPM",
        "Nonsynonymous_mutation_burden",
        "MANTIS_score"
      ),
    n_PCD_C1 =
      c(
        sum(
          luad_cd274_subtype$PCD_cluster ==
            "PCD_C1"
        ),
        sum(
          luad_mutation_subtype$PCD_cluster ==
            "PCD_C1"
        ),
        sum(
          luad_msi_subtype$PCD_cluster ==
            "PCD_C1"
        )
      ),
    n_PCD_C2 =
      c(
        sum(
          luad_cd274_subtype$PCD_cluster ==
            "PCD_C2"
        ),
        sum(
          luad_mutation_subtype$PCD_cluster ==
            "PCD_C2"
        ),
        sum(
          luad_msi_subtype$PCD_cluster ==
            "PCD_C2"
        )
      ),
    median_PCD_C1 =
      c(
        median(
          luad_cd274_subtype$CD274_log2TPM[
            luad_cd274_subtype$PCD_cluster ==
              "PCD_C1"
          ]
        ),
        median(
          luad_mutation_subtype$log2_nonsyn_mutations[
            luad_mutation_subtype$PCD_cluster ==
              "PCD_C1"
          ]
        ),
        median(
          luad_msi_subtype$MANTIS_score[
            luad_msi_subtype$PCD_cluster ==
              "PCD_C1"
          ]
        )
      ),
    median_PCD_C2 =
      c(
        median(
          luad_cd274_subtype$CD274_log2TPM[
            luad_cd274_subtype$PCD_cluster ==
              "PCD_C2"
          ]
        ),
        median(
          luad_mutation_subtype$log2_nonsyn_mutations[
            luad_mutation_subtype$PCD_cluster ==
              "PCD_C2"
          ]
        ),
        median(
          luad_msi_subtype$MANTIS_score[
            luad_msi_subtype$PCD_cluster ==
              "PCD_C2"
          ]
        )
      ),
    W =
      c(
        unname(
          cd274_subtype_wilcox$statistic
        ),
        unname(
          mutation_subtype_wilcox$statistic
        ),
        unname(
          mantis_subtype_wilcox$statistic
        )
      ),
    p_value =
      c(
        cd274_subtype_wilcox$p.value,
        mutation_subtype_wilcox$p.value,
        mantis_subtype_wilcox$p.value
      ),
    stringsAsFactors =
      FALSE
  )


cat(
  "\n========================================\n",
  "PCD SUBTYPE INDICATOR COMPARISONS\n",
  "========================================\n"
)


print(
  subtype_indicator_comparisons,
  row.names =
    FALSE
)


# ============================================================
# 12. Save result tables
# ============================================================


write.csv(
  continuous_mpcds_associations,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_MPCDS_immunotherapy_indicator_correlations.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  subtype_indicator_comparisons,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_PCD_subtype_immunotherapy_indicator_comparisons.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  cd274_subtype_summary,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_CD274_subtype_descriptive_statistics.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  mutation_subtype_summary,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_mutation_burden_subtype_descriptive_statistics.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  mantis_subtype_summary,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_MANTIS_subtype_descriptive_statistics.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  luad_cd274_subtype,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_patient_CD274_with_PCD_subtype.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  luad_mutation_subtype,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_patient_mutation_burden_with_PCD_subtype.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  luad_msi_subtype,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_patient_MANTIS_with_PCD_subtype.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  luad_cd274_mpcds,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_patient_CD274_with_MPCDS.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  luad_mutation_mpcds,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_patient_mutation_burden_with_MPCDS.csv"
  ),
  row.names =
    FALSE
)


write.csv(
  luad_msi_mpcds,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_patient_MANTIS_with_MPCDS.csv"
  ),
  row.names =
    FALSE
)


msi_class_by_subtype <-
  as.data.frame(
    table(
      luad_msi_subtype$PCD_cluster,
      luad_msi_subtype$MSI_class
    )
  )


colnames(
  msi_class_by_subtype
) <-
  c(
    "PCD_cluster",
    "MSI_class",
    "n"
  )


write.csv(
  msi_class_by_subtype,
  file.path(
    immunotherapy_results_dir,
    "TCGA_LUAD_MSI_class_by_PCD_subtype_descriptive.csv"
  ),
  row.names =
    FALSE
)


saved_immunotherapy_files <-
  list.files(
    immunotherapy_results_dir,
    full.names =
      FALSE
  )


cat(
  "\n========================================\n",
  "SAVED IMMUNOTHERAPY-INDICATOR FILES\n",
  "========================================\n"
)


print(
  saved_immunotherapy_files
)


cat(
  "\nNumber of files saved:\n"
)


print(
  length(
    saved_immunotherapy_files
  )
)


# ============================================================
# 13. Figure generation
# ============================================================

# Final manuscript figures are intentionally deferred to Script 32.
# Script 26 produces validated analysis tables only.


# ============================================================
# 14. Final validation
# ============================================================


stopifnot(
  nrow(
    luad_cd274_biological
  ) ==
    517
)


stopifnot(
  nrow(
    luad_cd274_mpcds
  ) ==
    504
)


stopifnot(
  nrow(
    luad_mutation_burden
  ) ==
    505
)


stopifnot(
  nrow(
    luad_mutation_mpcds
  ) ==
    492
)


stopifnot(
  nrow(
    luad_msi_core
  ) ==
    569
)


stopifnot(
  nrow(
    luad_msi_subtype
  ) ==
    514
)


stopifnot(
  nrow(
    luad_msi_mpcds
  ) ==
    501
)


stopifnot(
  sum(
    luad_msi_core$MSI_class ==
      "MSI-H"
  ) ==
    3
)




# ------------------------------------------------------------
# Regression checks against established immunotherapy-indicator results
# ------------------------------------------------------------

stopifnot(
  abs(
    unname(
      cd274_mpcds_cor$estimate
    ) -
      0.1281195
  ) <
    1e-5,
  abs(
    cd274_mpcds_cor$p.value -
      0.003963948
  ) <
    1e-6,
  abs(
    unname(
      mutation_mpcds_cor$estimate
    ) -
      0.26767
  ) <
    1e-4,
  abs(
    mutation_mpcds_cor$p.value -
      1.615e-09
  ) <
    1e-10,
  abs(
    cd274_subtype_test_summary$median_PCD_C1 -
      3.4932
  ) <
    1e-3,
  abs(
    cd274_subtype_test_summary$median_PCD_C2 -
      2.4696
  ) <
    1e-3,
  abs(
    cd274_subtype_test_summary$p_value -
      4.8267e-20
  ) <
    1e-21,
  nrow(
    luad_msi_core
  ) ==
    569,
  sum(
    luad_msi_core$MSI_class ==
      "MSI-H"
  ) ==
    3,
  sum(
    luad_msi_core$MSI_class ==
      "MSS"
  ) ==
    566
)

cat(
  "\\nEstablished immunotherapy-indicator regression checks: PASSED\\n"
)


cat(
  "\n========================================\n",
  "FINAL KEY RESULTS\n",
  "========================================\n"
)


cat(
  "\nContinuous MPCDS associations:\n"
)


print(
  continuous_mpcds_associations,
  row.names =
    FALSE
)


cat(
  "\nPCD subtype comparisons:\n"
)


print(
  subtype_indicator_comparisons,
  row.names =
    FALSE
)


cat(
  "\nMSI class by PCD subtype:\n"
)


print(
  msi_class_by_subtype,
  row.names =
    FALSE
)


cat(
  "\n========================================\n",
  "SCRIPT 26 IMMUNOTHERAPY-RELATED INDICATORS COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)