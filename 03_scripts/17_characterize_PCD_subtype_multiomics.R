# ============================================================
# Script 17: Characterize PCD subtype multi-omics
# ============================================================
#
# Thesis:
# Exploring a Specialized Programmed Cell Death Pattern to
# Predict Prognosis and Treatment Sensitivity in Lung Cancer
# by Machine Learning and Multi-Omics Analysis
#
# Cohorts:
#   TCGA-LUAD
#   TCGA-LUSC
#
# Multi-omics modalities:
#   1. Somatic mutation
#   2. Copy-number variation (CNV)
#   3. DNA methylation
#
# PCD clusters:
#   PCD_C1
#   PCD_C2
#
# Important:
#   - LUAD and LUSC are analysed separately.
#   - Cluster labels remain biologically neutral.
#   - All comparisons are associative, not causal.
#
# ============================================================


# ============================================================
# 0. Shared project configuration
# ============================================================

source("03_scripts/00_project_config.R")


# ============================================================
# 1. Required packages and canonical directories
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "limma"
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

mutation_dir <- file.path(
  processed_dir,
  "mutation"
)

cnv_data_dir <- file.path(
  processed_dir,
  "cnv"
)

methylation_data_dir <- file.path(
  processed_dir,
  "methylation"
)

expression_dir <- file.path(
  processed_dir,
  "expression"
)

clustering_dir <- file.path(
  results_dir,
  "clustering"
)

mutation_results_dir <- file.path(
  results_dir,
  "mutation"
)

cnv_results_dir <- file.path(
  results_dir,
  "cnv"
)

methylation_results_dir <- file.path(
  results_dir,
  "methylation"
)

ensure_dir(mutation_results_dir)
ensure_dir(cnv_results_dir)
ensure_dir(methylation_results_dir)


# ============================================================
# 2. Required input files
# ============================================================

luad_cluster_file <- file.path(
  clustering_dir,
  "TCGA_LUAD_PCD_cluster_assignments.csv"
)

lusc_cluster_file <- file.path(
  clustering_dir,
  "TCGA_LUSC_PCD_cluster_assignments.csv"
)

luad_counts_file <- file.path(
  expression_dir,
  "TCGA_LUAD_raw_counts.rds"
)

lusc_counts_file <- file.path(
  expression_dir,
  "TCGA_LUSC_raw_counts.rds"
)

luad_metadata_file <- file.path(
  expression_dir,
  "TCGA_LUAD_sample_metadata.rds"
)

lusc_metadata_file <- file.path(
  expression_dir,
  "TCGA_LUSC_sample_metadata.rds"
)

luad_mutation_file <- file.path(
  mutation_dir,
  "TCGA_LUAD_mutation_raw_combined.rds"
)

lusc_mutation_file <- file.path(
  mutation_dir,
  "TCGA_LUSC_mutation_raw_combined.rds"
)

luad_cnv_file <- file.path(
  cnv_data_dir,
  "TCGA_LUAD_CNV_raw_combined.rds"
)

lusc_cnv_file <- file.path(
  cnv_data_dir,
  "TCGA_LUSC_CNV_raw_combined.rds"
)

luad_cnv_meta_file <- file.path(
  cnv_data_dir,
  "TCGA_LUAD_CNV_sample_metadata.csv"
)

lusc_cnv_meta_file <- file.path(
  cnv_data_dir,
  "TCGA_LUSC_CNV_sample_metadata.csv"
)

luad_meth_file <- file.path(
  methylation_data_dir,
  "TCGA_LUAD_HM450_beta_matrix.rds"
)

lusc_meth_file <- file.path(
  methylation_data_dir,
  "TCGA_LUSC_HM450_beta_matrix.rds"
)

luad_meth_meta_file <- file.path(
  methylation_data_dir,
  "TCGA_LUAD_HM450_matrix_column_metadata.csv"
)

lusc_meth_meta_file <- file.path(
  methylation_data_dir,
  "TCGA_LUSC_HM450_matrix_column_metadata.csv"
)

required_files <- c(
  luad_cluster_file,
  lusc_cluster_file,
  luad_counts_file,
  lusc_counts_file,
  luad_metadata_file,
  lusc_metadata_file,
  luad_mutation_file,
  lusc_mutation_file,
  luad_cnv_file,
  lusc_cnv_file,
  luad_cnv_meta_file,
  lusc_cnv_meta_file,
  luad_meth_file,
  lusc_meth_file,
  luad_meth_meta_file,
  lusc_meth_meta_file
)

check_files_exist(required_files)


# ============================================================
# 2. Load PCD cluster assignments
# ============================================================

luad_cluster_assignment <- read.csv(
  luad_cluster_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lusc_cluster_assignment <- read.csv(
  lusc_cluster_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster"
    ) %in%
      colnames(luad_cluster_assignment)
  )
)

stopifnot(
  all(
    c(
      "patient_id",
      "PCD_cluster"
    ) %in%
      colnames(lusc_cluster_assignment)
  )
)


# ============================================================
# 3. Load expression data
#
# Required for cross-modality biological-sample matching.
# ============================================================

luad_counts <- readRDS(
  luad_counts_file
)

luad_metadata <- readRDS(
  luad_metadata_file
)

lusc_counts <- readRDS(
  lusc_counts_file
)

lusc_metadata <- readRDS(
  lusc_metadata_file
)


stopifnot(
  identical(
    colnames(luad_counts),
    luad_metadata$cases
  )
)

stopifnot(
  identical(
    colnames(lusc_counts),
    lusc_metadata$cases
  )
)


# ============================================================
# 4. Recover expression-selected Primary Tumor sample
#
# This follows the same rule previously used for patient-level
# expression:
#
#   - Primary Tumor only
#   - one profile per patient
#   - highest raw-count library size
# ============================================================

get_expression_selected_samples <- function(
    counts,
    metadata,
    cluster_assignment
) {
  
  metadata |>
    dplyr::mutate(
      
      patient_id =
        substr(
          cases,
          1,
          12
        ),
      
      library_size =
        colSums(counts)
      
    ) |>
    dplyr::filter(
      
      sample_type ==
        "Primary Tumor",
      
      patient_id %in%
        cluster_assignment$patient_id
      
    ) |>
    dplyr::arrange(
      
      patient_id,
      
      dplyr::desc(
        library_size
      )
      
    ) |>
    dplyr::group_by(
      patient_id
    ) |>
    dplyr::slice_head(
      n = 1
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      
      patient_id,
      
      expression_case =
        cases,
      
      expression_sample_id =
        substr(
          cases,
          1,
          16
        )
      
    )
}


luad_expression_selected <-
  get_expression_selected_samples(
    
    counts =
      luad_counts,
    
    metadata =
      luad_metadata,
    
    cluster_assignment =
      luad_cluster_assignment
    
  )


lusc_expression_selected <-
  get_expression_selected_samples(
    
    counts =
      lusc_counts,
    
    metadata =
      lusc_metadata,
    
    cluster_assignment =
      lusc_cluster_assignment
    
  )


# ============================================================
# PART A
# SOMATIC MUTATION CHARACTERIZATION
# ============================================================


# ============================================================
# 5. Load mutation data
# ============================================================

luad_mutation_raw <- readRDS(
  luad_mutation_file
)

lusc_mutation_raw <- readRDS(
  lusc_mutation_file
)


# ============================================================
# 6. Define nonsynonymous mutation classes
# ============================================================

nonsynonymous_classes <- c(
  
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Translation_Start_Site",
  "Nonstop_Mutation"
  
)


# ============================================================
# 7. Prepare patient-level somatic variants
#
# Rules:
#
#   - Primary Tumor sample code = 01
#   - clustered patients only
#   - somatic variants only
#   - predefined nonsynonymous classes
#   - valid gene symbol required
#   - identical genomic variants occurring in multiple
#     primary aliquots are counted once per patient
#
# GDC_FILTER is retained as annotation and is NOT used as an
# additional blanket exclusion.
# ============================================================

prepare_mutations <- function(
    mutation_data,
    cluster_assignment
) {
  
  mutation_data |>
    dplyr::mutate(
      
      sample_type_code =
        substr(
          Tumor_Sample_Barcode,
          14,
          15
        )
      
    ) |>
    dplyr::filter(
      
      sample_type_code == "01",
      
      patient_id %in%
        cluster_assignment$patient_id,
      
      Mutation_Status ==
        "Somatic",
      
      Variant_Classification %in%
        nonsynonymous_classes,
      
      !is.na(
        Hugo_Symbol
      ),
      
      Hugo_Symbol != ""
      
    ) |>
    dplyr::distinct(
      
      patient_id,
      Hugo_Symbol,
      Chromosome,
      Start_Position,
      End_Position,
      Reference_Allele,
      Tumor_Seq_Allele2,
      Variant_Classification,
      
      .keep_all = TRUE
      
    )
}


luad_mut_patient <-
  prepare_mutations(
    
    mutation_data =
      luad_mutation_raw,
    
    cluster_assignment =
      luad_cluster_assignment
    
  )


lusc_mut_patient <-
  prepare_mutations(
    
    mutation_data =
      lusc_mutation_raw,
    
    cluster_assignment =
      lusc_cluster_assignment
    
  )


# ============================================================
# 8. Determine patients with Primary Tumor mutation coverage
#
# Absence from the MAF is not automatically interpreted as
# zero mutation unless Primary Tumor mutation data exist.
# ============================================================

get_mutation_coverage <- function(
    mutation_data,
    cluster_assignment
) {
  
  mutation_data |>
    dplyr::mutate(
      
      sample_type_code =
        substr(
          Tumor_Sample_Barcode,
          14,
          15
        )
      
    ) |>
    dplyr::filter(
      
      sample_type_code == "01",
      
      patient_id %in%
        cluster_assignment$patient_id
      
    ) |>
    dplyr::distinct(
      patient_id
    ) |>
    dplyr::left_join(
      
      cluster_assignment,
      
      by =
        "patient_id"
      
    )
}


luad_mutation_coverage <-
  get_mutation_coverage(
    
    mutation_data =
      luad_mutation_raw,
    
    cluster_assignment =
      luad_cluster_assignment
    
  )


lusc_mutation_coverage <-
  get_mutation_coverage(
    
    mutation_data =
      lusc_mutation_raw,
    
    cluster_assignment =
      lusc_cluster_assignment
    
  )


# ============================================================
# 9. Patient-gene mutation presence
# ============================================================

luad_patient_gene_presence <-
  luad_mut_patient |>
  dplyr::distinct(
    patient_id,
    Hugo_Symbol
  )


lusc_patient_gene_presence <-
  lusc_mut_patient |>
  dplyr::distinct(
    patient_id,
    Hugo_Symbol
  )


# ============================================================
# 10. Overall gene mutation frequencies
# ============================================================

calculate_gene_frequency <- function(
    patient_gene_presence,
    mutation_coverage
) {
  
  denominator <-
    dplyr::n_distinct(
      mutation_coverage$patient_id
    )
  
  patient_gene_presence |>
    dplyr::count(
      
      Hugo_Symbol,
      
      name =
        "mutated_patients"
      
    ) |>
    dplyr::mutate(
      
      covered_patients =
        denominator,
      
      mutation_frequency =
        mutated_patients /
        covered_patients
      
    ) |>
    dplyr::arrange(
      
      dplyr::desc(
        mutation_frequency
      )
      
    )
}


luad_mutation_frequency <-
  calculate_gene_frequency(
    
    patient_gene_presence =
      luad_patient_gene_presence,
    
    mutation_coverage =
      luad_mutation_coverage
    
  )


lusc_mutation_frequency <-
  calculate_gene_frequency(
    
    patient_gene_presence =
      lusc_patient_gene_presence,
    
    mutation_coverage =
      lusc_mutation_coverage
    
  )


# ============================================================
# 11. Gene-wise Fisher tests
#
# Odds ratio:
#
#   odds of mutation in PCD_C1 /
#   odds of mutation in PCD_C2
#
# OR > 1 -> enriched in C1
# OR < 1 -> enriched in C2
# ============================================================

run_gene_fisher <- function(
    patient_gene_presence,
    mutation_coverage
) {
  
  genes <-
    sort(
      unique(
        patient_gene_presence$Hugo_Symbol
      )
    )
  
  c1_patients <-
    mutation_coverage$patient_id[
      mutation_coverage$PCD_cluster ==
        "PCD_C1"
    ]
  
  c2_patients <-
    mutation_coverage$patient_id[
      mutation_coverage$PCD_cluster ==
        "PCD_C2"
    ]
  
  n_c1 <-
    length(
      unique(
        c1_patients
      )
    )
  
  n_c2 <-
    length(
      unique(
        c2_patients
      )
    )
  
  
  results <- lapply(
    genes,
    function(gene) {
      
      mutated <-
        patient_gene_presence$patient_id[
          patient_gene_presence$Hugo_Symbol ==
            gene
        ]
      
      c1_mut <-
        sum(
          c1_patients %in%
            mutated
        )
      
      c2_mut <-
        sum(
          c2_patients %in%
            mutated
        )
      
      
      contingency <- matrix(
        
        c(
          c1_mut,
          n_c1 - c1_mut,
          c2_mut,
          n_c2 - c2_mut
        ),
        
        nrow = 2,
        
        byrow = TRUE
        
      )
      
      
      fisher_result <-
        fisher.test(
          contingency
        )
      
      
      data.frame(
        
        gene =
          gene,
        
        C1_mutated =
          c1_mut,
        
        C1_total =
          n_c1,
        
        C1_frequency =
          c1_mut /
          n_c1,
        
        C2_mutated =
          c2_mut,
        
        C2_total =
          n_c2,
        
        C2_frequency =
          c2_mut /
          n_c2,
        
        odds_ratio_C1_vs_C2 =
          unname(
            fisher_result$estimate
          ),
        
        pvalue =
          fisher_result$p.value,
        
        stringsAsFactors =
          FALSE
        
      )
      
    }
  )
  
  
  dplyr::bind_rows(
    results
  ) |>
    dplyr::mutate(
      
      padj =
        p.adjust(
          pvalue,
          method = "BH"
        )
      
    ) |>
    dplyr::arrange(
      padj,
      pvalue
    )
}


luad_mutation_fisher <-
  run_gene_fisher(
    
    patient_gene_presence =
      luad_patient_gene_presence,
    
    mutation_coverage =
      luad_mutation_coverage
    
  )


lusc_mutation_fisher <-
  run_gene_fisher(
    
    patient_gene_presence =
      lusc_patient_gene_presence,
    
    mutation_coverage =
      lusc_mutation_coverage
    
  )


# ============================================================
# 12. Cross-histology significant mutation associations
# ============================================================

shared_significant_mutation_genes <-
  dplyr::inner_join(
    
    luad_mutation_fisher |>
      dplyr::filter(
        padj < 0.05
      ) |>
      dplyr::select(
        
        gene,
        
        LUAD_C1_frequency =
          C1_frequency,
        
        LUAD_C2_frequency =
          C2_frequency,
        
        LUAD_OR =
          odds_ratio_C1_vs_C2,
        
        LUAD_padj =
          padj
        
      ),
    
    lusc_mutation_fisher |>
      dplyr::filter(
        padj < 0.05
      ) |>
      dplyr::select(
        
        gene,
        
        LUSC_C1_frequency =
          C1_frequency,
        
        LUSC_C2_frequency =
          C2_frequency,
        
        LUSC_OR =
          odds_ratio_C1_vs_C2,
        
        LUSC_padj =
          padj
        
      ),
    
    by =
      "gene"
    
  )


# ============================================================
# 13. Patient-level nonsynonymous mutation burden
#
# This is a nonsynonymous mutation count/burden.
#
# It is NOT called TMB because no callable-Mb denominator is
# available.
# Do not relabel this variable as clinical TMB.
# ============================================================

calculate_mutation_burden <- function(
    mutation_patient,
    mutation_coverage
) {
  
  counts <-
    mutation_patient |>
    dplyr::count(
      
      patient_id,
      
      name =
        "n_nonsyn_mutations"
      
    )
  
  
  mutation_coverage |>
    dplyr::left_join(
      
      counts,
      
      by =
        "patient_id"
      
    ) |>
    dplyr::mutate(
      
      n_nonsyn_mutations =
        tidyr::replace_na(
          n_nonsyn_mutations,
          0L
        )
      
    )
}


luad_mutation_burden <-
  calculate_mutation_burden(
    
    mutation_patient =
      luad_mut_patient,
    
    mutation_coverage =
      luad_mutation_coverage
    
  )


lusc_mutation_burden <-
  calculate_mutation_burden(
    
    mutation_patient =
      lusc_mut_patient,
    
    mutation_coverage =
      lusc_mutation_coverage
    
  )


# ============================================================
# 14. Mutation burden summaries
# ============================================================

summarize_mutation_burden <- function(
    burden_data
) {
  
  burden_data |>
    dplyr::group_by(
      PCD_cluster
    ) |>
    dplyr::summarise(
      
      n =
        dplyr::n(),
      
      median =
        median(
          n_nonsyn_mutations
        ),
      
      mean =
        mean(
          n_nonsyn_mutations
        ),
      
      IQR =
        IQR(
          n_nonsyn_mutations
        ),
      
      .groups =
        "drop"
      
    )
}


luad_mutation_burden_summary <-
  summarize_mutation_burden(
    luad_mutation_burden
  )


lusc_mutation_burden_summary <-
  summarize_mutation_burden(
    lusc_mutation_burden
  )


luad_mutation_burden_test <-
  wilcox.test(
    
    n_nonsyn_mutations ~ PCD_cluster,
    
    data =
      luad_mutation_burden,
    
    exact =
      FALSE
    
  )


lusc_mutation_burden_test <-
  wilcox.test(
    
    n_nonsyn_mutations ~ PCD_cluster,
    
    data =
      lusc_mutation_burden,
    
    exact =
      FALSE
    
  )


mutation_burden_tests <- data.frame(
  
  histology =
    c(
      "LUAD",
      "LUSC"
    ),
  
  wilcoxon_W =
    c(
      unname(
        luad_mutation_burden_test$statistic
      ),
      
      unname(
        lusc_mutation_burden_test$statistic
      )
    ),
  
  pvalue =
    c(
      luad_mutation_burden_test$p.value,
      lusc_mutation_burden_test$p.value
    )
  
)


# ============================================================
# 15. Mutation burden sensitivity analysis
#
# Restrict to patients having exactly one Primary Tumor
# mutation sample.
# ============================================================

count_primary_mutation_samples <- function(
    mutation_data,
    cluster_assignment
) {
  
  mutation_data |>
    dplyr::mutate(
      
      sample_type_code =
        substr(
          Tumor_Sample_Barcode,
          14,
          15
        )
      
    ) |>
    dplyr::filter(
      
      sample_type_code ==
        "01",
      
      patient_id %in%
        cluster_assignment$patient_id
      
    ) |>
    dplyr::distinct(
      
      patient_id,
      Tumor_Sample_Barcode
      
    ) |>
    dplyr::count(
      
      patient_id,
      
      name =
        "n_primary_samples"
      
    )
}


luad_primary_per_patient <-
  count_primary_mutation_samples(
    
    mutation_data =
      luad_mutation_raw,
    
    cluster_assignment =
      luad_cluster_assignment
    
  )


lusc_primary_per_patient <-
  count_primary_mutation_samples(
    
    mutation_data =
      lusc_mutation_raw,
    
    cluster_assignment =
      lusc_cluster_assignment
    
  )


luad_burden_check <-
  luad_mutation_burden |>
  dplyr::left_join(
    
    luad_primary_per_patient,
    
    by =
      "patient_id"
    
  )


lusc_burden_check <-
  lusc_mutation_burden |>
  dplyr::left_join(
    
    lusc_primary_per_patient,
    
    by =
      "patient_id"
    
  )


luad_burden_single <-
  luad_burden_check |>
  dplyr::filter(
    n_primary_samples == 1
  )


lusc_burden_single <-
  lusc_burden_check |>
  dplyr::filter(
    n_primary_samples == 1
  )


luad_burden_single_summary <-
  summarize_mutation_burden(
    luad_burden_single
  )


lusc_burden_single_summary <-
  summarize_mutation_burden(
    lusc_burden_single
  )


luad_burden_single_test <-
  wilcox.test(
    
    n_nonsyn_mutations ~ PCD_cluster,
    
    data =
      luad_burden_single,
    
    exact =
      FALSE
    
  )


lusc_burden_single_test <-
  wilcox.test(
    
    n_nonsyn_mutations ~ PCD_cluster,
    
    data =
      lusc_burden_single,
    
    exact =
      FALSE
    
  )


mutation_burden_sensitivity_tests <-
  data.frame(
    
    histology =
      c(
        "LUAD",
        "LUSC"
      ),
    
    n_single_primary =
      c(
        nrow(
          luad_burden_single
        ),
        
        nrow(
          lusc_burden_single
        )
      ),
    
    wilcoxon_W =
      c(
        unname(
          luad_burden_single_test$statistic
        ),
        
        unname(
          lusc_burden_single_test$statistic
        )
      ),
    
    pvalue =
      c(
        luad_burden_single_test$p.value,
        lusc_burden_single_test$p.value
      )
    
  )


# ============================================================
# 16. Mutation QC summaries
# ============================================================

luad_gdc_filter_qc <-
  as.data.frame(
    table(
      luad_mut_patient$GDC_FILTER,
      useNA = "ifany"
    )
  )


lusc_gdc_filter_qc <-
  as.data.frame(
    table(
      lusc_mut_patient$GDC_FILTER,
      useNA = "ifany"
    )
  )


luad_mutation_status_qc <-
  as.data.frame(
    table(
      luad_mut_patient$Mutation_Status,
      useNA = "ifany"
    )
  )


lusc_mutation_status_qc <-
  as.data.frame(
    table(
      lusc_mut_patient$Mutation_Status,
      useNA = "ifany"
    )
  )


# ============================================================
# 17. Save mutation results
# ============================================================

write.csv(
  
  luad_mutation_frequency,
  
  file.path(mutation_results_dir, "TCGA_LUAD_gene_mutation_frequency.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_mutation_frequency,
  
  file.path(mutation_results_dir, "TCGA_LUSC_gene_mutation_frequency.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  luad_mutation_fisher,
  
  file.path(mutation_results_dir, "TCGA_LUAD_PCD_cluster_mutation_fisher.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_mutation_fisher,
  
  file.path(mutation_results_dir, "TCGA_LUSC_PCD_cluster_mutation_fisher.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  shared_significant_mutation_genes,
  
  file.path(mutation_results_dir, "TCGA_LUAD_LUSC_shared_significant_mutation_genes.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  luad_mutation_burden,
  
  file.path(mutation_results_dir, "TCGA_LUAD_nonsynonymous_mutation_burden.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_mutation_burden,
  
  file.path(mutation_results_dir, "TCGA_LUSC_nonsynonymous_mutation_burden.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  mutation_burden_tests,
  
  file.path(mutation_results_dir, "TCGA_mutation_burden_cluster_tests.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  mutation_burden_sensitivity_tests,
  
  file.path(mutation_results_dir, "TCGA_mutation_burden_single_primary_sensitivity.csv"),
  
  row.names = FALSE
  
)


# ============================================================
# PART B
# COPY-NUMBER VARIATION CHARACTERIZATION
# ============================================================


# ============================================================
# 18. Load CNV data
# ============================================================

luad_cnv <- readRDS(
  luad_cnv_file
)

lusc_cnv <- readRDS(
  lusc_cnv_file
)


luad_cnv_meta <- read.csv(
  luad_cnv_meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lusc_cnv_meta <- read.csv(
  lusc_cnv_meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 19. Restrict CNV metadata to clustered Primary Tumors
# ============================================================

luad_cnv_meta_cluster <-
  luad_cnv_meta |>
  dplyr::filter(
    
    sample_type ==
      "Primary Tumor",
    
    patient_id %in%
      luad_cluster_assignment$patient_id
    
  )


lusc_cnv_meta_cluster <-
  lusc_cnv_meta |>
  dplyr::filter(
    
    sample_type ==
      "Primary Tumor",
    
    patient_id %in%
      lusc_cluster_assignment$patient_id
    
  )


# ============================================================
# 20. CNV profile quality metrics
#
# Segment_Count is descriptive only.
#
# Profile selection does NOT prioritize Segment_Count.
# ============================================================

calculate_cnv_file_quality <- function(
    cnv_data,
    cluster_assignment
) {
  
  cnv_data |>
    dplyr::filter(
      
      sample_type ==
        "Primary Tumor",
      
      patient_id %in%
        cluster_assignment$patient_id
      
    ) |>
    dplyr::group_by(
      
      patient_id,
      sample.submitter_id,
      cases,
      source_file,
      GDC_Aliquot
      
    ) |>
    dplyr::summarise(
      
      Segment_Count =
        dplyr::n(),
      
      Total_Probes =
        sum(
          Num_Probes,
          na.rm = TRUE
        ),
      
      Genome_Length =
        sum(
          End - Start + 1,
          na.rm = TRUE
        ),
      
      .groups =
        "drop"
      
    )
}


luad_cnv_file_quality <-
  calculate_cnv_file_quality(
    
    cnv_data =
      luad_cnv,
    
    cluster_assignment =
      luad_cluster_assignment
    
  )


lusc_cnv_file_quality <-
  calculate_cnv_file_quality(
    
    cnv_data =
      lusc_cnv,
    
    cluster_assignment =
      lusc_cluster_assignment
    
  )


# ============================================================
# 21. Select one CNV profile per patient
#
# Selection hierarchy:
#
#   1. Prefer CNV biological sample matching expression sample.
#   2. Highest Total_Probes.
#   3. source_file alphabetically as deterministic tie-breaker.
# ============================================================

select_one_cnv_profile <- function(
    cnv_file_quality,
    expression_selected
) {
  
  cnv_file_quality |>
    dplyr::left_join(
      
      expression_selected,
      
      by =
        "patient_id"
      
    ) |>
    dplyr::mutate(
      
      expression_match =
        sample.submitter_id ==
        expression_sample_id
      
    ) |>
    dplyr::group_by(
      patient_id
    ) |>
    dplyr::arrange(
      
      dplyr::desc(
        expression_match
      ),
      
      dplyr::desc(
        Total_Probes
      ),
      
      source_file,
      
      .by_group =
        TRUE
      
    ) |>
    dplyr::slice_head(
      n = 1
    ) |>
    dplyr::ungroup()
}


luad_cnv_selected <-
  select_one_cnv_profile(
    
    cnv_file_quality =
      luad_cnv_file_quality,
    
    expression_selected =
      luad_expression_selected
    
  )


lusc_cnv_selected <-
  select_one_cnv_profile(
    
    cnv_file_quality =
      lusc_cnv_file_quality,
    
    expression_selected =
      lusc_expression_selected
    
  )


stopifnot(
  !any(
    duplicated(
      luad_cnv_selected$patient_id
    )
  )
)

stopifnot(
  !any(
    duplicated(
      lusc_cnv_selected$patient_id
    )
  )
)


# ============================================================
# 22. Extract selected CNV segments
# ============================================================

luad_cnv_final <-
  luad_cnv |>
  dplyr::semi_join(
    
    luad_cnv_selected |>
      dplyr::select(
        
        patient_id,
        source_file,
        GDC_Aliquot
        
      ),
    
    by =
      c(
        "patient_id",
        "source_file",
        "GDC_Aliquot"
      )
    
  )


lusc_cnv_final <-
  lusc_cnv |>
  dplyr::semi_join(
    
    lusc_cnv_selected |>
      dplyr::select(
        
        patient_id,
        source_file,
        GDC_Aliquot
        
      ),
    
    by =
      c(
        "patient_id",
        "source_file",
        "GDC_Aliquot"
      )
    
  )


# ============================================================
# 23. Continuous genome-wide CNV burden
#
# Autosomes only: chromosomes 1-22.
#
# No amplification/deletion threshold is imposed.
#
# absolute_CNV_burden:
#   length-weighted mean |Segment_Mean|
#
# gain_burden:
#   length-weighted positive Segment_Mean component
#
# loss_burden:
#   length-weighted magnitude of negative Segment_Mean
# ============================================================

prepare_cnv_burden <- function(
    cnv_segments,
    cluster_assignment
) {
  
  cnv_segments |>
    dplyr::mutate(
      
      Chromosome =
        as.character(
          Chromosome
        ),
      
      segment_length =
        End - Start + 1
      
    ) |>
    dplyr::filter(
      
      Chromosome %in%
        as.character(
          1:22
        ),
      
      !is.na(
        Segment_Mean
      ),
      
      !is.na(
        segment_length
      ),
      
      segment_length > 0
      
    ) |>
    dplyr::mutate(
      
      absolute_component =
        abs(
          Segment_Mean
        ),
      
      gain_component =
        pmax(
          Segment_Mean,
          0
        ),
      
      loss_component =
        pmax(
          -Segment_Mean,
          0
        )
      
    ) |>
    dplyr::group_by(
      patient_id
    ) |>
    dplyr::summarise(
      
      covered_length =
        sum(
          segment_length
        ),
      
      absolute_CNV_burden =
        weighted.mean(
          
          absolute_component,
          
          w =
            segment_length
          
        ),
      
      gain_burden =
        weighted.mean(
          
          gain_component,
          
          w =
            segment_length
          
        ),
      
      loss_burden =
        weighted.mean(
          
          loss_component,
          
          w =
            segment_length
          
        ),
      
      weighted_signed_mean =
        weighted.mean(
          
          Segment_Mean,
          
          w =
            segment_length
          
        ),
      
      .groups =
        "drop"
      
    ) |>
    dplyr::left_join(
      
      cluster_assignment |>
        dplyr::select(
          
          patient_id,
          PCD_cluster
          
        ),
      
      by =
        "patient_id"
      
    )
}


luad_cnv_burden <-
  prepare_cnv_burden(
    
    cnv_segments =
      luad_cnv_final,
    
    cluster_assignment =
      luad_cluster_assignment
    
  )


lusc_cnv_burden <-
  prepare_cnv_burden(
    
    cnv_segments =
      lusc_cnv_final,
    
    cluster_assignment =
      lusc_cluster_assignment
    
  )


# ============================================================
# 24. CNV cluster summaries
# ============================================================

summarize_cnv_burden <- function(
    cnv_burden
) {
  
  cnv_burden |>
    dplyr::group_by(
      PCD_cluster
    ) |>
    dplyr::summarise(
      
      n =
        dplyr::n(),
      
      median_absolute =
        median(
          absolute_CNV_burden
        ),
      
      mean_absolute =
        mean(
          absolute_CNV_burden
        ),
      
      IQR_absolute =
        IQR(
          absolute_CNV_burden
        ),
      
      median_gain =
        median(
          gain_burden
        ),
      
      median_loss =
        median(
          loss_burden
        ),
      
      .groups =
        "drop"
      
    )
}


luad_cnv_burden_summary <-
  summarize_cnv_burden(
    luad_cnv_burden
  )


lusc_cnv_burden_summary <-
  summarize_cnv_burden(
    lusc_cnv_burden
  )


# ============================================================
# 25. CNV Wilcoxon tests
# ============================================================

run_cnv_tests <- function(
    cnv_burden
) {
  
  absolute_test <-
    wilcox.test(
      
      absolute_CNV_burden ~
        PCD_cluster,
      
      data =
        cnv_burden,
      
      exact =
        FALSE
      
    )
  
  
  gain_test <-
    wilcox.test(
      
      gain_burden ~
        PCD_cluster,
      
      data =
        cnv_burden,
      
      exact =
        FALSE
      
    )
  
  
  loss_test <-
    wilcox.test(
      
      loss_burden ~
        PCD_cluster,
      
      data =
        cnv_burden,
      
      exact =
        FALSE
      
    )
  
  
  data.frame(
    
    metric =
      c(
        "absolute_CNV_burden",
        "gain_burden",
        "loss_burden"
      ),
    
    pvalue =
      c(
        absolute_test$p.value,
        gain_test$p.value,
        loss_test$p.value
      )
    
  ) |>
    dplyr::mutate(
      
      padj =
        p.adjust(
          pvalue,
          method = "BH"
        )
      
    )
}


luad_cnv_test_summary <-
  run_cnv_tests(
    luad_cnv_burden
  )


lusc_cnv_test_summary <-
  run_cnv_tests(
    lusc_cnv_burden
  )


# ============================================================
# 26. Save CNV results
# ============================================================

write.csv(
  
  luad_cnv_selected,
  
  file.path(cnv_results_dir, "TCGA_LUAD_selected_CNV_profiles.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_cnv_selected,
  
  file.path(cnv_results_dir, "TCGA_LUSC_selected_CNV_profiles.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  luad_cnv_burden,
  
  file.path(cnv_results_dir, "TCGA_LUAD_continuous_CNV_burden.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_cnv_burden,
  
  file.path(cnv_results_dir, "TCGA_LUSC_continuous_CNV_burden.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  luad_cnv_burden_summary,
  
  file.path(cnv_results_dir, "TCGA_LUAD_CNV_burden_summary.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_cnv_burden_summary,
  
  file.path(cnv_results_dir, "TCGA_LUSC_CNV_burden_summary.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  luad_cnv_test_summary,
  
  file.path(cnv_results_dir, "TCGA_LUAD_CNV_cluster_tests.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_cnv_test_summary,
  
  file.path(cnv_results_dir, "TCGA_LUSC_CNV_cluster_tests.csv"),
  
  row.names = FALSE
  
)


# ============================================================
# PART C
# DNA METHYLATION CHARACTERIZATION
# ============================================================


# ============================================================
# 27. Load HM450 beta matrices
# ============================================================

luad_meth <- readRDS(
  luad_meth_file
)

lusc_meth <- readRDS(
  lusc_meth_file
)


luad_meth_colmeta <- read.csv(
  
  luad_meth_meta_file,
  
  stringsAsFactors =
    FALSE,
  
  check.names =
    FALSE
  
)


lusc_meth_colmeta <- read.csv(
  
  lusc_meth_meta_file,
  
  stringsAsFactors =
    FALSE,
  
  check.names =
    FALSE
  
)


stopifnot(
  identical(
    
    colnames(
      luad_meth
    ),
    
    luad_meth_colmeta$file_name
    
  )
)


stopifnot(
  identical(
    
    colnames(
      lusc_meth
    ),
    
    lusc_meth_colmeta$file_name
    
  )
)


# ============================================================
# 28. Select one expression-matched methylation profile
#
# Exact biological-sample match is required.
#
# When multiple methylation files exist for the same matched
# biological sample, file_name provides a deterministic
# tie-breaker.
# ============================================================

select_one_methylation_sample <- function(
    meth_colmeta,
    expression_selected,
    cluster_assignment
) {
  
  meth_colmeta |>
    dplyr::filter(
      
      patient_id %in%
        cluster_assignment$patient_id
      
    ) |>
    dplyr::left_join(
      
      expression_selected,
      
      by =
        "patient_id"
      
    ) |>
    dplyr::mutate(
      
      expression_match =
        sample.submitter_id ==
        expression_sample_id
      
    ) |>
    dplyr::filter(
      expression_match
    ) |>
    dplyr::arrange(
      
      patient_id,
      file_name
      
    ) |>
    dplyr::group_by(
      patient_id
    ) |>
    dplyr::slice_head(
      n = 1
    ) |>
    dplyr::ungroup()
}


luad_meth_selected <-
  select_one_methylation_sample(
    
    meth_colmeta =
      luad_meth_colmeta,
    
    expression_selected =
      luad_expression_selected,
    
    cluster_assignment =
      luad_cluster_assignment
    
  )


lusc_meth_selected <-
  select_one_methylation_sample(
    
    meth_colmeta =
      lusc_meth_colmeta,
    
    expression_selected =
      lusc_expression_selected,
    
    cluster_assignment =
      lusc_cluster_assignment
    
  )


# ============================================================
# 29. Extract selected methylation columns
# ============================================================

luad_meth_final <-
  luad_meth[
    ,
    match(
      
      luad_meth_selected$file_name,
      
      colnames(
        luad_meth
      )
      
    ),
    drop = FALSE
  ]


lusc_meth_final <-
  lusc_meth[
    ,
    match(
      
      lusc_meth_selected$file_name,
      
      colnames(
        lusc_meth
      )
      
    ),
    drop = FALSE
  ]


stopifnot(
  identical(
    
    colnames(
      luad_meth_final
    ),
    
    luad_meth_selected$file_name
    
  )
)


stopifnot(
  identical(
    
    colnames(
      lusc_meth_final
    ),
    
    lusc_meth_selected$file_name
    
  )
)


# ============================================================
# 30. Attach PCD clusters
# ============================================================

luad_meth_selected <-
  luad_meth_selected |>
  dplyr::left_join(
    
    luad_cluster_assignment |>
      dplyr::select(
        
        patient_id,
        PCD_cluster
        
      ),
    
    by =
      "patient_id"
    
  )


lusc_meth_selected <-
  lusc_meth_selected |>
  dplyr::left_join(
    
    lusc_cluster_assignment |>
      dplyr::select(
        
        patient_id,
        PCD_cluster
        
      ),
    
    by =
      "patient_id"
    
  )


# ============================================================
# 31. Final missingness QC
#
# Reapply <=10% probe missingness after restricting to the
# final cluster-matched cohort.
#
# Remaining values are NOT imputed.
# ============================================================

luad_missing_final <-
  rowMeans(
    is.na(
      luad_meth_final
    )
  )


lusc_missing_final <-
  rowMeans(
    is.na(
      lusc_meth_final
    )
  )


luad_meth_analysis <-
  luad_meth_final[
    luad_missing_final <= 0.10,
    ,
    drop = FALSE
  ]


lusc_meth_analysis <-
  lusc_meth_final[
    lusc_missing_final <= 0.10,
    ,
    drop = FALSE
  ]


# ============================================================
# 32. Convert beta values to M-values
#
# M-values are used for statistical modelling.
#
# Beta values are retained for biological interpretation.
# ============================================================

luad_m_values <-
  log2(
    
    luad_meth_analysis /
      (
        1 -
          luad_meth_analysis
      )
    
  )


lusc_m_values <-
  log2(
    
    lusc_meth_analysis /
      (
        1 -
          lusc_meth_analysis
      )
    
  )


stopifnot(
  sum(
    is.infinite(
      luad_m_values
    )
  ) == 0
)


stopifnot(
  sum(
    is.infinite(
      lusc_m_values
    )
  ) == 0
)


# ============================================================
# 33. Differential methylation design
# ============================================================

luad_meth_group <- factor(
  
  luad_meth_selected$PCD_cluster,
  
  levels =
    c(
      "PCD_C1",
      "PCD_C2"
    )
  
)


lusc_meth_group <- factor(
  
  lusc_meth_selected$PCD_cluster,
  
  levels =
    c(
      "PCD_C1",
      "PCD_C2"
    )
  
)


luad_meth_design <-
  model.matrix(
    ~ 0 + luad_meth_group
  )


lusc_meth_design <-
  model.matrix(
    ~ 0 + lusc_meth_group
  )


colnames(
  luad_meth_design
) <-
  c(
    "PCD_C1",
    "PCD_C2"
  )


colnames(
  lusc_meth_design
) <-
  c(
    "PCD_C1",
    "PCD_C2"
  )


# ============================================================
# 34. limma differential methylation
#
# Contrast:
#
#   PCD_C2 - PCD_C1
#
# Positive effect:
#   higher methylation in C2
#
# Negative effect:
#   lower methylation in C2
# ============================================================

meth_contrast <-
  limma::makeContrasts(
    
    PCD_C2_vs_PCD_C1 =
      PCD_C2 - PCD_C1,
    
    levels =
      luad_meth_design
    
  )


luad_meth_fit <-
  limma::lmFit(
    
    luad_m_values,
    
    luad_meth_design
    
  )


luad_meth_fit <-
  limma::contrasts.fit(
    
    luad_meth_fit,
    
    meth_contrast
    
  )


luad_meth_fit <-
  limma::eBayes(
    luad_meth_fit
  )


lusc_meth_fit <-
  limma::lmFit(
    
    lusc_m_values,
    
    lusc_meth_design
    
  )


lusc_meth_fit <-
  limma::contrasts.fit(
    
    lusc_meth_fit,
    
    meth_contrast
    
  )


lusc_meth_fit <-
  limma::eBayes(
    lusc_meth_fit
  )


# ============================================================
# 35. Extract complete probe-level statistics
# ============================================================

luad_meth_results <-
  limma::topTable(
    
    luad_meth_fit,
    
    coef =
      "PCD_C2_vs_PCD_C1",
    
    number =
      Inf,
    
    adjust.method =
      "BH",
    
    sort.by =
      "none"
    
  )


lusc_meth_results <-
  limma::topTable(
    
    lusc_meth_fit,
    
    coef =
      "PCD_C2_vs_PCD_C1",
    
    number =
      Inf,
    
    adjust.method =
      "BH",
    
    sort.by =
      "none"
    
  )


luad_meth_results$probe_id <-
  rownames(
    luad_meth_results
  )


lusc_meth_results$probe_id <-
  rownames(
    lusc_meth_results
  )


# ============================================================
# 36. Mean beta values and delta-beta
#
# delta_beta = C2 - C1
# ============================================================

luad_c1_beta <-
  rowMeans(
    
    luad_meth_analysis[
      ,
      luad_meth_group ==
        "PCD_C1",
      drop = FALSE
    ],
    
    na.rm =
      TRUE
    
  )


luad_c2_beta <-
  rowMeans(
    
    luad_meth_analysis[
      ,
      luad_meth_group ==
        "PCD_C2",
      drop = FALSE
    ],
    
    na.rm =
      TRUE
    
  )


lusc_c1_beta <-
  rowMeans(
    
    lusc_meth_analysis[
      ,
      lusc_meth_group ==
        "PCD_C1",
      drop = FALSE
    ],
    
    na.rm =
      TRUE
    
  )


lusc_c2_beta <-
  rowMeans(
    
    lusc_meth_analysis[
      ,
      lusc_meth_group ==
        "PCD_C2",
      drop = FALSE
    ],
    
    na.rm =
      TRUE
    
  )


luad_meth_results$mean_beta_C1 <-
  luad_c1_beta[
    luad_meth_results$probe_id
  ]


luad_meth_results$mean_beta_C2 <-
  luad_c2_beta[
    luad_meth_results$probe_id
  ]


luad_meth_results$delta_beta <-
  luad_meth_results$mean_beta_C2 -
  luad_meth_results$mean_beta_C1


lusc_meth_results$mean_beta_C1 <-
  lusc_c1_beta[
    lusc_meth_results$probe_id
  ]


lusc_meth_results$mean_beta_C2 <-
  lusc_c2_beta[
    lusc_meth_results$probe_id
  ]


lusc_meth_results$delta_beta <-
  lusc_meth_results$mean_beta_C2 -
  lusc_meth_results$mean_beta_C1


# ============================================================
# 37. Primary differentially methylated probes
#
# Criteria:
#
#   BH FDR < 0.05
#   absolute delta-beta >= 0.10
# ============================================================

luad_dmps <-
  luad_meth_results |>
  dplyr::filter(
    
    adj.P.Val < 0.05,
    
    abs(
      delta_beta
    ) >= 0.10
    
  ) |>
  dplyr::arrange(
    
    adj.P.Val,
    
    dplyr::desc(
      abs(
        delta_beta
      )
    )
    
  )


lusc_dmps <-
  lusc_meth_results |>
  dplyr::filter(
    
    adj.P.Val < 0.05,
    
    abs(
      delta_beta
    ) >= 0.10
    
  ) |>
  dplyr::arrange(
    
    adj.P.Val,
    
    dplyr::desc(
      abs(
        delta_beta
      )
    )
    
  )


# ============================================================
# 38. DMP summaries
# ============================================================

luad_dmp_summary <- data.frame(
  
  histology =
    "LUAD",
  
  total =
    nrow(
      luad_dmps
    ),
  
  hypermethylated_C2 =
    sum(
      luad_dmps$delta_beta >=
        0.10
    ),
  
  hypomethylated_C2 =
    sum(
      luad_dmps$delta_beta <=
        -0.10
    )
  
)


lusc_dmp_summary <- data.frame(
  
  histology =
    "LUSC",
  
  total =
    nrow(
      lusc_dmps
    ),
  
  hypermethylated_C2 =
    sum(
      lusc_dmps$delta_beta >=
        0.10
    ),
  
  hypomethylated_C2 =
    sum(
      lusc_dmps$delta_beta <=
        -0.10
    )
  
)


dmp_summary <-
  dplyr::bind_rows(
    
    luad_dmp_summary,
    lusc_dmp_summary
    
  )


# ============================================================
# 39. Cross-histology DMP overlap
# ============================================================

shared_dmp_comparison <-
  dplyr::inner_join(
    
    luad_dmps |>
      dplyr::select(
        
        probe_id,
        
        LUAD_mean_beta_C1 =
          mean_beta_C1,
        
        LUAD_mean_beta_C2 =
          mean_beta_C2,
        
        LUAD_delta_beta =
          delta_beta,
        
        LUAD_padj =
          adj.P.Val
        
      ),
    
    lusc_dmps |>
      dplyr::select(
        
        probe_id,
        
        LUSC_mean_beta_C1 =
          mean_beta_C1,
        
        LUSC_mean_beta_C2 =
          mean_beta_C2,
        
        LUSC_delta_beta =
          delta_beta,
        
        LUSC_padj =
          adj.P.Val
        
      ),
    
    by =
      "probe_id"
    
  ) |>
  dplyr::mutate(
    
    same_direction =
      
      sign(
        LUAD_delta_beta
      ) ==
      
      sign(
        LUSC_delta_beta
      )
    
  )


shared_concordant_dmps <-
  shared_dmp_comparison |>
  dplyr::filter(
    same_direction
  )


# ============================================================
# 40. Save methylation results
# ============================================================

write.csv(
  
  luad_meth_selected,
  
  file.path(methylation_results_dir, "TCGA_LUAD_selected_HM450_samples.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_meth_selected,
  
  file.path(methylation_results_dir, "TCGA_LUSC_selected_HM450_samples.csv"),
  
  row.names = FALSE
  
)


# Full probe-level results are stored as RDS because they are
# large.

saveRDS(
  
  luad_meth_results,
  
  file.path(methylation_results_dir, "TCGA_LUAD_full_differential_methylation.rds")
  
)


saveRDS(
  
  lusc_meth_results,
  
  file.path(methylation_results_dir, "TCGA_LUSC_full_differential_methylation.rds")
  
)


# Primary biologically meaningful DMPs are also exported as CSV.

write.csv(
  
  luad_dmps,
  
  file.path(methylation_results_dir, "TCGA_LUAD_significant_DMPs.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  lusc_dmps,
  
  file.path(methylation_results_dir, "TCGA_LUSC_significant_DMPs.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  dmp_summary,
  
  file.path(methylation_results_dir, "TCGA_DMP_summary.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  shared_dmp_comparison,
  
  file.path(methylation_results_dir, "TCGA_LUAD_LUSC_shared_DMPs.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  shared_concordant_dmps,
  
  file.path(methylation_results_dir, "TCGA_LUAD_LUSC_concordant_DMPs.csv"),
  
  row.names = FALSE
  
)


# ============================================================
# 41. Final validation summary
# ============================================================

cat(
  "\n",
  "============================================================\n",
  "SCRIPT 17 COMPLETE\n",
  "============================================================\n"
)


cat(
  "\nMutation coverage:\n"
)

cat(
  "LUAD:",
  nrow(
    luad_mutation_coverage
  ),
  "patients\n"
)

cat(
  "LUSC:",
  nrow(
    lusc_mutation_coverage
  ),
  "patients\n"
)


cat(
  "\nSignificant mutation-associated genes (FDR < 0.05):\n"
)

cat(
  "LUAD:",
  sum(
    luad_mutation_fisher$padj <
      0.05
  ),
  "\n"
)

cat(
  "LUSC:",
  sum(
    lusc_mutation_fisher$padj <
      0.05
  ),
  "\n"
)


cat(
  "\nCNV patients:\n"
)

cat(
  "LUAD:",
  nrow(
    luad_cnv_burden
  ),
  "\n"
)

cat(
  "LUSC:",
  nrow(
    lusc_cnv_burden
  ),
  "\n"
)


cat(
  "\nMethylation analysis dimensions:\n"
)

cat(
  "LUAD:",
  nrow(
    luad_meth_analysis
  ),
  "probes x",
  ncol(
    luad_meth_analysis
  ),
  "patients\n"
)

cat(
  "LUSC:",
  nrow(
    lusc_meth_analysis
  ),
  "probes x",
  ncol(
    lusc_meth_analysis
  ),
  "patients\n"
)


cat(
  "\nDifferentially methylated probes:\n"
)

cat(
  "LUAD:",
  nrow(
    luad_dmps
  ),
  "\n"
)

cat(
  "LUSC:",
  nrow(
    lusc_dmps
  ),
  "\n"
)

cat(
  "Shared:",
  nrow(
    shared_dmp_comparison
  ),
  "\n"
)

cat(
  "Shared and concordant:",
  nrow(
    shared_concordant_dmps
  ),
  "\n"
)


cat(
  "\nAll Script 17 results saved successfully.\n"
)

cat(
  "\n========================================\n",
  "SCRIPT 17 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "Mutation outputs: ", mutation_results_dir, "\n",
  "CNV outputs: ", cnv_results_dir, "\n",
  "Methylation outputs: ", methylation_results_dir, "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================