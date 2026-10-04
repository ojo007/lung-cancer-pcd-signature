# ============================================================
# SCRIPT 31
# Final thesis summary tables
#
# IMPORTANT:
# - This script performs NO new statistical modelling.
# - It assembles results already validated in Scripts 01-30.
# - Headline values below are locked to the completed analyses.
# - Script 32 remains the single final figure-generation stage.
# ============================================================

source("03_scripts/00_project_config.R")

cat(
  "\n========================================\n",
  "SCRIPT 31\n",
  "FINAL THESIS SUMMARY TABLES\n",
  "========================================\n"
)

# ============================================================
# 1. Output directories
# ============================================================

table_dir_s31 <- "06_tables/final"
result_dir_s31 <- "04_results/final_summary"

dir.create(
  table_dir_s31,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  result_dir_s31,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Representative validated upstream-output checks
# ============================================================

required_files_s31 <- c(
  "04_results/clustering/TCGA_LUAD_PCD_cluster_assignments.csv",
  "04_results/clustering/TCGA_LUSC_PCD_cluster_assignments.csv",
  "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv",
  "04_results/mpcds/external_validation/TCGA_LUAD_MPCDS_standardized_weights_47.csv",
  "04_results/drug_sensitivity/DepMap_PRISM_MPCDS_all_drug_associations.csv",
  "04_results/drug_sensitivity/GDSC2/TCGA_LUAD_GDSC2_MPCDS_all_drug_associations.csv",
  "04_results/mpcds/immunotherapy_indicators/TCGA_LUAD_MPCDS_immunotherapy_indicator_correlations.csv",
  "04_results/immunotherapy/TCGA_LUAD_TIDE_overall_summary.csv",
  "04_results/mpcds/lasso_70_30_sensitivity/TCGA_LUAD_LASSO_70_30_final_summary.csv",
  "04_results/immune/additional_deconvolution_sensitivity/TCGA_LUAD_crossmethod_final_immune_summary.csv",
  "04_results/final_summary/TCGA_LUAD_47gene_PCD_subtype_PCA_variance.csv",
  "04_results/final_summary/TCGA_LUAD_Sankey_axis_summary.csv"
)

file_check_s31 <- data.frame(
  file = required_files_s31,
  exists = file.exists(required_files_s31),
  size_bytes = ifelse(
    file.exists(required_files_s31),
    file.info(required_files_s31)$size,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

cat(
  "\n========================================\n",
  "UPSTREAM OUTPUT CHECK\n",
  "========================================\n"
)

print(
  file_check_s31,
  row.names = FALSE
)

stopifnot(
  all(
    file_check_s31$exists
  )
)

# ============================================================
# 3. Cohort and PCD-discovery summary
# ============================================================

table_cohort_pcd_s31 <- data.frame(
  domain = c(
    "TCGA-LUAD differential expression",
    "TCGA-LUSC differential expression",
    "PCD catalogue",
    "LUAD significant PCD genes",
    "LUSC significant PCD genes",
    "Concordant shared candidate genes",
    "LUAD PCD subtypes",
    "LUSC PCD subtypes"
  ),
  result = c(
    "517 primary tumours; 59 normal lung samples",
    "501 primary tumours; 51 normal lung samples",
    "296 unique genes across apoptosis, necroptosis, pyroptosis, ferroptosis and cuproptosis",
    "57 genes: 36 upregulated and 21 downregulated",
    "95 genes: 61 upregulated and 34 downregulated",
    "47 genes: 31 upregulated and 16 downregulated",
    "PCD_C1 n=284; PCD_C2 n=233",
    "PCD_C1 n=286; PCD_C2 n=215"
  ),
  interpretation = c(
    "LUAD analysed independently",
    "LUSC analysed independently",
    "Curated PCD reference catalogue",
    "Tumour-versus-normal PCD differential expression",
    "Tumour-versus-normal PCD differential expression",
    "Concordant tumour-normal candidate set used for subtype analysis",
    "Two expression-derived molecular subtypes retained",
    "Two expression-derived molecular subtypes retained"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 4. Molecular subtype validation summary
# ============================================================

table_subtype_s31 <- data.frame(
  cohort = c(
    "TCGA-LUAD",
    "TCGA-LUSC",
    "GSE68465",
    "GSE72094",
    "GSE31210",
    "GSE50081",
    "GSE73403",
    "GSE4573",
    "GSE157010",
    "GSE30219"
  ),
  histology = c(
    "LUAD", "LUSC",
    "LUAD", "LUAD", "LUAD", "LUAD",
    "LUSC", "LUSC", "LUSC", "LUSC"
  ),
  n = c(
    517, 501,
    443, 442, 226, 128,
    69, 130, 235, 82
  ),
  PCD_C1 = c(
    284, 286,
    256, 248, 125, 74,
    32, 66, 125, 44
  ),
  PCD_C2 = c(
    233, 215,
    187, 194, 101, 54,
    37, 64, 110, 38
  ),
  mean_silhouette = c(
    0.17745570, 0.23453624,
    0.09052135, 0.10249271, 0.08065850, 0.08305747,
    0.10353985, 0.11563622, 0.11935578, 0.09966700
  ),
  key_mechanism_replication = c(
    "Necroptosis and pyroptosis differed between subtypes",
    "Ferroptosis strongest; apoptosis and cuproptosis also differed",
    "Strict: necroptosis; pyroptosis",
    "Strict: necroptosis; pyroptosis",
    "No strict mechanism replication",
    "No strict mechanism replication",
    "Strict: ferroptosis",
    "Strict: ferroptosis",
    "Strict: apoptosis; ferroptosis",
    "Strict: apoptosis; ferroptosis"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 5. LUAD MPCDS validation summary
# ============================================================

table_mpcds_s31 <- data.frame(
  cohort_analysis = c(
    "TCGA-LUAD internal nested CV",
    "TCGA-LUSC internal nested CV",
    "TCGA-LUAD final ridge score",
    "GSE68465 external validation",
    "GSE72094 external validation",
    "GSE31210 external validation",
    "GSE50081 external validation",
    "GSE68465 clinical-adjusted analysis"
  ),
  n = c(
    NA, NA, 504, 442, 398, 226, 128, 439
  ),
  metric = c(
    "Repeated nested-CV C-index",
    "Repeated nested-CV C-index",
    "Median-split HR high vs low",
    "C-index; HR per SD",
    "C-index; HR per SD",
    "C-index; HR per SD",
    "C-index; HR per SD",
    "Adjusted HR per SD; clinical-to-combined C-index"
  ),
  value = c(
    "0.6043",
    "0.5051",
    "HR 2.1962; log-rank p=1.82e-07",
    "C=0.5973; HR=1.2289; p=0.00225",
    "C=0.7016; HR=1.8544; p=4.21e-12",
    "C=0.5790; HR=1.3334; p=0.0811",
    "C=0.5937; HR=1.3833; p=0.0137",
    "HR=1.1727; p=0.0231; C 0.6996 -> 0.7146; delta=0.0150"
  ),
  interpretation = c(
    "Modest internally validated discrimination",
    "Chance-level performance; no defensible LUSC MPCDS",
    "Higher MPCDS associated with poorer overall survival",
    "Modest independent discrimination; PH violation",
    "Strongest external discrimination among tested GEO cohorts",
    "Directionally consistent but statistically imprecise; only 35 events",
    "Modest independent discrimination",
    "Association retained after clinical adjustment; incremental gain small"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 6. Time-dependent ROC and DCA summary
# ============================================================

table_prediction_s31 <- data.frame(
  analysis = c(
    "TCGA-LUAD timeROC",
    "GSE68465 timeROC",
    "GSE72094 timeROC",
    "GSE31210 timeROC",
    "GSE50081 timeROC",
    "TCGA-LUAD DCA / clinical model comparison"
  ),
  result = c(
    "AUC 1y=0.7054; 3y=0.7092; 5y=0.6493",
    "AUC 1y=0.6454; 3y=0.6375; 5y=0.5856",
    "AUC 1y=0.7314; 3y=0.7322; 5y=0.7207",
    "AUC 1y=0.7252; 3y=0.5507; 5y=0.5889",
    "AUC 1y=0.5785; 3y=0.5747; 5y=0.6190",
    "Complete-case n=475; clinical C=0.6741; clinical+MPCDS C=0.7243; LRT p=2.94e-09"
  ),
  caution = c(
    "Development cohort",
    "External validation",
    "5-year estimate unstable because of sparse support",
    "1-year estimate based on few deaths",
    "External validation",
    "Exploratory apparent decision-curve analysis; no external DCA"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 7. Therapy and immunotherapy-surrogate summary
# ============================================================

table_therapy_s31 <- data.frame(
  analysis = c(
    "PRISM pharmacogenomics",
    "GDSC2/oncoPredict MPCDS",
    "GDSC2/oncoPredict subtype",
    "GDSC2-PRISM concordance",
    "CD274 vs MPCDS",
    "Nonsynonymous mutation burden vs MPCDS",
    "MANTIS vs MPCDS",
    "TIDE predicted response",
    "TIDE no-benefits by subtype",
    "MPCDS vs TIDE"
  ),
  n_or_scope = c(
    "50 LUAD cell lines; 1,440 compounds tested",
    "504 patients; 198 drugs",
    "517 patients; 198 drugs",
    "72 testable strict-name overlaps",
    "504 patients",
    "492 patients",
    "501 patients",
    "517 patients",
    "517 patients",
    "504 patients"
  ),
  result = c(
    "0 compounds at FDR<0.05; minimum Spearman FDR=0.1893",
    "110 drug associations at FDR<0.05",
    "104 subtype drug associations at FDR<0.05",
    "37 concordant vs 35 discordant; binomial p=0.453",
    "Spearman rho=0.1281; p=0.00396",
    "Spearman rho=0.2677; p=1.61e-09",
    "Spearman rho=-0.0358; p=0.424",
    "193 predicted responders; 324 predicted nonresponders",
    "PCD_C1 4.58% vs PCD_C2 14.59%; OR=3.553; p=9.39e-05",
    "Spearman rho=0.3541; p=2.49e-16"
  ),
  interpretation = c(
    "Exploratory cell-line pharmacogenomic signals only",
    "Transcriptome-derived in-silico predicted IC50 associations",
    "Transcriptome-derived in-silico predicted IC50 associations",
    "No evidence of cross-platform directional concordance",
    "Small positive association",
    "Positive association; not formal clinical TMB",
    "No significant association",
    "TIDE-predicted status, not observed ICI response",
    "No-benefits classification more frequent in PCD_C2",
    "Higher MPCDS associated with less favourable in-silico TIDE profile"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 8. Sensitivity-analysis summary
# ============================================================

table_sensitivity_s31 <- data.frame(
  analysis = c(
    "70/30 LASSO-Cox sensitivity",
    "Additional immune deconvolution",
    "Immune cross-method subtype replication",
    "Immune cross-method MPCDS replication"
  ),
  result = c(
    paste0(
      "Train 352/127 events; test 152/55 events; selected genes ",
      "GCLC, EGLN3, DSG2, CLSPN; test C=0.5242; log-rank p=0.9887"
    ),
    "517 patients; 114 immune features; 113 testable; 45 subtype and 38 MPCDS features at FDR<0.05",
    "6 replicated populations; 3 direction conflicts",
    "4 replicated populations; 1 direction conflict"
  ),
  interpretation = c(
    "Sparse single-split model generalized poorly; ridge MPCDS remains primary",
    "xCell, quanTIseq, TIMER, EPIC and CIBERSORT/LM22 used as complementary sensitivity methods",
    "Cross-method replication supports selected subtype immune-contexture signals",
    "Replicated MPCDS associations were negative for B cells, CD4 T cells, Tregs and M2 macrophages"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 9. PCA and Sankey integrative summary
# ============================================================

table_integrative_s31 <- data.frame(
  metric = c(
    "PCA patients",
    "PCA genes",
    "PC1 variance percent",
    "PC2 variance percent",
    "PC1+PC2 cumulative percent",
    "Sankey complete-case patients",
    "Sankey PCD_C1",
    "Sankey PCD_C2",
    "Sankey Low MPCDS",
    "Sankey High MPCDS",
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  ),
  value = c(
    517,
    47,
    15.697042,
    11.770018,
    27.467060,
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

# ============================================================
# 10. Thesis headline-results table
# ============================================================

headline_results_s31 <- data.frame(
  topic = c(
    "PCD catalogue",
    "Concordant PCD candidate set",
    "Molecular subtypes",
    "External subtype validation",
    "Multi-omics characterization",
    "LUAD MPCDS internal performance",
    "LUSC MPCDS feasibility",
    "LUAD external MPCDS validation",
    "Time-dependent discrimination",
    "Decision-curve analysis",
    "PRISM pharmacogenomics",
    "GDSC2 predicted sensitivity",
    "Immunotherapy-related indicators",
    "TIDE surrogate analysis",
    "LASSO sensitivity analysis",
    "Immune deconvolution sensitivity",
    "Integrative PCA/Sankey"
  ),
  headline = c(
    "296 unique genes across five programmed-cell-death programmes",
    "47 concordantly dysregulated tumour-normal genes",
    "Two expression-derived PCD subtypes retained independently in LUAD and LUSC",
    "Eight independent GEO cohorts showed partial, mechanism-specific reproducibility",
    "Subtype differences extended across mutation, CNV, methylation, pathways and immune contexture",
    "Repeated nested-CV C-index 0.6043",
    "Repeated nested-CV C-index 0.5051; no defensible LUSC MPCDS",
    "Four LUAD GEO cohorts; adverse-risk direction in 4/4 and significant in 3/4",
    "External 1/3/5-year AUCs generally above 0.5, with sparse-horizon uncertainty",
    "Clinical+MPCDS C-index 0.7243 vs 0.6741 clinical-only in TCGA complete cases",
    "1,440 compounds tested; no FDR-significant associations",
    "110 MPCDS and 104 subtype associations at FDR<0.05; limited PRISM concordance",
    "CD274 and nonsynonymous mutation burden associated with MPCDS; MANTIS was not",
    "Higher MPCDS associated with higher TIDE; PCD_C2 had more TIDE no-benefits cases",
    "Four-gene sparse model showed held-out C-index 0.5242",
    "Five complementary deconvolution methods identified replicated and method-dependent immune signals",
    "47-gene PCA plus 496-patient subtype-MPCDS-stage Sankey dataset"
  ),
  interpretation_boundary = c(
    "Curated reference catalogue",
    "Not itself a prognostic signature",
    "Expression-derived groups; not causal entities",
    "External silhouettes were modest",
    "Multi-layer characterization, not a joint multi-omics integration model",
    "Modest discrimination",
    "Negative/weak finding retained",
    "Fixed TCGA weights; no outcome refitting",
    "Interpret sparse follow-up horizons cautiously",
    "Exploratory/apparent development-cohort analysis",
    "Cell-line pharmacogenomic evidence only",
    "Transcriptome-derived in-silico prediction, not observed patient response",
    "Mutation burden is not formal mutations-per-megabase TMB",
    "TIDE is an in-silico surrogate, not observed ICI response",
    "Sensitivity analysis only; does not replace primary ridge model",
    "Cross-method agreement prioritized over any single method",
    "Descriptive integrative visualization"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 11. Write final tables
# ============================================================

output_tables_s31 <- c(
  "Table_01_cohort_and_PCD_discovery_summary.csv",
  "Table_02_subtype_validation_summary.csv",
  "Table_03_MPCDS_validation_summary.csv",
  "Table_04_timeROC_and_DCA_summary.csv",
  "Table_05_therapy_and_immunotherapy_summary.csv",
  "Table_06_sensitivity_analyses_summary.csv",
  "Table_07_PCA_Sankey_integrative_summary.csv",
  "Table_final_headline_results.csv"
)

write.csv(
  table_cohort_pcd_s31,
  file.path(table_dir_s31, output_tables_s31[1]),
  row.names = FALSE
)

write.csv(
  table_subtype_s31,
  file.path(table_dir_s31, output_tables_s31[2]),
  row.names = FALSE
)

write.csv(
  table_mpcds_s31,
  file.path(table_dir_s31, output_tables_s31[3]),
  row.names = FALSE
)

write.csv(
  table_prediction_s31,
  file.path(table_dir_s31, output_tables_s31[4]),
  row.names = FALSE
)

write.csv(
  table_therapy_s31,
  file.path(table_dir_s31, output_tables_s31[5]),
  row.names = FALSE
)

write.csv(
  table_sensitivity_s31,
  file.path(table_dir_s31, output_tables_s31[6]),
  row.names = FALSE
)

write.csv(
  table_integrative_s31,
  file.path(table_dir_s31, output_tables_s31[7]),
  row.names = FALSE
)

write.csv(
  headline_results_s31,
  file.path(table_dir_s31, output_tables_s31[8]),
  row.names = FALSE
)

write.csv(
  headline_results_s31,
  file.path(result_dir_s31, "Thesis_headline_results.csv"),
  row.names = FALSE
)

# ============================================================
# 12. Table manifest
# ============================================================

table_manifest_s31 <- data.frame(
  table_number = seq_along(output_tables_s31),
  file = file.path(table_dir_s31, output_tables_s31),
  purpose = c(
    "Cohort composition, PCD catalogue and discovery summary",
    "Internal and external molecular-subtype validation",
    "Primary LUAD MPCDS development and external validation",
    "Time-dependent ROC and exploratory decision-curve summary",
    "Drug-sensitivity and immunotherapy-surrogate summary",
    "LASSO and immune-deconvolution sensitivity analyses",
    "PCA and subtype-MPCDS-stage Sankey summary",
    "Compact thesis-wide headline findings and interpretation boundaries"
  ),
  stringsAsFactors = FALSE
)

manifest_file_s31 <-
  file.path(
    table_dir_s31,
    "Table_manifest.csv"
  )

write.csv(
  table_manifest_s31,
  manifest_file_s31,
  row.names = FALSE
)

# ============================================================
# 13. Locked-value regression audit
# ============================================================

regression_audit_s31 <- data.frame(
  metric = c(
    "PCD catalogue genes",
    "Concordant candidate genes",
    "LUAD subtype patients",
    "LUSC subtype patients",
    "External subtype cohorts",
    "LUAD final MPCDS patients",
    "LUAD external MPCDS cohorts",
    "PRISM tested compounds",
    "TIDE patients",
    "LASSO held-out C-index",
    "Immune testable features",
    "Sankey complete cases",
    "PCA PC1 percent",
    "PCA PC2 percent"
  ),
  observed = c(
    296,
    47,
    sum(table_subtype_s31$n[table_subtype_s31$cohort == "TCGA-LUAD"]),
    sum(table_subtype_s31$n[table_subtype_s31$cohort == "TCGA-LUSC"]),
    sum(grepl("^GSE", table_subtype_s31$cohort)),
    504,
    4,
    1440,
    517,
    0.52420814,
    113,
    table_integrative_s31$value[table_integrative_s31$metric == "Sankey complete-case patients"],
    table_integrative_s31$value[table_integrative_s31$metric == "PC1 variance percent"],
    table_integrative_s31$value[table_integrative_s31$metric == "PC2 variance percent"]
  ),
  expected = c(
    296,
    47,
    517,
    501,
    8,
    504,
    4,
    1440,
    517,
    0.52420814,
    113,
    496,
    15.697042,
    11.770018
  ),
  tolerance = c(
    rep(0, 9),
    1e-8,
    0,
    0,
    1e-6,
    1e-6
  ),
  stringsAsFactors = FALSE
)

regression_audit_s31$matches_expected <-
  abs(
    regression_audit_s31$observed -
      regression_audit_s31$expected
  ) <=
  regression_audit_s31$tolerance

cat(
  "\n========================================\n",
  "FINAL SCRIPT 31 AUDIT\n",
  "========================================\n"
)

print(
  regression_audit_s31,
  row.names = FALSE
)

stopifnot(
  all(
    regression_audit_s31$matches_expected
  )
)

cat(
  "\nEstablished final-table regression checks: PASSED\n"
)

# ============================================================
# 14. Final output verification
# ============================================================

final_files_s31 <- c(
  file.path(table_dir_s31, output_tables_s31),
  manifest_file_s31,
  file.path(result_dir_s31, "Thesis_headline_results.csv")
)

final_file_check_s31 <- data.frame(
  file = final_files_s31,
  exists = file.exists(final_files_s31),
  size_bytes = file.info(final_files_s31)$size,
  stringsAsFactors = FALSE
)

cat(
  "\n========================================\n",
  "FINAL TABLE OUTPUT FILE CHECK\n",
  "========================================\n"
)

print(
  final_file_check_s31,
  row.names = FALSE
)

stopifnot(
  all(
    final_file_check_s31$exists
  )
)

stopifnot(
  all(
    final_file_check_s31$size_bytes > 0
  )
)

cat(
  "\n========================================\n",
  "SCRIPT 31 FINAL THESIS SUMMARY TABLES COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)