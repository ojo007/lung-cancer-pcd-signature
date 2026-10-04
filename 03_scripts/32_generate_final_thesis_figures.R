# ============================================================
# SCRIPT 32
# Final manuscript figure generation
#
# Purpose:
#   Reformat finalized analyses into consistent, colored,
#   publication-style multi-panel figures.
#
# IMPORTANT:
#   - No statistical models are refitted.
#   - No tests are rerun.
#   - Existing finalized outputs from Scripts 01-31 are used.
#   - PCD_C1 and PCD_C2 retain neutral names.
# ============================================================


# ------------------------------------------------------------
# 1. Project setup
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")


# ------------------------------------------------------------
# 2. Required packages
# ------------------------------------------------------------

required_packages <-
  c(
    "ggplot2",
    "dplyr",
    "patchwork",
    "ggalluvial"
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
      "Missing package(s):",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


library(
  ggplot2
)

library(
  dplyr
)

library(
  patchwork
)


# ------------------------------------------------------------
# 3. Canonical figure directories
# ------------------------------------------------------------

manuscript_figure_dir <-
  "05_figures/manuscript/main"

supplementary_figure_dir <-
  "05_figures/manuscript/supplementary"

individual_figure_dir <-
  "05_figures/manuscript/individual"

individual_category_dirs_s32 <-
  file.path(
    individual_figure_dir,
    c(
      "clustering",
      "immune",
      "enrichment",
      "mutation",
      "cnv",
      "methylation",
      "mpcds",
      "drug_sensitivity"
    )
  )


for (
  dir_s32 in c(
    manuscript_figure_dir,
    supplementary_figure_dir,
    individual_figure_dir,
    individual_category_dirs_s32
  )
) {
  dir.create(
    dir_s32,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 4. Consistent manuscript color system
# ------------------------------------------------------------

cluster_colors <-
  c(
    "PCD_C1" = "#2F6BFF",
    "PCD_C2" = "#E34A4A"
  )


histology_colors <-
  c(
    "LUAD" = "#2C7FB8",
    "LUSC" = "#F28E2B"
  )


mpcds_colors <-
  c(
    "Low" = "#2F6BFF",
    "High" = "#E34A4A"
  )


methylation_colors <-
  c(
    "Concordant:\nC2 hypomethylated" = "#2F6BFF",
    "Concordant:\nC2 hypermethylated" = "#E34A4A",
    "Discordant" = "#9E9E9E"
  )


# ------------------------------------------------------------
# 5. Shared manuscript theme
# ------------------------------------------------------------

theme_manuscript <-
  function(
    base_size = 11
  ) {
    
    theme_classic(
      base_size = base_size
    ) +
      
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = base_size + 2,
            hjust = 0
          ),
        
        plot.subtitle =
          element_text(
            size = base_size,
            margin =
              margin(
                b = 8
              )
          ),
        
        axis.title =
          element_text(
            face = "bold"
          ),
        
        axis.text =
          element_text(
            colour = "black"
          ),
        
        strip.background =
          element_blank(),
        
        strip.text =
          element_text(
            face = "bold",
            size = base_size + 1
          ),
        
        legend.title =
          element_text(
            face = "bold"
          ),
        
        legend.position =
          "top",
        
        plot.margin =
          margin(
            10,
            14,
            10,
            10
          )
      )
  }


# ------------------------------------------------------------
# 6. Load finalized biological data
# ------------------------------------------------------------

immune_robust <-
  read.csv(
    "04_results/immune/TCGA_LUAD_LUSC_MCPcounter_robust_populations.csv",
    check.names = FALSE
  )


hallmark_robust <-
  read.csv(
    "04_results/enrichment/TCGA_LUAD_LUSC_Hallmark_robust_pathways.csv",
    check.names = FALSE
  )


luad_mut_burden <-
  read.csv(
    "04_results/mutation/TCGA_LUAD_nonsynonymous_mutation_burden.csv",
    check.names = FALSE
  )


lusc_mut_burden <-
  read.csv(
    "04_results/mutation/TCGA_LUSC_nonsynonymous_mutation_burden.csv",
    check.names = FALSE
  )


mutation_tests <-
  read.csv(
    "04_results/mutation/TCGA_mutation_burden_cluster_tests.csv",
    check.names = FALSE
  )


mutation_shared <-
  read.csv(
    "04_results/mutation/TCGA_LUAD_LUSC_shared_significant_mutation_genes.csv",
    check.names = FALSE
  )


luad_cnv_summary <-
  read.csv(
    "04_results/cnv/TCGA_LUAD_CNV_burden_summary.csv",
    check.names = FALSE
  )


lusc_cnv_summary <-
  read.csv(
    "04_results/cnv/TCGA_LUSC_CNV_burden_summary.csv",
    check.names = FALSE
  )


luad_cnv_tests <-
  read.csv(
    "04_results/cnv/TCGA_LUAD_CNV_cluster_tests.csv",
    check.names = FALSE
  )


lusc_cnv_tests <-
  read.csv(
    "04_results/cnv/TCGA_LUSC_CNV_cluster_tests.csv",
    check.names = FALSE
  )


shared_dmps <-
  read.csv(
    "04_results/methylation/TCGA_LUAD_LUSC_shared_DMPs.csv",
    check.names = FALSE
  )


concordant_dmps <-
  read.csv(
    "04_results/methylation/TCGA_LUAD_LUSC_concordant_DMPs.csv",
    check.names = FALSE
  )


# ============================================================
# FIGURE 3
# Immune and functional characterization
# ============================================================


# ------------------------------------------------------------
# 7. Figure 3A - immune/stromal differences
# ------------------------------------------------------------

immune_plot_data <-
  data.frame(
    cell_population =
      rep(
        immune_robust$cell_population,
        2
      ),
    
    histology =
      rep(
        c(
          "LUAD",
          "LUSC"
        ),
        each =
          nrow(
            immune_robust
          )
      ),
    
    difference =
      c(
        immune_robust$LUAD_difference,
        immune_robust$LUSC_difference
      ),
    
    stringsAsFactors =
      FALSE
  )


immune_order <-
  immune_robust |>
  dplyr::mutate(
    mean_difference =
      (
        LUAD_difference +
          LUSC_difference
      ) / 2
  ) |>
  dplyr::arrange(
    mean_difference
  ) |>
  dplyr::pull(
    cell_population
  )


immune_plot_data$cell_population <-
  factor(
    immune_plot_data$cell_population,
    levels = immune_order
  )


plot_immune <-
  ggplot(
    immune_plot_data,
    aes(
      x = difference,
      y = cell_population,
      colour = histology,
      shape = histology
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 3.4,
    stroke = 0.8
  ) +
  
  scale_colour_manual(
    values = histology_colors
  ) +
  
  labs(
    title =
      "A  Immune and stromal infiltration",
    
    subtitle =
      "Negative values indicate greater abundance in PCD_C1",
    
    x =
      "MCP-counter difference (PCD_C2 - PCD_C1)",
    
    y =
      NULL,
    
    colour =
      "Histology",
    
    shape =
      "Histology"
  ) +
  
  theme_manuscript(
    base_size = 11
  )


# ------------------------------------------------------------
# 8. Figure 3B - Hallmark pathway differences
# ------------------------------------------------------------

hallmark_robust$mean_abs_NES <-
  rowMeans(
    abs(
      hallmark_robust[
        ,
        c(
          "LUAD_NES",
          "LUSC_NES"
        )
      ]
    )
  )


hallmark_top_C1 <-
  hallmark_robust |>
  dplyr::filter(
    enriched_cluster ==
      "PCD_C1"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      mean_abs_NES
    )
  ) |>
  dplyr::slice_head(
    n = 8
  )


hallmark_top_C2 <-
  hallmark_robust |>
  dplyr::filter(
    enriched_cluster ==
      "PCD_C2"
  ) |>
  dplyr::arrange(
    dplyr::desc(
      mean_abs_NES
    )
  ) |>
  dplyr::slice_head(
    n = 8
  )


hallmark_plot_data <-
  bind_rows(
    hallmark_top_C1,
    hallmark_top_C2
  )


hallmark_plot_data$pathway_label <-
  gsub(
    "^HALLMARK_",
    "",
    hallmark_plot_data$pathway
  )


hallmark_plot_data$pathway_label <-
  gsub(
    "_",
    " ",
    hallmark_plot_data$pathway_label
  )


hallmark_long <-
  data.frame(
    pathway =
      rep(
        hallmark_plot_data$pathway_label,
        2
      ),
    
    histology =
      rep(
        c(
          "LUAD",
          "LUSC"
        ),
        each =
          nrow(
            hallmark_plot_data
          )
      ),
    
    NES =
      c(
        hallmark_plot_data$LUAD_NES,
        hallmark_plot_data$LUSC_NES
      ),
    
    enriched_cluster =
      rep(
        hallmark_plot_data$enriched_cluster,
        2
      ),
    
    stringsAsFactors =
      FALSE
  )


hallmark_order <-
  hallmark_plot_data |>
  dplyr::mutate(
    mean_NES =
      (
        LUAD_NES +
          LUSC_NES
      ) / 2
  ) |>
  dplyr::arrange(
    mean_NES
  ) |>
  dplyr::pull(
    pathway_label
  )


hallmark_long$pathway <-
  factor(
    hallmark_long$pathway,
    levels = hallmark_order
  )


plot_hallmark <-
  ggplot(
    hallmark_long,
    aes(
      x = NES,
      y = pathway,
      colour = enriched_cluster,
      shape = histology
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 3,
    stroke = 0.8
  ) +
  
  scale_colour_manual(
    values = cluster_colors
  ) +
  
  labs(
    title =
      "B  Hallmark pathway enrichment",
    
    subtitle =
      "Negative NES: PCD_C1 enrichment; positive NES: PCD_C2 enrichment",
    
    x =
      "Normalized enrichment score (NES)",
    
    y =
      NULL,
    
    colour =
      "Enriched cluster",
    
    shape =
      "Histology"
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 9. Assemble Figure 3
# ------------------------------------------------------------

figure_3 <-
  plot_immune /
  plot_hallmark +
  
  plot_layout(
    heights =
      c(
        0.85,
        1.35
      )
  )


figure_3


ggsave(
  "05_figures/manuscript/main/Figure_03_immune_and_functional_characterization.png",
  figure_3,
  width = 10,
  height = 12,
  dpi = 300
)


ggsave(
  "05_figures/manuscript/main/Figure_03_immune_and_functional_characterization.pdf",
  figure_3,
  width = 10,
  height = 12
)


# ============================================================
# FIGURE 4
# Multi-omics characterization
# ============================================================


# ------------------------------------------------------------
# 10. Figure 4A - mutation burden
# ------------------------------------------------------------

mutation_plot_data <-
  bind_rows(
    luad_mut_burden |>
      mutate(
        histology = "LUAD"
      ),
    
    lusc_mut_burden |>
      mutate(
        histology = "LUSC"
      )
  )


mutation_plot_data$PCD_cluster <-
  factor(
    mutation_plot_data$PCD_cluster,
    levels =
      c(
        "PCD_C1",
        "PCD_C2"
      )
  )


luad_mutation_p <-
  mutation_tests$pvalue[
    mutation_tests$histology ==
      "LUAD"
  ]


lusc_mutation_p <-
  mutation_tests$pvalue[
    mutation_tests$histology ==
      "LUSC"
  ]


plot_mutation <-
  ggplot(
    mutation_plot_data,
    aes(
      x = PCD_cluster,
      y = log1p(
        n_nonsyn_mutations
      ),
      fill = PCD_cluster
    )
  ) +
  
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  
  geom_jitter(
    aes(
      colour = PCD_cluster
    ),
    width = 0.14,
    alpha = 0.20,
    size = 0.7,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~ histology,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = cluster_colors
  ) +
  
  scale_colour_manual(
    values = cluster_colors
  ) +
  
  labs(
    title =
      "A  Nonsynonymous mutation burden",
    
    subtitle =
      paste0(
        "LUAD p = ",
        format(
          luad_mutation_p,
          scientific = TRUE,
          digits = 2
        ),
        "; LUSC p = ",
        format(
          lusc_mutation_p,
          scientific = TRUE,
          digits = 2
        )
      ),
    
    x =
      NULL,
    
    y =
      "log(1 + mutation count)",
    
    fill =
      "PCD cluster"
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 11. Figure 4B - KEAP1 mutation frequency
# ------------------------------------------------------------

keap1 <-
  mutation_shared[
    mutation_shared$gene ==
      "KEAP1",
    ,
    drop = FALSE
  ]


keap1_plot_data <-
  data.frame(
    histology =
      rep(
        c(
          "LUAD",
          "LUSC"
        ),
        each = 2
      ),
    
    PCD_cluster =
      rep(
        c(
          "PCD_C1",
          "PCD_C2"
        ),
        2
      ),
    
    frequency =
      c(
        keap1$LUAD_C1_frequency,
        keap1$LUAD_C2_frequency,
        keap1$LUSC_C1_frequency,
        keap1$LUSC_C2_frequency
      ),
    
    stringsAsFactors =
      FALSE
  )


plot_keap1 <-
  ggplot(
    keap1_plot_data,
    aes(
      x = PCD_cluster,
      y = frequency,
      fill = PCD_cluster
    )
  ) +
  
  geom_col(
    width = 0.62
  ) +
  
  facet_wrap(
    ~ histology
  ) +
  
  scale_fill_manual(
    values = cluster_colors
  ) +
  
  scale_y_continuous(
    labels =
      function(
    x
      ) {
        
        paste0(
          round(
            x * 100
          ),
          "%"
        )
      }
  ) +
  
  labs(
    title =
      "B  KEAP1 mutation frequency",
    
    x =
      NULL,
    
    y =
      "Mutation frequency",
    
    fill =
      "PCD cluster"
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 12. Figure 4C - CNV burden
# ------------------------------------------------------------

cnv_plot_data <-
  bind_rows(
    luad_cnv_summary |>
      mutate(
        histology = "LUAD"
      ),
    
    lusc_cnv_summary |>
      mutate(
        histology = "LUSC"
      )
  )


plot_cnv <-
  ggplot(
    cnv_plot_data,
    aes(
      x = PCD_cluster,
      y = median_absolute,
      fill = PCD_cluster
    )
  ) +
  
  geom_col(
    width = 0.62
  ) +
  
  facet_wrap(
    ~ histology
  ) +
  
  scale_fill_manual(
    values = cluster_colors
  ) +
  
  labs(
    title =
      "C  Copy-number alteration burden",
    
    subtitle =
      "PCD_C2 shows greater absolute CNV burden in both histologies",
    
    x =
      NULL,
    
    y =
      "Median absolute CNV burden",
    
    fill =
      "PCD cluster"
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 13. Figure 4D - methylation concordance
# ------------------------------------------------------------

concordant_hyper <-
  sum(
    concordant_dmps$LUAD_delta_beta > 0 &
      concordant_dmps$LUSC_delta_beta > 0
  )


concordant_hypo <-
  sum(
    concordant_dmps$LUAD_delta_beta < 0 &
      concordant_dmps$LUSC_delta_beta < 0
  )


discordant <-
  sum(
    !shared_dmps$same_direction
  )


methylation_plot_data <-
  data.frame(
    category =
      c(
        "Concordant:\nC2 hypomethylated",
        "Concordant:\nC2 hypermethylated",
        "Discordant"
      ),
    
    n =
      c(
        concordant_hypo,
        concordant_hyper,
        discordant
      ),
    
    stringsAsFactors =
      FALSE
  )


methylation_plot_data$category <-
  factor(
    methylation_plot_data$category,
    levels =
      c(
        "Concordant:\nC2 hypomethylated",
        "Concordant:\nC2 hypermethylated",
        "Discordant"
      )
  )


plot_methylation <-
  ggplot(
    methylation_plot_data,
    aes(
      x = category,
      y = n,
      fill = category
    )
  ) +
  
  geom_col(
    width = 0.62
  ) +
  
  scale_fill_manual(
    values = methylation_colors
  ) +
  
  labs(
    title =
      "D  Cross-cohort methylation concordance",
    
    subtitle =
      paste0(
        sum(
          shared_dmps$same_direction
        ),
        " of ",
        nrow(
          shared_dmps
        ),
        " shared DMPs are directionally concordant"
      ),
    
    x =
      NULL,
    
    y =
      "Number of shared DMPs",
    
    fill =
      NULL
  ) +
  
  theme_manuscript(
    base_size = 10
  ) +
  
  theme(
    legend.position =
      "none",
    
    axis.text.x =
      element_text(
        size = 8.5
      )
  )


# ------------------------------------------------------------
# 14. Assemble Figure 4
# ------------------------------------------------------------

plot_mutation <-
  plot_mutation +
  theme(
    legend.position = "none"
  )


plot_keap1 <-
  plot_keap1 +
  theme(
    legend.position = "none"
  )


plot_cnv <-
  plot_cnv +
  theme(
    legend.position = "none"
  )


plot_methylation <-
  plot_methylation +
  theme(
    legend.position = "none"
  )


figure_4 <-
  (
    plot_mutation |
      plot_keap1
  ) /
  (
    plot_cnv |
      plot_methylation
  )


figure_4


ggsave(
  "05_figures/manuscript/main/Figure_04_multiomics_characterization.png",
  figure_4,
  width = 12,
  height = 9,
  dpi = 300
)


ggsave(
  "05_figures/manuscript/main/Figure_04_multiomics_characterization.pdf",
  figure_4,
  width = 12,
  height = 9
)


# ------------------------------------------------------------
# 15. Verify first manuscript-style figures
# ------------------------------------------------------------

first_manuscript_figures <-
  c(
    "05_figures/manuscript/main/Figure_03_immune_and_functional_characterization.png",
    "05_figures/manuscript/main/Figure_03_immune_and_functional_characterization.pdf",
    "05_figures/manuscript/main/Figure_04_multiomics_characterization.png",
    "05_figures/manuscript/main/Figure_04_multiomics_characterization.pdf"
  )


cat(
  "\nFirst manuscript-style figures created:\n"
)


print(
  file.exists(
    first_manuscript_figures
  )
)


# ============================================================
# FIGURE 5
# Construction and internal validation of LUAD MPCDS
# ============================================================


# ------------------------------------------------------------
# 16. Additional package for Kaplan-Meier estimation
# ------------------------------------------------------------

if (
  !requireNamespace(
    "survival",
    quietly = TRUE
  )
) {
  
  stop(
    "Package 'survival' is required."
  )
}


library(
  survival
)


# ------------------------------------------------------------
# 17. Load finalized MPCDS data
# ------------------------------------------------------------

luad_repeat_performance <-
  read.csv(
    "04_results/mpcds/TCGA_LUAD_MPCDS_repeat_performance.csv",
    check.names = FALSE
  )


lusc_repeat_performance <-
  read.csv(
    "04_results/mpcds/TCGA_LUSC_MPCDS_repeat_performance.csv",
    check.names = FALSE
  )


luad_coefficients <-
  read.csv(
    "04_results/mpcds/TCGA_LUAD_MPCDS_ridge_coefficients.csv",
    check.names = FALSE
  )


clinical_repeat_summary <-
  read.csv(
    "04_results/mpcds/TCGA_LUAD_clinical_MPCDS_repeat_summary.csv",
    check.names = FALSE
  )


incremental_summary <-
  read.csv(
    "04_results/mpcds/TCGA_LUAD_MPCDS_incremental_prediction_summary.csv",
    check.names = FALSE
  )


tcga_scores <-
  read.csv(
    "04_results/mpcds/TCGA_LUAD_MPCDS_patient_scores.csv",
    check.names = FALSE
  )


tcga_km_summary <-
  read.csv(
    "04_results/mpcds/TCGA_LUAD_MPCDS_KM_summary.csv",
    check.names = FALSE
  )


external_scores <-
  read.csv(
    paste0(
      "04_results/mpcds/external_validation/",
      "GSE68465_MPCDS_patient_scores.csv"
    ),
    check.names = FALSE
  )


external_km_summary <-
  read.csv(
    paste0(
      "04_results/mpcds/external_validation/",
      "GSE68465_KM_summary.csv"
    ),
    check.names = FALSE
  )


external_primary <-
  read.csv(
    paste0(
      "04_results/mpcds/external_validation/",
      "GSE68465_primary_validation_summary.csv"
    ),
    check.names = FALSE
  )


external_adjusted <-
  read.csv(
    paste0(
      "04_results/mpcds/external_validation/",
      "GSE68465_clinical_adjusted_MPCDS_summary.csv"
    ),
    check.names = FALSE
  )


external_time_hr <-
  read.csv(
    paste0(
      "04_results/mpcds/external_validation/",
      "GSE68465_time_specific_HR_adjusted.csv"
    ),
    check.names = FALSE
  )


# Four independent LUAD GEO cohorts from Script 21. These files are
# distinct from the legacy GSE68465-only clinical-adjustment outputs below.
geo_accessions_s32 <-
  c("GSE68465", "GSE72094", "GSE31210", "GSE50081")

geo_validation_dir_s32 <-
  "04_results/mpcds/external_validation"

geo_validation_summary_file_s32 <-
  file.path(
    geo_validation_dir_s32,
    "LUAD_MPCDS_external_GEO_validation_summary.csv"
  )

geo_score_files_s32 <-
  file.path(
    geo_validation_dir_s32,
    geo_accessions_s32,
    paste0(geo_accessions_s32, "_MPCDS_patient_scores.csv")
  )

stopifnot(
  file.exists(geo_validation_summary_file_s32),
  all(file.exists(geo_score_files_s32))
)

geo_validation_summary_s32 <-
  read.csv(
    geo_validation_summary_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

required_geo_summary_columns_s32 <-
  c(
    "accession", "validation_status", "n_OS_samples", "n_events",
    "C_index", "C_index_lower95", "C_index_upper95",
    "continuous_HR_per_SD", "continuous_HR_lower95",
    "continuous_HR_upper95", "continuous_p"
  )

stopifnot(
  all(required_geo_summary_columns_s32 %in%
        names(geo_validation_summary_s32)),
  nrow(geo_validation_summary_s32) == 4L,
  setequal(geo_validation_summary_s32$accession, geo_accessions_s32),
  all(geo_validation_summary_s32$validation_status == "SUCCESS")
)

geo_validation_summary_s32 <-
  geo_validation_summary_s32[
    match(geo_accessions_s32, geo_validation_summary_s32$accession),
    ,
    drop = FALSE
  ]

geo_scores_s32 <-
  setNames(
    lapply(
      geo_score_files_s32,
      function(score_file_s32) {
        read.csv(
          score_file_s32,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    ),
    geo_accessions_s32
  )

for (accession_s32 in geo_accessions_s32) {
  scores_s32 <- geo_scores_s32[[accession_s32]]
  summary_s32 <- geo_validation_summary_s32[
    geo_validation_summary_s32$accession == accession_s32,
    ,
    drop = FALSE
  ]
  stopifnot(
    all(c("OS_time", "OS_status", "risk_group") %in%
          names(scores_s32)),
    nrow(scores_s32) == summary_s32$n_OS_samples,
    sum(scores_s32$OS_status == 1L) == summary_s32$n_events,
    all(is.finite(scores_s32$OS_time)),
    all(scores_s32$OS_time >= 0),
    all(scores_s32$OS_status %in% c(0L, 1L)),
    all(scores_s32$risk_group %in%
          c("Low MPCDS", "High MPCDS"))
  )
}


# ------------------------------------------------------------
# 18. Helper to convert survfit object into ggplot data
# ------------------------------------------------------------

survfit_to_df <-
  function(
    fit
  ) {
    
    fit_summary <-
      summary(
        fit
      )
    
    
    strata_names <-
      if (
        is.null(
          fit_summary$strata
        )
      ) {
        
        rep(
          "All",
          length(
            fit_summary$time
          )
        )
        
      } else {
        
        as.character(
          fit_summary$strata
        )
      }
    
    
    data.frame(
      time =
        fit_summary$time,
      
      survival =
        fit_summary$surv,
      
      lower =
        fit_summary$lower,
      
      upper =
        fit_summary$upper,
      
      strata =
        strata_names,
      
      stringsAsFactors =
        FALSE
    )
  }


# ------------------------------------------------------------
# 19. Figure 5A - repeated nested CV performance
# ------------------------------------------------------------

cv_plot_data <-
  bind_rows(
    luad_repeat_performance |>
      mutate(
        histology =
          "LUAD"
      ),
    
    lusc_repeat_performance |>
      mutate(
        histology =
          "LUSC"
      )
  )


cv_plot_data$histology <-
  factor(
    cv_plot_data$histology,
    levels =
      c(
        "LUAD",
        "LUSC"
      )
  )


plot_cv <-
  ggplot(
    cv_plot_data,
    aes(
      x = histology,
      y = mean_C,
      fill = histology
    )
  ) +
  
  geom_hline(
    yintercept = 0.5,
    linetype = 2,
    linewidth = 0.5
  ) +
  
  geom_boxplot(
    width = 0.55,
    alpha = 0.75,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    aes(
      colour = histology
    ),
    width = 0.10,
    size = 2,
    alpha = 0.80,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = histology_colors
  ) +
  
  scale_colour_manual(
    values = histology_colors
  ) +
  
  scale_y_continuous(
    limits =
      c(
        0.45,
        0.65
      )
  ) +
  
  labs(
    title =
      "A  Repeated nested-CV performance",
    
    subtitle =
      paste0(
        "LUAD mean C = ",
        sprintf(
          "%.3f",
          mean(
            luad_repeat_performance$mean_C
          )
        ),
        "; LUSC mean C = ",
        sprintf(
          "%.3f",
          mean(
            lusc_repeat_performance$mean_C
          )
        )
      ),
    
    x =
      NULL,
    
    y =
      "Mean C-index",
    
    fill =
      "Histology"
  ) +
  
  theme_manuscript(
    base_size = 10
  ) +
  
  theme(
    legend.position =
      "none"
  )


# ------------------------------------------------------------
# 20. Figure 5B - top ridge coefficients
# ------------------------------------------------------------

coefficient_plot_data <-
  luad_coefficients |>
  mutate(
    abs_coefficient =
      abs(
        coefficient
      ),
    
    direction =
      ifelse(
        coefficient > 0,
        "Positive",
        "Negative"
      )
  ) |>
  arrange(
    desc(
      abs_coefficient
    )
  ) |>
  slice_head(
    n = 20
  )


coefficient_plot_data$gene <-
  factor(
    coefficient_plot_data$gene,
    levels =
      rev(
        coefficient_plot_data$gene
      )
  )


coefficient_colors <-
  c(
    "Negative" = "#2F6BFF",
    "Positive" = "#E34A4A"
  )


plot_coefficients <-
  ggplot(
    coefficient_plot_data,
    aes(
      x = coefficient,
      y = gene,
      fill = direction
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +
  
  geom_col(
    width = 0.68
  ) +
  
  scale_fill_manual(
    values = coefficient_colors
  ) +
  
  labs(
    title =
      "B  Final LUAD ridge coefficients",
    
    subtitle =
      "Top 20 absolute coefficients shown; final model contains all 47 PCD genes",
    
    x =
      "Ridge Cox coefficient",
    
    y =
      NULL,
    
    fill =
      "Coefficient"
  ) +
  
  theme_manuscript(
    base_size = 9
  ) +
  
  theme(
    legend.position =
      "top"
  )


# ------------------------------------------------------------
# 21. Figure 5C - TCGA-LUAD Kaplan-Meier
# ------------------------------------------------------------

tcga_scores$risk_group <-
  factor(
    tcga_scores$risk_group,
    levels =
      c(
        "Low MPCDS",
        "High MPCDS"
      )
  )


tcga_survfit <-
  survival::survfit(
    survival::Surv(
      OS_time,
      OS_status
    ) ~ risk_group,
    data = tcga_scores
  )


tcga_km_data <-
  survfit_to_df(
    tcga_survfit
  )


tcga_km_data$risk_group <-
  sub(
    "^risk_group=",
    "",
    tcga_km_data$strata
  )


tcga_km_data$risk_group <-
  factor(
    tcga_km_data$risk_group,
    levels =
      c(
        "Low MPCDS",
        "High MPCDS"
      )
  )


risk_colors <-
  c(
    "Low MPCDS" = "#2F6BFF",
    "High MPCDS" = "#E34A4A"
  )


tcga_logrank_p <-
  unique(
    tcga_km_summary$logrank_p
  )[1]


plot_tcga_km <-
  ggplot(
    tcga_km_data,
    aes(
      x = time / 365.25,
      y = survival,
      colour = risk_group
    )
  ) +
  
  geom_step(
    linewidth = 1
  ) +
  
  scale_colour_manual(
    values = risk_colors
  ) +
  
  scale_y_continuous(
    limits =
      c(
        0,
        1
      ),
    labels =
      function(
    x
      ) {
        
        paste0(
          round(
            100 * x
          ),
          "%"
        )
      }
  ) +
  
  labs(
    title =
      "C  TCGA-LUAD survival by MPCDS",
    
    subtitle =
      paste0(
        "Development-cohort median split; log-rank p = ",
        format(
          tcga_logrank_p,
          scientific = TRUE,
          digits = 2
        )
      ),
    
    x =
      "Overall survival (years)",
    
    y =
      "Survival probability",
    
    colour =
      "MPCDS group"
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 22. Figure 5D - incremental prediction across repeats
# ------------------------------------------------------------

incremental_long <-
  bind_rows(
    clinical_repeat_summary |>
      transmute(
        repeat_id =
          repeat_id,
        
        model =
          "Clinical",
        
        C_index =
          mean_clinical_C
      ),
    
    clinical_repeat_summary |>
      transmute(
        repeat_id =
          repeat_id,
        
        model =
          "Clinical + MPCDS",
        
        C_index =
          mean_combined_C
      )
  )


incremental_long$model <-
  factor(
    incremental_long$model,
    levels =
      c(
        "Clinical",
        "Clinical + MPCDS"
      )
  )


incremental_model_colors <-
  c(
    "Clinical" = "#8A8A8A",
    "Clinical + MPCDS" = "#E34A4A"
  )


plot_incremental <-
  ggplot(
    incremental_long,
    aes(
      x = model,
      y = C_index,
      group = repeat_id
    )
  ) +
  
  geom_line(
    colour = "#B0B0B0",
    alpha = 0.7,
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      colour = model
    ),
    size = 2.6
  ) +
  
  scale_colour_manual(
    values = incremental_model_colors
  ) +
  
  labs(
    title =
      "D  Incremental prognostic performance",
    
    subtitle =
      paste0(
        incremental_summary$repeats_positive,
        "/",
        incremental_summary$repeats_total,
        " repeats improved; mean delta C = ",
        sprintf(
          "%.3f",
          incremental_summary$overall_delta_C
        )
      ),
    
    x =
      NULL,
    
    y =
      "Mean C-index",
    
    colour =
      NULL
  ) +
  
  theme_manuscript(
    base_size = 10
  ) +
  
  theme(
    legend.position =
      "top"
  )


# ------------------------------------------------------------
# 23. Assemble Figure 5
# ------------------------------------------------------------

figure_5 <-
  (
    plot_cv |
      plot_coefficients
  ) /
  (
    plot_tcga_km |
      plot_incremental
  )


figure_5


ggsave(
  "05_figures/manuscript/main/Figure_05_MPCDS_construction_internal_validation.png",
  figure_5,
  width = 12,
  height = 10,
  dpi = 300
)


ggsave(
  "05_figures/manuscript/main/Figure_05_MPCDS_construction_internal_validation.pdf",
  figure_5,
  width = 12,
  height = 10
)


# ============================================================
# FIGURE 6
# Independent external validation of MPCDS
# ============================================================


# ------------------------------------------------------------
# 24. Figure 6A - GSE68465 Kaplan-Meier
# ------------------------------------------------------------

external_scores$risk_group <-
  factor(
    external_scores$risk_group,
    levels =
      c(
        "Low MPCDS",
        "High MPCDS"
      )
  )


external_survfit <-
  survival::survfit(
    survival::Surv(
      OS_months,
      OS_status
    ) ~ risk_group,
    data = external_scores
  )


external_km_data <-
  survfit_to_df(
    external_survfit
  )


external_km_data$risk_group <-
  sub(
    "^risk_group=",
    "",
    external_km_data$strata
  )


external_km_data$risk_group <-
  factor(
    external_km_data$risk_group,
    levels =
      c(
        "Low MPCDS",
        "High MPCDS"
      )
  )


external_logrank_p <-
  unique(
    external_km_summary$logrank_p
  )[1]


plot_external_km <-
  ggplot(
    external_km_data,
    aes(
      x = time / 12,
      y = survival,
      colour = risk_group
    )
  ) +
  
  geom_step(
    linewidth = 1
  ) +
  
  scale_colour_manual(
    values = risk_colors
  ) +
  
  scale_y_continuous(
    limits =
      c(
        0,
        1
      ),
    labels =
      function(
    x
      ) {
        
        paste0(
          round(
            100 * x
          ),
          "%"
        )
      }
  ) +
  
  labs(
    title =
      "GSE68465 survival by MPCDS",
    
    subtitle =
      paste0(
        "Median-group log-rank p = ",
        format(
          external_logrank_p,
          scientific = TRUE,
          digits = 2
        ),
        "; continuous C-index = ",
        sprintf(
          "%.3f",
          external_primary$C_index
        )
      ),
    
    x =
      "Overall survival (years)",
    
    y =
      "Survival probability",
    
    colour =
      "MPCDS group"
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 25. Figure 6A - one Kaplan-Meier panel for each GEO cohort
# ------------------------------------------------------------

make_geo_km_s32 <-
  function(accession_s32, panel_letter_s32) {
    scores_s32 <- geo_scores_s32[[accession_s32]]
    info_s32 <- geo_validation_summary_s32[
      geo_validation_summary_s32$accession == accession_s32,
      ,
      drop = FALSE
    ]
    scores_s32$risk_group <-
      factor(
        scores_s32$risk_group,
        levels = c("Low MPCDS", "High MPCDS")
      )
    
    fit_s32 <-
      survival::survfit(
        survival::Surv(OS_time, OS_status) ~ risk_group,
        data = scores_s32
      )
    km_s32 <- survfit_to_df(fit_s32)
    km_s32$risk_group <-
      factor(
        sub("^risk_group=", "", km_s32$strata),
        levels = c("Low MPCDS", "High MPCDS")
      )
    
    ggplot(
      km_s32,
      aes(x = time / 12, y = survival, colour = risk_group)
    ) +
      geom_step(linewidth = 0.85) +
      scale_colour_manual(values = risk_colors, drop = FALSE) +
      scale_y_continuous(limits = c(0, 1), labels =
                           function(x) paste0(round(100 * x), "%")) +
      labs(
        title = paste0(panel_letter_s32, "  ", accession_s32),
        subtitle = paste0(
          "n = ", info_s32$n_OS_samples,
          "; deaths = ", info_s32$n_events,
          "; C-index = ", sprintf("%.3f", info_s32$C_index)
        ),
        x = "Overall survival (years)",
        y = "Survival probability",
        colour = "MPCDS group"
      ) +
      theme_manuscript(base_size = 9)
  }

geo_km_panels_s32 <-
  Map(make_geo_km_s32, geo_accessions_s32, c("A", "B", "C", "D"))

# ------------------------------------------------------------
# 25B. Figure 6E - continuous MPCDS associations in all four
#      independent GEO cohorts; no refitting occurs here.
# ------------------------------------------------------------

geo_forest_s32 <-
  geo_validation_summary_s32

geo_forest_s32$accession <-
  factor(
    geo_forest_s32$accession,
    levels = rev(geo_accessions_s32)
  )

plot_forest <-
  ggplot(
    geo_forest_s32,
    aes(x = continuous_HR_per_SD, y = accession)
  ) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.5) +
  geom_errorbar(
    aes(
      xmin = continuous_HR_lower95,
      xmax = continuous_HR_upper95
    ),
    width = 0.18,
    linewidth = 0.8,
    orientation = "y",
    colour = "#2C7FB8"
  ) +
  geom_point(size = 3, colour = "#2C7FB8") +
  labs(
    title = "E  Four independent GEO cohorts",
    subtitle = "Continuous MPCDS HR per 1-SD increase (95% CI)",
    x = "Hazard ratio (95% CI)",
    y = NULL
  ) +
  theme_manuscript(base_size = 10)


# ------------------------------------------------------------
# 26. Figure 6C - external incremental C-index
# ------------------------------------------------------------

external_cindex_data <-
  data.frame(
    model =
      factor(
        c(
          "Clinical",
          "Clinical + MPCDS"
        ),
        levels =
          c(
            "Clinical",
            "Clinical + MPCDS"
          )
      ),
    
    C_index =
      c(
        external_adjusted$clinical_C,
        external_adjusted$clinical_plus_MPCDS_C
      )
  )


plot_external_cindex <-
  ggplot(
    external_cindex_data,
    aes(
      x = model,
      y = C_index,
      fill = model
    )
  ) +
  
  geom_col(
    width = 0.60
  ) +
  
  geom_text(
    aes(
      label =
        sprintf(
          "%.3f",
          C_index
        )
    ),
    vjust = -0.5,
    size = 4
  ) +
  
  scale_fill_manual(
    values = incremental_model_colors
  ) +
  
  scale_y_continuous(
    limits =
      c(
        0,
        0.80
      )
  ) +
  
  labs(
    title =
      "F  GSE68465 incremental discrimination",
    
    subtitle =
      paste0(
        "Delta C = ",
        sprintf(
          "%.3f",
          external_adjusted$delta_C
        ),
        "; LRT p = ",
        sprintf(
          "%.3f",
          external_adjusted$LRT_p_value
        )
      ),
    
    x =
      NULL,
    
    y =
      "C-index",
    
    fill =
      NULL
  ) +
  
  theme_manuscript(
    base_size = 10
  ) +
  
  theme(
    legend.position =
      "none"
  )


# ------------------------------------------------------------
# 27. Figure 6D - time-specific adjusted HR
# ------------------------------------------------------------

external_time_hr$time_label <-
  paste0(
    external_time_hr$years,
    "-year"
  )


external_time_hr$time_label <-
  factor(
    external_time_hr$time_label,
    levels =
      c(
        "1-year",
        "3-year",
        "5-year"
      )
  )


plot_time_hr <-
  ggplot(
    external_time_hr,
    aes(
      x = time_label,
      y = HR
    )
  ) +
  
  geom_hline(
    yintercept = 1,
    linetype = 2,
    linewidth = 0.5
  ) +
  
  geom_errorbar(
    aes(
      ymin = CI_lower,
      ymax = CI_upper
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 3.3,
    colour = "#E34A4A"
  ) +
  
  labs(
    title =
      "G  GSE68465 time-varying MPCDS effect",
    
    subtitle =
      "Adjusted HR/SD is strongest early and attenuates over time",
    
    x =
      "Time after diagnosis",
    
    y =
      "Adjusted hazard ratio"
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 28. Assemble Figure 6
# ------------------------------------------------------------

figure_6 <-
  patchwork::wrap_plots(
    A = geo_km_panels_s32[[1]],
    B = geo_km_panels_s32[[2]],
    C = geo_km_panels_s32[[3]],
    D = geo_km_panels_s32[[4]],
    E = plot_forest,
    F = plot_external_cindex,
    G = plot_time_hr,
    design = "AB\nCD\nEE\nFG"
  ) +
  patchwork::plot_layout(heights = c(1.2, 1.2, 0.75, 1))


figure_6


ggsave(
  "05_figures/manuscript/main/Figure_06_external_MPCDS_validation.png",
  figure_6,
  width = 13,
  height = 14,
  dpi = 300
)


ggsave(
  "05_figures/manuscript/main/Figure_06_external_MPCDS_validation.pdf",
  figure_6,
  width = 13,
  height = 14
)


# ------------------------------------------------------------
# 29. Verify Figures 5 and 6
# ------------------------------------------------------------

mpcds_manuscript_figures <-
  c(
    "05_figures/manuscript/main/Figure_05_MPCDS_construction_internal_validation.png",
    "05_figures/manuscript/main/Figure_05_MPCDS_construction_internal_validation.pdf",
    "05_figures/manuscript/main/Figure_06_external_MPCDS_validation.png",
    "05_figures/manuscript/main/Figure_06_external_MPCDS_validation.pdf"
  )


cat(
  "\nMPCDS manuscript-style figures created:\n"
)


print(
  file.exists(
    mpcds_manuscript_figures
  )
)


# ============================================================
# FIGURE 7
# PRISM pharmacogenomic associations
# ============================================================


# ------------------------------------------------------------
# 30. Load finalized PRISM outputs
# ------------------------------------------------------------

prism_all <-
  read.csv(
    "04_results/drug_sensitivity/DepMap_PRISM_MPCDS_all_drug_associations.csv",
    check.names = FALSE
  )


prism_nominal_0001 <-
  read.csv(
    "04_results/drug_sensitivity/DepMap_PRISM_MPCDS_nominal_p_lt_0.001.csv",
    check.names = FALSE
  )


prism_summary <-
  read.csv(
    "04_results/drug_sensitivity/DepMap_PRISM_MPCDS_analysis_summary.csv",
    check.names = FALSE
  )


prism_sensitivity_top <-
  read.csv(
    "04_results/drug_sensitivity/DepMap_PRISM_top_20_sensitivity_associations.csv",
    check.names = FALSE
  )


prism_resistance_top <-
  read.csv(
    "04_results/drug_sensitivity/DepMap_PRISM_top_20_resistance_associations.csv",
    check.names = FALSE
  )


get_prism_metric <-
  function(
    metric_name
  ) {
    
    x <-
      prism_summary$value[
        prism_summary$metric ==
          metric_name
      ]
    
    if (
      length(
        x
      ) != 1
    ) {
      
      stop(
        paste(
          "PRISM metric not uniquely found:",
          metric_name
        )
      )
    }
    
    x
  }


# ------------------------------------------------------------
# 31. PRISM headline counts
# ------------------------------------------------------------

prism_tested <-
  get_prism_metric(
    "Spearman-tested compounds"
  )


prism_nominal_005 <-
  get_prism_metric(
    "Spearman nominal p < 0.05"
  )


prism_nominal_001 <-
  get_prism_metric(
    "Spearman nominal p < 0.01"
  )


prism_nominal_0001_n <-
  get_prism_metric(
    "Spearman nominal p < 0.001"
  )


prism_fdr_005 <-
  get_prism_metric(
    "Spearman FDR < 0.05"
  )


prism_min_fdr <-
  get_prism_metric(
    "Minimum Spearman FDR"
  )


# ------------------------------------------------------------
# 32. Filter actually tested compounds
# ------------------------------------------------------------

prism_tested_data <-
  prism_all |>
  dplyr::filter(
    !is.na(
      spearman_rho
    ),
    !is.na(
      p_value
    ),
    !is.na(
      FDR
    )
  )


stopifnot(
  nrow(
    prism_tested_data
  ) ==
    prism_tested
)


# ------------------------------------------------------------
# 33. Association direction colors
# ------------------------------------------------------------

prism_direction_colors <-
  c(
    "Greater sensitivity" = "#2F6BFF",
    "Greater resistance" = "#E34A4A"
  )


prism_tested_data$direction_short <-
  ifelse(
    prism_tested_data$spearman_rho < 0,
    "Greater sensitivity",
    "Greater resistance"
  )


# ------------------------------------------------------------
# 34. Figure 7A - global PRISM association landscape
# ------------------------------------------------------------

prism_tested_data$neg_log10_p <-
  -log10(
    prism_tested_data$p_value
  )


plot_prism_global <-
  ggplot(
    prism_tested_data,
    aes(
      x = spearman_rho,
      y = neg_log10_p,
      colour = direction_short
    )
  ) +
  
  geom_hline(
    yintercept =
      -log10(
        0.05
      ),
    linetype = 2,
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +
  
  geom_point(
    alpha = 0.55,
    size = 1.6
  ) +
  
  scale_colour_manual(
    values = prism_direction_colors
  ) +
  
  labs(
    title =
      "A  Global PRISM association landscape",
    
    subtitle =
      paste0(
        as.integer(
          prism_tested
        ),
        " compounds tested; no association reached BH FDR < 0.05"
      ),
    
    x =
      "Spearman correlation: MPCDS vs PRISM AUC",
    
    y =
      "-log10(nominal p-value)",
    
    colour =
      "Higher MPCDS associated with"
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 35. Figure 7B - strongest nominal associations
# ------------------------------------------------------------

top_nominal_plot_data <-
  prism_nominal_0001 |>
  dplyr::arrange(
    spearman_rho
  )


top_nominal_plot_data$drug_name <-
  factor(
    top_nominal_plot_data$drug_name,
    levels =
      top_nominal_plot_data$drug_name
  )


plot_prism_top <-
  ggplot(
    top_nominal_plot_data,
    aes(
      x = spearman_rho,
      y = drug_name
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = spearman_rho,
      y = drug_name,
      yend = drug_name
    ),
    linewidth = 1,
    colour = "#2F6BFF"
  ) +
  
  geom_point(
    size = 3.5,
    colour = "#2F6BFF"
  ) +
  
  labs(
    title =
      "B  Strongest nominal MPCDS-drug associations",
    
    subtitle =
      "All five nominal p < 0.001 signals indicate greater relative sensitivity",
    
    x =
      "Spearman rho",
    
    y =
      NULL
  ) +
  
  theme_manuscript(
    base_size = 10
  )


# ------------------------------------------------------------
# 36. Figure 7C - strongest sensitivity vs resistance signals
# ------------------------------------------------------------

prism_extremes <-
  bind_rows(
    prism_sensitivity_top |>
      dplyr::slice_head(
        n = 8
      ) |>
      dplyr::mutate(
        direction =
          "Greater sensitivity"
      ),
    
    prism_resistance_top |>
      dplyr::slice_head(
        n = 8
      ) |>
      dplyr::mutate(
        direction =
          "Greater resistance"
      )
  )


prism_extremes$label <-
  paste0(
    prism_extremes$drug_name,
    "  "
  )


prism_extremes$label <-
  factor(
    prism_extremes$label,
    levels =
      prism_extremes |>
      dplyr::arrange(
        spearman_rho
      ) |>
      dplyr::pull(
        label
      )
  )


plot_prism_extremes <-
  ggplot(
    prism_extremes,
    aes(
      x = spearman_rho,
      y = label,
      fill = direction
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.5
  ) +
  
  geom_col(
    width = 0.65
  ) +
  
  scale_fill_manual(
    values = prism_direction_colors
  ) +
  
  labs(
    title =
      "C  Strongest directional associations",
    
    subtitle =
      "Top eight sensitivity-side and resistance-side correlations",
    
    x =
      "Spearman rho",
    
    y =
      NULL,
    
    fill =
      "Higher MPCDS associated with"
  ) +
  
  theme_manuscript(
    base_size = 8.8
  ) +
  
  theme(
    legend.position =
      "top"
  )


# ------------------------------------------------------------
# 37. Figure 7D - multiple-testing summary
# ------------------------------------------------------------

prism_testing_summary <-
  data.frame(
    threshold =
      factor(
        c(
          "Tested",
          "p < 0.05",
          "p < 0.01",
          "p < 0.001",
          "FDR < 0.05"
        ),
        levels =
          c(
            "Tested",
            "p < 0.05",
            "p < 0.01",
            "p < 0.001",
            "FDR < 0.05"
          )
      ),
    
    n =
      c(
        prism_tested,
        prism_nominal_005,
        prism_nominal_001,
        prism_nominal_0001_n,
        prism_fdr_005
      )
  )


prism_testing_summary$category <-
  ifelse(
    prism_testing_summary$threshold ==
      "FDR < 0.05",
    "Adjusted significance",
    "Nominal/testing"
  )


prism_summary_colors <-
  c(
    "Nominal/testing" = "#2F6BFF",
    "Adjusted significance" = "#E34A4A"
  )


plot_prism_testing <-
  ggplot(
    prism_testing_summary,
    aes(
      x = threshold,
      y = n,
      fill = category
    )
  ) +
  
  geom_col(
    width = 0.65
  ) +
  
  geom_text(
    aes(
      label =
        format(
          n,
          big.mark = ",",
          scientific = FALSE
        )
    ),
    vjust = -0.5,
    size = 3.8
  ) +
  
  scale_fill_manual(
    values = prism_summary_colors
  ) +
  
  scale_y_continuous(
    expand =
      expansion(
        mult =
          c(
            0,
            0.10
          )
      )
  ) +
  
  labs(
    title =
      "D  Multiple-testing summary",
    
    subtitle =
      paste0(
        "Minimum BH FDR = ",
        sprintf(
          "%.3f",
          prism_min_fdr
        )
      ),
    
    x =
      NULL,
    
    y =
      "Number of compounds",
    
    fill =
      NULL
  ) +
  
  theme_manuscript(
    base_size = 10
  ) +
  
  theme(
    legend.position =
      "none",
    
    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1
      )
  )


# ------------------------------------------------------------
# 38. Assemble Figure 7
# ------------------------------------------------------------

figure_7 <-
  (
    plot_prism_global |
      plot_prism_top
  ) /
  (
    plot_prism_extremes |
      plot_prism_testing
  )


figure_7


ggsave(
  "05_figures/manuscript/main/Figure_07_PRISM_pharmacogenomic_associations.png",
  figure_7,
  width = 12,
  height = 10,
  dpi = 300
)


ggsave(
  "05_figures/manuscript/main/Figure_07_PRISM_pharmacogenomic_associations.pdf",
  figure_7,
  width = 12,
  height = 10
)


# ------------------------------------------------------------
# 39. Verify Figure 7
# ------------------------------------------------------------

prism_manuscript_figures <-
  c(
    "05_figures/manuscript/main/Figure_07_PRISM_pharmacogenomic_associations.png",
    "05_figures/manuscript/main/Figure_07_PRISM_pharmacogenomic_associations.pdf"
  )


cat(
  "\nPRISM manuscript-style figure created:\n"
)


print(
  file.exists(
    prism_manuscript_figures
  )
)


# ============================================================
# FIGURE 1
# Overall study workflow
# ============================================================


# ------------------------------------------------------------
# 40. Complete study-workflow data
#
# This figure summarizes the full scientific workflow rather
# than mirroring script numbers one-by-one.
# ------------------------------------------------------------

workflow_data <-
  data.frame(
    step =
      factor(
        paste(
          "Step",
          1:6
        ),
        levels =
          paste(
            "Step",
            1:6
          )
      ),
    
    x =
      1:6,
    
    title =
      c(
        "Data acquisition and\nPCD gene discovery",
        
        "PCD subtype discovery\nand characterization",
        
        "External GEO\nsubtype validation",
        
        "LUAD MPCDS construction\nand prognostic validation",
        
        "Therapeutic and\nimmunotherapy characterization",
        
        "Sensitivity and\nintegrative analyses"
      ),
    
    detail =
      c(
        paste0(
          "TCGA-LUAD + TCGA-LUSC\n",
          "RNA-seq • Clinical • Mutation\n",
          "CNV • DNA methylation\n",
          "296 curated PCD genes\n",
          "Tumour-vs-normal DE analysis\n",
          "47 concordant PCD candidates"
        ),
        
        paste0(
          "Consensus clustering (k = 2)\n",
          "CPI • Gap • PAC • Silhouette\n",
          "PCD mechanism scoring\n",
          "Immune/stromal infiltration\n",
          "Hallmark pathway enrichment\n",
          "Mutation • CNV • Methylation"
        ),
        
        paste0(
          "Fixed TCGA-centroid projection\n",
          "LUAD: GSE68465 • GSE72094\n",
          "GSE31210 • GSE50081\n",
          "LUSC: GSE73403 • GSE4573\n",
          "GSE157010 • GSE30219\n",
          "Mechanism-specific reproducibility"
        ),
        
        paste0(
          "47-gene ridge-Cox MPCDS\n",
          "Repeated nested cross-validation\n",
          "TCGA-LUAD internal validation\n",
          "4 independent LUAD GEO cohorts\n",
          "Time-dependent ROC\n",
          "Decision-curve analysis"
        ),
        
        paste0(
          "DepMap PRISM screening\n",
          "GDSC2 / oncoPredict\n",
          "CD274 expression\n",
          "Nonsynonymous mutation burden\n",
          "MANTIS MSI\n",
          "TIDE immunotherapy surrogates"
        ),
        
        paste0(
          "70/30 LASSO-Cox sensitivity\n",
          "xCell • quanTIseq • TIMER\n",
          "EPIC • CIBERSORT/LM22\n",
          "47-gene PCA\n",
          "Subtype–MPCDS–stage Sankey\n",
          "Final tables and figures"
        )
      ),
    
    stringsAsFactors =
      FALSE
  )


workflow_colors <-
  c(
    "Step 1" = "#4E79A7",
    "Step 2" = "#59A14F",
    "Step 3" = "#76B7B2",
    "Step 4" = "#E15759",
    "Step 5" = "#B07AA1",
    "Step 6" = "#F28E2B"
  )


# ------------------------------------------------------------
# 41. Figure 1: complete overall study workflow
# ------------------------------------------------------------

workflow_arrows <-
  data.frame(
    x =
      seq(
        1.43,
        5.43,
        by = 1
      ),
    
    xend =
      seq(
        1.57,
        5.57,
        by = 1
      ),
    
    y =
      rep(
        0.50,
        5
      ),
    
    yend =
      rep(
        0.50,
        5
      )
  )


figure_1 <-
  ggplot() +
  
  geom_segment(
    data =
      workflow_arrows,
    
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    
    arrow =
      arrow(
        length =
          unit(
            0.16,
            "inches"
          ),
        type =
          "closed"
      ),
    
    linewidth =
      0.85,
    
    colour =
      "#555555"
  ) +
  
  geom_rect(
    data =
      workflow_data,
    
    aes(
      xmin =
        x - 0.40,
      
      xmax =
        x + 0.40,
      
      ymin =
        0.08,
      
      ymax =
        0.92,
      
      fill =
        step
    ),
    
    colour =
      "white",
    
    linewidth =
      1.1,
    
    alpha =
      0.96
  ) +
  
  geom_text(
    data =
      workflow_data,
    
    aes(
      x =
        x,
      
      y =
        0.82,
      
      label =
        step
    ),
    
    colour =
      "white",
    
    fontface =
      "bold",
    
    size =
      4.0
  ) +
  
  geom_text(
    data =
      workflow_data,
    
    aes(
      x =
        x,
      
      y =
        0.66,
      
      label =
        title
    ),
    
    colour =
      "white",
    
    fontface =
      "bold",
    
    size =
      3.45,
    
    lineheight =
      0.94
  ) +
  
  geom_text(
    data =
      workflow_data,
    
    aes(
      x =
        x,
      
      y =
        0.36,
      
      label =
        detail
    ),
    
    colour =
      "white",
    
    size =
      2.62,
    
    lineheight =
      1.10
  ) +
  
  scale_fill_manual(
    values =
      workflow_colors
  ) +
  
  coord_cartesian(
    xlim =
      c(
        0.45,
        6.55
      ),
    
    ylim =
      c(
        0,
        1
      ),
    
    clip =
      "off"
  ) +
  
  labs(
    title =
      "Overall study workflow",
    
    subtitle =
      paste0(
        "Programmed cell death profiling, molecular subtyping, external validation, ",
        "prognostic modelling, therapeutic characterization and integrative analysis"
      )
  ) +
  
  theme_void(
    base_size =
      12
  ) +
  
  theme(
    legend.position =
      "none",
    
    plot.title =
      element_text(
        face =
          "bold",
        
        size =
          18,
        
        hjust =
          0.5,
        
        margin =
          margin(
            b = 6
          )
      ),
    
    plot.subtitle =
      element_text(
        size =
          10.8,
        
        hjust =
          0.5,
        
        margin =
          margin(
            b = 18
          )
      ),
    
    plot.margin =
      margin(
        20,
        25,
        20,
        25
      )
  )


figure_1


ggsave(
  "05_figures/manuscript/main/Figure_01_overall_study_workflow.png",
  figure_1,
  width = 19,
  height = 6.5,
  dpi = 300,
  bg = "white"
)


ggsave(
  "05_figures/manuscript/main/Figure_01_overall_study_workflow.pdf",
  figure_1,
  width = 19,
  height = 6.5,
  bg = "white"
)


workflow_files <-
  c(
    "05_figures/manuscript/main/Figure_01_overall_study_workflow.png",
    "05_figures/manuscript/main/Figure_01_overall_study_workflow.pdf"
  )


cat(
  "\nComplete workflow figure created:\n"
)


print(
  file.exists(
    workflow_files
  )
)


cat(
  "External GEO subtype cohorts shown: 8\n"
)


cat(
  "LUAD MPCDS external-validation cohorts shown: 4\n"
)


# ============================================================
# FIGURE 2
# Identification and characterization of expression-derived
# PCD molecular subtypes
# ============================================================


# ------------------------------------------------------------
# 43. Additional package for reading consensus PNG files
# ------------------------------------------------------------

if (
  !requireNamespace(
    "png",
    quietly = TRUE
  )
) {
  
  stop(
    paste0(
      "Package 'png' is required for Figure 2. ",
      "Install it with install.packages('png')."
    )
  )
}


# ------------------------------------------------------------
# 44. Candidate 47-gene PCD set
# ------------------------------------------------------------

pcd_47_genes <-
  c(
    "SLC39A8",
    "MAP1LC3C",
    "DAPK2",
    "TLR4",
    "BMX",
    "SEPTIN4",
    "CYBB",
    "DAPK1",
    "PRKCQ",
    "CAMK2A",
    "VIM",
    "KAT2B",
    "CASP5",
    "ACSL4",
    "UACA",
    "IL1B",
    
    "DSG3",
    "AKR1C2",
    "PKP1",
    "PSAT1",
    "H1-5",
    "CDKN2A",
    "AKR1C1",
    "H1-3",
    "SLC7A11",
    "DSG1",
    "EGLN3",
    "H1-4",
    "DSP",
    "NQO1",
    "CLSPN",
    "GCLC",
    "AKR1C3",
    "H1-1",
    "TF",
    "CP",
    "LMNB1",
    "PMAIP1",
    "DSG2",
    "GCLM",
    "E2F1",
    "H1-2",
    "TXNRD1",
    "ADGRG1",
    "CBS",
    "NOX4",
    "H1-0"
  )


stopifnot(
  length(
    pcd_47_genes
  ) == 47
)


# ------------------------------------------------------------
# 45. Load clustering assignments
# ------------------------------------------------------------

luad_clusters <-
  read.csv(
    "04_results/clustering/TCGA_LUAD_PCD_cluster_assignments.csv",
    check.names = FALSE
  )


lusc_clusters <-
  read.csv(
    "04_results/clustering/TCGA_LUSC_PCD_cluster_assignments.csv",
    check.names = FALSE
  )


# ------------------------------------------------------------
# 46. Load PCD expression matrices
# ------------------------------------------------------------

luad_pcd_tpm <-
  readRDS(
    "02_processed_data/expression/TCGA_LUAD_PCD_TPM.rds"
  )


lusc_pcd_tpm <-
  readRDS(
    "02_processed_data/expression/TCGA_LUSC_PCD_TPM.rds"
  )


# ------------------------------------------------------------
# 47. Confirm availability of all 47 genes
# ------------------------------------------------------------

luad_genes_available <-
  intersect(
    pcd_47_genes,
    rownames(
      luad_pcd_tpm
    )
  )


lusc_genes_available <-
  intersect(
    pcd_47_genes,
    rownames(
      lusc_pcd_tpm
    )
  )


cat(
  "\nFigure 2 gene availability:\n"
)


cat(
  "LUAD:",
  length(
    luad_genes_available
  ),
  "/ 47\n"
)


cat(
  "LUSC:",
  length(
    lusc_genes_available
  ),
  "/ 47\n"
)


if (
  length(
    luad_genes_available
  ) != 47
) {
  
  stop(
    paste(
      "LUAD missing candidate genes:",
      paste(
        setdiff(
          pcd_47_genes,
          luad_genes_available
        ),
        collapse = ", "
      )
    )
  )
}


if (
  length(
    lusc_genes_available
  ) != 47
) {
  
  stop(
    paste(
      "LUSC missing candidate genes:",
      paste(
        setdiff(
          pcd_47_genes,
          lusc_genes_available
        ),
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 48. Helper: prepare expression heatmap data
# ------------------------------------------------------------

prepare_pcd_heatmap <-
  function(
    tpm_matrix,
    cluster_assignments,
    gene_set
  ) {
    
    sample_barcodes <-
      colnames(
        tpm_matrix
      )
    
    
    patient_ids <-
      substr(
        sample_barcodes,
        1,
        12
      )
    
    
    sample_type <-
      substr(
        sample_barcodes,
        14,
        15
      )
    
    
    sample_info <-
      data.frame(
        sample_barcode =
          sample_barcodes,
        
        patient_id =
          patient_ids,
        
        sample_type =
          sample_type,
        
        stringsAsFactors =
          FALSE
      )
    
    
    # Primary tumour samples only
    sample_info <-
      sample_info[
        sample_info$sample_type ==
          "01",
        ,
        drop = FALSE
      ]
    
    
    # Retain only patients with finalized cluster assignments
    sample_info <-
      merge(
        sample_info,
        cluster_assignments,
        by =
          "patient_id",
        all = FALSE,
        sort = FALSE
      )
    
    
    # Ensure one tumour sample per patient for visualization
    sample_info <-
      sample_info[
        !duplicated(
          sample_info$patient_id
        ),
        ,
        drop = FALSE
      ]
    
    
    expression <-
      tpm_matrix[
        gene_set,
        sample_info$sample_barcode,
        drop = FALSE
      ]
    
    
    # log2 TPM+1
    expression_log <-
      log2(
        expression + 1
      )
    
    
    # Gene-wise standardization
    expression_z <-
      t(
        scale(
          t(
            expression_log
          )
        )
      )
    
    
    expression_z[
      is.na(
        expression_z
      )
    ] <-
      0
    
    
    # Cap extremes for manuscript visualization
    expression_z[
      expression_z > 2
    ] <-
      2
    
    
    expression_z[
      expression_z < -2
    ] <-
      -2
    
    
    # ----------------------------------------------------
    # Order patients by cluster, then by within-cluster
    # expression similarity
    # ----------------------------------------------------
    
    ordered_samples <-
      character(
        0
      )
    
    
    for (
      cluster_name in
      c(
        "PCD_C1",
        "PCD_C2"
      )
    ) {
      
      cluster_samples <-
        sample_info$sample_barcode[
          sample_info$PCD_cluster ==
            cluster_name
        ]
      
      
      if (
        length(
          cluster_samples
        ) > 2
      ) {
        
        hc <-
          hclust(
            dist(
              t(
                expression_z[
                  ,
                  cluster_samples,
                  drop = FALSE
                ]
              )
            ),
            method =
              "ward.D2"
          )
        
        
        cluster_samples <-
          cluster_samples[
            hc$order
          ]
      }
      
      
      ordered_samples <-
        c(
          ordered_samples,
          cluster_samples
        )
    }
    
    
    expression_z <-
      expression_z[
        ,
        ordered_samples,
        drop = FALSE
      ]
    
    
    sample_info <-
      sample_info[
        match(
          ordered_samples,
          sample_info$sample_barcode
        ),
        ,
        drop = FALSE
      ]
    
    
    # ----------------------------------------------------
    # Long plotting table
    # ----------------------------------------------------
    
    heatmap_long <-
      as.data.frame(
        as.table(
          expression_z
        ),
        stringsAsFactors = FALSE
      )
    
    
    colnames(
      heatmap_long
    ) <-
      c(
        "gene",
        "sample_barcode",
        "z_score"
      )
    
    
    heatmap_long$PCD_cluster <-
      sample_info$PCD_cluster[
        match(
          heatmap_long$sample_barcode,
          sample_info$sample_barcode
        )
      ]
    
    
    heatmap_long$sample_order <-
      match(
        heatmap_long$sample_barcode,
        ordered_samples
      )
    
    
    heatmap_long$sample_order <-
      factor(
        heatmap_long$sample_order,
        levels =
          seq_along(
            ordered_samples
          )
      )
    
    
    # Preserve same gene order in both cohorts
    heatmap_long$gene <-
      factor(
        heatmap_long$gene,
        levels =
          rev(
            gene_set
          )
      )
    
    
    heatmap_long$PCD_cluster <-
      factor(
        heatmap_long$PCD_cluster,
        levels =
          c(
            "PCD_C1",
            "PCD_C2"
          )
      )
    
    
    list(
      data =
        heatmap_long,
      
      sample_info =
        sample_info,
      
      expression =
        expression_z
    )
  }


# ------------------------------------------------------------
# 49. Prepare LUAD and LUSC heatmap data
# ------------------------------------------------------------

luad_heatmap_object <-
  prepare_pcd_heatmap(
    tpm_matrix =
      luad_pcd_tpm,
    
    cluster_assignments =
      luad_clusters,
    
    gene_set =
      pcd_47_genes
  )


lusc_heatmap_object <-
  prepare_pcd_heatmap(
    tpm_matrix =
      lusc_pcd_tpm,
    
    cluster_assignments =
      lusc_clusters,
    
    gene_set =
      pcd_47_genes
  )


luad_heatmap_data <-
  luad_heatmap_object$data


lusc_heatmap_data <-
  lusc_heatmap_object$data


cat(
  "\nFigure 2 tumour samples represented:\n"
)


cat(
  "LUAD:",
  nrow(
    luad_heatmap_object$sample_info
  ),
  "\n"
)


cat(
  "LUSC:",
  nrow(
    lusc_heatmap_object$sample_info
  ),
  "\n"
)


cat(
  "\nLUAD cluster distribution:\n"
)


print(
  table(
    luad_heatmap_object$sample_info$PCD_cluster
  )
)


cat(
  "\nLUSC cluster distribution:\n"
)


print(
  table(
    lusc_heatmap_object$sample_info$PCD_cluster
  )
)


# ------------------------------------------------------------
# 50. Figure 2C - LUAD expression heatmap
# ------------------------------------------------------------

plot_luad_heatmap <-
  ggplot(
    luad_heatmap_data,
    aes(
      x =
        sample_order,
      
      y =
        gene,
      
      fill =
        z_score
    )
  ) +
  
  geom_raster() +
  
  facet_grid(
    ~ PCD_cluster,
    scales =
      "free_x",
    
    space =
      "free_x"
  ) +
  
  scale_fill_gradient2(
    low =
      "#2F6BFF",
    
    mid =
      "white",
    
    high =
      "#E34A4A",
    
    midpoint =
      0,
    
    limits =
      c(
        -2,
        2
      ),
    
    name =
      "Expression\nZ-score"
  ) +
  
  labs(
    title =
      "C  TCGA-LUAD PCD expression pattern",
    
    subtitle =
      "47 candidate PCD genes; primary tumour samples",
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  theme_minimal(
    base_size =
      7
  ) +
  
  theme(
    plot.title =
      element_text(
        face =
          "bold",
        
        size =
          11
      ),
    
    plot.subtitle =
      element_text(
        size =
          8
      ),
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank(),
    
    axis.text.y =
      element_text(
        size =
          5.5,
        colour =
          "black"
      ),
    
    strip.text =
      element_text(
        face =
          "bold",
        size =
          8
      ),
    
    strip.background =
      element_rect(
        fill =
          "grey95",
        colour =
          NA
      ),
    
    legend.position =
      "right"
  )


# ------------------------------------------------------------
# 51. Figure 2D - LUSC expression heatmap
# ------------------------------------------------------------

plot_lusc_heatmap <-
  ggplot(
    lusc_heatmap_data,
    aes(
      x =
        sample_order,
      
      y =
        gene,
      
      fill =
        z_score
    )
  ) +
  
  geom_raster() +
  
  facet_grid(
    ~ PCD_cluster,
    scales =
      "free_x",
    
    space =
      "free_x"
  ) +
  
  scale_fill_gradient2(
    low =
      "#2F6BFF",
    
    mid =
      "white",
    
    high =
      "#E34A4A",
    
    midpoint =
      0,
    
    limits =
      c(
        -2,
        2
      ),
    
    name =
      "Expression\nZ-score"
  ) +
  
  labs(
    title =
      "D  TCGA-LUSC PCD expression pattern",
    
    subtitle =
      "47 candidate PCD genes; primary tumour samples",
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  theme_minimal(
    base_size =
      7
  ) +
  
  theme(
    plot.title =
      element_text(
        face =
          "bold",
        
        size =
          11
      ),
    
    plot.subtitle =
      element_text(
        size =
          8
      ),
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank(),
    
    axis.text.y =
      element_text(
        size =
          5.5,
        colour =
          "black"
      ),
    
    strip.text =
      element_text(
        face =
          "bold",
        size =
          8
      ),
    
    strip.background =
      element_rect(
        fill =
          "grey95",
        colour =
          NA
      ),
    
    legend.position =
      "right"
  )


# ------------------------------------------------------------
# 52. Figure 2E - PCD mechanism differences
# ------------------------------------------------------------

luad_mechanism_tests <-
  read.csv(
    "04_results/clustering/TCGA_LUAD_PCD_mechanism_cluster_tests.csv",
    check.names = FALSE
  )


lusc_mechanism_tests <-
  read.csv(
    "04_results/clustering/TCGA_LUSC_PCD_mechanism_cluster_tests.csv",
    check.names = FALSE
  )


mechanism_plot_data <-
  bind_rows(
    luad_mechanism_tests |>
      mutate(
        histology =
          "LUAD"
      ),
    
    lusc_mechanism_tests |>
      mutate(
        histology =
          "LUSC"
      )
  )


mechanism_plot_data$enriched_cluster <-
  ifelse(
    mechanism_plot_data$C2_minus_C1 > 0,
    "PCD_C2",
    "PCD_C1"
  )


mechanism_plot_data$significance <-
  ifelse(
    mechanism_plot_data$padj < 0.05,
    "FDR < 0.05",
    "Not significant"
  )


mechanism_plot_data$PCD_type <-
  factor(
    mechanism_plot_data$PCD_type,
    levels =
      c(
        "Necroptosis",
        "Pyroptosis",
        "Apoptosis",
        "Cuproptosis",
        "Ferroptosis"
      )
  )


mechanism_plot_data$histology <-
  factor(
    mechanism_plot_data$histology,
    levels =
      c(
        "LUAD",
        "LUSC"
      )
  )


plot_mechanisms <-
  ggplot(
    mechanism_plot_data,
    aes(
      x =
        C2_minus_C1,
      
      y =
        PCD_type,
      
      colour =
        enriched_cluster,
      
      shape =
        significance
    )
  ) +
  
  geom_vline(
    xintercept =
      0,
    
    linetype =
      2,
    
    linewidth =
      0.5
  ) +
  
  geom_segment(
    aes(
      x =
        0,
      
      xend =
        C2_minus_C1,
      
      yend =
        PCD_type
    ),
    
    linewidth =
      0.8,
    
    alpha =
      0.65
  ) +
  
  geom_point(
    size =
      3.4,
    
    stroke =
      1
  ) +
  
  facet_wrap(
    ~ histology,
    nrow =
      1
  ) +
  
  scale_colour_manual(
    values =
      cluster_colors
  ) +
  
  scale_shape_manual(
    values =
      c(
        "FDR < 0.05" = 16,
        "Not significant" = 1
      )
  ) +
  
  labs(
    title =
      "E  PCD-mechanism differences between molecular subtypes",
    
    subtitle =
      "Difference = PCD_C2 - PCD_C1; filled symbols indicate FDR < 0.05",
    
    x =
      "Difference in mechanism score",
    
    y =
      NULL,
    
    colour =
      "Higher in",
    
    shape =
      "Statistical evidence"
  ) +
  
  theme_manuscript(
    base_size =
      9
  )


# ------------------------------------------------------------
# 53. Read and format consensus k = 2 images
# ------------------------------------------------------------

luad_consensus_image <-
  png::readPNG(
    "04_results/clustering/LUAD_consensus/consensus002.png"
  )


lusc_consensus_image <-
  png::readPNG(
    "04_results/clustering/LUSC_consensus/consensus002.png"
  )


luad_consensus_grob <-
  grid::rasterGrob(
    luad_consensus_image,
    interpolate = TRUE
  )


lusc_consensus_grob <-
  grid::rasterGrob(
    lusc_consensus_image,
    interpolate = TRUE
  )


plot_luad_consensus <-
  ggplot() +
  
  annotation_custom(
    grob = luad_consensus_grob,
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +
  
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  
  labs(
    title =
      "A  TCGA-LUAD consensus clustering",
    
    subtitle =
      "Consensus matrix at k = 2"
  ) +
  
  theme_void() +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11,
        hjust = 0
      ),
    
    plot.subtitle =
      element_text(
        size = 8.5,
        hjust = 0,
        margin = margin(
          b = 5
        )
      ),
    
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
  )


plot_lusc_consensus <-
  ggplot() +
  
  annotation_custom(
    grob = lusc_consensus_grob,
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +
  
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  
  labs(
    title =
      "B  TCGA-LUSC consensus clustering",
    
    subtitle =
      "Consensus matrix at k = 2"
  ) +
  
  theme_void() +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11,
        hjust = 0
      ),
    
    plot.subtitle =
      element_text(
        size = 8.5,
        hjust = 0,
        margin = margin(
          b = 5
        )
      ),
    
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
  )


# ------------------------------------------------------------
# 54. Assemble Figure 2
# ------------------------------------------------------------

figure_2 <-
  (
    plot_luad_consensus |
      plot_lusc_consensus
  ) /
  (
    plot_luad_heatmap |
      plot_lusc_heatmap
  ) /
  plot_mechanisms +
  
  plot_layout(
    heights =
      c(
        0.95,
        1.25,
        0.65
      )
  ) +
  
  plot_annotation(
    title =
      "Identification and characterization of expression-derived PCD molecular subtypes",
    
    subtitle =
      paste0(
        "Consensus clustering of the 47-gene candidate PCD expression landscape ",
        "in TCGA-LUAD and TCGA-LUSC"
      )
  ) &
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 15
      ),
    
    plot.subtitle =
      element_text(
        size = 10
      )
  )


figure_2


ggsave(
  "05_figures/manuscript/main/Figure_02_PCD_subtype_identification.png",
  figure_2,
  width = 13,
  height = 15,
  dpi = 300,
  bg = "white"
)


ggsave(
  "05_figures/manuscript/main/Figure_02_PCD_subtype_identification.pdf",
  figure_2,
  width = 13,
  height = 15,
  bg = "white"
)


# ------------------------------------------------------------
# 56. Verify Figure 2
# ------------------------------------------------------------

subtype_manuscript_files <-
  c(
    "05_figures/manuscript/main/Figure_02_PCD_subtype_identification.png",
    "05_figures/manuscript/main/Figure_02_PCD_subtype_identification.pdf"
  )


cat(
  "\nFigure 2 manuscript-style subtype figure created:\n"
)


print(
  file.exists(
    subtype_manuscript_files
  )
)


# ============================================================
# FINAL SCRIPT 32 VERIFICATION
# ============================================================

manuscript_final_files <-
  c(
    "05_figures/manuscript/main/Figure_01_overall_study_workflow.png",
    "05_figures/manuscript/main/Figure_01_overall_study_workflow.pdf",
    
    "05_figures/manuscript/main/Figure_02_PCD_subtype_identification.png",
    "05_figures/manuscript/main/Figure_02_PCD_subtype_identification.pdf",
    
    "05_figures/manuscript/main/Figure_03_immune_and_functional_characterization.png",
    "05_figures/manuscript/main/Figure_03_immune_and_functional_characterization.pdf",
    
    "05_figures/manuscript/main/Figure_04_multiomics_characterization.png",
    "05_figures/manuscript/main/Figure_04_multiomics_characterization.pdf",
    
    "05_figures/manuscript/main/Figure_05_MPCDS_construction_internal_validation.png",
    "05_figures/manuscript/main/Figure_05_MPCDS_construction_internal_validation.pdf",
    
    "05_figures/manuscript/main/Figure_06_external_MPCDS_validation.png",
    "05_figures/manuscript/main/Figure_06_external_MPCDS_validation.pdf",
    
    "05_figures/manuscript/main/Figure_07_PRISM_pharmacogenomic_associations.png",
    "05_figures/manuscript/main/Figure_07_PRISM_pharmacogenomic_associations.pdf"
  )


manuscript_verified <-
  file.exists(
    manuscript_final_files
  )


cat(
  "\n",
  paste(
    rep(
      "=",
      70
    ),
    collapse = ""
  ),
  "\n"
)


cat(
  "SCRIPT 32 COMPLETE\n"
)


cat(
  paste(
    rep(
      "=",
      70
    ),
    collapse = ""
  ),
  "\n"
)


cat(
  "Manuscript figures:",
  7,
  "\n"
)


cat(
  "Expected figure files:",
  length(
    manuscript_final_files
  ),
  "\n"
)


cat(
  "LUAD Figure 2 genes:",
  length(
    luad_genes_available
  ),
  "/ 47\n"
)


cat(
  "LUSC Figure 2 genes:",
  length(
    lusc_genes_available
  ),
  "/ 47\n"
)


cat(
  "LUAD clustered tumours:",
  nrow(
    luad_heatmap_object$sample_info
  ),
  "\n"
)


cat(
  "LUSC clustered tumours:",
  nrow(
    lusc_heatmap_object$sample_info
  ),
  "\n"
)


cat(
  "PRISM compounds tested:",
  as.integer(
    prism_tested
  ),
  "\n"
)


cat(
  "PRISM FDR < 0.05:",
  as.integer(
    prism_fdr_005
  ),
  "\n"
)


cat(
  "All manuscript figure files verified:",
  all(
    manuscript_verified
  ),
  "\n"
)


cat(
  paste(
    rep(
      "=",
      70
    ),
    collapse = ""
  ),
  "\n"
)



# ============================================================
# CURRENT PIPELINE SUPPLEMENTARY FIGURES
# ============================================================

# No model is refitted below. These plots use locked results from
# Scripts 14, 20, 28, 29, 30 and 31.


# ------------------------------------------------------------
# A. Cluster-number diagnostic summary
# ------------------------------------------------------------

cluster_diag_s32 <-
  data.frame(
    histology = rep(c("LUAD", "LUSC"), each = 5),
    criterion = rep(
      c(
        "PAC minimum",
        "CPI peak",
        "Gap global maximum",
        "Gap 1-SE",
        "Silhouette maximum"
      ),
      2
    ),
    selected_k = c(
      8, 2, 7, 3, 2,
      8, 2, 3, 2, 2
    ),
    stringsAsFactors = FALSE
  )


silhouette_s32 <-
  rbind(
    data.frame(
      histology = "LUAD",
      k = 2:8,
      mean_silhouette = c(
        0.17745570, 0.16381477, 0.12700184,
        0.14104407, 0.11683610, 0.11188014,
        0.09340932
      )
    ),
    data.frame(
      histology = "LUSC",
      k = 2:8,
      mean_silhouette = c(
        0.23453624, 0.16458417, 0.13414381,
        0.10908706, 0.11361394, 0.08278907,
        0.07297321
      )
    )
  )


plot_cluster_choice_s32 <-
  ggplot(
    cluster_diag_s32,
    aes(x = selected_k, y = criterion)
  ) +
  geom_vline(
    xintercept = 2,
    linetype = 2,
    colour = "grey45"
  ) +
  geom_point(size = 3) +
  facet_wrap(~ histology, ncol = 1) +
  scale_x_continuous(
    breaks = 2:8,
    limits = c(1.7, 8.3)
  ) +
  labs(
    title = "A. Cluster-number criteria",
    subtitle =
      paste0(
        "k=2 was retained as the parsimonious solution; ",
        "not all diagnostics selected k=2."
      ),
    x = "Selected cluster number (k)",
    y = NULL
  ) +
  theme_manuscript(base_size = 10)


plot_silhouette_s32 <-
  ggplot(
    silhouette_s32,
    aes(
      x = k,
      y = mean_silhouette,
      colour = histology,
      group = histology
    )
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_wrap(~ histology, ncol = 1) +
  scale_x_continuous(breaks = 2:8) +
  scale_colour_manual(values = histology_colors) +
  labs(
    title = "B. Mean silhouette width",
    x = "Cluster number (k)",
    y = "Mean silhouette width"
  ) +
  theme_manuscript(base_size = 10) +
  theme(legend.position = "none")


figure_cluster_diagnostics_s32 <-
  plot_cluster_choice_s32 |
  plot_silhouette_s32


# ------------------------------------------------------------
# B. Eight-cohort external subtype overview
# ------------------------------------------------------------

external_subtype_file_s32 <-
  "06_tables/final/Table_02_subtype_validation_summary.csv"

stopifnot(file.exists(external_subtype_file_s32))

external_subtype_s32 <-
  read.csv(
    external_subtype_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

external_subtype_s32 <-
  external_subtype_s32[
    grepl("^GSE", external_subtype_s32$cohort),
    ,
    drop = FALSE
  ]

external_subtype_s32$cohort <-
  factor(
    external_subtype_s32$cohort,
    levels = rev(external_subtype_s32$cohort)
  )

plot_external_silhouette_s32 <-
  ggplot(
    external_subtype_s32,
    aes(
      x = mean_silhouette,
      y = cohort,
      shape = histology
    )
  ) +
  geom_point(size = 3) +
  labs(
    title = "A. External cohort silhouette",
    x = "Mean silhouette width",
    y = NULL
  ) +
  theme_manuscript(base_size = 10)

external_counts_s32 <-
  rbind(
    data.frame(
      cohort = external_subtype_s32$cohort,
      subtype = "PCD_C1",
      n = external_subtype_s32$PCD_C1
    ),
    data.frame(
      cohort = external_subtype_s32$cohort,
      subtype = "PCD_C2",
      n = external_subtype_s32$PCD_C2
    )
  )

plot_external_counts_s32 <-
  ggplot(
    external_counts_s32,
    aes(
      x = n,
      y = cohort,
      fill = subtype
    )
  ) +
  geom_col() +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "B. Projected subtype composition",
    x = "Patients",
    y = NULL,
    fill = "Projected subtype"
  ) +
  theme_manuscript(base_size = 10)

figure_external_subtypes_s32 <-
  plot_external_silhouette_s32 |
  plot_external_counts_s32


# ------------------------------------------------------------
# C. LASSO 70/30 sensitivity
# ------------------------------------------------------------

lasso_summary_file_s32 <-
  paste0(
    "04_results/mpcds/lasso_70_30_sensitivity/",
    "TCGA_LUAD_LASSO_70_30_final_summary.csv"
  )

stopifnot(file.exists(lasso_summary_file_s32))

lasso_summary_s32 <-
  read.csv(
    lasso_summary_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

get_lasso_metric_s32 <-
  function(metric_s32) {
    as.numeric(
      lasso_summary_s32$value[
        lasso_summary_s32$metric == metric_s32
      ][1]
    )
  }

lasso_perf_s32 <-
  data.frame(
    measure = c(
      "Training C-index",
      "Held-out C-index",
      "1-year AUC",
      "3-year AUC",
      "5-year AUC"
    ),
    value = c(
      get_lasso_metric_s32("Apparent training C-index"),
      get_lasso_metric_s32("Held-out test C-index"),
      get_lasso_metric_s32("Held-out 1-year AUC"),
      get_lasso_metric_s32("Held-out 3-year AUC"),
      get_lasso_metric_s32("Held-out 5-year AUC")
    ),
    stringsAsFactors = FALSE
  )

lasso_perf_s32$measure <-
  factor(
    lasso_perf_s32$measure,
    levels = rev(lasso_perf_s32$measure)
  )

figure_lasso_s32 <-
  ggplot(
    lasso_perf_s32,
    aes(x = value, y = measure)
  ) +
  geom_vline(
    xintercept = 0.5,
    linetype = 2,
    colour = "grey50"
  ) +
  geom_point(size = 3) +
  coord_cartesian(xlim = c(0.40, 0.72)) +
  labs(
    title = "TCGA-LUAD 70/30 LASSO-Cox sensitivity analysis",
    subtitle =
      "The sparse four-gene model showed weak held-out generalization.",
    x = "Discrimination",
    y = NULL
  ) +
  theme_manuscript(base_size = 11)


# ------------------------------------------------------------
# D. Multi-method immune replication
# ------------------------------------------------------------

immune_rep_file_s32 <-
  paste0(
    "04_results/immune/additional_deconvolution_sensitivity/",
    "TCGA_LUAD_crossmethod_final_immune_summary.csv"
  )

stopifnot(file.exists(immune_rep_file_s32))

immune_rep_s32 <-
  read.csv(
    immune_rep_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

immune_rep_plot_s32 <-
  immune_rep_s32[
    immune_rep_s32$subtype_replicated_FDR05 |
      immune_rep_s32$MPCDS_replicated_FDR05 |
      immune_rep_s32$subtype_significant_conflict |
      immune_rep_s32$MPCDS_significant_conflict,
    ,
    drop = FALSE
  ]

immune_rep_long_s32 <-
  rbind(
    data.frame(
      biological_population =
        immune_rep_plot_s32$biological_population,
      analysis = "PCD subtype",
      n_significant_methods =
        immune_rep_plot_s32$subtype_n_FDR05,
      replicated =
        immune_rep_plot_s32$subtype_replicated_FDR05,
      conflict =
        immune_rep_plot_s32$subtype_significant_conflict
    ),
    data.frame(
      biological_population =
        immune_rep_plot_s32$biological_population,
      analysis = "MPCDS",
      n_significant_methods =
        immune_rep_plot_s32$MPCDS_n_FDR05,
      replicated =
        immune_rep_plot_s32$MPCDS_replicated_FDR05,
      conflict =
        immune_rep_plot_s32$MPCDS_significant_conflict
    )
  )

immune_rep_long_s32$status <-
  ifelse(
    immune_rep_long_s32$conflict,
    "Conflicting",
    ifelse(
      immune_rep_long_s32$replicated,
      "Replicated",
      "Other"
    )
  )

figure_immune_replication_s32 <-
  ggplot(
    immune_rep_long_s32,
    aes(
      x = n_significant_methods,
      y = biological_population,
      shape = status
    )
  ) +
  geom_point(size = 3) +
  facet_wrap(~ analysis, ncol = 2) +
  scale_x_continuous(breaks = 0:5) +
  labs(
    title = "Cross-method immune-deconvolution sensitivity",
    x = "Methods significant at FDR < 0.05",
    y = NULL,
    shape = "Status"
  ) +
  theme_manuscript(base_size = 10)


# ------------------------------------------------------------
# E. 47-gene PCA
# ------------------------------------------------------------

pca_scores_file_s32 <-
  "04_results/final_summary/TCGA_LUAD_47gene_PCD_subtype_PCA_scores.csv"

stopifnot(file.exists(pca_scores_file_s32))

pca_scores_current_s32 <-
  read.csv(
    pca_scores_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

stopifnot(
  all(
    c("PC1", "PC2", "PCD_cluster") %in%
      colnames(pca_scores_current_s32)
  )
)

figure_pca_s32 <-
  ggplot(
    pca_scores_current_s32,
    aes(
      x = PC1,
      y = PC2,
      colour = PCD_cluster
    )
  ) +
  geom_point(
    alpha = 0.75,
    size = 1.8
  ) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    title = "TCGA-LUAD 47-gene PCD expression space",
    subtitle = "PC1 = 15.70%; PC2 = 11.77%",
    colour = "PCD subtype"
  ) +
  theme_manuscript(base_size = 11)


# ------------------------------------------------------------
# F. 496-patient Sankey/alluvial
# ------------------------------------------------------------

sankey_file_s32 <-
  "04_results/final_summary/TCGA_LUAD_Sankey_complete_case_patients.csv"

stopifnot(file.exists(sankey_file_s32))

sankey_current_s32 <-
  read.csv(
    sankey_file_s32,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

stopifnot(nrow(sankey_current_s32) == 496L)

sankey_current_s32$PCD_cluster <-
  factor(
    sankey_current_s32$PCD_cluster,
    levels = c("PCD_C1", "PCD_C2")
  )

sankey_current_s32$risk_group <-
  factor(
    sankey_current_s32$risk_group,
    levels = c("Low MPCDS", "High MPCDS")
  )

sankey_current_s32$stage_group <-
  factor(
    sankey_current_s32$stage_group,
    levels = c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  )

figure_sankey_s32 <-
  ggplot(
    sankey_current_s32,
    aes(
      axis1 = PCD_cluster,
      axis2 = risk_group,
      axis3 = stage_group,
      y = 1
    )
  ) +
  ggalluvial::geom_alluvium(
    aes(fill = PCD_cluster),
    width = 0.12,
    alpha = 0.65
  ) +
  ggalluvial::geom_stratum(
    width = 0.12,
    fill = "grey95",
    colour = "grey45"
  ) +
  ggalluvial::stat_stratum(
    geom = "text",
    aes(label = after_stat(stratum)),
    size = 3.2
  ) +
  scale_x_discrete(
    limits = c(
      "PCD subtype",
      "MPCDS group",
      "AJCC stage"
    ),
    expand = c(0.08, 0.05)
  ) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "TCGA-LUAD PCD subtype–MPCDS–stage integration",
    subtitle = "Complete-case cohort: n = 496",
    x = NULL,
    y = "Patients",
    fill = "PCD subtype"
  ) +
  theme_manuscript(base_size = 11)


# ------------------------------------------------------------
# Save helpers and final files
# ------------------------------------------------------------

save_all_formats_s32 <-
  function(
    plot_s32,
    base_path_s32,
    width_s32,
    height_s32
  ) {
    
    ggplot2::ggsave(
      paste0(base_path_s32, ".png"),
      plot_s32,
      width = width_s32,
      height = height_s32,
      dpi = 300,
      bg = "white"
    )
    
    ggplot2::ggsave(
      paste0(base_path_s32, ".pdf"),
      plot_s32,
      width = width_s32,
      height = height_s32,
      bg = "white"
    )
    
    ggplot2::ggsave(
      paste0(base_path_s32, ".tiff"),
      plot_s32,
      width = width_s32,
      height = height_s32,
      dpi = 300,
      compression = "lzw",
      bg = "white"
    )
  }


main_plot_list_s32 <-
  list(
    Figure_01_overall_study_workflow = figure_1,
    Figure_02_PCD_subtype_identification = figure_2,
    Figure_03_immune_and_functional_characterization = figure_3,
    Figure_04_multiomics_characterization = figure_4,
    Figure_05_MPCDS_construction_internal_validation = figure_5,
    Figure_06_external_MPCDS_validation = figure_6,
    Figure_07_PRISM_pharmacogenomic_associations = figure_7
  )

main_dims_s32 <-
  data.frame(
    figure = names(main_plot_list_s32),
    width = c(19, 14, 13, 13, 13, 13, 12),
    height = c(6.5, 11, 9, 10, 10, 14, 10),
    stringsAsFactors = FALSE
  )

for (i_s32 in seq_len(nrow(main_dims_s32))) {
  fig_s32 <- main_dims_s32$figure[i_s32]
  
  # PNG and PDF were already saved by the original manuscript
  # figure sections. Add required TIFF output here.
  ggplot2::ggsave(
    filename =
      file.path(
        manuscript_figure_dir,
        paste0(fig_s32, ".tiff")
      ),
    plot = main_plot_list_s32[[fig_s32]],
    width = main_dims_s32$width[i_s32],
    height = main_dims_s32$height[i_s32],
    dpi = 300,
    compression = "lzw",
    bg = "white"
  )
}


# ------------------------------------------------------------
# Save the 23 individual manuscript panels
# ------------------------------------------------------------

individual_plot_specs_s32 <-
  list(
    list(
      plot = plot_luad_consensus,
      relative = "clustering/Figure_02A_LUAD_consensus_k2",
      width = 6,
      height = 5.5
    ),
    list(
      plot = plot_lusc_consensus,
      relative = "clustering/Figure_02B_LUSC_consensus_k2",
      width = 6,
      height = 5.5
    ),
    list(
      plot = plot_luad_heatmap,
      relative = "clustering/Figure_02C_LUAD_47gene_expression_heatmap",
      width = 8,
      height = 8
    ),
    list(
      plot = plot_lusc_heatmap,
      relative = "clustering/Figure_02D_LUSC_47gene_expression_heatmap",
      width = 8,
      height = 8
    ),
    list(
      plot = plot_mechanisms,
      relative = "clustering/Figure_02E_PCD_mechanism_differences",
      width = 10,
      height = 5
    ),
    list(
      plot = plot_immune,
      relative = "immune/Figure_03A_immune_stromal_infiltration",
      width = 8,
      height = 5.5
    ),
    list(
      plot = plot_hallmark,
      relative = "enrichment/Figure_03B_Hallmark_pathway_enrichment",
      width = 9,
      height = 7
    ),
    list(
      plot = plot_mutation,
      relative = "mutation/Figure_04A_nonsynonymous_mutation_burden",
      width = 8,
      height = 5.5
    ),
    list(
      plot = plot_keap1,
      relative = "mutation/Figure_04B_KEAP1_mutation_frequency",
      width = 7,
      height = 5
    ),
    list(
      plot = plot_cnv,
      relative = "cnv/Figure_04C_CNV_burden",
      width = 7,
      height = 5
    ),
    list(
      plot = plot_methylation,
      relative = "methylation/Figure_04D_methylation_concordance",
      width = 7,
      height = 5
    ),
    list(
      plot = plot_cv,
      relative = "mpcds/Figure_05A_repeated_nested_CV",
      width = 6.5,
      height = 5
    ),
    list(
      plot = plot_coefficients,
      relative = "mpcds/Figure_05B_LUAD_ridge_coefficients",
      width = 7,
      height = 6
    ),
    list(
      plot = plot_tcga_km,
      relative = "mpcds/Figure_05C_TCGA_LUAD_MPCDS_KM",
      width = 7,
      height = 5.5
    ),
    list(
      plot = plot_incremental,
      relative = "mpcds/Figure_05D_internal_incremental_prediction",
      width = 7,
      height = 5.5
    ),
    list(
      plot = plot_external_km,
      relative = "mpcds/Figure_06A_GSE68465_MPCDS_KM",
      width = 7,
      height = 5.5
    ),
    list(
      plot = plot_forest,
      relative = "mpcds/Figure_06B_MPCDS_prognostic_forest",
      width = 8,
      height = 5.5
    ),
    list(
      plot = plot_external_cindex,
      relative = "mpcds/Figure_06C_external_incremental_Cindex",
      width = 6.5,
      height = 5
    ),
    list(
      plot = plot_time_hr,
      relative = "mpcds/Figure_06D_external_time_varying_HR",
      width = 6.5,
      height = 5
    ),
    list(
      plot = plot_prism_global,
      relative = "drug_sensitivity/Figure_07A_PRISM_global_landscape",
      width = 7,
      height = 5.5
    ),
    list(
      plot = plot_prism_top,
      relative = "drug_sensitivity/Figure_07B_PRISM_top_nominal_associations",
      width = 7,
      height = 5
    ),
    list(
      plot = plot_prism_extremes,
      relative = "drug_sensitivity/Figure_07C_PRISM_directional_extremes",
      width = 8,
      height = 6
    ),
    list(
      plot = plot_prism_testing,
      relative = "drug_sensitivity/Figure_07D_PRISM_multiple_testing_summary",
      width = 7,
      height = 5
    )
  )


stopifnot(
  length(
    individual_plot_specs_s32
  ) ==
    23L
)


for (
  spec_s32 in individual_plot_specs_s32
) {
  
  save_all_formats_s32(
    plot_s32 =
      spec_s32$plot,
    base_path_s32 =
      file.path(
        individual_figure_dir,
        spec_s32$relative
      ),
    width_s32 =
      spec_s32$width,
    height_s32 =
      spec_s32$height
  )
}


supplementary_plots_s32 <-
  list(
    Cluster_number_diagnostics =
      figure_cluster_diagnostics_s32,
    External_subtype_validation_overview =
      figure_external_subtypes_s32,
    LASSO_70_30_sensitivity =
      figure_lasso_s32,
    Immune_deconvolution_replication =
      figure_immune_replication_s32,
    LUAD_47gene_PCA =
      figure_pca_s32,
    LUAD_PCD_MPCDS_stage_Sankey =
      figure_sankey_s32
  )

supp_dims_s32 <-
  data.frame(
    figure = names(supplementary_plots_s32),
    width = c(11, 12, 9, 12, 8.5, 11),
    height = c(8, 7, 6, 8, 7, 7),
    stringsAsFactors = FALSE
  )

for (i_s32 in seq_len(nrow(supp_dims_s32))) {
  fig_s32 <- supp_dims_s32$figure[i_s32]
  
  save_all_formats_s32(
    plot_s32 = supplementary_plots_s32[[fig_s32]],
    base_path_s32 =
      file.path(
        supplementary_figure_dir,
        paste0("Supplementary_", fig_s32)
      ),
    width_s32 = supp_dims_s32$width[i_s32],
    height_s32 = supp_dims_s32$height[i_s32]
  )
}



# ============================================================
# LUAD MPCDS time-dependent ROC: TCGA plus four GEO cohorts
# Source: Script 22 ROC objects and summary; no models refitted.
# Export combined and standalone external validation figures.
# ============================================================

roc_dir_s32roc <- "04_results/mpcds/time_dependent_ROC"
object_file_s32roc <- file.path(roc_dir_s32roc, "LUAD_MPCDS_timeROC_objects.rds")
summary_file_s32roc <- file.path(
  roc_dir_s32roc, "LUAD_MPCDS_time_dependent_ROC_all_cohorts.csv"
)
if (!file.exists(object_file_s32roc) || !file.exists(summary_file_s32roc)) {
  stop("Script 22 ROC outputs are missing. Run Script 22 before Script 32.")
}

roc_list_s32roc <- readRDS(object_file_s32roc)
summary_s32roc <- read.csv(summary_file_s32roc, stringsAsFactors = FALSE,
                           check.names = FALSE)
cohorts_s32roc <- c("TCGA-LUAD", "GSE68465", "GSE72094", "GSE31210", "GSE50081")
keys_s32roc <- c("TCGA_LUAD", "GSE68465", "GSE72094", "GSE31210", "GSE50081")
horizons_s32roc <- c("1 year", "3 years", "5 years")
cols_s32roc <- c("#2874B2", "#D48116", "#30916B")
stopifnot(
  all(keys_s32roc %in% names(roc_list_s32roc)),
  all(c("cohort", "timepoint", "time_value", "n", "events",
        "AUC", "CI_lower", "CI_upper", "known_survivors_beyond_time") %in%
        names(summary_s32roc)),
  nrow(summary_s32roc) == 15L,
  setequal(summary_s32roc$cohort, cohorts_s32roc),
  all(is.finite(summary_s32roc$AUC)),
  !anyDuplicated(summary_s32roc[c("cohort", "timepoint")])
)

ordered_rows_s32roc <- do.call(rbind, lapply(cohorts_s32roc, function(cohort) {
  d <- summary_s32roc[summary_s32roc$cohort == cohort, , drop = FALSE]
  d <- d[match(horizons_s32roc, d$timepoint), , drop = FALSE]
  stopifnot(nrow(d) == 3L, !anyNA(d$timepoint))
  d
}))
rownames(ordered_rows_s32roc) <- NULL

# Confirm the ROC objects and printed AUC values come from the same run.
for (i in seq_along(cohorts_s32roc)) {
  d <- ordered_rows_s32roc[ordered_rows_s32roc$cohort == cohorts_s32roc[i], ]
  roc <- roc_list_s32roc[[keys_s32roc[i]]]
  stopifnot(is.matrix(roc$FP), is.matrix(roc$TP),
            ncol(roc$FP) == ncol(roc$TP))
  for (j in seq_len(3L)) {
    pos <- which(abs(roc$times - d$time_value[j]) < 1e-7)
    stopifnot(length(pos) == 1L,
              abs(roc$AUC[pos] - d$AUC[j]) < 1e-6)
  }
}

draw_cohort_s32roc <- function(i, standalone = FALSE) {
  d <- ordered_rows_s32roc[ordered_rows_s32roc$cohort == cohorts_s32roc[i], ]
  roc <- roc_list_s32roc[[keys_s32roc[i]]]
  graphics::plot(NA_real_, NA_real_, xlim = c(0, 1), ylim = c(0, 1),
                 asp = 1, xlab = "1 - specificity", ylab = "Sensitivity",
                 main = if (standalone) {
                   paste0(cohorts_s32roc[i], " LUAD MPCDS")
                 } else {
                   paste0(LETTERS[i], "  ", cohorts_s32roc[i])
                 },
                 cex.main = 1.07, cex.axis = 0.83, cex.lab = 0.93)
  graphics::abline(0, 1, lty = 3, col = "#9BA4AE")
  for (j in seq_len(3L)) {
    pos <- which(abs(roc$times - d$time_value[j]) < 1e-7)
    x <- as.numeric(roc$FP[, pos])
    y <- as.numeric(roc$TP[, pos])
    ok <- is.finite(x) & is.finite(y)
    stopifnot(any(ok))
    ord <- order(x[ok], y[ok])
    graphics::lines(c(0, x[ok][ord], 1), c(0, y[ok][ord], 1),
                    col = cols_s32roc[j], lwd = 2.4)
  }
  graphics::legend(
    "bottomright",
    legend = paste0(c("1 year", "3 years", "5 years"), ": AUC ",
                    sprintf("%.3f", d$AUC)),
    col = cols_s32roc, lwd = 2.4, cex = 0.77, bty = "n"
  )
  graphics::mtext(sprintf("n = %d; deaths = %d", d$n[1], d$events[1]),
                  side = 3, line = 0.2, cex = 0.76)
}

draw_summary_s32roc <- function() {
  d <- ordered_rows_s32roc[ordered_rows_s32roc$cohort != "TCGA-LUAD", ]
  cohort_levels <- cohorts_s32roc[-1]
  y <- 5 - match(d$cohort, cohort_levels) + rep(c(0.18, 0, -0.18), 4L)
  graphics::plot(NA_real_, NA_real_, xlim = c(0.35, 1.03),
                 ylim = c(0.5, 4.65), yaxt = "n", ylab = "",
                 xlab = "Time-dependent AUC (95% CI)",
                 main = "F  Four external LUAD cohorts",
                 cex.main = 1.03, cex.axis = 0.83, cex.lab = 0.91)
  graphics::axis(2, at = 4:1, labels = cohort_levels, las = 1, cex.axis = 0.79)
  graphics::abline(v = 0.5, lty = 3, col = "#8E99A5")
  graphics::segments(pmax(d$CI_lower, 0.35), y,
                     pmin(d$CI_upper, 1.03), y,
                     col = cols_s32roc[match(d$timepoint, horizons_s32roc)],
                     lwd = 1.4)
  graphics::points(d$AUC, y, pch = 19, cex = 0.8,
                   col = cols_s32roc[match(d$timepoint, horizons_s32roc)])
  graphics::legend("topright", legend = horizons_s32roc, col = cols_s32roc,
                   pch = 19, cex = 0.73, bty = "n", horiz = TRUE)
  graphics::mtext("Interpret sparse follow-up horizons cautiously",
                  side = 1, line = 3.2, cex = 0.71)
}

draw_figure_s32roc <- function() {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mfrow = c(2, 3), mar = c(4.3, 4.6, 3.9, 1.5),
                oma = c(0.9, 0.4, 2, 0.4))
  for (i in seq_along(cohorts_s32roc)) draw_cohort_s32roc(i)
  draw_summary_s32roc()
  graphics::mtext(
    "LUAD MPCDS: 1-, 3- and 5-year time-dependent ROC in TCGA and four GEO cohorts",
    side = 3, outer = TRUE, line = 0.25, font = 2, cex = 1.15
  )
}

out_dir_s32roc <- "05_figures/manuscript/supplementary"
dir.create(out_dir_s32roc, recursive = TRUE, showWarnings = FALSE)
stem_s32roc <- file.path(out_dir_s32roc, "Figure_S11_LUAD_MPCDS_time_dependent_ROC")

grDevices::png(paste0(stem_s32roc, ".png"), width = 4800, height = 3200,
               res = 300, bg = "white")
tryCatch(draw_figure_s32roc(), finally = grDevices::dev.off())

grDevices::pdf(paste0(stem_s32roc, ".pdf"), width = 16, height = 10.67,
               useDingbats = FALSE, bg = "white")
tryCatch(draw_figure_s32roc(), finally = grDevices::dev.off())

grDevices::tiff(paste0(stem_s32roc, ".tiff"), width = 4800, height = 3200,
                res = 300, compression = "lzw", bg = "white")
tryCatch(draw_figure_s32roc(), finally = grDevices::dev.off())

output_s32roc <- paste0(stem_s32roc, c(".png", ".pdf", ".tiff"))
stopifnot(all(file.exists(output_s32roc)),
          all(file.info(output_s32roc)$size > 0))

# Export each external cohort panel separately as a thesis-ready figure.
individual_dir_s32roc <- "05_figures/manuscript/individual/external_validation"
dir.create(individual_dir_s32roc, recursive = TRUE, showWarnings = FALSE)
draw_individual_s32roc <- function(i) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(5, 5, 4, 2))
  draw_cohort_s32roc(i, standalone = TRUE)
  graphics::mtext("IPCW ROC at 1, 3 and 5 years; continuous MPCDS",
                  side = 1, line = 3.4, cex = 0.82)
}
individual_files_s32roc <- character()
for (i in 2:5) {
  base <- file.path(individual_dir_s32roc,
                    paste0("Figure_EV_", cohorts_s32roc[i], "_MPCDS_timeROC"))
  grDevices::png(paste0(base, ".png"), width = 2400, height = 2100,
                 res = 300, bg = "white")
  tryCatch(draw_individual_s32roc(i), finally = grDevices::dev.off())
  grDevices::pdf(paste0(base, ".pdf"), width = 8, height = 7,
                 useDingbats = FALSE, bg = "white")
  tryCatch(draw_individual_s32roc(i), finally = grDevices::dev.off())
  grDevices::tiff(paste0(base, ".tiff"), width = 2400, height = 2100,
                  res = 300, compression = "lzw", bg = "white")
  tryCatch(draw_individual_s32roc(i), finally = grDevices::dev.off())
  individual_files_s32roc <- c(individual_files_s32roc,
                               paste0(base, c(".png", ".pdf", ".tiff")))
}
stopifnot(length(individual_files_s32roc) == 12L,
          all(file.exists(individual_files_s32roc)),
          all(file.info(individual_files_s32roc)$size > 0))


# ============================================================
# FINAL MANIFEST AND AUDIT
# ============================================================

main_expected_s32 <-
  unlist(
    lapply(
      names(main_plot_list_s32),
      function(name_s32) {
        file.path(
          manuscript_figure_dir,
          paste0(
            name_s32,
            c(".png", ".pdf", ".tiff")
          )
        )
      }
    ),
    use.names = FALSE
  )


individual_expected_s32 <-
  unlist(
    lapply(
      individual_plot_specs_s32,
      function(spec_s32) {
        paste0(
          file.path(
            individual_figure_dir,
            spec_s32$relative
          ),
          c(".png", ".pdf", ".tiff")
        )
      }
    ),
    use.names = FALSE
  )


supp_expected_s32 <-
  unlist(
    lapply(
      names(supplementary_plots_s32),
      function(name_s32) {
        file.path(
          supplementary_figure_dir,
          paste0(
            "Supplementary_",
            name_s32,
            c(".png", ".pdf", ".tiff")
          )
        )
      }
    ),
    use.names = FALSE
  )

individual_expected_s32 <- c(individual_expected_s32, individual_files_s32roc)
supp_expected_s32 <- c(supp_expected_s32, output_s32roc)

all_expected_figures_s32 <-
  c(
    main_expected_s32,
    individual_expected_s32,
    supp_expected_s32
  )


figure_manifest_s32 <-
  data.frame(
    file =
      all_expected_figures_s32,
    category =
      c(
        rep(
          "Main",
          length(
            main_expected_s32
          )
        ),
        rep(
          "Individual",
          length(
            individual_expected_s32
          )
        ),
        rep(
          "Supplementary",
          length(
            supp_expected_s32
          )
        )
      ),
    exists =
      file.exists(
        all_expected_figures_s32
      ),
    size_bytes =
      file.info(
        all_expected_figures_s32
      )$size,
    stringsAsFactors = FALSE
  )

write.csv(
  figure_manifest_s32,
  "05_figures/manuscript/Publication_figure_manifest.csv",
  row.names = FALSE
)

cat(
  "\n========================================\n",
  "FINAL SCRIPT 32 FIGURE AUDIT\n",
  "========================================\n"
)

cat(
  "Canonical main figures:",
  length(main_plot_list_s32),
  "\n"
)


cat(
  "Workflow steps:",
  nrow(workflow_data),
  "\n"
)


cat(
  "Workflow external GEO subtype cohorts: 8\n"
)


cat(
  "Workflow LUAD MPCDS external cohorts: 4\n"
)

cat(
  "Individual manuscript panels:",
  length(individual_expected_s32) / 3L,
  "\n"
)


cat(
  "Additional supplementary figures:",
  length(supp_expected_s32) / 3L,
  "\n"
)


cat(
  "Canonical structure: main / individual / supplementary\n"
)


cat(
  "Formats: PNG + PDF + TIFF\n"
)

print(
  figure_manifest_s32,
  row.names = FALSE
)

stopifnot(
  all(figure_manifest_s32$exists),
  all(figure_manifest_s32$size_bytes > 0)
)

cat(
  "\nEstablished final-figure regression checks: PASSED\n"
)

cat(
  "\n========================================\n",
  "SCRIPT 32 FINAL FIGURE GENERATION COMPLETED SUCCESSFULLY\n",
  "========================================\n"
)