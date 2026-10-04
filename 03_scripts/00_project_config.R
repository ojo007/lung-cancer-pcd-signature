# ============================================================
# MSc Lung Cancer PCD Project
# Script: 00_project_config.R
#
# Purpose:
#   Define the canonical project root, shared directories,
#   external staging paths, and small reusable helper functions.
#
# IMPORTANT:
#   - This script performs no statistical analysis.
#   - All downstream scripts should source this file first.
# ============================================================

project_root <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"

if (!dir.exists(project_root)) {
  stop(
    "Project root does not exist:\n",
    project_root,
    "\nUpdate project_root in 03_scripts/00_project_config.R."
  )
}

raw_dir        <- file.path(project_root, "01_raw_data")
processed_dir  <- file.path(project_root, "02_processed_data")
scripts_dir    <- file.path(project_root, "03_scripts")
results_dir    <- file.path(project_root, "04_results")
figures_dir    <- file.path(project_root, "05_figures")
tables_dir     <- file.path(project_root, "06_tables")
manuscript_dir <- file.path(project_root, "07_manuscript")

project_dirs <- c(
  raw_dir,
  processed_dir,
  scripts_dir,
  results_dir,
  figures_dir,
  tables_dir,
  manuscript_dir
)

for (d in project_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

manuscript_figure_dir <- file.path(figures_dir, "manuscript")
main_figure_dir        <- file.path(manuscript_figure_dir, "main")
supp_figure_dir        <- file.path(manuscript_figure_dir, "supplementary")
individual_figure_dir  <- file.path(manuscript_figure_dir, "individual")

figure_subdirs <- c(
  main_figure_dir,
  supp_figure_dir,
  file.path(individual_figure_dir, "clustering"),
  file.path(individual_figure_dir, "immune"),
  file.path(individual_figure_dir, "enrichment"),
  file.path(individual_figure_dir, "mutation"),
  file.path(individual_figure_dir, "cnv"),
  file.path(individual_figure_dir, "methylation"),
  file.path(individual_figure_dir, "mpcds"),
  file.path(individual_figure_dir, "external_validation"),
  file.path(individual_figure_dir, "drug_sensitivity"),
  file.path(individual_figure_dir, "immunotherapy"),
  file.path(individual_figure_dir, "pcd_pathways")
)

for (d in figure_subdirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

gdc_client <- "C:/Users/ibrah/gdc-client/gdc-client.exe"
gdc_staging_root <- "C:/Users/ibrah/GDC"

gdc_luad_rnaseq_stage <- file.path(gdc_staging_root, "GDC_LUAD")
gdc_lusc_rnaseq_stage <- file.path(gdc_staging_root, "GDC_LUSC")
gdc_luad_mut_stage    <- file.path(gdc_staging_root, "GDC_LUAD_MUT")
gdc_lusc_mut_stage    <- file.path(gdc_staging_root, "GDC_LUSC_MUT")
gdc_luad_cnv_stage    <- file.path(gdc_staging_root, "GDC_LUAD_CNV")
gdc_lusc_cnv_stage    <- file.path(gdc_staging_root, "GDC_LUSC_CNV")
gdc_luad_meth_stage   <- file.path(gdc_staging_root, "GDC_LUAD_METH")
gdc_lusc_meth_stage   <- file.path(gdc_staging_root, "GDC_LUSC_METH")

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) {
    stop("Could not create directory:\n", path)
  }
  invisible(path)
}

check_files_exist <- function(paths, label = "Required file(s)") {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      label,
      " missing:\n",
      paste(missing, collapse = "\n")
    )
  }
  invisible(TRUE)
}

check_dirs_exist <- function(paths, label = "Required directorie(s)") {
  missing <- paths[!dir.exists(paths)]
  if (length(missing) > 0) {
    stop(
      label,
      " missing:\n",
      paste(missing, collapse = "\n")
    )
  }
  invisible(TRUE)
}

setwd(project_root)

cat(
  "\n========================================\n",
  "PROJECT CONFIGURATION LOADED\n",
  "========================================\n",
  "Project root: ", project_root, "\n",
  sep = ""
)

# ============================================================
# End of script
# ============================================================