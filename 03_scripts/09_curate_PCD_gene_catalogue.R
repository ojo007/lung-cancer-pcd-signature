# ============================================================
# MSc Lung Cancer PCD Project
# Script: 09_curate_PCD_gene_catalogue.R
#
# Purpose:
#   Curate a reproducible programmed cell death (PCD)
#   gene reference for:
#
#   - Apoptosis
#   - Necroptosis
#   - Pyroptosis
#   - Ferroptosis
#   - Cuproptosis
#
# Sources:
#   Apoptosis:
#     REACTOME_APOPTOSIS
#
#   Necroptosis:
#     REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS
#
#   Pyroptosis:
#     REACTOME_PYROPTOSIS
#
#   Ferroptosis:
#     GOBP_FERROPTOSIS
#     WP_FERROPTOSIS
#
#   Cuproptosis:
#     Literature-curated 12-gene core set
#
# Final curated set:
#   296 unique PCD genes
#
# Important:
#   Genes belonging to multiple PCD mechanisms are retained
#   in all applicable mechanisms.
# ============================================================


# ============================================================
# 0. LOAD SHARED PROJECT CONFIGURATION
# ============================================================

source("03_scripts/00_project_config.R")


# ============================================================
# 1. REQUIRED PACKAGES AND PROJECT PATHS
# ============================================================

required_packages <- c(
  "msigdbr",
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
  library(msigdbr)
  library(dplyr)
})

pcd_raw_dir <- file.path(
  raw_dir,
  "PCD_genes"
)

pcd_processed_dir <- file.path(
  processed_dir,
  "PCD_genes"
)

ensure_dir(pcd_raw_dir)
ensure_dir(pcd_processed_dir)


# ============================================================
# 2. LOAD HUMAN MSigDB
# ============================================================

msig_human <- msigdbr(
  species = "Homo sapiens"
)

# Record the software/database context used for this exact curation run.
# The fixed count checks below intentionally stop the script if a future
# MSigDB release changes the selected gene-set membership.
pcd_curation_provenance <- data.frame(
  item = c(
    "R_version",
    "msigdbr_version",
    "dplyr_version",
    "species",
    "curation_date"
  ),
  value = c(
    R.version.string,
    as.character(packageVersion("msigdbr")),
    as.character(packageVersion("dplyr")),
    "Homo sapiens",
    as.character(Sys.Date())
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 3. SELECT PCD PATHWAY GENE SETS
# ============================================================

selected_msig_sets <- c(
  "REACTOME_APOPTOSIS",
  "REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS",
  "REACTOME_PYROPTOSIS",
  "GOBP_FERROPTOSIS",
  "WP_FERROPTOSIS"
)


selected_set_check <- msig_human |>
  dplyr::filter(
    gs_name %in% selected_msig_sets
  ) |>
  dplyr::distinct(
    gs_name,
    gs_collection,
    gs_subcollection
  )


stopifnot(
  nrow(
    selected_set_check
  ) == 5
)


selected_set_sizes_final <- msig_human |>
  dplyr::filter(
    gs_name %in% selected_msig_sets
  ) |>
  dplyr::group_by(
    gs_name
  ) |>
  dplyr::summarise(
    n_genes =
      dplyr::n_distinct(
        gene_symbol
      ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    gs_name
  )


# ============================================================
# 4. EXTRACT MSigDB PCD GENES
# ============================================================

pcd_msig_long <- msig_human |>
  dplyr::filter(
    gs_name %in% selected_msig_sets
  ) |>
  dplyr::select(
    gene_symbol,
    gs_name,
    gs_collection,
    gs_subcollection
  ) |>
  dplyr::distinct() |>
  dplyr::mutate(
    PCD_type =
      dplyr::case_when(
        
        gs_name ==
          "REACTOME_APOPTOSIS" ~
          "Apoptosis",
        
        gs_name ==
          "REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS" ~
          "Necroptosis",
        
        gs_name ==
          "REACTOME_PYROPTOSIS" ~
          "Pyroptosis",
        
        gs_name %in% c(
          "GOBP_FERROPTOSIS",
          "WP_FERROPTOSIS"
        ) ~
          "Ferroptosis",
        
        TRUE ~
          NA_character_
      )
  )


stopifnot(
  !anyNA(
    pcd_msig_long$PCD_type
  )
)


pcd_msig_long$role <-
  NA_character_


# ============================================================
# 5. CUPROPTOSIS CORE GENES
# ============================================================

cuproptosis_genes <- c(
  "FDX1",
  "LIAS",
  "LIPT1",
  "DLD",
  "DLAT",
  "PDHA1",
  "PDHB",
  "MTF1",
  "GLS",
  "CDKN2A",
  "SLC31A1",
  "ATP7B"
)


stopifnot(
  length(
    unique(
      cuproptosis_genes
    )
  ) == 12
)


cuproptosis_table <- data.frame(
  
  gene_symbol =
    cuproptosis_genes,
  
  role = c(
    rep(
      "Pro-cuproptosis",
      7
    ),
    rep(
      "Anti-cuproptosis",
      3
    ),
    rep(
      "Copper transporter",
      2
    )
  ),
  
  PCD_type =
    "Cuproptosis",
  
  source_gene_set =
    "Literature-curated cuproptosis core genes",
  
  stringsAsFactors = FALSE
)


cuproptosis_long <- cuproptosis_table |>
  dplyr::transmute(
    
    gene_symbol =
      gene_symbol,
    
    gs_name =
      source_gene_set,
    
    gs_collection =
      "Literature",
    
    gs_subcollection =
      "Cuproptosis core genes",
    
    PCD_type =
      PCD_type,
    
    role =
      role
  )


# ============================================================
# 6. COMBINE ALL PCD SOURCES
# ============================================================

pcd_master_long <- dplyr::bind_rows(
  pcd_msig_long,
  cuproptosis_long
)


pcd_master_long_final <- pcd_master_long |>
  dplyr::rename(
    
    source_gene_set =
      gs_name,
    
    source_collection =
      gs_collection,
    
    source_subcollection =
      gs_subcollection
  ) |>
  dplyr::arrange(
    PCD_type,
    gene_symbol,
    source_gene_set
  )


stopifnot(
  nrow(
    pcd_master_long_final
  ) == 332
)


# ============================================================
# 7. UNIQUE GENE COUNTS BY PCD TYPE
# ============================================================

pcd_unique_counts <- pcd_master_long_final |>
  dplyr::group_by(
    PCD_type
  ) |>
  dplyr::summarise(
    
    unique_genes =
      dplyr::n_distinct(
        gene_symbol
      ),
    
    .groups =
      "drop"
  ) |>
  dplyr::arrange(
    PCD_type
  )


# Expected:
#
# Apoptosis       166
# Cuproptosis      12
# Ferroptosis      89
# Necroptosis      32
# Pyroptosis       27


# ============================================================
# 8. ONE-ROW-PER-GENE MASTER TABLE
# ============================================================

pcd_master_gene <- pcd_master_long_final |>
  dplyr::group_by(
    gene_symbol
  ) |>
  dplyr::summarise(
    
    n_PCD_types =
      dplyr::n_distinct(
        PCD_type
      ),
    
    PCD_types =
      paste(
        sort(
          unique(
            PCD_type
          )
        ),
        collapse = "; "
      ),
    
    source_gene_sets =
      paste(
        sort(
          unique(
            source_gene_set
          )
        ),
        collapse = "; "
      ),
    
    source_collections =
      paste(
        sort(
          unique(
            source_collection
          )
        ),
        collapse = "; "
      ),
    
    cuproptosis_role =
      paste(
        unique(
          stats::na.omit(
            role
          )
        ),
        collapse = "; "
      ),
    
    .groups =
      "drop"
  ) |>
  
  dplyr::mutate(
    
    cuproptosis_role =
      dplyr::if_else(
        cuproptosis_role == "",
        NA_character_,
        cuproptosis_role
      )
  ) |>
  
  dplyr::arrange(
    gene_symbol
  )


stopifnot(
  nrow(
    pcd_master_gene
  ) == 296
)


stopifnot(
  sum(
    duplicated(
      pcd_master_gene$gene_symbol
    )
  ) == 0
)


stopifnot(
  !anyNA(
    pcd_master_gene$gene_symbol
  )
)


# ============================================================
# 9. OVERLAPPING PCD GENES
# ============================================================

pcd_overlap_genes <- pcd_master_gene |>
  dplyr::filter(
    n_PCD_types > 1
  )


stopifnot(
  nrow(
    pcd_overlap_genes
  ) == 29
)


# Expected membership distribution:
#
# 267 genes -> one PCD type
#  28 genes -> two PCD types
#   1 gene  -> three PCD types


stopifnot(
  identical(
    as.integer(
      table(
        pcd_master_gene$n_PCD_types
      )
    ),
    c(
      267L,
      28L,
      1L
    )
  )
)


# ============================================================
# 10. SAVE RAW SOURCE RECORDS
# ============================================================

write.csv(
  pcd_master_long_final,
  file.path(
    pcd_raw_dir,
    "PCD_gene_sources_master_long.csv"
  ),
  row.names = FALSE
)


write.csv(
  selected_set_sizes_final,
  file.path(
    pcd_raw_dir,
    "PCD_selected_MSigDB_gene_sets.csv"
  ),
  row.names = FALSE
)


write.csv(
  cuproptosis_table,
  file.path(
    pcd_raw_dir,
    "Cuproptosis_literature_curated_genes.csv"
  ),
  row.names = FALSE
)


write.csv(
  pcd_curation_provenance,
  file.path(
    pcd_raw_dir,
    "PCD_curation_provenance.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. SAVE PROCESSED PCD REFERENCES
# ============================================================

write.csv(
  pcd_master_gene,
  file.path(
    pcd_processed_dir,
    "PCD_master_unique_genes.csv"
  ),
  row.names = FALSE
)


write.csv(
  pcd_overlap_genes,
  file.path(
    pcd_processed_dir,
    "PCD_overlapping_genes.csv"
  ),
  row.names = FALSE
)


write.csv(
  pcd_unique_counts,
  file.path(
    pcd_processed_dir,
    "PCD_gene_counts_by_type.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. FINAL SUMMARY
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "PCD gene curation complete.\n"
)

cat(
  "========================================\n"
)


cat(
  "\nUnique genes by PCD mechanism:\n"
)

print(
  pcd_unique_counts
)


cat(
  "\nTotal unique PCD genes:\n"
)

print(
  nrow(
    pcd_master_gene
  )
)


cat(
  "\nGenes belonging to >1 PCD mechanism:\n"
)

print(
  nrow(
    pcd_overlap_genes
  )
)


cat(
  "\nPCD membership distribution:\n"
)

print(
  table(
    pcd_master_gene$n_PCD_types
  )
)


cat(
  "\nGene belonging to >=3 mechanisms:\n"
)

print(
  pcd_master_gene |>
    dplyr::filter(
      n_PCD_types >= 3
    )
)


cat(
  "\n========================================\n",
  "SCRIPT 09 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  sep = ""
)


# ============================================================
# End of script
# ============================================================