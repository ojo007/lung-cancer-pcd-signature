# ============================================================
# Script 19: Prepare external GEO cohorts
#
# Purpose:
#   Prepare multiple independent GEO cohorts for external
#   validation of TCGA-derived PCD molecular subtypes and,
#   where appropriate, the LUAD MPCDS.
#
# Cohorts selected:
#
#   LUAD:
#     GSE68465
#     GSE72094
#     GSE31210
#     GSE50081 (adenocarcinoma subset only)
#
#   LUSC:
#     GSE73403
#     GSE4573
#     GSE157010
#     GSE30219 (squamous subset only)
#
# Important:
#   - Expression preprocessing is performed independently
#     within each GEO cohort.
#   - No ComBat is applied across GEO and TCGA.
#   - No survival outcome is used to select samples, genes,
#     probes, or scaling parameters.
#   - Probe-level values are collapsed to genes using the
#     median across all probes mapped to each gene, consistent
#     with the existing GSE68465 workflow.
#   - Gene-wise z-scores are computed within each GEO cohort
#     for later fixed-centroid subtype projection.
# ============================================================


# ------------------------------------------------------------
# 0. Shared project configuration
# ------------------------------------------------------------

source("03_scripts/00_project_config.R")

# Large GEO series matrices can exceed GEOquery's default
# 120-second transfer window. Use a 30-minute timeout for
# NCBI downloads while preserving the rest of the workflow.
options(
  timeout = max(
    1800,
    getOption("timeout")
  ),
  GEOquery.download.timeout = 1800
)



# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_packages <- c(
  "GEOquery",
  "Biobase",
  "dplyr",
  "tidyr"
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

geo_raw_dir <- file.path(
  raw_dir,
  "GEO"
)

geo_processed_dir <- file.path(
  processed_dir,
  "GEO",
  "external_validation"
)

ensure_dir(geo_raw_dir)
ensure_dir(geo_processed_dir)


# ------------------------------------------------------------
# 3. Required TCGA-derived gene lists
# ------------------------------------------------------------

pcd_master_file <- file.path(
  processed_dir,
  "PCD_genes",
  "PCD_master_unique_genes.csv"
)

candidate_gene_file <- file.path(
  results_dir,
  "differential_expression",
  "TCGA_concordant_PCD_candidates_LUAD_LUSC.csv"
)

check_files_exist(
  c(
    pcd_master_file,
    candidate_gene_file
  )
)

pcd_master <- read.csv(
  pcd_master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_genes <- read.csv(
  candidate_gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  "gene_symbol" %in% colnames(pcd_master),
  "gene_name" %in% colnames(candidate_genes)
)

pcd_genes <- unique(
  pcd_master$gene_symbol
)

model_genes <- unique(
  candidate_genes$gene_name
)

stopifnot(
  length(pcd_genes) == 296,
  length(model_genes) == 47
)


# ------------------------------------------------------------
# 4. Historical gene-symbol rescue
# ------------------------------------------------------------

legacy_to_current <- c(
  "SEPT4" = "SEPTIN4",
  "GPR56" = "ADGRG1",
  "HIST1H1A" = "H1-1",
  "HIST1H1C" = "H1-2",
  "HIST1H1D" = "H1-3",
  "HIST1H1E" = "H1-4",
  "HIST1H1B" = "H1-5",
  "H1F0" = "H1-0",
  "NUCLING" = "UACA"
)


# ------------------------------------------------------------
# 5. External-cohort manifest
#
# selection_mode:
#   "target_regex" = require phenotype evidence of the target
#                    histology and subset accordingly.
#   "histology_specific" = cohort is histology-specific;
#                          phenotype target matches are used
#                          when available, otherwise obvious
#                          normal/control samples are excluded.
# ------------------------------------------------------------

cohort_manifest <- data.frame(
  accession = c(
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
    "LUAD",
    "LUAD",
    "LUAD",
    "LUAD",
    "LUSC",
    "LUSC",
    "LUSC",
    "LUSC"
  ),
  selection_mode = c(
    "histology_specific",
    "histology_specific",
    "histology_specific",
    "target_regex",
    "histology_specific",
    "histology_specific",
    "histology_specific",
    "target_regex"
  ),
  target_regex = c(
    "adenocarcinoma",
    "adenocarcinoma",
    "adenocarcinoma|\\badc\\b",
    "adenocarcinoma",
    "squamous|lsqc|lusc",
    "squamous|scc",
    "squamous|scc",
    "squamous|scc"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cohort_manifest,
  file.path(
    geo_processed_dir,
    "external_GEO_cohort_manifest.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Utility: GEO download with retry
# ------------------------------------------------------------

get_geo_with_retry <- function(
    accession,
    raw_cohort_dir,
    max_attempts = 3
) {
  
  last_error <- NULL
  
  for (
    attempt in seq_len(
      max_attempts
    )
  ) {
    
    cat(
      "GEO download attempt ",
      attempt,
      "/",
      max_attempts,
      " for ",
      accession,
      "\n",
      sep = ""
    )
    
    result <- tryCatch(
      {
        GEOquery::getGEO(
          accession,
          GSEMatrix = TRUE,
          getGPL = TRUE,
          destdir = raw_cohort_dir
        )
      },
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    
    if (!is.null(result)) {
      return(result)
    }
    
    if (attempt < max_attempts) {
      Sys.sleep(
        5 * attempt
      )
    }
  }
  
  stop(
    accession,
    ": GEO download failed after ",
    max_attempts,
    " attempts. Last error: ",
    conditionMessage(
      last_error
    )
  )
}


# ------------------------------------------------------------
# 7. Utility: choose one ExpressionSet
# ------------------------------------------------------------

choose_eset <- function(
    gse_object,
    accession
) {
  
  if (
    methods::is(
      gse_object,
      "ExpressionSet"
    )
  ) {
    return(gse_object)
  }
  
  if (!is.list(gse_object)) {
    stop(
      accession,
      ": GEOquery returned an unsupported object."
    )
  }
  
  if (length(gse_object) == 1) {
    return(gse_object[[1]])
  }
  
  n_samples <- vapply(
    gse_object,
    function(x) {
      if (
        methods::is(
          x,
          "ExpressionSet"
        )
      ) {
        ncol(
          Biobase::exprs(x)
        )
      } else {
        0L
      }
    },
    integer(1)
  )
  
  if (all(n_samples == 0)) {
    stop(
      accession,
      ": no usable ExpressionSet was found."
    )
  }
  
  chosen <- which.max(
    n_samples
  )
  
  message(
    accession,
    ": multiple ExpressionSets found; using the largest (",
    n_samples[chosen],
    " samples)."
  )
  
  gse_object[[chosen]]
}


# ------------------------------------------------------------
# 7. Utility: phenotype text for robust cohort filtering
# ------------------------------------------------------------

collapse_pheno_text <- function(
    pheno
) {
  
  apply(
    pheno,
    1,
    function(x) {
      paste(
        x,
        collapse = " | "
      )
    }
  )
}


# ------------------------------------------------------------
# 8. Utility: select target histology samples
# ------------------------------------------------------------

select_target_samples <- function(
    pheno,
    target_regex,
    selection_mode,
    accession
) {
  
  if (!"geo_accession" %in% colnames(pheno)) {
    stop(
      accession,
      ": phenotype table has no geo_accession column."
    )
  }
  
  combined_text <- tolower(
    collapse_pheno_text(
      pheno
    )
  )
  
  target_hit <- grepl(
    target_regex,
    combined_text,
    ignore.case = TRUE,
    perl = TRUE
  )
  
  obvious_non_tumor <- grepl(
    paste0(
      "\\bnormal\\b|",
      "non[- ]?tumou?r|",
      "adjacent normal|",
      "healthy|",
      "control lung"
    ),
    combined_text,
    ignore.case = TRUE,
    perl = TRUE
  )
  
  if (
    selection_mode ==
    "target_regex"
  ) {
    
    if (sum(target_hit) < 20) {
      stop(
        accession,
        ": fewer than 20 samples matched target histology regex. ",
        "Matched = ",
        sum(target_hit),
        ". Inspect the saved phenotype table before proceeding."
      )
    }
    
    keep <- target_hit
    
  } else {
    
    if (sum(target_hit) >= 20) {
      
      keep <- target_hit
      
    } else {
      
      keep <- !obvious_non_tumor
      
      message(
        accession,
        ": phenotype did not contain >=20 explicit target-histology ",
        "matches; using histology-specific cohort definition and ",
        "excluding obvious normal/control samples."
      )
    }
  }
  
  selected_ids <- pheno$geo_accession[
    keep
  ]
  
  if (length(selected_ids) < 20) {
    stop(
      accession,
      ": fewer than 20 target samples remained after filtering."
    )
  }
  
  selected_ids
}


# ------------------------------------------------------------
# 9. Utility: identify gene-symbol annotation column
# ------------------------------------------------------------

find_symbol_column <- function(
    feature_table
) {
  
  nms <- colnames(
    feature_table
  )
  
  priority_exact <- c(
    "Gene Symbol",
    "Gene symbol",
    "GENE_SYMBOL",
    "GENE SYMBOL",
    "GeneSymbol",
    "gene_symbol",
    "Symbol",
    "SYMBOL"
  )
  
  exact_hit <- priority_exact[
    priority_exact %in%
      nms
  ]
  
  if (length(exact_hit) > 0) {
    return(exact_hit[1])
  }
  
  regex_hit <- grep(
    "gene.*symbol|symbol.*gene|(^|_)symbol($|_)",
    nms,
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(regex_hit) > 0) {
    return(regex_hit[1])
  }
  
  assignment_hit <- grep(
    "gene.assignment|gene_assignment",
    nms,
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(assignment_hit) > 0) {
    return(assignment_hit[1])
  }
  
  NA_character_
}


# ------------------------------------------------------------
# 10. Utility: obtain a platform annotation table
# ------------------------------------------------------------

get_feature_annotation <- function(
    eset,
    accession,
    raw_cohort_dir
) {
  
  feature_table <- Biobase::fData(
    eset
  )
  
  symbol_column <- find_symbol_column(
    feature_table
  )
  
  if (!is.na(symbol_column)) {
    
    feature_table$probe_id <- rownames(
      feature_table
    )
    
    return(
      list(
        feature_table = feature_table,
        symbol_column = symbol_column,
        annotation_source = "ExpressionSet_fData"
      )
    )
  }
  
  platform_id <- Biobase::annotation(
    eset
  )
  
  if (
    is.na(platform_id) ||
    platform_id == ""
  ) {
    stop(
      accession,
      ": no usable gene-symbol annotation in fData and no GPL ID."
    )
  }
  
  message(
    accession,
    ": downloading platform annotation ",
    platform_id,
    "."
  )
  
  gpl <- GEOquery::getGEO(
    platform_id,
    destdir = raw_cohort_dir
  )
  
  platform_table <- GEOquery::Table(
    gpl
  )
  
  symbol_column <- find_symbol_column(
    platform_table
  )
  
  if (is.na(symbol_column)) {
    stop(
      accession,
      ": could not identify a gene-symbol column in ",
      platform_id,
      ". Available columns: ",
      paste(
        colnames(platform_table),
        collapse = ", "
      )
    )
  }
  
  id_candidates <- c(
    "ID",
    "ID_REF",
    "ProbeID",
    "PROBE_ID",
    "SPOT_ID"
  )
  
  id_column <- id_candidates[
    id_candidates %in%
      colnames(platform_table)
  ]
  
  if (length(id_column) == 0) {
    id_column <- colnames(
      platform_table
    )[1]
  } else {
    id_column <- id_column[1]
  }
  
  platform_table$probe_id <- as.character(
    platform_table[[id_column]]
  )
  
  list(
    feature_table = platform_table,
    symbol_column = symbol_column,
    annotation_source = paste0(
      "GPL:",
      platform_id
    )
  )
}


# ------------------------------------------------------------
# 11. Utility: clean and expand probe-to-symbol mappings
# ------------------------------------------------------------

build_probe_map <- function(
    feature_table,
    symbol_column
) {
  
  map <- data.frame(
    probe_id = as.character(
      feature_table$probe_id
    ),
    gene_symbol = as.character(
      feature_table[[symbol_column]]
    ),
    stringsAsFactors = FALSE
  )
  
  map <- map |>
    dplyr::filter(
      !is.na(gene_symbol),
      gene_symbol != "",
      gene_symbol != "---"
    )
  
  # Standard Affymetrix and many GEO annotations use "///"
  # for multiple mapped symbols. Some use "//" or ";".
  map <- map |>
    tidyr::separate_rows(
      gene_symbol,
      sep = "\\s*///\\s*|\\s*//\\s*|\\s*;\\s*"
    ) |>
    dplyr::mutate(
      gene_symbol = trimws(
        gene_symbol
      )
    ) |>
    dplyr::filter(
      gene_symbol != "",
      gene_symbol != "---"
    )
  
  # For assignment-style annotations, remove accession /
  # description fragments and retain tokens resembling symbols.
  map$gene_symbol <- sub(
    "^.*?\\|",
    "",
    map$gene_symbol
  )
  
  map$gene_symbol <- trimws(
    map$gene_symbol
  )
  
  for (
    old_symbol in names(
      legacy_to_current
    )
  ) {
    
    map$gene_symbol[
      map$gene_symbol ==
        old_symbol
    ] <- legacy_to_current[[old_symbol]]
  }
  
  map |>
    dplyr::distinct(
      probe_id,
      gene_symbol
    )
}


# ------------------------------------------------------------
# 12. Utility: determine whether log2 transformation is needed
# ------------------------------------------------------------

ensure_log2_expression <- function(
    expr,
    accession
) {
  
  finite_values <- expr[
    is.finite(expr)
  ]
  
  if (length(finite_values) == 0) {
    stop(
      accession,
      ": expression matrix contains no finite values."
    )
  }
  
  q99 <- as.numeric(
    stats::quantile(
      finite_values,
      probs = 0.99,
      na.rm = TRUE
    )
  )
  
  max_value <- max(
    finite_values,
    na.rm = TRUE
  )
  
  min_value <- min(
    finite_values,
    na.rm = TRUE
  )
  
  # Typical MAS5 / linear-scale arrays have large positive
  # intensities. Already normalized log-expression matrices
  # generally occupy a much smaller range.
  needs_log2 <-
    min_value >= 0 &&
    (
      q99 > 100 ||
        max_value > 1000
    )
  
  if (needs_log2) {
    
    message(
      accession,
      ": applying log2(x + 1) transformation."
    )
    
    expr <- log2(
      expr + 1
    )
    
  } else {
    
    message(
      accession,
      ": expression values appear already log-scaled; ",
      "no additional log2 transformation applied."
    )
  }
  
  attr(
    expr,
    "log2_transform_applied"
  ) <- needs_log2
  
  expr
}


# ------------------------------------------------------------
# 13. Utility: collapse mapped probes to gene-level median
# ------------------------------------------------------------

collapse_probes_to_genes <- function(
    expr,
    probe_map,
    target_genes
) {
  
  probe_map <- probe_map |>
    dplyr::filter(
      gene_symbol %in%
        target_genes,
      probe_id %in%
        rownames(expr)
    ) |>
    dplyr::distinct(
      probe_id,
      gene_symbol
    )
  
  genes_present <- intersect(
    target_genes,
    unique(
      probe_map$gene_symbol
    )
  )
  
  if (length(genes_present) == 0) {
    stop(
      "No target genes could be mapped on this platform."
    )
  }
  
  gene_matrix <- matrix(
    NA_real_,
    nrow = length(
      genes_present
    ),
    ncol = ncol(
      expr
    ),
    dimnames = list(
      genes_present,
      colnames(expr)
    )
  )
  
  for (
    gene in genes_present
  ) {
    
    probes <- probe_map$probe_id[
      probe_map$gene_symbol ==
        gene
    ]
    
    probe_expr <- expr[
      probes,
      ,
      drop = FALSE
    ]
    
    if (nrow(probe_expr) == 1) {
      
      gene_matrix[
        gene,
      ] <- probe_expr[
        1,
      ]
      
    } else {
      
      gene_matrix[
        gene,
      ] <- apply(
        probe_expr,
        2,
        stats::median,
        na.rm = TRUE
      )
    }
  }
  
  gene_matrix
}


# ------------------------------------------------------------
# 14. Utility: gene-wise z-score
# ------------------------------------------------------------

gene_zscore <- function(
    gene_matrix
) {
  
  row_mean <- rowMeans(
    gene_matrix,
    na.rm = TRUE
  )
  
  row_sd <- apply(
    gene_matrix,
    1,
    stats::sd,
    na.rm = TRUE
  )
  
  keep <- is.finite(
    row_sd
  ) &
    row_sd >
    0
  
  z <- sweep(
    gene_matrix[
      keep,
      ,
      drop = FALSE
    ],
    1,
    row_mean[
      keep
    ],
    "-"
  )
  
  z <- sweep(
    z,
    1,
    row_sd[
      keep
    ],
    "/"
  )
  
  z
}


# ------------------------------------------------------------
# 15. Process one GEO cohort
# ------------------------------------------------------------

prepare_geo_cohort <- function(
    accession,
    histology,
    selection_mode,
    target_regex
) {
  
  cat(
    "\n============================================================\n",
    "PREPARING ",
    accession,
    " (",
    histology,
    ")\n",
    "============================================================\n",
    sep = ""
  )
  
  raw_cohort_dir <- file.path(
    geo_raw_dir,
    accession
  )
  
  processed_cohort_dir <- file.path(
    geo_processed_dir,
    accession
  )
  
  ensure_dir(
    raw_cohort_dir
  )
  
  ensure_dir(
    processed_cohort_dir
  )
  
  gse_object <- get_geo_with_retry(
    accession = accession,
    raw_cohort_dir = raw_cohort_dir,
    max_attempts = 3
  )
  
  eset <- choose_eset(
    gse_object,
    accession
  )
  
  expr <- Biobase::exprs(
    eset
  )
  
  pheno <- Biobase::pData(
    eset
  )
  
  if (
    !"geo_accession" %in%
    colnames(pheno)
  ) {
    pheno$geo_accession <- colnames(
      expr
    )
  }
  
  pheno <- pheno[
    match(
      colnames(expr),
      pheno$geo_accession
    ),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      colnames(expr),
      pheno$geo_accession
    )
  )
  
  write.csv(
    pheno,
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_phenotype_full.csv"
      )
    ),
    row.names = FALSE
  )
  
  selected_ids <- select_target_samples(
    pheno = pheno,
    target_regex = target_regex,
    selection_mode = selection_mode,
    accession = accession
  )
  
  expr_selected <- expr[
    ,
    selected_ids,
    drop = FALSE
  ]
  
  pheno_selected <- pheno[
    match(
      selected_ids,
      pheno$geo_accession
    ),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      colnames(expr_selected),
      pheno_selected$geo_accession
    )
  )
  
  expr_log2 <- ensure_log2_expression(
    expr_selected,
    accession
  )
  
  log2_applied <- isTRUE(
    attr(
      expr_log2,
      "log2_transform_applied"
    )
  )
  
  annotation_result <- get_feature_annotation(
    eset = eset,
    accession = accession,
    raw_cohort_dir = raw_cohort_dir
  )
  
  feature_table <- annotation_result$feature_table
  
  # Ensure probe IDs are available for fData-derived annotation.
  if (!"probe_id" %in% colnames(feature_table)) {
    feature_table$probe_id <- rownames(
      feature_table
    )
  }
  
  probe_map <- build_probe_map(
    feature_table = feature_table,
    symbol_column = annotation_result$symbol_column
  )
  
  pcd_expr <- collapse_probes_to_genes(
    expr = expr_log2,
    probe_map = probe_map,
    target_genes = pcd_genes
  )
  
  model_expr <- pcd_expr[
    intersect(
      model_genes,
      rownames(pcd_expr)
    ),
    ,
    drop = FALSE
  ]
  
  pcd_z <- gene_zscore(
    pcd_expr
  )
  
  model_z <- pcd_z[
    intersect(
      model_genes,
      rownames(pcd_z)
    ),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      colnames(pcd_expr),
      selected_ids
    ),
    identical(
      colnames(model_expr),
      selected_ids
    ),
    identical(
      colnames(pcd_z),
      selected_ids
    ),
    identical(
      colnames(model_z),
      selected_ids
    )
  )
  
  pcd_present <- intersect(
    pcd_genes,
    rownames(pcd_expr)
  )
  
  model_present <- intersect(
    model_genes,
    rownames(model_expr)
  )
  
  pcd_missing <- setdiff(
    pcd_genes,
    pcd_present
  )
  
  model_missing <- setdiff(
    model_genes,
    model_present
  )
  
  # A cohort can still be useful with incomplete gene coverage,
  # but low coverage is flagged for review rather than hidden.
  coverage_status <- dplyr::case_when(
    length(model_present) >= 44 ~ "PASS_HIGH",
    length(model_present) >= 38 ~ "PASS_MODERATE",
    TRUE ~ "REVIEW_LOW_COVERAGE"
  )
  
  qc <- data.frame(
    accession = accession,
    histology = histology,
    platform = Biobase::annotation(
      eset
    ),
    n_series_samples = ncol(
      expr
    ),
    n_selected_samples = ncol(
      expr_selected
    ),
    n_expression_features = nrow(
      expr
    ),
    log2_transform_applied = log2_applied,
    annotation_source = annotation_result$annotation_source,
    symbol_column = annotation_result$symbol_column,
    n_PCD_genes_present = length(
      pcd_present
    ),
    PCD_coverage_percent = 100 *
      length(
        pcd_present
      ) /
      length(
        pcd_genes
      ),
    n_candidate47_present = length(
      model_present
    ),
    candidate47_coverage_percent = 100 *
      length(
        model_present
      ) /
      length(
        model_genes
      ),
    coverage_status = coverage_status,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    pheno_selected,
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_phenotype_selected.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    probe_map,
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_probe_gene_mapping.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      gene = pcd_missing,
      stringsAsFactors = FALSE
    ),
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_missing_PCD_genes.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      gene = model_missing,
      stringsAsFactors = FALSE
    ),
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_missing_candidate47_genes.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    qc,
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_preparation_QC.csv"
      )
    ),
    row.names = FALSE
  )
  
  saveRDS(
    pcd_expr,
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_PCD296_expression_log2.rds"
      )
    )
  )
  
  saveRDS(
    pcd_z,
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_PCD296_expression_z.rds"
      )
    )
  )
  
  saveRDS(
    model_expr,
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_candidate47_expression_log2.rds"
      )
    )
  )
  
  saveRDS(
    model_z,
    file.path(
      processed_cohort_dir,
      paste0(
        accession,
        "_candidate47_expression_z.rds"
      )
    )
  )
  
  cat(
    "Selected samples: ",
    ncol(expr_selected),
    "\n",
    "PCD genes mapped: ",
    length(pcd_present),
    "/296\n",
    "Candidate genes mapped: ",
    length(model_present),
    "/47\n",
    "Coverage status: ",
    coverage_status,
    "\n",
    sep = ""
  )
  
  qc
}


# ------------------------------------------------------------
# 16. Prepare all eight cohorts
# ------------------------------------------------------------

qc_results <- vector(
  "list",
  nrow(
    cohort_manifest
  )
)

for (
  i in seq_len(
    nrow(
      cohort_manifest
    )
  )
) {
  
  cohort_row <- cohort_manifest[
    i,
    ,
    drop = FALSE
  ]
  
  qc_results[[i]] <- tryCatch(
    {
      prepare_geo_cohort(
        accession = cohort_row$accession,
        histology = cohort_row$histology,
        selection_mode = cohort_row$selection_mode,
        target_regex = cohort_row$target_regex
      )
    },
    error = function(e) {
      
      warning(
        cohort_row$accession,
        " preparation failed: ",
        conditionMessage(e)
      )
      
      data.frame(
        accession = cohort_row$accession,
        histology = cohort_row$histology,
        platform = NA_character_,
        n_series_samples = NA_integer_,
        n_selected_samples = NA_integer_,
        n_expression_features = NA_integer_,
        log2_transform_applied = NA,
        annotation_source = NA_character_,
        symbol_column = NA_character_,
        n_PCD_genes_present = NA_integer_,
        PCD_coverage_percent = NA_real_,
        n_candidate47_present = NA_integer_,
        candidate47_coverage_percent = NA_real_,
        coverage_status = paste0(
          "FAILED: ",
          conditionMessage(e)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
}

qc_summary <- dplyr::bind_rows(
  qc_results
)

write.csv(
  qc_summary,
  file.path(
    geo_processed_dir,
    "external_GEO_preparation_QC_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Final validation and console summary
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "SCRIPT 19 EXTERNAL GEO PREPARATION SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  qc_summary
)

n_success <- sum(
  !grepl(
    "^FAILED:",
    qc_summary$coverage_status
  )
)

n_high_or_moderate <- sum(
  qc_summary$coverage_status %in%
    c(
      "PASS_HIGH",
      "PASS_MODERATE"
    )
)

cat(
  "\nCohorts requested: ",
  nrow(cohort_manifest),
  "\n",
  "Cohorts prepared successfully: ",
  n_success,
  "\n",
  "Cohorts with >=38/47 candidate genes: ",
  n_high_or_moderate,
  "\n",
  "Outputs: ",
  geo_processed_dir,
  "\n",
  sep = ""
)

if (
  n_success <
  nrow(cohort_manifest)
) {
  
  cat(
    "\nIMPORTANT: One or more cohorts require inspection. ",
    "Do not proceed to Script 20 until failed cohorts are resolved.\n",
    sep = ""
  )
  
} else {
  
  cat(
    "\nAll eight external GEO cohorts were prepared successfully.\n"
  )
}

cat(
  "============================================================\n"
)

# ============================================================
# End of script
# ============================================================