## ============================================================================
## Both_CosMXNew_GeoMX_Pathway_validation.R
## ============================================================================
## Purpose: Pathway validation comparing New CosMX (DFCI + Rochester) with
##          GeoMX (DFCI + Rochester) using comprehensive GSVA heatmaps
##
## Input:
##   - New_CosMX_Top1_BatchCorrected.qs (from CosMX4TMANew_1_3_EBV_celltype_correlation.R)
##   - GeoMX_DFCI_Rochester_top1.qs (batch-corrected GeoMX data)
##   - Customized pathway gene sets (from 5-Preload_customized_gene_list.R)
##
## Output:
##   - Four-level GSVA heatmaps with pooling strategies
##   - Lineage GSVA confidence stratification plots
##   - Combined pathway validation results
##
## Author: New CosMX analysis pipeline
## Date: 2025
## ============================================================================

## ---------------------------------------------------------------------------
## Libraries
## ---------------------------------------------------------------------------
library(Seurat)
library(SpatialExperiment)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(GSVA)
library(enrichR)
library(qs)
library(knitr)
library(DT)
library(grid)
library(RColorBrewer)

## ---------------------------------------------------------------------------
## Setup: Base paths and dependencies
## ---------------------------------------------------------------------------

base_root <- "/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/CosMXNew_code_data_publish/"
standr_path <- file.path(base_root, "standR_covariate")
output_base <- file.path(base_root, "Data/CosMXNew/Exampleoutpu")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
# Path for loading functions
loadingfunction_wdpath <- file.path(base_root, "Code/src")
# Source pathway and gene list definitions
source(file.path(loadingfunction_wdpath, "5-Preload_customized_gene_list.R"))
# Download standR_covariate package to local directory if not present
if (!dir.exists(standr_path)) {
  cat("Downloading standR_covariate from GitHub to local directory...\n")
  system(paste0("git clone https://github.com/BMEngineeR/standR_covariate.git ", standr_path))
}

# Load the package from local directory
devtools::load_all(standr_path)

## ---------------------------------------------------------------------------
## Paths and loading of gene sets / helper functions
## ---------------------------------------------------------------------------

# Path to saved top1 objects
geomx_dir_objects <- file.path(base_root, "Data/Three_Cohorts_Top_Objects/")

# New CosMX results directory
cosmx_dir_objects <- file.path(output_base, "CosMX4TMANew_1_3_EBV_celltype_correlation_finer")

# Output directory for pathway validation
results_dir <- file.path(output_base, "Pathway_Validation")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}






## Helper to rename pathways with gene counts
rename_pathways_with_counts <- function(pathway_list, dataset_genes) {
  dataset_genes <- as.character(dataset_genes)

  new_names <- sapply(names(pathway_list), function(pathway_name) {
    pathway_genes <- pathway_list[[pathway_name]]
    overlap_genes <- intersect(pathway_genes, dataset_genes)
    n_overlap <- length(overlap_genes)
    n_total <- length(pathway_genes)
    paste0("(", n_overlap, "/", n_total, ")_", pathway_name)
  })

  names(pathway_list) <- new_names
  pathway_list
}

## ---------------------------------------------------------------------------
## Load New CosMX (DFCI + Rochester)
## ---------------------------------------------------------------------------

cat("=== Loading New CosMX batch-corrected data ===\n")

cosmx_new_path <- file.path(cosmx_dir_objects, "New_CosMX_Top1_BatchCorrected_finer_subsetTMacroTumor.qs")
if (!file.exists(cosmx_new_path)) {
  stop("New CosMX batch-corrected file not found: ", cosmx_new_path)
}

cosmx_new_spe_obj <- qread(cosmx_new_path)
cat("New CosMX object loaded:\n")
cat(" - Dimensions:", nrow(cosmx_new_spe_obj), "genes x", ncol(cosmx_new_spe_obj), "samples\n")
cat(" - Sites:", paste(unique(colData(cosmx_new_spe_obj)$TMA), collapse = ", "), "\n")

# Standardize metadata to match template structure
cosmx_new_spe_obj$shared_EBV <- cosmx_new_spe_obj$EBV_status
cosmx_new_spe_obj$shared_site <- cosmx_new_spe_obj$TMA
cosmx_new_spe_obj$shared_tech <- "CosMX_New"
cosmx_new_spe_obj$shared_FOV <- paste0("FOV_", cosmx_new_spe_obj$coreName)

# Map cell types: "B cell" -> "Tumor", "T cell" -> "Tcell", "Macrophage" -> "Macrophage"
cosmx_new_spe_obj$shared_CT <- dplyr::recode(
  cosmx_new_spe_obj$cell_type,
  "B cell" = "Tumor",
  "CD4 T" = "Tcell",
  "CD8 T" = "Tcell",
  "T cell" = "Tcell",
  "Macrophage" = "Macrophage"
)

# Filter to DFCI and Rochester only
keep_sites <- c("DFCI", "Rochester")
cosmx_new_keep <- cosmx_new_spe_obj$shared_site %in% keep_sites
cosmx_new_spe_obj <- cosmx_new_spe_obj[, cosmx_new_keep]

cosmx_new_spe_obj$shared_ruv_W1 <- cosmx_new_spe_obj$ruv_W1
cosmx_new_spe_obj$shared_cell_counts <- cosmx_new_spe_obj$cell_count
cosmx_new_spe_obj$shared_librarysize <- colSums(cosmx_new_spe_obj@assays@data$counts)

cat("After filtering to DFCI/Rochester:\n")
cat(" - Samples:", ncol(cosmx_new_spe_obj), "\n")
cat(" - Sites:", paste(unique(cosmx_new_spe_obj$shared_site), collapse = ", "), "\n")
cat(" - Cell types:", paste(unique(cosmx_new_spe_obj$shared_CT), collapse = ", "), "\n")
cat(" - EBV status:\n")
print(table(cosmx_new_spe_obj$shared_EBV))

## ---------------------------------------------------------------------------
## Load GeoMX (DFCI + Rochester)
## ---------------------------------------------------------------------------

cat("\n=== Loading GeoMX batch-corrected data ===\n")
# this code refer to GeoMX_1_3_object_producing.R
geomx_save <- qread(file.path(geomx_dir_objects, "GeoMX_DFCI_Rochester_top1.qs"))
geomx_spe_obj <- geomx_save$spe_obj
geomx_sample_metadata <- geomx_save$sample_metadata
geomx_top5 <- "seurat"

geomx_spe_obj$shared_EBV <- geomx_spe_obj$EBV_Indicator
geomx_spe_obj$shared_EBV <- dplyr::recode(
  geomx_spe_obj$shared_EBV,
  "yes" = "EBV+",
  "no" = "EBV-"
)

geomx_spe_obj$shared_site <- geomx_spe_obj$Cohort
geomx_spe_obj$shared_tech <- "GeoMX"
geomx_spe_obj$shared_FOV <- geomx_spe_obj$ROI

geomx_spe_obj$shared_CT <- geomx_spe_obj$MergedLabel
geomx_spe_obj$shared_CT <- dplyr::recode(
  geomx_spe_obj$shared_CT,
  "Tumor" = "Tumor",
  "CD4T" = "Tcell",
  "CD8T" = "Tcell",
  "Macro" = "Macrophage",
  .default = NA_character_
)

geomx_keep <- !is.na(geomx_spe_obj$shared_CT)
geomx_spe_obj <- geomx_spe_obj[, geomx_keep]

geomx_cd <- as.data.frame(colData(geomx_spe_obj))

geomx_group <- interaction(
  geomx_cd$shared_EBV,
  geomx_cd$shared_site,
  geomx_cd$shared_tech,
  geomx_cd$shared_FOV,
  geomx_cd$shared_CT,
  drop = TRUE
)

stopifnot(all(tapply(geomx_cd$shared_EBV, geomx_group, function(x) length(unique(x))) == 1))
stopifnot(all(tapply(geomx_cd$shared_site, geomx_group, function(x) length(unique(x))) == 1))
stopifnot(all(tapply(geomx_cd$shared_tech, geomx_group, function(x) length(unique(x))) == 1))

geomx_counts <- as.matrix(assay(geomx_spe_obj, "counts"))
geomx_logcounts <- as.matrix(assay(geomx_spe_obj, "logcounts"))

geomx_pooled_counts <- t(rowsum(t(geomx_counts), geomx_group))

geomx_sum_logcounts <- t(rowsum(t(geomx_logcounts), geomx_group))
geomx_group_sizes <- as.numeric(table(geomx_group))
names(geomx_group_sizes) <- names(table(geomx_group))
geomx_group_sizes <- geomx_group_sizes[colnames(geomx_pooled_counts)]
geomx_pooled_logcounts <- sweep(geomx_sum_logcounts, 2, geomx_group_sizes, "/")

geomx_grp_levels <- colnames(geomx_pooled_counts)

g_EBV <- tapply(geomx_cd$shared_EBV, geomx_group, function(x) x[1])
g_site <- tapply(geomx_cd$shared_site, geomx_group, function(x) x[1])
g_tech <- tapply(geomx_cd$shared_tech, geomx_group, function(x) x[1])
g_FOV <- tapply(geomx_cd$shared_FOV, geomx_group, function(x) x[1])
g_CT <- tapply(geomx_cd$shared_CT, geomx_group, function(x) x[1])

g_ruv <- tapply(geomx_spe_obj$ruv_W1, geomx_group, mean, na.rm = TRUE)
# g_cells <- tapply(geomx_cd$cell_count, geomx_group, sum, na.rm = TRUE)

geomx_meta <- data.frame(
  shared_EBV = as.character(g_EBV[geomx_grp_levels]),
  shared_site = as.character(g_site[geomx_grp_levels]),
  shared_tech = as.character(g_tech[geomx_grp_levels]),
  shared_FOV = as.character(g_FOV[geomx_grp_levels]),
  shared_CT = as.character(g_CT[geomx_grp_levels]),
  shared_ruv_W1 = as.numeric(g_ruv[geomx_grp_levels]),
 # shared_cell_counts = as.numeric(g_cells[geomx_grp_levels]),
  shared_librarysize = as.numeric(colSums(geomx_pooled_counts)),
  row.names = geomx_grp_levels,
  stringsAsFactors = FALSE
)

geomx_spe_obj <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = geomx_pooled_counts,
                logcounts = geomx_pooled_logcounts),
  colData = S4Vectors::DataFrame(geomx_meta)
)

cat("GeoMX object prepared:\n")
cat(" - Samples:", ncol(geomx_spe_obj), "\n")
cat(" - Sites:", paste(unique(geomx_spe_obj$shared_site), collapse = ", "), "\n")
cat(" - Cell types:", paste(unique(geomx_spe_obj$shared_CT), collapse = ", "), "\n")
cat(" - EBV status:\n")
print(table(geomx_spe_obj$shared_EBV))

## ---------------------------------------------------------------------------
## Combine New CosMX and GeoMX datasets
## ---------------------------------------------------------------------------

cat("\n=== Combining New CosMX and GeoMX datasets ===\n")

# Find common genes between the two cohorts
common_genes <- intersect(rownames(cosmx_new_spe_obj), rownames(geomx_spe_obj))
cat("Common genes:", length(common_genes), "\n")

# Subset both to common genes
cosmx_new_spe_obj <- cosmx_new_spe_obj[common_genes, ]
geomx_spe_obj <- geomx_spe_obj[common_genes, ]

# Extract logcounts (batch-corrected expression)
cosmx_new_logcounts <- as.matrix(assay(cosmx_new_spe_obj, "logcounts"))
geomx_logcounts <- as.matrix(assay(geomx_spe_obj, "logcounts"))

# Combine expression matrices
full_expr <- cbind(cosmx_new_logcounts, geomx_logcounts)
cat("Combined expression matrix:", nrow(full_expr), "genes x", ncol(full_expr), "samples\n")

# Combine metadata
meta_cols <- c("shared_EBV", "shared_site", "shared_tech", "shared_FOV", "shared_CT",
               "shared_ruv_W1", "shared_librarysize")

cosmx_new_meta <- as.data.frame(colData(cosmx_new_spe_obj))[, meta_cols]
geomx_meta_use <- as.data.frame(colData(geomx_spe_obj))[, meta_cols]

combined_meta <- rbind(cosmx_new_meta, geomx_meta_use)
rownames(combined_meta) <- colnames(full_expr)

cat("Combined metadata:", nrow(combined_meta), "samples\n")
cat("Metadata summary:\n")
cat(" - Technologies:", paste(unique(combined_meta$shared_tech), collapse = ", "), "\n")
cat(" - Sites:", paste(unique(combined_meta$shared_site), collapse = ", "), "\n")
cat(" - Cell types:", paste(unique(combined_meta$shared_CT), collapse = ", "), "\n")
cat(" - EBV status:\n")
print(table(combined_meta$shared_EBV))

# Summary table by technology
cat("\nSamples by technology and EBV status:\n")
print(table(combined_meta$shared_tech, combined_meta$shared_EBV))

cat("\nSamples by technology and site:\n")
print(table(combined_meta$shared_tech, combined_meta$shared_site))

cat("\nSamples by technology and cell type:\n")
print(table(combined_meta$shared_tech, combined_meta$shared_CT))

## ---------------------------------------------------------------------------
## Create annotation bars for heatmaps
## ---------------------------------------------------------------------------

create_simple_annotation <- function(metadata,
                                     barplot_vars = c("total_counts"),
                                     barplot_height = 1) {
  ann_list <- list()
  col_maps <- list()

  ## Categorical annotations --------------------------------------------------

  # EBV status (force order: EBV+ then EBV- then unknown)
  if ("shared_EBV" %in% colnames(metadata)) {
    ebv_vals <- as.character(metadata$shared_EBV)
    ebv_levels <- c("EBV+", "EBV-", "unknown")
    # keep only levels that are present
    ebv_levels <- ebv_levels[ebv_levels %in% unique(ebv_vals)]
    ann_list$EBV <- factor(ebv_vals, levels = ebv_levels)
    col_maps$EBV <- c(
      "EBV+"   = "#d95f02",
      "EBV-"   = "#1b9e77",
      "unknown" = "grey80"
    )[ebv_levels]
  }

  # Cell type (Tumor / Tcell / Macrophage)
  if ("shared_CT" %in% colnames(metadata)) {
    ct_vals <- as.character(metadata$shared_CT)
    ann_list$Cell_Type <- ct_vals
    ct_palette <- c(
      "Tumor"      = "#7aadff",
      "Tcell"      = "#7dc79b",
      "Macrophage" = "#cfa7cd",
      "Others" = "grey80",
      "Stromal" = "grey80",
      "Treg" = "grey80",
      "pooled"     = "grey80"
    )
    ct_levels <- intersect(names(ct_palette), unique(ct_vals))
    col_maps$Cell_Type <- ct_palette[ct_levels]
  }

  # TMA / site information
  if ("shared_site" %in% colnames(metadata)) {
    site_palette <- c(
      "DFCI"           = "#e7298a",
      "Rochester"      = "#66a61e",
      "pooled"         = "grey80"
    )

    site_levels <- intersect(names(site_palette),
                             unique(as.character(metadata$shared_site)))
    metadata$shared_site <- factor(metadata$shared_site, levels = site_levels)

    ann_list$Site <- metadata$shared_site
    col_maps$Site <- site_palette[site_levels]
  }

  # Technology (CosMX_New / GeoMX)
  if ("shared_tech" %in% colnames(metadata)) {
    tech_palette <- c(
      "CosMX_New"  = "#377eb8",
      "GeoMX"      = "#4daf4a",
      "pooled"     = "grey80"
    )
    tech_levels <- intersect(names(tech_palette),
                             unique(as.character(metadata$shared_tech)))
    metadata$shared_tech <- factor(metadata$shared_tech, levels = tech_levels)

    ann_list$Tech <- metadata$shared_tech
    col_maps$Tech <- tech_palette[tech_levels]
  }

  # FOV annotation (kept simple; pooled FOV will appear as grey)
  if ("shared_FOV" %in% colnames(metadata)) {
    fov_vals <- as.character(metadata$shared_FOV)
    fov_levels <- sort(unique(fov_vals))
    fov_palette <- setNames(rep("#bdbdbd", length(fov_levels)), fov_levels)
    if ("pooled" %in% fov_levels) {
      fov_palette["pooled"] <- "grey80"
    }
    metadata$shared_FOV <- factor(fov_vals, levels = fov_levels)

    ann_list$FOV <- metadata$shared_FOV
    col_maps$FOV <- fov_palette
  }

  ## Barplot annotations (total_counts, cell_count) ---------------------------

  barplot_colors <- c(
    total_counts = "#e78ac3",
    cell_count   = "#fc8d62"
  )

  for (var in barplot_vars) {
    if (var %in% colnames(metadata)) {
      values <- metadata[[var]]
      if (!all(is.na(values))) {
        var_name <- gsub("_", " ", var)
        var_name <- paste0(toupper(substring(var_name, 1, 1)), substring(var_name, 2))
        var_name <- gsub(" ", "_", var_name)

        ann_list[[var_name]] <- anno_barplot(
          values,
          gp = gpar(fill = barplot_colors[var]),
          axis_param = list(gp = gpar(fontsize = 8)),
          height = unit(barplot_height, "cm")
        )
      }
    }
  }

  do.call(HeatmapAnnotation, c(ann_list, list(col = col_maps)))
}
## ---------------------------------------------------------------------------
## Load marker gene lists (for lineage marker GSVA) – borrowed from CosMX4TMA_1_3
## ---------------------------------------------------------------------------

marker_genes_list <- qread(marker_genes_path <- file.path(base_root,"Data","marker_genes_list.qs"))

CD4T_markers         <- marker_genes_list$CD4T_markers
CD8T_markers         <- marker_genes_list$CD8T_markers
B_cellmarkers        <- marker_genes_list$B_cellmarkers
Macrophage_marker_list <- marker_genes_list$Macrophage_marker_list
TAM_gene_clusters    <- marker_genes_list$TAM_gene_clusters
CD4T_panelgene       <- marker_genes_list$CD4T_panelgene
CD8T_panelgene       <- marker_genes_list$CD8T_panelgene
Tumor_genes_list     <- marker_genes_list$Tumor_genes_list

# Overlapped DEG gene lists (line 347–358 in CosMX4TMA_1_3_EBV_celltype_correlation.R)
Overlapped_DEG_Macrophage <- c(
  "LAMP3", "C1QA", "C1QC", "APOC1", "C1QB", "CCL18", "CD163", "CTSB",
  "GPNMB", "MMP9", "MRC1", "NR1H3", "THBS1", "HLA.DRA", "IL4R", "MMP14"
)

Overlapped_DEG_Bcell <- c(
  "REL", "SSR4", "HSPA6", "BTG2", "CD79A", "GSTM4", "NFKBIA", "PPIB",
  "UBE2S", "EGR1", "LIMD2", "SPIB", "EI24", "SSR2", "IKZF1", "SELL",
  "TNFSF9", "SPCS2", "TRAF4", "LRPAP1", "SPINT2", "CENPM", "UPF2",
  "UBALD2", "CD19", "PAX5", "IRF4", "CR2", "FCRL2", "FCRL3"
)

Overlapped_DEG_CD8T <- c(
  "FOXO1", "NKG7", "LAG3", "CD3D", "CD3E", "PAG1", "GZMA", "GZMK",
  "IL2RG", "CXCR3", "CDC26", "PSMD3", "VCP", "MAPK14", "TRIM14",
  "ALDH16A1", "BCL11B", "BAG4", "CD8A", "CD8B"
)

## ============================================================================
## EBV gene list (shared across datasets) – used for all heatmaps
## ============================================================================

cat("=== EBV genes shared across CosMX1 and GeoMX ===\n")

expr_ebv_cosmx1 <- assay(cosmx_new_spe_obj,  "logcounts")
expr_ebv_geomx  <- assay(geomx_spe_obj,  "logcounts")

virus_genes_cosmx1 <- c(
  grep("^LMP",   rownames(expr_ebv_cosmx1), value = TRUE),
  grep("^BCRF1", rownames(expr_ebv_cosmx1), value = TRUE),
  grep("^BNLF",  rownames(expr_ebv_cosmx1), value = TRUE),
  grep("^BLLF1", rownames(expr_ebv_cosmx1), value = TRUE),
  grep("^BZL",   rownames(expr_ebv_cosmx1), value = TRUE),
  grep("^EBNA",  rownames(expr_ebv_cosmx1), value = TRUE),
  grep("^RPM",   rownames(expr_ebv_cosmx1), value = TRUE)
)

virus_genes_geomx <- c(
  grep("^LMP",   rownames(expr_ebv_geomx), value = TRUE),
  grep("^BCRF1", rownames(expr_ebv_geomx), value = TRUE),
  grep("^BNLF",  rownames(expr_ebv_geomx), value = TRUE),
  grep("^BLLF1", rownames(expr_ebv_geomx), value = TRUE),
  grep("^BZL",   rownames(expr_ebv_geomx), value = TRUE),
  grep("^EBNA",  rownames(expr_ebv_geomx), value = TRUE),
  grep("^RPM",   rownames(expr_ebv_geomx), value = TRUE)
)

virus_genes_cosmx1 <- unique(virus_genes_cosmx1)
virus_genes_geomx  <- unique(virus_genes_geomx)

virus_genes_overlap <- sort(
  Reduce(
    intersect,
    list(virus_genes_cosmx1, virus_genes_geomx)
  )
)

cat("Shared EBV genes across all 3 datasets:", length(virus_genes_overlap), "\n")
if (length(virus_genes_overlap) > 0) {
  cat("EBV genes:", paste(virus_genes_overlap, collapse = ", "), "\n")
}

ebv_gene_candidates <- virus_genes_overlap

## ---------------------------------------------------------------------------
## Compute combined GSVA scores (functional + lineage + EBV)
## ---------------------------------------------------------------------------

compute_combined_scores <- function(expr_mat) {
  dataset_genes <- rownames(expr_mat)

  # Information 1: functional GSVA (Macro / Tumor / CD4T)
  Macro_pathways <- functionalenrichment[grep("^Macro", names(functionalenrichment))]
  Tumor_pathways <- functionalenrichment[grep("^Tumor", names(functionalenrichment))]
  CD4T_pathways  <- functionalenrichment[grep("^CD4T",  names(functionalenrichment))]

  Macro_pathways_renamed <- rename_pathways_with_counts(Macro_pathways, dataset_genes)
  Tumor_pathways_renamed <- rename_pathways_with_counts(Tumor_pathways, dataset_genes)
  CD4T_pathways_renamed  <- rename_pathways_with_counts(CD4T_pathways,  dataset_genes)

  function_gene_sets <- c(
    Macro_pathways_renamed,
    Tumor_pathways_renamed,
    CD4T_pathways_renamed
  )
  function_gene_sets <- function_gene_sets[sapply(function_gene_sets, length) > 0]

  gsva_results_function <- NULL
  if (length(function_gene_sets) > 0) {
    SetGSVAPar_function <- gsvaParam(
      exprData = expr_mat,
      geneSets = function_gene_sets,
      minSize  = 3,
      maxSize  = 400,
      kcdf     = "Gaussian"
    )
    gsva_results_function <- gsva(SetGSVAPar_function, verbose = TRUE)
  }

  CD4T_marker_tonsil <- c("TCF7","FKBP5","AC004585.1","ITM2A","PASK","CD40LG","THEMIS","ID3","NFIA","PDCD1","PTPN13","NUCB2",
                          "TOX2","COTL1","DRAIC","CXCR4","SORL1","NR3C1","ASAP1","KSR2","SEMA4D",
                          "TMEM123","LEF1","AFF3","TRAT1","ST8SIA1","IL21","ZFP36L2",
                          "NIN","SCGB3A1","ACTN1","ARMH1","GNG4","STK17A","SCML4","TC2N",
                          "TNFSF8","CDK5R1","IFITM1","CD226","PPP1CC","ATM","CCDC50","LRMP",
                          "RBMS1","NCOA7","AC068587.4","CPM","H2AFZ","DHRS7","PRRC2B","BTLA",
                          "IL6ST","PLAC8","PDE3B","AC012645.3","TRERF1","AHI1","MCUB","P2RY8",
                          "CD84","TRIM8","CXCR5","SMCO4","CORO1B","MAML2","CCND3","ANK3","ITPKB",
                          "LINC01934","CXCL13","CAMK4","AL450352.1","CEP128","TOX")
  CD8T_marker_tonsil <- c(
    "CCL5","GZMK","NKG7","GZMA","CST7","CCL4","KLRK1","CD8A","KLRG1",
    "SAMD3","TRGC2","CLDND1","CTSW","CMC1","APOBEC3G","KLRD1","EOMES",
    "RUNX3","SLAMF7","PLEK","LYST","GNLY","HLA.DPB1","CD8B","AOAH","GZMH",
    "ITGA1","ITM2C","PRF1","CRTAM","LYAR","IRF1","IFNG.AS1","AH K",
    "HLA.DPA1","CCL4L2","HLA.DRB1","CXCR6","CXCR3","MYO1F","CLEC2B","FCRL3",
    "GPR171","HCST","ANXA1","SLF1","CD63","CD99","ARHGAP26","CD74","PLAAT4",
    "STOM","DTHD1","A2M.AS1","TRDC","ZEB2","LINC01871","RNF213","XCL2",
    "PTGDR","MATK","A2M","ANXA2","PPP2R2B","APMAP","TRGC1","XCL1","IL32",
    "F2R","PRR5L","CCL3L1","CD160","GZMB","MYBL1","PLAAT3","PECAM1","CCR5",
    "GZMM","CHST12","PDCD4","MIAT","NIBAN1","PSMB9","CYBA","LCP1","ADGRE5",
    "MBP","PARP8","CD96","PILRB","DUSP2","IL10RA","TRG.AS1","GSTP1","ITGA4",
    "VIM","CLIC1","SH2D1A","HLA.F","C5orf56","CYTOR","S100A6","PREX1","LY6E",
    "ACTN4","BTN3A1","PIP4K2A","TAPBP","BTN3A2","D JC1","EMP3","ARPC5L",
    "APOBEC3C","STAT4","SLFN12L","PTPN22","CD81","GUK1","PYHIN1","MT2A",
    "CTSC","GPR174","PIK3R1","NELL2","HLA-DRA","S100A10"
  )


  # Information 2: lineage marker GSVA
  gene_sets_lineage <- list(
    CD4T_panelgene_pancancerT = intersect(unique(unlist(CD4T_panelgene)),         dataset_genes),
    CD4T_marker_tonsil = intersect(unique(CD4T_marker_tonsil),         dataset_genes),
    CD4T_markers              = intersect(CD4T_markers,                           dataset_genes),
    CD8T_panelgene_pancancerT = intersect(unique(unlist(CD8T_panelgene)),         dataset_genes),
    CD8T_marker_tonsil = intersect(unique(CD8T_marker_tonsil),         dataset_genes),
    CD8T_markers              = intersect(CD8T_markers,                           dataset_genes),
    TAM_gene_clusters         = intersect(unique(unlist(TAM_gene_clusters)),      dataset_genes),
    Macrophage_markers        = intersect(unique(unlist(Macrophage_marker_list)), dataset_genes),
    Tumor_genes               = intersect(unique(unlist(Tumor_genes_list)),       dataset_genes),
    B_cellmarkers             = intersect(B_cellmarkers,                          dataset_genes),
    Overlapped_DEG_Macrophage = intersect(Overlapped_DEG_Macrophage,             dataset_genes),
    Overlapped_DEG_Bcell      = intersect(Overlapped_DEG_Bcell,                  dataset_genes),
    Overlapped_DEG_CD8T       = intersect(Overlapped_DEG_CD8T,                   dataset_genes)
  )
  gene_sets_lineage <- gene_sets_lineage[sapply(gene_sets_lineage, length) > 0]

  gsva_results_lineage <- NULL
  if (length(gene_sets_lineage) > 0) {
    SetGSVAPar_lineage <- gsvaParam(
      exprData = expr_mat,
      geneSets = gene_sets_lineage,
      minSize  = 1,
      maxSize  = 400,
      kcdf     = "Gaussian"
    )
    gsva_results_lineage <- gsva(SetGSVAPar_lineage, verbose = TRUE)
  }

  score_mats <- list()
  if (!is.null(gsva_results_function)) score_mats$gsva_function <- gsva_results_function
  if (!is.null(gsva_results_lineage))  score_mats$gsva_lineage  <- gsva_results_lineage

  # Information 3: EBV gene block (GSVA enrichment)
  if (!is.null(ebv_gene_candidates) && length(ebv_gene_candidates) > 0) {
    ebv_genes_use <- intersect(ebv_gene_candidates, dataset_genes)
    if (length(ebv_genes_use) > 0) {
      ebv_gene_sets <- list(EBV_genes = ebv_genes_use)
      SetGSVAPar_ebv <- gsvaParam(
        exprData = expr_mat,
        geneSets = ebv_gene_sets,
        minSize  = 1,
        maxSize  = max(1, length(ebv_genes_use)),
        kcdf     = "Gaussian"
      )
      ebv_gsva <- gsva(SetGSVAPar_ebv, verbose = TRUE)
      # Single EBV GSVA score row
      rownames(ebv_gsva) <- "EBV_GSVA"
      score_mats$EBV <- ebv_gsva
    }
  }

  do.call(rbind, score_mats)
}

## ============================================================================
## HELPER: color scale for a score matrix
## ============================================================================

build_color_fun <- function(score_mat) {
  score_vals <- as.numeric(score_mat)
  score_vals <- score_vals[is.finite(score_vals)]

  if (length(score_vals) > 0) {
    score_min <- quantile(score_vals, 0.02, na.rm = TRUE)
    score_mid <- median(score_vals, na.rm = TRUE)
    score_max <- quantile(score_vals, 0.98, na.rm = TRUE)
  } else {
    score_min <- -2
    score_mid <- 0
    score_max <- 2
  }

  colorRamp2(
    c(score_min, score_mid, score_max),
    c("#b953a0", "black", "#f4ed17")
  )
}

## ============================================================================
## HELPER: pooling by grouping variables
## ============================================================================

pool_by_group <- function(expr_mat, meta, group_cols) {
  group_df <- meta[, group_cols, drop = FALSE]
  group <- interaction(group_df, drop = TRUE)
  levels_g <- levels(group)

  # pooled expression: mean across original samples in each group
  pooled_expr <- sapply(levels_g, function(g) {
    idx <- which(group == g)
    if (length(idx) == 1L) {
      as.numeric(expr_mat[, idx, drop = FALSE])
    } else {
      rowMeans(expr_mat[, idx, drop = FALSE])
    }
  })
  if (is.vector(pooled_expr)) {
    pooled_expr <- matrix(
      pooled_expr,
      nrow = nrow(expr_mat),
      dimnames = list(rownames(expr_mat), levels_g)
    )
  } else {
    rownames(pooled_expr) <- rownames(expr_mat)
    colnames(pooled_expr) <- levels_g
  }

  # pooled metadata: first row for factors, sum counts/librarysize
  pooled_meta <- do.call(rbind, lapply(levels_g, function(g) {
    idx <- which(group == g)
    subm <- meta[idx, , drop = FALSE]
    out <- subm[1, , drop = FALSE]

    if ("shared_librarysize" %in% colnames(subm)) {
      out$shared_librarysize <- sum(subm$shared_librarysize, na.rm = TRUE)
    }
    if ("shared_cell_counts" %in% colnames(subm)) {
      out$shared_cell_counts <- sum(subm$shared_cell_counts, na.rm = TRUE)
    }

    rownames(out) <- g
    out
  }))

  pooled_meta$total_counts <- pooled_meta$shared_librarysize
  pooled_meta$cell_count <- pooled_meta$shared_cell_counts

  list(expr = pooled_expr, meta = pooled_meta)
}

## ============================================================================
## Heatmap 2: EBV+/− x technology x site, pooled across FOV/CT
## ============================================================================

cat("\n=== Heatmap 2: pooled by shared_tech, shared_site, shared_EBV (EBV last) ===\n")



pool2 <- pool_by_group(
  expr_mat = full_expr,
  meta = combined_meta,
  group_cols = c("shared_EBV", "shared_tech", "shared_site")
)

expr_h2 <- pool2$expr
meta_h2 <- pool2$meta

meta_h2$shared_FOV <- "pooled"
meta_h2$shared_CT <- "pooled"

scores_h2 <- compute_combined_scores(expr_h2)
col_fun_h2 <- build_color_fun(scores_h2)
ann_h2 <- create_simple_annotation(
  metadata = meta_h2,
  barplot_vars = c("total_counts", "cell_count"),
  barplot_height = 1
)

ht_h2 <- Heatmap(
  scores_h2,
  name = "Score",
  col = col_fun_h2,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = ann_h2,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_title = "Heatmap 2: EBV+/− x technology x site (pooled FOV/CT)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Score",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  )
)

pdf(file.path(results_dir,"Supplementary Figure.9F.pdf"), width = 16, height = 10)
draw(ht_h2, heatmap_legend_side = "right", show_annotation_legend = TRUE)
dev.off()
