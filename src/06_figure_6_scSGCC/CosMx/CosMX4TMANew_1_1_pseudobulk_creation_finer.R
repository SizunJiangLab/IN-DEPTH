# New CosMX Pseudobulk Data Creation for GeoMX Reproduction (Finer Cell Type)
# Script: New_CosMX4TMA_1_1_pseudobulk_creation_finer.R
# Create pseudobulk expression data from New_CosMX data using FINER cell type annotations

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(edgeR)
library(readxl)
library(SpatialExperiment)
base_root <- "/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/CosMXNew_code_data_publish/"
standr_path <- file.path(base_root, "standR_covariate")
output_base <- file.path(base_root, "Data/CosMXNew/Exampleoutpu")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

# Download standR_covariate package to local directory if not present
if (!dir.exists(standr_path)) {
  cat("Downloading standR_covariate from GitHub to local directory...\n")
  system(paste0("git clone https://github.com/BMEngineeR/standR_covariate.git ", standr_path))
}

# Load the package from local directory
devtools::load_all(standr_path)
library(qs)
library(networkD3)

# Load modular functions
source(file.path(base_root, "Code/src/6-Pathway_validation_DLBCL_preload_function-sc.R"))
source(file.path(base_root, "Code/src/utils.R"))

# Set paths
data_path <- file.path(base_root,"Data/New_CosMX")
results_path <- file.path(output_base, "CosMX4TMANew_1_1_pseudobulk_creation_finer")

# Create results directory if it doesn't exist
if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
}

setwd(data_path)

# Load data
cat("Loading New CosMX data...\n")
cat("NOTE: This script uses FINER cell type annotations for pseudobulk creation\n\n")

# 1. Load and standardize metadata from core_meta.csv
cat("Step 1: Loading and filtering core metadata...\n")
core_meta <- read.csv("core_meta.csv", stringsAsFactors = FALSE)

# Filter core_meta based on specified criteria
# Remove: TMAmap == NA, CosMXprofile == No or "No, but need to keep for protein", ExcludedOnBoard == Yes
core_meta_filtered <- core_meta %>%
  filter(!is.na(TMAmap)) %>%
  filter(CosMXProfile != "No") %>%
  filter(CosMXProfile != "No, but need to keep for protein") %>%
  filter(ExcludedOnBoard != "Yes") %>%
  # Standardize coreName (handle case inconsistency like C06 vs c06)
  # mutate(coreName = tolower(coreName)) %>%
  # Add EBV_status column based on TMAmap
  mutate(EBV_status = ifelse(grepl("EBV\\+", TMAmap, ignore.case = TRUE), "EBV+",
                              ifelse(grepl("EBV-", TMAmap, ignore.case = TRUE), "EBV-", "Unknown"))) %>%
  # Create TMA_core identifier for merging
  mutate(TMA_core = paste0(TMA, "_", coreName))

cat("Filtered core metadata:\n")
print(table(core_meta_filtered$EBV_status,core_meta_filtered$TMA))
cat("Total cores to process:", nrow(core_meta_filtered), "\n")

# 2. Load coarse annotation data
cat("\nStep 2: Loading coarse annotation data...\n")
mymeta_lineage <- read.csv("corase_annotation.csv", stringsAsFactors = FALSE)

# Standardize coreName in annotation data
mymeta_lineage <- mymeta_lineage %>%
  mutate(coreName = coreName) %>%
  mutate(cropName = coreName) %>%
  mutate(TMA_core = paste0(TMA, "_", coreName))

cat("Annotation data dimensions:", dim(mymeta_lineage), "\n")
cat("Cell types in annotation:\n")
print(table(mymeta_lineage$Annotation))

# 3. Merge core_meta with annotation data to create cell metadata
cat("\nStep 3: Merging core metadata with annotation data...\n")
cell_meta <- mymeta_lineage %>%
  right_join(core_meta_filtered %>%
               select(TMA_core, TMA, coreName, EBV_status, TMAmap, CosMXProfile, runID),
             by = c("TMA_core","TMA","coreName"))

cat("Cell metadata after merging:\n")
cat("Total cells:", nrow(cell_meta), "\n")
cat("Cores represented:", length(unique(cell_meta$TMA_core)), "\n")
cat("EBV status distribution:\n")
print(table(cell_meta$EBV_status))

# 3.1 Load finer cell type annotation and subset cell_meta
cat("\nStep 3.1: Loading finer cell type annotation...\n")
mymeta_finer <- read.csv("FINAL_DF_1219.csv")
mymeta_finer2 <- read.csv("FINAL_DF_with_scaled_fxn.csv")
mymeta_finer_backup <- mymeta_finer2 %>%
  mutate(unique_id = paste0(TMA,"_",cropName,"_",cellLabel))

cat("Finer annotation dimensions:", dim(mymeta_finer_backup), "\n")
cat("Unique cells in finer annotation:", length(unique(mymeta_finer_backup$unique_id)), "\n")
cat("Cell types in finer annotation:\n")
print(table(mymeta_finer_backup$Annotation))

# Rename Annotation column in cell_meta to avoid conflict
# mymeta_finer_backup has more detailed cell types
cat("\nRenaming cell_meta$Annotation to Annotation_coarse...\n")
cell_meta <- cell_meta %>%
  dplyr::rename(Annotation_coarse = Annotation)

# Store original cell_meta before subsetting to track dropped cells
cell_meta_before_subset <- cell_meta
cells_before_subset <- nrow(cell_meta_before_subset)
cat("\nCells in cell_meta before subsetting:", cells_before_subset, "\n")

# Find intersection
overlapping_cells <- intersect(mymeta_finer_backup$unique_id, cell_meta$unique_id)
cat("Overlapping cells:", length(overlapping_cells), "\n")

# Identify dropped cells before subsetting
dropped_cells <- cell_meta_before_subset %>%
  filter(!unique_id %in% mymeta_finer_backup$unique_id) %>%
  select(unique_id, Annotation_coarse, TMA_core, EBV_status)

# Subset cell_meta
cell_meta <- cell_meta %>%
  filter(unique_id %in% mymeta_finer_backup$unique_id)

cells_after_subset <- nrow(cell_meta)
cells_dropped <- cells_before_subset - cells_after_subset

cat("\n=== Cell Subsetting Summary ===\n")
cat("Cells before subsetting:", cells_before_subset, "\n")
cat("Cells after subsetting:", cells_after_subset, "\n")
cat("Cells dropped:", cells_dropped, "\n")
cat("Percentage retained:", round(cells_after_subset / cells_before_subset * 100, 2), "%\n")

# Show dropped cells information
if (nrow(dropped_cells) > 0) {
  cat("\n=== Dropped Cells Information ===\n")
  cat("Total dropped cells:", nrow(dropped_cells), "\n")
  cat("Dropped cells by Annotation_coarse:\n")
  print(table(dropped_cells$Annotation_coarse))
  cat("\nDropped cells by EBV status:\n")
  print(table(dropped_cells$EBV_status))
  cat("\nFirst 10 dropped cell IDs:\n")
  print(head(dropped_cells, 10))
}

# Add finer cell type annotation to cell_meta from mymeta_finer_backup
cat("\n=== Adding Finer Cell Type Annotation ===\n")
cell_meta <- cell_meta %>%
  left_join(mymeta_finer_backup %>% select(unique_id, Annotation,LMP1,LAG3,C1q,FoxP3,Myc),#select(unique_id, Annotation,LMP1,LAG3,C1q,COL1A1,FoxP3,BCL6),
            by = "unique_id",
            suffix = c("", "_finer"))

# Rename the new Annotation column to Annotation_finer for clarity
if ("Annotation" %in% colnames(cell_meta)) {
  cell_meta <- cell_meta %>%
    dplyr::rename(Annotation_finer = Annotation)
}

cat("Finer annotation added successfully!\n")
cat("Cell types in Annotation_finer:\n")
print(table(cell_meta$Annotation_finer))

# Show comparison of coarse vs finer annotations
cat("\n=== Annotation Comparison (Coarse vs Finer) ===\n")
annotation_comparison <- cell_meta %>%
  group_by(Annotation_coarse, Annotation_finer) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Annotation_coarse, desc(count))
print(annotation_comparison)

cat("\nRetained cells by Annotation_coarse:\n")
print(table(cell_meta$Annotation_coarse))

# 4. Read expression matrices from 06_tx_inframe folder
cat("\nStep 4: Reading expression matrices from 06_tx_inframe...\n")

# Get list of available expression files
expr_files <- list.files(file.path(data_path, "06_tx_inframe"),
                         pattern = "\\.csv\\.gz$",
                         full.names = TRUE)

cat("Found", length(expr_files), "expression files\n")

# Extract TMA and coreName from filenames
expr_file_info <- data.frame(
  file_path = expr_files,
  file_name = basename(expr_files)
) %>%
  mutate(
    file_base = gsub("\\.csv\\.gz$", "", file_name),
    TMA = sub("_.*", "", file_base),
    coreName = sub(".*_", "", file_base),
    TMA_core = paste0(TMA, "_", coreName)
  )

# Filter expression files based on cores in filtered metadata
expr_files_to_load <- expr_file_info %>%
  filter(TMA_core %in% core_meta_filtered$TMA_core)

cat("Expression files to load based on core_meta filtering:", nrow(expr_files_to_load), "\n")
cat("Cores to load:\n")
print(expr_files_to_load$TMA_core)

# Load and combine expression matrices
expr_list <- list()
for (i in 1:nrow(expr_files_to_load)) {
  cat("Loading", expr_files_to_load$file_name[i], "...\n")
  expr_data <- read.csv(gzfile(expr_files_to_load$file_path[i]), stringsAsFactors = FALSE)
  expr_list[[i]] <- expr_data
}

# Combine all expression data
cat("\nCombining all expression matrices...\n")
mydata <- do.call(rbind, expr_list)

cat("Combined expression matrix dimensions:", dim(mydata), "\n")

# Create unique_id from data_tag, region_name, and cell_id
mydata <- mydata %>%
  mutate(unique_id = paste0(data_tag, "_", region_name, "_", cell_id))

# Check overlap with cell_meta
overlap_cells <- intersect(mydata$unique_id, cell_meta$unique_id)
cat("Overlapping cells between expression and metadata:", length(overlap_cells), "\n")

# Filter to keep only overlapping cells
mydata_filtered <- mydata %>%
  filter(unique_id %in% overlap_cells)

cell_meta_filtered <- cell_meta %>%
  filter(unique_id %in% overlap_cells)

# Set unique_id as rownames and remove metadata columns from expression matrix
rownames(mydata_filtered) <- mydata_filtered$unique_id
mydata_expr <- mydata_filtered %>%
  select(-data_tag, -region_name, -cell_id, -unique_id)

# Set unique_id as rownames for cell metadata
rownames(cell_meta_filtered) <- cell_meta_filtered$unique_id

cat("\nFinal data dimensions:\n")
cat("Expression matrix (cells x genes):", dim(mydata_expr), "\n")
cat("Cell metadata:", dim(cell_meta_filtered), "\n")
rownames(cell_meta_filtered) <- cell_meta_filtered$unique_id
# Save single-cell level SpatialExperiment object
cat("\n=== Creating and Saving Single-Cell SpatialExperiment Object ===\n")

# Transpose expression matrix to genes x cells for SpatialExperiment
expr_matrix_for_spe <- t(as.matrix(mydata_expr))
expr_matrix_for_spe <- expr_matrix_for_spe[,rownames(cell_meta_filtered)]
# Create SpatialExperiment object
fov_mapping_table <- read.csv(file.path(data_path, "06_tx_inframe","mapping_coor","cell_fov_map.csv.gz"))
fov_mapping_table$unique_id <- paste0(fov_mapping_table$data_tag,"_",fov_mapping_table$region_name,"_",fov_mapping_table$cell_id)
intersect_cells <- intersect(fov_mapping_table$unique_id,cell_meta_filtered$unique_id)
fov_mapping_table <- fov_mapping_table[match(intersect_cells,fov_mapping_table$unique_id),]
fov_mapping_table <- fov_mapping_table[match(cell_meta_filtered$unique_id,fov_mapping_table$unique_id),]
cell_meta_filtered_2 <- fov_mapping_table %>% select(unique_id,fov_id) %>%  right_join(cell_meta_filtered,by = "unique_id")
spe_single_cell <- SingleCellExperiment(
  assays = list(counts = expr_matrix_for_spe),
  colData = cell_meta_filtered_2
)

cat("Single-cell SpatialExperiment object created:\n")
cat("- Dimensions:", nrow(spe_single_cell), "genes x", ncol(spe_single_cell), "cells\n")
cat("- Metadata columns:", ncol(colData(spe_single_cell)), "\n")

# Save the single-cell SPE object
qsave(spe_single_cell, file = file.path(results_path, "CosMX4TMANew_SPE_object.qs"))
cat("\nSingle-cell SPE object saved to: CosMX4TMANew_SPE_object.qs\n")

# 5. Create ROI-based pseudobulk using FINER cell type annotations
cat("\n=== Creating ROI-based Pseudobulk Data with FINER Cell Types ===\n")

# Add ROI identifier using TMA_core
cell_meta_filtered$ROI_ID <- cell_meta_filtered$TMA_core

cat("Available ROIs (TMA_cores):\n")
print(table(cell_meta_filtered$ROI_ID))

# Define cell types of interest for pseudobulk creation using FINER annotations
unique_cell_types_finer <- unique(cell_meta_filtered$Annotation_finer)
cat("\nUnique FINER cell types found:\n")
print(unique_cell_types_finer)
cat("\nTotal unique finer cell types:", length(unique_cell_types_finer), "\n")

# Show cell count for each finer cell type
cat("\nCell count by finer cell type:\n")
finer_cell_type_counts <- table(cell_meta_filtered$Annotation_finer)
print(sort(finer_cell_type_counts, decreasing = TRUE))

# Filter out finer cell types with too few cells (e.g., < 10 cells total)
min_cells_threshold <- 1
cell_types_with_enough_cells <- names(finer_cell_type_counts[finer_cell_type_counts >= min_cells_threshold])

cat("\nFiner cell types with at least", min_cells_threshold, "cells:\n")
for (ct in cell_types_with_enough_cells) {
  count <- sum(cell_meta_filtered$Annotation_finer == ct, na.rm = TRUE)
  cat("-", ct, ":", count, "cells\n")
}

# Use all finer cell types with enough cells for pseudobulk creation
cell_types_of_interest_finer <- cell_types_with_enough_cells

cat("\nUsing", length(cell_types_of_interest_finer), "finer cell types for pseudobulk creation\n")

# Create modified metadata for pseudobulk function
# The CreatePseudoBulk function expects 'window_id' column
cell_meta_for_pseudobulk <- cell_meta_filtered
cell_meta_for_pseudobulk$window_id <- cell_meta_for_pseudobulk$ROI_ID

# Filter for cell types of interest (finer annotations)
cell_meta_for_pb <- cell_meta_for_pseudobulk %>%
  filter(Annotation_finer %in% cell_types_of_interest_finer)

# Filter expression data accordingly
mydata_for_pb <- mydata_expr[cell_meta_for_pb$unique_id, ]

cat("\nFiltered data for pseudobulk:\n")
cat("- Metadata:", dim(cell_meta_for_pb), "\n")
cat("- Expression data:", dim(mydata_for_pb), "\n")

# Create pseudobulk expression data using the modular function with FINER annotations
cat("\nCreating pseudobulk expression data for all finer cell types...\n")
cell_meta_for_pb$Annotation_pooled <- cell_meta_for_pb$Annotation_finer
rownames(cell_meta_for_pb) <- cell_meta_for_pb$unique_id
cell_meta_for_pb$cell_id <- cell_meta_for_pb$unique_id
pseudobulk_result <- CreatePseudoBulk(
  mymeta = cell_meta_for_pb,
  mydata = mydata_for_pb,
  cell_types_for_pseudobulk = cell_types_of_interest_finer,
  aggregate_meta_cols = c("CD3", "Pax5","CD68","CD4","CD8","CD163","CD206","LMP1","LAG3","C1q","FoxP3","Myc"), #c("LMP1","CD3","Pax5","CD68","CD4","CD8","CD163","CD206","LAG3","C1q","COL1A1","FoxP3","BCL6"),
  normalize_by_cell_count = FALSE,
  min_cells_per_sample = 5,  # Require at least 5 cells per pseudobulk sample
  verbose = TRUE
)

# Extract results
pseudobulk_expr <- pseudobulk_result$pseudobulk_expr
sample_metadata <- pseudobulk_result$sample_metadata
cell_type_summary <- pseudobulk_result$cell_type_summary

# Remove negative and system control genes
pseudobulk_expr <- pseudobulk_expr[!grepl("Negative|system", rownames(pseudobulk_expr), ignore.case = TRUE), ]

cat("Pseudobulk creation completed:\n")
cat("- Expression matrix dimensions:", dim(pseudobulk_expr), "\n")
cat("- Sample metadata dimensions:", dim(sample_metadata), "\n")

# Add EBV status and core information to sample metadata
sample_metadata <- sample_metadata %>%
  left_join(core_meta_filtered %>%
              select(TMA_core, EBV_status, TMAmap, CosMXProfile, runID,TMA),
            by = c("window_id" = "TMA_core"))

cat("\nSample metadata with EBV status:\n")
print(head(sample_metadata))

# Create summary statistics
cat("\n=== Pseudobulk Summary Statistics (Finer Cell Types) ===\n")
cat("Samples per finer cell type:\n")
print(table(sample_metadata$cell_type))

cat("\nSamples per EBV status:\n")
print(table(sample_metadata$EBV_status))

cat("\nSamples per finer cell type and EBV status:\n")
print(table(sample_metadata$cell_type, sample_metadata$EBV_status))

# Filter out genes with very low expression
gene_sums <- rowSums(pseudobulk_expr)
expressed_genes <- names(gene_sums[gene_sums > 1])
pseudobulk_expr_filtered <- pseudobulk_expr[expressed_genes, ]

cat("\nGene filtering completed:\n")
cat("- Original genes:", nrow(pseudobulk_expr), "\n")
cat("- Filtered genes:", nrow(pseudobulk_expr_filtered), "\n")

# Save intermediate results
cat("\n=== Saving Results (Finer Cell Type Pseudobulk) ===\n")

# Save pseudobulk results with "_finer" suffix
qsave(pseudobulk_expr_filtered, file = file.path(results_path, "New_pseudobulk_expr_filtered_finer.qs"))
qsave(sample_metadata, file = file.path(results_path, "New_sample_metadata_finer.qs"))
qsave(cell_type_summary, file = file.path(results_path, "New_cell_type_summary_finer.qs"))

# Also save cell metadata for reference
qsave(cell_meta_filtered, file = file.path(results_path, "New_cell_meta_filtered_finer.qs"))

cat("\nResults saved to:", results_path, "\n")
cat("Files created:\n")
cat("- New_pseudobulk_expr_filtered_finer.qs\n")
cat("- New_sample_metadata_finer.qs\n")
cat("- New_cell_type_summary_finer.qs\n")
cat("- New_cell_meta_filtered_finer.qs\n")

cat("\n=== Script completed successfully ===\n")
cat("NOTE: This pseudobulk data was created using FINER cell type annotations\n")
