# SGCC and DEG Analysis for Cell-Cell Interactions (New CosMX TMA)
# Analysis: High SGCC vs Low SGCC differential expression, stratified by EBV status
# WITH CELL FILTERING: Only include windows where BOTH cell types have >= MIN_CELLS_PER_WINDOW

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("BioGSP", quietly = TRUE)) {
  devtools::install_github("BMEngineeR/BioGSP")
}else{
  library(BioGSP)
}

library(edgeR)
library(limma)
library(Matrix.utils)
library(ComplexHeatmap)
library(circlize)
library(qs)
library(enrichR)
library(GSVA)
library(gridExtra)
# Load batch correction functions
base_root <- "/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/CosMXNew_code_data_publish"
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
source(file.path(base_root, "/Code/src/utils.R"))

## ============================================================================
## CONFIGURATION
## ============================================================================
CELL_TYPE_1 <- "B cell"# CD8 T; CD4 T; B cell
CELL_TYPE_2 <- "Macrophage"
EXPRESSION_CELL_TYPE <- "Macrophage"

# Cell count threshold for filtering windows and pseudobulk
MIN_CELLS_PER_WINDOW <- 5

analysis_label <- paste0(gsub(" ", "_", CELL_TYPE_1), "_", gsub(" ", "_", CELL_TYPE_2))
expression_label <- gsub(" ", "_", EXPRESSION_CELL_TYPE)

cat("\n=== SGCC DEG Analysis Configuration (WITH CELL FILTERING) ===\n")
cat("Cell Type 1:", CELL_TYPE_1, "\n")
cat("Cell Type 2:", CELL_TYPE_2, "\n")
cat("Expression Cell Type:", EXPRESSION_CELL_TYPE, "\n")
cat("Min Cells Per Window (BOTH cell types):", MIN_CELLS_PER_WINDOW, "\n")
cat("Comparison: High SGCC (top 1/3) vs Low SGCC (bottom 1/3)\n")
cat("Stratification: EBV+ and EBV- separately\n\n")

# Set paths
data_path <- file.path(output_base, "CosMX4TMANew_1_1_pseudobulk_creation_finer")
base_results_path <- file.path(output_base, "CosMX4TMANew_3_SGCC_DEG_analysis_finer-filtercell")

# Create subfolder for this cell type pair
results_path <- file.path(base_results_path, analysis_label)

if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
  cat("Created results folder:", results_path, "\n")
}

setwd(data_path)

## ============================================================================
## LOAD DATA AND CREATE SLIDING WINDOWS
## ============================================================================
cat("\n=== Loading Data ===\n")
spe_single_cell <- qread(file.path(data_path, "CosMX4TMANew_SPE_object.qs"))
mymeta <- colData(spe_single_cell)
# mymeta$fov <- paste0(mymeta$TMA, "_", mymeta$coreName)
mymeta$fov <- paste0(mymeta$TMA,"_",mymeta$coreName,"_",mymeta$fov_id)
mymeta$unique_id_fov <- paste0(mymeta$TMA,"_",mymeta$coreName,"_","fov",mymeta$fov_id,"_",mymeta$cellLabel)
mymeta$cell_id <- mymeta$unique_id_fov
mymeta$Annotation_pooled <- mymeta$Annotation_finer

cat("\n=== Creating Sliding Windows ===\n")
window_result <- CreateSlidingWindow(
  mymeta = as.data.frame(mymeta),
  window_size =500,
  stride = 400,
  verbose = TRUE
)

mymeta <- window_result$mymeta_windowed
all_windows <- window_result$all_windows

## ============================================================================
## FILTER WINDOWS BY CELL COUNT
## ============================================================================
cat("\n=== Filtering Windows by Cell Count ===\n")
cat("Filtering windows where BOTH cell types have >=", MIN_CELLS_PER_WINDOW, "cells\n")

# Count cells of each type in each window
window_cell_counts <- mymeta %>%
  filter(Annotation_pooled %in% c(CELL_TYPE_1, CELL_TYPE_2)) %>%
  group_by(window_id, Annotation_pooled) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Annotation_pooled, values_from = cell_count, values_fill = 0)

# Filter to windows where both cell types have >= MIN_CELLS_PER_WINDOW
valid_windows <- window_cell_counts %>%
  filter(.data[[CELL_TYPE_1]] >= MIN_CELLS_PER_WINDOW & .data[[CELL_TYPE_2]] >= MIN_CELLS_PER_WINDOW)

cat("Total windows before filtering:", length(unique(mymeta$window_id)), "\n")
cat("Windows passing filter (both cell types >=", MIN_CELLS_PER_WINDOW, "):", nrow(valid_windows), "\n")
cat("Windows removed:", length(unique(mymeta$window_id)) - nrow(valid_windows), "\n")

# Keep only valid windows in mymeta and all_windows
mymeta <- mymeta %>% filter(window_id %in% valid_windows$window_id)
all_windows <- all_windows %>% filter(window_id %in% valid_windows$window_id)

## ============================================================================
## CREATE SPATIAL SIGNALS AND CALCULATE SGCC
## ============================================================================
cat("\n=== Creating Spatial Signals ===\n")
grid_size <- 50
base_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)
cell_types_of_interest <- c(CELL_TYPE_1, CELL_TYPE_2)

all_signals <- list()
signal_names <- c()

for (window_id in unique(mymeta$window_id)) {
  for (cell_type in cell_types_of_interest) {
    cell_count <- sum(mymeta$window_id == window_id & mymeta$Annotation_pooled == cell_type)
    if (cell_count > 0) {
      signal_name <- paste(cell_type, window_id, sep = "_")
      signal_data <- CreateSCSignalData(mymeta, window_id, cell_type, all_windows, grid_size)
      if (sum(signal_data) > 0) {
        all_signals[[signal_name]] <- signal_data
        signal_names <- c(signal_names, signal_name)
      }
    }
  }
}

combined_signal_data <- base_grid
for (signal_name in signal_names) {
  combined_signal_data[[signal_name]] <- all_signals[[signal_name]]
}

cat("Total signals created:", length(signal_names), "\n")

cat("\n=== Running SGWT Analysis ===\n")
SG_combined <- initSGWT(
  data.in = combined_signal_data,
  x_col = "x", y_col = "y",
  signals = signal_names,
  J = 4, scaling_factor = 2,
  kernel_type = "mexican_hat"
)

SG_combined <- runSpecGraph(SG_combined, k = 8, laplacian_type = "normalized",
                            length_eigenvalue = grid_size^2, verbose = FALSE)
SG_combined <- runSGWT(SG_combined, verbose = FALSE)
plot_sgwt_decomposition(SG_combined, signal_name = SG_combined$Data$signals[1])
cat("\n=== Calculating SGCC ===\n")
sgcc_scores_df <- data.frame()

for (window_id in unique(mymeta$window_id)) {
  signal_1 <- paste(CELL_TYPE_1, window_id, sep = "_")
  signal_2 <- paste(CELL_TYPE_2, window_id, sep = "_")

  if (signal_1 %in% signal_names && signal_2 %in% signal_names) {
    sgcc_result <- runSGCC(signal_1, signal_2, SG = SG_combined, return_parts = TRUE)
    sgcc_row <- data.frame(
      window_id = window_id,
      SGCC_score = sgcc_result$S,
      stringsAsFactors = FALSE
    )
    sgcc_scores_df <- rbind(sgcc_scores_df, sgcc_row)
  }
}

cat("Total SGCC scores:", nrow(sgcc_scores_df), "\n")
sgcc_ranked <- sgcc_scores_df[order(sgcc_scores_df$SGCC_score), ]

## ============================================================================
## TRACK WINDOWS FOR VISUALIZATION (PLACEHOLDER - MOVED LATER)
## ============================================================================
# Visualization will be done after DEG analysis to show only windows used in DEG
sgcc_ranked_for_viz <- sgcc_ranked

# This section moved after DEG analysis

## ============================================================================
## CREATE PSEUDOBULK EXPRESSION
## ============================================================================
cat("\n=== Creating Pseudobulk Expression ===\n")
mydata_sub <- t(assay(spe_single_cell))
pseudobulk_result <- CreatePseudoBulk(
  mymeta = mymeta,
  mydata = mydata_sub,
  cell_types_for_pseudobulk = unique(mymeta$Annotation_pooled),
  min_cells_per_sample = MIN_CELLS_PER_WINDOW,
  verbose = TRUE
)

pseudobulk_expr <- pseudobulk_result$pseudobulk_expr
pseudobulk_expr <- pseudobulk_expr[!grepl("Negative|system", rownames(pseudobulk_expr), ignore.case = TRUE), ]

# Filter low expression genes
gene_sums <- rowSums(pseudobulk_expr)
expressed_genes <- names(gene_sums[gene_sums > 1])
pseudobulk_expr_filtered <- pseudobulk_expr[expressed_genes, ]

# Match with SGCC scores
expr_ct_samples <- colnames(pseudobulk_expr_filtered)
expr_ct_window_ids <- gsub(paste0("^", EXPRESSION_CELL_TYPE, "_"), "", expr_ct_samples)
available_window_ids <- intersect(expr_ct_window_ids, sgcc_ranked$window_id)
available_samples <- paste0(EXPRESSION_CELL_TYPE, "_", available_window_ids)

pseudobulk_expr_final <- pseudobulk_expr_filtered[, available_samples]
sample_sgcc_scores <- sgcc_ranked$SGCC_score[match(available_window_ids, sgcc_ranked$window_id)]
names(sample_sgcc_scores) <- available_samples

cat("Final pseudobulk samples:", length(available_samples), "\n")

# Save windows used in pseudobulk for later visualization
windows_used_in_pseudobulk <- available_window_ids

## ============================================================================
## CREATE ANNOTATION WITH SGCC GROUPS
## ============================================================================
cat("\n=== Creating Sample Annotations ===\n")

# Create annotation DataFrame
dfAnnotation <- data.frame(
  Sample = names(sample_sgcc_scores),
  SGCC = as.numeric(sample_sgcc_scores),
  Batch = ifelse(grepl("Rochester", names(sample_sgcc_scores)), "Rochester", "DFCI"),
  Batch2 = sub(".*(FOV[A-Za-z]+_[Cc][0-9]{2}).*", "\\1", names(sample_sgcc_scores)),
  stringsAsFactors = FALSE
)

# Add EBV status
mymeta_statusmapping <- mymeta %>%
  select(window_id, EBV_status) %>%
  distinct() %>%
  mutate(Sample = paste0(EXPRESSION_CELL_TYPE, "_", window_id))

dfAnnotation <- dfAnnotation %>% left_join(mymeta_statusmapping, by = "Sample")
dfAnnotation <- dfAnnotation[order(dfAnnotation$SGCC), ]
pseudobulk_expr_final <- pseudobulk_expr_final[, dfAnnotation$Sample]

# Define High and Low SGCC groups (top 1/3 and bottom 1/3)
n_samples <- nrow(dfAnnotation)
n_third <- floor(n_samples /4)

dfAnnotation$SGCC_Group <- NA
dfAnnotation$SGCC_Group[1:n_third] <- "Low"
dfAnnotation$SGCC_Group[(n_samples - n_third + 1):n_samples] <- "High"

# Remove middle third
dfAnnotation <- dfAnnotation[!is.na(dfAnnotation$SGCC_Group), ]
pseudobulk_expr_final <- pseudobulk_expr_final[, dfAnnotation$Sample]

## ============================================================================
## CALCULATE MACROPHAGE C1Q EXPRESSION PER SLIDING WINDOW
## ============================================================================
cat("\n=== Calculating Macrophage C1Q Expression per Window ===\n")

c1q_col <- intersect(c("C1q", "C1Q"), colnames(mymeta))

if (length(c1q_col) == 0) {
  warning("No C1q/C1Q column found in mymeta; Macrophage_C1Q annotation will be NA")
  dfAnnotation$Macrophage_C1Q <- NA_real_
} else {
  macrophage_c1q_by_window <- mymeta %>%
    as.data.frame() %>%
    filter(Annotation_pooled == "Macrophage") %>%
    group_by(window_id) %>%
    summarise(
      Macrophage_C1Q = mean(.data[[c1q_col[1]]], na.rm = TRUE),
      .groups = "drop"
    )

  dfAnnotation$Macrophage_C1Q <- macrophage_c1q_by_window$Macrophage_C1Q[
    match(
      gsub(paste0("^", EXPRESSION_CELL_TYPE, "_"), "", dfAnnotation$Sample),
      macrophage_c1q_by_window$window_id
    )
  ]
}

cat("Samples in Low SGCC:", sum(dfAnnotation$SGCC_Group == "Low"), "\n")
cat("Samples in High SGCC:", sum(dfAnnotation$SGCC_Group == "High"), "\n")
cat("EBV+ samples:", sum(dfAnnotation$EBV_status == "EBV+"), "\n")
cat("EBV- samples:", sum(dfAnnotation$EBV_status == "EBV-"), "\n")
cat("SGCC groups: ",print(table(dfAnnotation$SGCC_Group,dfAnnotation$EBV_status)))

## ============================================================================
## VISUALIZE EXTREME SAMPLES IN EACH SGCC GROUP (BY EBV STATUS)
## ============================================================================
cat("\n=== Creating Cell Type Layout Plots for Extreme Samples ===\n")

# Function to get extreme samples for a given EBV status and SGCC group
get_extreme_samples_by_ebv <- function(ebv_status_val, sgcc_group_val) {
  subset_data <- dfAnnotation[dfAnnotation$EBV_status == ebv_status_val &
                               dfAnnotation$SGCC_Group == sgcc_group_val, ]
  subset_data <- subset_data[order(subset_data$SGCC), ]

  if (nrow(subset_data) == 0) {
    return(list(lowest = NULL, highest = NULL))
  }

  return(list(
    lowest = subset_data$Sample[1],
    highest = subset_data$Sample[nrow(subset_data)]
  ))
}

# EBV+ samples
cat("\n--- EBV+ Group ---\n")
ebv_pos_low <- get_extreme_samples_by_ebv("EBV+", "Low")
ebv_pos_high <- get_extreme_samples_by_ebv("EBV+", "High")

sample_ebvpos_low_lowest <- ebv_pos_low$lowest
sample_ebvpos_low_highest <- ebv_pos_low$highest
sample_ebvpos_high_lowest <- ebv_pos_high$lowest
sample_ebvpos_high_highest <- ebv_pos_high$highest

# EBV- samples
cat("\n--- EBV- Group ---\n")
ebv_neg_low <- get_extreme_samples_by_ebv("EBV-", "Low")
ebv_neg_high <- get_extreme_samples_by_ebv("EBV-", "High")

sample_ebvneg_low_lowest <- ebv_neg_low$lowest
sample_ebvneg_low_highest <- ebv_neg_low$highest
sample_ebvneg_high_lowest <- ebv_neg_high$lowest
sample_ebvneg_high_highest <- ebv_neg_high$highest

# Extract window IDs and get info
get_sample_info <- function(sample_name) {
  if (is.null(sample_name)) return(NULL)

  window_id <- gsub(paste0("^", EXPRESSION_CELL_TYPE, "_"), "", sample_name)
  sgcc_score <- dfAnnotation[dfAnnotation$Sample == sample_name, "SGCC"]

  info <- mymeta[mymeta$window_id == window_id, c("TMA_core", "EBV_status")] %>%
    distinct()

  return(list(
    window_id = window_id,
    sgcc_score = sgcc_score,
    tma_core = info$TMA_core[1],
    ebv_status = info$EBV_status[1]
  ))
}

# Get info for all 8 samples
info_ebvpos_low_lowest <- get_sample_info(sample_ebvpos_low_lowest)
info_ebvpos_low_highest <- get_sample_info(sample_ebvpos_low_highest)
info_ebvpos_high_lowest <- get_sample_info(sample_ebvpos_high_lowest)
info_ebvpos_high_highest <- get_sample_info(sample_ebvpos_high_highest)

info_ebvneg_low_lowest <- get_sample_info(sample_ebvneg_low_lowest)
info_ebvneg_low_highest <- get_sample_info(sample_ebvneg_low_highest)
info_ebvneg_high_lowest <- get_sample_info(sample_ebvneg_high_lowest)
info_ebvneg_high_highest <- get_sample_info(sample_ebvneg_high_highest)

# Print sample information
cat("\nEBV+ Low SGCC - Lowest:", sample_ebvpos_low_lowest, "(SGCC =", round(info_ebvpos_low_lowest$sgcc_score, 4), ")\n")
cat("EBV+ Low SGCC - Highest:", sample_ebvpos_low_highest, "(SGCC =", round(info_ebvpos_low_highest$sgcc_score, 4), ")\n")
cat("EBV+ High SGCC - Lowest:", sample_ebvpos_high_lowest, "(SGCC =", round(info_ebvpos_high_lowest$sgcc_score, 4), ")\n")
cat("EBV+ High SGCC - Highest:", sample_ebvpos_high_highest, "(SGCC =", round(info_ebvpos_high_highest$sgcc_score, 4), ")\n")

cat("\nEBV- Low SGCC - Lowest:", sample_ebvneg_low_lowest, "(SGCC =", round(info_ebvneg_low_lowest$sgcc_score, 4), ")\n")
cat("EBV- Low SGCC - Highest:", sample_ebvneg_low_highest, "(SGCC =", round(info_ebvneg_low_highest$sgcc_score, 4), ")\n")
cat("EBV- High SGCC - Lowest:", sample_ebvneg_high_lowest, "(SGCC =", round(info_ebvneg_high_lowest$sgcc_score, 4), ")\n")
cat("EBV- High SGCC - Highest:", sample_ebvneg_high_highest, "(SGCC =", round(info_ebvneg_high_highest$sgcc_score, 4), ")\n")

# Function to create cell type layout plot
create_cell_layout_plot <- function(window_id, sgcc_score, tma_core, ebv_status, sgcc_group, position) {
  # Get only the two cell types of interest from this window
  window_cells <- mymeta[mymeta$window_id == window_id &
                         mymeta$Annotation_pooled %in% c(CELL_TYPE_1, CELL_TYPE_2), ]

  if (nrow(window_cells) == 0) {
    return(NULL)
  }

  # Create color palette for the two cell types
  cell_type_colors <- c(
     "CD4 T" = "#3cb44e",
      "Macrophage" = "#8e52a1",
      "B cell" = "#ec1e2f",
      "CD8 T" = "#1f77b4"
  )

  p <- ggplot(window_cells, aes(x = X_cent, y = Y_cent, color = Annotation_pooled)) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_manual(values = cell_type_colors, name = "Cell Type") +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 9, hjust = 0.5)
    ) +
    scale_y_reverse() +
    labs(
      title = paste0(ebv_status, " | ", sgcc_group, " SGCC - ", position),
      subtitle = paste0(tma_core, " | SGCC: ", round(sgcc_score, 4)),
      x = "X Position", y = "Y Position"
    )

  return(p)
}

# Create the 8 plots (4 for EBV+, 4 for EBV-)
plot1 <- create_cell_layout_plot(info_ebvpos_low_lowest$window_id, info_ebvpos_low_lowest$sgcc_score,
                                 info_ebvpos_low_lowest$tma_core, info_ebvpos_low_lowest$ebv_status,
                                 "Low", "Lowest")
plot2 <- create_cell_layout_plot(info_ebvpos_low_highest$window_id, info_ebvpos_low_highest$sgcc_score,
                                 info_ebvpos_low_highest$tma_core, info_ebvpos_low_highest$ebv_status,
                                 "Low", "Highest")
plot3 <- create_cell_layout_plot(info_ebvpos_high_lowest$window_id, info_ebvpos_high_lowest$sgcc_score,
                                 info_ebvpos_high_lowest$tma_core, info_ebvpos_high_lowest$ebv_status,
                                 "High", "Lowest")
plot4 <- create_cell_layout_plot(info_ebvpos_high_highest$window_id, info_ebvpos_high_highest$sgcc_score,
                                 info_ebvpos_high_highest$tma_core, info_ebvpos_high_highest$ebv_status,
                                 "High", "Highest")

plot5 <- create_cell_layout_plot(info_ebvneg_low_lowest$window_id, info_ebvneg_low_lowest$sgcc_score,
                                 info_ebvneg_low_lowest$tma_core, info_ebvneg_low_lowest$ebv_status,
                                 "Low", "Lowest")
plot6 <- create_cell_layout_plot(info_ebvneg_low_highest$window_id, info_ebvneg_low_highest$sgcc_score,
                                 info_ebvneg_low_highest$tma_core, info_ebvneg_low_highest$ebv_status,
                                 "Low", "Highest")
plot7 <- create_cell_layout_plot(info_ebvneg_high_lowest$window_id, info_ebvneg_high_lowest$sgcc_score,
                                 info_ebvneg_high_lowest$tma_core, info_ebvneg_high_lowest$ebv_status,
                                 "High", "Lowest")
plot8 <- create_cell_layout_plot(info_ebvneg_high_highest$window_id, info_ebvneg_high_highest$sgcc_score,
                                 info_ebvneg_high_highest$tma_core, info_ebvneg_high_highest$ebv_status,
                                 "High", "Highest")

# Save the 8 plots in a 4x2 grid
extreme_samples_file <- file.path(results_path, paste0("CosMX4TMANew_", analysis_label, "_SGCC_extreme_samples_cell_layout.pdf"))
pdf(extreme_samples_file, width = 25, height = 10)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol = 4,
             top = paste0("Cell Type Layout for Extreme SGCC Samples by EBV Status (",
                         CELL_TYPE_1, " - ", CELL_TYPE_2, ")"))
dev.off()
cat("Extreme samples cell layout plots saved to:", extreme_samples_file, "\n")

## ============================================================================
## BATCH CORRECTION
## ============================================================================
cat("\n=== Performing Batch Correction ===\n")

spe_for_bc <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = pseudobulk_expr_final),
  colData = S4Vectors::DataFrame(dfAnnotation)
)

spe_for_bc$total_counts <- colSums(assay(spe_for_bc, "counts"))
spe_normalized <- geomxNorm(spe_for_bc, method = "CPM", log = TRUE)
spe_with_ncgs <- findNCGs(spe_normalized, batch_name = "Batch2", top_n = 1500)
spe_batch_corrected <- geomxBatchCorrection(
  spe_with_ncgs,
  factors = c("EBV_status", "SGCC_Group"),
  NCGs = S4Vectors::metadata(spe_with_ncgs)$NCGs,
  k = 2,
  covariates = NULL,
  isLog = TRUE
)

norm_expr_corrected <- assay(spe_batch_corrected, "logcounts")
cat("Batch correction completed.\n")

## ============================================================================
## DEG ANALYSIS: HIGH vs LOW SGCC, STRATIFIED BY EBV STATUS
## ============================================================================
cat("\n=== Running DEG Analysis (High vs Low SGCC) ===\n")

# Function to run limma DEG analysis on batch-corrected expression
run_deg_analysis <- function(expr_matrix, annotation_df, group_col, comparison_name) {
  cat("\n--- DEG Analysis:", comparison_name, "---\n")

  # Design matrix with SGCC group and Batch
  design <- model.matrix(~ 0 + annotation_df[[group_col]] + annotation_df$Batch)
  colnames(design) <- gsub("annotation_df\\[\\[group_col\\]\\]", "", colnames(design))
  colnames(design) <- gsub("annotation_df\\$Batch", "Batch", colnames(design))

  # Fit linear model using limma on batch-corrected log-expression
  fit <- lmFit(expr_matrix, design)

  # Contrast: High vs Low
  contrast_matrix <- makeContrasts(High - Low, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # Get results
  results <- topTable(fit2, number = Inf, sort.by = "none")
  results$Gene <- rownames(results)
  colnames(results)[colnames(results) == "adj.P.Val"] <- "FDR"

  # Add column indicating which group has higher expression
  results$Higher_in <- ifelse(results$logFC > 0, "High_SGCC", "Low_SGCC")

  # Significant genes
  sig_genes <- rownames(results[results$P.Value < 0.01 & abs(results$logFC) > 0.0, ])
  sig_genes_high <- rownames(results[results$P.Value < 0.01 & results$logFC > 0.0, ])
  sig_genes_low <- rownames(results[results$P.Value < 0.01 & results$logFC < 0, ])

  cat("Total DEGs (P.Value < 0.005, |logFC| > 0.5):", length(sig_genes), "\n")
  cat("  Upregulated in High SGCC:", length(sig_genes_high), "\n")
  cat("  Downregulated in High SGCC (Upregulated in Low SGCC):", length(sig_genes_low), "\n")

  return(list(results = results, sig_genes = sig_genes,
              sig_genes_high = sig_genes_high, sig_genes_low = sig_genes_low))
}

# Split by EBV status
dfAnnotation_ebv_pos <- dfAnnotation[dfAnnotation$EBV_status == "EBV+", ]
dfAnnotation_ebv_neg <- dfAnnotation[dfAnnotation$EBV_status == "EBV-", ]

# DEG for EBV+
deg_ebv_pos <- run_deg_analysis(
  expr_matrix = norm_expr_corrected[, dfAnnotation_ebv_pos$Sample],
  annotation_df = dfAnnotation_ebv_pos,
  group_col = "SGCC_Group",
  comparison_name = "EBV+ (High vs Low SGCC)"
)

# DEG for EBV-
deg_ebv_neg <- run_deg_analysis(
  expr_matrix = norm_expr_corrected[, dfAnnotation_ebv_neg$Sample],
  annotation_df = dfAnnotation_ebv_neg,
  group_col = "SGCC_Group",
  comparison_name = "EBV- (High vs Low SGCC)"
)
View(deg_ebv_pos$results[deg_ebv_pos$results$P.Value < 0.01,])
View(deg_ebv_neg$results[deg_ebv_neg$results$P.Value < 0.01,])
# Save DEG results
write.csv(deg_ebv_pos$results,
          file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_EBVpos_HighVsLow.csv")),
          row.names = FALSE)
write.csv(deg_ebv_neg$results,
          file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_EBVneg_HighVsLow.csv")),
          row.names = FALSE)

# Aggregate results
all_sig_genes <- unique(c(deg_ebv_pos$sig_genes, deg_ebv_neg$sig_genes))
cat("\nTotal unique DEGs across both groups:", length(all_sig_genes), "\n")
cat("Overlap between EBV+ and EBV-:", length(intersect(deg_ebv_pos$sig_genes, deg_ebv_neg$sig_genes)), "genes\n")

## ============================================================================
## CREATE HEATMAPS
## ============================================================================
cat("\n=== Creating Heatmaps ===\n")

create_deg_heatmap <- function(expr_matrix, annotation_df, sig_genes, title) {
  if (length(sig_genes) == 0) return(NULL)

  # Order by SGCC
  anno_ordered <- annotation_df[order(annotation_df$SGCC), ]
  expr_subset <- expr_matrix[sig_genes, anno_ordered$Sample, drop = FALSE]
  expr_scaled <- t(scale(t(expr_subset)))

  # Column annotation
  col_annotation <- HeatmapAnnotation(
    SGCC_Score = anno_lines(anno_ordered$SGCC, which = "column",
                            gp = gpar(col = "red", lwd = 2)),
    annotation_name_gp = gpar(fontsize = 8)
  )

  ht <- Heatmap(
    expr_scaled,
    name = "Expression\n(Z-score)",
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    top_annotation = col_annotation,
    col = colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b")),
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE,
    column_title = title,
    column_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
  return(ht)
}

# EBV+ heatmap
if (length(deg_ebv_pos$sig_genes) > 0) {
  ht_ebv_pos <- create_deg_heatmap(
    norm_expr_corrected, dfAnnotation_ebv_pos,
    deg_ebv_pos$sig_genes,
    paste0(EXPRESSION_CELL_TYPE, " DEGs (EBV+: High vs Low SGCC)")
  )

  pdf(file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_heatmap_EBVpos.pdf")),
      width = 12, height = 10)
  draw(ht_ebv_pos, heatmap_legend_side = "right")
  dev.off()
}

# EBV- heatmap
if (length(deg_ebv_neg$sig_genes) > 0) {
  ht_ebv_neg <- create_deg_heatmap(
    norm_expr_corrected, dfAnnotation_ebv_neg,
    deg_ebv_neg$sig_genes,
    paste0(EXPRESSION_CELL_TYPE, " DEGs (EBV-: High vs Low SGCC)")
  )

  pdf(file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_heatmap_EBVneg.pdf")),
      width = 12, height = 10)
  draw(ht_ebv_neg, heatmap_legend_side = "right")
  dev.off()
}

# Combined heatmap (EBV+ and EBV- together, ranked by SGCC)
all_sig_genes_combined <- unique(c(deg_ebv_pos$sig_genes, deg_ebv_neg$sig_genes))
if (length(all_sig_genes_combined) > 0) {
  # Use all samples (both EBV+ and EBV-), ordered by SGCC
  dfAnnotation_all <- dfAnnotation[order(dfAnnotation$SGCC), ]
  expr_subset_all <- norm_expr_corrected[all_sig_genes_combined, dfAnnotation_all$Sample, drop = FALSE]
  expr_scaled_all <- t(scale(t(expr_subset_all)))

  # Column annotation with EBV status
  col_annotation_combined <- HeatmapAnnotation(
    EBV_Status = dfAnnotation_all$EBV_status,
    SGCC_Score = anno_lines(dfAnnotation_all$SGCC, which = "column",
                            gp = gpar(col = "red", lwd = 2)),
    col = list(
      EBV_Status = c("EBV+" = "#d95f02", "EBV-" = "#1b9e77")
    ),
    annotation_name_gp = gpar(fontsize = 8)
  )

  ht_combined <- Heatmap(
    expr_scaled_all,
    name = "Expression\n(Z-score)",
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    top_annotation = col_annotation_combined,
    col = colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b")),
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE,
    column_title = paste0(EXPRESSION_CELL_TYPE, " DEGs (EBV+ & EBV- Combined, Ranked by SGCC)"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold")
  )

  pdf(file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_heatmap_Combined.pdf")),
      width = 14, height = 10)
  draw(ht_combined, heatmap_legend_side = "right")
  dev.off()
}

## ============================================================================
## VISUALIZE SLIDING WINDOWS USED IN DEG ANALYSIS
## ============================================================================
cat("\n=== Creating Spatial Plots for Windows Used in DEG Analysis ===\n")

# Get windows that were actually used in DEG analysis
deg_window_ids <- gsub(paste0("^", EXPRESSION_CELL_TYPE, "_"), "", dfAnnotation$Sample)

# Add EBV status, TMA core info, and SGCC group
sgcc_ranked_viz <- sgcc_ranked %>%
  filter(window_id %in% deg_window_ids) %>%
  left_join(mymeta %>% select(window_id, EBV_status, TMA_core) %>% distinct(),
            by = "window_id") %>%
  left_join(dfAnnotation %>%
              mutate(window_id = gsub(paste0("^", EXPRESSION_CELL_TYPE, "_"), "", Sample)) %>%
              select(window_id, SGCC_Group),
            by = "window_id")

cat("Total windows used in DEG analysis:", nrow(sgcc_ranked_viz), "\n")
cat("  High SGCC:", sum(sgcc_ranked_viz$SGCC_Group == "High", na.rm = TRUE), "\n")
cat("  Low SGCC:", sum(sgcc_ranked_viz$SGCC_Group == "Low", na.rm = TRUE), "\n")

# Select windows for visualization: 10 from High SGCC, 10 from Medium, 10 from Low SGCC
# High and Low come from DEG-filtered windows, Medium comes from all windows
high_sgcc_windows <- sgcc_ranked_viz %>%
  filter(SGCC_Group == "High") %>%
  arrange(desc(SGCC_score)) %>%
  slice_head(n = 10)

low_sgcc_windows <- sgcc_ranked_viz %>%
  filter(SGCC_Group == "Low") %>%
  arrange(SGCC_score) %>%
  slice_head(n = 10)

# For medium windows, get from the middle range of ALL windows (not just DEG-filtered)
n_total_all <- nrow(sgcc_ranked)
middle_start <- round(n_total_all/2) - 4
middle_end <- round(n_total_all/2) + 5
medium_sgcc_windows <- sgcc_ranked[middle_start:middle_end, ] %>%
  left_join(mymeta %>% select(window_id, EBV_status, TMA_core) %>% distinct(),
            by = "window_id") %>%
  mutate(SGCC_Group = "Medium")

viz_data <- bind_rows(
  low_sgcc_windows %>% mutate(Category = "Low SGCC"),
  medium_sgcc_windows %>% mutate(Category = "Medium SGCC"),
  high_sgcc_windows %>% mutate(Category = "High SGCC")
)

cat("Selected", nrow(viz_data), "windows for visualization\n")
cat("  Low SGCC:", sum(viz_data$Category == "Low SGCC"), "\n")
cat("  Medium SGCC:", sum(viz_data$Category == "Medium SGCC"), "\n")
cat("  High SGCC:", sum(viz_data$Category == "High SGCC"), "\n")

# Create spatial plots for each window
plot_list <- list()

for (i in seq_len(nrow(viz_data))) {
  window_id <- viz_data$window_id[i]
  sgcc_score <- round(viz_data$SGCC_score[i], 4)
  tma_core <- viz_data$TMA_core[i]
  ebv_status <- viz_data$EBV_status[i]
  sgcc_group <- viz_data$SGCC_Group[i]

  # Get cells from this window - only the two cell types of interest
  window_cells <- mymeta[mymeta$window_id == window_id &
                           mymeta$Annotation_pooled %in% c(CELL_TYPE_1, CELL_TYPE_2), ]

  if (nrow(window_cells) > 0) {
    # Create color palette for cell types
    cell_type_colors <- c(
      "CD4 T" = "#3cb44e",
      "Macrophage" = "#8e52a1",
      "B cell" = "#ec1e2f",
      "CD8 T" = "#1f77b4"
    )
    color_palette <- cell_type_colors[c(CELL_TYPE_1, CELL_TYPE_2)]
    names(color_palette) <- c(CELL_TYPE_1, CELL_TYPE_2)

    p <- ggplot(window_cells, aes(x = X_cent, y = Y_cent, color = Annotation_pooled)) +
      geom_point(size = 1.2, alpha = 0.8) +
      scale_color_manual(values = color_palette) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, hjust = 0.5),
        plot.subtitle = element_text(size = 7, hjust = 0.5)
      ) +
      labs(
        title = paste0(tma_core, " (", ebv_status, ") - ", sgcc_group),
        subtitle = paste0("SGCC: ", sgcc_score),
        x = "X", y = "Y"
      )

    plot_list[[i]] <- p
  }
}

# Arrange all plots in a grid (3 rows: Low, Medium, High)
cat("Creating composite plot with", length(plot_list), "windows...\n")

# Save combined plot
phenotype_file <- file.path(results_path, paste0("CosMX4TMANew_", analysis_label, "_SGCC_DEG_windows_spatial.pdf"))
pdf(phenotype_file, width = 20, height = 20)
grid.arrange(grobs = plot_list, ncol = 5,
             top = paste0("Sliding Windows Ranked by SGCC (",
                         CELL_TYPE_1, " - ", CELL_TYPE_2, ")\n",
                         "Row 1: Low SGCC (n=10) | Row 2: Medium SGCC (n=10) | Row 3: High SGCC (n=10)\n",
                         "DEG Analysis Used: High=", sum(sgcc_ranked_viz$SGCC_Group == "High", na.rm = TRUE),
                         ", Low=", sum(sgcc_ranked_viz$SGCC_Group == "Low", na.rm = TRUE)))
dev.off()
cat("Spatial window plots saved to:", phenotype_file, "\n")

## ============================================================================
## ENRICHR PATHWAY ENRICHMENT
## ============================================================================
cat("\n=== Running Enrichr Pathway Enrichment ===\n")

run_enrichr_pathway <- function(gene_vec, label, databases = c("Reactome_2022", "GO_Biological_Process_2025"),
                                results_dir = ".", top_n = 20) {
  if (length(gene_vec) < 3) {
    cat("Not enough genes for Enrichr (n < 3)\n")
    return(NULL)
  }

  enr_list <- enrichR::enrichr(gene_vec, databases)

  for (db in databases) {
    if (!is.null(enr_list[[db]]) && nrow(enr_list[[db]]) > 0) {
      csv_file <- file.path(results_dir, paste0("Enrichr_", label, "_", db, "_full.csv"))
      write.csv(enr_list[[db]], csv_file, row.names = FALSE)
    }
  }

  db_to_plot <- NULL
  for (db in databases) {
    if (!is.null(enr_list[[db]]) && nrow(enr_list[[db]]) > 0) {
      db_to_plot <- db
      break
    }
  }

  if (is.null(db_to_plot) || !"Adjusted.P.value" %in% colnames(enr_list[[db_to_plot]])) {
    return(list(results = enr_list, plot = NULL))
  }

  res <- enr_list[[db_to_plot]]
  res$log10AdjP <- -log10(res$Adjusted.P.value + 1e-16)
  res_sig <- res[res$Adjusted.P.value < 0.05, , drop = FALSE]

  if (nrow(res_sig) > 0) {
    res_sig <- res_sig[order(res_sig$Adjusted.P.value), , drop = FALSE]
    top_res <- res_sig[seq_len(min(top_n, nrow(res_sig))), , drop = FALSE]
    top_res$Term <- factor(top_res$Term, levels = rev(top_res$Term))

    p <- ggplot(top_res, aes(x = log10AdjP, y = Term)) +
      geom_col(fill = "#1b9e77") +
      labs(title = paste0(label, " - ", db_to_plot),
           x = "-log10(adj. P-value)", y = "Pathway") +
      theme_bw(base_size = 10)

    plot_file <- file.path(results_dir, paste0("Enrichr_", label, "_", db_to_plot, "_top", top_n, ".pdf"))
    ggsave(plot_file, plot = p, width = 8, height = 6)
  }

  return(list(results = enr_list, plot = NULL))
}


# Enrichr for EBV+ - High SGCC DEGs
enrichr_ebv_pos_high <- run_enrichr_pathway(
  gene_vec = deg_ebv_pos$sig_genes_high,
  label = paste0(expression_label, "_", analysis_label, "_DEG_EBVpos_HighSGCC"),
  results_dir = results_path
)

# Enrichr for EBV+ - Low SGCC DEGs
enrichr_ebv_pos_low <- run_enrichr_pathway(
  gene_vec = deg_ebv_pos$sig_genes_low,
  label = paste0(expression_label, "_", analysis_label, "_DEG_EBVpos_LowSGCC"),
  results_dir = results_path
)


# Enrichr for EBV- - High SGCC DEGs
enrichr_ebv_neg_high <- run_enrichr_pathway(
  gene_vec = deg_ebv_neg$sig_genes_high,
  label = paste0(expression_label, "_", analysis_label, "_DEG_EBVneg_HighSGCC"),
  results_dir = results_path
)

# Enrichr for EBV- - Low SGCC DEGs
enrichr_ebv_neg_low <- run_enrichr_pathway(
  gene_vec = deg_ebv_neg$sig_genes_low,
  label = paste0(expression_label, "_", analysis_label, "_DEG_EBVneg_LowSGCC"),
  results_dir = results_path
)

cat("\n=== Pathway Enrichment Summary ===\n")
cat("EBV+ High SGCC genes:", length(deg_ebv_pos$sig_genes_high), "\n")
cat("EBV+ Low SGCC genes:", length(deg_ebv_pos$sig_genes_low), "\n")
cat("EBV- High SGCC genes:", length(deg_ebv_neg$sig_genes_high), "\n")
cat("EBV- Low SGCC genes:", length(deg_ebv_neg$sig_genes_low), "\n")

## ============================================================================
## SAVE SEPARATE DEG LISTS BY SGCC GROUP
## ============================================================================
# Save High SGCC and Low SGCC gene lists separately
write.csv(deg_ebv_pos$results[deg_ebv_pos$results$P.Value < 0.01 & deg_ebv_pos$results$logFC > 0.25, ],
          file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_EBVpos_HighSGCC.csv")),
          row.names = FALSE)
write.csv(deg_ebv_pos$results[deg_ebv_pos$results$P.Value < 0.01 & deg_ebv_pos$results$logFC < -0.25, ],
          file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_EBVpos_LowSGCC.csv")),
          row.names = FALSE)

write.csv(deg_ebv_neg$results[deg_ebv_neg$results$P.Value < 0.01 & deg_ebv_neg$results$logFC > 0.25, ],
          file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_EBVneg_HighSGCC.csv")),
          row.names = FALSE)
write.csv(deg_ebv_neg$results[deg_ebv_neg$results$P.Value < 0.01 & deg_ebv_neg$results$logFC < -0.25, ],
          file.path(results_path, paste0("CosMX4TMANew_", expression_label, "_", analysis_label, "_DEG_EBVneg_LowSGCC.csv")),
          row.names = FALSE)

## ============================================================================
## GSVA ANALYSIS
## ============================================================================
cat("\n=== Running GSVA Analysis ===\n")

run_gsva_analysis <- function(enrichr_results, expr_matrix, annotation_df, label, results_dir) {
  if (is.null(enrichr_results)) return(NULL)

  databases <- c("Reactome_2022", "GO_Biological_Process_2025")

  for (db_name in databases) {
    if (is.null(enrichr_results$results[[db_name]])) next

    db_results <- enrichr_results$results[[db_name]]
    sig_pathways <- db_results[db_results$Adjusted.P.value < 0.05, ]

    if (nrow(sig_pathways) == 0) next

    top_pathways <- sig_pathways[order(sig_pathways$Adjusted.P.value), ]
    top_pathways <- top_pathways[seq_len(min(30, nrow(top_pathways))), ]

    gsva_gene_sets <- list()
    for (i in seq_len(nrow(top_pathways))) {
      pathway_name <- top_pathways$Term[i]
      genes_vec <- unlist(strsplit(top_pathways$Genes[i], ";"))

      if (db_name == "Reactome_2022") {
        pathway_display <- sub(" R-HSA-[0-9]+$", "", pathway_name)
      } else if (db_name == "GO_Biological_Process_2025") {
        pathway_display <- sub(" \\(GO:[0-9]+\\)$", "", pathway_name)
      } else {
        pathway_display <- pathway_name
      }

      gsva_gene_sets[[pathway_display]] <- genes_vec
    }

    SetGSVAPar <- gsvaParam(
      exprData = as.matrix(expr_matrix[, annotation_df$Sample]),
      geneSets = gsva_gene_sets,
      minSize = 3, maxSize = 500,
      kcdf = "Gaussian"
    )
    gsva_scores <- gsva(SetGSVAPar, verbose = FALSE)

    anno_ordered <- annotation_df[order(annotation_df$SGCC), ]
    gsva_scores_ordered <- gsva_scores[, anno_ordered$Sample, drop = FALSE]
    gsva_scores_scaled <- t(scale(t(gsva_scores_ordered)))

    col_annotation <- HeatmapAnnotation(
      SGCC_Group = anno_ordered$SGCC_Group,
      SGCC_Score = anno_lines(anno_ordered$SGCC, which = "column",
                              gp = gpar(col = "red", lwd = 2)),
      col = list(SGCC_Group = c("High" = "#d95f02", "Low" = "#1b9e77")),
      annotation_name_gp = gpar(fontsize = 8)
    )

    ht_gsva <- Heatmap(
      gsva_scores_scaled,
      name = "GSVA\n(Z-score)",
      cluster_columns = FALSE,
      cluster_rows = TRUE,
      top_annotation = col_annotation,
      col = colorRamp2(c(-2, 0, 2), c("#b953a0", "#0f0b0e", "#efea1f")),
      row_names_gp = gpar(fontsize = 8),
      show_column_names = FALSE,
      column_title = paste0("GSVA: ", db_name, " (", label, ")"),
      column_title_gp = gpar(fontsize = 10, fontface = "bold")
    )

    pdf(file.path(results_dir, paste0("CosMX4TMANew_", label, "_GSVA_", db_name, "_heatmap.pdf")),
        width = 14, height = 12)
    draw(ht_gsva, heatmap_legend_side = "right")
    dev.off()

    write.csv(gsva_scores,
              file.path(results_dir, paste0("CosMX4TMANew_", label, "_GSVA_", db_name, "_scores.csv")))
  }
}

## ============================================================================
## COMBINED GSVA ANALYSIS BY EBV STATUS (High and Low SGCC pathways together)
## ============================================================================
cat("\n=== Running Combined GSVA Analysis (High & Low SGCC pathways, Split by EBV status) ===\n")

run_gsva_combined_by_ebv <- function(enrichr_high, enrichr_low, expr_matrix,
                                     annotation_df, label, results_dir) {
  databases <- c("Reactome_2022", "GO_Biological_Process_2025")

  for (db_name in databases) {
    all_pathways <- list()

    # Collect pathways from High SGCC
    if (!is.null(enrichr_high) && !is.null(enrichr_high$results[[db_name]])) {
      db_results_high <- enrichr_high$results[[db_name]]
      sig_pathways_high <- db_results_high[db_results_high$Adjusted.P.value < 0.05, ]

      if (nrow(sig_pathways_high) > 0) {
        top_pathways_high <- sig_pathways_high[order(sig_pathways_high$Adjusted.P.value), ]
        top_pathways_high <- top_pathways_high[seq_len(min(15, nrow(top_pathways_high))), ]
        all_pathways[["High_SGCC"]] <- top_pathways_high
      }
    }

    # Collect pathways from Low SGCC
    if (!is.null(enrichr_low) && !is.null(enrichr_low$results[[db_name]])) {
      db_results_low <- enrichr_low$results[[db_name]]
      sig_pathways_low <- db_results_low[db_results_low$Adjusted.P.value < 0.05, ]

      if (nrow(sig_pathways_low) > 0) {
        top_pathways_low <- sig_pathways_low[order(sig_pathways_low$Adjusted.P.value), ]
        top_pathways_low <- top_pathways_low[seq_len(min(15, nrow(top_pathways_low))), ]
        all_pathways[["Low_SGCC"]] <- top_pathways_low
      }
    }

    if (length(all_pathways) == 0) {
      cat("No significant pathways found for", label, db_name, "\n")
      next
    }

    # Combine unique pathways
    combined_pathways <- do.call(rbind, all_pathways)
    combined_pathways <- combined_pathways[!duplicated(combined_pathways$Term), ]
    combined_pathways <- combined_pathways[order(combined_pathways$Adjusted.P.value), ]
    combined_pathways <- combined_pathways[seq_len(min(30, nrow(combined_pathways))), ]

    # Create gene sets
    gsva_gene_sets <- list()
    for (i in seq_len(nrow(combined_pathways))) {
      pathway_name <- combined_pathways$Term[i]
      genes_vec <- unlist(strsplit(combined_pathways$Genes[i], ";"))

      if (db_name == "Reactome_2022") {
        pathway_display <- sub(" R-HSA-[0-9]+$", "", pathway_name)
      } else if (db_name == "GO_Biological_Process_2025") {
        pathway_display <- sub(" \\(GO:[0-9]+\\)$", "", pathway_name)
      } else {
        pathway_display <- pathway_name
      }

      gsva_gene_sets[[pathway_display]] <- genes_vec
    }

    # Run GSVA on samples from this EBV group
    SetGSVAPar <- gsvaParam(
      exprData = as.matrix(expr_matrix[, annotation_df$Sample]),
      geneSets = gsva_gene_sets,
      minSize = 3, maxSize = 500,
      kcdf = "Gaussian"
    )
    gsva_scores <- gsva(SetGSVAPar, verbose = FALSE)

    # Order by SGCC within this EBV group
    anno_ordered <- annotation_df[order(annotation_df$SGCC), ]
    gsva_scores_ordered <- gsva_scores[, anno_ordered$Sample, drop = FALSE]
    gsva_scores_scaled <- t(scale(t(gsva_scores_ordered)))

    # Column annotation
    col_annotation <- HeatmapAnnotation(
     # SGCC_Group = anno_ordered$SGCC_Group,
      SGCC_Score = anno_lines(anno_ordered$SGCC, which = "column",
                              gp = gpar(col = "red", lwd = 2)),
      col = list(
        SGCC_Group = c("High" = "#d95f02", "Low" = "#1b9e77")
      ),
      annotation_name_gp = gpar(fontsize = 8)
    )

    ht_gsva <- Heatmap(
      gsva_scores_scaled,
      name = "GSVA\n(Z-score)",
      cluster_columns = FALSE,
      cluster_rows = TRUE,
      top_annotation = col_annotation,
      col = colorRamp2(c(-2, 0, 2), c("#b953a0", "#0f0b0e", "#efea1f")),
      row_names_gp = gpar(fontsize = 8),
      show_column_names = FALSE,
      column_title = paste0("GSVA: ", db_name, " (", label, ", Ranked by SGCC)"),
      column_title_gp = gpar(fontsize = 10, fontface = "bold")
    )

    pdf(file.path(results_dir, paste0("CosMX4TMANew_", label, "_GSVA_", db_name, "_heatmap.pdf")),
        width = 14, height = 12)
    draw(ht_gsva, heatmap_legend_side = "right")
    dev.off()

    write.csv(gsva_scores,
              file.path(results_dir, paste0("CosMX4TMANew_", label, "_GSVA_", db_name, "_scores.csv")))

    cat("GSVA heatmap for", label, db_name, "saved\n")
  }
}

# Run GSVA for EBV+ (combining High and Low SGCC pathways)
run_gsva_combined_by_ebv(
  enrichr_ebv_pos_high, enrichr_ebv_pos_low,
  norm_expr_corrected, dfAnnotation_ebv_pos,
  paste0(expression_label, "_", analysis_label, "_EBVpos"),
  results_path
)

# Run GSVA for EBV- (combining High and Low SGCC pathways)
run_gsva_combined_by_ebv(
  enrichr_ebv_neg_high, enrichr_ebv_neg_low,
  norm_expr_corrected, dfAnnotation_ebv_neg,
  paste0(expression_label, "_", analysis_label, "_EBVneg"),
  results_path
)

## ============================================================================
## COMBINED GSVA ANALYSIS - ALL SAMPLES (EBV+ and EBV- on same heatmap)
## ============================================================================
cat("\n=== Running Combined GSVA Analysis (All Samples, EBV+ & EBV- together) ===\n")

run_gsva_all_samples <- function(enrichr_pos_high, enrichr_pos_low,
                                 enrichr_neg_high, enrichr_neg_low,
                                 expr_matrix, annotation_df, label, results_dir) {
  databases <- c("Reactome_2022", "GO_Biological_Process_2025")

  for (db_name in databases) {
    all_pathways <- list()

    # Collect pathways from all four groups
    for (enrichr_result in list(enrichr_pos_high, enrichr_pos_low, enrichr_neg_high, enrichr_neg_low)) {
      if (!is.null(enrichr_result) && !is.null(enrichr_result$results[[db_name]])) {
        db_results <- enrichr_result$results[[db_name]]
        sig_pathways <- db_results[db_results$Adjusted.P.value < 0.05, ]

        if (nrow(sig_pathways) > 0) {
          all_pathways[[length(all_pathways) + 1]] <- sig_pathways
        }
      }
    }

    if (length(all_pathways) == 0) {
      cat("No significant pathways found for", label, db_name, "\n")
      next
    }

    # Combine unique pathways
    combined_pathways <- do.call(rbind, all_pathways)
    combined_pathways <- combined_pathways[!duplicated(combined_pathways$Term), ]
    combined_pathways <- combined_pathways[order(combined_pathways$Adjusted.P.value), ]
    combined_pathways <- combined_pathways[seq_len(min(40, nrow(combined_pathways))), ]

    # Create gene sets
    gsva_gene_sets <- list()
    for (i in seq_len(nrow(combined_pathways))) {
      pathway_name <- combined_pathways$Term[i]
      genes_vec <- unlist(strsplit(combined_pathways$Genes[i], ";"))

      if (db_name == "Reactome_2022") {
        pathway_display <- sub(" R-HSA-[0-9]+$", "", pathway_name)
      } else if (db_name == "GO_Biological_Process_2025") {
        pathway_display <- sub(" \\(GO:[0-9]+\\)$", "", pathway_name)
      } else {
        pathway_display <- pathway_name
      }

      gsva_gene_sets[[pathway_display]] <- genes_vec
    }

    # Run GSVA on all samples
    SetGSVAPar <- gsvaParam(
      exprData = as.matrix(expr_matrix[, annotation_df$Sample]),
      geneSets = gsva_gene_sets,
      minSize = 1, maxSize = 500,
      kcdf = "Gaussian"
    )
    gsva_scores <- gsva(SetGSVAPar, verbose = FALSE)

    # Order samples: EBV+ by SGCC, then EBV- by SGCC
    anno_ebv_pos <- annotation_df[annotation_df$EBV_status == "EBV+", ]
    anno_ebv_pos <- anno_ebv_pos[order(anno_ebv_pos$SGCC), ]

    anno_ebv_neg <- annotation_df[annotation_df$EBV_status == "EBV-", ]
    anno_ebv_neg <- anno_ebv_neg[order(anno_ebv_neg$SGCC), ]

    anno_ordered <- rbind(anno_ebv_pos, anno_ebv_neg)

    if (!"Macrophage_C1Q" %in% colnames(anno_ordered)) {
      anno_ordered$Macrophage_C1Q <- NA_real_
    }

    c1q_vals <- anno_ordered$Macrophage_C1Q
    c1q_rng <- range(c1q_vals, na.rm = TRUE)
    if (length(c1q_rng) != 2 || any(!is.finite(c1q_rng)) || c1q_rng[1] == c1q_rng[2]) {
      c1q_rng <- c(0, 1)
    }
    c1q_col_fun <- circlize::colorRamp2(c1q_rng, c("#f7f7f7", "#d95f02"))

    gsva_scores_ordered <- gsva_scores[, anno_ordered$Sample, drop = FALSE]
    gsva_scores_scaled <- t(scale(t(gsva_scores_ordered)))

    # Column annotation with split by EBV status
    col_annotation <- HeatmapAnnotation(
      EBV_Status = anno_ordered$EBV_status,
      Macrophage_C1Q = c1q_vals,
      SGCC_Score = anno_lines(anno_ordered$SGCC, which = "column",
                              gp = gpar(col = "red", lwd = 2)),
      col = list(
        EBV_Status = c("EBV+" = "#d95f02", "EBV-" = "#1b9e77"),
        Macrophage_C1Q = c1q_col_fun
      ),
      annotation_name_gp = gpar(fontsize = 8)
    )

    # Create column split based on EBV status
    column_split <- factor(anno_ordered$EBV_status, levels = c("EBV+", "EBV-"))

    ht_gsva <- Heatmap(
      gsva_scores_scaled,
      name = "GSVA\n(Z-score)",
      cluster_columns = FALSE,
      cluster_rows = TRUE,
      column_split = column_split,
      top_annotation = col_annotation,
      col = colorRamp2(c(-2, 0, 2), c("#b953a0", "#0f0b0e", "#efea1f")),
      row_names_gp = gpar(fontsize = 8),
      show_column_names = FALSE,
      column_title = paste0("GSVA: ", db_name, " (", label, " - All Samples, Ranked by SGCC within EBV Status)"),
      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
      column_gap = unit(3, "mm")
    )

    pdf(file.path(results_dir, paste0("CosMX4TMANew_", label, "_GSVA_", db_name, "_heatmap_AllSamples.pdf")),
        width = 16, height = 12)
    draw(ht_gsva, heatmap_legend_side = "right")
    dev.off()

    write.csv(gsva_scores,
              file.path(results_dir, paste0("CosMX4TMANew_", label, "_GSVA_", db_name, "_scores_AllSamples.csv")))

    cat("All-samples GSVA heatmap for", db_name, "saved\n")
  }
}

# Run GSVA for all samples together
run_gsva_all_samples(
  enrichr_ebv_pos_high, enrichr_ebv_pos_low,
  enrichr_ebv_neg_high, enrichr_ebv_neg_low,
  norm_expr_corrected, dfAnnotation,
  paste0(expression_label, "_", analysis_label),
  results_path
)

cat("\n=== DEG Analysis Complete ===\n")
cat("Results saved to:", results_path, "\n")

## ============================================================================
## EXPORT CELL-LEVEL METADATA WITH SGCC
## ============================================================================
cat("\n=== Exporting Cell-level Metadata with SGCC ===\n")

# Map window_id -> SGCC group (only extreme windows used in DEG)
deg_window_groups <- dfAnnotation %>%
  mutate(window_id = gsub(paste0("^", EXPRESSION_CELL_TYPE, "_"), "", Sample)) %>%
  select(window_id, SGCC_Group) %>%
  distinct()

cell_level_export <- mymeta %>%
  as.data.frame() %>%
  select(
    unique_id_fov,
    TMA,
    EBV_status,
    coreName,
    cellLabel,
    fov,
    Annotation_pooled,
    X_cent,
    Y_cent,
    CD3,Pax5,CD68,CD4,CD8,CD163,CD206, LMP1, LAG3,C1q,FoxP3,Myc,
    window_id
  ) %>%
  left_join(sgcc_scores_df %>% select(window_id, SGCC_score), by = "window_id") %>%
  left_join(deg_window_groups, by = "window_id") %>%
  mutate(
    SGCC_category = dplyr::case_when(
      SGCC_Group == "Low" ~ "SGCC low",
      SGCC_Group == "High" ~ "SGCC high",
      TRUE ~ "others"
    )
  ) %>%
  select(
    unique_id_fov,
    EBV_status,
    fov,
    coreName,
    cellLabel,
    Annotation_pooled,
    X_cent,
    Y_cent,
    SGCC_score,
    SGCC_category,
    CD3,Pax5,CD68,CD4,CD8,CD163,CD206, LMP1, LAG3,C1q,FoxP3,Myc,
    window_id
  )

cell_export_file <- file.path(results_path, paste0("CosMX4TMANew_", analysis_label, "_cell_metadata_with_SGCC.csv"))
write.csv(cell_level_export, cell_export_file, row.names = FALSE)
cat("Cell-level metadata saved to:", cell_export_file, "\n")

