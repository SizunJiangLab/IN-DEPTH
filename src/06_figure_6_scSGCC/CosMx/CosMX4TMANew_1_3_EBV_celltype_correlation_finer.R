# New CosMX EBV, Cell Type, and Correlation Analysis (Finer Cell Type)
# Script 1_3: Check EBV gene expression, cell type specific genes, and Pearson correlation
# Load processed data from New_CosMX4TMA scripts 1_1_finer and 1_2_2_finer

# Load required libraries
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(qs)
library(GSVA)

# Set paths
base_root <- "/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/CosMXNew_code_data_publish/"
output_base <- file.path(base_root, "Data/CosMXNew/Exampleoutpu")
input_intermediate_dir <- file.path(output_base, "CosMX4TMANew_1_2_2_batch_correction_benchmarking_finer")
results_path <- file.path(output_base, "CosMX4TMANew_1_3_EBV_celltype_correlation_finer")

cat("=== Script New_1_3: EBV, Cell Type, and Correlation Analysis (FINER Cell Types) ===\n")

# Load batch-corrected results from CosMX4TMANew_1_2_2_finer script
cat("Loading batch-corrected results from CosMX4TMANew_1_2_2_batch_correction_benchmarking_finer.R...\n")

# These files are created by CosMX4TMANew_1_2_2_batch_correction_benchmarking_finer.R
New_scores_path <- file.path(input_intermediate_dir, "New_parameter_search_scores_TMA-DFCI_Rochester_finer_subsetTMacroTumor.csv")
New_spe_path <- file.path(input_intermediate_dir, "New_parameter_search_spe_lists_TMA-DFCI_Rochester_finer_subsetTMacroTumor.qs")

# Check if files exist
if (!file.exists(New_scores_path)) {
  stop("Batch correction scores not found: ", New_scores_path,
       "\nPlease run CosMX4TMANew_1_2_2_batch_correction_benchmarking_finer.R first.")
}
if (!file.exists(New_spe_path)) {
  stop("Batch correction SPE lists not found: ", New_spe_path,
       "\nPlease run CosMX4TMANew_1_2_2_batch_correction_benchmarking_finer.R first.")
}
# Create results directory if it doesn't exist
if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
}

# Load parameter scores and SPE lists from batch correction benchmarking
New_param_scores <- read.csv(New_scores_path)
New_spe_lists <- qread(New_spe_path)

cat("Loaded batch correction results:\n")
cat("- Number of parameter combinations:", nrow(New_param_scores), "\n")
cat("- Number of SPE objects:", length(New_spe_lists), "\n")

# Select top-performing batch correction based on overall score
New_param_scores <- New_param_scores[order(-New_param_scores$overall), ]
New_top5 <- head(New_param_scores$name, 5)
New_top5 <- New_top5[New_top5 %in% names(New_spe_lists)]
New_top5_spe <- New_spe_lists[New_top5]
spe_obj <- New_top5_spe[[1]]

cat("Using top batch-corrected dataset:", names(New_top5_spe)[1], "\n")
cat("Top 5 parameter combinations:\n")
print(head(New_param_scores[, c("name", "overall", "bio_mean", "batch_mean")], 5))
cat("SPE object dimensions:", nrow(spe_obj), "genes x", ncol(spe_obj), "samples\n")

# Extract metadata
sample_metadata <- as.data.frame(colData(spe_obj))
sample_metadata$sample_name <- rownames(sample_metadata)

cat("\nFiner cell types in dataset:\n")
print(table(sample_metadata$cell_type))

# Load marker genes list (assuming it exists from the original pipeline)
marker_genes_path <- file.path(base_root,"Data","marker_genes_list.qs")

if (file.exists(marker_genes_path)) {
  marker_genes_list <- qread(marker_genes_path)

  # Extract marker gene lists
  CD4T_markers <- marker_genes_list$CD4T_markers
  CD8T_markers <- marker_genes_list$CD8T_markers
  B_cellmarkers <- marker_genes_list$B_cellmarkers
  Macrophage_marker_list <- marker_genes_list$Macrophage_marker_list
  TAM_gene_clusters <- marker_genes_list$TAM_gene_clusters
  CD4T_panelgene <- marker_genes_list$CD4T_panelgene
  CD8T_panelgene <- marker_genes_list$CD8T_panelgene
  Tumor_genes_list <- marker_genes_list$Tumor_genes_list

  cat("Marker genes loaded successfully\n")
} else {
  cat("Warning: marker_genes_list.qs not found. Using default marker genes.\n")
  # Define default marker genes
  CD4T_markers <- c("CD3D", "CD3E", "CD3G", "CD4", "CD28")
  CD8T_markers <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B")
  B_cellmarkers <- c("MS4A1", "CD19", "CD79A", "PAX5")
  Macrophage_marker_list <- list(c("CD68", "CD163", "CD14"))
  TAM_gene_clusters <- list(c("CD68", "CD163", "MMP9"))
  CD4T_panelgene <- list(CD4T_markers)
  CD8T_panelgene <- list(CD8T_markers)
  Tumor_genes_list <- list(B_cellmarkers)
}

cat("Data loaded successfully:\n")
cat("- SPE object:", ncol(spe_obj), "samples,", nrow(spe_obj), "genes\n")
cat("- Sample metadata:", nrow(sample_metadata), "samples\n")

# Save the batch-corrected SPE object for future use
qsave(spe_obj, file.path(results_path, "New_CosMX_Top1_BatchCorrected_finer_subsetTMacroTumor.qs"))
cat("Saved batch-corrected SPE object to: New_CosMX_Top1_BatchCorrected_finer.qs\n")

# ============================================================================
# HELPER FUNCTION: Create comprehensive heatmap annotation
# ============================================================================

create_heatmap_annotation <- function(metadata,
                                      include_barplots = TRUE,
                                      barplot_vars = c("total_counts", "cell_count"),
                                      lineage_markers = NULL,
                                      ebv_col = "EBV_status",
                                      celltype_col = "cell_type",
                                      barplot_height = 1) {

  # Initialize annotation list
  ann_list <- list()
  col_maps <- list()

  # Add EBV status
  if (ebv_col %in% colnames(metadata)) {
    ann_list$EBV_Status <- metadata[[ebv_col]]
    col_maps$EBV_Status <- c("EBV+" = "#d95f02", "EBV-" = "#1b9e77", "unknown" = "grey")
  }

  # Add cell type (finer annotations - use automatic colors for many categories)
  if (celltype_col %in% colnames(metadata)) {
    ann_list$Cell_Type <- metadata[[celltype_col]]
    # For finer cell types, we may have many categories, so we'll let ComplexHeatmap auto-assign colors
    # Or create a custom palette if needed
    unique_celltypes <- unique(metadata[[celltype_col]])
    if (length(unique_celltypes) <= 20) {
      # Use a colorblind-friendly palette for up to 20 categories
      library(RColorBrewer)
      n_colors <- length(unique_celltypes)
      if (n_colors <= 8) {
        colors <- brewer.pal(max(3, n_colors), "Set2")
      } else if (n_colors <= 12) {
        colors <- brewer.pal(n_colors, "Set3")
      } else {
        colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
      }
      col_maps$Cell_Type <- setNames(colors[1:n_colors], unique_celltypes)
    }
    # If more than 20, let ComplexHeatmap auto-assign
  }

  # Add bar plots if requested
  if (include_barplots) {
    barplot_colors <- c(
      total_counts = "#e78ac3",
      cell_count = "#fc8d62",
      Virus_count = "#8da0cb"
    )

    for (var in barplot_vars) {
      if (var %in% colnames(metadata)) {
        values <- metadata[[var]]
        # Skip if all NA
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
  }

  # Add lineage markers as continuous scales (with scaling)
  if (!is.null(lineage_markers)) {
    present_markers <- intersect(lineage_markers, colnames(metadata))
    for (mk in present_markers) {
      marker_vals <- as.numeric(metadata[[mk]])
      if (!all(is.na(marker_vals))) {
        # Scale the marker values (z-score)
        marker_vals_scaled <- scale(marker_vals)
        marker_vals_scaled[is.na(marker_vals_scaled)] <- 0  # Handle NAs after scaling

        ann_list[[mk]] <- as.numeric(marker_vals_scaled)

        # Use fixed scale for all markers: -2 to +2 (typical z-score range)
        col_maps[[mk]] <- circlize::colorRamp2(
          c(-2, 0, 2),
          c("#2166ac", "#f7f7f7", "#b2182b")
        )
      }
    }
  }

  # Create the annotation object
  col_annotation <- do.call(HeatmapAnnotation,
                            c(ann_list,
                              list(col = col_maps,
                                   annotation_name_gp = gpar(fontsize = 10))))

  return(col_annotation)
}

# ============================================================================
# CHECK EBV GENE EXPRESSION (focus on B cell-related finer types)
# ============================================================================

cat("\n=== Analyzing EBV Gene Expression (Finer Cell Types) ===\n")

virus_genes_cosmx <- c(grep("^LMP", rownames(spe_obj), value = TRUE),
                       grep("^BCRF1", rownames(spe_obj), value = TRUE),
                       grep("^BNLF", rownames(spe_obj), value = TRUE),
                       grep("^BLLF1", rownames(spe_obj), value = TRUE),
                       grep("^BZL", rownames(spe_obj), value = TRUE),
                       grep("^EBNA", rownames(spe_obj), value = TRUE),
                       grep("^RPM", rownames(spe_obj), value = TRUE))

cat("EBV genes found:", length(virus_genes_cosmx), "\n")
if (length(virus_genes_cosmx) > 0) {
  cat("EBV genes:", paste(virus_genes_cosmx, collapse = ", "), "\n")
}

if (length(virus_genes_cosmx) > 0) {
  # For finer cell types, identify B cell-related subtypes
  # This might include "B cell", "Germinal center B cell", "Memory B cell", etc.
  all_cell_types <- unique(sample_metadata$cell_type)
  b_cell_related <- grep("B|Plasma", all_cell_types, value = TRUE, ignore.case = TRUE)

  cat("B cell-related finer types found:", paste(b_cell_related, collapse = ", "), "\n")

  if (length(b_cell_related) > 0) {
    b_cell_samples <- sample_metadata$cell_type %in% b_cell_related

    EBV_counts <- colSums(assay(spe_obj, "counts")[virus_genes_cosmx, b_cell_samples])
    EBV_expr <- assay(spe_obj, "logcounts")[virus_genes_cosmx, b_cell_samples]
    EBV_meta <- sample_metadata[b_cell_samples, ]
    EBV_meta$Virus_count <- EBV_counts

    # Add total counts for bar plot
    raw_counts_ebv <- assay(spe_obj, 1)[, rownames(EBV_meta), drop = FALSE]
    EBV_meta$total_counts <- colSums(as.matrix(raw_counts_ebv))

    cat("EBV expression analysis:\n")
    cat("- B cell-related samples:", ncol(EBV_expr), "\n")
    cat("- EBV genes with expression:", sum(rowSums(EBV_expr) > 0), "\n")

    # Create column annotation for EBV heatmap
    col_annotation_ebv <- create_heatmap_annotation(
      metadata = EBV_meta,
      include_barplots = TRUE,
      barplot_vars = c("total_counts", "cell_count", "Virus_count"),
      lineage_markers = NULL,
      ebv_col = "EBV_status",
      celltype_col = "cell_type",
      barplot_height = 1
    )

    # Create column split based on EBV status (EBV+ on left, EBV- on right)
    column_split_ebv <- factor(EBV_meta$EBV_status, levels = c("EBV+", "EBV-"))

    # Create heatmap
    EBV_heatmap_pathway <- Heatmap(
      t(scale(t(EBV_expr))),
      name = "EBV gene",
      col = colorRamp2(c(-2, 0, 2), c("#b953a0", "black", "#f4ed17")),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      column_split = column_split_ebv,
      top_annotation = col_annotation_ebv,
      show_column_names = TRUE,
      show_row_names = TRUE,
      column_names_rot = 45,
      row_names_gp = gpar(fontsize = 8),
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      heatmap_legend_param = list(
        title = "Scaled Expression",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
      )
    )

    # Display heatmap
    cat("Displaying EBV gene expression heatmap (Finer B cell types)...\n")
    draw(EBV_heatmap_pathway, padding = unit(c(35, 10, 5, 45), "mm"))
  } else {
    cat("No B cell-related finer cell types found\n")
  }
} else {
  cat("No EBV genes found in the dataset\n")
}

# ============================================================================
# CHECK CELL TYPE SPECIFIC GENES (using coarse markers with finer types)
# ============================================================================

cat("\n=== Analyzing Cell Type Specific Genes (Finer Cell Types) ===\n")

CD4Tcell_genelist <- intersect(unique(unlist(c(CD4T_panelgene, CD4T_markers))), rownames(spe_obj))
CD8Tcell_genelist <- intersect(unique(unlist(c(CD8T_panelgene, CD8T_markers))), rownames(spe_obj))
Macrophage_genelist <- intersect(unique(unlist(c(TAM_gene_clusters, unlist(Macrophage_marker_list)))), rownames(spe_obj))
Tumor_genelist <- intersect(union(unique(unlist(Tumor_genes_list)), B_cellmarkers), rownames(spe_obj))

cat("Cell type specific gene lists:\n")
cat("- CD4T genes:", length(CD4Tcell_genelist), "\n")
cat("- CD8T genes:", length(CD8Tcell_genelist), "\n")
cat("- Macrophage genes:", length(Macrophage_genelist), "\n")
cat("- B cell/Tumor genes:", length(Tumor_genelist), "\n")

# Use CD4T markers for demonstration (can be changed to other cell types)
demo_genes <- c(intersect(CD4T_markers, rownames(spe_obj)))
if (length(demo_genes) > 0) {
  all_expr <- assay(spe_obj, "logcounts")[demo_genes, ]
  all_meta <- sample_metadata

  # Add total counts for bar plot
  raw_counts_all <- assay(spe_obj, 1)[, rownames(all_meta), drop = FALSE]
  all_meta$total_counts <- colSums(as.matrix(raw_counts_all))

  cat("Cell type specific gene expression analysis:\n")
  cat("- Genes analyzed:", nrow(all_expr), "\n")
  cat("- Samples analyzed:", ncol(all_expr), "\n")

  # Create column annotation for cell type heatmap
  col_annotation_all <- create_heatmap_annotation(
    metadata = all_meta,
    include_barplots = TRUE,
    barplot_vars = c("total_counts", "cell_count"),
    lineage_markers = NULL,
    ebv_col = "EBV_status",
    celltype_col = "cell_type",
    barplot_height = 1
  )

  # Create column split based on EBV status
  column_split_all <- factor(all_meta$EBV_status, levels = c("EBV+", "EBV-"))

  # Create heatmap
  all_heatmap_pathway <- Heatmap(
    t(scale(t(all_expr))),
    name = "CTS gene",
    col = colorRamp2(c(-2, 0, 2), c("#b953a0", "black", "#f4ed17")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_split = column_split_all,
    top_annotation = col_annotation_all,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Scaled Expression",
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
  )

  # Display heatmap
  cat("Displaying cell type specific gene expression heatmap (Finer types)...\n")
  draw(all_heatmap_pathway, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
}


