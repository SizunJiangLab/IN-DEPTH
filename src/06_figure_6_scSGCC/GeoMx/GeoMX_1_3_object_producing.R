# Load required libraries
library(Seurat)  # kept for downstream compatibility if needed
library(SpatialExperiment)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(qs)

# ============================================================================
# GeoMX-only object generation (Top 1)
# ============================================================================

# Input paths (aligned to GeoMX_1_2_batch_correction_benchmarking.R)
wdpath <- "/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/GeoMXNew_code_data_publish/Data/"
LogCPM_seurat_list <- qread(file.path(wdpath,"GeoMX_batchcorrection","LogCPM-Seurat_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated-update.qs"))
raw_obj <- qread(file.path(wdpath,"GeoMx_count_table/spe_for_DLBCL_mergeTumorTmacro_update.qs"))
cat("Loading GeoMX intermediate results...\n")
geomx_spe_obj <- LogCPM_seurat_list$Seurat_top1000_k3
raw_obj <- raw_obj[,colnames(geomx_spe_obj)]
exp_matrix <- as.matrix( raw_obj@assays@data$counts)
colnames(exp_matrix) <- colnames(geomx_spe_obj)
rownames(exp_matrix) <- rownames(geomx_spe_obj)
identical(colnames(exp_matrix),colnames(raw_obj))
geomx_sample_metadata <- as.data.frame(geomx_spe_obj@meta.data)
geomx_sample_metadata$sample_name <- rownames(geomx_sample_metadata)
geomx_sample_metadata$total_counts <- colSums(as.matrix(raw_obj@assays@data$counts))
geomx_sample_metadata$cell_type_clean <- geomx_sample_metadata$MergedLabel

colnames(raw_obj)[1:5]
colnames(geomx_spe_obj)[1:5]
# Output path (GeoMX only)
output_dir_objects <- file.path(wdpath, "Three_Cohorts_Top_Objects")
if (!dir.exists(output_dir_objects)) {
  dir.create(output_dir_objects, recursive = TRUE)
}

geomx_spe_obj <- SpatialExperiment(
  assay = list(counts = raw_obj@assays@data$counts,
               logcounts = geomx_spe_obj@assays$RNA@layers$counts),
  colData = geomx_spe_obj@meta.data)

geomx_save <- list(
  spe_obj = geomx_spe_obj,
  sample_metadata = geomx_sample_metadata,
  top_parameter_name = "Seurat_top1000_k3"
)
geomx_file <- file.path(output_dir_objects, "GeoMX_DFCI_Rochester_top1.qs")
qsave(geomx_save, file = geomx_file)
cat("GeoMX top 1 object saved to:", geomx_file, "\n")

