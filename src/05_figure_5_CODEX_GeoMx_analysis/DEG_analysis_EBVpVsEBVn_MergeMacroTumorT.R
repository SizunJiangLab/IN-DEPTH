# set global path
wdpath <- "./data/DLBCL_run/"
loadingfunction_wdpath <- c("./src/")
# set env
library(dplyr)
library(Seurat)
library(readr)
library(tidyr)
library(tidyverse)
library(deldir)
library(igraph)
library(progress)
library(qs)
library(pheatmap)
library(epitools)  # For odds ratio calculation
library(SingleCellExperiment)
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","2-rank_SGCC_Impulse_preload_function.R"))
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","3-WithinCrossDomainTableAndVisualization_preload_function.R"))

# read in SGCC matrix 
SGCC_matrix <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","60_SGCC_matrix_DLBCL_mergeMacro_mergeTumor_mergeT.csv"),row.names = 1)

#### DEG analysis
# Load the ROI annotation data
annotation <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","DLBCL_ROIlevel_annotation.csv"))
# Standardize the annotations for CD8, CD4, and Tumor cells
annotation <- annotation %>%
  mutate(Annotation = case_when(
    grepl("CD8", Annotation) ~ "CD8T",
    grepl("CD4", Annotation) ~ "CD4T",
    grepl("Tumor", Annotation) ~ "Tumor",
    grepl("M1", Annotation) ~ "Macro",
    grepl("M2", Annotation) ~ "Macro",
    TRUE ~ Annotation
  ))

#### analysis pipeline
# Convert the matrix to a dataframe
SGCC_df <- as.data.frame(SGCC_matrix)
# Add a column for the cell type pairs
SGCC_df$CellTypePair <- rownames(SGCC_matrix)
# Exclude "Other" and "Neutrophil" annotations
annotation <- annotation %>%
  filter(!Annotation %in% c("Other", "Neutrophil"))
# Exclude cores containing "Tonsil" in their names
annotation <- annotation %>%
  filter(!grepl("Tonsil", coreName))
# Specify the selected cores
selected_cores <- c("Rochester_4", "Rochester_6",
                    "Rochester_7", "Rochester_9", "Rochester_11", "Rochester_12",
                    "Rochester_13", "Rochester_14",
                    "Rochester_15", "Rochester_16", "Rochester_17", "Rochester_18",
                    "Rochester_19", "Rochester_21", "Rochester_23",
                    "Rochester_25", "Rochester_TonsilA", "DFCI_2.2", "DFCI_3.2",
                    "DFCI_4.1", "DFCI_7.1", "DFCI_8.1",
                    "DFCI_12.1", "DFCI_13.2", "DFCI_14.1", "DFCI_15.2", "DFCI_17.1",
                    "DFCI_18.2", "DFCI_19.2", "DFCI_22.2", "DFCI_23.2")

# Filter for selected cores
annotation <- annotation %>%
  filter(coreName %in% selected_cores)
# summarize cell number
cellnumber_percore <- annotation %>%
  # Filter only the relevant cell types
  filter(Annotation %in% c("CD4T", "CD8T", "Macro", "Endothelial")) %>%
  # Group by coreName and MergedTumor to count the number of occurrences
  group_by(coreName, Annotation) %>%
  summarise(cell_count = n(), .groups = 'drop') %>% 
  mutate(ROI_CT = paste0(coreName,"_",Annotation))

################ Macro
# CD4 T new LogCPM_Seurat_top5000_k2_CD4T_vs_Endothelial_ruv_W1_ruv_W2_EBV_yes-no_CTnumber_TRUE
# macro: LogCPM_Seurat_top1000_k3_Macro_vs_Endothelial_ruv_W1_ruv_W2_ruv_W3_EBV_no_CTnumber_FALSE
# Tumor: LogCPM_Seurat_top1000_k3_Tumor_vs_Endothelial_ruv_W1_EBV_no_CTnumber_FALSE
# changed items for generating different results
my.seurat <- qread(file.path(wdpath,"GeoMX_batchcorrection","LogCPM-Seurat_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated.qs"))
spe_obj <- as.SingleCellExperiment(my.seurat$Seurat_top1000_k3)
spe_obj <- spe_obj[, spe_obj$ROI_rename %in% selected_cores]
convariate_component <- c("ruv_W1","ruv_W2","ruv_W3")
consider_CTnumber = F
covariate_matrix <- spe_obj@colData[,grep("ruv",colnames(spe_obj@colData))]
if(is.null(dim(covariate_matrix))){
  names(covariate_matrix) <- colnames(spe_obj)
}
# load previous calculated SGCC score and calculate gene in CT1
prepare_meta_output_CT1 <- prepare_meta_SGCC(spe_obj = spe_obj,
                                             SGCC_df = SGCC_matrix,
                                             cell_pair = "CD4T_Macro",
                                             CT.use = "Macro",
                                             condition_name = "EBV_Indicator",sample.drop = c("DFCI_9.1","DFCI_21.1","DFCI_11.1"),batch.factor = "Cohort")
dim(prepare_meta_output_CT1[[1]])
prepare_meta_output_CT1[[2]][1:5,1:5]
#"ruv_W1","ruv_W2","ruv_W3"
my.results_Fix_Time3 <- runDEGs_limma_fix_Time(meta_data = prepare_meta_output_CT1[[1]], 
                                               expression_matrix = prepare_meta_output_CT1[[2]], 
                                               consider_CTnumber = consider_CTnumber, 
                                               cellnumber_percore = cellnumber_percore,
                                               time_column = "Time", 
                                               EBV_column = "Condition",  # Column for EBV status
                                               time_fix = c(1,2,3), 
                                               logfc_threshold = 0.05, 
                                               p_val_threshold = 0.01,
                                               covariate_ID = convariate_component,
                                               covariate_mat = covariate_matrix)
rownames(my.results_Fix_Time3)[my.results_Fix_Time3$logFC < 0]
write.csv(my.results_Fix_Time3,
          file = paste0("/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/INDEPTH_data/DLBCL//pattern_analysis/Limma_DEG_overallEBV/",CT.use,"_EBV+ vs EBV-.csv"))

#
################ Macro
# based on the benchmarking results macro: LogCPM_Seurat_top1000_k3_Macro_vs_Endothelial_ruv_W1_ruv_W2_ruv_W3_EBV_no_CTnumber_FALSE
# changed items for generating different results
my.seurat <- qread(file.path(wdpath,"GeoMX_batchcorrection","LogCPM-Seurat_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated.qs"))
spe_obj <- as.SingleCellExperiment(my.seurat$Seurat_top1000_k3)
spe_obj <- spe_obj[, spe_obj$ROI_rename %in% selected_cores]
convariate_component <- c("ruv_W1","ruv_W2","ruv_W3")
consider_CTnumber = F
covariate_matrix <- spe_obj@colData[,grep("ruv",colnames(spe_obj@colData))]
if(is.null(dim(covariate_matrix))){
  names(covariate_matrix) <- colnames(spe_obj)
}
# load previous calculated SGCC score and calculate gene in CT1
prepare_meta_output_CT1 <- prepare_meta_SGCC(spe_obj = spe_obj,
                                             SGCC_df = SGCC_matrix,
                                             cell_pair = "CD4T_Macro",
                                             CT.use = "Macro",
                                             condition_name = "EBV_Indicator",sample.drop = c("DFCI_9.1","DFCI_21.1","DFCI_11.1"),batch.factor = "Cohort")
#"ruv_W1","ruv_W2","ruv_W3"
my.results <- runDEGs_limma_fix_Time(meta_data = prepare_meta_output_CT1[[1]], 
                                     expression_matrix = prepare_meta_output_CT1[[2]], 
                                     consider_CTnumber = consider_CTnumber, 
                                     cellnumber_percore = cellnumber_percore,
                                     time_column = "Time", 
                                     EBV_column = "Condition",  # Column for EBV status
                                     time_fix = c(1,2,3), 
                                     logfc_threshold = 0.05, 
                                     p_val_threshold = 0.01,
                                     covariate_ID = convariate_component,
                                     covariate_mat = covariate_matrix)
write.csv(my.results,
          file = file.path(wdpath,"GeoMX_differential_expressed_analysis","Macro_EBV+ vs EBV-.csv"))


################ CD4T
# based on the benchmarking results CD4 T:  LogCPM_Seurat_top5000_k2_CD4T_vs_Endothelial_ruv_W1_ruv_W2_EBV_yes-no_CTnumber_TRUE
# changed items for generating different results
spe_obj <- as.SingleCellExperiment(my.seurat$Seurat_top5000_k2)
spe_obj <- spe_obj[, spe_obj$ROI_rename %in% selected_cores]
convariate_component <- c("ruv_W1","ruv_W2")
consider_CTnumber = T
covariate_matrix <- spe_obj@colData[,grep("ruv",colnames(spe_obj@colData))]
if(is.null(dim(covariate_matrix))){
  names(covariate_matrix) <- colnames(spe_obj)
}
# load previous calculated SGCC score and calculate gene in CT1
prepare_meta_output_CT1 <- prepare_meta_SGCC(spe_obj = spe_obj,
                                             SGCC_df = SGCC_matrix,
                                             cell_pair = "CD4T_Macro",
                                             CT.use = "CD4T",
                                             condition_name = "EBV_Indicator",sample.drop = c("DFCI_9.1","DFCI_21.1","DFCI_11.1"),batch.factor = "Cohort")
#"ruv_W1","ruv_W2","ruv_W3"
my.results <- runDEGs_limma_fix_Time(meta_data = prepare_meta_output_CT1[[1]], 
                                     expression_matrix = prepare_meta_output_CT1[[2]], 
                                     consider_CTnumber = consider_CTnumber, 
                                     cellnumber_percore = cellnumber_percore,
                                     time_column = "Time", 
                                     EBV_column = "Condition",  # Column for EBV status
                                     time_fix = c(1,2,3), 
                                     logfc_threshold = 0.05, 
                                     p_val_threshold = 0.01,
                                     covariate_ID = convariate_component,
                                     covariate_mat = covariate_matrix)
write.csv(my.results,
          file = file.path(wdpath,"GeoMX_differential_expressed_analysis","CD4T_EBV+ vs EBV-.csv"))
################ Tumor
# Tumor: LogCPM_Seurat_top1000_k3_Tumor_vs_Endothelial_ruv_W1_EBV_no_CTnumber_FALSE
# changed items for generating different results
spe_obj <- as.SingleCellExperiment(my.seurat$Seurat_top1000_k3)
spe_obj <- spe_obj[, spe_obj$ROI_rename %in% selected_cores]
convariate_component <- c("ruv_W1")
consider_CTnumber = F
covariate_matrix <- spe_obj@colData[,grep("ruv",colnames(spe_obj@colData))]
if(is.null(dim(covariate_matrix))){
  names(covariate_matrix) <- colnames(spe_obj)
}
# load previous calculated SGCC score and calculate gene in CT1
prepare_meta_output_CT1 <- prepare_meta_SGCC(spe_obj = spe_obj,
                                             SGCC_df = SGCC_matrix,
                                             cell_pair = "Macro_Tumor",
                                             CT.use = "Tumor",
                                             condition_name = "EBV_Indicator",sample.drop = c("DFCI_9.1","DFCI_21.1","DFCI_11.1"),batch.factor = "Cohort")
my.results <- runDEGs_limma_fix_Time(meta_data = prepare_meta_output_CT1[[1]], 
                                     expression_matrix = prepare_meta_output_CT1[[2]], 
                                     consider_CTnumber = consider_CTnumber, 
                                     cellnumber_percore = cellnumber_percore,
                                     time_column = "Time", 
                                     EBV_column = "Condition",  # Column for EBV status
                                     time_fix = c(1,2,3), 
                                     logfc_threshold = 0.05, 
                                     p_val_threshold = 0.01,
                                     covariate_ID = convariate_component,
                                     covariate_mat = covariate_matrix)
write.csv(my.results,
          file = file.path(wdpath,"GeoMX_differential_expressed_analysis","Tumor_EBV+ vs EBV-.csv"))

