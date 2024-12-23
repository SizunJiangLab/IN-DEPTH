wdpath <- "./data/DLBCL_run/"
loadingfunction_wdpath <- c("./src/")
# set env
library(ggtern)
library(dplyr)
library(Seurat)
library(GSVA)
library(readr)
library(tidyr)
library(tidyverse)
library(deldir)
library(igraph)
library(progress)
library(qs)
library(ggpubr)
library(SingleCellExperiment)
library(ggrepel)
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","5-Preload_customized_gene_list.R"))
# when running 4-post_DEG_heatmap_visualization_preloadfunction-mergeMacroTumorT.R, you need to setwd() to the directory contain data downloaded from zenodo.
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","4-post_DEG_heatmap_visualization_preloadfunction-mergeMacroTumorT.R"))
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","2-rank_SGCC_Impulse_preload_function.R"))
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","3-WithinCrossDomainTableAndVisualization_preload_function.R"))### changed variables
# claim your cell type pairs
Celltypepairs <- c("Tumor","Macro")
SGCC_pair <- "Macro_Tumor"
# assign cell types
CT.1 <-Celltypepairs[[1]]
CT.2 <-Celltypepairs[[2]]
CT.1.exp <- list_exp_matrix_metadata[[CT.1]]
CT.2.exp <- list_exp_matrix_metadata[[CT.2]]

spe_obj_CT1 <- SingleCellExperiment(assays = list(counts = as.matrix(CT.1.exp),
                                                  logcounts = as.matrix(CT.1.exp)),
                                    colData = list_exp_matrix_metadata$reformedMeta)
spe_obj_CT2 <- SingleCellExperiment(assays = list(counts = as.matrix(CT.2.exp),
                                                  logcounts = as.matrix(CT.2.exp)),
                                    colData = list_exp_matrix_metadata$reformedMeta)


#
SGCC_matrix <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","60_SGCC_matrix_DLBCL_mergeMacro_mergeTumor_mergeT.csv"),row.names = 1)
recalled_CTs <- combn(c("CD4T","Macro","Tumor"), 2, FUN = function(x) paste(x, collapse = "_"))
SGCC_submatrix <- SGCC_matrix[recalled_CTs,]
SGCC_submatrix
# Load the ROI annotation data
annotation <- list_exp_matrix_metadata$annotation
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

# Exclude "Other" and "Neutrophil" annotations
annotation <- annotation %>%
  filter(!Annotation %in% c("Other", "Neutrophil"))

# Exclude cores containing "Tonsil" in their names
annotation <- annotation %>%
  filter(!grepl("Tonsil", coreName))
unique(annotation$coreName)
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

# select by cell number difference
filtered_summary <- filter_samples_summary(annotation, 
                                           max_diff = 1, 
                                           cell_types = c(Celltypepairs[[1]], Celltypepairs[[2]]), 
                                           tumor_annotation = "Tumor",  # Specify the annotation used to identify tumor cells
                                           tumor_upper = 1, tumor_lower =0 ,logfc =10)

drop.samples <- intersect(selected_cores, filtered_summary$coreName[filtered_summary$filtered_out ==T])
drop.samples
# load previous calculated SGCC score and calculate gene in CT1
# CD4T_Macro, Macro_Tumor
prepare_meta_output_CT1 <- prepare_meta_SGCC(spe_obj = spe_obj_CT1,
                                             SGCC_df = SGCC_matrix,
                                             # cell_pair = paste0(Celltypepairs[[1]],"_",Celltypepairs[[2]]),
                                             cell_pair = c(SGCC_pair),
                                             CT.use = CT.1,
                                             condition_name = "EBV_Indicator",sample.drop = c(drop.samples),batch.factor = "Cohort")
prepare_meta_output_CT2 <- prepare_meta_SGCC(spe_obj = spe_obj_CT2,
                                             SGCC_df = SGCC_matrix,
                                            # cell_pair = paste0(Celltypepairs[[1]],"_",Celltypepairs[[2]]),
                                            cell_pair = c(SGCC_pair),
                                             CT.use = CT.2,
                                             condition_name = "EBV_Indicator",sample.drop = c(drop.samples),batch.factor = "Cohort")


metaoutput_CT1 <- spe_obj_CT1@colData %>% 
  as.data.frame() %>% 
  select(ROI_CT, virus_loading, virus_loading_InTumor, total_cells, total_tumor_cells, 
         Tumor_other_cell, Tumor_BCL2_cell, Tumor_BCL6_cell, LMP1_g_tumor,
         Tumor_Myc_cell,LMP1_filtered_mean,LMP1_positive_numebr, LMP1_pct, 
         m1_pct,m2_pct,m1_odds, m2_odds,PD1_CD4T, PDL1_macro, Tox_CD4T, LAG3_CD4T, HLADR_macro,dysfunction_score_protein) %>% 
  right_join(prepare_meta_output_CT1[[1]], by = "ROI_CT") %>% 
  mutate(SGCC_cat = case_when(
    Time == 1 ~ "SGCC low",
    Time == 2 ~ "SGCC mediate",
    Time == 3 ~ "SGCC high",
    TRUE ~ NA_character_  # Set default to NA if Time doesn't match any specified value
  ),
  EBV_status = case_when(
    Condition == "case" ~ "EBV+",
    Condition == "control" ~ "EBV-"
  )
)

metaoutput_CT2 <- spe_obj_CT2@colData %>% 
  as.data.frame() %>% 
  select(ROI_CT, virus_loading, virus_loading_InTumor, total_cells, total_tumor_cells, 
         Tumor_other_cell, Tumor_BCL2_cell, Tumor_BCL6_cell, LMP1_g_tumor,
         Tumor_Myc_cell,LMP1_filtered_mean,LMP1_positive_numebr, LMP1_pct, 
         m1_pct,m2_pct,m1_odds,m2_odds, PD1_CD4T, PDL1_macro, Tox_CD4T, LAG3_CD4T, HLADR_macro,dysfunction_score_protein) %>% 
  right_join(prepare_meta_output_CT2[[1]], by = "ROI_CT") %>% 
  mutate(SGCC_cat = case_when(
    Time == 1 ~ "SGCC low",
    Time == 2 ~ "SGCC mediate",
    Time == 3 ~ "SGCC high",
    TRUE ~ NA_character_  # Set default to NA if Time doesn't match any specified value
  ),
  EBV_status = case_when(
    Condition == "case" ~ "EBV+",
    Condition == "control" ~ "EBV-"
  ))

# read in DEG
Celltypepairs_con <- paste0(Celltypepairs,collapse = "-")
# run GSVA
SetGSVAPar_CT1 <- gsvaParam(exprData = prepare_meta_output_CT1[[2]], 
                            geneSets = functionalenrichment[grep(paste0("^",Celltypepairs[[1]]),(names(functionalenrichment)))],
                        minSize=3, maxSize=400,kcdf = "Gaussian")
gsva_results_CT1 <- gsva(SetGSVAPar_CT1, verbose = TRUE)
GSVA_pathwayname_CT1 <- rownames(gsva_results_CT1)
SetGSVAPar_CT2 <- gsvaParam(exprData = prepare_meta_output_CT2[[2]], 
                            geneSets = functionalenrichment[grep(paste0("^",Celltypepairs[[2]]),(names(functionalenrichment)))],
                            minSize=3, maxSize=100,kcdf = "Gaussian")
gsva_results_CT2 <- gsva(SetGSVAPar_CT2, verbose = TRUE)
GSVA_pathwayname_CT2 <- rownames(gsva_results_CT2)
# 
# Genelist_CT2 <- grep(paste0("^",Celltypepairs[[2]]), names(functionalenrichment)[!names(functionalenrichment) %in% GSVA_pathwayname_CT2],value = T)
Genelist_CT2 <- intersect(Macro_customized_M1M2,rownames(prepare_meta_output_CT2[[2]]))

Macropage_heatmap_pathway <- create_complex_heatmap_FixEBV(gene_list = NA,scale_color_ht = 2,
                                                        pathway_matrix = gsva_results_CT2,
                                                        expr_data = prepare_meta_output_CT2[[2]],
                                                        meta_data = metaoutput_CT2,
                                                        cluster_rows = T,cluster_columns = F
)

# Genelist_CT1 <- grep(paste0("^",Celltypepairs[[1]]), names(functionalenrichment)[!names(functionalenrichment) %in% GSVA_pathwayname_CT1],value = T)
Genelist_CT1 <- intersect(Tumor_customized_M1M2_formation,rownames(prepare_meta_output_CT1[[2]]))
Tumor_heatmap_gene <- create_complex_heatmap_FixEBV(gene_list = Genelist_CT1,scale_color_ht = 2,
                                                   pathway_matrix = NA,
                                                   expr_data = prepare_meta_output_CT1[[2]],
                                                   meta_data = metaoutput_CT1,
                                                   cluster_rows = T,cluster_columns = F
)

Tumor_heatmap_pathway <- create_complex_heatmap_FixEBV(gene_list = NA,scale_color_ht = 2,
                                                   pathway_matrix = gsva_results_CT1,
                                                   expr_data = prepare_meta_output_CT1[[2]],
                                                   meta_data = metaoutput_CT1,
                                                   cluster_rows = T,cluster_columns = F
)
# save plot
pdf("Figures 5A and 5B and Supplementary Figure 7A heatmap.pdf", width = 10,height = 9)
draw(Macropage_heatmap_pathway)
dev.off()

# save plot
pdf("Supplementary Figure 5A heatmap.pdf", width = 10,height = 9)
draw(Tumor_heatmap_EBV_gene)
dev.off()
# save plot
pdf("Figures 5A and 5B and Supplementary Figure 7A heatmap.pdf", width = 10,height = 9)
draw(Tumor_heatmap_pathway)
dev.off()





