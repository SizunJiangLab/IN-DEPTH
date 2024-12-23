# set global path
wdpath <- "./data/DLBCL_run/"
loadingfunction_wdpath <- c("./src/")
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
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","3-WithinCrossDomainTableAndVisualization_preload_function.R"))
### changed variables
# claim your cell type pairs
Celltypepairs <- c("CD4T","Macro")
SGCC_pair <- "CD4T_Macro"
###
figure_save_path <- file.path("/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/INDEPTH_data/DLBCL/pattern_analysis/Limma_DEG_overallEBV//pathwayenrichment/",paste(Celltypepairs,collapse = '_'))

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
                    "Rochester_13", "Rochester_14", "Rochester_25",#"DFCI_14.1", "Rochester_TonsilA",
                    "Rochester_15", "Rochester_16", "Rochester_17", "Rochester_18",
                    "Rochester_19", "Rochester_21", "Rochester_23",
                     "DFCI_2.2", "DFCI_3.2",
                    "DFCI_4.1", "DFCI_7.1", "DFCI_8.1",
                    "DFCI_12.1", "DFCI_13.2",  "DFCI_15.2", "DFCI_17.1",
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

drop.samples <- setdiff(spe_obj_CT1$ROI_rename,selected_cores)
drop.samples
# load previous calculated SGCC score and calculate gene in CT1
# CD4T_Macro, Macro_Tumor
prepare_meta_output_CT1 <- prepare_meta_SGCC(spe_obj = spe_obj_CT1,
                                             SGCC_df = SGCC_matrix,
                                             cell_pair = c(SGCC_pair),
                                             CT.use = CT.1,
                                             condition_name = "EBV_Indicator",sample.drop = c(drop.samples),batch.factor = "Cohort")
prepare_meta_output_CT2 <- prepare_meta_SGCC(spe_obj = spe_obj_CT2,
                                             SGCC_df = SGCC_matrix,
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
Genelist_CT2 <- intersect(Macro_customized_M1M2,rownames(prepare_meta_output_CT2[[2]]))
Macropage_heatmap_pathway <- create_complex_heatmap_FixEBV(gene_list = NA,
                                                           pathway_matrix = gsva_results_CT2,scale_color_ht = 2,
                                                           expr_data = prepare_meta_output_CT2[[2]],
                                                           meta_data = metaoutput_CT2,
                                                           cluster_rows = T,cluster_columns = F
)

Genelist_CT1 <- intersect(CD4T_customized_Tcell_exhaustion, rownames(prepare_meta_output_CT1[[2]]))
CD4T_heatmap_pathway <- create_complex_heatmap_FixEBV(gene_list = NA,scale_color_ht = 2,
                                                      pathway_matrix = gsva_results_CT1,
                                                      expr_data = prepare_meta_output_CT1[[2]],
                                                      meta_data = metaoutput_CT1,
                                                      cluster_rows = T,cluster_columns = F
)
CD4T_heatmap_annotation_exhaustion_score <- create_complex_heatmap_annotation_from_pathway(pathway_matrix = gsva_results_CT1,scale_color_ht = 2,
                                                                                           meta_data = metaoutput_CT1,pathwayname = "CD4T_customized_Tcell_exhaustion")
# save plot
pdf("Figure5B-T dysfunction RNA annotation bar.pdf", width = 10,height = 9)
draw(CD4T_heatmap_annotation_exhaustion_score)
dev.off()
# save plot
pdf("Figures 5A and 5B and Supplementary Figure 7A heatmap.pdf", width = 10,height = 9)
draw(Macropage_heatmap_pathway)
dev.off()
# save plot
pdf("Figure 5B and Supplementary Figure 7B heatmap.pdf", width = 10,height = 9)
draw(CD4T_heatmap_pathway)
dev.off()

# Ternary plot
minmax_rescale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# scale data
transformed <- apply(as.matrix(SGCC_submatrix),1,FUN = minmax_rescale)
transformed <- transformed/rowSums(transformed)*100
transformed <- cbind.data.frame(transformed, t(SGCC_submatrix))
# rename
colnames(transformed) <- c("T_M","T_Tu","M_Tu","SGCC_T_M","SGCC_T_Tu","SGCC_M_Tu")
transformed$Core <- rownames(transformed)
transformed$ROI_CT = paste0(rownames(transformed),"_CD4T")
transformed <- transformed %>% 
  right_join(metaoutput_CT1,by = "ROI_CT")
# 
p1.merge <- ggtern(data =transformed , aes(T_M,T_Tu,M_Tu))+
  geom_point(aes(fill = EBV_status,shape = EBV_status),size = 3)+
  geom_density_tern(aes(color=EBV_status, alpha=..level..),bdl=0.02,bins=5)+
  scale_fill_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  )  +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  )+scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))


p1.merge.withname <- ggtern(data =transformed , aes(T_M,T_Tu,M_Tu))+
  geom_point(aes(fill = EBV_status,shape = EBV_status),size = 3)+
  geom_density_tern(aes(color=EBV_status, alpha=..level..),bdl=0.02,bins=5)+
  geom_text(aes(label = ROI_CT), size = 1) + # Add text labels
  scale_fill_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  )  +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  )+scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))





p_LAG3_CD4T <- ggtern(data = transformed, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = LAG3_CD4T,shape = EBV_status), 
    size = 3,     # Adjust point size as needed
    data = transformed
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))


p_Tox_CD4T <- ggtern(data = transformed, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = Tox_CD4T,shape = EBV_status), 
    size = 3,     # Adjust point size as needed
    data = transformed
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))

p_HLADR_macro <- ggtern(data = transformed, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = HLADR_macro,shape = EBV_status), 
    size = 3,     # Adjust point size as needed
    data = transformed
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))


p_PDL1_macro <- ggtern(data = transformed, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = PDL1_macro,shape = EBV_status), 
    size = 3,     # Adjust point size as needed
    data = transformed
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))


# read in AES results
CD4T_Macro <- qread(file.path(wdpath,"SGCC_relevant_analysis","Macro_CD4Tstat.result.qs"))
Tumor_Macro <- qread(file.path(wdpath,"SGCC_relevant_analysis","Macro_Tumorstat.result.qs"))
CD4T_Tumor <- qread(file.path(wdpath,"SGCC_relevant_analysis","CD4T_Tumorstat.result.qs"))

# rename
colnames(CD4T_Tumor) <- paste0(colnames(CD4T_Tumor),"_T_Tu")
CD4T_Tumor$Core <- rownames(CD4T_Tumor)
colnames(Tumor_Macro) <- paste0(colnames(Tumor_Macro),"_M_Tu")
Tumor_Macro$Core <- rownames(Tumor_Macro)
colnames(CD4T_Macro) <- paste0(colnames(CD4T_Macro),"_T_M")
CD4T_Macro$Core <- rownames(CD4T_Macro)
# add exhaustion_score 
sorted_CD4T_customized_Tcell_exhaustion <-  gsva_results_CT1["CD4T_customized_Tcell_exhaustion",][transformed$Sample]
# visualize AES score for CD4T Macro
transformed_CD4T_Macro_Tumor <- Tumor_Macro %>% 
  filter(Core %in% selected_cores) %>% 
  left_join(CD4T_Macro, by = "Core") %>% 
  left_join(CD4T_Tumor, by = "Core") %>% 
  left_join(transformed, by = "Core") 
  

transformed_CD4T_Macro_Tumor$CD4T_exhaustion_gene <- sorted_CD4T_customized_Tcell_exhaustion

p_CD4T_LMP1_filtered_mean <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = LMP1_filtered_mean,shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))

# plot 
p_CD4T_pd1 <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = PD1_CD4T,shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))

# plot 
p_CD4T_exhaustion_gene <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = CD4T_exhaustion_gene,shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))

#
p_CD4T_dysfunciton_protein <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = dysfunction_score_protein,shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))
# plot 
p_AES_M_Tu <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = log1p(AES_M_Tu),shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#2424db",high = "yellow")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))


p_AES_T_M <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = log1p(AES_T_M),shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#2424db",high = "yellow")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))


p_AES_CD4T_T <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = log1p(AES_T_Tu),shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#2424db",high = "yellow")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))

# plot 
p_LMPtumorpct <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = LMP1_pct,shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))

# plot 
p_M2PCT <- ggtern(data = transformed_CD4T_Macro_Tumor, aes(T_M, T_Tu, M_Tu)) +
  geom_density_tern(
    aes(x = T_M, y = T_Tu, z = M_Tu, color = EBV_status, alpha = ..level..),
    bdl = 0.02,
    bins = 5,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_color_manual(
    values = c("#ca6938", "#4a9d7a"),
    breaks = c("EBV+", "EBV-")
  ) +
  geom_point(
    mapping = aes(fill = m2_pct,shape = EBV_status), size = 3,
    data = transformed_CD4T_Macro_Tumor
  ) +
  scale_fill_gradient(low = "#fafafa",high = "#fc7051")+
  scale_shape_manual(values = c(21,24),breaks = c("EBV+", "EBV-"))

pdf("Figure 5C, 5D, supplementary Figure 7C.pdf", width = 30,height = 30)
grid.arrange( p1.merge,
              p1.merge.withname, 
              p_CD4T_exhaustion_gene,
              p_CD4T_dysfunciton_protein,
              p_CD4T_LMP1_filtered_mean,
              p_AES_T_M,
              p_AES_M_Tu,
              p_AES_CD4T_T,
              p_HLADR_macro,
              p_PDL1_macro,
              p_Tox_CD4T,
              p_LAG3_CD4T,
              p_CD4T_pd1,
              p_LMPtumorpct,
              p_M2PCT)
dev.off()


colnames(transformed_CD4T_Macro_Tumor)
functional_test_for_LMP1 <- c("dysfunction_score_protein","LMP1_filtered_mean","LMP1_g_tumor","LMP1_positive_numebr","m2_pct","m1_pct","CD4T_exhaustion_gene")

selected_data <- transformed_CD4T_Macro_Tumor %>% 
#  filter(EBV_status== "EBV-") %>% 
  select(functional_test_for_LMP1) 
correlation.tst <- 
  cor(selected_data, use = "pairwise.complete.obs", method = "spearman")
correlation.tst <- ifelse(is.na(correlation.tst), 0, correlation.tst)
diag(correlation.tst) <- 0
# If p-values are needed
correlation_pvalues <- matrix(NA, ncol = ncol(selected_data), nrow = ncol(selected_data))
rownames(correlation_pvalues) <- colnames(selected_data)
colnames(correlation_pvalues) <- colnames(selected_data)

# Perform pairwise correlation and generate scatter plots
for (i in 1:ncol(selected_data)) {
  for (j in 1:ncol(selected_data)) {
    if (i != j) {
      # Perform Spearman correlation test
      test_result <- cor.test(selected_data[[i]], selected_data[[j]], method = "spearman")
      correlation_pvalues[i, j] <- test_result$p.value
      if(correlation_pvalues[i, j] <=0.05){
        # Generate scatter plot with ggplot2
        plot_data <- data.frame(x = selected_data[[i]], y = selected_data[[j]], EBV = transformed_CD4T_Macro_Tumor$EBV_status)
        plot <- ggplot(plot_data, aes(x = x, y = y, color = EBV)) +
          geom_point(size =3) + # Semi-transparent points
          geom_smooth(method = "lm", color = "blue", se = FALSE) + # Line of best fit
          scale_color_manual(values = c("#ca6938", "#4a9d7a"),
                             breaks = c("EBV+", "EBV-"))+
          theme_classic() + # Minimal theme
          labs(
            title = paste(colnames(selected_data)[i], "vs", colnames(selected_data)[j]),
            x = NULL,
            y = NULL
          ) +
          theme(
            panel.background = element_rect(fill = "transparent", color = NA), # Transparent background
            plot.background = element_rect(fill = "transparent", color = NA),
            legend.background = element_rect(fill = "transparent", color = NA)
          )
        
        ggsave(
          filename = file.path(paste0("Supplementary Figure 7_scatterplot_", colnames(selected_data)[i], "_vs_", colnames(selected_data)[j], ".pdf")),
          plot = plot,device = "pdf",
          bg = "transparent", # Transparent background
          width = 8,
          height = 8,
          units = "in",
          dpi = 300
        )
      }
    }
  }
}


indicator_matrix <- ifelse(correlation_pvalues < 0.05 | is.na(correlation_pvalues), 1, 0.3)
heatmap_mat_indicator <- correlation.tst * indicator_matrix

alph_non_signaficant <- pheatmap(heatmap_mat_indicator,
         color = colorRamp2(c(-1,0, 1), c("#528ad0", "white", "#cb438d")),na_col = "#ffffff")
pheatmap(correlation_pvalues,display_numbers = TRUE,    number_format = "%.4f",
                                 color = colorRamp2(c(0, 1), c( "#cb438d","white")),na_col = "#ffffff")
pdf(file = "Supplementary Figure 7_LMP1_correlation_heatmap.pdf",width = 10,height = 10)
print(alph_non_signaficant)
dev.off()


