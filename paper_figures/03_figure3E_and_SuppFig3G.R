# set env
library(ComplexHeatmap)
library(circlize)
library(qs)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
## set working directory
wdpath <- "./data/Tonsil_run/"
# L24
L24_norm_mat_B <- qread(file.path(wdpath,"SGCC_relevant_analysis","norm_mat_B_L24.qs"))
L24_norm_mat_T <- qread(file.path(wdpath,"SGCC_relevant_analysis","norm_mat_T_L24.qs"))
L24_matrix_merge <- qread(file.path(wdpath,"SGCC_relevant_analysis","matrix_merge_L24.qs"))
L24_SGCC_matrix <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","SGCC_matrix_L24.csv"),row.names = 1)
L24_CODEX <- read.csv(file.path(wdpath,"Annotation/L24_final_annotation_ROIlevel.csv"))
# change names
mapping_table <- data.frame(
  SegmentLabel = c("BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "DC", "Endo", "Full ROI", "M1", "M2", "Myeloid", "Other"),
  Annotation5 = c("BCL6- B Cell", "BCL6+ B Cell", "CD4 T", "CD4 Treg", "CD8 T", "DC", "Endothelial", "Full ROI", "M1", "M2", "Myeloid", "Other")
)
L24_CODEX <- L24_CODEX %>%
  left_join(mapping_table, by = c("Annotation5" = "Annotation5"))

# L25
L25_norm_mat_B <- qread(file.path(wdpath,"SGCC_relevant_analysis","norm_mat_B_L25.qs"))
L25_norm_mat_T <- qread(file.path(wdpath,"SGCC_relevant_analysis","norm_mat_T_L25.qs"))
L25_matrix_merge <- qread(file.path(wdpath,"SGCC_relevant_analysis","matrix_merge_L25.qs"))
L25_SGCC_matrix <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","SGCC_matrix_L25.csv"),row.names = 1)
L25_CODEX <- read.csv(file.path(wdpath,"Annotation/L25_final_annotation_ROIlevel.csv"))
# change names
mapping_table <- data.frame(
  SegmentLabel = c("BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "DC", "Endo", "Full ROI", "M1", "M2", "Myeloid", "Other"),
  Annotation7 = c("BCL6- B Cell", "BCL6+ B Cell", "CD4 T", "CD4 Treg", "CD8 T", "DC", "Endothelial", "Full ROI", "M1", "M2", "Myeloid", "Other")
)
L25_CODEX <- L25_CODEX %>%
  left_join(mapping_table, by = c("Annotation7" = "Annotation7"))

# rename 
colnames(L24_norm_mat_B) <- paste0("L24_",colnames(L24_norm_mat_B))
colnames(L24_norm_mat_T) <- paste0("L24_",colnames(L24_norm_mat_T))
colnames(L24_matrix_merge) <- paste0("L24_",colnames(L24_matrix_merge))
colnames(L24_SGCC_matrix) <- paste0("L24_",colnames(L24_SGCC_matrix))
colnames(L25_norm_mat_B) <- paste0("L25_",colnames(L25_norm_mat_B))
colnames(L25_norm_mat_T) <- paste0("L25_",colnames(L25_norm_mat_T))
colnames(L25_matrix_merge) <- paste0("L25_",colnames(L25_matrix_merge))
colnames(L25_SGCC_matrix) <- paste0("L25_",colnames(L25_SGCC_matrix))
# plot select heatmap 
gene_selected <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","BCL6pB_CD4T_Gene_select_from_Impulse_prediction.csv"))
L24_L25_SGCC <- c(L24_SGCC_matrix["BCL6pB_CD4T",],L25_SGCC_matrix["BCL6pB_CD4T",])
SGCC_reorder.name <- names(L24_L25_SGCC)
SGCC_reorder.name <- SGCC_reorder.name[order(as.numeric(L24_L25_SGCC))]
L24_L25_SGCC_order <- as.numeric(L24_L25_SGCC[SGCC_reorder.name])
L24_L25_matrix_merge <- cbind(L25_matrix_merge,L24_matrix_merge)
L24_L25_matrix_merge <- L24_L25_matrix_merge[,SGCC_reorder.name]
#
# Use grid graphics to create density plots for the rows and columns
column_density_plot <- HeatmapAnnotation(
  SGCC = anno_lines(L24_L25_SGCC_order,which = "column"),
  gp = gpar(lwd = 10),
  height = unit(1.3, "cm")
)
`%notin%` <- Negate(`%in%`)
scale_mat <- t(scale(t(as.matrix(L24_L25_matrix_merge))))
ht_list <- Heatmap(scale_mat[rownames(scale_mat)%notin% c("MKI67","CDK1","CDT1","CCNB1") ,],
                   cluster_columns = F,
                   cluster_rows = F, 
                   row_split = c(rep("A",3),
                                 rep("B",1),
                                 rep("C",2),
                                 rep("D",6),
                                 rep("E",2),
                                 # rep("F",4),
                                 rep("G",2)),
                   row_names_gp = gpar(fontsize = 10),top_annotation = column_density_plot,
                   col = colorRamp2(c(min(t(scale(t(as.matrix(scale_mat)))))+1, 
                                      0, max(t(scale(t(as.matrix(scale_mat)))))-1), 
                                    c("#b9529f", "#0e0a0c", "#f3ec19")))
svg("Figure E gene expression.svg", width = 10,height = 30)
draw(ht_list)
dev.off()
# Vector for B cell pathways
msigdb_species <- "Homo sapiens"
b_cell_gene_sets <- msigdbr(species = msigdb_species, category = "C5", subcategory = "GO:BP") %>%
  dplyr::filter(gs_name %in% c(
    "GOBP_B_CELL_PROLIFERATION",
    "GOBP_B_CELL_PROLIFERATION_INVOLVED_IN_IMMUNE_RESPONSE",
    "GOBP_GERMINAL_CENTER_B_CELL_DIFFERENTIATION",
    "GOBP_MATURE_B_CELL_DIFFERENTIATION",
    "GOBP_B_CELL_DIFFERENTIATION",
    "GOBP_B_CELL_ACTIVATION",
    "GOBP_INTERLEUKIN_6_MEDIATED_SIGNALING_PATHWAY",
    "GOBP_MHC_CLASS_II_BIOSYNTHETIC_PROCESS",
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
    "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX",
    "GOBP_DNA_REPLICATION",
    "GOBP_DNA_TOPOLOGICAL_CHANGE",
    "GOBP_DNA_LIGATION"
  ))

# Vector for T cell pathways
t_cell_gene_sets <- msigdbr(species = msigdb_species, category = "C5", subcategory = "GO:BP") %>%
  dplyr::filter(gs_name %in% c(
    "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_INTERLEUKIN_21_PRODUCTION",
    "GOBP_INTERLEUKIN_4_PRODUCTION",
    "GOBP_TYPE_II_INTERFERON_PRODUCTION",
    "GOBP_LYMPHOCYTE_CHEMOTAXIS",
    "GOBP_TRANSLATIONAL_ELONGATION",
    "GOBP_TRANSLATIONAL_INITIATION"
  ))

# Combine the gene sets
b_gene_sets <- dplyr::bind_rows(b_cell_gene_sets)
t_gene_sets <- dplyr::bind_rows( t_cell_gene_sets)
# Create a list of gene sets
b_gene_sets <- split(b_gene_sets$human_gene_symbol, b_gene_sets$gs_name)
t_gene_sets <- split(t_gene_sets$human_gene_symbol, t_gene_sets$gs_name)
# try downloaded version comment above code
b_cell_gene_sets <- c( "GOBP_RESPONSE_TO_CYTOKINE",
                       "GSE23925_DARK_ZONE_VS_NAIVE_BCELL_UP",
                       "CELL_PROLIFERATION_GO_0008283",
                       "GOBP_REGULATION_OF_B_CELL_PROLIFERATION",
                       "GOBP_B_CELL_PROLIFERATION",
                       "GOBP_B_CELL_PROLIFERATION_INVOLVED_IN_IMMUNE_RESPONSE",
                       "GOBP_GERMINAL_CENTER_B_CELL_DIFFERENTIATION",
                       "GOBP_MATURE_B_CELL_DIFFERENTIATION",
                       "GOBP_B_CELL_DIFFERENTIATION",
                       "GOBP_B_CELL_ACTIVATION",
                       "GOBP_INTERLEUKIN_6_MEDIATED_SIGNALING_PATHWAY",
                       "GOBP_MHC_CLASS_II_BIOSYNTHETIC_PROCESS",
                       "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
                       "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX",
                       "GOBP_DNA_REPLICATION",
                       "GOBP_DNA_TOPOLOGICAL_CHANGE",
                       "GOBP_DNA_LIGATION")
t_cell_gene_sets <- c("GOBP_CYTOKINE_PRODUCTION",
                      "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                      "GOBP_INTERLEUKIN_21_PRODUCTION",
                      "GOBP_INTERLEUKIN_4_PRODUCTION",
                      "GOBP_TYPE_II_INTERFERON_PRODUCTION",
                      "GOBP_LYMPHOCYTE_CHEMOTAXIS",
                      "GOBP_TRANSLATIONAL_ELONGATION",
                      "GOBP_TRANSLATIONAL_INITIATION")
b_gene_sets <- list()
for (i in 1:length(b_cell_gene_sets)){
  tmp.file.read <- read.table(file.path(wdpath,"SGCC_relevant_analysis","GSEA_signature",
                                        paste0(b_cell_gene_sets[i],".v2023.2.Hs.grp")))$V1[-1]
  b_gene_sets[[i]] <- tmp.file.read
  names(b_gene_sets)[i] <- b_cell_gene_sets[i]
}
t_gene_sets <- list()
for (i in 1:length(t_cell_gene_sets)){
  tmp.file.read <- read.table(file.path(wdpath,"SGCC_relevant_analysis","GSEA_signature",
                                        paste0(t_cell_gene_sets[i],".v2023.2.Hs.grp")))$V1[-1]
  t_gene_sets[[i]] <- tmp.file.read
  names(t_gene_sets)[i] <- t_cell_gene_sets[i]
}
# Perform GSVA
L24_B_gsva <- gsvaParam(exprData = as.matrix(L24_norm_mat_B), geneSets = b_gene_sets,
                      minSize=1, maxSize=1000)
L24_T_gsva <- gsvaParam(exprData = as.matrix(L24_norm_mat_T), geneSets = t_gene_sets,
                        minSize=1, maxSize=1000)
L25_B_gsva <- gsvaParam(exprData = as.matrix(L25_norm_mat_B), geneSets = b_gene_sets,
                        minSize=1, maxSize=1000)
L25_T_gsva <- gsvaParam(exprData = as.matrix(L25_norm_mat_T), geneSets = t_gene_sets,
                        minSize=1, maxSize=1000)

L24_B_gsva_results <- gsva(L24_B_gsva)
L24_T_gsva_results <- gsva(L24_T_gsva)
L25_B_gsva_results <- gsva(L25_B_gsva)
L25_T_gsva_results <- gsva(L25_T_gsva)

L24_B_T_gsva <- rbind(L24_B_gsva_results,L24_T_gsva_results)
L25_B_T_gsva <- rbind(L25_B_gsva_results,L25_T_gsva_results)

merged_gsva <- cbind(L24_B_T_gsva,L25_B_T_gsva)
merged_gsva <- merged_gsva[,SGCC_reorder.name]

scale_mat <- t(scale(t(as.matrix(merged_gsva))))

rowname.order <- c( "GOBP_RESPONSE_TO_CYTOKINE",
                   "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX",
                   "GOBP_DNA_TOPOLOGICAL_CHANGE",
                   "GOBP_INTERLEUKIN_6_MEDIATED_SIGNALING_PATHWAY",
                   "GOBP_INTERLEUKIN_21_PRODUCTION",
                   "GOBP_INTERLEUKIN_4_PRODUCTION",
                   "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                   "GOBP_TRANSLATIONAL_INITIATION",
                   "GOBP_CYTOKINE_PRODUCTION")
ht_list <- Heatmap(scale_mat[rowname.order,],
                   cluster_columns = F,
                   cluster_rows = F, 
                   row_split = c(rep("A",4),
                                 rep("B",5)),
                   row_names_gp = gpar(fontsize = 10),top_annotation = column_density_plot,
                   col = colorRamp2(c(min(t(scale(t(as.matrix(scale_mat))))+0.5), 
                                      0, max(t(scale(t(as.matrix(scale_mat)))))-0.5), 
                                    c("#b9529f", "#0e0a0c", "#f3ec19")))
svg("Figure E pathway.svg", width = 10,height = 30)
draw(ht_list)
dev.off()
# plot annotation bar
extract_protein_cell <- function(CODEX_Data){
  # 1. Summarize total cell number for each core
  total_cell_count <- CODEX_Data %>%
    group_by(ROI_rename) %>%
    summarise(Total_Cells = n(), .groups = "drop")
  
  # 2. Summarize the number of CD4T, BCL6nB, and BCL6pB cells in each core and their proportions
  specific_cell_counts <- CODEX_Data %>%
    filter(SegmentLabel %in% c("CD4T", "BCL6nB", "BCL6pB")) %>%
    group_by(ROI_rename, SegmentLabel) %>%
    summarise(CellCount = n(), .groups = "drop") %>%
    pivot_wider(names_from = SegmentLabel, values_from = CellCount, values_fill = 0) %>%
    left_join(total_cell_count, by = "ROI_rename") %>%
    mutate(
      CD4T_Proportion = CD4T / Total_Cells,
      BCL6nB_Proportion = BCL6nB / Total_Cells,
      BCL6pB_Proportion = BCL6pB / Total_Cells
    ) %>% 
    dplyr::select(ROI_rename,BCL6nB, BCL6pB,  CD4T, CD4T_Proportion, BCL6nB_Proportion, BCL6pB_Proportion)
  
  # 3. Summarize M1 and M2 cell numbers, their proportions, and the M1/M2 ratio in each core
  m1_m2_summary <- CODEX_Data %>%
    filter(SegmentLabel %in% c("M1", "M2")) %>%
    group_by(ROI_rename, SegmentLabel) %>%
    summarise(CellCount = n(), .groups = "drop") %>%
    pivot_wider(names_from = SegmentLabel, values_from = CellCount, values_fill = 0) %>%
    mutate(
      M1_Proportion = M1 / (M1+M2),
      M2_Proportion = M2 / (M1+M2),
      M1_M2_Ratio = ifelse(M2 > 0, M1 / M2, NA)
    )
  
  # 4. Summarize mean expression of PAX5 and BCL6 for BCL6nB in each core
  bcl6nb_expression_summary <- CODEX_Data %>%
    filter(SegmentLabel == "BCL6nB") %>%
    group_by(ROI_rename) %>%
    summarise(
      Mean_PAX5 = mean(Pax5, na.rm = TRUE),
      Mean_BCL6 = mean(BCL6, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 5. Summarize mean expression of CD163 protein for M1 and M2 in each core
  cd163_expression_summary <- CODEX_Data %>%
    filter(SegmentLabel %in% c("M1", "M2")) %>%
    group_by(ROI_rename, SegmentLabel) %>%
    summarise(Mean_CD163 = mean(CD163, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = SegmentLabel, values_from = Mean_CD163, names_prefix = "Mean_CD163_")
  
  
  # Combine summaries as needed
  data_summary <- total_cell_count %>%
    left_join(specific_cell_counts, by = "ROI_rename") %>%
    left_join(m1_m2_summary, by = "ROI_rename") %>%
    left_join(bcl6nb_expression_summary, by = "ROI_rename") %>%
    left_join(cd163_expression_summary, by = "ROI_rename") %>%
    dplyr::select(
      ROI_rename,
      Total_Cells,
      CD4T, CD4T_Proportion,
      BCL6nB, BCL6nB_Proportion,
      BCL6pB, BCL6pB_Proportion,
      M1, M1_Proportion,
      M2, M2_Proportion,
      M1_M2_Ratio,
      Mean_PAX5, Mean_BCL6,
      Mean_CD163_M1, Mean_CD163_M2
    )
  return(data_summary)
}


L24_CODEX$ROI_rename <- paste0("L24_",L24_CODEX$ROI_num)
L25_CODEX$ROI_rename <- paste0("L25_",L25_CODEX$ROI_num)

L24_codex_summary <- extract_protein_cell(L24_CODEX)
L25_codex_summary <- extract_protein_cell(L25_CODEX)
L24_L25_codex_summary <- rbind(L24_codex_summary, L25_codex_summary)
L24_L25_codex_summary <- L24_L25_codex_summary[match(SGCC_reorder.name,L24_L25_codex_summary$ROI_rename),]

# Define the column annotation
# Define color functions for continuous variables
color_fun_pax5 <- colorRamp2(c(min(L24_L25_codex_summary$Mean_PAX5, na.rm = TRUE), 
                               mean(L24_L25_codex_summary$Mean_PAX5, na.rm = TRUE), 
                               max(L24_L25_codex_summary$Mean_PAX5, na.rm = TRUE)),
                             c("blue", "white", "red"))

color_fun_bcl6 <- colorRamp2(c(min(L24_L25_codex_summary$Mean_BCL6, na.rm = TRUE), 
                               mean(L24_L25_codex_summary$Mean_BCL6, na.rm = TRUE), 
                               max(L24_L25_codex_summary$Mean_BCL6, na.rm = TRUE)),
                             c("blue", "white", "red"))

color_fun_cd163_m1 <- colorRamp2(c(min(L24_L25_codex_summary$Mean_CD163_M1, na.rm = TRUE), 
                                   mean(L24_L25_codex_summary$Mean_CD163_M1, na.rm = TRUE), 
                                   max(L24_L25_codex_summary$Mean_CD163_M1, na.rm = TRUE)),
                                 c("blue", "white", "red"))

color_fun_cd163_m2 <- colorRamp2(c(min(L24_L25_codex_summary$Mean_CD163_M2, na.rm = TRUE), 
                                   mean(L24_L25_codex_summary$Mean_CD163_M2, na.rm = TRUE), 
                                   max(L24_L25_codex_summary$Mean_CD163_M2, na.rm = TRUE)),
                                 c("blue", "white", "red"))

# Define column annotation
col_annotation <- HeatmapAnnotation(
  Total_Cells = anno_barplot(L24_L25_codex_summary$Total_Cells, gp = gpar(fill = "black")),
 # CD4T = anno_barplot(L24_L25_codex_summary$CD4T, gp = gpar(fill = "#2958a7")),
  CD4T_Proportion = anno_barplot(L24_L25_codex_summary$CD4T_Proportion, gp = gpar(fill = "#2958a7")),
 # BCL6nB = anno_barplot(L24_L25_codex_summary$BCL6nB, gp = gpar(fill = "#e07229")),
  BCL6nB_Proportion = anno_barplot(L24_L25_codex_summary$BCL6nB_Proportion, gp = gpar(fill = "#e07229")),
#  BCL6pB = anno_barplot(L24_L25_codex_summary$BCL6pB, gp = gpar(fill = "#12964b")),
  BCL6pB_Proportion = anno_barplot(L24_L25_codex_summary$BCL6pB_Proportion, gp = gpar(fill = "#12964b")),
#  M1 = anno_barplot(L24_L25_codex_summary$M1, gp = gpar(fill = "#0796b1")),
  M1_Proportion = anno_barplot(L24_L25_codex_summary$M1_Proportion, gp = gpar(fill = "#0796b1")),
#  M2 = anno_barplot(L24_L25_codex_summary$M2, gp = gpar(fill = "#b8e4f2")),
  M2_Proportion = anno_barplot(L24_L25_codex_summary$M2_Proportion, gp = gpar(fill = "#b8e4f2")),
  M1_M2_Ratio = anno_barplot(L24_L25_codex_summary$M1_M2_Ratio, gp = gpar(fill = "red")),
  M1_M2_stackBar = anno_barplot(cbind(L24_L25_codex_summary$M1_Proportion,
                                      L24_L25_codex_summary$M2_Proportion),
                                gp = gpar(fill = c("#0796b1", "#b8e4f2"), 
                                          col = c("#0796b1", "#b8e4f2"))),
  Mean_PAX5 = anno_simple(L24_L25_codex_summary$Mean_PAX5, col = color_fun_pax5),
  Mean_BCL6 = anno_simple(L24_L25_codex_summary$Mean_BCL6, col = color_fun_bcl6),
  Mean_CD163_M1 = anno_simple(L24_L25_codex_summary$Mean_CD163_M1, col = color_fun_cd163_m1),
  Mean_CD163_M2 = anno_simple(L24_L25_codex_summary$Mean_CD163_M2, col = color_fun_cd163_m2),
  annotation_name_gp = gpar(fontsize = 10)
)

# Display the column annotation

pdf("Figure 3E and Figure 3G annotation bar.pdf", width = 10,height = 9)
draw(Heatmap(matrix(nrow = 0, ncol = nrow(L24_L25_codex_summary)), 
             top_annotation = col_annotation, 
             show_heatmap_legend = FALSE, 
             show_row_names = FALSE, 
             show_column_names = FALSE))
dev.off()
