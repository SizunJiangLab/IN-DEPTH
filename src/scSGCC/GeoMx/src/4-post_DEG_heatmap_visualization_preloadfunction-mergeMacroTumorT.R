library(dplyr)
library(ComplexHeatmap)
library(qs)
library(Seurat)
library(circlize)
library(tidyr)
library(tidyverse)
# get structral gene list
library(msigdbr)
#
# Retrieve GO Cellular Component gene sets for Homo sapiens (human)
go_cc_genes <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")

go_cc_genes_list <- unique(go_cc_genes$gene_symbol[go_cc_genes$gs_name %in% unique(grep("house",go_cc_genes$gs_name,ignore.case = T,value = T))])
# prepare the EBV genes
# Creating an R vector with EBV gene names
EBV_genes <- c("EBNA1", "EBNA1.1", "EBNA1.2",
               "EBNA2", "EBNA2.1", "EBNA2.2",
               "EBNA1BP2",
               "EBNALP", "EBNALP.1", "EBNALP.2",
               "LMP1", "LMP1.1", "LMP1.2",
               "BZLF1", "BZLF1.1", "BZLF1.2",
               "BHRF1", "BHRF1.1", "BHRF1.2",
               "RPMS1", "RPMS1.1", "RPMS1.2",
               "BALF1", "BALF1.1", "BALF1.2",
               "BCRF1", "BCRF1.1", "BCRF1.2",
               "BNLF2A", "BNLF2A.1", "BNLF2A.2",
               "BNRF1", "BNRF1.1", "BNRF1.2",
               "EBER1", "EBER1.1", "EBER2", "EBER2.1")

#mt.seurat_log_mergermacro1 <- qread(file.path(wdpath,"GeoMX_batchcorrection","LogCPM-SPE_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated.qs"))
mt.seurat_log_mergermacro <- qread(file.path(wdpath,"GeoMX_batchcorrection","LogCPM-Seurat_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated-update.qs"))
# CD4T_expr_mat <- mt.seurat_log_mergermacro1$`LogCPM_RUV4_top1000_k1_bn_Cohort_wf_MergedLabel+EBV_Indicator`@assays@data$logcounts
CD4T_expr_mat <- mt.seurat_log_mergermacro$Seurat_top1000_k3@assays$RNA@layers$counts
rownames(CD4T_expr_mat) <- rownames(mt.seurat_log_mergermacro$Seurat_top1000_k3)
colnames(CD4T_expr_mat) <- colnames(mt.seurat_log_mergermacro$Seurat_top1000_k3)
Macro_expr_mat <- mt.seurat_log_mergermacro$Seurat_top1000_k3@assays$RNA@layers$counts
rownames(Macro_expr_mat) <- rownames(mt.seurat_log_mergermacro$Seurat_top1000_k3)
colnames(Macro_expr_mat) <- colnames(mt.seurat_log_mergermacro$Seurat_top1000_k3)
Tumor_expr_mat <- mt.seurat_log_mergermacro$Seurat_top1000_k3@assays$RNA@layers$counts
rownames(Tumor_expr_mat) <- rownames(mt.seurat_log_mergermacro$Seurat_top1000_k3)
colnames(Tumor_expr_mat) <- colnames(mt.seurat_log_mergermacro$Seurat_top1000_k3)
#
Seurat_virus <- mt.seurat_log_mergermacro$Seurat_top1000_k3
rownames(Seurat_virus@assays$RNA@layers$counts) <- rownames(Seurat_virus)
EBV_expr <- Seurat_virus@assays$RNA@layers$counts[EBV_genes, ]
total_viral_expr <- colSums(EBV_expr)
total_gene_expr <- colSums(Seurat_virus@assays$RNA@layers$counts)
virus_loading <- total_viral_expr / total_gene_expr
Seurat_virus <- AddMetaData(Seurat_virus, metadata = data.frame(row.names = colnames(Seurat_virus), virus_loading = virus_loading))
Seurat_virus$CXCL9 <- Seurat_virus@assays$RNA@layers$counts["CXCL9",]
Seurat_virus$SPP1 <- Seurat_virus@assays$RNA@layers$counts["SPP1",]
Seurat_virus$CXCL9_SPP1 <- Seurat_virus$CXCL9/Seurat_virus$SPP1
Seurat_virus$LMP1_g <-Seurat_virus@assays$RNA@layers$counts["LMP1",]
# subset M1 and M2
CXCL9_SPP1_score <- Seurat_virus@meta.data %>%
  filter(MergedLabel == "Macro" ) %>%
  group_by(ROI_rename) %>%
  summarise(CXCL9_SPP1_mean = mean(CXCL9_SPP1) )
Overall_LMP1Expression <- Seurat_virus@meta.data %>%
  filter(MergedLabel == "Tumor" ) %>%
  group_by(ROI_rename) %>%
  summarise(LMP1_g_tumor = mean(LMP1_g))
# Seurat_virus@meta.data[1:5,]
# calculate the tumor virus loading
tumor_meta <- subset(Seurat_virus@meta.data, grepl("Tumor", MergedLabel))
tumor_meta <- tumor_meta[,c("ROI_rename","virus_loading")]
colnames(tumor_meta)[2] <- "virus_loading_InTumor"
tumor_meta$virus_loading_InTumor <- scale(tumor_meta$virus_loading_InTumor)

reformedMeta <- Seurat_virus@meta.data %>%
  left_join(y = tumor_meta,by = "ROI_rename") %>%
  left_join(y=CXCL9_SPP1_score, by = "ROI_rename") %>%
  left_join(y = Overall_LMP1Expression, by = "ROI_rename")

# calculate tumor proportion
annotation <- read_csv(file.path(wdpath,"SGCC_relevant_analysis","DLBCL_ROIlevel_annotation.csv"))
annotation_old <- read_csv(file.path(wdpath,"SGCC_relevant_analysis","DLBCL_ROIlevel_markers.csv"))
annotation$LMP1<- annotation_old$LMP1
table(annotation$Annotation)
# Define tumor-related annotations
tumor_annotations <- c("Other Tumor", "Tumor BCL2", "Tumor BCL6", "Tumor Myc")

# Calculate total tumor cells and tumor proportion for each coreName
tumor_summary <- annotation %>%
  group_by(coreName) %>%
  summarise(
    total_cells = n(),  # Total number of cells in each coreName
    total_tumor_cells = sum(Annotation %in% tumor_annotations),  # Total tumor cells in each coreName
    Tumor_other_cell = sum(Annotation == "Other Tumor"),  # Number of "Tumor" cells
    Tumor_BCL2_cell = sum(Annotation == "Tumor BCL2"),  # Number of "Tumor BCL2" cells
    Tumor_BCL6_cell = sum(Annotation == "Tumor BCL6"),  # Number of "Tumor BCL6" cells
    Tumor_Myc_cell = sum(Annotation == "Tumor Myc")  # Number of "Tumor Myc" cells
  )
tumor_summary <- tumor_summary %>%
  dplyr::rename(ROI_rename = coreName)
# calculate M1 ratio
M1_ratio <- annotation %>%
  group_by(coreName) %>%
  summarise(
    m1_cells = sum(Annotation == "M1"),
    m2_cells = sum(Annotation == "M2"),
    macrocell = m1_cells+ m2_cells,
    m1_pct = m1_cells/macrocell,
    m2_pct = m2_cells/macrocell,
    m1_odds =  m1_cells/m2_cells,
    m2_odds =  m2_cells/m1_cells,
  )
M1_ratio <- M1_ratio %>%
  dplyr::rename(ROI_rename = coreName)
###
## read in LMP1 table
LMP1_table <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","DLBCL_Tumor_LMP1.csv"))
LMP1_table$CellID <- paste0(LMP1_table$cellLabel,"_", LMP1_table$coreName)

# read in tumor cell inforamtion
tumor_LMP1 <- annotation %>%
  mutate(CellID = paste0(cellLabel, "_", coreName)) %>%
  filter(Annotation %in% tumor_annotations) %>%
  left_join(LMP1_table %>% select(CellID, LMP1_thresholded, LMP1_Positive), by = "CellID") %>%
  group_by(coreName) %>%
  summarise(
    LMP1_filtered_mean = mean(LMP1, na.rm = TRUE),
    LMP1_positive_numebr = sum(LMP1_Positive, na.rm = TRUE)
  )
tumor_LMP1$LMP1_filtered_mean[tumor_LMP1$coreName %in% "DFCI_14.1"] <- NA

tumor_LMP1_summary <- tumor_LMP1 %>%
  dplyr::rename(ROI_rename = coreName)
# other protein information
annotation_merged <- annotation %>%
  mutate(Annotation = case_when(
    grepl("CD8", Annotation) ~ "CD8T",
    grepl("CD4", Annotation) ~ "CD4T",
    grepl("Tumor", Annotation) ~ "Tumor",
    grepl("M1", Annotation) ~ "Macro",
    grepl("M2", Annotation) ~ "Macro",
    TRUE ~ Annotation
  ))
#
CD4T_function <- annotation_merged %>%
  filter(Annotation %in% c("CD4T")) %>%
  group_by(coreName) %>%
  summarise(
    PD1_CD4T = mean(`PD-1`),
    Tox_CD4T = mean(Tox),
    LAG3_CD4T = mean(LAG3),
    CD45RA_CD4T = mean(CD45RA),
    CD45RO_CD4T = mean(CD45RO),
    Ki67_CD4T = mean(Ki67),
    GZMB_CD4T = mean(GZMB)
  ) %>%
  mutate(dysfunction_score_protein = LAG3_CD4T +
           CD45RO_CD4T - CD45RA_CD4T + Tox_CD4T - Ki67_CD4T - GZMB_CD4T)
CD4T_function <- CD4T_function %>%
  dplyr::rename(ROI_rename = coreName)
Macro_function <- annotation_merged %>%
  filter(Annotation %in% c("Macro")) %>%
  group_by(coreName) %>%
  summarise(
    HLADR_macro = mean(`HLA-DR`),
    PDL1_macro = mean(`PD-L1`)
  )
Macro_function <- Macro_function %>%
  dplyr::rename(ROI_rename = coreName)


reformedMeta <- reformedMeta %>%
  left_join(y = tumor_summary, by = "ROI_rename") %>%
  left_join(y = tumor_LMP1_summary, by = "ROI_rename") %>%
  left_join(y = M1_ratio,by = "ROI_rename" ) %>%
  left_join( y = CD4T_function,by = "ROI_rename") %>%
  left_join( y = Macro_function,by = "ROI_rename") %>%
  mutate(LMP1_pct = LMP1_positive_numebr/total_tumor_cells)

reformedMeta$LMP1_pct = reformedMeta$LMP1_positive_numebr/ reformedMeta$total_tumor_cells

reformedMeta$ROI_CT <- paste0(reformedMeta$ROI_rename,"_",reformedMeta$MergedLabel)

list_exp_matrix_metadata <- list(CD4T = CD4T_expr_mat,
                                 #CD8T = CD8T_expr_mat,
                                 #M1 = M1_expr_mat,
                                 #M2 = M2_expr_mat,
                                 Macro = Macro_expr_mat,
                                 Tumor = Tumor_expr_mat,
                                 reformedMeta = reformedMeta,
                                 annotation = annotation)


create_complex_heatmap_FixEBV <- function(gene_list,pathway_matrix, expr_data, meta_data,scale_color_ht = 3,cluster_rows = T, cluster_columns = T) {
  meta_data <- meta_data %>%
    arrange(EBV_status, SGCC)
  # Define color functions
  scale_virus_loading_InTumor <- scale(meta_data$virus_loading_InTumor)
  scale_LMP1_g_tumor <- scale(meta_data$LMP1_g_tumor)
  scale_LMP1_filtered_mean <-scale(meta_data$LMP1_filtered_mean)
  scale_LMP1_positive_numebr <- scale(meta_data$LMP1_positive_numebr)
  scale_PD1_CD4T <- scale(meta_data$PD1_CD4T)
  scale_Tox_CD4T <- scale(meta_data$Tox_CD4T)
  scale_LAG3_CD4T <- scale(meta_data$LAG3_CD4T)
  scale_HLADR_macro <- scale(meta_data$HLADR_macro)
  scale_PDL1_macro <- scale(meta_data$PDL1_macro)
  scale_fun_CD4TDysfunction <- scale(meta_data$dysfunction_score_protein)

  global_max <- max(scale_LMP1_g_tumor,scale_LMP1_filtered_mean,scale_LMP1_positive_numebr,
                    scale_PD1_CD4T, scale_Tox_CD4T, scale_LAG3_CD4T,scale_HLADR_macro,scale_PDL1_macro,
                    scale_fun_CD4TDysfunction)
  col_global <- colorRamp2(c(-2,0, 2), c("blue","white", "red"))
  col_fun_LMP1_pct <- colorRamp2(c(0,0.5,1), c("blue","white", "red"))
  col_fun_m1_ratio <- colorRamp2(c(0,1), c("#f6f7f5", "#df536b"))
  col_fun_m2_ratio <- colorRamp2(c(0,1), c("#f6f7f5", "#62cf4f"))
  # 2. Create top annotation (stacked barplot) for tumor proportions
  tumor_data <- meta_data %>%
    select(Tumor_other_cell, Tumor_BCL2_cell, Tumor_BCL6_cell, Tumor_Myc_cell) %>%
    mutate(across(everything(), ~ . / meta_data$total_tumor_cells)) %>%
    as.matrix()
  rownames(tumor_data) <- meta_data$Sample

  meta_data[,colnames(tumor_data)] <- tumor_data
  meta_data$tumorProportion <- meta_data$total_tumor_cells/meta_data$total_cells
  tumor_proportion_anno <- HeatmapAnnotation(
    EBV_Status = meta_data$EBV_status,
    SGCC_Category = meta_data$SGCC_cat,
    SGCC_vaue = anno_lines(meta_data$SGCC),
    Tumor_Proportion = anno_barplot(meta_data[,colnames(tumor_data)],
                                    border = TRUE,
                                    gp = gpar(fill = c("Tumor_other_cell" = "#66c2a5",
                                                       "Tumor_BCL2_cell" = "#fc8d62",
                                                       "Tumor_BCL6_cell"="#8da0cb",
                                                       "Tumor_Myc_cell" = "#e78ac3"))),
    TumorAbundance_ROI = anno_barplot(meta_data$tumorProportion, gp = gpar(fill = "red")),
    EBV_Loading = anno_simple(scale_virus_loading_InTumor, col = col_global),
    LMP1_pos_tumor_number = anno_simple(scale_LMP1_positive_numebr, col = col_global),
    LMP1_pct_over_tumor = anno_barplot(meta_data$LMP1_pct, gp = gpar(fill = "#32a8cc")),
    LMP1_gene = anno_simple(scale_LMP1_g_tumor,col = col_global),
    LMP1_protein = anno_simple(scale_LMP1_filtered_mean,col = col_global),
    m1_m2 = anno_barplot(cbind(meta_data$m1_pct, meta_data$m2_pct), gp = gpar(fill = 2:3, col = 2:3)),
    m1_pct = anno_simple(meta_data$m1_pct,col = col_fun_m1_ratio),
    m2_pct = anno_simple(meta_data$m2_pct,col = col_fun_m2_ratio),
    PD1_CD4T = anno_simple(scale_PD1_CD4T,col = col_global),
    Tox_CD4T = anno_simple(scale_Tox_CD4T,col = col_global),
    LAG3_CD4T = anno_simple(scale_LAG3_CD4T,col = col_global),
    HLADR_macro = anno_simple(scale_HLADR_macro,col = col_global),
    PDL1_macro = anno_simple(scale_PDL1_macro,col = col_global),
    dysfunction_score_protein = anno_simple(scale_fun_CD4TDysfunction,col = col_global),
    col = list( LMP1_g = col_global,
                SGCC_Category = c("SGCC low" = "#ffcad4", "SGCC mediate" = "#f3abb6", "SGCC high" = "#9f8189"),
                EBV_Status = c("EBV+" = "#d95f02", "EBV-" = "#1b9e77"),
                Tumor_Proportion = c("Tumor_other_cell" = "#66c2a5",
                                    "Tumor_BCL2_cell" = "#fc8d62",
                                    "Tumor_BCL6_cell"="#8da0cb",
                                    "Tumor_Myc_cell" = "#e78ac3")),
    show_annotation_name = TRUE,show_legend = T,
    which = "column",
    height = unit(8, "cm"),
    gap = unit(1, "mm")
  )

  lgd_list <- list(
    Legend(labels = c("Tumor_other_cell", "Tumor_BCL2_cell", "Tumor_BCL6_cell", "Tumor_Myc_cell"),
           title = "Tumor Proportions",
           type = "points",
           pch = 20,
           legend_gp = gpar(col =  c("#66c2a5",
                                     "#fc8d62",
                                     "#8da0cb",
                                     "#e78ac3"))))
  # Create heatmap object
  if(length(gene_list) < 1|all(is.na(gene_list))){
    pathway_matrix_subset <- pathway_matrix[,meta_data$Sample]
    heatmap_mat <- t(scale(t(pathway_matrix_subset)))
  }else if(is.na(pathway_matrix)){
    expr_data_subset <- expr_data[gene_list,]
    expr_data_subset <- expr_data_subset[, meta_data$Sample]
    heatmap_mat <- t(scale(t(expr_data_subset)))
  }else{
    expr_data_subset <- expr_data[gene_list,]
    expr_data_subset <- expr_data_subset[, meta_data$Sample]
    scaled_mat <- t(scale(t(expr_data_subset)))
    pathway_matrix_subset <- pathway_matrix[,meta_data$Sample]
    scaled_pathway_matrix_subset <- t(scale(t(pathway_matrix_subset)))
    heatmap_mat <- rbind(scaled_mat, scaled_pathway_matrix_subset)
  }
  uni_value <- max(abs(scale_color_ht))
  heatmap <- Heatmap(
    heatmap_mat,
    name = "Expression",
    row_title = "features",
    column_title = "Samples",
    top_annotation = tumor_proportion_anno,
    show_column_names = T,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    col = colorRamp2(c(-uni_value,0, uni_value), c("#b953a0", "black", "#f4ed17"))
  )
  # Draw heatmap
  draw(heatmap, heatmap_legend_side = "right",annotation_legend_list = lgd_list)
}


create_complex_heatmap_annotation_from_pathway <- function(pathway_matrix,scale_color_ht = 2, meta_data, pathwayname = "CD4T_customized_Tcell_exhaustion") {
  meta_data <- meta_data %>%
    arrange(EBV_status, SGCC)
  pathway_score <- scale(pathway_matrix[pathwayname,meta_data$Sample])[,1]
  col_global <- colorRamp2(c(-2,0, 2), c("blue","white", "red"))
  tumor_proportion_anno <- HeatmapAnnotation(
  EBV_Status = meta_data$EBV_status,
  SGCC_Category = meta_data$SGCC_cat,
  dysfunction_gene_CD4T =  anno_simple(pathway_score,col = col_global),
  col = list( LMP1_g = col_global,
              SGCC_Category = c("SGCC low" = "#ffcad4", "SGCC mediate" = "#f3abb6", "SGCC high" = "#9f8189"),
              EBV_Status = c("EBV+" = "#d95f02", "EBV-" = "#1b9e77"),
              Tumor_Proportion = c("Tumor_other_cell" = "#66c2a5",
                                   "Tumor_BCL2_cell" = "#fc8d62",
                                   "Tumor_BCL6_cell"="#8da0cb",
                                   "Tumor_Myc_cell" = "#e78ac3")),
  show_annotation_name = TRUE,show_legend = T,
  which = "column",
  height = unit(8, "cm"),
  gap = unit(1, "mm")
  #  CXCL9_SPP1 = anno_simple(meta_data$CXCL9_SPP1_mean,col = col_fun_CXCL9),
  )
  pathway_matrix_ht <- t(scale(t(pathway_matrix[,meta_data$Sample])))
  heatmap <- Heatmap(
    pathway_matrix_ht,
    name = "Expression",
    row_title = "features",
    column_title = "Samples",
    top_annotation = tumor_proportion_anno,
    show_column_names = T,
    cluster_rows = T,
    cluster_columns = F,
    col = colorRamp2(c(-scale_color_ht,0, scale_color_ht), c("#b953a0", "black", "#f4ed17"))
  )
  draw(heatmap)
}


