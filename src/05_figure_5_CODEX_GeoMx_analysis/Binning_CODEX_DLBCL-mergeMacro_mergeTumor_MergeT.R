## set working directory
loadingfunction_wdpath <- c("./src/")
wdpath <- "./data/DLBCL_run/"
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(standR)
library(SingleCellExperiment)
library(SpatialExperiment)
library(GSVA)
library(readxl)
library(utils)
library(reshape2)
library(Matrix.utils)
library(ComplexHeatmap)
library(circlize)
library(edgeR)
source(file = file.path(loadingfunction_wdpath,"SGCC_code","1-SpaGFT_preload_function.R"))
#remotes::install_github("YosefLab/ImpulseDE2")
# Set dimensions
x_border <- 60
# Create a base scatter plot
base_data <- expand.grid(x = 1:x_border, y = 1:x_border)
output.res <- Cal_Eigen(data.in = base_data,k = 400,k_fold = 10)
My_CODEXposition <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","DLBCL_ROIlevel_annotation.csv"))
dim(My_CODEXposition)
dim(annotation)
My_CODEXposition[1:5,]
table(My_CODEXposition$Annotation)
names(table(My_CODEXposition$Annotation))
table(My_CODEXposition$coreName)
## remove samples
My_CODEXposition <- My_CODEXposition[My_CODEXposition$coreName != "Rochester_TonsilA" & My_CODEXposition$coreName != "DFCI_Tonsil1",]
# change names
mapping_table <- data.frame(
  SegmentLabel = c("CD4mem", "CD4naive", "CD8mem", "CD8naive", "DC", "Endothelial",
                   "M1", "M2", "Neutrophil", "Other", "Treg", "Tumor", "TumorBCL2",
                   "TumorBCL6", "TumorMyc", "TumorOther"),
  Annotation = c("CD4mem", "CD4naive", "CD8mem", "CD8naive", "DC", "Endothelial",
                 "M1", "M2", "Neutrophil", "Other", "Treg", "Other Tumor", "Tumor BCL2",
                 "Tumor BCL6", "Tumor Myc", "Other Tumor"),
  MergedTumor = c("CD4T", "CD4T", "CD8T", "CD8T", "DC", "Endothelial",
                  "Macro", "Macro", "Neutrophil", "Other", "Treg", "Tumor", "Tumor",
                  "Tumor", "Tumor", "Tumor")
)
My_CODEXposition_merge <- left_join(My_CODEXposition,y = mapping_table, by = "Annotation")
table(My_CODEXposition_merge$MergedTumor,My_CODEXposition_merge$Annotation)
My_CODEXposition_list <- split(My_CODEXposition_merge, ~coreName)
My_CODEXposition_list[[1]]
sapply(My_CODEXposition_list, function(x) table(x$Annotation))
sapply(My_CODEXposition_list, dim)
# plot
# ggplot(My_CODEXAnnotation_list$DFCI_1.2, aes(x= X_cent, y = Y_cent,color = Annotation))+
#   geom_point()+
#   scale_y_reverse()


cellpairs <- combn(c("CD4T","CD8T","DC", "Endothelial",
                     "Macro", "Neutrophil", "Other", "Treg", "Tumor"), 2)
color_mapping <- data.frame(
  MergedTumor = c("CD4T","CD8T" , "DC","Endothelial","Macro",
                  "Neutrophil","Treg","Other","Tumor"),
  Color = c("#5cb259","#FF7600", "#ab3fa0", "#5A5AB7","#8c549c",
            "#2182C2","#AD7442","#E92C7C","#da3a38"))
cellcombination.name <- apply(cellpairs, 2, function(x) paste(x, collapse = "_"))
# LSCI_matrix <- matrix(NA, nrow = dim(cellpairs)[2], ncol = length(My_CODEXposition_list))
# rownames(LSCI_matrix) <- cellcombination.name
# colnames(LSCI_matrix) <- names(My_CODEXposition_list)
SGCC_matrix <- matrix(NA, nrow = dim(cellpairs)[2], ncol = length(My_CODEXposition_list))
rownames(SGCC_matrix) <- cellcombination.name
colnames(SGCC_matrix) <- names(My_CODEXposition_list)
# run scores
for (i in names(My_CODEXposition_list)){
  ROI_number <- i
  data <- My_CODEXposition_list[[ROI_number]]
  # Create 3600 bins (60x60 grid)
  x_bins <- x_border
  y_bins <- x_border
  data <- data %>%
    mutate(
      x_bin = cut(X_cent_fusion, breaks = x_bins, labels = FALSE),
      y_bin = cut(Y_cent_fusion, breaks = y_bins, labels = FALSE)
    )

  # Count the number of cells of each type in each bin
  bin_counts <- data %>%
    group_by(x_bin, y_bin, MergedTumor) %>%
    summarise(cell_count = n(), .groups = 'drop')

  # Calculate the total number of cells of each type
  total_cells_per_type <- data %>%
    group_by(MergedTumor) %>%
    summarise(total_count = n(), .groups = 'drop')

  # Merge total counts with bin counts
  bin_counts <- bin_counts %>%
    left_join(total_cells_per_type, by = "MergedTumor") %>%
    mutate(cell_proportion = cell_count / total_count)
  ###
  for (k in 1:length(cellcombination.name)){
    # Filter for a specific cell type, for example, "CD4T"
    tmp_cellcombination.name <- cellcombination.name[k]
    CT1 <- strsplit(tmp_cellcombination.name,"_")[[1]][1]

    CT2 <- strsplit(tmp_cellcombination.name,"_")[[1]][2]
    #

    color.use <- color_mapping$Color[color_mapping$MergedTumor==CT1 | color_mapping$MergedTumor==CT2]
    # Custom color scale based on cell proportion
    plot_data <- bin_counts %>% filter(MergedTumor == CT1 | MergedTumor == CT2)

    plot_data <- plot_data %>%
      mutate(Color = ifelse(MergedTumor == CT1, color_mapping$Color[color_mapping$MergedTumor == CT1],
                            color_mapping$Color[color_mapping$MergedTumor == CT2]))
    # Plot the data with custom colors
    p.out <- ggplot(plot_data, aes(x = as.numeric(x_bin), y = as.numeric(y_bin), fill = Color)) +
      geom_tile(aes(alpha = cell_proportion)) +
      scale_fill_identity() +  # Use identity scale to apply colors directly
      scale_alpha_continuous(range = c(0.2, 1), guide = "none") +  # Adjust alpha for transparency
      # labs(
      #   title = paste("Cell Proportion in 3600 Bins for", CT1, "and", CT2),
      #   x = "X Bin",
      #   y = "Y Bin",
      #   fill = "Cell Type"
      # ) +
      scale_y_reverse() +
      theme_void()
    p.out
    save.name <- paste0(ROI_number,"_",CT1,"_",CT2,".tiff")
    ggsave(plot = p.out,filename = file.path(wdpath,"SGCC_relevant_analysis","Bin_phenotype_map",save.name),
           device = "tiff",dpi = 200,width = 2,height = 2)
    if(!all(c(CT1,CT2)%in% unique(plot_data$MergedTumor))){
      next()
    }
    # Set dimensions
    # Create a base scatter plot
    base_data <- expand.grid(x = 1:x_border, y = 1:x_border)
    new_data <- base_data
    # Join plot_data with base_data
    new_data <- as_tibble(new_data) %>%
      left_join(plot_data %>% filter(MergedTumor == CT1) %>% select(x_bin, y_bin, cell_proportion),
                by = c("x" = "x_bin", "y" = "y_bin")) %>%
      mutate(!!CT1:= ifelse(is.na(cell_proportion), 0, cell_proportion)) %>%
      select(-cell_proportion)
    #
    new_data <- as_tibble(new_data) %>%
      left_join(plot_data %>% filter(MergedTumor == CT2) %>% select(x_bin, y_bin, cell_proportion),
                by = c("x" = "x_bin", "y" = "y_bin")) %>%
      mutate(!!CT2:= ifelse(is.na(cell_proportion), 0, cell_proportion)) %>%
      select(-cell_proportion)

    consine_sim <- Cal_GCC(data.in = new_data,knee = output.res[[1]][1], eigenvector = output.res[[2]],signal1 = CT1 , signal2 = CT2)

    SGCC_matrix[tmp_cellcombination.name,ROI_number] <- consine_sim
  }
}

write.csv(SGCC_matrix,file = file.path(wdpath,"SGCC_relevant_analysis", paste0("60_SGCC_matrix_DLBCL_mergeMacro_mergeTumor_mergeT.csv")),quote = F)





