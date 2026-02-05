## set working directory
loadingfunction_wdpath <- c("./src/")
# set env
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
library(qs)
source(file = file.path(loadingfunction_wdpath,"SGCC_code","1-SpaGFT_preload_function.R"))
############
wdpath <- "./data/Tonsil_run/"
cell_data <- read.csv(file.path(wdpath,"Annotation/L25_final_annotation_ROIlevel.csv"))
# change names
mapping_table <- data.frame(
  SegmentLabel = c("BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "DC", "Endo", "Full ROI", "M1", "M2", "Myeloid", "Other"),
  Annotation7 = c("BCL6- B Cell", "BCL6+ B Cell", "CD4 T", "CD4 Treg", "CD8 T", "DC", "Endothelial", "Full ROI", "M1", "M2", "Myeloid", "Other")
)
cell_data <- cell_data %>%
  left_join(mapping_table, by = c("Annotation7" = "Annotation7"))
# Filter for a specific cell type, for example, "CD4 T"
# "BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "DC", "Endo", "Full ROI", "M1", "M2", "Myeloid", "Other"
# Create a color mapping table
fill_colors <- c("#e07329", 
                 "#16964a",  
                 "#2958a8",  
                 "#dd2246",  
                 "#703594",  
                 "#d0bd2a",
                 "#69bd48",  
                 "#2297b1",  
                 "#b9e4f3", 
                 "#ae5c71",
                 "#ab2170")

color_mapping <- data.frame(
  SegmentLabel = c("BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "DC", "Endo", "M1", "M2", "Myeloid", "Other"),
  Color = fill_colors
)
# Create a named vector for the color mapping
color_mapping_vector <- setNames(color_mapping$Color, color_mapping$SegmentLabel)
# Create a summary data frame with ROI center coordinates
roi_centers_clean <- cell_data %>%
  group_by(ROI_num) %>%
  summarize(center_x = median(X_cent),
            center_y = median(Y_cent)) %>% 
  ungroup() %>% 
  mutate(ROILabel = sprintf("%03d",as.numeric(gsub("ROI_","",ROI_num))))


roi_centers_clean_list <- split(cell_data, cell_data$ROI_num)
# Set dimensions
x_border <- 60
# Create a base grid
base_data <- expand.grid(x = 1:x_border, y = 1:x_border)
output.res <- Cal_Eigen(data.in = base_data,k = 400,k_fold = 10)
#
cellpairs <- combn(c("BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "DC", "Endo", "M1", "M2", "Myeloid", "Other"), 2)
cellcombination.name <- apply(cellpairs, 2, function(x) paste(x, collapse = "_"))
# create SGCC matrix
SGCC_matrix <- matrix(NA, nrow = dim(cellpairs)[2], ncol = length(roi_centers_clean_list))
rownames(SGCC_matrix) <- cellcombination.name
colnames(SGCC_matrix) <- paste0("ROI_",0:15)
for (i in paste0("ROI_",0:15)){
  ROI_number <- i
  data <- roi_centers_clean_list[[ROI_number]]
  
  # Create 3600 bins (60x60 grid)
  x_bins <- x_border
  y_bins <- x_border
  
  data <- data %>%
    mutate(
      x_bin = cut(X_cent, breaks = x_bins, labels = FALSE),
      y_bin = cut(Y_cent, breaks = y_bins, labels = FALSE)
    )
  
  # Count the number of cells of each type in each bin
  bin_counts <- data %>%
    group_by(x_bin, y_bin, SegmentLabel) %>%
    summarise(cell_count = n(), .groups = 'drop')
  
  # Calculate the total number of cells of each type
  total_cells_per_type <- data %>%
    group_by(SegmentLabel) %>%
    summarise(total_count = n(), .groups = 'drop')
  
  # Merge total counts with bin counts
  bin_counts <- bin_counts %>%
    left_join(total_cells_per_type, by = "SegmentLabel") %>%
    mutate(cell_proportion = cell_count / total_count)
  
  # Function to adjust color transparency based on proportion
  for (k in 1:length(cellcombination.name)){
    # Filter for a specific cell type, for example, "CD4T"
    tmp_cellcombination.name <- cellcombination.name[k]
    CT1 <- strsplit(tmp_cellcombination.name,"_")[[1]][1]
    
    CT2 <- strsplit(tmp_cellcombination.name,"_")[[1]][2]
    
    color.use <- color_mapping$Color[color_mapping$SegmentLabel==CT1 | color_mapping$SegmentLabel==CT2]
    # Custom color scale based on cell proportion
    plot_data <- bin_counts %>% filter(SegmentLabel == CT1 | SegmentLabel == CT2)
    plot_data <- plot_data %>%
      mutate(Color = ifelse(SegmentLabel == CT1, color_mapping$Color[color_mapping$SegmentLabel == CT1], 
                            color_mapping$Color[color_mapping$SegmentLabel == CT2]))
    # Plot the data with custom colors
    p.out <- ggplot(plot_data, aes(x = as.numeric(x_bin), y = as.numeric(y_bin), fill = Color)) +
      geom_tile(aes(alpha = cell_proportion)) +
      scale_fill_identity() +  # Use identity scale to apply colors directly
      scale_alpha_continuous(range = c(0.2, 1), guide = "none") +  # Adjust alpha for transparency
      scale_y_reverse() +
      theme_void()
    p.out
    save.name <- paste("L25",ROI_number,CT1,CT2,".tiff",sep = "_")
    ggsave(plot = p.out,filename = file.path(wdpath,"SGCC_relevant_analysis","Bin_phenotype_map",save.name),
           device = "tiff",dpi = 150,width = 2,height = 2)
    # 
    if(!all(c(CT1,CT2)%in% unique(plot_data$SegmentLabel))){
      next()
    }
    # Create a base scatter plot
    base_data <- expand.grid(x = 1:x_border, y = 1:x_border)
    new_data <- base_data
    # Join plot_data with base_data
    new_data <- as_tibble(new_data) %>%
      left_join(plot_data %>% filter(SegmentLabel == CT1) %>% dplyr::select(x_bin, y_bin, cell_proportion), 
                by = c("x" = "x_bin", "y" = "y_bin")) %>%
      mutate(!!CT1:= ifelse(is.na(cell_proportion), 0, cell_proportion)) %>%
      dplyr::select(-cell_proportion)
    # 
    new_data <- as_tibble(new_data) %>%
      left_join(plot_data %>% filter(SegmentLabel == CT2) %>% dplyr::select(x_bin, y_bin, cell_proportion), 
                by = c("x" = "x_bin", "y" = "y_bin")) %>%
      mutate(!!CT2:= ifelse(is.na(cell_proportion), 0, cell_proportion)) %>%
      dplyr::select(-cell_proportion)
    
    consine_sim <- Cal_GCC(data.in = new_data,knee = output.res[[1]][1], eigenvector = output.res[[2]],signal1 = CT1 , signal2 = CT2) 
    
    SGCC_matrix[tmp_cellcombination.name,ROI_number] <- consine_sim
  }
}
SGCC_matrix <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","SGCC_matrix_L25.csv"),row.names = 1)
write.csv(SGCC_matrix,file = file.path(wdpath,"SGCC_relevant_analysis","SGCC_matrix_L25.csv"),quote = F)

# one way anova
# Convert matrix to data frame
SGCC_df <- as.data.frame(SGCC_matrix)

# Create ROI factor
roi_factor <- as.factor(colnames(SGCC_df))

####### gene expression data SGCC analysis
# read in spe object 
spe <- qread(file.path(wdpath,"GeoMX_L24_L25_batchcorrection/L25_spe_raw_counts.qs"))
## add ROI number 
# subset ROI_num and ROILabel
roi_subset <- select(roi_centers_clean, ROI_num, ROILabel)
# Convert colData to a data.frame for manipulation
current_metadata <- as.data.frame(colData(spe))
# Perform the left join to update metadata
updated_metadata <- left_join(current_metadata, roi_subset, by = "ROILabel")
# Update the colData in spe
colData(spe) <- DataFrame(updated_metadata)
#select BCL6pB_CD4T
group.name <- "BCL6pB_CD4T"
CT1 = strsplit(group.name,"_")[[1]][1]
CT2 =strsplit(group.name,"_")[[1]][2]

# Select a specific cell type (CT1) for analysis
CT.use <- CT1
# Subset the SingleCellExperiment object to include only the samples with the specified cell type
spe_sub <- spe[, colData(spe)$SegmentLabel == CT.use]
# Extract ROI_num from the column data of the subsetted object
groups <- colData(spe_sub)[, c("ROI_num")]

# rename matrix
pb <- t(aggregate.Matrix(t(counts(spe_sub)), groupings = groups, fun = "sum"))

# Calculate scores, assuming SGCC_df is previously defined and contains scores for each group
SGCC_score = as.numeric(SGCC_df[group.name, colnames(pb)])

# Create an annotation DataFrame with metadata for ImpulseDE2 analysis
dfAnnotation <- data.frame(
  row.names = colnames(pb),
  Sample = colnames(pb),
  Condition = "case",
  SGCC = SGCC_score,
  Time = as.numeric(cut(SGCC_score, 
                        breaks = quantile(SGCC_score, probs = seq(0, 1, by = 1/3)), 
                        labels = 1:3, 
                        include.lowest = TRUE))
)
# Order the DataFrame by SGCC scores
dfAnnotation <- dfAnnotation[order(dfAnnotation$SGCC),]

# Reorder pb matrix columns based on the ordered annotation
pb <- pb[, dfAnnotation$Sample]

# Run ImpulseDE2 to identify differentially expressed genes over time or condition
objectImpulseDE2 <- runImpulseDE2(
  matCountData = as.matrix(pb), 
  dfAnnotation = dfAnnotation,
  scaQThres = 0.1,
  boolCaseCtrl = FALSE,
  boolIdentifyTransients = FALSE,
  scaNProc = 10
)


# Create density plots for the rows and columns of the heatmap
column_density_plot <- HeatmapAnnotation(
  SGCC = anno_lines(dfAnnotation$SGCC, which = "column"),
  gp = gpar(lwd = 10),
  height = unit(1.3, "cm")
)

# Define a function to evaluate Impulse model
evalImpulse <- function(vecImpulseParam, vecTimepoints) {
  vecImpulseValue <- sapply(vecTimepoints, function(t) {
    (1/vecImpulseParam[3]) * 
      (vecImpulseParam[2] + (vecImpulseParam[3] - vecImpulseParam[2]) *
         (1/(1 + exp(-vecImpulseParam[1] * (t - vecImpulseParam[5]))))) *
      (vecImpulseParam[4] + (vecImpulseParam[3] - vecImpulseParam[4]) *
         (1/(1 + exp(vecImpulseParam[1] * (t - vecImpulseParam[6])))))
  })
  vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
  
  return(vecImpulseValue)
}
# Compile the evalImpulse function to improve performance
evalImpulse_comp <- compiler::cmpfun(evalImpulse)

# Retrieve processed annotation and count data
dfAnnot <- get_dfAnnotationProc(obj=objectImpulseDE2)
scaNGenes <- dim(get_matCountDataProc(obj=objectImpulseDE2))[1]

# Order genes by time of extremum (peak/valley) occurrence
vecSignificantIDs <- rownames(objectImpulseDE2$dfImpulseDE2Results[
  !is.na(objectImpulseDE2$dfImpulseDE2Results$padj) & 
    objectImpulseDE2$dfImpulseDE2Results$p < 0.01, ])
vecTimePointsToEval <- sort(unique(dfAnnot$Time), decreasing = FALSE)
scaNTPtoEvaluate <- length(vecTimePointsToEval)

# Calculate Impulse values across significant genes
matImpulseValue <- do.call(rbind, lapply(
  vecSignificantIDs, function(x) {
    evalImpulse_comp(
      vecImpulseParam = 
        get_lsModelFits(obj=objectImpulseDE2)[["case"]][[x]]$
        lsImpulseFit$vecImpulseParam, 
      vecTimepoints = vecTimePointsToEval)
  }))
rownames(matImpulseValue) <- vecSignificantIDs

# Sort genes by time points
matidxMaxTimeSort <- t(apply(matImpulseValue, 1, function(genevalues) {
  sort(genevalues, decreasing = TRUE, index.return = TRUE)$ix
}))
vecMaxTime <- vecTimePointsToEval[matidxMaxTimeSort[, 1]]
matidxMinTimeSort <- t(apply(matImpulseValue, 1, function(genevalues) {
  sort(genevalues, decreasing = FALSE, index.return = TRUE)$ix
}))

# Analyze time differences to identify trends in gene expression
Time_DIFF <- matidxMinTimeSort[,3]-matidxMinTimeSort[,1]
increase_gene <- rownames(matidxMinTimeSort)[Time_DIFF>0]
Decrease_gene <- rownames(matidxMinTimeSort)[Time_DIFF<0]

# Create output data frame with results
Impulse_res_out <- data.frame(
  T1 = matidxMinTimeSort[,1],
  T2 = matidxMinTimeSort[,2],
  T3 = matidxMinTimeSort[,3],
  Time_Difference = Time_DIFF,
  monotone = ifelse(Time_DIFF > 0, "Increase", ifelse(Time_DIFF < 0, "Decrease", "NA"))
)

# Write results to CSV file
write.csv(Impulse_res_out,
          file = file.path(wdpath, "SGCC_relevant_analysis", paste0("L25", "--", CT.use, "--", group.name, "_Gene.csv")),
          quote = F)

# Generate log-transformed counts per million matrix
norm_mat_B <- edgeR::cpm(pb, log = TRUE)

# Select a specific cell type (CT2) for analysis
CT.use <- CT2
# Subset the SingleCellExperiment object to include only the samples with the specified cell type
spe_sub <- spe[, colData(spe)$SegmentLabel == CT.use]
# Extract ROI_num from the column data of the subsetted object
groups <- colData(spe_sub)[, c("ROI_num")]

# Rename matrix
pb <- t(aggregate.Matrix(t(counts(spe_sub)), groupings = groups, fun = "sum"))

# Calculate scores, assuming SGCC_df is previously defined and contains scores for each group
SGCC_score = as.numeric(SGCC_df[group.name, colnames(pb)])

# Create an annotation DataFrame with metadata for ImpulseDE2 analysis
dfAnnotation <- data.frame(
  row.names = colnames(pb),
  Sample = colnames(pb),
  Condition = "case",
  SGCC = SGCC_score,
  Time = as.numeric(cut(SGCC_score, 
                        breaks = quantile(SGCC_score, probs = seq(0, 1, by = 1/3)), 
                        labels = 1:3, 
                        include.lowest = TRUE))
)
# Order the DataFrame by SGCC scores
dfAnnotation <- dfAnnotation[order(dfAnnotation$SGCC),]

# Reorder pb matrix columns based on the ordered annotation
pb <- pb[, dfAnnotation$Sample]

# Run ImpulseDE2 to identify differentially expressed genes over time or condition
objectImpulseDE2 <- runImpulseDE2(
  matCountData = as.matrix(pb), 
  dfAnnotation = dfAnnotation,
  scaQThres = 0.1,
  boolCaseCtrl = FALSE,
  boolIdentifyTransients = FALSE,
  scaNProc = 10
)


# Create density plots for the rows and columns of the heatmap
column_density_plot <- HeatmapAnnotation(
  SGCC = anno_lines(dfAnnotation$SGCC, which = "column"),
  gp = gpar(lwd = 10),
  height = unit(1.3, "cm")
)

# Define a function to evaluate Impulse model
evalImpulse <- function(vecImpulseParam, vecTimepoints) {
  vecImpulseValue <- sapply(vecTimepoints, function(t) {
    (1/vecImpulseParam[3]) * 
      (vecImpulseParam[2] + (vecImpulseParam[3] - vecImpulseParam[2]) *
         (1/(1 + exp(-vecImpulseParam[1] * (t - vecImpulseParam[5]))))) *
      (vecImpulseParam[4] + (vecImpulseParam[3] - vecImpulseParam[4]) *
         (1/(1 + exp(vecImpulseParam[1] * (t - vecImpulseParam[6])))))
  })
  vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
  
  return(vecImpulseValue)
}
# Compile the evalImpulse function to improve performance
evalImpulse_comp <- compiler::cmpfun(evalImpulse)

# Retrieve processed annotation and count data
dfAnnot <- get_dfAnnotationProc(obj=objectImpulseDE2)
scaNGenes <- dim(get_matCountDataProc(obj=objectImpulseDE2))[1]

# Order genes by time of extremum (peak/valley) occurrence
vecSignificantIDs <- rownames(objectImpulseDE2$dfImpulseDE2Results[
  !is.na(objectImpulseDE2$dfImpulseDE2Results$padj) & 
    objectImpulseDE2$dfImpulseDE2Results$p < 0.01, ])
vecTimePointsToEval <- sort(unique(dfAnnot$Time), decreasing = FALSE)
scaNTPtoEvaluate <- length(vecTimePointsToEval)

# Calculate Impulse values across significant genes
matImpulseValue <- do.call(rbind, lapply(
  vecSignificantIDs, function(x) {
    evalImpulse_comp(
      vecImpulseParam = 
        get_lsModelFits(obj=objectImpulseDE2)[["case"]][[x]]$
        lsImpulseFit$vecImpulseParam, 
      vecTimepoints = vecTimePointsToEval)
  }))
rownames(matImpulseValue) <- vecSignificantIDs

# Sort genes by time points
matidxMaxTimeSort <- t(apply(matImpulseValue, 1, function(genevalues) {
  sort(genevalues, decreasing = TRUE, index.return = TRUE)$ix
}))
vecMaxTime <- vecTimePointsToEval[matidxMaxTimeSort[, 1]]
matidxMinTimeSort <- t(apply(matImpulseValue, 1, function(genevalues) {
  sort(genevalues, decreasing = FALSE, index.return = TRUE)$ix
}))

# Analyze time differences to identify trends in gene expression
Time_DIFF <- matidxMinTimeSort[,3]-matidxMinTimeSort[,1]
increase_gene <- rownames(matidxMinTimeSort)[Time_DIFF>0]
Decrease_gene <- rownames(matidxMinTimeSort)[Time_DIFF<0]

# Create output data frame with results
Impulse_res_out <- data.frame(
  T1 = matidxMinTimeSort[,1],
  T2 = matidxMinTimeSort[,2],
  T3 = matidxMinTimeSort[,3],
  Time_Difference = Time_DIFF,
  monotone = ifelse(Time_DIFF > 0, "Increase", ifelse(Time_DIFF < 0, "Decrease", "NA"))
)

# Write results to CSV file
write.csv(Impulse_res_out,
          file = file.path(wdpath, "SGCC_relevant_analysis", paste0("L25", "--", CT.use, "--", group.name, "_Gene.csv")),
          quote = F)

# Generate log-transformed counts per million matrix
norm_mat_T <- edgeR::cpm(pb, log = TRUE)

# plot select heatmap 
gene_selected <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","BCL6pB_CD4T_Gene_select_from_Impulse_prediction.csv"))

matrix_merge <- rbind(norm_mat_T[gene_selected$Gene[gene_selected$CellType=="CD4T"],],
                      norm_mat_B[gene_selected$Gene[gene_selected$CellType=="BCL6pB"],])

qs::qsave(norm_mat_B, file = file.path(wdpath,"SGCC_relevant_analysis","norm_mat_B_L25.qs"))
qs::qsave(norm_mat_T, file = file.path(wdpath,"SGCC_relevant_analysis","norm_mat_T_L25.qs"))
qs::qsave(matrix_merge, file = file.path(wdpath,"SGCC_relevant_analysis","matrix_merge_L25.qs"))










