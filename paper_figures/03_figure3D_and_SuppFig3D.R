## set working directory
loadingfunction_wdpath <- c("./src/")
wdpath <- c("./data/SGCC_data/")
# set env
source(file = file.path(loadingfunction_wdpath,"SGCC_code","1-SpaGFT_preload_function.R"))
library(qs)
library(ComplexHeatmap)
library(circlize)
# Figure 3C
filelists <- list.files(file.path(wdpath,"Simulation_Data","RingPattern"),pattern = ".qs")
list_of_movements <- qread(file.path(wdpath,"Simulation_Data","RingPattern",filelists[1]))
## GFT method set k=100, and k_fold = 10 (number of eigen vector)
output.res <- Cal_Eigen(data.in = list_of_movements$`Movement 1`,k = 400,k_fold = 10)
GSCC_df <- matrix(NA,nrow = 10,ncol = length(filelists))
rownames(GSCC_df) <- paste0("Movement ",1:10)
colnames(GSCC_df) <- paste0("radius-",gsub(".*_radius-([0-9.]+)\\.qs", "\\1", filelists))
for (i in 1:length(filelists)){
  tmp.filelist <- filelists[i]
  tmp.colname <- paste0("radius-",gsub(".*_radius-([0-9.]+)\\.qs", "\\1", tmp.filelist))
  #
  list_of_movements <- qread(file.path(wdpath,"Simulation_Data","RingPattern",tmp.filelist))
  for (k in 1:length(list_of_movements)){
    name_movement <- names(list_of_movements)[k]
    test.data <- list_of_movements[[name_movement]]
    tmp_sim <- Cal_GCC(data.in = test.data,knee = output.res[[1]][1], eigenvector = output.res[[2]],
                       signal1 = "inside_circle1",
                       signal2 = "inside_ring")
    GSCC_df[name_movement,tmp.colname] <- tmp_sim
  }
  
}
# heatmap
GSCC_df <- GSCC_df[,paste0("radius-",seq(from = 2.5 ,to= 20, by = 2.5))]
value_to_color <- function(x) {
  colors <- colorRamp2(c(-1, 0, 1), c("#4276b6", "#fcfdfe", "#d93027"))
  colors(x)
}
# Use grid graphics to create density plots for the rows and columns
row_density_plot <- rowAnnotation(
  distribution = anno_boxplot(GSCC_df,which = "row",
                              gp =  gpar(fill = value_to_color(matrixStats::rowMedians(GSCC_df)))),
  width = unit(1.3, "cm")
)
#draw(row_density_plot)
column_density_plot <- HeatmapAnnotation(
  distribution = anno_boxplot(GSCC_df,which = "column",
                              gp =  gpar(fill = value_to_color(matrixStats::colMedians(GSCC_df)))),
  height = unit(1.3, "cm")
)
#draw(column_density_plot)

# Combine heatmap with side bars
ht_list <- Heatmap(as.matrix(GSCC_df),cluster_columns = F,cluster_rows = F,top_annotation = column_density_plot,
                   right_annotation = row_density_plot,
                   col = colorRamp2(c(-0.8, 0,1), c("#4276b6", "#fcfdfe", "#d93027")))
svg("Figure 3C.svg",width = 10,height = 10)
draw(ht_list)
dev.off()

# Supplementary Figure 3C
filelists <- list.files(file.path(wdpath,"Simulation_Data","MovingPattern"),pattern = ".qs")
list_of_movements <- qread(file.path(wdpath,"Simulation_Data","MovingPattern",filelists[1]))
## GFT method set k=100, and k_fold = 10 (number of eigen vector)
output.res <- Cal_Eigen(data.in = list_of_movements$`Movement 1`,k = 400,k_fold = 10)
GSCC_df <- matrix(NA,nrow = 10,ncol = length(filelists))
rownames(GSCC_df) <- paste0("Movement ",1:10)
colnames(GSCC_df) <- paste0("radius-",gsub(".*_radius-([0-9.]+)\\.qs", "\\1", filelists))
for (i in 1:length(filelists)){
  tmp.filelist <- filelists[i]
  tmp.colname <- paste0("radius-",gsub(".*_radius-([0-9.]+)\\.qs", "\\1", tmp.filelist))
  #
  list_of_movements <- qread(file.path(wdpath,"Simulation_Data","MovingPattern",tmp.filelist))
  for (k in 1:length(list_of_movements)){
    name_movement <- names(list_of_movements)[k]
    test.data <- list_of_movements[[name_movement]]
    tmp_sim <- Cal_GCC(data.in = test.data,knee = output.res[[1]][1], eigenvector = output.res[[2]],
                       signal1 = "inside_circle1",
                       signal2 = "inside_circle2")
    GSCC_df[name_movement,tmp.colname] <- tmp_sim
  }
  
}
# heatmap
GSCC_df <- GSCC_df[,paste0("radius-",seq(from =6 ,to= 14, by = 1))]
value_to_color <- function(x) {
  colors <- colorRamp2(c(-1, 0, 1), c("#4276b6", "#fcfdfe", "#d93027"))
  colors(x)
}
# Use grid graphics to create density plots for the rows and columns
row_density_plot <- rowAnnotation(
  distribution = anno_boxplot(GSCC_df,which = "row",
                              gp =  gpar(fill = value_to_color(matrixStats::rowMedians(GSCC_df)))),
  width = unit(1.3, "cm")
)
#draw(row_density_plot)
column_density_plot <- HeatmapAnnotation(
  distribution = anno_boxplot(GSCC_df,which = "column",
                              gp =  gpar(fill = value_to_color(matrixStats::colMedians(GSCC_df)))),
  height = unit(1.3, "cm")
)
#draw(column_density_plot)

# Combine heatmap with side bars
ht_list <- Heatmap(as.matrix(GSCC_df),cluster_columns = F,cluster_rows = F,top_annotation = column_density_plot,
                   right_annotation = row_density_plot,
                   col = colorRamp2(c(-0.8, 0,1), c("#4276b6", "#fcfdfe", "#d93027")))
svg("Supplementary Figure 3C.svg",width = 10,height = 10)
draw(ht_list)
dev.off()