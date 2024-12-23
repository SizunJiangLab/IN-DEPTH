loadingfunction_wdpath <- c("./src/")
wdpath <- "./data/DLBCL_run/"
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(deldir)
library(igraph)
library(progress)
library(qs)
library(pheatmap)
library(epitools)  # For odds ratio calculation
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","2-rank_SGCC_Impulse_preload_function.R"))
source(file.path(loadingfunction_wdpath,"CODEX_GeoMX analysis for DLBCL","3-WithinCrossDomainTableAndVisualization_preload_function.R"))
# setup your cell type color
cell_type_colors <- c("Tumor" = "#ed1c2c",   # Red
                      "CD8T" = "#0082bf",   # Green
                      "CD4T" = "#2db44e",   # Blue
                      "Treg" = "#85c540",    # Light Green
                      "Macro" = "#8f52a0", # purple
                      "DC" = "#66cbdb",      # Blue-Violet
                      "Endothelial" = "#a87a66") # Brown
# claim your cell type pairs
# Celltypepairs <- c("Macro","CD4T")
# Celltypepairs <- c("Macro","Tumor")
Celltypepairslist <- list(c("CD4T","Tumor"),c("Macro","Tumor"),c("Macro","CD4T"))
for (i in 1:3){
  Celltypepairs <- Celltypepairslist[[i]]
  # Load the ROI annotation data
  annotation <- read_csv(file.path(wdpath,"SGCC_relevant_analysis","DLBCL_ROIlevel_annotation.csv"))
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
  
  ### run SGCC score vs cell edges
  stat.result <- list()
  for (core in unique(annotation$coreName)) {
    print(paste("Working on", core, "..."))
    
    # Subset data for the current core
    sub_df <- annotation %>%
      filter(coreName == core) 
    result <- generate_edge_data(annotation = sub_df,cell_types = c(Celltypepairs[1], Celltypepairs[2]),pixel_to_um = 0.5)
    stat.result.tmp <- (calculate_odds(cell_data = result$cell_data,
                                       edge_data =result$edge_data ,
                                       ct1_name = Celltypepairs[1],
                                       ct2_name = Celltypepairs[2]))
    stat.result[[core]] <- stat.result.tmp
  }
  
  # Create a data frame where each column corresponds to a different element of stat.result (e.g., "DFCI_4.1")
  stat.result.df <- do.call(cbind, lapply(stat.result, function(x) {
    c(unlist(x))  # Convert each list component into a vector (unlisting if necessary)
  }))
  
  # Convert the matrix to a data frame if not automatically converted
  stat.result.df <- as.data.frame(stat.result.df)
  
  # Set row names to reflect the features from any single element of stat.result
  row.names(stat.result.df) <- names(stat.result$DFCI_4.1)
  ##
  stat.result.df.t <- as.data.frame(t(stat.result.df))
  
  stat.result.df.t$AES <- stat.result.df.t$edge_across/(2* stat.result.df.t$CT1_n*stat.result.df.t$CT2_n*stat.result.df.t$total_edges/(stat.result.df.t$total_n)^2)-1
  
  qsave(stat.result.df.t,
        file = file.path(wdpath,"SGCC_relevant_analysis",
                      paste0(Celltypepairs[1],"_", Celltypepairs[2],"stat.result.qs")))
}


