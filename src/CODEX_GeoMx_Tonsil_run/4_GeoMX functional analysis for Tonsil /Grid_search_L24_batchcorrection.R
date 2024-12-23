# set env
library(qs)
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
library(SeuratObject)
library(Seurat)
library(knitr)
library(FNN)
library(limma)
library(scater)
library(ExperimentHub)
library(ggalluvial)
library(kBET)
# set up working directory
wdpath <- "./data/Tonsil_run/"
# read cell data
cell_data <- read.csv(file.path(wdpath,"Annotation/L24_final_annotation_ROIlevel.csv"))
roi_centers <- cell_data %>%
  group_by(ROI_num) %>%
  summarize(center_x = median(X_cent),
            center_y = median(Y_cent))
roi_centers_clean <- cell_data %>%
  group_by(ROI_num) %>%
  summarize(center_x = median(X_cent),
            center_y = median(Y_cent)) %>% 
  ungroup() %>% 
  mutate(ROI_indx = sprintf("%03d",as.numeric(gsub("ROI_","",ROI_num))))
# summarize cell data
summarized_data <- cell_data %>%
  group_by(ROI_num, Annotation5) %>%
  summarize(
    cell_count_CODEX = n(),
    median_cell_size = median(cellSize, na.rm = TRUE),
    median_dapi = median(DAPI, na.rm = TRUE),
    median_cd11b = median(CD11b, na.rm = TRUE),
    median_cd68 = median(CD68, na.rm = TRUE),
    median_cd11c = median(CD11c, na.rm = TRUE),
    median_pax5 = median(Pax5, na.rm = TRUE),
    median_cd20 = median(CD20, na.rm = TRUE),
    median_bcl6 = median(BCL6, na.rm = TRUE),
    median_cd8 = median(CD8, na.rm = TRUE),
    median_cd163 = median(CD163, na.rm = TRUE),
    median_cd4 = median(CD4, na.rm = TRUE),
    median_foxp3 = median(FoxP3, na.rm = TRUE),
    median_cd31 = median(CD31, na.rm = TRUE),
    median_cd3 = median(CD3, na.rm = TRUE)
  ) %>% 
  ungroup() %>% group_by(ROI_num) %>% 
  mutate(TotalCellPerROI = sum(cell_count_CODEX))
summarized_data <- full_join(summarized_data,roi_centers_clean, by = "ROI_num")


#read data
mysegmentproperty <- read_excel(file.path(wdpath,"GeoMx_count_table/L24_GeoMXMeta.xlsx"),sheet = "SegmentProperties")
mycountmat <- read_excel(file.path(wdpath,"GeoMx_count_table/L24_GeoMXRawData.xlsx"),sheet = "BioProbeCountMatrix")
# move neagtive gene to meta negative to meta data of ROI
mycountmat_duplicated <- mycountmat[duplicated(mycountmat$TargetName),] 
mycountmat_clean <- mycountmat[!duplicated(mycountmat$TargetName),] 
#####
mygeneannotation <- mycountmat_clean[,1:12]
mycountmat_clean <- mycountmat_clean[,c(13:ncol(mycountmat_clean))]
# map cell annotaiton from CODEX and GeoMX
mapping_table <- data.frame(
  SegmentLabel = c("BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "DC", "Endo", "Full ROI", "M1", "M2", "Myeloid", "Other"),
  Annotation5 = c("BCL6- B Cell", "BCL6+ B Cell", "CD4 T", "CD4 Treg", "CD8 T", "DC", "Endothelial", "Full ROI", "M1", "M2", "Myeloid", "Other")
)
mysegmentproperty_mapped <- mysegmentproperty %>%
  left_join(mapping_table, by = "SegmentLabel")
# add codex data summarie to AOI
merged_table <- mysegmentproperty_mapped %>%
  full_join(summarized_data, by = c("ROILabel" = "ROI_indx", "Annotation5" = "Annotation5"))
merged_table <- merged_table[!is.na(merged_table$SlideName),]
# create spatial object based on SpatialExperiment
SegmentMeta <- merged_table[,c("SegmentDisplayName","SegmentLabel","AOISurfaceArea","AOINucleiCount",
                               "RawReads","AlignedReads","DeduplicatedReads", "TrimmedReads","StitchedReads" , "SequencingSaturation",
                               "Annotation5","cell_count_CODEX", "median_cell_size" ,"median_dapi", "median_cd11b" ,"median_cd68" ,"median_cd11c",
                               "median_pax5", "median_cd20","median_bcl6" ,"median_cd8" , "median_cd163",
                               "median_cd4" ,"median_foxp3", "median_cd31","median_cd3","TotalCellPerROI", "center_x","center_y","ROILabel")]

# Predefined some functions for calculating silhoutte and Kbet score (Run whole chunk)
# Silhoutte
cal_sil <- function(se.object = spe.list[[names(spe.list)[1]]],
                    #obj.name = names(spe.list)[1],
                    batch.name = "patient"){
  data <- assay(se.object,i = 2)
  # find HVG
  dec <-scran:: modelGeneVar(data)
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes,])
  #
  batch <- colData(se.object)[,batch.name]
  pca.data.f <- gmodels::fast.prcomp(data_hvg, center=TRUE)
  dd <- as.matrix(dist(pca.data.f$x[, 1:10]))
  batch.silhouette <- summary(cluster::silhouette(as.numeric(factor(batch,
                                                                    levels = sort(unique(batch)),
                                                                    labels = 1:length(sort(unique(batch))))), dd))$avg.width
  return(batch.silhouette)
}

# kBET
cal_kbet <- function(se.object = spe.list[[names(spe.list)[1]]],
                     # obj.name = names(spe.list)[1],
                     batch.name = "Annotation5"){
  data <- assay(se.object,i = 2)
  # find HVG
  dec <-scran:: modelGeneVar(data)
  top_genes <- scran::getTopHVGs(dec, n = 1000)
  data_hvg <- t(data[rownames(data) %in% top_genes,])
  #
  batch <- colData(se.object)[,batch.name]
  k0=floor(mean(table(batch))) #neighbourhood size: mean batch size 
  knn <- get.knn(data_hvg, k=k0, algorithm = 'cover_tree')
  batch.estimate <- kBET(data_hvg, batch, k = k0, knn = knn,plot = F)
  return(batch.estimate)
}


# Set up benchmarking parameter for RUV4 functions (Action required for user)

# To accelerate the execution speed of the tutorial, we only consider four types of normalization here. Users will need to set the remaining parameters themselves, including wanted variable, unwanted variable, number of negative control genes, and k from RUV4.
# 
# * factors: it indicate wanted factors and unwanted factors, i.e. the biological variation to keep;
# 
# * NCGs: the list of negative control genes detected using the function findNCGs;
# 
# * k: the number of unwanted factors to use. Based on RUV’s documentation, it is suggest to use the smallest k possible where the observed technical variation is no longer observed.
# define unwanted factor
unwantedVariable <- c("ROILabel")
# define wanted factor
WantedVariable <- c("SegmentLabel")
normalization_methods <- c("TMM", "CPM", "upperquartile", "sizefactor")
# number of negative control gene
findNCG_topn <- seq(from = 1000, to = 3000, by = 500)
# k of RUV4 (The number of unwanted factors to use. Can be 0. This is required for the RUV4 method.)
BatchCorrection_k <- seq(from = 3, to = 9, by = 2)
# show the number of combinations
print(paste0("Combination Parameter Number: ",length(BatchCorrection_k)*length(findNCG_topn)*length(normalization_methods)))

# run batch correction
# create spatial Experiment obj
mycountmat_clean <- as.matrix(mycountmat_clean)
rownames(mycountmat_clean) <- mygeneannotation$TargetName
spe <- SpatialExperiment(
  assay = list(counts = as.matrix(mycountmat_clean)),
  rowData = mygeneannotation,
  colData = SegmentMeta,
  spatialCoords = as.matrix(data.frame(X = SegmentMeta$center_x, Y = SegmentMeta$center_y))
)
# save processed raw counts object 
# qsave(spe, file.path(wdpath,"GeoMX_L24_L25_batchcorrection/L24_spe_raw_counts.qs"))

spe.list <- list()
# Total number of iterations
total_iterations <- length(normalization_methods) * length(findNCG_topn) * length(BatchCorrection_k)

# Initialize progress bar
progress <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration_count <- 0

for (n_m in 1:length(normalization_methods)) {
  tmp_n_m <- normalization_methods[n_m]
  tmp_spe <- geomxNorm(spe, method = tmp_n_m)
  if(tmp_n_m == "CPM"){tmp_spe <- scater::logNormCounts(spe)}
  
  for (topn in 1:length(findNCG_topn)) {
    tmp_topn <- findNCG_topn[topn]
    #### here need to set up the batch factor based on meta data
    tmp_spe <- findNCGs(tmp_spe, batch_name = unwantedVariable, top_n = tmp_topn)
    
    for (k in 1:length(BatchCorrection_k)) {
      tmp_k <- BatchCorrection_k[k]
      # [to do] change factors to biological variation remain based on meta data
      # here to setup wanted biology variable
      tmp_spe <- geomxBatchCorrection(tmp_spe, factors = WantedVariable,
                                      NCGs = metadata(tmp_spe)$NCGs, k = tmp_k)
      tmp_name <- paste0(tmp_n_m, "-nNCG_", tmp_topn, "-k_", tmp_k)
      spe.list <- c(spe.list, tmp_spe)
      names(spe.list)[length(spe.list)] <- tmp_name
      
      # Update progress bar
      iteration_count <- iteration_count + 1
      setTxtProgressBar(progress, iteration_count)
    }
  }
}

# Close progress bar
close(progress)



# evaluation using silouette and Kbet (Run whole chunk)
#
# run for each of factor, including patientID (batch confounder), celltype (biological variation), and Disease (biological variation)
print("calculated batch score")
batch.name = c(unwantedVariable, WantedVariable)
sil_score_list <- list()
k_bet_list <- list()
# progress bar
total_iterations <- length(spe.list)*length(batch.name)

# Initialize progress bar
progress <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration_count <- 0
for (i in 1:length(spe.list)){
  se.object = spe.list[[names(spe.list)[i]]]
  obj.name = names(spe.list)[i]
  for (j in 1:length(batch.name)){
    sil_score <- cal_sil(se.object = se.object, batch.name = batch.name[j] )
    k_bet.obj <- cal_kbet(se.object = se.object, batch.name = batch.name[j] )
    sil_score_list <- c(sil_score_list, sil_score)
    k_bet_list <- c(k_bet_list, mean(k_bet.obj$stats$kBET.observed))
    names(sil_score_list)[length(sil_score_list)] <-  names(k_bet_list)[length(k_bet_list)] <- paste0(batch.name[j],"_",obj.name)
    # Update progress bar
    iteration_count <- iteration_count + 1
    setTxtProgressBar(progress, iteration_count)
  }
}


# evaluation of all factors
# minimize unwanted variance
unwantedVariable_df <- matrix(rep(NA,length(unwantedVariable)*length(spe.list)*5),ncol =length(unwantedVariable)*5 )

unwantedVariable_combinations <- expand.grid(unwantedVariable,  c("_sil","_kbet","_Silrank","_kbetrank","_MeanRank"))

colnames(unwantedVariable_df) <- paste0(unwantedVariable_combinations$Var1, unwantedVariable_combinations$Var2)

rownames(unwantedVariable_df) <-  names(spe.list)
for (i in 1:length(unwantedVariable)){
  tmp_unwantedVariable <- unwantedVariable[i]
  tmp_kbet_name_vec <- grep(tmp_unwantedVariable,(names(k_bet_list)),value = F)
  tmp_kbet <-  unlist(k_bet_list[tmp_kbet_name_vec])
  tmp_sil_name_vec <- grep(tmp_unwantedVariable,(names(sil_score_list)),value = F)
  tmp_sil <-  unlist(sil_score_list[tmp_sil_name_vec])
  
  #
  names(tmp_sil) <- names(tmp_kbet) <- gsub(paste0("^",tmp_unwantedVariable,"_"),"",names(tmp_sil))
  tmp_unwantedVariable_kbet_rank <- rank(tmp_kbet)
  tmp_unwantedVariable_sil_rank <- rank(tmp_sil)
  # mean rank of batch
  tmp_unwantedVariable_kbet_rank <- tmp_unwantedVariable_kbet_rank[rownames(unwantedVariable_df)]
  tmp_unwantedVariable_sil_rank <- tmp_unwantedVariable_sil_rank[rownames(unwantedVariable_df)]
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_sil"))] <- as.numeric(tmp_sil)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_kbet"))] <- as.numeric(tmp_kbet)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_Silrank"))] <- as.numeric(tmp_unwantedVariable_sil_rank)
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_kbetrank"))] <- as.numeric(tmp_unwantedVariable_kbet_rank)
  # cal mean rank of sil and kbet
  unwantedVariable_df[,paste0(tmp_unwantedVariable,c("_MeanRank"))] <- rowSums(unwantedVariable_df[,c(paste0(tmp_unwantedVariable,c("_Silrank")),
                                                                                                      paste0(tmp_unwantedVariable,c("_kbetrank")))])/2
}

# Maximize wanted variance
WantedVariable_df <- matrix(rep(NA,length(WantedVariable)*length(spe.list)*5),ncol =length(WantedVariable)*5 )
WantedVariable_combinations <- expand.grid(WantedVariable, c("_sil","_kbet","_Silrank","_kbetrank","_MeanRank"))

colnames(WantedVariable_df) <- paste0(WantedVariable_combinations$Var1, WantedVariable_combinations$Var2)
rownames(WantedVariable_df) <-  names(spe.list)
for (i in 1:length(WantedVariable)){
  tmp_WantedVariable <- WantedVariable[i]
  tmp_kbet_name_vec <- grep(tmp_WantedVariable,(names(k_bet_list)),value = F)
  tmp_kbet <-  unlist(k_bet_list[tmp_kbet_name_vec])
  tmp_sil_name_vec <- grep(tmp_WantedVariable,(names(sil_score_list)),value = F)
  tmp_sil <-  unlist(sil_score_list[tmp_sil_name_vec])
  
  #
  names(tmp_sil) <- names(tmp_kbet) <- gsub(paste0("^",tmp_WantedVariable,"_"),"",names(tmp_sil))
  tmp_WantedVariable_kbet_rank <- rank(-tmp_kbet)
  tmp_WantedVariable_sil_rank <- rank(-tmp_sil)
  # mean rank of batch
  tmp_WantedVariable_kbet_rank <- tmp_WantedVariable_kbet_rank[rownames(WantedVariable_df)]
  tmp_WantedVariable_sil_rank <- tmp_WantedVariable_sil_rank[rownames(WantedVariable_df)]
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_sil"))] <- as.numeric(tmp_sil)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_kbet"))] <- as.numeric(tmp_kbet)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_Silrank"))] <- as.numeric(tmp_WantedVariable_sil_rank)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_kbetrank"))] <- as.numeric(tmp_WantedVariable_kbet_rank)
  WantedVariable_df[,paste0(tmp_WantedVariable,c("_MeanRank"))] <- rowSums(WantedVariable_df[,c(paste0(tmp_WantedVariable,c("_Silrank")),
                                                                                                paste0(tmp_WantedVariable,c("_kbetrank")))])/2
}


# select your parameters (Action required for user)
# merge results of wanted variable and unwanted variable
merge_res <- cbind.data.frame(unwantedVariable_df,WantedVariable_df)
merge_res <- merge_res[,grep("_MeanRank",colnames(merge_res))]
merge_res$overall_rank <- rowSums(merge_res[,grep("_MeanRank",colnames(merge_res))])/length(merge_res[,grep("_MeanRank",colnames(merge_res))]) 
merge_res <- merge_res[order(merge_res$overall_rank),]
Save.name_par <- file.path(wdpath,"GeoMX_L24_L25_batchcorrection","L24_para.csv")
# Save parameters
# write.csv(merge_res,file = Save.name_par,quote = F)












