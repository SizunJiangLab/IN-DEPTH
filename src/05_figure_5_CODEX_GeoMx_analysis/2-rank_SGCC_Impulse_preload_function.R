library(ImpulseDE2)
library(ComplexHeatmap)
library(circlize)
library(limma)
library(edgeR)
# this function need to be customized based on the user' object name
prepare_meta_SGCC <- function(spe_obj = NULL,SGCC_df = NULL, 
                              cell_pair = "CD8T_M1",CT.use = "CD8T",condition_name = "EBV_Indicator",
                              sample.drop = c("DFCI_9.1","DFCI_21.1","DFCI_11.1"),batch.factor = "Cohort"){
  # filter out sgcc score = na
  sample.sgcc.na <- colnames(SGCC_df)[is.na(SGCC_df[cell_pair,])]
  # filter out samples
  spe_sub <- spe_obj[, colData(spe_obj)$MergedLabel  == CT.use & 
                   ! colData(spe_obj)$ROI_rename  %in%sample.drop &
                   ! colData(spe_obj)$ROI_rename %in% sample.sgcc.na]
  # input spe_obj object
  pb <- counts(spe_sub)
  norm.pb <- spe_obj@assays@data$logcounts
  SGCC_score = as.numeric(SGCC_df[cell_pair,spe_sub$ROI_rename])
  
  dfAnnotation <- data.frame(row.names = colnames(pb),
                             Sample = colnames(pb),
                             Condition = ifelse(colData(spe_sub)[,condition_name] == "yes", "case","control"),
                             SGCC = SGCC_score,
                             Time = as.numeric(cut(SGCC_score, 
                                                   breaks = quantile(SGCC_score, probs = seq(0, 1, by = 1/3)), 
                                                   labels = 1:3, 
                                                   include.lowest = TRUE)),
                             Batch = colData(spe_sub)[,batch.factor],
                             ROI_CT = paste0(colData(spe_sub)$ROI_rename,"_",colData(spe_sub)$MergedLabel))
  
  dfAnnotation <- dfAnnotation[order(dfAnnotation$SGCC),]
  pb <- pb[,dfAnnotation$Sample]
  return(list(dfAnnotation, as.matrix(pb),as.matrix(norm.pb)))
} 



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
evalImpulse_comp <- compiler::cmpfun(evalImpulse)

rank_DEG <- function(objectImpulseDE2 = objectImpulseDE2, pvalue.cut= 0.05){
  pvalue.cut = pvalue.cut
  dfAnnot <- get_dfAnnotationProc(obj=objectImpulseDE2)
  
  scaNGenes <- dim(get_matCountDataProc(obj=objectImpulseDE2))[1]
  # Order genes by time of extremum (peak/valley)
  vecSignificantIDs <- rownames(objectImpulseDE2$dfImpulseDE2Results[
    !is.na(objectImpulseDE2$dfImpulseDE2Results$padj) & 
      objectImpulseDE2$dfImpulseDE2Results$p < pvalue.cut, ])
  vecTimePointsToEval <- sort(unique(dfAnnot$Time), 
                              decreasing = FALSE)
  scaNTPtoEvaluate <- length(vecTimePointsToEval)
  matImpulseValue <- do.call(rbind, lapply(
    vecSignificantIDs, function(x) {
      evalImpulse_comp(
        vecImpulseParam = 
          get_lsModelFits(obj=objectImpulseDE2)[["case"]][[x]]$
          lsImpulseFit$vecImpulseParam, 
        vecTimepoints = vecTimePointsToEval)
    }))
  rownames(matImpulseValue) <- vecSignificantIDs
  matidxMaxTimeSort <- t(apply(matImpulseValue, 1, function(genevalues) {
    sort(genevalues, decreasing = TRUE, index.return = TRUE)$ix
  }))
  vecMaxTime <- vecTimePointsToEval[matidxMaxTimeSort[, 1]]
  matidxMinTimeSort <- t(apply(matImpulseValue, 1, function(genevalues) {
    sort(genevalues, decreasing = FALSE, index.return = TRUE)$ix
  }))
  
  # Time_DIFF <- rowSums(matidxMinTimeSort[,c(3,4)])-rowSums(matidxMinTimeSort[,c(1,2)])
  # Time_DIFF <- rowSums(matidxMinTimeSort[,c(3,4)])-rowSums(matidxMinTimeSort[,c(1,2)])
  Time_DIFF <- matidxMinTimeSort[,3]-matidxMinTimeSort[,1]
  increase_gene <- rownames(matidxMinTimeSort)[Time_DIFF>0]
  Decrease_gene <- rownames(matidxMinTimeSort)[Time_DIFF<0]
  Impulse_res_out <- data.frame(T1 = matidxMinTimeSort[,1],
                                T2 = matidxMinTimeSort[,2],
                                T3 = matidxMinTimeSort[,3],
                                #  T4 = matidxMinTimeSort[,4],
                                Time_Difference = Time_DIFF,
                                monotone = ifelse( Time_DIFF > 0, "Increase",
                                                   ifelse(Time_DIFF < 0,"Decrease","NA")))
  Impulse_res_out$p <- objectImpulseDE2$dfImpulseDE2Results[rownames(Impulse_res_out),]$p
  return(Impulse_res_out)
  # lsHeatmaps <- plotHeatmap(
  #   objectImpulseDE2       = objectImpulseDE2,
  #   strCondition           = "case",
  #   boolIdentifyTransients = TRUE,
  #   scaQThres              = 0.9)
  # Combine heatmap with side bars
  
}

## Heatmap function
PlotGeneHeatmap <- function(ExpMat = NULL, gene = NULL, SampleAnnotation = NULL,cluster_rows = TRUE,... ){
  SGCC = SampleAnnotation$SGCC
  column_density <- HeatmapAnnotation(
    SGCC = anno_lines(SGCC,which = "column"),
    gp = gpar(lwd = 10),
    height = unit(1.3, "cm")
  )
  interested.gene <- intersect(gene,rownames(ExpMat))
  ExpMat_sub <- ExpMat[interested.gene,SampleAnnotation$Sample]
  scale.ExpMat_sub <- t(scale(t(as.matrix(ExpMat_sub))))
  upper.value <- ceiling(max(c(abs(min(scale.ExpMat_sub)),max(scale.ExpMat_sub))))
  ht_list <- Heatmap(scale.ExpMat_sub,
                     cluster_columns = F,
                     cluster_rows = T, 
                     row_names_gp = gpar(fontsize = 10),top_annotation = column_density,
                     col = colorRamp2(c(-upper.value, 0, upper.value), 
                                      c("#4276b6", "#fcfdfe", "#d93027")))
  return(ht_list)
}
# Load necessary libraries


runDEGs_limma_fix_EBV<- function(meta_data = prepare_meta_output_CT1[[1]], 
                                    expression_matrix = prepare_meta_output_CT1[[2]], 
                                    consider_CTnumber = T, 
                                    cellnumber_percore = cellnumber_percore,
                                    time_column = "Time", 
                                    EBV_column = "Condition",  # Column for EBV status
                                    EBV_status = "case",            # Filter for specific EBV status
                                    time1 = 1, 
                                    time3 = 3, 
                                    logfc_threshold = 0.05, 
                                    p_val_threshold = 0.05,
                                    covariate_ID = c("ruv_W1"),
                                    covariate_mat = covariate_matrix) {
  
  # Step 1: Filter metadata for the specified Time points (Time 1 and Time 3) and EBV status
  meta_data_filtered <- meta_data[meta_data[[time_column]] %in% c(time1, time3) & meta_data[[EBV_column]] %in% EBV_status, ]
  if(length(covariate_ID)==1){
    covariate_df <- data.frame(row.names = rownames(meta_data_filtered),
                               ruv_covariate = covariate_mat[rownames(meta_data_filtered)])
    meta_data_filtered <- cbind(meta_data_filtered,covariate_df)
  }
  if(length(covariate_ID)>1){
    covariate_df <- data.frame(row.names = rownames(meta_data_filtered),
                               covariate_mat[rownames(meta_data_filtered),covariate_ID])
    meta_data_filtered <- cbind(meta_data_filtered,covariate_df)
  }
  if(length(covariate_ID)<1){
    meta_data_filtered <- cbind(meta_data_filtered)
  }
  
  
  # Step 2: Subset the expression matrix to include only the filtered samples
  samples_to_keep <- rownames(meta_data_filtered)
  expression_matrix_filtered <- expression_matrix[, samples_to_keep]
  
  # Step 3: Create the design matrix
  # Convert Time to factor and ensure Time 1 is the baseline level
  meta_data_filtered[[time_column]] <- as.factor(meta_data_filtered[[time_column]])
  meta_data_filtered[[time_column]] <- relevel(meta_data_filtered[[time_column]], ref = as.character(time1))
  rownames(cellnumber_percore) <- cellnumber_percore$ROI_CT
  cellnumber = cellnumber_percore[meta_data_filtered$ROI_CT,]$cell_count
  if( consider_CTnumber == T){
    design_data <- data.frame(Time = meta_data_filtered[[time_column]], 
                              meta_data_filtered[,grep("ruv",colnames(meta_data_filtered))],
                              cellnumber = cellnumber)
  }else{
    design_data <- data.frame(Time = meta_data_filtered[[time_column]], 
                              meta_data_filtered[,grep("ruv",colnames(meta_data_filtered))])
  }
  
  design <- model.matrix(~ Time + ., data = design_data)
  print(design)
  # Step 4: Set up DGEList object for edgeR
  # dge <- DGEList(counts = (exp(expression_matrix_filtered)-1))
  dge <- DGEList(counts = expression_matrix_filtered)
  # Step 5: Normalize the data using voom (limma-voom)
  v <- voom(dge, design, plot = TRUE)
  
  # Step 6: Fit the linear model using limma
  fit <- lmFit(v, design)
  
  # Step 7: Apply empirical Bayes moderation
  fit <- eBayes(fit)
  
  # Step 8: Extract DEGs between Time 1 and Time 3 (Time3 - Time1)
  DEG_result <- topTable(fit, coef = "Time3", number = Inf, adjust.method = "fdr", sort.by = "P")
  DEG_result <- DEG_result %>%
    mutate(positive_in = ifelse(logFC > 0, paste0(time_column,time3), paste0(time_column,time1))) # Positive logFC means upregulated in ident2
  # Step 9: Filter DEGs based on thresholds
  filtered_DEGs <- DEG_result[DEG_result$P.Value < p_val_threshold & abs(DEG_result$logFC) > logfc_threshold, ]

  
  filtered_DEGs <- DEG_result[order(filtered_DEGs$logFC),]
  # Return the filtered DEGs
  return(filtered_DEGs)
}


runDEGs_limma_fix_Time<- function(meta_data = prepare_meta_output_CT1[[1]], 
                                 expression_matrix = prepare_meta_output_CT1[[2]], 
                                 consider_CTnumber = T, 
                                 cellnumber_percore = cellnumber_percore,
                                 time_column = "Time", 
                                 EBV_column = "Condition",  # Column for EBV status
                                 time_fix = 1, 
                                 logfc_threshold = 0.05, 
                                 p_val_threshold = 0.05,
                                 covariate_ID = c("ruv_W1"),
                                 covariate_mat = covariate_matrix) {
  
  # Step 1: Filter metadata for the specified Time points (Time 1 and Time 3) and EBV status
  meta_data_filtered <- meta_data[meta_data[[time_column]] %in% c(time_fix) ,]
  if(length(covariate_ID)==1){
    covariate_df <- data.frame(row.names = rownames(meta_data_filtered),
                               ruv_covariate = covariate_mat[rownames(meta_data_filtered),covariate_ID])
    meta_data_filtered <- cbind(meta_data_filtered,covariate_df)
  }
  if(length(covariate_ID)>1){
    covariate_df <- data.frame(row.names = rownames(meta_data_filtered),
                               covariate_matrix[rownames(meta_data_filtered),covariate_ID])
    meta_data_filtered <- cbind(meta_data_filtered,covariate_df)
  }
  if(length(covariate_ID)<1){
    meta_data_filtered <- cbind(meta_data_filtered)
  }
  
  
  # Step 2: Subset the expression matrix to include only the filtered samples
  samples_to_keep <- rownames(meta_data_filtered)
  expression_matrix_filtered <- expression_matrix[, samples_to_keep]
  
  # Step 3: Create the design matrix
  rownames(cellnumber_percore) <- cellnumber_percore$ROI_CT
  cellnumber = cellnumber_percore[meta_data_filtered$ROI_CT,]$cell_count
  if( consider_CTnumber == T){
    design_data <- data.frame(EBV = meta_data_filtered[[EBV_column]], 
                              meta_data_filtered[,grep("ruv",colnames(meta_data_filtered))],
                              cellnumber = cellnumber)
  }else{
    design_data <- data.frame(EBV = meta_data_filtered[[EBV_column]], 
                              meta_data_filtered[,grep("ruv",colnames(meta_data_filtered))])
  }
  
  design <- model.matrix(~ EBV + ., data = design_data)
  print(design[1:3,])
  # Step 4: Set up DGEList object for edgeR
  dge <- DGEList(counts = expression_matrix_filtered)
  
  # Step 5: Normalize the data using voom (limma-voom)
  v <- voom(dge, design, plot = TRUE)
  
  # Step 6: Fit the linear model using limma
  fit <- lmFit(v, design)
  
  # Step 7: Apply empirical Bayes moderation
  fit <- eBayes(fit)
  
  # Step 8: Extract DEGs 
  DEG_result <- topTable(fit, coef = "EBVcontrol", number = Inf, adjust.method = "fdr", sort.by = "P")
  DEG_result <- DEG_result %>%
    mutate(positive_in = ifelse(logFC > 0, "EBV-", "EBV+")) # Positive logFC means upregulated in ident2
  # Step 9: Filter DEGs based on thresholds
  filtered_DEGs <- DEG_result[DEG_result$P.Value < p_val_threshold & abs(DEG_result$logFC) > logfc_threshold, ]
  
  
  filtered_DEGs <- DEG_result[order(filtered_DEGs$logFC),]
  # Return the filtered DEGs
  return(filtered_DEGs)
}


