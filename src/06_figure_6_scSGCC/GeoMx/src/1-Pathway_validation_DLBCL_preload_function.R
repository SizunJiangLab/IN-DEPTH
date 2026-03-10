# Load necessary libraries
library(Seurat)
library(dplyr)
library(limma)
library(edgeR)
library(progress)
perform_DEA_timepoint <- function(expr_data, meta_data, timepoint_var = "Time", condition_var = "Condition", batch_var = "Batch", time1 = 1, time3 = 3) {
  # Ensure Time and Condition are factors
  meta_data[[timepoint_var]] <- factor(meta_data[[timepoint_var]], levels = c(time1, time3))  # Time 1 is the reference
  meta_data[[condition_var]] <- factor(meta_data[[condition_var]])

  # Function to perform limma analysis for a specific condition (case or control)
  run_DEA <- function(condition_value) {
    # Subset metadata and expression data for the specific condition
    meta_subset <- meta_data[meta_data[[condition_var]] == condition_value & meta_data[[timepoint_var]] %in% c(time1, time3), ]
    expr_subset <- expr_data[, rownames(meta_subset)]
    meta_subset$Time <- as.factor(meta_subset$Time)
    dge <- DGEList(expr_subset,group = meta_subset$Time)
    # Create the design matrix
    design <- model.matrix(~ 0+ Time,data = meta_subset)
    contr.matrix <- makeContrasts(
      Highvlow = Time3 - Time1,
      levels = colnames(design))
    dge_all <- estimateDisp(dge, design = design, robust = TRUE)
    bcv_df <- data.frame(
      'BCV' = sqrt(dge_all$tagwise.dispersion),
      'AveLogCPM' = dge_all$AveLogCPM,
      'gene_id' = rownames(dge_all)
    )

    # Fit the linear model using limma
    fit <- lmFit(expr_subset, design)
    fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
    # Apply empirical Bayes moderation to the linear model
    fit <- eBayes(fit)
    results_efit<- decideTests(fit, p.value = 0.05)
    summary_efit <- summary(results_efit)
    # Extract the differential expression results for Time 3 vs Time 1
    results <- topTable(fit, coef = 1,number = Inf,sort.by = "P")

    return(results)
  }

  # Perform DEA for case group
  dea_case <- run_DEA("case")

  # Perform DEA for control group
  dea_control <- run_DEA("control")

  return(list(case_results = dea_case, control_results = dea_control))
}


perform_DEA_limma_edgeR <- function(seurat_obj = my.seurat,
                                   CT_number = NULL,
                                   consider_CTnumber = TRUE,
                                   group_column = "MergedLabel",
                                   EBV_column = "EBV_Indicator",
                                   ident1 = "M1",
                                   ident2 = "Endothelial",
                                   EBV_status = "yes",
                                   confounder_columns = c("ruv_W1", "ruv_W2", "ruv_W3"),
                                   logfc_threshold = 0.05,
                                   p_val_threshold = 0.01,
                                   assay = "data") {

  # Use counts layer (normalized counts are commonly stored here in these pipelines)
  expr_matrix <- GetAssayData(object = seurat_obj, layer = "counts")

  # Extract metadata
  metadata <- seurat_obj@meta.data

  # Optionally merge cell counts if provided via CT_number (expects ROI_CT and cell_count)
  if (!is.null(CT_number) && all(c("ROI_rename", group_column) %in% colnames(metadata))) {
    suppressWarnings({
      metadata <- metadata %>%
        mutate(ROI_CT = paste0(ROI_rename, "_", .data[[group_column]])) %>%
        left_join(CT_number, by = "ROI_CT")
    })
  }

  # Subset by EBV status and build grouping (support target vs REST)
  group <- metadata[[group_column]]
  EBV <- metadata[[EBV_column]]

  if (ident2 == "REST") {
    # Keep only EBV-matching samples; group target vs REST
    cells_to_keep <- EBV %in% EBV_status
    expr_matrix <- expr_matrix[, cells_to_keep]
    group <- group[cells_to_keep]
    metadata <- metadata[cells_to_keep, ]

    group_binary <- ifelse(group == ident1, ident1, "REST")
    group <- factor(group_binary, levels = c("REST", ident1))  # REST as baseline
  } else {
    # Pairwise comparison
    cells_to_keep <- (group %in% c(ident1, ident2)) & (EBV %in% EBV_status)
    expr_matrix <- expr_matrix[, cells_to_keep]
    group <- group[cells_to_keep]
    metadata <- metadata[cells_to_keep, ]
    group <- factor(group, levels = c(ident1, ident2))
  }

  # Build design matrix with optional confounders and CT count
  if (!is.null(confounder_columns)) {
    confounders <- metadata[, confounder_columns, drop = FALSE]
    if (consider_CTnumber) {
      design_data <- data.frame(group = as.factor(group), confounders, CT_count = metadata$cell_count)
    } else {
      design_data <- data.frame(group = as.factor(group), confounders)
    }
    design <- model.matrix(~ group + ., data = design_data)
  } else {
    if (consider_CTnumber) {
      design_data <- data.frame(group = as.factor(group), CT_count = metadata$cell_count)
    } else {
      design_data <- data.frame(group = as.factor(group))
    }
    design <- model.matrix(~ group + ., data = design_data)
  }

  # Fit model depending on data scale
  if (tolower(assay) == "counts") {
    dge <- DGEList(counts = expr_matrix, group = group)
    v <- voom(dge, design, plot = TRUE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
  } else {
    fit <- lmFit(expr_matrix, design)
    fit <- eBayes(fit, trend = TRUE)
  }

  # Coefficient selection by name and labeling
  baseline_level <- levels(group)[1]
  contrast_level <- levels(group)[2]
  group_coef_name <- grep("^group", colnames(design), value = TRUE)[1]
  DEG_result <- topTable(fit, coef = group_coef_name, number = Inf, adjust.method = "fdr", sort.by = "P")

  DEG_result <- DEG_result %>%
    mutate(positive_in = ifelse(logFC > 0, contrast_level, baseline_level))

  filtered_DEGs <- DEG_result %>%
    filter(P.Value < p_val_threshold & abs(logFC) > logfc_threshold) %>%
    arrange(desc(logFC))

  return(filtered_DEGs)
}

hypergeo_test <- function(predicted_list = cd4t_markers,
                          groundtruthlist = genes_list){
  # Initialize a list to store hypergeometric p-values
  hypergeo_pvals <- list()

  # Iterate over each element in cd4t_markers
  for (element_name in names(predicted_list)) {
    # Extract the DEG list for the current element
    deg_list <- predicted_list[[element_name]]

    # Get the DEGs in the current element
    deg_genes <- rownames(deg_list)

    # Perform hypergeometric test for the overlap of CD4T markers and DEGs
    overlap_genes <- intersect(groundtruthlist, deg_genes)  # Overlap between CD4T markers and DEGs
    k <- length(overlap_genes)  # Number of overlapping genes
    m <- length(CD4T_markers)   # Number of CD4T markers
    n <- 18000 - m  # Total genes not in CD4T_markers
    N <- length(deg_genes)  # Number of genes in the DEG list

    # Hypergeometric test (phyper is used for cumulative probability)
    p_value <- phyper(k - 1, m, n, N, lower.tail = FALSE)

    # Store the p-value
    hypergeo_pvals[[element_name]] <- list(pval = p_value, intersect_name = overlap_genes)
  }
  # Create a dataframe to store the results
  hypergeo_results_df <- data.frame(
    parameter = character(),  # Stores the parameter name (element_name)
    p_value = numeric(),      # Stores the hypergeometric p-value
    intersect_genes = character(),  # Stores the intersected genes as a string
    stringsAsFactors = FALSE
  )

  # Iterate over the hypergeo_pvals list to extract p-value and intersected genes
  for (element_name in names(hypergeo_pvals)) {
    # Extract p-value and intersected genes
    p_value <- hypergeo_pvals[[element_name]]$pval
    intersect_genes <- hypergeo_pvals[[element_name]]$intersect_name

    # Convert the intersected genes to a comma-separated string
    intersect_genes_str <- paste(intersect_genes, collapse = ", ")

    # Append a new row to the dataframe
    hypergeo_results_df <- rbind(hypergeo_results_df, data.frame(
      parameter = element_name,
      p_value = p_value,
      intersect_genes = intersect_genes_str,
      stringsAsFactors = FALSE
    ))
  }
  hypergeo_results_df <- hypergeo_results_df[order(hypergeo_results_df$p_value),]
  return(hypergeo_results_df)
}

# Define the benchmarking function for DEG analysis with progress bar
benchmark_DEA_pipeline <- function(seurat_obj,
                                   comparisons,               # List of ident1 vs ident2 comparisons
                                   CT_number= cellnumber_percore,
                                   consider_CTnumberList = c(TRUE, FALSE),
                                   group_column = "MergedLabel",
                                   EBV_column = "EBV_Indicator",
                                   EBV_status = list(c("yes"),
                                                     c("no"),
                                                     c("yes", "no")),  # Handle yes/no options
                                   confounder_combinations = list(NULL,
                                                                  c("ruv_W1"),
                                                                  c("ruv_W1", "ruv_W2"),
                                                                  c("ruv_W1", "ruv_W2", "ruv_W3")),
                                   logfc_threshold = 0.05,
                                   p_val_threshold = 0.001,
                                   assay = "data") {

  # List to store the results
  results_list <- list()


  # Calculate the total number of iterations for the progress bar
  total_iterations <- length(EBV_status) * length(comparisons) * length(confounder_combinations)*length(consider_CTnumberList)

  # Initialize the progress bar
  pb <- progress_bar$new(
    format = "  Benchmarking [:bar] :percent ETA: :eta",
    total = total_iterations, clear = FALSE, width = 60
  )

  # Loop over each EBV status
  for (ebv in EBV_status) {
    # Loop over each comparison (ident1 vs ident2)
    for (comparison in comparisons) {
      ident1 <- comparison[1]
      ident2 <- comparison[2]

      # Loop over each confounder combination
      for (confounder_set in confounder_combinations) {

        # Check if confounders exist in metadata or proceed if NULL
        if (!is.null(confounder_set) && any(!confounder_set %in% colnames(seurat_obj@meta.data))) {
          next()  # Skip invalid confounder sets
        }
        for (consider_CTnumber in consider_CTnumberList ){
          # Perform DEA with the current confounder combination, EBV status, and cell type comparison
          degs <- perform_DEA_limma_edgeR(
            seurat_obj = seurat_obj,
            group_column = group_column,
            CT_number= cellnumber_percore,
            consider_CTnumber = consider_CTnumber,
            EBV_column = EBV_column,
            ident1 = ident1,
            ident2 = ident2,
            EBV_status = ebv,  # Pass the current EBV status
            confounder_columns = confounder_set,
            logfc_threshold = logfc_threshold,
            p_val_threshold = p_val_threshold,
            assay = assay
          )

          # Store the result in the list, including all parameters
          comparison_name <- paste0(ident1, "_vs_", ident2)
          confounder_name <- if (is.null(confounder_set)) "No_Confounders" else paste(confounder_set, collapse = "_")
          ebv_name <- paste0("EBV_", paste(ebv,collapse = "-",sep = ""))
          result_key <- paste(comparison_name, confounder_name, ebv_name, "CTnumber",consider_CTnumber,sep = "_")
          results_list[[result_key]] <- list(
            "DEG_list" = degs,
            "comparison" = comparison,
            "confounders" = confounder_set,
            "EBV_status" = paste0(ebv,"_"),
            "considerCTnumber" = consider_CTnumber
          )

          # Update the progress bar
          pb$tick()
        }
      }
    }
  }

  # Return the benchmark results
  return(results_list)
}
