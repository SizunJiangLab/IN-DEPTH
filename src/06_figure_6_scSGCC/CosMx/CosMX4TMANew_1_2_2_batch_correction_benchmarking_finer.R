# New CosMX Batch Effect Correction and Benchmarking (Finer Cell Type)
# Script 1_2: Batch effect correction, benchmarking parameter, and seurat validation
# Load pseudobulk data from New_CosMX4TMA_1_1_pseudobulk_creation_finer.R and perform batch correction and benchmarking

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(edgeR)
library(readxl)
library(SpatialExperiment)
library(qs)
library(cowplot)
library(GSVA)
library(patchwork)

# Load modular functions
base_root <- "/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/CosMXNew_code_data_publish/"
standr_path <- file.path(base_root, "standR_covariate")
output_base <- file.path(base_root, "Data/CosMXNew/Exampleoutpu")
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

# Download standR_covariate package to local directory if not present
if (!dir.exists(standr_path)) {
  cat("Downloading standR_covariate from GitHub to local directory...\n")
  system(paste0("git clone https://github.com/BMEngineeR/standR_covariate.git ", standr_path))
}

# Load the package from local directory
devtools::load_all(standr_path)
source(file.path(base_root, "Code/src/6-Pathway_validation_DLBCL_preload_function-sc.R"))
source(file.path(base_root, "Code/src/utils.R"))
source(file.path(base_root, "Code/src/5-Preload_customized_gene_list.R"))


# Set paths
# wdpath <- file.path(base_root, "Data")
# data_path <- file.path(wdpath, "New_CosMX")
input_intermediate_dir <- file.path(output_base, "CosMX4TMANew_1_1_pseudobulk_creation_finer")
results_path <- file.path(output_base, "CosMX4TMANew_1_2_2_batch_correction_benchmarking_finer")
# Create results directory if it doesn't exist
if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
}


# setwd(data_path)

cat("=== Script New_1_2: Batch Effect Correction and Benchmarking (FINER Cell Types) ===\n")

## Load intermediate results from New_CosMX4TMA_1_1_pseudobulk_creation_finer.R
cat("Loading intermediate results from New_CosMX4TMA_1_1_pseudobulk_creation_finer.R...\n")
pseudobulk_expr_filtered <- qread(file.path(input_intermediate_dir, "New_pseudobulk_expr_filtered_finer.qs"))
sample_metadata <- qread(file.path(input_intermediate_dir, "New_sample_metadata_finer.qs"))
cell_type_summary <- qread(file.path(input_intermediate_dir, "New_cell_type_summary_finer.qs"))

cat("Data loaded successfully:\n")
cat("- Expression matrix dimensions:", dim(pseudobulk_expr_filtered), "\n")
cat("- Sample metadata dimensions:", dim(sample_metadata), "\n")

# ============================================================================
# PARAMETER SEARCH WITH standR-covariate
# - Batch factor: TMA (batch across different TMA sources)
# - Biological factors: EBV_status, cell_type (using FINER annotations)
# - Covariate: cell number (if available)
# - Outputs: parameter_search_scores_finer.csv, parameter_search_spe_lists_finer.qs -> results_path
# ============================================================================

cat("\n=== Parameter Search (standR-covariate) with FINER Cell Types ===\n")

# Create SPE object
spe_obj <- SpatialExperiment(
  assays = list(counts = pseudobulk_expr_filtered),
  colData = sample_metadata
)
spe_obj$total_counts <- colSums(spe_obj@assays@data$counts)

# ============================================================================
# REMOVE VIRUS GENES FROM EXPRESSION MATRIX (Optional)
# ============================================================================
cat("\n=== Checking for Virus Genes ===\n")

# Identify virus genes (EBV-related)
virus_genes_cosmx <- c(grep("^LMP", rownames(spe_obj), value = TRUE),
                       grep("^BCRF1", rownames(spe_obj), value = TRUE),
                       grep("^BNLF", rownames(spe_obj), value = TRUE),
                       grep("^BLLF1", rownames(spe_obj), value = TRUE),
                       grep("^BZL", rownames(spe_obj), value = TRUE),
                       grep("^EBNA", rownames(spe_obj), value = TRUE),
                       grep("^RPM", rownames(spe_obj), value = TRUE))

cat("Total genes before filtering:", nrow(spe_obj), "\n")
cat("Virus genes identified:", length(virus_genes_cosmx), "\n")
if (length(virus_genes_cosmx) > 0) {
  cat("Virus genes:", paste(virus_genes_cosmx, collapse = ", "), "\n")

  # Option to remove virus genes (uncomment if needed)
  # non_virus_genes <- setdiff(rownames(spe_obj), virus_genes_cosmx)
  # spe_obj <- spe_obj[non_virus_genes, ]
  # cat("Genes retained after virus gene removal:", nrow(spe_obj), "\n")
} else {
  cat("No virus genes found in the dataset.\n")
}

# ============================================================================
# OPTIONAL: SUBSET TO SELECTED TMAs (CORES) FOR BATCH CORRECTION
# ============================================================================
# Set the TMA(s) of interest here, e.g. c("DFCI","Rochester"); use character(0)
# or leave as NULL to keep all TMAs.
selected_TMA <- c("DFCI", "Rochester")  # Set to NULL to use all TMAs, or specify specific TMAs like c("DFCI", "Rochester")

# For finer cell types, we'll select the most abundant ones across all samples
# Get unique cell types and their counts
finer_cell_type_counts <- table(sample_metadata$cell_type)
cat("\nFiner cell type distribution in pseudobulk samples:\n")
print(sort(finer_cell_type_counts, decreasing = TRUE))

# Select finer cell types with sufficient samples (e.g., at least 10 samples)
min_samples_per_cell_type <- 5
wanted_groups <- names(finer_cell_type_counts[finer_cell_type_counts >= min_samples_per_cell_type])

cat("\nSelected finer cell types with at least", min_samples_per_cell_type, "samples:\n")
print(wanted_groups)

# Subset to desired cell groups if present
subset_spe_obj <- if ("cell_type" %in% colnames(colData(spe_obj))) {
  spe_obj[, spe_obj$cell_type %in% wanted_groups]
} else {
  spe_obj
}

subset_spe_obj$cell_count <- colSums(subset_spe_obj@assays@data$counts)

if (!is.null(selected_TMA) && length(selected_TMA) > 0 && "TMA" %in% colnames(colData(subset_spe_obj))) {
  keep_tma <- subset_spe_obj$TMA %in% selected_TMA
  cat("Subsetting SPE object to selected TMA(s):", paste(selected_TMA, collapse = ", "), "\n")
  cat("Samples before TMA subset:", ncol(subset_spe_obj), "\n")
  subset_spe_obj <- subset_spe_obj[, keep_tma]
  cat("Samples after TMA subset:", ncol(subset_spe_obj), "\n")
} else {
  cat("No TMA subset specified or 'TMA' column not present; using all samples.\n")
}

# Attach/derive cell count covariate if possible
cat("Attaching cell count covariate if available...\n")
cell_count_vec <- NULL

# Try common column names in sample metadata
for (cand in c("cell_count", "cell_number", "n_cells", "num_cells", "cells")) {
  if (cand %in% colnames(colData(subset_spe_obj))) {
    cell_count_vec <- colData(subset_spe_obj)[[cand]]
    break
  } else if (cand %in% colnames(sample_metadata)) {
    cell_count_vec <- sample_metadata[[cand]]
    break
  }
}

# Try deriving from cell_type_summary if needed
if (is.null(cell_count_vec) && is.data.frame(cell_type_summary)) {
  id_col <- intersect(c("sample_name", "sample", "sample_id", "Sample", "SampleID", "name", "window_id"),
                      colnames(cell_type_summary))
  cnt_col <- intersect(c("total_cells", "cell_count", "n_cells", "num_cells", "cells", "n"),
                       colnames(cell_type_summary))

  if (length(id_col) >= 1 && length(cnt_col) >= 1) {
    agg_df <- cell_type_summary %>%
      dplyr::group_by(.data[[id_col[1]]]) %>%
      dplyr::summarise(cell_count = sum(.data[[cnt_col[1]]], na.rm = TRUE), .groups = "drop")

    # Map onto column order of subset_spe_obj using sample_name if present, else colnames
    sample_id_col <- intersect(c("sample_name", "sample", "sample_id", "Sample", "SampleID", "name", "window_id"),
                               colnames(colData(subset_spe_obj)))
    if (length(sample_id_col) >= 1) {
      key_vec <- colData(subset_spe_obj)[[sample_id_col[1]]]
      m <- match(key_vec, agg_df[[id_col[1]]])
      cell_count_vec <- agg_df$cell_count[m]
    }
  }
}

if (!is.null(cell_count_vec)) {
  colData(subset_spe_obj)$cell_count <- as.numeric(cell_count_vec)
  cat("Cell count covariate attached.\n")
} else {
  cat("Cell count covariate not found; will run without covariate where required.\n")
}

# Ensure EBV_status is properly formatted
if (!"EBV_status" %in% colnames(colData(subset_spe_obj)) && "EBV_status.x" %in% colnames(colData(subset_spe_obj))) {
  colData(subset_spe_obj)$EBV_status <- colData(subset_spe_obj)$EBV_status.x
}

# Verify biological and batch labels
bio_labels <- intersect(c("cell_type", "EBV_status"), colnames(colData(subset_spe_obj)))
batch_label <- intersect(c("TMA"), colnames(colData(subset_spe_obj)))

cat("\nBiological labels found:", paste(bio_labels, collapse = ", "), "\n")
cat("Batch labels found:", paste(batch_label, collapse = ", "), "\n")

# Ensure we have at least one batch label
if (length(batch_label) == 0) {
  stop("No batch label (TMA, Source, or runID) found in the data!")
}

# Scoring helpers (ARI)
compute_ari_scores <- function(spe_in, bio_labels, batch_labels, resolution = 0.5) {
  mx <- assay(spe_in, "logcounts")
  mx <- mx + matrix(rnorm(length(mx), mean = 0, sd = 1e-6),
                    nrow = nrow(mx), ncol = ncol(mx))
  seu <- CreateSeuratObject(mx, meta.data = as.data.frame(colData(spe_in)))
  seu[["RNA"]]$data <- seu[["RNA"]]$counts
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 100, assay = "RNA")
  seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE, npcs = 15)
  seu <- FindNeighbors(seu, dims = 1:10, verbose = FALSE)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  clust <- as.factor(Idents(seu))
  md <- as.data.frame(seu@meta.data)

  bio_ari <- setNames(numeric(0), character(0))
  if (length(bio_labels) > 0) {
    for (lb in bio_labels) {
      if (!is.null(md[[lb]])) {
        bio_ari[lb] <- mclust::adjustedRandIndex(clust, as.factor(md[[lb]]))
      }
    }
  }

  batch_ari <- setNames(numeric(0), character(0))
  if (length(batch_labels) > 0) {
    for (lb in batch_labels) {
      if (!is.null(md[[lb]])) {
        batch_ari[lb] <- mclust::adjustedRandIndex(clust, as.factor(md[[lb]]))
      }
    }
  }

  list(bio = bio_ari, batch = batch_ari)
}

score_from_ari <- function(bio_ari, batch_ari) {
  bmean <- if (length(bio_ari)) mean(bio_ari, na.rm = TRUE) else NA_real_
  tmean <- if (length(batch_ari)) mean(batch_ari, na.rm = TRUE) else NA_real_
  overall <- bmean - tmean
  list(bio_mean = bmean, batch_mean = tmean, overall = overall)
}

# Search space
norm_methods <- c("CPM", "TMM")
norm_prefix <- function(m) if (m == "upperquartile") "Upper" else if (m == "CPM") "LogCPM" else "TMM"
covariate_options <- c(TRUE, FALSE)
topn_values <- seq(1000, 3000, by = 1000)
k_values <- 1:3

# Results
spe_lists <- list()
param_scores <- data.frame()
subset_spe_obj <- subset_spe_obj[, subset_spe_obj$cell_type %in% c("B cell","CD4 T", "CD8 T", "Macrophage")]
for (nmeth in norm_methods) {
  cat("Normalization:", nmeth, "\n")
  normalized_spe <- geomxNorm(subset_spe_obj, method = nmeth, log = TRUE)
  prefix <- norm_prefix(nmeth)

  for (top_n in topn_values) {
    for (k in k_values) {
      set.seed(123)

      # Discover NCGs using the specified batch factor
      tmp_spe_bn <- findNCGs(normalized_spe, batch_name = batch_label[1], top_n = top_n)

      for (use_cov in covariate_options) {
        cov_vec <- if (use_cov && "cell_count" %in% colnames(colData(tmp_spe_bn))) {
          colData(tmp_spe_bn)$cell_count
        } else {
          NULL
        }

        if (use_cov && is.null(cov_vec)) {
          # Skip covariate-on run if covariate is unavailable
          next()
        }

        tmp_spe <- geomxBatchCorrection(
          tmp_spe_bn,
          factors = as.vector(bio_labels),
          NCGs = S4Vectors::metadata(tmp_spe_bn)$NCGs,
          k = k,
          covariates = cov_vec,
          isLog = TRUE
        )

        # Tag parameters
        S4Vectors::metadata(tmp_spe)$param_roles <- list(
          method = "RUV4",
          normalization = nmeth,
          biological_factors = bio_labels,
          batch_factors = batch_label,
          covariates = if (!is.null(cov_vec)) "cell_count" else character(0),
          ruv_k = k,
          top_n = top_n
        )

        wf_str <- paste(bio_labels, collapse = "+")
        spe_name <- paste0(prefix, "_RUV4_top", top_n, "_k", k, "_bnset_", batch_label[1],
                          "_cov_", use_cov, "_wf_", wf_str)
        spe_lists[[spe_name]] <- tmp_spe

        # ARI scoring
        ari <- compute_ari_scores(tmp_spe, bio_labels, batch_label)
        sc <- score_from_ari(ari$bio, ari$batch)

        bio_str <- if (length(ari$bio)) paste0(names(ari$bio), "=", sprintf("%.4f", ari$bio), collapse = ";") else ""
        batch_str <- if (length(ari$batch)) paste0(names(ari$batch), "=", sprintf("%.4f", ari$batch), collapse = ";") else ""

        param_scores <- rbind(param_scores, data.frame(
          name = spe_name,
          method = "RUV4",
          normalization = nmeth,
          top_n = top_n,
          k = k,
          batch_set = batch_label[1],
          covariate_on = use_cov,
          wanted = wf_str,
          bio_mean = sc$bio_mean,
          batch_mean = sc$batch_mean,
          overall = sc$overall,
          bio_ari = bio_str,
          batch_ari = batch_str,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Save outputs (tagged by selected TMA if any) with "_finer" suffix
cat("\nSaving parameter search outputs...\n")

if (!is.null(selected_TMA) && length(selected_TMA) > 0) {
  tma_tag <- paste(selected_TMA, collapse = "_")
  scores_path <- file.path(results_path, paste0("New_parameter_search_scores_TMA-", tma_tag, "_finer_subsetTMacroTumor.csv"))
  spe_path    <- file.path(results_path, paste0("New_parameter_search_spe_lists_TMA-", tma_tag, "_finer_subsetTMacroTumor.qs"))
} else {
  scores_path <- file.path(results_path, "New_parameter_search_scores_finer_subsetTMacroTumor.csv")
  spe_path    <- file.path(results_path, "New_parameter_search_spe_lists_finer_subsetTMacroTumor.qs")
}

write.csv(param_scores, scores_path, row.names = FALSE)
qsave(spe_lists, file = spe_path)

cat("Parameter search completed.\n")
cat("Output files:\n")
cat("- Scores:", scores_path, "\n")
cat("- SPE lists:", spe_path, "\n")

# ============================================================================
# TOP-5 SELECTION AND UMAP PLOTS
# ============================================================================
cat("\n=== Selecting Top 5 and Rendering UMAPs ===\n")

# Reload to ensure consistency
param_scores <- read.csv(scores_path)
spe_lists <- qread(spe_path)

# Select top 5 by ARI-based overall score
param_scores <- param_scores[order(-param_scores$overall), ]
top5 <- head(param_scores$name, 5)
top5 <- top5[top5 %in% names(spe_lists)]
top5_spe <- spe_lists[top5]

cat("Top 5 parameter combinations:\n")
print(param_scores[1:min(5, nrow(param_scores)), c("name", "bio_mean", "batch_mean", "overall")])

# Helper function to parse parameter name
parse_param_name <- function(nm) {
  parts <- strsplit(nm, "_")[[1]]
  norm <- parts[1]
  topN <- sub("^top", "", parts[grep("^top[0-9]+$", parts, perl = TRUE)][1])
  klev <- sub("^k", "", parts[grep("^k[0-9]+$", parts, perl = TRUE)][1])
  cov_flag <- {
    i <- which(parts == "cov")
    if (length(i) == 1 && i < length(parts)) parts[i + 1] else NA_character_
  }
  wf <- {
    i <- which(parts == "wf")
    if (length(i) == 1 && i < length(parts)) paste(parts[(i + 1):length(parts)], collapse = "_") else NA_character_
  }
  list(normalization = norm, top_n = topN, k = klev, covariate = cov_flag, wanted = wf)
}

# UMAP plotting function
plot_umap_grid <- function(spe_in, spe_name = NULL, bio_cols = bio_labels, batch_cols = batch_label) {
  mx <- assay(spe_in, "logcounts")

  ## EBV GSVA score per sample (CosMX EBV gene set)
  virus_genes <- c(
    grep("^LMP",   rownames(mx), value = TRUE),
    grep("^BCRF1", rownames(mx), value = TRUE),
    grep("^BNLF",  rownames(mx), value = TRUE),
    grep("^BLLF1", rownames(mx), value = TRUE),
    grep("^BZL",   rownames(mx), value = TRUE),
    grep("^EBNA",  rownames(mx), value = TRUE),
    grep("^RPM",   rownames(mx), value = TRUE)
  )

  ebv_score <- rep(NA_real_, ncol(mx))
  if (length(virus_genes) > 0) {
    ebv_gene_sets <- list(EBV_genes = virus_genes)
    SetGSVAPar_function <- gsvaParam(
      exprData = as.matrix(mx),
      geneSets = ebv_gene_sets,
      minSize  = 1,
      maxSize  = max(1L, length(virus_genes)),
      kcdf     = "Gaussian"
    )
    ebv_gsva <- gsva(SetGSVAPar_function, verbose = TRUE)
    ebv_score <- as.numeric(ebv_gsva["EBV_genes", ])
  }

  seu <- CreateSeuratObject(mx, meta.data = as.data.frame(colData(spe_in)))
  seu$EBV_GSVA <- ebv_score
  seu[["RNA"]]$data <- seu[["RNA"]]$counts
  DefaultAssay(seu) <- "RNA"
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 200, assay = "RNA")
  seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE, npcs = 14)
  seu <- RunUMAP(seu, dims = 1:10, verbose = FALSE, n.neighbors = 10)

  emb <- as.data.frame(Seurat::Embeddings(seu, "umap"))
  emb$cell_id <- rownames(emb)
  md <- seu@meta.data
  md$cell_id <- rownames(md)
  merged <- dplyr::left_join(emb, md, by = "cell_id")

  cols_requested <- c(bio_cols, batch_cols)
  cols_requested <- cols_requested[!is.null(cols_requested)]
  cols_present <- cols_requested[cols_requested %in% colnames(merged)]

  if (length(cols_present) == 0) stop("No valid factors found to plot.")

  # Keep only factors that have at least one non-NA value
  cols_valid <- vapply(cols_present, function(fn) any(!is.na(merged[[fn]])), logical(1))
  cols_to_plot <- cols_present[cols_valid]

  # Also add continuous covariates if available
  cont_candidates <- c("total_counts", "cell_count", "EBV_GSVA")
  cont_present <- cont_candidates[cont_candidates %in% colnames(merged)]
  cols_to_plot <- unique(c(cols_to_plot, cont_present))

  if (length(cols_to_plot) == 0) stop("All selected factors have only NA values; nothing to plot.")

  # Robust UMAP column selection
  umap_cols <- colnames(emb)
  if (length(umap_cols) < 2) stop("UMAP embedding not found.")
  xcol <- umap_cols[1]; ycol <- umap_cols[2]

  # Build per-factor plots
  mk_plot <- function(fn) {
    fv <- merged[[fn]]
    keep <- !is.na(fv)
    df <- data.frame(
      UMAP_1 = merged[[xcol]][keep],
      UMAP_2 = merged[[ycol]][keep],
      value  = fv[keep],
      stringsAsFactors = FALSE
    )

    if (is.numeric(df$value)) {
      # Continuous visualization
      ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = log1p(value))) +
        geom_point(size = 1, alpha = 0.9) +
        labs(title = paste0("UMAP: log(", fn, ")"), color = fn) +
        scale_color_viridis_c(option = "plasma") +
        theme_classic() +
        theme(legend.position = "bottom")
    } else {
      # Discrete factor visualization
      ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = as.factor(value))) +
        geom_point(size = 1, alpha = 0.9) +
        labs(title = paste0("UMAP: ", fn), color = fn) +
        theme_classic() +
        theme(legend.position = "bottom")
    }
  }

  p_list <- lapply(cols_to_plot, mk_plot)

  # Compose title safely
  title_str <- {
    if (!is.null(spe_name)) {
      info <- parse_param_name(spe_name)
      paste0(
        "Set: ", spe_name, " | Norm=", info$normalization,
        " | top=", info$top_n, " | k=", info$k,
        " | cov=", info$covariate
      )
    } else {
      "UMAP facets"
    }
  }

  patchwork::wrap_plots(plotlist = p_list, nrow = 2) +
    patchwork::plot_annotation(title = title_str)
}

# Plot raw (no batch correction) UMAP
plot_umap_grid_raw <- function(spe_raw, bio_cols = bio_labels, batch_cols = batch_label, normalization = "CPM") {
  cat("Rendering raw (no batch correction) UMAP with normalization:", normalization, "\n")
  normalized <- geomxNorm(spe_raw, method = normalization, log = TRUE)
  plot_umap_grid(normalized, spe_name = paste0("RAW_", normalization), bio_cols = bio_cols, batch_cols = batch_cols)
}

cat("\nPlotting raw (uncorrected) data UMAP...\n")
plot_umap_grid_raw(subset_spe_obj)

if (length(top5_spe) > 0) {
  cat("\nPlotting top batch-corrected parameter set UMAP...\n")
  plot_umap_grid(top5_spe[[1]])
}

cat("\n=== Script completed successfully ===\n")
cat("NOTE: This batch correction was performed using FINER cell type annotations\n")
cat("Results saved to:", results_path, "\n")
