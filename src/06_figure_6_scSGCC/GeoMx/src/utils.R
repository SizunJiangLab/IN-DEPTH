#' Create Sliding Windows for Spatial Analysis
#'
#' This function creates sliding windows across all FOVs in the metadata for spatial analysis.
#' It divides each FOV into overlapping or non-overlapping windows based on the specified
#' window size and stride parameters.
#'
#' @param mymeta Data frame containing cell metadata with columns: fov, X_cent, Y_cent, cell_id
#' @param window_size Numeric, spatial units for window size (adjustable)
#' @param stride Numeric, step size for sliding window (adjustable)
#' @param min_cells_per_window Numeric, minimum number of cells required for a window to be retained (default: 1)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return List containing:
#'   - mymeta_windowed: Updated metadata with window_id assignments
#'   - all_windows: Data frame with window information (window_id, fov_id, x_start, y_start, x_end, y_end, cell_count)
#'
#' @examples
#' result <- CreateSlidingWindow(mymeta, window_size = 150, stride = 100, verbose = TRUE)
#' mymeta_windowed <- result$mymeta_windowed
#' all_windows <- result$all_windows
#'
#' Create Sliding Windows for Spatial Analysis
#'
#' This function creates sliding windows across all FOVs in the metadata for spatial analysis.
#' It treats the spatial coordinate space as an image, creating a grid of windows based on
#' the overall spatial boundaries (xmin, ymin, xmax, ymax), then assigns cells to windows.
#'
#' @param mymeta Data frame containing cell metadata with columns: fov, X_cent, Y_cent, cell_id
#' @param window_size Numeric, spatial units for window size (adjustable)
#' @param stride Numeric, step size for sliding window (adjustable)
#' @param min_cells_per_window Numeric, minimum number of cells required for a window to be retained (default: 1)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return List containing:
#'   - mymeta_windowed: Updated metadata with window_id assignments
#'   - all_windows: Data frame with window information (window_id, fov_id, x_start, y_start, x_end, y_end, cell_count)
#'
#' @examples
#' result <- CreateSlidingWindow(mymeta, window_size = 150, stride = 100, verbose = TRUE)
#' mymeta_windowed <- result$mymeta_windowed
#' all_windows <- result$all_windows
#'
CreateSlidingWindow <- function(mymeta, window_size = 150, stride = 100, min_cells_per_window = 1, verbose = TRUE) {

  # Input validation
  if (!is.data.frame(mymeta)) {
    stop("mymeta must be a data frame")
  }

  required_cols <- c("fov", "X_cent", "Y_cent", "cell_id")
  missing_cols <- setdiff(required_cols, colnames(mymeta))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  if (window_size <= 0 || stride <= 0) {
    stop("window_size and stride must be positive numbers")
  }
  
  if (min_cells_per_window < 1) {
    stop("min_cells_per_window must be at least 1")
  }

  if (verbose) {
    cat("\n=== Creating Sliding Windows for Spatial Analysis ===\n")
    cat("Window parameters:\n")
    cat("- Window size:", window_size, "spatial units\n")
    cat("- Stride:", stride, "spatial units\n")
    cat("Available FOVs:\n")
    print(table(mymeta$fov))
  }

  # Initialize results
  mymeta$window_id <- NA
  all_windows <- data.frame()

  # Process each FOV separately
  for (fov_id in unique(mymeta$fov)) {
    fov_mask <- mymeta$fov == fov_id
    fov_data <- mymeta[fov_mask, ]

    # Step 1: Get spatial boundaries (treat as image)
    xmin <- min(fov_data$X_cent)
    ymin <- min(fov_data$Y_cent)
    xmax <- max(fov_data$X_cent)
    ymax <- max(fov_data$Y_cent)

    if (verbose) {
      cat("\nFOV", fov_id, "image borders - X: [", round(xmin, 2), ",", round(xmax, 2), "], Y: [", round(ymin, 2), ",", round(ymax, 2), "]\n")
      cat("  Image dimensions:", round(xmax - xmin, 2), "x", round(ymax - ymin, 2), "\n")
    }

    # Step 2: Create sliding window grid on the image space
    # Handle cases where FOV dimension is smaller than window_size
    if (xmax - xmin >= window_size) {
      x_starts <- seq(xmin, xmax - window_size, by = stride)
    } else {
      # FOV width is smaller than window_size - create one window covering entire width
      x_starts <- xmin
      if (verbose) {
        cat("  Note: FOV width (", round(xmax - xmin, 2), ") < window_size (", window_size, ") - using single window for X\n")
      }
    }
    
    if (ymax - ymin >= window_size) {
      y_starts <- seq(ymin, ymax - window_size, by = stride)
    } else {
      # FOV height is smaller than window_size - create one window covering entire height
      y_starts <- ymin
      if (verbose) {
        cat("  Note: FOV height (", round(ymax - ymin, 2), ") < window_size (", window_size, ") - using single window for Y\n")
      }
    }

    if (verbose) {
      cat("  Grid dimensions:", length(x_starts), "x", length(y_starts), "windows\n")
    }

    window_count <- 0

    # Step 3: Create all windows and assign cells
    for (x_start in x_starts) {
      for (y_start in y_starts) {
        window_count <- window_count + 1

        # Define window boundaries
        x_end <- min(x_start + window_size, xmax)
        y_end <- min(y_start + window_size, ymax)

        # Create window ID
        window_id <- paste0("FOV", fov_id, "_window_", window_count)

        # Step 4: Assign cells to this window
        cells_in_window <- fov_data$X_cent >= x_start & fov_data$X_cent <= x_end &
                          fov_data$Y_cent >= y_start & fov_data$Y_cent <= y_end

        cell_count_in_window <- sum(cells_in_window)

        # Only keep windows with sufficient cells
        if (cell_count_in_window >= min_cells_per_window) {
          # Assign window ID to cells
          mymeta$window_id[fov_mask][cells_in_window] <- window_id

          # Store window information
          window_info <- data.frame(
            window_id = window_id,
            fov_id = fov_id,
            x_start = x_start,
            y_start = y_start,
            x_end = x_end,
            y_end = y_end,
            cell_count = cell_count_in_window,
            stringsAsFactors = FALSE
          )
          all_windows <- rbind(all_windows, window_info)
        } else if (verbose) {
          cat("  Skipping window", window_id, "- only", cell_count_in_window, "cell(s)\n")
        }
      }
    }

    if (verbose) {
      cat("FOV", fov_id, "- Created", nrow(all_windows[all_windows$fov_id == fov_id, ]), "valid windows\n")
    }
  }

  # Step 5: Return windowed metadata (only cells assigned to windows)
  mymeta_windowed <- mymeta[!is.na(mymeta$window_id), ]

  if (verbose) {
    cat("\n=== Summary ===\n")
    cat("Total cells after window assignment:", nrow(mymeta_windowed), "\n")
    cat("Total unique windows created:", nrow(all_windows), "\n")
    cat("Windows per FOV:\n")
    print(table(all_windows$fov_id))
    cat("Average cells per window:", round(mean(all_windows$cell_count), 2), "\n")
  }

  # Return results as a list
  return(list(
    mymeta_windowed = mymeta_windowed,
    all_windows = all_windows
  ))
}

#' Create Single Cell Signal Data for SGCC Analysis
#'
#' This function creates spatial signal data for a specific window and cell type using
#' pre-calculated window boundaries. It bins cells into a spatial grid and calculates
#' cell proportions for SGCC analysis.
#'
#' @param meta_data Data frame containing cell metadata with columns: window_id, Annotation_pooled, X_cent, Y_cent
#' @param window_id Character, specific window identifier to process
#' @param cell_type Character, cell type annotation to filter for
#' @param window_info_df Data frame containing window boundary information with columns: window_id, x_start, y_start, x_end, y_end
#' @param grid_size Numeric, size of the spatial grid (default: 10, creates 10x10 grid)
#'
#' @return Numeric vector of length grid_size^2 containing cell proportions for each grid bin
#'
#' @examples
#' signal_data <- CreateSCSignalData(mymeta, "FOV3_window_1", "Macrophage", all_windows, grid_size = 10)
#'
CreateSCSignalData <- function(meta_data, window_id, cell_type, window_info_df, grid_size = 10) {

  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }

  # Import pipe operator
  `%>%` <- dplyr::`%>%`

  # Input validation
  if (!is.data.frame(meta_data)) {
    stop("meta_data must be a data frame")
  }

  if (!is.data.frame(window_info_df)) {
    stop("window_info_df must be a data frame")
  }

  required_meta_cols <- c("window_id", "Annotation_pooled", "X_cent", "Y_cent")
  missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_data))
  if (length(missing_meta_cols) > 0) {
    stop(paste("Missing required columns in meta_data:", paste(missing_meta_cols, collapse = ", ")))
  }

  required_window_cols <- c("window_id", "x_start", "y_start", "x_end", "y_end")
  missing_window_cols <- setdiff(required_window_cols, colnames(window_info_df))
  if (length(missing_window_cols) > 0) {
    stop(paste("Missing required columns in window_info_df:", paste(missing_window_cols, collapse = ", ")))
  }

  if (grid_size <= 0) {
    stop("grid_size must be a positive integer")
  }

  # Create base grid for consistent output structure
  base_grid <- expand.grid(x = 1:grid_size, y = 1:grid_size)

  # Filter data for specific window and cell type
  window_data <- meta_data[meta_data$window_id == window_id &
                             meta_data$Annotation_pooled == cell_type, ]

  if (nrow(window_data) == 0) {
    return(rep(0, grid_size^2))
  }

  # Get pre-calculated spatial boundaries for this window from window_info_df
  window_info <- window_info_df[window_info_df$window_id == window_id, ]

  if (nrow(window_info) == 0) {
    # Fallback to original method if window info not found
    warning(paste("Window info not found for", window_id, "- using fallback method"))
    window_x_range <- range(meta_data$X_cent[meta_data$window_id == window_id])
    window_y_range <- range(meta_data$Y_cent[meta_data$window_id == window_id])
  } else {
    # Use pre-calculated boundaries from window_info_df
    window_x_range <- c(window_info$x_start[1], window_info$x_end[1])
    window_y_range <- c(window_info$y_start[1], window_info$y_end[1])
  }

  # Create bins within this window using the exact boundaries
  x_bin_breaks <- seq(window_x_range[1], window_x_range[2], length.out = grid_size + 1)
  y_bin_breaks <- seq(window_y_range[1], window_y_range[2], length.out = grid_size + 1)

  # Assign cells to bins
  window_data$x_bin <- cut(window_data$X_cent, breaks = x_bin_breaks,
                           labels = 1:grid_size, include.lowest = TRUE)
  window_data$y_bin <- cut(window_data$Y_cent, breaks = y_bin_breaks,
                           labels = 1:grid_size, include.lowest = TRUE)

  # Count cells per bin
  # Use explicit variable references to avoid R CMD check warnings
  x_bin <- y_bin <- cell_proportion <- signal_value <- NULL

  bin_counts <- window_data %>%
    dplyr::filter(!is.na(x_bin) & !is.na(y_bin)) %>%
    dplyr::group_by(x_bin, y_bin) %>%
    dplyr::summarise(cell_count = dplyr::n(), .groups = 'drop') %>%
    dplyr::mutate(x_bin = as.numeric(as.character(x_bin)),
                  y_bin = as.numeric(as.character(y_bin)))

  # Calculate cell proportion (relative to total cells of this type in this window)
  total_cells <- nrow(window_data)
  if (total_cells > 0) {
    bin_counts$cell_proportion <- bin_counts$cell_count / total_cells
  } else {
    bin_counts$cell_proportion <- 0
  }

  # Create full grid and merge with bin counts
  full_signal <- base_grid %>%
    dplyr::left_join(bin_counts, by = c("x" = "x_bin", "y" = "y_bin")) %>%
    dplyr::mutate(signal_value = ifelse(is.na(cell_proportion), 0, cell_proportion)) %>%
    dplyr::select(signal_value)

  return(full_signal$signal_value)
}
#' Create Pseudobulk Expression Data
#'
#' This function creates pseudobulk expression data by aggregating single-cell expression
#' profiles for all cell phenotypes within spatial windows. It can either sum expression across
#' cells or calculate average expression per cell within each window.
#'
#' @param mymeta Data frame containing cell metadata with columns: cell_id, window_id, Annotation_pooled
#' @param mydata Data frame or matrix containing expression data with cells as rows and genes as columns
#' @param cell_types_for_pseudobulk Character vector of cell types to include in pseudobulk creation (default: NULL, uses all available cell types)
#' @param min_cells_per_sample Numeric, minimum number of cells required to create a pseudobulk sample (default: 1)
#' @param normalize_by_cell_count Logical, whether to normalize expression by cell count to get average per cell (default: FALSE, returns sum)
#' @param cell_factor Numeric, if provided, scales expression to what it would be with this many cells (default: NULL, no scaling)
#' @param aggregate_meta_cols Character vector, names of metadata columns to aggregate (calculate mean) for each pseudobulk sample (default: NULL)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return List containing:
#'   - pseudobulk_expr: Matrix with genes as rows and pseudobulk samples as columns
#'   - sample_metadata: Data frame with metadata for each pseudobulk sample (includes aggregated metadata if specified)
#'   - cell_type_summary: Summary of cell counts per cell type and window
#'
#' @examples
#' # Get total expression (sum across cells)
#' result <- CreatePseudoBulk(mymeta, mydata, normalize_by_cell_count = FALSE)
#' 
#' # Get average expression per cell
#' result <- CreatePseudoBulk(mymeta, mydata, normalize_by_cell_count = TRUE)
#' 
#' # Scale expression to what it would be with 1000 cells
#' result <- CreatePseudoBulk(mymeta, mydata, cell_factor = 1000)
#' 
#' # Aggregate specific metadata columns (calculate mean per pseudobulk sample)
#' result <- CreatePseudoBulk(mymeta, mydata, aggregate_meta_cols = c("SGCC", "score", "intensity"))
#'
CreatePseudoBulk <- function(mymeta, mydata, cell_types_for_pseudobulk = NULL, min_cells_per_sample = 1, normalize_by_cell_count = FALSE, cell_factor = NULL, aggregate_meta_cols = NULL, verbose = TRUE) {

  # Input validation
  if (!is.data.frame(mymeta)) {
    stop("mymeta must be a data frame")
  }

  if (!is.data.frame(mydata) && !is.matrix(mydata)) {
    stop("mydata must be a data frame or matrix")
  }

  required_meta_cols <- c("cell_id", "window_id", "Annotation_pooled")
  missing_meta_cols <- setdiff(required_meta_cols, colnames(mymeta))
  if (length(missing_meta_cols) > 0) {
    stop(paste("Missing required columns in mymeta:", paste(missing_meta_cols, collapse = ", ")))
  }

  if (min_cells_per_sample < 1) {
    stop("min_cells_per_sample must be at least 1")
  }
  
  if (!is.null(cell_factor) && (cell_factor <= 0)) {
    stop("cell_factor must be a positive number")
  }
  
  # Validate aggregate_meta_cols
  if (!is.null(aggregate_meta_cols)) {
    if (!is.character(aggregate_meta_cols)) {
      stop("aggregate_meta_cols must be a character vector")
    }
    missing_cols <- setdiff(aggregate_meta_cols, colnames(mymeta))
    if (length(missing_cols) > 0) {
      stop(paste("aggregate_meta_cols contains columns not found in mymeta:", paste(missing_cols, collapse = ", ")))
    }
    # Check if columns are numeric
    non_numeric_cols <- aggregate_meta_cols[!sapply(aggregate_meta_cols, function(col) is.numeric(mymeta[[col]]))]
    if (length(non_numeric_cols) > 0) {
      warning(paste("Non-numeric columns in aggregate_meta_cols will be converted to numeric:", paste(non_numeric_cols, collapse = ", ")))
    }
  }

  if (verbose) {
    cat("\n=== Creating Pseudobulk Expression Data ===\n")
    cat("Input data dimensions:\n")
    cat("- Metadata:", dim(mymeta), "\n")
    cat("- Expression data:", dim(mydata), "\n")
    
    if (!is.null(cell_factor)) {
      cat("- Expression scaling: Normalized to", cell_factor, "cells\n")
    } else if (normalize_by_cell_count) {
      cat("- Normalization method: Average per cell\n")
    } else {
      cat("- Normalization method: Total sum\n")
    }
    
    if (!is.null(aggregate_meta_cols)) {
      cat("- Metadata aggregation: Calculating mean for", length(aggregate_meta_cols), "columns:", paste(aggregate_meta_cols, collapse = ", "), "\n")
    }
  }

  # Match cell IDs between metadata and expression data
  common_cells <- intersect(rownames(mymeta), rownames(mydata))
  if (length(common_cells) == 0) {
    # Try using cell_id column if rownames don't match
    if ("cell_id" %in% colnames(mymeta)) {
      common_cells <- intersect(mymeta$cell_id, rownames(mydata))
      if (length(common_cells) > 0) {
        # Set rownames for mymeta to match
        rownames(mymeta) <- mymeta$cell_id
        mymeta_matched <- mymeta[common_cells, ]
        mydata_matched <- mydata[common_cells, ]
      } else {
        stop("No matching cell IDs found between mymeta and mydata")
      }
    } else {
      stop("No matching cell IDs found between mymeta and mydata")
    }
  } else {
    mymeta_matched <- mymeta[common_cells, ]
    mydata_matched <- mydata[common_cells, ]
  }

  if (verbose) {
    cat("Matched cells between metadata and expression:", length(common_cells), "\n")
  }

  # Determine cell types for pseudobulk creation
  if (is.null(cell_types_for_pseudobulk)) {
    cell_types_for_pseudobulk <- unique(mymeta_matched$Annotation_pooled)
    cell_types_for_pseudobulk <- cell_types_for_pseudobulk[!is.na(cell_types_for_pseudobulk)]
  }

  if (verbose) {
    cat("Cell types for pseudobulk creation:\n")
    print(cell_types_for_pseudobulk)
    cat("Available windows:", length(unique(mymeta_matched$window_id)), "\n")
  }

  # Get all unique windows
  window_ids <- unique(mymeta_matched$window_id)
  window_ids <- window_ids[!is.na(window_ids)]

  # Initialize pseudobulk data structures
  pseudobulk_expr_list <- list()
  sample_metadata_list <- list()
  cell_type_summary <- data.frame()

  # Create pseudobulk samples
  total_samples_created <- 0

  for (cell_type in cell_types_for_pseudobulk) {
    for (window_id in window_ids) {
      # Get cells of this type in this window
      type_window_cells <- mymeta_matched[mymeta_matched$Annotation_pooled == cell_type &
                                            mymeta_matched$window_id == window_id, ]

      if (nrow(type_window_cells) >= min_cells_per_sample) {
        cell_ids <- rownames(type_window_cells)
        sample_name <- paste(cell_type, window_id, sep = "_")

        # Sum expression across all cells of this type in this window
        if (length(cell_ids) == 1) {
          pseudobulk_expr_sample <- as.numeric(mydata_matched[cell_ids, ])
        } else {
          pseudobulk_expr_sample <- colSums(mydata_matched[cell_ids, , drop = FALSE])
        }
        
        # Apply normalization/scaling based on parameters
        if (!is.null(cell_factor)) {
          # Scale expression to what it would be with cell_factor cells
          # Formula: (total_expression / actual_cell_count) * cell_factor
          pseudobulk_expr_sample <- (pseudobulk_expr_sample / length(cell_ids)) * cell_factor
        } else if (normalize_by_cell_count) {
          # Normalize by cell count to get average expression per cell
          pseudobulk_expr_sample <- pseudobulk_expr_sample / length(cell_ids)
        }
        # If neither parameter is set, keep total sum (original behavior)

        # Store pseudobulk expression
        pseudobulk_expr_list[[sample_name]] <- pseudobulk_expr_sample

        # Create sample metadata
        sample_meta <- data.frame(
          sample_name = sample_name,
          cell_type = cell_type,
          window_id = window_id,
          cell_count = nrow(type_window_cells),
          stringsAsFactors = FALSE
        )

        # Extract FOV information from window_id if possible
        if (grepl("FOV", window_id)) {
          fov_part <- strsplit(window_id, "_")[[1]][1]
          sample_meta$fov_id <- as.numeric(gsub("FOV", "", fov_part))
        } else {
          sample_meta$fov_id <- NA
        }
        
        # Aggregate specified metadata columns (calculate mean)
        if (!is.null(aggregate_meta_cols)) {
          for (meta_col in aggregate_meta_cols) {
            # Get values for this metadata column from cells in this pseudobulk sample
            meta_values <- type_window_cells[[meta_col]]
            
            # Convert to numeric if needed and calculate mean
            if (!is.numeric(meta_values)) {
              meta_values <- as.numeric(as.character(meta_values))
            }
            
            # Calculate mean, handling NA values
            mean_value <- mean(meta_values, na.rm = TRUE)
            
            # Add to sample metadata with descriptive column name
            sample_meta[[paste0("mean_", meta_col)]] <- mean_value
          }
        }

        sample_metadata_list[[sample_name]] <- sample_meta

        # Update summary
        summary_row <- data.frame(
          cell_type = cell_type,
          window_id = window_id,
          cell_count = nrow(type_window_cells),
          stringsAsFactors = FALSE
        )
        cell_type_summary <- rbind(cell_type_summary, summary_row)

        total_samples_created <- total_samples_created + 1

        if (verbose && total_samples_created %% 100 == 0) {
          cat("Created", total_samples_created, "pseudobulk samples so far...\n")
        }
      }
    }
  }

  if (length(pseudobulk_expr_list) == 0) {
    stop("No pseudobulk samples could be created with the given parameters")
  }

  # Convert lists to matrices/data frames
  pseudobulk_expr <- do.call(cbind, pseudobulk_expr_list)
  rownames(pseudobulk_expr) <- colnames(mydata_matched)  # genes
  colnames(pseudobulk_expr) <- names(pseudobulk_expr_list)

  sample_metadata <- do.call(rbind, sample_metadata_list)
  rownames(sample_metadata) <- sample_metadata$sample_name

  if (verbose) {
    cat("Pseudobulk creation completed:\n")
    cat("- Total samples created:", ncol(pseudobulk_expr), "\n")
    cat("- Genes:", nrow(pseudobulk_expr), "\n")
    cat("- Samples per cell type:\n")
    print(table(sample_metadata$cell_type))
    
    if (!is.null(aggregate_meta_cols)) {
      cat("- Aggregated metadata columns added:\n")
      for (meta_col in aggregate_meta_cols) {
        new_col_name <- paste0("mean_", meta_col)
        if (new_col_name %in% colnames(sample_metadata)) {
          col_summary <- summary(sample_metadata[[new_col_name]])
          cat("  ", new_col_name, ": Min =", round(col_summary[1], 3), 
              ", Mean =", round(col_summary[4], 3), 
              ", Max =", round(col_summary[6], 3), "\n")
        }
      }
    }
  }

  # Return results as a list
  return(list(
    pseudobulk_expr = pseudobulk_expr,
    sample_metadata = sample_metadata,
    cell_type_summary = cell_type_summary
  ))
}

#' Plot Sliding Windows Overlay
#'
#' This function creates an overlay visualization showing cell distribution (bottom layer)
#' with fixed-size sliding window boundaries (top layer) as rectangles. Uses the actual
#' window boundaries from CreateSlidingWindow output to show true fixed-size windows.
#'
#' @param mymeta Data frame containing cell metadata with columns: fov, X_cent, Y_cent, Annotation_pooled, window_id
#' @param window_info_df Data frame containing window boundary information from CreateSlidingWindow (columns: window_id, fov_id, x_start, y_start, x_end, y_end)
#' @param fov_id Character or numeric, specific FOV identifier to visualize (default: NULL, plots all FOVs)
#' @param cell_size Numeric, size of cell points in the plot (default: 0.5)
#' @param window_color Character, color for window rectangles (default: "black")
#' @param window_linewidth Numeric, line width for window rectangles (default: 0.8)
#' @param window_alpha Numeric, transparency of window rectangles (0-1, default: 1)
#' @param show_legend Logical, whether to show cell type legend (default: TRUE)
#' @param title Character, custom plot title (default: NULL, auto-generated)
#'
#' @return ggplot object showing cell distribution with fixed-size window overlay
#'
#' @examples
#' # Get sliding windows
#' result <- CreateSlidingWindow(mymeta, window_size = 150, stride = 100)
#' 
#' # Plot single FOV with fixed-size windows
#' p <- PlotSlidingWindowsOverlay(result$mymeta_windowed, result$all_windows, fov_id = "DFCI_c02")
#' print(p)
#' 
#' # Plot all FOVs
#' p <- PlotSlidingWindowsOverlay(result$mymeta_windowed, result$all_windows)
#' 
#' # Customize appearance
#' p <- PlotSlidingWindowsOverlay(result$mymeta_windowed, result$all_windows, 
#'                                fov_id = "DFCI_c02",
#'                                cell_size = 1, window_color = "red", 
#'                                window_linewidth = 1.2)
#'
PlotSlidingWindowsOverlay <- function(mymeta, 
                                      window_info_df = NULL,
                                      fov_id = NULL, 
                                      cell_size = 0.5,
                                      window_color = "black",
                                      window_linewidth = 0.8,
                                      window_alpha = 1,
                                      show_legend = TRUE,
                                      title = NULL) {
  
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  
  # Input validation
  if (!is.data.frame(mymeta)) {
    stop("mymeta must be a data frame")
  }
  
  required_cols <- c("fov", "X_cent", "Y_cent", "Annotation_pooled")
  missing_cols <- setdiff(required_cols, colnames(mymeta))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in mymeta:", paste(missing_cols, collapse = ", ")))
  }
  
  # Filter for specific FOV if requested
  if (!is.null(fov_id)) {
    plot_data <- mymeta[mymeta$fov == fov_id, ]
    if (nrow(plot_data) == 0) {
      stop(paste("No cells found for FOV:", fov_id))
    }
  } else {
    plot_data <- mymeta
  }
  
  if (nrow(plot_data) == 0) {
    stop("No cells found")
  }
  
  # Get window boundaries: either from window_info_df or extract from cell positions
  if (!is.null(window_info_df)) {
    # Use fixed-size windows from CreateSlidingWindow output
    if (!is.data.frame(window_info_df)) {
      stop("window_info_df must be a data frame")
    }
    
    required_window_cols <- c("window_id", "fov_id", "x_start", "y_start", "x_end", "y_end")
    missing_window_cols <- setdiff(required_window_cols, colnames(window_info_df))
    if (length(missing_window_cols) > 0) {
      stop(paste("Missing required columns in window_info_df:", paste(missing_window_cols, collapse = ", ")))
    }
    
    # Filter windows for FOV
    if (!is.null(fov_id)) {
      window_boundaries <- window_info_df[window_info_df$fov_id == fov_id, ]
    } else {
      window_boundaries <- window_info_df
    }
    
    if (nrow(window_boundaries) == 0) {
      stop("No windows found in window_info_df")
    }
    
    # Calculate uniform window size from the largest window (to ensure all windows same size)
    # This handles cases where edge windows might be truncated in the data
    window_sizes_x <- window_boundaries$x_end - window_boundaries$x_start
    window_sizes_y <- window_boundaries$y_end - window_boundaries$y_start
    uniform_window_size_x <- max(window_sizes_x)
    uniform_window_size_y <- max(window_sizes_y)
    
    # Create uniform fixed-size rectangles using x_start/y_start + uniform size
    # This ensures all windows are displayed with the same dimensions
    window_boundaries$x_min <- window_boundaries$x_start
    window_boundaries$x_max <- window_boundaries$x_start + uniform_window_size_x
    window_boundaries$y_min <- window_boundaries$y_start
    window_boundaries$y_max <- window_boundaries$y_start + uniform_window_size_y
    
    cat("Using", nrow(window_boundaries), "fixed-size windows (", 
        round(uniform_window_size_x, 1), "x", round(uniform_window_size_y, 1), 
        ") from window_info_df\n")
    
  } else {
    # Extract boundaries from cell assignments (old behavior)
    if (!"window_id" %in% colnames(mymeta)) {
      stop("window_id column required in mymeta when window_info_df is not provided")
    }
    
    # Load dplyr for grouping
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("dplyr package is required for this function")
    }
    
    `%>%` <- dplyr::`%>%`
    X_cent <- Y_cent <- window_id <- x_min <- x_max <- y_min <- y_max <- NULL  # Avoid R CMD check warnings
    
    # Remove cells without window assignment
    plot_data <- plot_data[!is.na(plot_data$window_id), ]
    
    if (nrow(plot_data) == 0) {
      stop("No cells with window assignments found")
    }
    
    # Extract window boundaries from cell positions
    window_boundaries <- plot_data %>%
      dplyr::group_by(window_id) %>%
      dplyr::summarise(
        x_min = min(X_cent),
        x_max = max(X_cent),
        y_min = min(Y_cent),
        y_max = max(Y_cent),
        cell_count = dplyr::n(),
        .groups = 'drop'
      )
    
    cat("Extracted", nrow(window_boundaries), "window boundaries from cell assignments\n")
  }
  
  # Create base plot - Layer 1: Cell distribution
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = X_cent, y = Y_cent, color = Annotation_pooled)) +
    ggplot2::geom_point(size = cell_size, alpha = 0.7) +
    ggplot2::scale_color_discrete(name = "Cell Type")
  
  # Add Layer 2: Fixed-size window rectangles on top
  p <- p + ggplot2::geom_rect(
    data = window_boundaries,
    ggplot2::aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
    fill = NA, 
    color = window_color, 
    linewidth = window_linewidth, 
    alpha = window_alpha,
    inherit.aes = FALSE
  )
  
  # Customize plot appearance
  if (is.null(title)) {
    if (!is.null(fov_id)) {
      title <- paste0("Fixed-Size Sliding Windows for ", fov_id, 
                     "\n", nrow(window_boundaries), " windows, ", 
                     nrow(plot_data), " cells")
    } else {
      n_fovs <- length(unique(plot_data$fov))
      title <- paste0("Fixed-Size Sliding Windows\n", 
                     n_fovs, " FOVs, ",
                     nrow(window_boundaries), " windows, ", 
                     nrow(plot_data), " cells")
    }
  }
  
  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = if (show_legend) "right" else "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}

#' Plot All FOVs with Sliding Windows Overlay
#'
#' This function creates separate overlay plots for each FOV in the dataset,
#' showing cell distribution with fixed-size sliding window boundaries. Returns either
#' a combined multi-panel plot or a list of individual plots.
#'
#' @param mymeta Data frame containing cell metadata with columns: fov, X_cent, Y_cent, Annotation_pooled, window_id
#' @param window_info_df Data frame containing window boundary information from CreateSlidingWindow (columns: window_id, fov_id, x_start, y_start, x_end, y_end)
#' @param cell_size Numeric, size of cell points in the plot (default: 0.5)
#' @param window_color Character, color for window rectangles (default: "black")
#' @param window_linewidth Numeric, line width for window rectangles (default: 0.8)
#' @param show_legend Logical, whether to show cell type legend (default: TRUE)
#' @param ncol Numeric, number of columns in multi-panel layout (default: NULL, auto-determined)
#' @param return_list Logical, whether to return a list of individual plots (default: FALSE)
#'
#' @return If return_list=FALSE, returns a combined ggplot object with all FOVs.
#'         If return_list=TRUE, returns a list of ggplot objects, one per FOV.
#'
#' @examples
#' # Get sliding windows
#' result <- CreateSlidingWindow(mymeta, window_size = 150, stride = 100)
#' 
#' # Create combined multi-panel plot
#' p <- PlotAllFOVsOverlay(result$mymeta_windowed, result$all_windows)
#' print(p)
#' 
#' # Get individual plots
#' plot_list <- PlotAllFOVsOverlay(result$mymeta_windowed, result$all_windows, return_list = TRUE)
#' print(plot_list[[1]])
#'
PlotAllFOVsOverlay <- function(mymeta,
                               window_info_df = NULL,
                               cell_size = 0.5,
                               window_color = "black",
                               window_linewidth = 0.8,
                               show_legend = TRUE,
                               ncol = NULL,
                               return_list = FALSE) {
  
  # Load required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  
  # Input validation
  if (!is.data.frame(mymeta)) {
    stop("mymeta must be a data frame")
  }
  
  # Get all unique FOVs
  all_fovs <- sort(unique(mymeta$fov))
  
  if (length(all_fovs) == 0) {
    stop("No FOVs found in mymeta")
  }
  
  cat("Creating overlay plots for", length(all_fovs), "FOVs...\n")
  
  # Create individual plots for each FOV
  plot_list <- list()
  
  for (i in seq_along(all_fovs)) {
    fov_id <- all_fovs[i]
    
    tryCatch({
      p <- PlotSlidingWindowsOverlay(
        mymeta = mymeta,
        window_info_df = window_info_df,
        fov_id = fov_id,
        cell_size = cell_size,
        window_color = window_color,
        window_linewidth = window_linewidth,
        show_legend = show_legend,
        title = NULL
      )
      
      plot_list[[as.character(fov_id)]] <- p
      
      if (i %% 5 == 0) {
        cat("  Created plots for", i, "/", length(all_fovs), "FOVs\n")
      }
    }, error = function(e) {
      warning(paste("Could not create plot for FOV", fov_id, ":", e$message))
    })
  }
  
  if (length(plot_list) == 0) {
    stop("No plots could be created")
  }
  
  cat("Successfully created", length(plot_list), "FOV plots\n")
  
  # Return list if requested
  if (return_list) {
    return(plot_list)
  }
  
  # Otherwise, combine plots into multi-panel figure
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork package not available. Returning list of plots instead.")
    warning("Install patchwork with: install.packages('patchwork')")
    return(plot_list)
  }
  
  # Determine number of columns for layout
  if (is.null(ncol)) {
    n_plots <- length(plot_list)
    if (n_plots <= 2) {
      ncol <- n_plots
    } else if (n_plots <= 4) {
      ncol <- 2
    } else if (n_plots <= 9) {
      ncol <- 3
    } else {
      ncol <- 4
    }
  }
  
  # Combine plots using patchwork
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol) +
    patchwork::plot_annotation(
      title = paste("Sliding Windows Overlay -", length(plot_list), "FOVs"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  return(combined_plot)
}
