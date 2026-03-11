#  Figure 6B and Figure S9E for codex_cosmx

# %% Library Loading =====

SEPARATOR <- paste(rep("=", 80), collapse = "")

cat("\n", SEPARATOR, "\n", sep = "")
cat("Library Loading\n")
cat(SEPARATOR, "\n", sep = "")

library(tidyverse)
library(reshape2)
library(corrplot)
library(patchwork)

cat("\nLibraries loaded successfully\n")

# %% Configurations =====

cat("\n", SEPARATOR, "\n", sep = "")
cat("Configurations\n")
cat(SEPARATOR, "\n", sep = "")

script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
data_root <- file.path(dirname(script_path), "src", "06_figure_6_correlation","codex_cosmx")
cat("\nData root:", data_root, "\n")

for (tag in c("core", "fov")) {
  cat("\n", SEPARATOR, "\n", sep = "")
  cat("Processing:", toupper(tag), "level analysis\n")
  cat(SEPARATOR, "\n", sep = "")
  
  # %% Load Data =====
  
  cat("\n--- Load Data ---\n")
  
  output_dir <- fs::path(data_root, paste0("02_correlation_", tag))
  df_score_f <- fs::path(data_root, paste0("01_scores_", tag, ".csv"))
  
  fs::dir_create(output_dir, recurse = TRUE)
  cat("Output directory:", as.character(output_dir), "\n")
  cat("Input file:", as.character(df_score_f), "\n")
  
  df_score <- read_csv(df_score_f, show_col_types = FALSE)
  df_score <- df_score %>%
    mutate(data_tag = str_extract(region_id, "^[^_]+"))
  
  cat("\nData loaded successfully\n")
  cat("   - Total regions:", nrow(df_score), "\n")
  cat("   - Total columns:", ncol(df_score), "\n")
  cat("   - Data tags:", paste(unique(df_score$data_tag), collapse = ", "), "\n")
  
  
  # %% Spearman Correlation Analysis =====
  
  cat("\n--- Spearman Correlation Analysis ---\n")
  
  score_columns <- c(
    "Dysfunction_CD4T", "Dysfunction_CD8T", "Dysfunction_T",
    "C1QC-Mac_Macrophage",
    "ebv_burden_B"
  )
  
  cat("\nScore columns for correlation analysis:\n")
  for (col in score_columns) {
    cat("   -", col, "\n")
  }
  
  df_score_selected <- df_score %>% select(all_of(score_columns))
  
  # Calculate Spearman correlation and p-values
  cat("\nComputing Spearman correlation...\n")
  cor_spearman <- cor(df_score_selected, use = "pairwise.complete.obs", method = "spearman")
  p_spearman <- cor.mtest(df_score_selected, method = "spearman")$p
  
  # Count significant correlations
  n_significant <- sum(p_spearman < 0.05, na.rm = TRUE) - nrow(p_spearman) # exclude diagonal
  cat("\nCorrelation analysis complete\n")
  cat("   - Total comparisons:", (nrow(cor_spearman)^2 - nrow(cor_spearman)) / 2, "\n")
  cat("   - Significant correlations (p < 0.05):", n_significant / 2, "\n")
  
  cat("\nSpearman Correlation Matrix:\n")
  print(cor_spearman)
  
  cat("\nSpearman p-values:\n")
  print(p_spearman)
  
  # %% Correlation Matrix Plot =====
  
  cat("\n--- Correlation Matrix Plot ---\n")
  
  cor_melt_sp <- melt(cor_spearman)
  p_melt_sp <- melt(p_spearman)
  colnames(cor_melt_sp) <- c("Var1", "Var2", "correlation")
  colnames(p_melt_sp) <- c("Var1", "Var2", "p_value")
  
  plot_data_sp <- cor_melt_sp %>%
    left_join(p_melt_sp, by = c("Var1", "Var2")) %>%
    mutate(significant = p_value < 0.05)
  
  p <- ggplot(plot_data_sp, aes(x = Var1, y = Var2, fill = correlation)) +
    geom_tile(color = "white") +
    geom_tile(
      data = filter(plot_data_sp, significant),
      aes(x = Var1, y = Var2),
      color = "black",
      linewidth = 1.5,
      fill = NA
    ) +
    geom_text(aes(label = sprintf("%.2f", correlation)),
              size = 3,
              color = "black"
    ) +
    scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      limits = c(-1, 1),
      name = "cor"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = "Spearman Correlation Matrix (black box: p < 0.05)") +
    coord_fixed()
  
  p
  
  output_file <- fs::path(output_dir, "01_correlation_matrix.pdf")
  ggsave(output_file, p, width = 6, height = 6)
  cat("\nCorrelation matrix plot saved to:", as.character(output_file), "\n")
  
  
  # %% Scatter Plots =====
  
  create_scatter_plot <- function(df, x_var, y_var, cor_mat, p_mat, color_by = NULL) {
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    
    # Data preparation
    if (is.null(color_by)) {
      df_plot <- df %>% select(x = all_of(x_var), y = all_of(y_var))
    } else {
      df_plot <- df %>% select(x = all_of(x_var), y = all_of(y_var), source = all_of(color_by))
    }
    
    # create figures
    df_plot <- df_plot %>% drop_na()
    p <- ggplot(df_plot, aes(x = x, y = y))
    
    if (!is.null(color_by)) {
      p <- p + geom_point(aes(color = source))
    } else {
      p <- p + geom_point()
    }
    
    p <- p +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
      labs(
        x = x_var,
        y = y_var,
        title = sprintf(
          "Spearman: cor = %.3f, p = %.3f",
          cor_mat[x_var, y_var],
          p_mat[x_var, y_var]
        )
      ) +
      theme_bw() +
      theme(legend.position = "bottom", legend.box = "horizontal")
    
    return(p)
  }
  
  xy <- list(
    c("C1QC-Mac_Macrophage", "Dysfunction_CD4T"),
    c("C1QC-Mac_Macrophage", "Dysfunction_CD8T"),
    c("ebv_burden_B", "C1QC-Mac_Macrophage")
  )
  
  cat("\nGenerating", length(xy), "scatter plots...\n")

  for (i in seq_along(xy)) {
    x_var <- xy[[i]][1]
    y_var <- xy[[i]][2]

    cat("   - Plot", i, ":", x_var, "vs", y_var, "\n")

    p <- create_scatter_plot(
      df = df_score,
      x_var = x_var,
      y_var = y_var,
      cor_mat = cor_spearman,
      p_mat = p_spearman, color_by = "data_tag"
    )

    ggsave(
      fs::path(output_dir, sprintf("03_scatter_plot_%s_vs_%s.pdf", x_var, y_var)),
      p,
      width = 5,
      height = 5
    )
  }

  cat("\nAll individual scatter plots saved\n")


  cat("\n", SEPARATOR, "\n", sep = "")
  cat(toupper(tag), "level analysis complete\n")
  cat(SEPARATOR, "\n", sep = "")
}

cat("\n", SEPARATOR, "\n", sep = "")
cat("ALL ANALYSES COMPLETE\n")
cat(SEPARATOR, "\n", sep = "")



#  Figure 6B for cosmx_only
# %% Library Loading =====

SEPARATOR <- paste(rep("=", 80), collapse = "")

cat("\n", SEPARATOR, "\n", sep = "")
cat("Library Loading\n")
cat(SEPARATOR, "\n", sep = "")

library(tidyverse)
library(reshape2)
library(corrplot)
library(patchwork)

cat("\nLibraries loaded successfully\n")

# %% Configurations =====

cat("\n", SEPARATOR, "\n", sep = "")
cat("Configurations\n")
cat(SEPARATOR, "\n", sep = "")

script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
data_root <- file.path(dirname(script_path), "src", "06_figure_6_correlation","cosmx_only")
cat("\nData root:", data_root, "\n")


for (tag in c("fov")) {
  cat("\n", SEPARATOR, "\n", sep = "")
  cat("Processing:", toupper(tag), "level analysis\n")
  cat(SEPARATOR, "\n", sep = "")
  
  # %% Load Data =====
  
  cat("\n--- Load Data ---\n")
  
  output_dir <- fs::path(data_root, paste0("02_correlation_", tag))
  df_score_f <- fs::path(data_root, paste0("01_scores_", tag, ".csv"))
  
  fs::dir_create(output_dir, recurse = TRUE)
  cat("Output directory:", as.character(output_dir), "\n")
  cat("Input file:", as.character(df_score_f), "\n")
  
  df_score <- read_csv(df_score_f, show_col_types = FALSE)
  df_score <- df_score %>%
    mutate(data_tag = str_extract(region_id, "^[^_]+"))
  
  cat("\nData loaded successfully\n")
  cat("   - Total regions:", nrow(df_score), "\n")
  cat("   - Total columns:", ncol(df_score), "\n")
  cat("   - Data tags:", paste(unique(df_score$data_tag), collapse = ", "), "\n")
  
  
  # %% Spearman Correlation Analysis =====
  
  cat("\n--- Spearman Correlation Analysis ---\n")
  
  score_columns <- c(
    "Dysfunction_T",
    "C1QC-Mac_Macrophage",
    "ebv_burden_B"
  )
  
  cat("\nScore columns for correlation analysis:\n")
  for (col in score_columns) {
    cat("   -", col, "\n")
  }
  
  df_score_selected <- df_score %>% select(all_of(score_columns))
  
  # Calculate Spearman correlation and p-values
  cat("\nComputing Spearman correlation...\n")
  cor_spearman <- cor(df_score_selected, use = "pairwise.complete.obs", method = "spearman")
  p_spearman <- cor.mtest(df_score_selected, method = "spearman")$p
  
  # Count significant correlations
  n_significant <- sum(p_spearman < 0.05, na.rm = TRUE) - nrow(p_spearman) # exclude diagonal
  cat("\nCorrelation analysis complete\n")
  cat("   - Total comparisons:", (nrow(cor_spearman)^2 - nrow(cor_spearman)) / 2, "\n")
  cat("   - Significant correlations (p < 0.05):", n_significant / 2, "\n")
  
  cat("\nSpearman Correlation Matrix:\n")
  print(cor_spearman)
  
  cat("\nSpearman p-values:\n")
  print(p_spearman)
  
  # %% Correlation Matrix Plot =====
  
  cat("\n--- Correlation Matrix Plot ---\n")
  
  cor_melt_sp <- melt(cor_spearman)
  p_melt_sp <- melt(p_spearman)
  colnames(cor_melt_sp) <- c("Var1", "Var2", "correlation")
  colnames(p_melt_sp) <- c("Var1", "Var2", "p_value")
  
  plot_data_sp <- cor_melt_sp %>%
    left_join(p_melt_sp, by = c("Var1", "Var2")) %>%
    mutate(significant = p_value < 0.05)
  
  p <- ggplot(plot_data_sp, aes(x = Var1, y = Var2, fill = correlation)) +
    geom_tile(color = "white") +
    geom_tile(
      data = filter(plot_data_sp, significant),
      aes(x = Var1, y = Var2),
      color = "black",
      linewidth = 1.5,
      fill = NA
    ) +
    geom_text(aes(label = sprintf("%.2f", correlation)),
              size = 3,
              color = "black"
    ) +
    scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      limits = c(-1, 1),
      name = "cor"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = "Spearman Correlation Matrix (black box: p < 0.05)") +
    coord_fixed()
  
  p
  
  output_file <- fs::path(output_dir, "01_correlation_matrix.pdf")
  ggsave(output_file, p, width = 6, height = 6)
  cat("\nCorrelation matrix plot saved to:", as.character(output_file), "\n")
  
  
  # %% Scatter Plots =====
  
  create_scatter_plot <- function(df, x_var, y_var, cor_mat, p_mat, color_by = NULL) {
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    
    # Data preparation
    if (is.null(color_by)) {
      df_plot <- df %>% select(x = all_of(x_var), y = all_of(y_var))
    } else {
      df_plot <- df %>% select(x = all_of(x_var), y = all_of(y_var), source = all_of(color_by))
    }
    
    #  create figures
    df_plot <- df_plot %>% drop_na()
    p <- ggplot(df_plot, aes(x = x, y = y))
    
    if (!is.null(color_by)) {
      p <- p + geom_point(aes(color = source))
    } else {
      p <- p + geom_point()
    }
    
    p <- p +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      scale_y_continuous(labels = function(x) sprintf("%.3f", x)) +
      labs(
        x = x_var,
        y = y_var,
        title = sprintf(
          "Spearman: cor = %.3f, p = %.3f",
          cor_mat[x_var, y_var],
          p_mat[x_var, y_var]
        )
      ) +
      theme_bw() +
      theme(legend.position = "bottom", legend.box = "horizontal")
    
    return(p)
  }
  
  xy <- list(
    c("C1QC-Mac_Macrophage", "Dysfunction_T"),
    c("ebv_burden_B", "C1QC-Mac_Macrophage")
  )
  
  cat("\nGenerating", length(xy), "scatter plots...\n")

  for (i in seq_along(xy)) {
    x_var <- xy[[i]][1]
    y_var <- xy[[i]][2]

    cat("   - Plot", i, ":", x_var, "vs", y_var, "\n")

    p <- create_scatter_plot(
      df = df_score,
      x_var = x_var,
      y_var = y_var,
      cor_mat = cor_spearman,
      p_mat = p_spearman, color_by = "data_tag"
    )

    ggsave(
      fs::path(output_dir, sprintf("03_scatter_plot_%s_vs_%s.pdf", x_var, y_var)),
      p,
      width = 5,
      height = 5
    )
  }

  cat("\nAll individual scatter plots saved\n")
  
  cat("\n", SEPARATOR, "\n", sep = "")
  cat(toupper(tag), "level analysis complete\n")
  cat(SEPARATOR, "\n", sep = "")
}

cat("\n", SEPARATOR, "\n", sep = "")
cat("ALL ANALYSES COMPLETE\n")
cat(SEPARATOR, "\n", sep = "")


