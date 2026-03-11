library(tidyverse)
library(arrow)
library(patchwork)

# Define ligand-receptor pairs
lr_pairs <- c(
  "CCL22::CCR4",
  "EBI3::IL6ST",
  "EBI3::IL27RA",
  "IL6::IL6ST",
  "EBI3::STAT3",
  "ICOSLG::ICOS",
  "ICAM3::CD209",
  "DCN::TLR4",
  "CCN1::TLR4"
)

# Set data directory
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
data_dir <- file.path(dirname(script_path), "src", "07_figure_7_ neighborhood_analysis")

# Read metadata
df_metadata <- arrow::read_feather(file.path(data_dir, "07_aucell_scores.feather")) %>%
  select(data_tag, core_id, cell_id, fov_id, smp_id) %>%
  mutate(
    cell_id = str_c(smp_id, cell_id, sep = "_"),
    smp_id_fov = str_c(smp_id, fov_id, sep = "_")
  ) %>%
  select(cell_id, smp_id_fov)

# Read macrophage L-R mean data
df_mac_lr_mean <- read_feather(file.path(data_dir, "06_macrophage_center_lr_mean.feather")) %>%
  select(cell_id_mac, lmp1, all_of(lr_pairs)) %>%
  filter(!is.na(lmp1)) %>%
  left_join(df_metadata, by = c("cell_id_mac" = "cell_id"))

# Read tumor-macrophage L-R mean data
df_lr_mean_tumor_mac <- read_feather(file.path(data_dir, "06_tumor_macrophage_lr_mean.feather")) %>%
  select(cell_id_tumor, cell_id_mac, lmp1, all_of(lr_pairs)) %>%
  filter(!is.na(lmp1)) %>%
  left_join(df_metadata, by = c("cell_id_mac" = "cell_id"))

# Function to format p-values
format_pval <- function(p) {
  if (p < 0.0001) {
    return(sprintf("%.2e", p))
  } else {
    return(sprintf("%.4f", p))
  }
}

# Function to create violin+boxplot for each L-R pair (single-cell level)
create_lr_plot <- function(data, lr_pair, title_prefix = "") {
  library(ggplot2)
  library(dplyr)
  
  # Prepare data for plotting
  plot_data <- data %>%
    select(!!sym("lmp1"), value = all_of(lr_pair)) %>%
    filter(!is.na(!!sym("value")))
  
  # Perform statistical tests
  lmp1p_values <- plot_data %>%
    filter(!!sym("lmp1") == "lmp1p") %>%
    pull(!!sym("value"))
  lmp1n_values <- plot_data %>%
    filter(!!sym("lmp1") == "lmp1n") %>%
    pull(!!sym("value"))
  
  
  # Wilcoxon test
  wilcox_result <- wilcox.test(lmp1p_values, lmp1n_values)
  wilcox_pval <- wilcox_result$p.value
  
  # Create title with p-values
  title_text <- sprintf(
    "%s%s\nWilcoxon p = %s",
    title_prefix,
    lr_pair,
    format_pval(wilcox_pval)
  )

  # Create plot
  p <- ggplot(plot_data, aes(x = !!sym("lmp1"), y = !!sym("value"), fill = !!sym("lmp1"))) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    scale_fill_manual(values = c("lmp1p" = "#E64B35", "lmp1n" = "#4DBBD5")) +
    theme_bw() +
    labs(
      title = title_text,
      x = "LMP1 Status",
      y = "Mean Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )
  
  return(p)
}

# Function to create pseudobulk violin+boxplot+scatter for each L-R pair
create_lr_plot_pseudobulk <- function(data, lr_pair, title_prefix = "", seed = 42) {
  library(ggplot2)
  library(dplyr)
  
  # Aggregate by pseudobulk (smp_id_fov)
  plot_data <- data %>%
    select(!!sym("lmp1"), !!sym("smp_id_fov"), value = all_of(lr_pair)) %>%
    filter(!is.na(!!sym("value"))) %>%
    group_by(!!sym("lmp1"), !!sym("smp_id_fov")) %>%
    summarise(mean_value = mean(!!sym("value"), na.rm = TRUE), .groups = "drop")
  
  # Perform statistical tests
  lmp1p_values <- plot_data %>%
    filter(!!sym("lmp1") == "lmp1p") %>%
    pull(!!sym("mean_value"))
  lmp1n_values <- plot_data %>%
    filter(!!sym("lmp1") == "lmp1n") %>%
    pull(!!sym("mean_value"))
  
  
  # Wilcoxon test
  wilcox_result <- wilcox.test(lmp1p_values, lmp1n_values)
  wilcox_pval <- wilcox_result$p.value
  
  # Create title with p-values
  title_text <- sprintf(
    "%s%s (Pseudobulk)\nWilcoxon p = %s",
    title_prefix,
    lr_pair,
    format_pval(wilcox_pval)
  )
  
  # Set seed for reproducible jitter
  set.seed(seed)
  
  # Create plot with scatter points
  p <- ggplot(plot_data, aes(x = !!sym("lmp1"), y = !!sym("mean_value"), fill = !!sym("lmp1"))) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, shape = 21, color = "black") +
    scale_fill_manual(values = c("lmp1p" = "#E64B35", "lmp1n" = "#4DBBD5")) +
    theme_bw() +
    labs(
      title = title_text,
      x = "LMP1 Status",
      y = "Mean Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )
  
  return(p)
}

# %% ========== Pseudobulk Analysis ==========

cat("\n=== Pseudobulk Analysis ===\n")

# Create pseudobulk plots for df_lr_mean_tumor_mac
cat("Creating pseudobulk plots for tumor-macrophage L-R pairs...\n")
plot_list_tumor_mac_pb <- list()
for (lr_pair in lr_pairs) {
  plot_list_tumor_mac_pb[[lr_pair]] <- create_lr_plot_pseudobulk(
    df_lr_mean_tumor_mac,
    lr_pair,
    title_prefix = "T-M: ",
    seed = 42
  )
}

# Combine pseudobulk plots using patchwork
combined_plot_tumor_mac_pb <- wrap_plots(plot_list_tumor_mac_pb, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Tumor-Macrophage L-R Pair Expression (Pseudobulk): LMP1+ vs LMP1-",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# Display pseudobulk plots
print(combined_plot_tumor_mac_pb)

# %% ========== Save Figures ==========

output_dir <- file.path(data_dir, "figures", "05_single_cell_lr_comparison")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


ggsave(
  filename = file.path(output_dir, "02_tumor_macrophage_lr_pair_pseudobulk.pdf"),
  plot = combined_plot_tumor_mac_pb,
  width = 15,
  height = 12,
)
