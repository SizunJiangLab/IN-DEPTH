# Squidpy Ligand-Receptor Analysis Visualization
#
# This script visualizes ligand-receptor interaction results from Squidpy analysis
# for tumor-macrophage interactions, comparing LMP1+ and LMP1- conditions.
#
# Output: Dot plot showing L-R pair mean expression levels colored by significance

# %% Load Libraries ==========
library(tidyverse)

# %% Data Paths ==========
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
data_dir <- file.path(dirname(script_path), "src", "07_figure_7_ neighborhood_analysis")

# %% Define L-R Pairs of Interest ==========

# Ligand-receptor pairs selected for tumor-macrophage interaction analysis
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


# %% Load and Filter Data ==========

# Load Squidpy results containing L-R interaction analysis
df_squidpy <- read_csv(fs::path(data_dir, "05_squidpy_res_df_merged.csv"))

df_squidpy <- df_squidpy %>%
  filter(interaction == "tumor::macrophage", lr_pair %in% lr_pairs) %>%
  arrange(desc(log2fc))

df_squidpy$lr_pair <- factor(df_squidpy$lr_pair, levels = unique(df_squidpy$lr_pair))


# %% Reshape Data for Visualization ==========

# Transform data from wide to long format for easier plotting
df_squidpy_long <- df_squidpy %>%
  select(lr_pair, contains("lr_means_"), contains("pvalue_")) %>%
  pivot_longer(-lr_pair) %>%
  mutate(lmp1 = ifelse(str_detect(name, "lmp1p"), "lmp1p", "lmp1n"))

# Extract p-values and determine significance (p < 0.05)
df_squidpy_pval <- df_squidpy_long %>%
  filter(name %>% str_detect("pvalue_")) %>%
  mutate(psig = value < 0.05) %>%
  rename(pval = value)

# Extract mean expression values
df_squidpy_mean <- df_squidpy_long %>%
  filter(name %>% str_detect("lr_means_")) %>%
  rename(mean = value)

# Combine mean expression and significance data
df_plot <- df_squidpy_mean %>%
  left_join(
    df_squidpy_pval %>%
      select(lr_pair, lmp1, pval, psig),
    by = c("lr_pair", "lmp1")
  )


# %% Create Visualization ==========

# Create dot plot with:
# - x-axis: LMP1 status (lmp1n vs lmp1p)
# - y-axis: L-R pairs (ordered by log2fc)
# - Fill color: Mean expression level (viridis colormap)
# - Point size: Mean expression level (larger = higher expression)
# - Border color: Statistical significance (red = p < 0.05, black = not significant)
p <- df_plot %>%
  ggplot(aes(x = lmp1, y = lr_pair)) +
  geom_point(aes(fill = mean, size = mean, color = psig), shape = 21) +
  scale_fill_viridis_c(name = "Mean\nExpression") +
  scale_size_continuous(name = "Mean\nExpression") +
  scale_color_manual(
    name = "Significant\n(p < 0.05)",
    values = c("TRUE" = "red", "FALSE" = "black")
  ) +
  labs(
    x = "LMP1 Status",
    y = "Ligand-Receptor Pair",
    title = "Tumor-Macrophage L-R Interaction Comparison"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# %% Save Output ==========

ggsave(
  fs::path(data_dir, "figures", "03_squidpy_lr_mean_comparison.pdf"),
  p,
  width = 4,
  height = 5
)

cat("Plot saved to:", fs::path(data_dir, "figures", "03_squidpy_lr_mean_comparison.pdf"), "\n")
