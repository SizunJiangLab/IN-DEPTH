#####################################################
# Title: Eigenvalue Spectra Analysis Across k Values
#
# Description: 
#   This script loads eigenvalue spectra data 
#   calculated for various k-nearest neighbor graphs,
#   reshapes the data, computes pairwise correlations, 
#   and visualizes the stability of the spectra 
#   as a function of k. It includes code for 
#   saving and loading data, and for generating 
#   summary statistics and plots.
#
# Author: Rongting Huang
# Date: 2025-06-25
#
#####################################################
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)

# set working directory
setwd("/Users/rongting/Documents/ResearchProject/IN-DEPTH/Figures/")

# Load the eigenvalue spectra data
eigvals_wide <- readRDS("/Users/rongting/Documents/ResearchProject/IN-DEPTH/data/eigvals_wide.rds")

# Convert to matrix (rows = k, columns = order)
eigvals_matrix <- as.matrix(eigvals_wide[,-1])
rownames(eigvals_matrix) <- eigvals_wide$name

# Calculate pairwise correlations
cor_matrix <- cor(t(eigvals_matrix))
pairwise_cor <- cor_matrix[lower.tri(cor_matrix)]


# Output results
mean(pairwise_cor)
summary(pairwise_cor)
hist(pairwise_cor, main = "Pairwise Correlations", xlab = "Correlation")


# 1. correlation Heatmap---------------------------------------------------------
# cor_matrix: correlation matrix (rows and columns are k values)
pheatmap(cor_matrix,
         cluster_rows = FALSE,  # Don't cluster, preserve order
         cluster_cols = FALSE,
         display_numbers = TRUE, # Show correlation values
         main = "Pairwise Correlation of Eigenvalue Spectra: Each k vs. Others",
         fontsize_number = 10)



## custom heatmap with ggplot2
# Suppose cor_matrix is your correlation matrix
# Set up k labels
k_labels <- rownames(cor_matrix)

# Convert matrix to tidy data frame
cor_df <- as.data.frame(as.table(cor_matrix)) %>%
  rename(k_row = Var1, k_col = Var2, correlation = Freq)

# Keep only the lower triangle (including or excluding diagonal)
cor_df <- cor_df %>%
  mutate(
    k_row_num = as.numeric(factor(k_row, levels = k_labels)),
    k_col_num = as.numeric(factor(k_col, levels = k_labels))
  ) %>%
  filter(k_row_num > k_col_num) # lower triangle, excluding diagonal
# Use >= if want to keep the diagonal: filter(k_row_num >= k_col_num)

# Plot
correlation_heatmap <- 
  ggplot(cor_df, aes(x = k_col, y = k_row, fill = correlation)) +
  geom_tile(color = NA) +  # No border
  geom_text(aes(label = sprintf("%.3f", correlation)), size = 3) + # Value labels
  scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
  # scale_fill_gradient2(low = "red", mid = "orange", high = "blue", na.value = "white") +
  scale_y_discrete(limits = rev(k_labels)) +  # Heatmap style: origin top-left
  labs(
    x = "",
    y = "",
    title = "Pairwise Correlation of Eigenvalue Spectra"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x labels
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "inside", # must set this
    legend.position.inside = c(0.75, 0.65), # (x, y) in [0, 1] relative to plot area
    legend.justification = c("left", "bottom"),
    # legend.background = element_rect(fill = "white", color = "gray80"),
    # legend.box.background = element_rect(color = "gray80"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave(
  filename = "Pairwise_Correlation_LowerTri_ggplot2.png",
  plot = correlation_heatmap,
  width = 7,
  height = 6
)





# ---------------------------------------------------------
# # 2. Correlation to Reference Spectrum (largest k)
# # Plot Mean Correlation to a Reference (e.g., Highest k)
# # Suppose eigvals_matrix rows = k (e.g., 'k_50', 'k_100', ...)
# # use the spectrum at the largest k as reference
# ref_k <- nrow(eigvals_matrix) # or match to 'k_600'
# ref_spectrum <- eigvals_matrix[ref_k, ]
# cor_to_ref <- apply(eigvals_matrix, 1, function(x) cor(x, ref_spectrum))
# 
# Extract numeric k for plotting
k_numeric <- as.numeric(sub("k_", "", rownames(eigvals_matrix)))
# 
# plot(k_numeric, cor_to_ref, type = "b", pch = 16,
#      xlab = "k", ylab = "Correlation to reference (k=600)",
#      main = "Correlation to Reference Spectrum Across k")
# abline(v = 400, col = "red", lty = 2)


# ---------------------------------------------------------
# 3. Stability of Eigenvalue Spectra

# Calculate SD for each eigenvalue rank across all k
sd_by_rank <- apply(eigvals_matrix, 2, sd)
plot(sd_by_rank, type="b", main = "SD for each eigenvalue rank across k",
     ylab="Standard Deviation", xlab="Eigenvalue Rank")

# Or: SD for each k group (low vs. high)

group_low <- which(k_numeric < 250)
group_high <- which(k_numeric >= 250)
sd_low1 <- apply(eigvals_matrix[group_low, ], 2, sd)
sd_high1 <- apply(eigvals_matrix[group_high, ], 2, sd)

group_low <- which(k_numeric < 400)
group_high <- which(k_numeric >= 400)
sd_low <- apply(eigvals_matrix[group_low, ], 2, sd)
sd_high <- apply(eigvals_matrix[group_high, ], 2, sd)

plot(sd_low, type="b", main = "SD for each eigenvalue rank across k",
     ylab="Standard Deviation", xlab="Eigenvalue Rank")

plot(sd_high, type="b", main = "SD for each eigenvalue rank across k",
     ylab="Standard Deviation", xlab="Eigenvalue Rank")

barplot(c(mean(sd_low), mean(sd_high)), names.arg = c("k < 400", "k >= 400"),
        main = "Mean SD of Eigenvalues: Low k vs. High k",
        ylab = "Mean SD across eigenvalue ranks", col = c("orange", "skyblue"))



## custom boxplot with ggplot2

# Prepare data
# sd_df <- tibble(
#   EigenvalueRank = rep(seq_along(sd_low), 2),
#   SD = c(sd_low, sd_high),
#   Group = rep(c("k < 400", "k ≥ 400"), each = length(sd_low))
# )

sd_df <- tibble(
  EigenvalueRank = rep(seq_along(sd_low), 4),
  SD = c(sd_low, sd_high, sd_low1, sd_high1),
  Group = rep(c("k < 400", "k ≥ 400", "k < 250", "k ≥ 250"), each = length(sd_low))
)

# Set up the order you want for the group labels
group_order <- c("k < 400", "k ≥ 400", "k < 250", "k ≥ 250")

# Make Group a factor with specified levels
sd_df <- sd_df %>%
  mutate(Group = factor(Group, levels = group_order))


# Violin + Boxplot
sd_boxplot <- ggplot(sd_df, aes(x = Group, y = SD, fill = Group)) +
  # geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  scale_fill_manual(values = c("orange", "skyblue", "gold", "dodgerblue")) +
  labs(
    title = "Distribution of SD across Eigenvalue Ranks",
    subtitle = "By (higher/lower) k group",
    x = "",
    y = "Standard Deviation"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "inside", # must set this
          legend.position.inside = c(0.75, 0.65), # (x, y) in [0, 1] relative to plot area
          legend.justification = c("left", "bottom"))

ggsave(
  filename = "Distribution_of_SD_across_Eigenvalue_Ranks.png",
  plot = sd_boxplot,
  width = 7,
  height = 6
)

# sd_by_rank_df <- tibble(
#   EigenvalueRank = seq_along(sd_low),
#   `k < 400` = sd_low,
#   `k ≥ 400` = sd_high
# ) %>%
#   pivot_longer(-EigenvalueRank, names_to = "Group", values_to = "SD")
# 
# ggplot(sd_by_rank_df, aes(x = EigenvalueRank, y = SD, color = Group)) +
#   geom_line(size = 1.2) +
#   geom_point(size = 1.8) +
#   scale_color_manual(values = c("orange", "skyblue")) +
#   labs(
#     title = "SD by Eigenvalue Rank",
#     subtitle = "Lines: Each k group",
#     x = "Eigenvalue Rank",
#     y = "Standard Deviation"
#   ) +
#   theme_minimal(base_size = 15)



sd_by_rank_df <- tibble(
  EigenvalueRank = seq_along(sd_low),
  `k < 400` = sd_low,
  `k ≥ 400` = sd_high,
  `k < 250` = sd_low1,
  `k ≥ 250` = sd_high1
) %>%
  pivot_longer(-EigenvalueRank, names_to = "Group", values_to = "SD")

ggplot(sd_by_rank_df, aes(x = EigenvalueRank, y = SD, color = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c("orange", "skyblue", "gold", "dodgerblue")) +
  labs(
    title = "SD by Eigenvalue Rank",
    subtitle = "Lines: Each k group",
    x = "Eigenvalue Rank",
    y = "Standard Deviation"
  ) +
  theme_minimal(base_size = 15)


# library(ggridges)
# 
# ggplot(sd_df, aes(x = SD, y = Group, fill = Group)) +
#   geom_density_ridges(alpha = 0.7) +
#   labs(
#     title = "SD Distribution Across Groups",
#     x = "Standard Deviation (SD)",
#     y = "k Group"
#   ) +
#   theme_minimal(base_size = 15) +
#   theme(legend.position = "none")



