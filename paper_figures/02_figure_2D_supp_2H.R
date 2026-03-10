library(tidyr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(colorRamp2)
library(gridExtra)  # For arranging the plots
library(ggpubr)
# set up working directory
wdpath <- "./data/Tonsil_run/"

#
cybersort_res <- read.csv(file.path(wdpath,"Deconvolution_data","cybersort_result.csv"),row.names = 1)
dtangle_res <- read.csv(file.path(wdpath,"Deconvolution_data","dtangle_result.csv"),row.names = 1)
music_res <- read.csv(file.path(wdpath,"Deconvolution_data","music_result.csv"),row.names = 1)
spatialdecon_res <- read.csv(file.path(wdpath,"Deconvolution_data","spatial_decon_result.csv"),row.names = 1)
# 
groundtruth <- read.delim(file.path(wdpath,"Deconvolution_data","ground_truth.txt"),header = F)
# data for plot pie chart for each ROI
# Extract and clean the cell type and score columns
groundtruth <- groundtruth %>%
  filter(!grepl("ROI", V1)) %>%
  separate(V1, into = c("cell_type", "score"), sep = ": ") %>%
  mutate(score = as.numeric(score)) %>%
  filter(cell_type %in% c("CD4Treg", "CD8T", "DC", "Myeloid", "CD4T", "M2", "MBC", "Tregs","BCL6pB","BCL6nB")) %>%
  group_by(cell_type) %>%
  summarize(average_score = mean(score, na.rm = TRUE))

# Rename cell types to match your data
groundtruth$cell_type <- recode(groundtruth$cell_type,
                                CD4T = "CD4.T",
                                CD8T = "CD8.T",
                                Myeloid = "myeloid",
                                `M2` = "M2.Macrophages",
                                CD4Treg = "Tregs",
                                BCL6pB = "GCBC",
                                BCL6nB = "MBC")
overall_PCC <- read.csv(file.path(wdpath,"Deconvolution_data","indepth_benchmark_result_pcc.csv"))
overall_JSD <- read.csv(file.path(wdpath,"Deconvolution_data","indepth_benchmark_result_JSD.csv"))
overall_RMSE <- read.csv(file.path(wdpath,"Deconvolution_data","indepth_benchmark_result_RMS.csv"))
# Remove ROI_Control
overall_PCC <- overall_PCC %>% filter(X != "ROI_Control")
overall_JSD <- overall_JSD %>% filter(X != "ROI_Control")
overall_RMSE <- overall_RMSE %>% filter(X != "ROI_Control")

# Rename the first column to 'ROI'
colnames(overall_PCC)[1] <- "ROI"
colnames(overall_JSD)[1] <- "ROI"
colnames(overall_RMSE)[1] <- "ROI"

# Add a column to indicate the metric
overall_PCC$metric <- "PCC"
overall_JSD$metric <- "JSD"
overall_RMSE$metric <- "RMSE"
# Melt the data frames to long format
pcc_long <- melt(overall_PCC, id.vars = c("ROI", "metric"), variable.name = "software", value.name = "score")
jsd_long <- melt(overall_JSD, id.vars = c("ROI", "metric"), variable.name = "software", value.name = "score")
rmse_long <- melt(overall_RMSE, id.vars = c("ROI", "metric"), variable.name = "software", value.name = "score")

# Combine all long data frames into one
combined_scores <- bind_rows(pcc_long, jsd_long, rmse_long)
y_range <- range(combined_scores$score, na.rm = TRUE)
combined_scores$metric <- factor(combined_scores$metric, level = c("PCC","JSD","RMSE"))
custom_colors <- c("PCC" = "#42b2d7", "JSD" = "#7ed0de", "RMSE" = "#d3ebed")
# Add a column for software name
cybersort_res$software <- 'cybersort'
dtangle_res$software <- 'dtangle'
music_res$software <- 'music'
spatialdecon_res$software <- 'spatialdecon'
# Melt the data frames to long format
cybersort_long <- melt(cybersort_res, id.vars = 'software', variable.name = 'cell_type', value.name = 'score')
dtangle_long <- melt(dtangle_res, id.vars = 'software', variable.name = 'cell_type', value.name = 'score')
music_long <- melt(music_res, id.vars = 'software', variable.name = 'cell_type', value.name = 'score')
spatialdecon_long <- melt(spatialdecon_res, id.vars = 'software', variable.name = 'cell_type', value.name = 'score')
# Combine all long data frames into one
combined_res <- rbind(cybersort_long, dtangle_long, music_long, spatialdecon_long)
# Define the cell types to subset
cell_types_to_subset <- c("CD4.T", "CD8.T", "DC", "GCBC", "M2.Macrophages", "MBC", "Tregs", "myeloid")
# Set the factor levels for cell types in combined result data and ground truth
cell_types_order <- c("GCBC", "MBC", "Tregs", "CD4.T", "CD8.T", "myeloid", "DC", "M2.Macrophages")
combined_res$cell_type <- factor(combined_res$cell_type, levels = cell_types_order)
groundtruth$cell_type <- factor(groundtruth$cell_type, levels = cell_types_order)
# Subset the data
subset_res <- combined_res %>% filter(cell_type %in% cell_types_to_subset)
#
# Define the cell types to subset and order
subset_res$cell_type <- factor(subset_res$cell_type, levels = c("MBC", "GCBC", "CD4.T", "Tregs", "CD8.T", "Endothelial", "DC", "M2.Macrophages", "myeloid"))
# plot benchmarking data
benchmarking_gg <- ggplot(subset_res, aes(x = software, y = score, fill = cell_type)) +
  # geom_boxplot() +
  geom_bar(position="fill", stat="identity")+
  # geom_hline(data = groundtruth, aes(yintercept = average_score), color = "red", linetype = "dashed") +
 # facet_wrap(~software, nrow = 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("MBC" = "#c76829", 
                               "GCBC" = "#16964a",
                               "CD4.T" = "#2958a8",
                               "Tregs" = "#dd2246",
                               "CD8.T" = "#663287",
                               "myeloid" = "#cf697f",
                               "DC" = "#b8a82c",
                               "M2.Macrophages" = "#b9e4f3"))
# plot ground truth data
groundtruth$cell_type <- factor(groundtruth$cell_type, levels = c("MBC", "GCBC", "CD4.T", "Tregs", "CD8.T", "Endothelial", "DC", "M2.Macrophages", "myeloid"))
gt_gg <- ggplot(groundtruth, aes(x = "cell_type", y = average_score,fill = cell_type)) +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("MBC" = "#c76829", 
                               "GCBC" = "#16964a",
                               "CD4.T" = "#2958a8",
                               "Tregs" = "#dd2246",
                               "CD8.T" = "#663287",
                               "myeloid" = "#cf697f",
                               "DC" = "#b8a82c",
                               "M2.Macrophages" = "#b9e4f3"))
# output 
ggsave("Figrue 2D-benchmarking-toppanel.svg", plot = boxplot_gg, device = "svg", width = 10, height = 6, units = "in")
ggsave("Figrue 2D-groundtruth-toppanel.svg", plot = barplot_gg, device = "svg", width = 10, height = 6, units = "in")


# plot supp figure 2H
# re-read in data
groundtruth <- read.delim(file.path(wdpath,"Deconvolution_data","ground_truth.txt"),header = F)
overall_PCC <- read.csv(file.path(wdpath,"Deconvolution_data","indepth_benchmark_result_pcc.csv"),row.names = 1)
# data for plot pie chart for each ROI
data <- groundtruth %>%
  mutate(ROI = cumsum(grepl("ROI:", V1)),  # Increment ROI number each time "ROI:" is found
         Entry = ifelse(grepl("ROI:", V1), NA, V1)) %>%  # Keep cell type and value string only when not "ROI:"
  filter(!is.na(Entry)) %>%  # Remove rows with NA in Entry
  separate(Entry, into = c("Type", "Value"), sep = ": ", convert = TRUE)
# Rename cell types to match your data
data$Type <- recode(data$Type,
                                CD4T = "CD4.T",
                                CD8T = "CD8.T",
                                Myeloid = "myeloid",
                                `M2` = "M2.Macrophages",
                                CD4Treg = "Tregs",
                                BCL6pB = "GCBC",
                                BCL6nB = "MBC")
# rank by entropy
entropy_scores <- data %>%
  group_by(ROI) %>%
  # Calculate the proportion of each cell type within the group
  mutate(Proportion = Value / sum(Value, na.rm = TRUE)) %>%
  # Calculate the entropy component for each row
  mutate(EntropyComponent = -Proportion * log(Proportion)) %>%
  # Sum up entropy components to get total entropy for each ROI
  summarise(Entropy = sum(EntropyComponent, na.rm = TRUE))
# gini index
gini_simpson_by_roi <- data %>%
  group_by(ROI) %>%
  summarise(Gini_Simpson = 1 - sum(Value^2))
entropy_scores$gini_simpson_by_roi <- gini_simpson_by_roi$Gini_Simpson
entropy_scores <- entropy_scores[entropy_scores$ROI!="17",]
entropy_scores$ROI <- sprintf("ROI_%03d", as.numeric(entropy_scores$ROI)-1)
# remove all. ROI
overall_PCC <- overall_PCC[-17,]
# Join entropy_scores with PCC data to sort them
sorted_data <- entropy_scores[order(entropy_scores$gini_simpson_by_roi),]
# Reorder rows of the PCC matrix based on sorted entropy
pcc_matrix <- overall_PCC[sorted_data$ROI, ]
# Prepare the annotation data
entropy_annotation <- rowAnnotation(df = data.frame(Entropy = sorted_data$gini_simpson_by_roi),
                                        col = list(Entropy = colorRamp2(c(min(sorted_data$gini_simpson_by_roi), max(sorted_data$gini_simpson_by_roi)), 
                                                                        c("#dcdcdc", "black"))))
# Draw the heatmap
p <- Heatmap(as.matrix(pcc_matrix[,c("SpatialDecon","Cybersort","MuSic","dtangle")]), 
        name = "PCC",
        left_annotation  = entropy_annotation,
        col = colorRampPalette(c("white","#5fa1d3", "#083776"))(100),
        na_col = "grey",  # Color for NA values
        show_row_names = TRUE,cluster_columns = F,cluster_rows = F,
        show_column_names = TRUE)

svg("Supplementary Figure 2H.svg", width = 8, height = 6)  # Set the size as needed
draw(p)
dev.off()  # Close the SVG device

### plot Figure 2D bottom panel
# Reshape the data into long format
pcc_df <- pcc_matrix
pcc_df$ROI <- rownames(pcc_matrix)
pcc_matrix_long <- pcc_df %>%
  pivot_longer(cols = -ROI, names_to = "Method", values_to = "PCC")

compare_means <- compare_means(PCC ~ Method, data = pcc_matrix_long, method = "wilcox.test",
                               comparisons = list(c("MuSic", "SpatialDecon"), 
                                                  c("MuSic", "Cybersort"), 
                                                  c("SpatialDecon", "Cybersort")))
pcc_matrix_long$Method <- factor(pcc_matrix_long$Method,levels = c("Cybersort","dtangle", "MuSic", "SpatialDecon") )
# Create the boxplot with significance annotations
p <- ggplot(pcc_matrix_long, aes(x = Method, y = PCC, fill = Method)) +
  geom_violin() +
  geom_boxplot(width=0.1)+
  stat_compare_means(comparisons = list(c("MuSic", "SpatialDecon"), 
                                        c("MuSic", "Cybersort"), 
                                        c("SpatialDecon", "Cybersort"))) +
  labs(title = "PCC Values by Method",
       x = "Method",
       y = "PCC") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = p,filename= "Figrue 2D-benchmarking-bottompanel.svg",device = "svg", width = 8, height = 6)  # Set the size as needed

