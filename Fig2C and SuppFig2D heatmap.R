library(Seurat)
library(qs)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
### set your working directory
wdpath <- c("./data/Tonsil_run/")
###
L24.seurat <- qread(file.path(wdpath,"GeoMx_count_table/L24_seurat_update.qs"))
L25.seurat <- qread(file.path(wdpath,"GeoMx_count_table/L25_seurat_update.qs"))

signatureName <- c("MBC", 
                  "GCBC", "CD4_CD4vsTregs", "Tregs_tregvscd4", 
                  "CD8 T_OneVsAll", "DC_DCvsAll", 
                  "M2_M2vsM1", "Myeloid_MyeloidvsAll", "Endo","epithelial")

# 
L24heatmap.df <- as.data.frame(L24.seurat@meta.data[,signatureName])
L25heatmap.df <- as.data.frame(L25.seurat@meta.data[,signatureName])

# Specify the column names to calculate mean
cols_to_mean <- signatureName

# Calculate the mean by SegmentLabel
L24_summary_cells <- L24.seurat@meta.data %>%
  group_by(SegmentLabel) %>%
  summarise(across(all_of(cols_to_mean), mean, na.rm = TRUE),
            total_cell_count_CODEX = mean(cell_count_CODEX, na.rm = TRUE))
# remove all cell type ROI
L24_summary_cells <- L24_summary_cells[-8,]
L24_mean_data_mt <- as.matrix(L24_summary_cells[,-c(1,12)])
rownames(L24_mean_data_mt) <- L24_summary_cells$SegmentLabel
# Create a row annotation barplot for total_cell_count_CODEX
row_annot <- rowAnnotation(
  Density = anno_barplot(L24_summary_cells$total_cell_count_CODEX,
                     bar_width = 0.8,
                     border = FALSE,
                     gp = gpar(fill = "skyblue"))
)
col_fun = colorRamp2(c(-1.5, 0, 2.3), c("#FF00FF", "black", "#FFFF00"))

ht <- Heatmap(t(scale(t(L24_mean_data_mt))),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 15),
              row_names_side = "left",
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              right_annotation = row_annot,
              heatmap_legend_param = list(title = "Value"),
              col = col_fun
)

ht
# L25
# Calculate the mean by SegmentLabel
L25_summary_cells <- L25.seurat@meta.data %>%
  group_by(SegmentLabel) %>%
  summarise(across(all_of(cols_to_mean), mean, na.rm = TRUE),
            total_cell_count_CODEX = mean(cell_count_CODEX, na.rm = TRUE))
L25_summary_cells <- L25_summary_cells[-8,]
L25_mean_data_mt <- as.matrix(L25_summary_cells[,-c(1,12)])
rownames(L25_mean_data_mt) <- L25_summary_cells$SegmentLabel
# Create a row annotation barplot for total_cell_count_CODEX
row_annot <- rowAnnotation(
  Density = anno_barplot(L25_summary_cells$total_cell_count_CODEX,
                         bar_width = 0.8,
                         border = FALSE,
                         gp = gpar(fill = "skyblue"))
)

ht <- Heatmap(t(scale(t(L25_mean_data_mt))),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 15),
              row_names_side = "left",
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              right_annotation = row_annot,
              heatmap_legend_param = list(title = "Value"),
              col = col_fun
)

ht

# merge L24 abd L25
L25_summary_cells_rename <- L25_summary_cells
L25_summary_cells_rename$SegmentLabel <- paste0(L25_summary_cells_rename$SegmentLabel,"_L25")
L24_summary_cells_rename <- L24_summary_cells
L24_summary_cells_rename$SegmentLabel <- paste0(L24_summary_cells_rename$SegmentLabel,"_L24")

merge_df <- rbind(L25_summary_cells_rename,L24_summary_cells_rename)
desired_order <- c("BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "DC", "M1", "M2", "Myeloid", "Endo","Other")
merge_df$segment_name <- sub("_(L24|L25)", "", merge_df$SegmentLabel)
merge_df$replicate <- sub("(BCL6nB|BCL6pB|CD4T|CD4Treg|CD8T|DC|Endo|M1|M2|Myeloid|Other)_(L24|L25)", "\\2", merge_df$SegmentLabel)
ordered_merge_df <- merge_df %>%
  arrange(match(segment_name, desired_order), replicate) %>%
  select(-segment_name, -replicate)

merge_mean_data_mt <- as.matrix(ordered_merge_df[,-c(1,12)])
rownames(merge_mean_data_mt) <- ordered_merge_df$SegmentLabel
# define group
# Assuming 'group' is the column in your dataframe with the group information

# Map the groups to the color vector
fill_colors <- c(rep("#f8cec1",2),
                 rep("#ef85b5",2),
                 rep("#4bb04a",2),
                 rep("#bfbbdc",2),
                 rep("#984f9e",2),
                 rep("#357fba",2),
                 rep("#70cddd",2),
                 rep("#d25229",2),
                 rep("#e51e26",2),
                 rep("#a65728",2),
                 rep("#bebebd",2))
# Create a row annotation barplot for total_cell_count_CODEX
row_annot <- rowAnnotation(
  Density = anno_barplot(ordered_merge_df$total_cell_count_CODEX,
                         bar_width = 0.8,
                         border = FALSE,
                         gp = gpar(fill = fill_colors)),width = unit(4, "cm")
)

colnames(merge_mean_data_mt) <- c("MBC geneset","GCBC geneset","CD4 geneset",
                                  "Treg geneset","CD8 geneset","DC geneset",
                                  "Macrophage geneset","Myeloid cell geneset","Endothelial geneset","Epithelial geneset")

# reorder the figure for visualization.
enrichment_order = c("MBC geneset","GCBC geneset","CD4 geneset",
                   "Treg geneset","CD8 geneset","Endothelial geneset","DC geneset",
                   "Macrophage geneset","Myeloid cell geneset","Epithelial geneset")
celltype_order = c("BCL6nB_L24", "BCL6nB_L25", "BCL6pB_L24", "BCL6pB_L25", "CD4T_L24", 
                     "CD4T_L25", "CD4Treg_L24", "CD4Treg_L25", "CD8T_L24", "CD8T_L25", "Endo_L24", "Endo_L25",
                     "DC_L24", "DC_L25", "M1_L24", "M1_L25", "M2_L24", "M2_L25", 
                     "Myeloid_L24", "Myeloid_L25",  "Other_L24", "Other_L25")
# plot heatmap
ht <- Heatmap(t(scale(t(merge_mean_data_mt[celltype_order, enrichment_order]))),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 15),column_names_rot =  45,
              row_names_side = "left",
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              right_annotation = row_annot,
              heatmap_legend_param = list(title = "Value"),
              col = col_fun
)
# save file
pdf(file = "Supplementary Figure 2D.pdf")
print(ht)
dev.off()
##############
# Step 1: Extract cell type and replicate information from SegmentLabel
ordered_merge_df <- ordered_merge_df %>%
  mutate(
    CellType = sub("_L[0-9]+", "", SegmentLabel),
    Replicate = sub(".*_L", "L", SegmentLabel)
  )

# Step 2: Calculate the mean for L24 and L25 replicates
mean_values <- ordered_merge_df %>%
  as.data.frame()%>%
  group_by(CellType) %>% 
  summarise(across(where(is.numeric), mean, .names = "{.col}")) %>%
  arrange(match(CellType, desired_order))
mean_vlaue.rownames <- mean_values$CellType

# define group
# Assuming 'group' is the column in your dataframe with the group information

# Map the groups to the color vector
fill_colors <- c(rep("#f8cec1",1),
                 rep("#ef85b5",1),
                 rep("#4bb04a",1),
                 rep("#bfbbdc",1),
                 rep("#984f9e",1),
                 rep("#357fba",1),
                 rep("#70cddd",1),
                 rep("#d25229",1),
                 rep("#e51e26",1),
                 rep("#a65728",1),
                 rep("#bebebd",1))
# Create a row annotation barplot for total_cell_count_CODEX
row_annot <- rowAnnotation(
  Density = anno_barplot(mean_values$total_cell_count_CODEX,
                         bar_width = 0.8,
                         border = FALSE,
                         gp = gpar(fill = fill_colors)),width = unit(4, "cm")
)

mean_values <- mean_values[,-c(1,12)]
mean_values <- as.matrix(mean_values)
rownames(mean_values) <- mean_vlaue.rownames
#col_fun = colorRamp2(c(-2, 0, 2), c("#557ebb", "#fffbf9", "#d93327"))
colnames(mean_values) <- c("MBC geneset","GCBC geneset","CD4 geneset",
                           "Treg geneset","CD8 geneset","DC geneset",
                           "Macrophage geneset","Myeloid cell geneset","Endothelial geneset","Epithelial geneset")
# reorder the figure for visualization.
enrichment_order = c("MBC geneset","GCBC geneset","CD4 geneset",
                   "Treg geneset","CD8 geneset","Endothelial geneset","DC geneset",
                   "Macrophage geneset","Myeloid cell geneset","Epithelial geneset")
celltype_order = c("BCL6nB", "BCL6pB", "CD4T", "CD4Treg", "CD8T", "Endo","DC", 
                     "M1", "M2", "Myeloid", "Other")
# plot heatmap
ht <- Heatmap(t(scale(t(mean_values[celltype_order, enrichment_order]))),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 15),column_names_rot =  45,
              row_names_side = "left",
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              right_annotation = row_annot,
              heatmap_legend_param = list(title = "Value"),
              col = col_fun
)
# save file
pdf(file = "Figure 2C.pdf")
print(ht)
dev.off()
