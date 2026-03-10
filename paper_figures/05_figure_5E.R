# Define working directory paths
wdpath <- "./data/DLBCL_run/"
loadingfunction_wdpath <- c("./src/")

# Load required libraries
library(Seurat)
library(ggpubr)
library(scGSVA)
library(qs)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(UCell)
library(tidyr)

# Source a script that presumably contains custom functions or gene lists for the analysis
source(file.path(loadingfunction_wdpath, "CODEX_GeoMX analysis for DLBCL", "5-Preload_customized_gene_list.R"))

# Load a Seurat object from a specified path
my.seurat <- readRDS(file.path(wdpath, "CosMX_data", "CosMXSeuratObject.rds"))

# Modify metadata within the Seurat object to classify samples based on field of view (fov) ranges
my.seurat@meta.data  <- my.seurat@meta.data %>% 
  mutate(
    ebv = case_when(
      fov %in% 3:10 ~ "+", 
      fov %in% 11:20 ~ "-"
    ))

# Store the modified metadata for easy access
my.meta <- my.seurat@meta.data

# Set identifiers for cells based on an 'annotation' field within the metadata
Idents(my.seurat) <- my.seurat$annotation

# Set the default assay in the Seurat object to 'RNA'
DefaultAssay(my.seurat) <- "RNA"

# Subset the Seurat object based on cell annotations to create separate objects for specific cell types
my.seurat.tumor <- subset(my.seurat, cells = colnames(my.seurat)[my.seurat$annotation == "Tumor"])
my.seurat.T <- subset(my.seurat, cells = colnames(my.seurat)[my.seurat$annotation == "T"])
my.seurat.Macrophage <- subset(my.seurat, cells = colnames(my.seurat)[my.seurat$annotation == "Macrophage"])

# Initialize an empty data frame for storing results from functional enrichment analysis
result_df <- data.frame(GeneID = character(), PATH = integer(), Annot = character(), stringsAsFactors = FALSE)

# Replace periods with hyphens in gene names within the functional enrichment results
functionalenrichment <- sapply(functionalenrichment, function(x) { gsub("\\.", "-", x) })

# Loop through the functional enrichment list to populate the result_df data frame
for (i in seq_along(functionalenrichment)) {
  # Create a temporary data frame for the current list element
  temp_df <- data.frame(
    GeneID = functionalenrichment[[i]],  # Access the genes in the ith list element
    PATH = i,                            # The index of the list element
    Annot = names(functionalenrichment)[i],  # The name of the list element
    stringsAsFactors = FALSE
  )
  
  # Append the temporary data frame to the result data frame
  result_df <- rbind(result_df, temp_df)
}

# Filter the result_df to include only genes that are present in the Seurat object
Pathway_result_df <- result_df[result_df$GeneID %in% intersect(unique(result_df$GeneID), rownames(my.seurat)),]

# Split the Seurat object by 'fov' to analyze data by specific fields of view
splited.objects <- SplitObject(my.seurat, split.by = "fov")

# Run gene set variation analysis (GSVA) using the 'ssgsea' method on the RNA assay counts
res <- scgsva(my.seurat, slot = "counts", assay = "RNA", Pathway_result_df, method="ssgsea", kcdf = "Poisson")

# Combine GSVA results with existing metadata
my.seurat@meta.data <- cbind.data.frame(my.seurat@meta.data, as.data.frame(res))

# Save the Seurat object with GSVA results
# qsave(my.seurat, file = file.path(wdpath, "CosMX_data", "my.seurat_pathway_calculated_Poisson.qs"))

# Load the saved Seurat object (this step would typically follow the save step)
my.seurat <- qread(file.path(wdpath, "CosMX_data", "my.seurat_pathway_calculated_Poisson.qs"))

# Subset the Seurat object to exclude certain fields of view
my.seurat_subset <- subset(my.seurat, subset = c(fov != 1 & fov != 2))
scmetadata <- my.seurat_subset@meta.data

# Identify pathways of interest based on their column names in metadata
pathway_cols <- names(scmetadata)[grep("^(CD4T_|Macro_|Tumor_)", names(scmetadata))]

# Calculate median values for pathways by annotation, FOV, and EBV status, then organize data
result_df <- scmetadata %>% 
  group_by(annotation, fov, ebv) %>%
  summarise(across(all_of(pathway_cols), median, na.rm = TRUE)) %>%
  ungroup()

# Prepare data for visualization
anntation_df <- result_df[, c("fov", "ebv", "annotation")]

# Long format transformation for T-cell specific pathways
long_df_Tcell <- result_df %>%
  filter(annotation == "T") %>% 
  pivot_longer(
    cols = starts_with("CD4T_"),  # Select T-cell columns for transformation
    names_to = "terminology",     # Define new column name for variable names
    values_to = "value"           # Define new column name for values
  )

# Generate a violin plot with statistical comparisons
T_cell_violin <- ggplot(long_df_Tcell, aes(x = ebv, y = value, fill = ebv)) +
  geom_violin(trim = FALSE, color = "black", width=1.4) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +
  facet_wrap(~ terminology, scales = "free_y") +  # Facet by pathway terminology
  stat_compare_means(
    method = "wilcox.test",  # Wilcoxon test for statistical comparison
    alternative = "both",
    aes(label = ..p.format..),
    label = "p.format",
    comparisons = list(c("+", "-"))  # Compare between EBV positive and negative samples
  ) +
  labs(
    title = "Violin Plot of EBV Status by Macro Terminology with P-values",
    x = "EBV Status",
    y = "Value"
  ) +
  scale_fill_manual(values = c("#d95f02", "#1b9e77"), breaks = c("+", "-")) + 
  geom_jitter(height = 0, width = 0.05) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the violin plot to a PDF file
ggsave(T_cell_violin, device = "pdf", filename = "Tcell_box_scGSVA.pdf", width = 8, height = 8)
