## set working directory
loadingfunction_wdpath <- c("/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/GeoMXNew_code_data_publish/Code/src/")
library(qs)
library(Seurat)
library(standR)
library(dplyr)
library(SpatialExperiment)
library(ggplot2)
# Create a list of CD4 T cell functional gene signatures
library(readxl)
source(file.path(loadingfunction_wdpath, "1-Pathway_validation_DLBCL_preload_function.R"))
### process codex data
wdpath <- "/bmbl_data/yuzhou/collaborative/Sizun_lab/INDEPTH/SGCC/SGWT_results/DLBCL_GeoMX/GeoMXNew_code_data_publish/Data/"
My_CODEXposition <- read.csv(file.path(wdpath,"SGCC_relevant_analysis","DLBCL_ROIlevel_annotation.csv"))
## remove samples
My_CODEXposition <- My_CODEXposition[My_CODEXposition$coreName != "Rochester_TonsilA" & My_CODEXposition$coreName != "DFCI_Tonsil1",]
# change names
mapping_table <- data.frame(
  SegmentLabel = c("CD4mem", "CD4naive", "CD8mem", "CD8naive", "DC", "Endothelial",
                   "M1", "M2", "Neutrophil", "Other", "Treg", "Tumor", "TumorBCL2",
                   "TumorBCL6", "TumorMyc", "TumorOther"),
  Annotation = c("CD4mem", "CD4naive", "CD8mem", "CD8naive", "DC", "Endothelial",
                 "M1", "M2", "Neutrophil", "Other", "Treg", "Other Tumor", "Tumor BCL2",
                 "Tumor BCL6", "Tumor Myc", "Other Tumor"),
  MergedTumor = c("CD4T", "CD4T", "CD8T", "CD8T", "DC", "Endothelial",
                  "Macro", "Macro", "Neutrophil", "Other", "Treg", "Tumor", "Tumor",
                  "Tumor", "Tumor", "Tumor")
)
My_CODEXposition_merge <- left_join(My_CODEXposition,y = mapping_table, by = "Annotation")
cellnumber_percore <- My_CODEXposition_merge %>%
  # Filter only the relevant cell types
  filter(MergedTumor %in% c("Tumor","CD4T", "CD8T", "Macro", "Endothelial")) %>%
  # Group by coreName and MergedTumor to count the number of occurrences
  group_by(coreName, MergedTumor) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(ROI_CT = paste0(coreName,"_",MergedTumor))
#### CD4 and CD8 are from PMID: 37248301
# tumor genes are from PMID: 35443163
CD4T_panelgene<- read_excel(file.path(wdpath,"GeoMX_batchcorrection","TumorAndTcell_MarkerGene","41591_2023_2371_MOESM3_ESM.xlsx"),sheet = "Table S6")
CD8T_panelgene<- read_excel(file.path(wdpath,"GeoMX_batchcorrection","TumorAndTcell_MarkerGene","41591_2023_2371_MOESM3_ESM.xlsx"),sheet = "Table S4")
Tumor_panelgene <- read_excel(file.path(wdpath,"GeoMX_batchcorrection","TumorAndTcell_MarkerGene","DLBCL-gene signature.xlsx"), sheet = "Top-30 genes for programs")
CD4T_genes_list <- as.list(CD4T_panelgene)
CD4T_genes_list <- lapply(CD4T_genes_list, function(x) x[!is.na(x)])
# Set names for each element in the list according to the column headers, with replacements for special characters or spaces
names(CD4T_genes_list) <- paste0("CD4_",names(CD4T_genes_list))
CD8T_panelgene <- as.list(CD8T_panelgene)
CD8T_panelgene <- lapply(CD8T_panelgene, function(x) x[!is.na(x)])
# Set names for each element in the list according to the column headers, with replacements for special characters or spaces
names(CD8T_panelgene) <- paste0("CD8_",names(CD8T_panelgene))
# Set names for each element in the list according to the column headers, with replacements for special characters or spaces
# A single-cell atlas of diffuse large B cell lymphoma paper for providing Tumor gene list
Tumor_genes_list <- as.list(Tumor_panelgene)
Tumor_genes_list <- lapply(Tumor_genes_list, function(x) {x[!is.na(x)]})
Tumor_genes_list <- lapply(Tumor_genes_list, function(x) {gsub("-",".",x)})

# Define the list with M1 and M2 markers from https://www.sciencedirect.com/science/article/pii/S0092867421000106#mmc2
Macrophage_marker_list <- list(
  M1_marker = c("IL23", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD86", "IL1A", "IL1B", "IL6", "CCL5", "IRF5", "IRF1", "CD40", "IDO1", "KYNU", "CCR7"),
  M2_marker = c("IL4R", "CCL4", "CCL13", "CCL20", "CCL17", "CCL18", "CCL22", "CCL24", "LYVE1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "EGF", "CTSA", "CTSB", "CTSC", "CTSD", "TGFB1", "TGFB2", "TGFB3", "MMP14", "MMP19", "MMP9", "CLEC7A", "WNT7B", "FASL", "TNFSF12", "TNFSF8", "CD276", "VTCN1", "MSR1", "FN1", "IRF4")
)
# Define the list of TAM gene clusters https://www.cell.com/trends/immunology/fulltext/S1471-4906(22)00094-1
TAM_gene_clusters <- list(
  IFN_TAMs = c("CASP1", "CASP4", "CCL2", "CCL3", "CCL4", "CCL7", "CCL8", "CD274", "CD40", "CXCL2", "CXCL3", "CXCL9", "CXCL10", "CXCL11", "IDO1", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", "IRF1", "IRF7", "ISG15", "LAMP3", "PDCD1LG2", "TNFSF10", "C1QA", "C1QC", "CD38", "IL4I1", "IFI44L"),
  Inflam_TAMs = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL20", "CCL3L1", "CCL3L3", "CCL4L2", "CCL4L4", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "G0S2", "IL1B", "IL1RN", "IL6", "INHBA", "KLF2", "KLF6", "NEDD9", "PMAIP1", "S100A8", "S100A9", "SPP1"),
  LA_TAMs = c("ACP5", "AOPE", "APOC1", "ATF1", "C1QA", "C1QB", "C1QC", "CCL18", "CD163", "CD36", "CD63", "CHI3L1", "CTSB", "CTSD", "CTSL", "F13A1", "FABP5", "FOLR2", "GPNMB", "IRF3", "LGALS3", "LIPA", "LPL", "MACRO", "MERTK", "MMP7", "MMP9", "MMP12", "MRC1", "NR1H3", "NRF1", "NUPR1", "PLA2G7", "RNASE1", "SPARC", "SPP1", "TFDP2", "TREM2", "ZEB1"),
  Angio_TAMs = c("ADAM8", "AREG", "BNIP3", "CCL2", "CCL4", "CCL20", "CD163", "CD300E", "CD44", "CD55", "CEBPB", "CLEC5A", "CTSB", "EREG", "FCN1", "FLT1", "FN1", "HES1", "IL1B", "IL1RN", "IL8", "MAF", "MIF", "NR1H3", "OLR1", "PPARG", "S100A8", "S100A9", "S100A12", "SERPINB2", "SLC2A1", "SPIC", "SPP1", "THBS1", "TIMP1", "VCAN", "VEGFA"),
  Reg_TAM = c("CCL2", "CD274", "CD40", "CD80", "CD86", "CHIT1", "CX3CR1", "HLA.A", "HLA.C", "HLA.DQA1", "HLA.DQB1", "HLA.DRA", "HLA.DRB1", "HLA.DRB5", "ICOSLG", "IL10", "ITGA4", "LGALS9", "MACRO", "MRC1", "TGFB2"),
  Prolif_TAM = c("CCNA2", "CDC45", "CDK1", "H2AFC", "HIST1H4C", "HMGB1", "HMGN2", "MKI67", "RRM2", "STMN1", "TOP2A", "TUBA1B", "TUBB", "TYMS")
)
# customized marker for CD8 T, CD4 T and B cell
CD8T_markers <- c("CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "TNF", "FASLG",
                  "CCR7", "CD3D", "CD3E", "CD3G", "KLRG1", "CXCR3", "CXCR6",
                  "EOMES", "TBX21", "LAG3", "PDCD1", "NKG7", "CXCL9", "CXCL10", "CXCL11")
CD4T_markers <- c("CD4", "CD3D", "CD3E", "CD3G",    # General T cell markers
                  "FOXP3", "IL2RA", "CTLA4",        # Treg markers
                  "TBX21", "IFNG", "STAT1",         # Th1 markers
                  "GATA3", "IL4", "IL5", "IL13",    # Th2 markers
                  "RORC", "IL17A", "IL17F",         # Th17 markers
                  "CXCR5", "ICOS", "PDCD1",         # T follicular helper (Tfh) markers
                  "CCR7", "SELL",                   # Naive/memory CD4 T cell markers
                  "STAT3", "IL21", "CCR6",          # Other Th markers
                  "CD40LG", "CD28")
B_cellmarkers <- c("CD19", "CD20", "MS4A1", "CD22", "CD79A", "CD79B",
                   "PAX5", "BLNK", "BCL6", "SPIB", "EBF1", "IRF4",
                   "IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2",
                   "IGHG3", "IGHG4", "IGHE", "IKZF1", "IKZF3",
                   "CR2", "FCRL1", "FCRL2", "FCRL3", "FCRL4", "FCRL5",
                   "CXCR5", "CXCR4", "CD38", "CD40", "CD72",
                   "CD83", "AIM2", "FOXP1", "MYC", "PRDM1")

# read in object
spe_obj <- qread(file.path(wdpath,"GeoMx_count_table","spe_for_DLBCL_mergeTumorTmacro_update.qs"))
dim(spe_obj)
selected_cores <- c("Rochester_4", "Rochester_6",
                    "Rochester_7", "Rochester_9", "Rochester_11", "Rochester_12",
                    "Rochester_13", "Rochester_14",
                    "Rochester_15", "Rochester_16", "Rochester_17", "Rochester_18",
                    "Rochester_19", "Rochester_21", "Rochester_23",
                    "Rochester_25", "Rochester_TonsilA", "DFCI_2.2", "DFCI_3.2",
                    "DFCI_4.1", "DFCI_7.1", "DFCI_8.1",
                    "DFCI_12.1", "DFCI_13.2", "DFCI_14.1", "DFCI_15.2", "DFCI_17.1",
                    "DFCI_18.2", "DFCI_19.2", "DFCI_22.2", "DFCI_23.2", "DFCI_Tonsil1")
# subset object
spe_obj <- spe_obj[, spe_obj$ROI_rename %in% selected_cores]

sample_annotation <- data.frame(MergedLabel = spe_obj$MergedLabel,
                                EBV = spe_obj$EBV_Indicator,
                                Cohort = spe_obj$Cohort,
                                ROI = spe_obj$ROI,
                                Sex = spe_obj$Sex,
                                Source2 = spe_obj$Source2,
                                Source2 = spe_obj$Source)
rownames(sample_annotation) <- colnames(spe_obj)
# subsample wanted samples
sample_annotation_sub <- sample_annotation %>% filter(MergedLabel %in% c("CD8T","CD4T","Macro","Tumor","Endothelial"))
expr_matrix_raw <- assay(spe_obj, "counts")
expr_matrix_raw_sub <- expr_matrix_raw[,rownames(sample_annotation_sub)]

# RUV4 batch correction
# Define the wanted groups and subset the spe object
wanted_groups <- c("CD8T", "CD4T", "Macro", "Endothelial","Tumor")
subset_spe_obj <- spe_obj[, spe_obj$MergedLabel %in% wanted_groups]

# Normalize the data
normalized_spe <- geomxNorm(subset_spe_obj, method = "upperquartile")
# Initialize a list to store the Seurat objects
Upper_seurat_list <- list()
# Loop over the ranges for top_n and k
for (top_n in seq(1000, 5000, by = 1000)) {
  for (k in 1:3) {
    set.seed(123)

    # Find NCGs based on the top_n parameter
    tmp_spe <- findNCGs(normalized_spe, batch_name = "Source2", top_n = top_n)

    # Apply batch correction using the k parameter
    tmp_spe <- geomxBatchCorrection(
      tmp_spe,
      factors = c("MergedLabel", "EBV_Indicator"),
      NCGs = metadata(tmp_spe)$NCGs,
      k = k
    )

    # Extract the adjusted counts
    adjusted_counts <- assay(tmp_spe, i = 2)

    # Create a Seurat object
    matrix.use <- adjusted_counts
    seurat_obj <- CreateSeuratObject(matrix.use, meta.data = as.data.frame(colData(tmp_spe)))

    # Name the Seurat object based on the parameters
    seurat_name <- paste0("Seurat_top", top_n, "_k", k)

    # Store the Seurat object in the list with the corresponding name
    Upper_seurat_list[[seurat_name]] <- seurat_obj
  }
}
qsave(Upper_seurat_list,file = file.path(wdpath,"GeoMX_batchcorrection","upperquartile-Seurat_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated-update.qs"))

#
normalized_spe <- scater::logNormCounts(subset_spe_obj)
LogCPM_seurat_list <- list()
# Loop over the ranges for top_n and k
for (top_n in seq(1000, 5000, by = 1000)) {
  for (k in 1:3) {
    set.seed(123)

    # Find NCGs based on the top_n parameter
    tmp_spe <- findNCGs(normalized_spe, batch_name = "Source2", top_n = top_n)

    # Apply batch correction using the k parameter
    tmp_spe <- geomxBatchCorrection(
      tmp_spe,
      factors = c("MergedLabel", "EBV_Indicator"),
      NCGs = metadata(tmp_spe)$NCGs,
      k = k
    )

    # Extract the adjusted counts
    adjusted_counts <- assay(tmp_spe, i = 2)

    # Create a Seurat object
    matrix.use <- adjusted_counts
    seurat_obj <- CreateSeuratObject(matrix.use, meta.data = as.data.frame(colData(tmp_spe)))

    # Name the Seurat object based on the parameters
    seurat_name <- paste0("Seurat_top", top_n, "_k", k)

    # Store the Seurat object in the list with the corresponding name
    LogCPM_seurat_list[[seurat_name]] <- seurat_obj
  }
}
qsave(LogCPM_seurat_list,file = file.path(wdpath,"GeoMX_batchcorrection","LogCPM-Seurat_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated-update.qs"))

LogCPM_seurat_list <- qread(file.path(wdpath,"GeoMX_batchcorrection","LogCPM-Seurat_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated-update.qs"))
Upper_seurat_list <- qread(file.path(wdpath,"GeoMX_batchcorrection","upperquartile-Seurat_for_DLBCL_includeTumorMergeMacro-batchcorrected-DEGvalidated-update.qs"))
names(LogCPM_seurat_list) <- paste0("LogCPM_",names(LogCPM_seurat_list))
names(Upper_seurat_list) <- paste0("Upper_",names(Upper_seurat_list))
seurat_list <- c(Upper_seurat_list, LogCPM_seurat_list)
# Define the comparisons to test
comparisons <- list(
  c("Macro", "Endothelial"),
  c("CD4T", "Endothelial"),
  c("CD8T", "Endothelial"),
  c("Tumor", "Endothelial")
)
# Apply the benchmarking function to each Seurat object in the list
all_benchmark_results <- list()

for (seurat_name in names(seurat_list)) {
  # Extract the Seurat object from the list
  seurat_obj <- seurat_list[[seurat_name]]
  if(min(seurat_obj@assays$RNA@layers$counts) < 0 ){next()}
  # Run the benchmarking pipeline on the current Seurat object
  benchmark_results <- benchmark_DEA_pipeline(
    seurat_obj = seurat_obj,
    comparisons = comparisons,
    CT_number= cellnumber_percore,
    consider_CTnumberList = c(TRUE, FALSE),
    group_column = "MergedLabel",
    EBV_column = "EBV_Indicator",
    EBV_status = list(c("yes"),
                      c("no"),
                      c("yes", "no")),    # Test both EBV statuses
    confounder_combinations = list(NULL,
                                   c("ruv_W1"),
                                   c("ruv_W1", "ruv_W2"),
                                   c("ruv_W1", "ruv_W2", "ruv_W3")),
    logfc_threshold = 0.05,
    p_val_threshold = 0.01,
    assay = "data"
  )

  # Store the results for this Seurat object
  all_benchmark_results[[seurat_name]] <- benchmark_results
  print(seurat_name)
}
qsave(all_benchmark_results, file = file.path(wdpath, "GeoMX_batchcorrection","Benchmark-includetumorMergeMacroTumorT.qs"))
# The variable `all_benchmark_results` contains the benchmark results for all Seurat objects and parameter combinations
all_benchmark_results <- qread(file.path(wdpath, "GeoMX_batchcorrection","Benchmark-includetumorMergeMacroTumorT.qs"))
#### summarize
# Initialize an empty list to store the flattened results
flattened_mergedbenchmark <- list()

# Iterate over each parameter set in all_benchmark_results
for (param_name in names(all_benchmark_results)) {
  # Iterate over the comparisons in each parameter set
  for (comparison_name in names(all_benchmark_results[[param_name]])) {
    # Create the combined name by concatenating the param_name and comparison_name
    combined_name <- paste(param_name, comparison_name, sep = "_")

    # Store the comparison results in the new list with the combined name
    flattened_mergedbenchmark[[combined_name]] <- all_benchmark_results[[param_name]][[comparison_name]]
  }
}
names(flattened_mergedbenchmark)[1:5]

# The resulting list `cell_type_specific_DEGs` will have DEGs specific for each cell type (CD4T, CD8T, M1, M2)
# Initialize lists to store markers for each cell type
cd4t_markers <- list()
cd8t_markers <- list()
macro_markers <- list()
tumor_markers <- list()
# Iterate through the elements of flattened_mergedbenchmark
for (element_name in names(flattened_mergedbenchmark)) {

  # Extract the cell type from the element name (CD4T, CD8T, M1, M2)
  if (grepl("CD4T_vs_Endothelial", element_name)) {
    cell_type <- "CD4T"
  } else if (grepl("CD8T_vs_Endothelial", element_name)) {
    cell_type <- "CD8T"
  } else if (grepl("Macro_vs_Endothelial", element_name)) {
    cell_type <- "Macro"
    # Tumor_vs_Endothelial， Tumor_vs_CD4T, Tumor_Macro
  } else if (grepl("Tumor_vs_Endothelial", element_name)) {
    cell_type <- "Tumor"
  }
  else {
    next()  # Skip if the cell type is not one of the four
  }

  # Extract the DEG list for the current element
  deg_list <- flattened_mergedbenchmark[[element_name]]$DEG_list

  # Filter DEGs where the positive_in column matches the cell type
  if (!is.null(deg_list)) {
    cell_type_deg <- deg_list[deg_list$positive_in == cell_type, ]

    # Append the filtered DEGs to the respective cell type list
    if (cell_type == "CD4T") {
      cd4t_markers[[element_name]] <- cell_type_deg
    } else if (cell_type == "CD8T") {
      cd8t_markers[[element_name]] <- cell_type_deg
    } else if (cell_type == "Macro") {
      macro_markers[[element_name]] <- cell_type_deg
    } else if (cell_type == "Tumor") {
      tumor_markers[[element_name]] <- cell_type_deg
    }
  }
}
cd4t_markers_results <- hypergeo_test(predicted_list = cd4t_markers,groundtruthlist = unique(unlist(c(CD4T_panelgene,CD4T_markers))))
cd8t_markers_results <- hypergeo_test(predicted_list = cd8t_markers,groundtruthlist = unique(unlist(c(CD8T_panelgene,CD8T_markers))))
macro_markers_results <- hypergeo_test(predicted_list = macro_markers,groundtruthlist = unique(unlist(c(TAM_gene_clusters,unlist(Macrophage_marker_list)))))
tumor_markers_results <- hypergeo_test(predicted_list = tumor_markers,groundtruthlist = union(unique(unlist(Tumor_genes_list)),B_cellmarkers))
# save files
# qsave(cd4t_markers_results, file = file.path(wdpath, "GeoMX_batchcorrection","cd4t_markers_results-includetumorMergeMacro.qs"))
# qsave(cd8t_markers_results, file = file.path(wdpath, "GeoMX_batchcorrection","cd8t_markers_results-includetumorMergeMacro.qs"))
# qsave(macro_markers_results, file = file.path(wdpath, "GeoMX_batchcorrection","macro_markers_results-includetumorMergeMacro.qs"))
# qsave(tumor_markers_results, file = file.path(wdpath, "GeoMX_batchcorrection","tumor_markers_results-includetumorMergeMacro.qs"))

cd4t_markers_results <- qread(file.path(wdpath, "GeoMX_batchcorrection","cd4t_markers_results-includetumorMergeMacro.qs"))
cd8t_markers_results<- qread(file.path(wdpath, "GeoMX_batchcorrection","cd8t_markers_results-includetumorMergeMacro.qs"))
macro_markers_results<- qread(file.path(wdpath, "GeoMX_batchcorrection","macro_markers_results-includetumorMergeMacro.qs"))
tumor_markers_results <- qread(file.path(wdpath, "GeoMX_batchcorrection","tumor_markers_results-includetumorMergeMacro.qs"))

# make supplementary table table
cd4t_markers_results$celltype <- "CD4"
cd8t_markers_results$celltype <- "CD8"
macro_markers_results$celltype <- "Macrophage"
tumor_markers_results$celltype <- "Tumor"


merge_results <- rbind(cd4t_markers_results, cd8t_markers_results, macro_markers_results, tumor_markers_results)
merge_results$normalization <- sapply(strsplit(merge_results$parameter,"_"),"[",1)
merge_results$NumberOfNCG <- sapply(strsplit(merge_results$parameter,"_"),"[",3)
merge_results$k <- sapply(strsplit(merge_results$parameter,"_"),"[",4)
merge_results$DEG_covarianceMat <- sapply(strsplit(merge_results$parameter,"_"), function(x) {
  ruv_index <- grep("^ruv", x) # Identify RUV components
  if (length(ruv_index) == 1) {
    "W1" # Combine RUV-related components
  } else if(length(ruv_index) == 2){
    "W1_W2"
  } else if(length(ruv_index) == 3){
    "W1_W2_W3"
  } else if ("No" %in% x && "Confounders" %in% x) {
    "No_Confounders" # Handle "No_Confounders" case
  } else {
    NA # If no RUV components
  }
})
#LogCPM_Seurat_top5000_k2_CD4T_vs_Endothelial_ruv_W1_ruv_W2_EBV_yes-no_CTnumber_TRUE
merge_results$EBV_status <- ifelse(grepl("_EBV_yes_",merge_results$parameter),"EBV+",
                                   ifelse(grepl("_EBV_no_",merge_results$parameter),"EBV-","both"))
merge_results$cellcount <- ifelse(grepl("CTnumber_TRUE",merge_results$parameter),"TRUE", "FALSE")
# write supplementary table 8
# openxlsx::write.xlsx(merge_results,file = "Supplementary Table 8.xlsx")

