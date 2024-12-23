library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(SpatialExperiment)
library(qs)
library(standR)
library(scran)
library(FNN)
library(kBET)
wdpath <- "./data/DLBCL_run/"
#read in meta data
mymeta <- read.csv(file.path(wdpath, "GeoMx_count_table","TMAmetadata.csv"))
mymeta
#show excel sheet name
mysegmentproperty_DFCI <- read_excel(file.path(wdpath, "GeoMx_count_table","DLBCL_DFCI.xlsx"),sheet = "SegmentProperties")
mycountmat_DFCI <- read_excel(file.path(wdpath, "GeoMx_count_table","DLBCL_DFCI.xlsx"),sheet = "BioProbeCountMatrix")
#
mysegmentproperty_Rochester <- read_excel(file.path(wdpath, "GeoMx_count_table","DLBCL_Rochester.xlsx"),sheet = "SegmentProperties")
mycountmat_Rochester <- read_excel(file.path(wdpath, "GeoMx_count_table","DLBCL_Rochester.xlsx"),sheet = "BioProbeCountMatrix")
# move negative gene to meta negative to meta data of ROI
mycountmat_clean_DFCI <- mycountmat_DFCI[-c(grep("NegProbe-WTX",mycountmat_DFCI$TargetName),
                                            grep("Custom Negative Set 1",mycountmat_DFCI$TargetName)),] 
mycountmat_clean_Rochester <- mycountmat_Rochester[-c(grep("NegProbe-WTX",mycountmat_DFCI$TargetName),
                                                      grep("Custom Negative Set 1",mycountmat_DFCI$TargetName)),] 
# merge expression matrix
mygeneannotation <- mycountmat_clean_DFCI[,1:12]
mycountmat_clean_DFCI <- mycountmat_clean_DFCI[,c(13:ncol(mycountmat_clean_DFCI))]
mycountmat_clean_Rochester <- mycountmat_clean_Rochester[,c(13:ncol(mycountmat_clean_Rochester))]
mycountmat_merge <- as.matrix(cbind(mycountmat_clean_DFCI,mycountmat_clean_Rochester))
rownames(mycountmat_merge) <- mygeneannotation$TargetName
# merge ROI annotation data
mysegmentproperty_Rochester$ROI <- paste0("Rochester_",mysegmentproperty_Rochester$ROILabel)
table(mysegmentproperty_Rochester$ROI)
mysegmentproperty_DFCI$ROI <- paste0("DFCI_",mysegmentproperty_DFCI$ROILabel)
table(mysegmentproperty_DFCI$ROI)
intersect_colname <- intersect(colnames(mysegmentproperty_Rochester),colnames(mysegmentproperty_DFCI))
setdiff(colnames(mysegmentproperty_Rochester),intersect_colname)
setdiff(colnames(mysegmentproperty_DFCI),intersect_colname)
mysegmentproperty_merge <- rbind.data.frame(mysegmentproperty_Rochester[,intersect_colname],mysegmentproperty_DFCI[,intersect_colname])
mysegmentproperty_merge_clean <- left_join(mysegmentproperty_merge,mymeta,by = "ROI")
mysegmentproperty_merge_clean$EBV_Indicator <- ifelse(mysegmentproperty_merge_clean$EBV == 1, "yes","no")
mysegmentproperty_merge_clean <- as.data.frame(mysegmentproperty_merge_clean)
rownames(mysegmentproperty_merge_clean) <- mysegmentproperty_merge_clean$SegmentDisplayName
keptsample <- mysegmentproperty_merge_clean$SegmentDisplayName[mysegmentproperty_merge_clean$ROI != "Rochester_TonsilA" & mysegmentproperty_merge_clean$ROI != "DFCI_Tonsil1" & mysegmentproperty_merge_clean$ROI != "Rochester_1"& mysegmentproperty_merge_clean$ROI != "Rochester_22"]
table(mysegmentproperty_merge_clean[keptsample,]$ROILabel)
mysegmentproperty_final <- mysegmentproperty_merge_clean[keptsample,]
mycountmat_final <- mycountmat_merge[,keptsample]
dim(mycountmat_final)
# merge tumor 
table(mysegmentproperty_final$ROI)
dim(mysegmentproperty_final)
mysegmentproperty_final$ROI_rename <- gsub("a",".1",mysegmentproperty_final$ROI)
mysegmentproperty_final$ROI_rename <- gsub("b",".2",mysegmentproperty_final$ROI_rename)
table(mysegmentproperty_final$ROI_rename)
table(mysegmentproperty_final$SegmentLabel)
names(table(mysegmentproperty_final$SegmentLabel))
dim(mysegmentproperty_final)

table(mysegmentproperty_final$EBV_Indicator)

# merge Tumor

cell_type_mapping <- c(
  "CD4mem" = "CD4T", 
  "CD4naive" = "CD4T", 
  "CD8mem" = "CD8T", 
  "CD8naive" = "CD8T",
  "M1" = "Macro",
  "M2" = "Macro",
  "Tumor"= "Tumor",
  "TumorBCL2"="Tumor",
  "TumorBCL6" = "Tumor",
  "TumorMyc" = "Tumor",
  "TumorOther" = "Tumor"
)

# Extract the segment labels, ROI labels, and core labels from the column names
extract_labels <- function(colname) {
  parts <- strsplit(colname, " \\| ")[[1]]
  list(segment_label = parts[3], roi_label = parts[2], core_label = parts[1])
}

# Identify unique segment, ROI, and core labels
col_labels <- lapply(colnames(mycountmat_final), extract_labels)
segment_labels <- sapply(col_labels, function(x) x$segment_label)
roi_labels <- sapply(col_labels, function(x) x$roi_label)
core_labels <- sapply(col_labels, function(x) x$core_label)

# Create an empty matrix to store the merged counts
unique_segments <- unique(segment_labels)
unique_rois <- unique(roi_labels)
unique_new_types <- unique(cell_type_mapping)
unique_cores <- unique(core_labels)

# Calculate the number of columns for the merged count matrix
num_cores <- length(unique_cores)
num_rois <- length(unique_rois)
num_new_types <- length(unique_new_types)

# Initialize the merged count matrix
merged_countmat <- matrix(0, nrow = nrow(mycountmat_final), ncol = num_cores * num_rois * num_new_types)
rownames(merged_countmat) <- rownames(mycountmat_final)

# Generate column names for the merged count matrix
col_names <- c()
for (core in unique_cores) {
  for (roi in unique_rois) {
    for (new_type in unique_new_types) {
      col_names <- c(col_names, paste(core, roi, new_type, sep = " | "))
    }
  }
}
colnames(merged_countmat) <- col_names

# Sum the counts for the merged cell types
for (core in unique_cores) {
  for (roi in unique_rois) {
    for (new_type in unique_new_types) {
      # Find the old types that map to the new type
      old_types <- names(cell_type_mapping)[cell_type_mapping == new_type]
      
      # Identify the columns corresponding to the old types, current ROI, and core
      cols_to_merge <- which(segment_labels %in% old_types & roi_labels == roi & core_labels == core)
      
      # Sum the counts for the identified columns
      if (length(cols_to_merge) > 0) {
        merged_countmat[, paste(core, roi, new_type, sep = " | ")] <- rowSums(mycountmat_final[, cols_to_merge, drop = FALSE])
      }
    }
  }
}

# Convert the merged count matrix to a data frame
merged_countmat_df <- as.data.frame(merged_countmat)
merged_countmat_df <- merged_countmat_df[,colSums(merged_countmat_df)!=0]
# Create a mapping from old cell types to new cell types
cell_type_mapping <- c(
  "CD4mem" = "CD4T", 
  "CD4naive" = "CD4T", 
  "CD8mem" = "CD8T", 
  "CD8naive" = "CD8T",
  "M1" = "Macro",
  "M2" = "Macro",
  "Tumor"= "Tumor",
  "TumorBCL2"="Tumor",
  "TumorBCL6" = "Tumor",
  "TumorMyc" = "Tumor",
  "TumorOther" = "Tumor"
)

# Extract the segment labels, ROI labels, and core labels from the row names
extract_labels <- function(rownames) {
  parts <- strsplit(rownames, " \\| ")[[1]]
  list(core_label = parts[1], roi_label = parts[2], segment_label = parts[3])
}

# Extract relevant rows from mysegmentproperty_final
segment_labels <- rownames(mysegmentproperty_final)
col_labels <- lapply(segment_labels, extract_labels)
core_labels <- sapply(col_labels, function(x) x$core_label)
roi_labels <- sapply(col_labels, function(x) x$roi_label)
segment_types <- sapply(col_labels, function(x) x$segment_label)

# Create a new column for the merged cell types in mysegmentproperty_final
mysegmentproperty_final$MergedLabel <- sapply(segment_types, function(x) {
  if (x %in% names(cell_type_mapping)) {
    cell_type_mapping[x]
  } else {
    x
  }
})

# Create a combined label for core, ROI, and merged cell type to identify duplicates
combined_label <- paste(core_labels, roi_labels, mysegmentproperty_final$MergedLabel, sep = " | ")

# Find unique rows based on the combined label
unique_indices <- !duplicated(combined_label)

# Subset the data frame to keep only unique rows
merged_segmentproperty_df <- mysegmentproperty_final[unique_indices, ]

# Ensure the row names are unique
rownames(merged_segmentproperty_df) <- NULL

# Convert the data frame columns to their appropriate types
merged_segmentproperty_df <- as.data.frame(merged_segmentproperty_df, stringsAsFactors = FALSE)

# Ensure the non-numeric columns are characters and numeric columns are numeric
merged_segmentproperty_df[] <- lapply(merged_segmentproperty_df, function(x) {
  if (is.character(x) || is.factor(x)) {
    as.character(x)
  } else if (is.numeric(x)) {
    as.numeric(x)
  } else {
    x
  }
})

# rename 
merged_segmentproperty_df$SegmentDisplayName_new <- paste(merged_segmentproperty_df$ScanLabel, merged_segmentproperty_df$ROILabel, merged_segmentproperty_df$MergedLabel, sep = " | ")
rownames(merged_segmentproperty_df) <- merged_segmentproperty_df$SegmentDisplayName_new
# megre expression with Tumor
mycountmat_final_sub <- mycountmat_final[,c(-grep("Tumor",colnames(mycountmat_final)),
                                            -grep("CD8",colnames(mycountmat_final)),
                                            -grep("CD4",colnames(mycountmat_final)),
                                            -grep("M1",colnames(mycountmat_final)),
                                            -grep("M2",colnames(mycountmat_final)))]
merge_final_GEM <- cbind.data.frame(merged_countmat_df,mycountmat_final_sub)
merge_final_GEM <- merge_final_GEM[,merged_segmentproperty_df$SegmentDisplayName_new]
#
mysegmentproperty_final <- merged_segmentproperty_df
mycountmat_final <- merge_final_GEM
# create spatial Experiment obj
rownames(mycountmat_final)[which(duplicated(rownames(mycountmat_final)))]
mycountmat_clean <- as.matrix(mycountmat_final)

# Convert counts to a data frame for aggregation
df_counts <- data.frame(gene=rownames(mycountmat_clean), 
                        counts=as.data.frame(mycountmat_clean),check.names = F)

# Aggregate by summing counts across duplicated gene names
agg_counts <- aggregate(. ~ gene, data=df_counts, FUN=sum)

# Convert back to matrix
new_counts <- as.matrix(agg_counts[,-1])
rownames(new_counts) <- agg_counts$gene
colnames(new_counts) <- gsub("^counts.","",colnames(new_counts))
# rowname  mygeneannotation
dim(new_counts)
dim(mygeneannotation)
dim(mysegmentproperty_final)
# identical(rownames(mygeneannotation),rownames(new_counts))
spe <- SpatialExperiment(
  assay = list(counts = new_counts),
  rowData = mygeneannotation,
  colData = mysegmentproperty_final)
qsave(spe,file = file.path(wdpath,"GeoMx_count_table","spe_for_DLBCL_mergeTumorTmacro.qs"))
