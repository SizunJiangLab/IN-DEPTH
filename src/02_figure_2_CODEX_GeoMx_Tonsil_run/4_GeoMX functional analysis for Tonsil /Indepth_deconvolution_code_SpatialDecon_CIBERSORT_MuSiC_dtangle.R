################################################################################################
#
# Merge L24 and L25 data
#
################################################################################################

# h5ad to Seurat object
# The following code can't works because of /layers problem. So I use anndata of python to get the expression matrix and save them in csv
# import anndata
# import pandas as pd
# import numpy as np
# from scipy.sparse import csr_matrix
# adata=anndata.read("L25_clean_GeoMX_rawcount.h5ad")
# adata.layers['log_transformed'] = adata.X
# data_new = adata.to_df(layer='log_transformed')
# data_new2 = data_new.T
# data_new2.to_csv("L25_clean_GeoMX_rawcount.csv")

library("Biobase")
library(qs)
library("Seurat")
library("TOAST")
library("SingleCellExperiment")

data_L24<-qread("L24_clean_GeoMX_rawcount_seurat.qs")
data_L250 <- read.csv("L25_clean_GeoMX_rawcount.csv")
rownames(data_L250) <- data_L250[,1]
data_L250 <- data_L250[,-1]
data_L250 <- as(as.matrix(data_L250), "dgCMatrix")
options(Seurat.object.assay.version = "v3")
data_L25 <- CreateSeuratObject(counts = data_L250)
#data_L25<-qread("L25_clean_GeoMX_rawcount.h5ad")
merged_24_25<-merge(x=data_L24, y=data_L25)
merged_24_25_mtx<-merged_24_25@assays$RNA@counts
write.csv(merged_24_25_mtx,"merged_24_25_raw.csv")


################################################################################################
#
#select celltype and marker
#
################################################################################################

file_path <- "Tonsil460k_reannotated_SeuratObj.qs"
sc_data <- qread(file_path)
sce.tmp = sc_data[, sample(1:ncol(sc_data),round(ncol(sc_data)/20)) ]
sce.tmp2 <- sce.tmp
sce.tmp2@meta.data[which(is.na(sce.tmp2@meta.data[,'newannotation'])),'newannotation'] = "other"
tmp_list <- c("M1 Macrophages","DC","myeloid","CD8 T","CD4 T","MBC","GCBC","Tregs")
tmp_list2 <- rownames(sce.tmp2@meta.data)[which(sce.tmp2@meta.data[,'newannotation'] == "M2 Macrophages")]
for (cell_type in tmp_list)
{
  tmp_list2 <- c(tmp_list2,rownames(sce.tmp2@meta.data)[which(sce.tmp2@meta.data[,'newannotation'] == cell_type)])
}
sce.tmp2@meta.data <- sce.tmp2@meta.data[tmp_list2,]
tmp_list3 <- read.table("deg_merker.txt", head=FALSE)
tmp_list3 <- list(tmp_list3[,1])
sce_gene<-rownames(sce.tmp2@assays$RNA@counts)
commongenes <- intersect (tmp_list3[[1]], sce_gene)
sce.tmp2@assays$RNA@counts <- sce.tmp2@assays$RNA@counts[commongenes,tmp_list2]
sce.tmp2@assays$RNA@meta.features <- sce.tmp2@assays$RNA@meta.features[commongenes,]
saveRDS(sce.tmp2, file = "simple_sampled_data2.rds")

################################################################################################
#
# SpatialDecon
#
################################################################################################

# Install the packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# install.packages("devtools")
# install.packages('e1071')
# install.packages('parallel')
# library(BiocManager)
# BiocManager::install('preprocessCore')
# # install.packages("pheatmap")
# # BiocManager::install('ComplexHeatmap')
# BiocManager::install("scRNAseq")
# devtools::install_github("Nanostring-Biostats/SpatialDecon")


library(e1071)
library(preprocessCore)
library(parallel)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(devtools)
library("Seurat")
library(SpatialDecon)

# Load data
data_merge <- read.csv("roi_merge_data_raw_clean.csv")
rownames(data_merge)=data_merge[,1]
data_merge=data_merge[,-1]
data_merge <- as(as.matrix(data_merge), "dgCMatrix")
options(Seurat.object.assay.version = "v3")
data_merge <- CreateSeuratObject(counts = data_merge)
data_merge <- NormalizeData(data_merge, normalization.method = "LogNormalize")
our_raw <- data_merge@assays$RNA@counts
our_norm <- data_merge@assays$RNA@data
our_norm <- as.matrix(our_norm)
our_raw <- as.matrix(our_raw)
our_norm<-rbind(our_norm,rep(1,len = length(colnames(our_norm))))#Add NegProb to our_norm
row.names(our_norm)[18678] <- "NegProbe"
per.observation.mean.neg = our_norm["NegProbe", ]
bg = sweep(our_norm * 0, 2, per.observation.mean.neg, "+")
bg2 = derive_GeoMx_background(norm = our_norm,
                              probepool = rep(1, nrow(our_norm)),
                              negnames = "NegProbe")
sampled_sc<-readRDS("simple_sampled_data2.rds")
Idents(sampled_sc) <- sampled_sc@meta.data$newannotation
annots <- data.frame(cbind(cellType=as.character(Idents(sampled_sc)), 
                           cellID=names(Idents(sampled_sc))))

# Get mixture reference
custom_mtx_seurat <- create_profile_matrix(mtx = Seurat::GetAssayData(object = sampled_sc, 
                                                                      assay = "RNA", 
                                                                      slot = "counts"), 
                                           cellAnnots = annots, 
                                           cellTypeCol = "cellType", 
                                           cellNameCol = "cellID", 
                                           matrixName = "custom_mini_colon",
                                           outDir = NULL, 
                                           normalize = FALSE, 
                                           minCellNum = 5, 
                                           minGenes = 10)

# Do SpatialDecon
res = spatialdecon(norm = our_norm, bg = bg, X = custom_mtx_seurat, align_genes = TRUE)
result_mtx = res$beta
result_mtx = t(result_mtx)

# Output result
write.csv(result_mtx,"spatial_decon_result.csv")

################################################################################################
#
# CIBERSORT
#
################################################################################################

# Install the packages
# if(!require(CIBERSORT))devtools::install_github("Moonerss/CIBERSORT")

library(CIBERSORT)


# Load data
data_merge <- read.csv("roi_merge_data_raw_clean.csv")

# Do clustering on sc_RNA data first
sampled_sc<-readRDS("simple_sampled_data2.rds")
counts <- GetAssayData(sampled_sc, slot = "counts")
metadata <- sampled_sc@meta.data
gene_annotation <- sampled_sc@assays$RNA@meta.features
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata,
  rowData = gene_annotation
)

# Do Cibersort
rownames(data_merge)=data_merge[,1]
data_merge=data_merge[,-1]
results <- cibersort(sig_matrix = custom_mtx_seurat, mixture_file = data_merge,QN=FALSE)
cibersort()
rowscale <- results[,1:ncol(custom_mtx_seurat)]#back up results

# Output result
write.csv(rowscale,"cybersort_result.csv")


################################################################################################
#
# MuSiC
#
################################################################################################

# Install the packages
# install.packages("tidyr")
# library("tidyr")
# library(qs)
# install.packages("TOAST_1.18.0.tar.gz", repos = NULL, type="source")
# BiocManager::install("TOAST")
# BiocManager::install("SingleCellExperiment")
# library("TOAST")
# library("SingleCellExperiment")
# devtools::install_github('xuranw/MuSiC')
# BiocManager::install("Biobase")

library(MuSiC)
library("Biobase")

# Load data
sampled_sc<-readRDS("simple_sampled_data2.rds")
counts <- GetAssayData(sampled_sc, slot = "counts")
metadata <- sampled_sc@meta.data
gene_annotation <- sampled_sc@assays$RNA@meta.features
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata,
  rowData = gene_annotation
)
data_merge <- read.csv("roi_merge_data_raw_clean.csv")
rownames(data_merge)=data_merge[,1]
data_merge=data_merge[,-1]
data_merge <- as.matrix(data_merge)

# Do MuSiC
Est.prop_24_25 = music_prop(bulk.mtx = data_merge, sc.sce = sce, 
                            clusters = 'newannotation', samples = 'gem_id',verbose = F)
prop_mat_24_25 <- Est.prop_24_25$Est.prop.weighted
rowscale <- prop_mat_24_25

# Output result
write.csv(rowscale,"music_result.csv")


################################################################################################
#
# dtangle
#
################################################################################################

# Install the packages
# BiocManager::install("GEOquery")
# install.packages("dtangle")

library(e1071)
library(preprocessCore)
library(parallel)
library(ggplot2)
library(pheatmap)
library(CIBERSORT)
library(scRNAseq)
library(GEOquery)
library(dtangle)
library(qs)
library("TOAST")
library("Seurat")
library("SingleCellExperiment")
library("Biobase")

# Load data
data_merge <- read.csv("roi_merge_data_raw_clean.csv")
rownames(data_merge)=data_merge[,1]
data_merge=data_merge[,-1]
data_merge <- as.matrix(data_merge)
data_merge <- as(data_merge, "dgCMatrix")
options(Seurat.object.assay.version = "v3")
data_merge <- CreateSeuratObject(counts = data_merge)
data_merge <- NormalizeData(data_merge, normalization.method = "LogNormalize")
me <- as.matrix(data_merge@assays$RNA@data)
sampled_sc<-readRDS("simple_sampled_data2.rds")
counts <- GetAssayData(sampled_sc, slot = "counts")
metadata <- sampled_sc@meta.data
gene_annotation <- sampled_sc@assays$RNA@meta.features
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata,
  rowData = gene_annotation
)


# Get the common gene an normalize
sce.mtx <- assay(sce,"counts")
sce.mtx <- as.matrix(sce.mtx)
sce_gene<-rownames(sce.mtx)
me_gene<-rownames(me)
commongenes <- intersect (me_gene, sce_gene)
me <- me[pmatch(commongenes, rownames(me)), ]
sce.mtx <- sce.mtx[pmatch(commongenes, rownames(sce)), ]
y <- cbind(sce.mtx, me)
y <- normalizeBetweenArrays(y)
y <- t(y)

# Record sample of each cell type
sce@colData <- sce@colData[!is.na(sce@colData$newannotation), ]
all_cell_type <- unique(sce@colData$newannotation)
pure_samples <- lapply(1:length(all_cell_type), function(i) {
  which(sce@colData$newannotation == all_cell_type[i])
})
names(pure_samples) = all_cell_type

# Choose the marker needed to use
marker_list = find_markers(y,pure_samples=pure_samples,data_type="rna-seq",marker_method='ratio')

q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_markers = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})
n_markers

# Do dtangle
marks = marker_list$L
dc <- dtangle(y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'microarray-gene', markers = marks)

# Get the proportion of the cell
final_est <- dc$estimates[(dim(sce)[2]+1):dim(y)[1],]
colnames(final_est) <-  all_cell_type
head(final_est)

# Output result
rowscale <- final_est
write.csv(rowscale,"dtangle_result.csv")



