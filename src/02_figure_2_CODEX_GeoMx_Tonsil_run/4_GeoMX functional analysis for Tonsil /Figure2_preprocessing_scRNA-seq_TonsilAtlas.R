
library(Seurat)
library(qs)
# set up working directory
wdpath <- c("./data/Tonsil_run/")
# read in data
tonsil_atlas_rna_seurat<- readRDS(file.path(wdpath,"Public_scRNASeq/20230911_tonsil_atlas_rna_seurat_obj.rds"))

# only use atlas RNA
# rename M1 and M2
tonsil_atlas_rna_seurat$newannotation <- tonsil_atlas_rna_seurat$annotation_figure_1
m1_name1 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$annotation_20220215 %in% "SELENOP FUCA1 PTGDS macrophages"]
m1_name2 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$annotation_20220215 %in% "C1Q HLA macrophages"]
m1_name3 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$annotation_20220215 %in% "ITGAX ZEB2 macrophages"]
m1_name4 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$annotation_20220215 %in% "IL7R MMP12 macrophages"]
m1_name5 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$annotation_20220215 %in% "M1 Macrophages"]

tonsil_atlas_rna_seurat$newannotation[names(tonsil_atlas_rna_seurat$newannotation) %in% c(m1_name1,m1_name2,m1_name3,m1_name4)] <- "M2 Macrophages"
tonsil_atlas_rna_seurat$newannotation[names(tonsil_atlas_rna_seurat$newannotation) %in% c(m1_name5)] <- "M1 Macrophages"


# rename Treg from CD4
grep("treg", unique(tonsil_atlas_rna_seurat$annotation_20220215),value = T,ignore.case = T)
Treg_name1 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$annotation_20220215 %in% "Eff-Tregs"]
tonsil_atlas_rna_seurat$newannotation[names(tonsil_atlas_rna_seurat$newannotation) %in% c(Treg_name1)] <- "Tregs"


# merge myeloid
grep("myeloid", unique(tonsil_atlas_rna_seurat$annotation_20220215),value = T,ignore.case = T)
tonsil_atlas_rna_seurat$newannotation[tonsil_atlas_rna_seurat$newannotation == "cycling myeloid"] <- "myeloid"
tonsil_atlas_rna_seurat$newannotation[tonsil_atlas_rna_seurat$newannotation == "Mono/Macro"] <- "myeloid"
table(tonsil_atlas_rna_seurat$newannotation)
# drop cell type
drop_name1 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$newannotation %in% "cycling FDC"]
drop_name3 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$newannotation %in% "cycling T"]
drop_name4 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$newannotation %in% "DN"]
drop_name5 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$newannotation %in% "Granulocytes"]
drop_name6 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$newannotation %in% "ILC"]
drop_name7 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$newannotation %in% "Mast"]
drop_name8 <- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$newannotation %in% "NK"]
drop_name9<- colnames(tonsil_atlas_rna_seurat)[tonsil_atlas_rna_seurat$newannotation %in% "preB/T"]

# define new function
`%notin%` <- Negate(`%in%`)
retain_name <- colnames(tonsil_atlas_rna_seurat)[colnames(tonsil_atlas_rna_seurat) %notin% c(drop_name1,drop_name3,
                                                                                             drop_name4,drop_name5,drop_name6,
                                                                                             drop_name7,drop_name8,drop_name9)]

tonsil_atlas_rna_seurat_sub <- subset(tonsil_atlas_rna_seurat,cells = retain_name)
# 
Idents(tonsil_atlas_rna_seurat_sub) <- tonsil_atlas_rna_seurat_sub$newannotation
table(Idents(tonsil_atlas_rna_seurat_sub))
# save object for benchmarking deconvolution
tonsil_atlas_rna_seurat_sub@meta.data <- tonsil_atlas_rna_seurat_sub@meta.data[,c(1:17,41)]
# data is saved as Tonsil460k_reannotated_SeuratObj.qs
# qsave(tonsil_atlas_rna_seurat_sub, file = file.path(wdpath,"Public_scRNASeq/Tonsil460k_reannotated_SeuratObj.qs"))

# calculate DEGs for all cell types
DefaultAssay(tonsil_atlas_rna_seurat_sub) <- "RNA"
# calculate all markers
my.markers <- FindAllMarkers(tonsil_atlas_rna_seurat_sub,slot = "data",only.pos = T)

## subset T cells and find markers
TcellTypes <- c("CD4 T","CD8 T","Naive CD4 T","Naive CD8 T","Tregs")
T_cell_seurat_names<- c()
for (i in 1:length(TcellTypes)){
  T_cell_seurat_names_tmp <- colnames(tonsil_atlas_rna_seurat_sub)[tonsil_atlas_rna_seurat_sub$newannotation %in% TcellTypes[i]]
  T_cell_seurat_names <- c(T_cell_seurat_names,T_cell_seurat_names_tmp)
}
TcellSeurat.obj <- subset(tonsil_atlas_rna_seurat_sub, cells = T_cell_seurat_names)
Idents(TcellSeurat.obj) <- TcellSeurat.obj$newannotation
DefaultAssay(TcellSeurat.obj) <- "RNA"
my.TcellmarkersOneVsAll <- FindAllMarkers(TcellSeurat.obj,slot = "data",only.pos = T)
my.TcellMarkerTregsvsCD4 <- FindMarkers(TcellSeurat.obj,ident.1 = "Tregs",ident.2 = "CD4 T" ,slot = "data", only.pos = T)
my.TcellMarkerCD4vsTregs <- FindMarkers(TcellSeurat.obj,ident.1 = "CD4 T",ident.2 = "Tregs" ,slot = "data", only.pos = T)

# subset myloid
## subset Myeloid cell
MyeloidcellTypes <- c("M1 Macrophages","M2 Macrophages","DC","PDC","myeloid")
Myeloid_seurat_names<- c()
for (i in 1:length(MyeloidcellTypes)){
  Myeloid_seurat_names_tmp <- colnames(tonsil_atlas_rna_seurat_sub)[tonsil_atlas_rna_seurat_sub$newannotation %in% MyeloidcellTypes[i]]
  Myeloid_seurat_names <- c(Myeloid_seurat_names,Myeloid_seurat_names_tmp)
}
MyeloidcellSeurat.obj <- subset(tonsil_atlas_rna_seurat_sub, cells = Myeloid_seurat_names)
Idents(MyeloidcellSeurat.obj) <- MyeloidcellSeurat.obj$newannotation
# set RNA as default assay
DefaultAssay(MyeloidcellSeurat.obj) <- "RNA"
my.MyeloidMarkerOnevsAll <- FindMarkers(MyeloidcellSeurat.obj,ident.1 = "myeloid",ident.2 = c("M1 Macrophages","M2 Macrophages","DC","PDC") ,slot = "data", only.pos = T)
my.DCMarkerOnevsAll <- FindMarkers(MyeloidcellSeurat.obj,ident.1 = "DC",ident.2 = c("M1 Macrophages","M2 Macrophages","PDC","myeloid") ,slot = "data", only.pos = T)
my.M2MarkerOnevsM1 <- FindMarkers(MyeloidcellSeurat.obj,ident.1 = "M2 Macrophages",ident.2 = c("M1 Macrophages") ,slot = "data", only.pos = T)

my.markerlist <- list(all_marker = DEGDF,
                      my.TcellmarkersOneVsAll = my.TcellmarkersOneVsAll, 
                      my.TcellMarkerTregsvsCD4 = my.TcellMarkerTregsvsCD4,
                      my.TcellMarkerCD4vsTregs= my.TcellMarkerCD4vsTregs,
                      my.MyeloidMarkerOnevsAll = my.MyeloidMarkerOnevsAll,
                      my.DCMarkerOnevsAll = my.DCMarkerOnevsAll,
                      my.M2MarkerOnevsM1 = my.M2MarkerOnevsM1)
# save marker for further analyses
# qsave(my.markerlist, file = file.path(wdpath,"Public_scRNASeq/calculated_marker_from_scRNA-seq.qs"))
