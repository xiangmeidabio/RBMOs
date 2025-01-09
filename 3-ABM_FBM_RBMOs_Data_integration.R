# R package---------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(devtools)
library(harmony)
library(ggplot2)
library(plyr)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(monocle3)
library(data.table)
library(readxl)
library(data.table)
library(ggpubr)
library(rstatix)
library(grid)
library(pheatmap)
library(scales)
library(org.Mm.eg.db)
library(clusterProfiler)
library(Hmisc)
library(VennDiagram)
library(destiny)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(SeuratDisk)

# RBMOs data--------------------------------------------------------------------
load("/media/xiangmeida/meida/8-Organoids/RDS/RBMOs.RData")
DimPlot(RBMOs, label = TRUE)
table(RBMOs@meta.data$celltype)

# FBM data----------------------------------------------------------------------
FBM <- LoadH5Seurat("/media/xiangmeida/meida/8-Organoids/data/FBM/fig1b_fbm_scaled_gex_updated_dr_20210104.h5seurat")

FBM <- NormalizeData(FBM)
FBM <- FindVariableFeatures(FBM, selection.method = "vst", nfeatures = 2000)
all.genes.FBM <- rownames(FBM)
FBM <- ScaleData(FBM, features = all.genes.FBM)

FBM <- RunPCA(FBM, features = VariableFeatures(object = FBM))
ElbowPlot(FBM)
FBM <- FindNeighbors(FBM, dims = 1:20)
FBM <- FindClusters(FBM, resolution = 0.5)
head(Idents(FBM), 5)
FBM_obj <- RunUMAP(FBM, dims = 1:20)
DimPlot(FBM, reduction = "umap",label = TRUE)

saveRDS(FBM, "/media/xiangmeida/meida/8-Organoids/RDS/FBM.rds")

# ABM data----------------------------------------------------------------------
ABM <- readRDS("/media/xiangmeida/meida/8-Organoids/data/ABM/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds")

ABM <- NormalizeData(ABM)
ABM <- FindVariableFeatures(ABM, selection.method = "vst", nfeatures = 2000)
all.genes.ABM <- rownames(ABM)
ABM <- ScaleData(ABM, features = all.genes.ABM)

ABM <- RunPCA(ABM, features = VariableFeatures(object = ABM))
ElbowPlot(ABM)
ABM <- FindNeighbors(ABM, dims = 1:20)
ABM <- FindClusters(ABM, resolution = 0.5)
head(Idents(ABM), 5)
ABM <- RunUMAP(ABM, dims = 1:20)
DimPlot(ABM, reduction = "umap",label = TRUE)

saveRDS(ABM, "/media/xiangmeida/meida/8-Organoids/RDS/ABM.rds")

#celltype-----------------------------------------------------------------------
#RBMOs
RBMOs@meta.data$celltype <- ""
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Fibroblasts")] <- "RBMOs-Fibroblasts"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Osteoclasts")] <- "RBMOs-Osteoclasts"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Osteoblasts")] <- "RBMOs-Osteoblasts"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Adipocytes")] <- "RBMOs-Adipocytes"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Chondrocytes")] <- "RBMOs-Chondrocytes"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("MSC proliferating")] <- "RBMOs-MSC proliferating"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("VSMCs")] <- "RBMOs-VSMCs"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Epithelial cells")] <- "RBMOs-Epithelial cells"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Erythroid")] <- "RBMOs-Erythroid"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Endothelial cells")] <- "RBMOs-Endothelial cells"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Monocytes")] <- "RBMOs-Monocytes"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Megakaryocytes")] <- "RBMOs-Megakaryocytes"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("HSCs")] <- "RBMOs-HSCs"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Macrophages")] <- "RBMOs-Macrophages"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Pericytes")] <- "RBMOs-Pericytes"
table(RBMOs@meta.data$celltype)

#FBM
FBM@meta.data$orig.ident <- "FBM"
table(FBM@meta.data[["broad_fig1_cell.labels"]])

FBM@meta.data$celltype <- ""
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("HSC/MPP and pro")] <- "FBM-HSC/MPP and pro"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("erythroid")] <- "FBM-erythroid"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("MK")] <- "FBM-MK"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("B_lineage")] <- "FBM-B_lineage"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("DC")] <- "FBM-DC"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("eo/baso/mast")] <- "FBM-eo/baso/mast"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("neutrophil")] <- "FBM-neutrophil"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("monocyte")] <- "FBM-monocyte"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("T_NK")] <- "FBM-T_NK"
FBM@meta.data$celltype[FBM@meta.data[["broad_fig1_cell.labels"]] %in% c("stroma")] <- "FBM-stroma"
table(FBM@meta.data$celltype)

#ABM
ABM@meta.data$orig.ident <- "ABM"

table(ABM@meta.data$celltype)
ABM@meta.data$celltype <- ""
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Adipo-MSC")] <- "ABM-Adipo-MSC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("AEC")] <- "ABM-AEC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("APOD+ MSC")] <- "ABM-APOD+ MSC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Ba/Eo/Ma")] <- "ABM-Ba/Eo/Ma"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("CD4+ T-Cell")] <- "ABM-CD4+ T-Cell"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("CD8+ T-Cell")] <- "ABM-CD8+ T-Cell"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("CLP")] <- "ABM-CLP"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Cycling DCs")] <- "ABM-Cycling DCs"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Cycling HSPC")] <- "ABM-Cycling HSPC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Early Myeloid Progenitor")] <- "ABM-Early Myeloid Progenitor"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Erythroblast")] <- "ABM-Erythroblast"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Fibro-MSC")] <- "ABM-Fibro-MSC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("GMP")] <- "ABM-GMP"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("HSC")] <- "ABM-HSC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Late Erythroid")] <- "ABM-Late Erythroid"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Late Myeloid")] <- "ABM-Late Myeloid"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Macrophages")] <- "ABM-Macrophages"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Mature B")] <- "ABM-Mature B"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Megakaryocyte")] <- "ABM-Megakaryocyte"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("MEP")] <- "ABM-MEP"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Monocyte")] <- "ABM-Monocyte"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("MPP")] <- "ABM-MPP"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Neutrophil")] <- "ABM-Neutrophil"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Osteo-MSC")] <- "ABM-Osteo-MSC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Osteoblast")] <- "ABM-Osteoblast"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("pDC")] <- "ABM-pDC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Plasma Cell")] <- "ABM-Plasma Cell"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Pre-B")] <- "ABM-Pre-B"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Pre-Pro B")] <- "ABM-Pre-Pro B"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("Pro-B")] <- "ABM-Pro-B"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("RBC")] <- "ABM-RBC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("RNAlo MSC")] <- "ABM-RNAlo MSC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("SEC")] <- "ABM-SEC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("THY1+ MSC")] <- "ABM-THY1+ MSC"
ABM@meta.data$celltype[ABM@meta.data[["cluster_anno_l2"]] %in% c("VSMC")] <- "ABM-VSMC"

# Data integration--------------------------------------------------------------
# RBMOs, FBM, ABM
BM <- list(RBMOs, FBM, ABM)
BM.1 <- c()
for (i in BM){
  i <- NormalizeData(i)
  i <- FindVariableFeatures(i, selection.method = "vst", nfeatures = 2000)
  BM.1 <- c(BM.1,i)
}

features <- SelectIntegrationFeatures(object.list = BM.1, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = BM.1, anchor.features = features)
BM.combined <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(BM.combined) <- "integrated"
BM.combined <- ScaleData(BM.combined, verbose = FALSE)
BM.combined <- RunPCA(BM.combined, npcs = 30, verbose = FALSE)
ElbowPlot(BM.combined)
BM.combined <- RunUMAP(BM.combined, reduction = "pca", dims = 1:20)
BM.combined <- FindNeighbors(BM.combined, reduction = "pca", dims = 1:20)
BM.combined <- FindClusters(BM.combined, resolution = 0.5)

DimPlot(BM.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(BM.combined, reduction = "umap")
DimPlot(BM.combined, reduction = "umap", split.by = "orig.ident")

saveRDS(BM.combined, "/media/xiangmeida/meida/8-Organoids/RDS/BM.combined.rds")

# RBMOs, FBM
RBMOs_FBM <- list(RBMOs, FBM)
RBMOs_FBM.1 <- c()
for (i in RBMOs_FBM){
  i <- NormalizeData(i)
  i <- FindVariableFeatures(i, selection.method = "vst", nfeatures = 2000)
  RBMOs_FBM.1 <- c(RBMOs_FBM.1,i)
}

features <- SelectIntegrationFeatures(object.list = RBMOs_FBM.1, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = RBMOs_FBM.1, anchor.features = features)
RBMOs_FBM.combined <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(RBMOs_FBM.combined) <- "integrated"
RBMOs_FBM.combined <- ScaleData(RBMOs_FBM.combined, verbose = FALSE)
RBMOs_FBM.combined <- RunPCA(RBMOs_FBM.combined, npcs = 30, verbose = FALSE)
ElbowPlot(RBMOs_FBM.combined)
RBMOs_FBM.combined <- RunUMAP(RBMOs_FBM.combined, reduction = "pca", dims = 1:20)
RBMOs_FBM.combined <- FindNeighbors(RBMOs_FBM.combined, reduction = "pca", dims = 1:20)
RBMOs_FBM.combined <- FindClusters(RBMOs_FBM.combined, resolution = 0.5)

DimPlot(RBMOs_FBM.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(RBMOs_FBM.combined, reduction = "umap")
DimPlot(RBMOs_FBM.combined, reduction = "umap", split.by = "orig.ident")

saveRDS(RBMOs_FBM.combined, "/media/xiangmeida/meida/8-Organoids/RDS/RBMOs_FBM.combined.rds")

# RBMOs, ABM
RBMOs_ABM <- list(RBMOs, ABM)
RBMOs_ABM.1 <- c()
for (i in RBMOs_ABM){
  i <- NormalizeData(i)
  i <- FindVariableFeatures(i, selection.method = "vst", nfeatures = 2000)
  RBMOs_ABM.1 <- c(RBMOs_ABM.1,i)
}

features <- SelectIntegrationFeatures(object.list = RBMOs_ABM.1, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = RBMOs_ABM.1, anchor.features = features)
RBMOs_ABM.combined <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(RBMOs_ABM.combined) <- "integrated"
RBMOs_ABM.combined <- ScaleData(RBMOs_ABM.combined, verbose = FALSE)
RBMOs_ABM.combined <- RunPCA(RBMOs_ABM.combined, npcs = 30, verbose = FALSE)
ElbowPlot(RBMOs_ABM.combined)
RBMOs_ABM.combined <- RunUMAP(RBMOs_ABM.combined, reduction = "pca", dims = 1:20)
RBMOs_ABM.combined <- FindNeighbors(RBMOs_ABM.combined, reduction = "pca", dims = 1:20)
RBMOs_ABM.combined <- FindClusters(RBMOs_ABM.combined, resolution = 0.5)

DimPlot(RBMOs_ABM.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(RBMOs_ABM.combined, reduction = "umap")
DimPlot(RBMOs_ABM.combined, reduction = "umap", split.by = "orig.ident")

saveRDS(RBMOs_ABM.combined, "/media/xiangmeida/meida/8-Organoids/RDS/RBMOs_ABM.combined.rds")

#correlation------------------------------------------------------------
BM.combined@meta.data$group[BM.combined@meta.data$orig.ident %in% c("RBMOs")] <- "RBMOs"
BM.combined@meta.data$group[BM.combined@meta.data$orig.ident %in% c("FBM")] <- "FBM"
BM.combined@meta.data$group[BM.combined@meta.data$orig.ident %in% c("ABM")] <- "ABM"

#Overall
av <- AverageExpression(BM.combined, group.by = "group", assays = "RNA")
av <- av[[1]]
dim(av)

cg <- names(tail(sort(apply(av,1,sd)),48133))
a <- as.matrix(av[cg,])
exp <- cor(x = a,y = NULL, method = "spearman")

pheatmap::pheatmap(exp, display_numbers = TRUE, number_format = "%.2f",
                   fontsize_number = 10,
                   color = colorRampPalette(c("#FBD37D","white","#4784C1"))(100))

#celltype: RBMOs FBM
av <- AverageExpression(RBMOs_FBM.combined, group.by = "celltype", assays = "RNA")
av <- av[[1]]
dim(av)

cg <- names(tail(sort(apply(av,1,sd)),2000))
a <- as.matrix(av[cg,])
exp <- cor(x = a,y = NULL, method = "spearman")

pheatmap::pheatmap(exp, display_numbers = TRUE, number_format = "%.2f",
                   fontsize_number = 10,
                   color = colorRampPalette(c("navy","white","firebrick3"))(100))

#celltype: RBMOs ABM
av <- AverageExpression(RBMOs_ABM.combined, group.by = "celltype", assays = "RNA")
av <- av[[1]]
dim(av)

cg <- names(tail(sort(apply(av,1,sd)),2000))
a <- as.matrix(av[cg,])
exp <- cor(x = a,y = NULL, method = "spearman")

pheatmap::pheatmap(exp, display_numbers = TRUE, number_format = "%.2f",
                   fontsize_number = 10,
                   color = colorRampPalette(c("navy","white","firebrick3"))(100))



