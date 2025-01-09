# R package---------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(xlsx)
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

# Read from 10X-----------------------------------------------------------------
bmo23_A.data <- Read10X(data.dir = "/media/xiangmeida/meida2/bmo/SRR18015620/VEGFA/outs/filtered_feature_bc_matrix")
bmo23_A <- CreateSeuratObject(counts = bmo23_A.data, project = "bmo23_A", min.cells = 3, min.features = 200)
bmo23_A

bmo23_AC.data <- Read10X(data.dir = "/media/xiangmeida/meida2/bmo/SRR18015621/VEGFAC/outs/filtered_feature_bc_matrix")
bmo23_AC <- CreateSeuratObject(counts = bmo23_AC.data, project = "bmo23_AC", min.cells = 3, min.features = 200)
bmo23_AC

bmo23_A[["percent.mt"]] <- PercentageFeatureSet(bmo23_A, pattern = "^MT-")
bmo23_AC[["percent.mt"]] <- PercentageFeatureSet(bmo23_AC, pattern = "^MT-")

# Data quality control----------------------------------------------------------
VlnPlot(bmo23_A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(bmo23_AC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

bmo23_A <- subset(bmo23_A, 
                  subset = nFeature_RNA > 300 & nFeature_RNA < 6000 &
                    nCount_RNA > 500 & nCount_RNA < 50000 &
                    percent.mt < 12)

bmo23_AC <- subset(bmo23_AC, 
                   subset = nFeature_RNA > 300 & nFeature_RNA < 6000 &
                     nCount_RNA > 500 & nCount_RNA < 50000 &
                     percent.mt < 12)

# data integration (bmo23_A and bmo23_AC)---------------------------------------
A_AC <- list(bmo23_A, bmo23_AC)
A_AC.1 <- c()
for (i in A_AC){
  i <- NormalizeData(i)
  i <- FindVariableFeatures(i, selection.method = "vst", nfeatures = 2000)
  A_AC.1 <- c(A_AC.1,i)
}

features <- SelectIntegrationFeatures(object.list = A_AC.1, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = A_AC.1, anchor.features = features)
A_AC.combined <- IntegrateData(anchorset = anchors, dims = 1:30)
A_AC.combined[["RNA"]] <- JoinLayers(A_AC.combined[["RNA"]])

DefaultAssay(A_AC.combined) <- "integrated"
A_AC.combined <- ScaleData(A_AC.combined, verbose = FALSE)
A_AC.combined <- RunPCA(A_AC.combined, npcs = 30, verbose = FALSE)
ElbowPlot(A_AC.combined)
A_AC.combined <- RunUMAP(A_AC.combined, reduction = "pca", dims = 1:20)
A_AC.combined <- FindNeighbors(A_AC.combined, reduction = "pca", dims = 1:20)
A_AC.combined <- FindClusters(A_AC.combined, resolution = 0.5)
DefaultAssay(A_AC.combined) <- "RNA"

DimPlot(A_AC.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(A_AC.combined, reduction = "umap",label = T)

save(A_AC.combined,file = "/media/xiangmeida/meida/8-leiqiguan/RData/A_AC.combined.RData")
load("/media/xiangmeida/meida2/bmo/A_AC.combined.RData")

#Cell annotation----------------------------------------------------------------
DefaultAssay(A_AC.combined) <- "RNA"
DimPlot(A_AC.combined, reduction = "umap", group.by = "seurat_clusters",label = T)
#MSCs
FeaturePlot(A_AC.combined, features = c("PDGFRB","PDGFRA","NES",
                                        "KITLG","ITGA4","EMCN",
                                        "CXCL12","CSPG4","COL3A1"),  max.cutoff = 3,  min.cutoff=0, cols = c("grey","red"), reduction = "umap")
#Endothelial cells
FeaturePlot(A_AC.combined, features = c("PECAM1","MCAM","FLT4",
                                        "ENG","CDH5","CD34"),  max.cutoff = 3,  min.cutoff=0, cols = c("grey","red"), reduction = "umap")

#HSPC
FeaturePlot(A_AC.combined, features = c("TRIM6","MYB","CD34"),  max.cutoff = 3,  min.cutoff=0, cols = c("grey","red"), reduction = "umap")

#Myeloid
FeaturePlot(A_AC.combined, features = c("TPSB2","PTPRC","PRSS57",
                                        "ITGAM","CD14"),  max.cutoff = 3,  min.cutoff=0, cols = c("grey","red"), reduction = "umap")

#MK
FeaturePlot(A_AC.combined, features = c("PPBP","PF4","GP9"),  max.cutoff = 3,  min.cutoff=0, cols = c("grey","red"), reduction = "umap")

#Erythroid
FeaturePlot(A_AC.combined, features = c("KLF1","GYPA","ALAS2"),  max.cutoff = 3,  min.cutoff=0, cols = c("grey","red"), reduction = "umap")

A_AC <- A_AC.combined[,A_AC.combined@meta.data$seurat_clusters %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)]
DimPlot(A_AC, reduction = "umap",label = T)

# data processing---------------------------------------------------------------
A_AC <- NormalizeData(A_AC)

A_AC <- FindVariableFeatures(A_AC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(A_AC), 10)
plot1 <- VariableFeaturePlot(A_AC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.A_AC <- rownames(A_AC)
A_AC <- ScaleData(A_AC, features = all.genes.A_AC)

A_AC <- RunPCA(A_AC, features = VariableFeatures(object = A_AC))
ElbowPlot(A_AC)
A_AC <- FindNeighbors(A_AC, dims = 1:20)
A_AC <- FindClusters(A_AC, resolution = 0.5)
head(Idents(A_AC), 5)

A_AC <- RunUMAP(A_AC, dims = 1:20)
DimPlot(A_AC, reduction = "umap",label = TRUE)

A_AC <- RenameIdents(A_AC, `0`="Fibroblasts",`5`="Fibroblasts",
                     `7`="MSCs", `13` ="MSCs",
                     `10` ="Endothelial cells",
                     `6`="Monocytes",`12`="Monocytes",
                     `1`="HSPC",
                     `9`="Monocyte/neutrophil prog.",
                     `8` ="Eosinophil/nasophil/mast prog.",
                     `3` ="Megakaryocytes",
                     `2`="Erythroid",`4`="Erythroid",`11`="Erythroid")
DimPlot(A_AC, label = TRUE)

table(A_AC@active.ident)
A_AC@meta.data$celltype <- ""
A_AC@meta.data$celltype[A_AC@active.ident %in% c("Fibroblasts")] <- "Fibroblasts"
A_AC@meta.data$celltype[A_AC@active.ident %in% c("MSCs")] <- "MSCs"
A_AC@meta.data$celltype[A_AC@active.ident %in% c("Endothelial cells")] <- "Endothelial cells"
A_AC@meta.data$celltype[A_AC@active.ident %in% c("Monocytes")] <- "Monocytes"
A_AC@meta.data$celltype[A_AC@active.ident %in% c("HSPC")] <- "HSPC"
A_AC@meta.data$celltype[A_AC@active.ident %in% c("Monocyte/neutrophil prog.")] <- "Monocyte/neutrophil prog."
A_AC@meta.data$celltype[A_AC@active.ident %in% c("Eosinophil/nasophil/mast prog.")] <- "Eosinophil/nasophil/mast prog."
A_AC@meta.data$celltype[A_AC@active.ident %in% c("Megakaryocytes")] <- "Megakaryocytes"
A_AC@meta.data$celltype[A_AC@active.ident %in% c("Erythroid")] <- "Erythroid"
table(A_AC@meta.data$celltype)

table(A_AC@meta.data[["orig.ident"]])
A_AC@meta.data$group <- ""
A_AC@meta.data$group[A_AC@meta.data[["orig.ident"]] %in% c("VEGFA")] <- "BMO-2023"
A_AC@meta.data$group[A_AC@meta.data[["orig.ident"]] %in% c("VEGFAC")] <- "BMO-2023"
table(A_AC@meta.data$group)

save(A_AC,file = "/media/xiangmeida/meida/8-Organoids/RData/A_AC.RData")

# data integration (bmo23_A , bmo23_AC, RBMOs, RBMOs_Blood)---------------------
load("/media/xiangmeida/meida/8-Organoids/RData/RBMOs.RData")
load("/media/xiangmeida/meida/8-Organoids/RData/RBMOs_Blood.RData")

A_AC_B_RBMO <- list(bmo23_A, bmo23_AC, RBMOs, RBMOs_Blood)
A_AC_B_RBMO.1 <- c()
for (i in A_AC_B_RBMO){
  A_AC_B_RBMO.1 <- c(A_AC_B_RBMO.1,i)
}

features <- SelectIntegrationFeatures(object.list = A_AC_B_RBMO.1, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = A_AC_B_RBMO.1, anchor.features = features)
A_AC_B_RBMO.combined <- IntegrateData(anchorset = anchors, dims = 1:30)
A_AC_B_RBMO.combined[["RNA"]] <- JoinLayers(A_AC_B_RBMO.combined[["RNA"]])

DefaultAssay(A_AC_B_RBMO.combined) <- "integrated"
A_AC_B_RBMO.combined <- ScaleData(A_AC_B_RBMO.combined, verbose = FALSE)
A_AC_B_RBMO.combined <- RunPCA(A_AC_B_RBMO.combined, npcs = 30, verbose = FALSE)
ElbowPlot(A_AC_B_RBMO.combined)
A_AC_B_RBMO.combined <- RunUMAP(A_AC_B_RBMO.combined, reduction = "pca", dims = 1:20)
A_AC_B_RBMO.combined <- FindNeighbors(A_AC_B_RBMO.combined, reduction = "pca", dims = 1:20)
A_AC_B_RBMO.combined <- FindClusters(A_AC_B_RBMO.combined, resolution = 4)

DefaultAssay(A_AC_B_RBMO.combined) <- "RNA"
DimPlot(A_AC_B_RBMO.combined, reduction = "umap", split.by = "orig.ident")
DimPlot(A_AC_B_RBMO.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(A_AC_B_RBMO.combined, reduction = "umap",label = T)

save(A_AC_B_RBMO.combined,file = "/media/xiangmeida/meida/8-8-Organoids/RData/A_AC_B_RBMO.combined.RData")

# Cell annotation---------------------------------------------------------------
source('/media/xiangmeida/meida/8-leiqiguan/画图/Dotplot_anno.R')

markers <- c("POSTN","ASPN",#Fibroblast
             "INSR","FGFR2",#Osteoblast
             "SOX9","HAPLN1","CCN1",#Chondrocyte
             "ACTA2","TAGLN",#Vascular smooth muscle cells
             "THY1","LPL",#Adipo-CAR
             "CCNE2","CDC6","MCM10",#MSC
             "CDH1","CLDN7","EPCAM",#Epithelial cells
             "TFRC","GYPA","SLC4A1",#Erythroid
             "PECAM1","CDH5",#Endothelial cells
             "CD14","MS4A7","MS4A6A",#Monocytes
             "THRA","PTH1R",#Osteoclast
             "ITGA2B","GP9","PPBP",#Megakaryocyte
             "CD68", "CD163",#Macrophage
             "PDGFRB","ANGPT1",#Pericytes
             "CCR3","ITGAM","ADGRE1",#Eo/baso/neutrophil
             "CD33","KIT","ENPP3",#Mast cells
             "CD34","ITGA6")#HSCs

Dotplot_anno(A_AC_B_RBMO.combined, features = markers, #celltype_color = dittoColors(),
             group = c( rep("Fibroblasts",2),rep("Osteoblasts",2),
                        rep("Chondrocytes",3),rep("VSMCs",2),
                        rep("Adipocytes",2),rep("MSC proliferating",3),
                        rep("Epithelial cells",3),rep("Erythroid",3),rep("Endothelial cells",2),
                        rep("Monocytes",3),rep("Osteoclasts",2),
                        rep("Megakaryocytes",3),rep("Macrophages",2),rep("Pericytes",2),
                        rep("Eo/baso/neutrophil",3),rep("Mast cells",3),rep("HSCs",2)),
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             celltype_color=c(
               '#F7E1ED','#E7DAD2',"#E7EFFA","#C497B2","#14517C",'#F3D266',"#82B0D2",
               "#96C37D","#9AC9D8",'#D8383A',"#FF8884","#FFBE7A","#FA7F6F","#2F7FC1",
               "#8ECFC9","#A9B8C6","#F8AC8C"),
             order = T)

A_AC_B_RBMO.combined <- RenameIdents(A_AC_B_RBMO.combined,
                                     `46` = "Erythroid cells",`36` = "Erythroid cells",
                                     `13` = "Erythroid cells",`1` = "Erythroid cells",
                                     `8` = "Erythroid cells",`7` = "Erythroid cells",
                                     `16` = "Erythroid cells",`19` = "Erythroid cells",
                                     `47` = "Erythroid cells",`29` = "Erythroid cells",
                                     `67` = "Erythroid cells",`24` = "Erythroid cells",
                                     `17` = "Erythroid cells", `52` = "Erythroid cells",
                                     `34` = "Erythroid cells", `27` = "Erythroid cells",
                                     `39` = "Erythroid cells",`54` = "Erythroid cells",
                                     `23` = "Erythroid cells",`15` = "Erythroid cells",
                                     `63` = "Erythroid cells",`55` = "Erythroid cells",
                                     `41` = "Erythroid cells",`0` = "Erythroid cells",
                                     `45` = "Erythroid cells",
                                     `3` = "Megakaryocytes", `25` = "Megakaryocytes",
                                     `28` = "Megakaryocytes", `33` = "Megakaryocytes" ,
                                     `14` = "Eo/baso/neutrophil", `66` = "Eo/baso/neutrophil",
                                     `5` = "Eo/baso/neutrophil", 
                                     `2` = "Mast cells" ,`62` = "Mast cells" ,`12` = "Mast cells" ,
                                     `51` = "Macrophages" ,`31` = "Macrophages" ,
                                     `37` = "Monocytes" , `61` = "Monocytes" ,`10` = "Monocytes" ,`65` = "Monocytes" ,
                                     `18` = "HSCs" ,`57` = "HSCs" ,`64` = "HSCs" ,
                                     `50` = "Endothelial cells" ,`35` = "Endothelial cells" ,
                                     `11` = "Osteoclasts" ,`49` = "Osteoclasts" ,`40` = "Osteoclasts" ,
                                     `21` = "Osteoblasts" ,
                                     `43` = "Chondrocytes" ,
                                     `38` = "Fibroblasts" ,`32` = "Fibroblasts" ,
                                     `58` = "Fibroblasts" ,`4` = "Fibroblasts" ,
                                     `9` = "Fibroblasts" ,`9` = "Fibroblasts" ,
                                     `44` = "Fibroblasts" ,`20` = "Fibroblasts" ,
                                     `60` = "Fibroblasts" ,`56` = "Fibroblasts" ,`30` = "Fibroblasts" ,
                                     `59` = "Epithelial cells" ,
                                     `26` = "Adipocytes" ,`22` = "Adipocytes" ,
                                     `6` = "Pericytes" ,
                                     `42` = "Pericytes" ,`48` = "VSMCs" ,
                                     `53` = "MSC proliferating")

DimPlot(A_AC_B_RBMO.combined, reduction = "umap",label = F)







