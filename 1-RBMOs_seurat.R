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
RBMOs.data <- Read10X(data.dir = "/media/xiangmeida/meida/8-Organoids/data/RBMOs")
RBMOs <- CreateSeuratObject(counts = RBMOs.data, project = "RBMOs", min.cells = 3, min.features = 200)
RBMOs[["percent.mt"]] <- PercentageFeatureSet(RBMOs, pattern = "^MT-")

# Data quality control----------------------------------------------------------
VlnPlot(RBMOs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RBMOs <- subset(RBMOs, subset =  nFeature_RNA < 7500 & nFeature_RNA > 400 & percent.mt < 15)

# data processing---------------------------------------------------------------
RBMOs <- NormalizeData(RBMOs)

RBMOs <- FindVariableFeatures(RBMOs, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(RBMOs), 10)
plot1 <- VariableFeaturePlot(RBMOs)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.RBMOs <- rownames(RBMOs)
RBMOs <- ScaleData(RBMOs, features = all.genes.RBMOs)

RBMOs <- RunPCA(RBMOs, features = VariableFeatures(object = RBMOs))

ElbowPlot(RBMOs)

RBMOs <- FindNeighbors(RBMOs, dims = 1:10)
RBMOs <- FindClusters(RBMOs, resolution = 10)
head(Idents(RBMOs), 5)
RBMOs <- RunUMAP(RBMOs, dims = 1:10)
DimPlot(RBMOs, reduction = "umap",label = TRUE)

# Cell annotation---------------------------------------------------------------
RBMOs <- RenameIdents(RBMOs, `17` = "Fibroblasts", `44` = "Fibroblasts",`23` = "Fibroblasts",`28` = "Fibroblasts",
                      `40` = "Fibroblasts",`85` = "Fibroblasts",`26` = "Fibroblasts",`32` = "Fibroblasts",
                      `21` = "Fibroblasts",`9` = "Fibroblasts",`37` = "Fibroblasts",`35` = "Fibroblasts",
                      `4` = "Fibroblasts",`62` = "Fibroblasts",`11` = "Fibroblasts",`66` = "Fibroblasts",
                      `30` = "Fibroblasts",`65` = "Fibroblasts",`56` = "Fibroblasts",
                      `55` = "Fibroblasts",`53` = "Fibroblasts",`19` = "Fibroblasts",`47` = "Fibroblasts",`24` = "Fibroblasts",
                      `1` = "Fibroblasts",`39` = "Fibroblasts",`42` = "Fibroblasts",`46` = "Fibroblasts",`5` = "Fibroblasts",`61` = "Fibroblasts",
                      
                      `12` = "Osteoclasts",`0` = "Osteoclasts",`10` = "Osteoclasts",`49` = "Osteoclasts",`22` = "Osteoclasts",
                      `86` = "Osteoclasts",`77` = "Osteoclasts",`7` = "Osteoclasts",`34` = "Osteoclasts",
                      
                      `15` = "Osteoblasts",`78` = "Osteoblasts",`51` = "Osteoblasts",`14` = "Osteoblasts",`54` = "Osteoblasts",
                      `18` = "Osteoblasts",`16` = "Osteoblasts",`69` = "Osteoblasts",`6` = "Osteoblasts",
                      
                      `82` = "Adipocytes",`33` = "Adipocytes",`13` = "Adipocytes",
                      
                      `71` = "Pericytes",`29` = "Pericytes",`38` = "Pericytes",
                      `70` = "Pericytes",`27` = "Pericytes",`74` = "Pericytes",`41` = "Pericytes",
                      
                      `43` = "Chondrocytes",`67` = "Chondrocytes",`84` = "Chondrocytes",`45` = "Chondrocytes",`31` = "Chondrocytes",
                      `36` = "Chondrocytes",`8` = "Chondrocytes",
                      `64` = "MSC proliferating",`58` = "MSC proliferating",`25` = "MSC proliferating",`89` = "MSC proliferating",`3` = "MSC proliferating",
                      `79` = "MSC proliferating",`20` = "MSC proliferating",
                      `59` = "VSMCs",`83` = "VSMCs",`88` = "VSMCs",`48` = "VSMCs",`2` = "VSMCs",
                      `76` = "Epithelial cells",`52` = "Epithelial cells",`60` = "Epithelial cells",`68` = "Epithelial cells",`63` = "Epithelial cells",
                      `57` = "Erythroid",`73` = "Erythroid",`75` = "Erythroid",
                      `50` = "Endothelial cells",
                      `72` = "Monocytes",`87` = "Macrophages",
                      `80` = "Megakaryocytes",`81` = "Megakaryocytes",
                      `90` = "HSCs")

DimPlot(RBMOs, label = TRUE)

source('/media/xiangmeida/meida/8-Organoids/Draw/Dotplot_anno.R')
markers <- c("POSTN","ASPN",#Fibroblast
             "INSR","FGFR2",#Osteoblast
             "SOX9","HAPLN1","CCN1",#Chondrocyte
             "ACTA2","TAGLN",#Vascular smooth muscle cells
             "THY1","LPL",#Adipo-CAR
             "CCNE2","CDC6","MCM10",#MSC
             "CDH1","CLDN7","EPCAM",#Epithelial cells
             "HBG1","HBA1","HBA2",#Erythroid
             "PECAM1","CDH5",#Endothelial cells
             "CD14","MS4A7","MS4A6A",#Monocytes
             "THRA","PTH1R",#Osteoclast
             "ITGA2B","GP9","PPBP",#Megakaryocyte
             "CD68", "CD163",#Macrophage
             "PDGFRB","ANGPT1",#Pericytes
             "CD34","ITGA6")#HSC

Dotplot_anno(RBMOs, features = markers,
             group = c( rep("Fibroblasts",2),rep("Osteoblasts",2),
                        rep("Chondrocytes",3),rep("VSMCs",2),
                        rep("Adipocytes",2),rep("MSC proliferating",3),
                        rep("Epithelial cells",3),rep("Erythroid",3),rep("Endothelial cells",2),
                        rep("Monocytes",3),rep("Osteoclasts",2),
                        rep("Megakaryocytes",3),rep("Macrophages",2),rep("Pericytes",2),rep("HSCs",2)),
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             celltype_color=c(
                              "#fb8d62","#FB9A99","#818181","#1F78B4",
                              "#E31A1C","#FF7F00","#9E0142","#EE1289",
                              "#B2DF8A", "#A6CEE3","#FDBF6F","#6A3D9A",
                              "#CAB2D6","#B15928","#33A02C"),
             order = T)

RBMOs@meta.data$celltype <- ""
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Fibroblasts")] <- "Fibroblasts"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Osteoclasts")] <- "Osteoclasts"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Osteoblasts")] <- "Osteoblasts"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Adipocytes")] <- "Adipocytes"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Chondrocytes")] <- "Chondrocytes"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("MSC proliferating")] <- "MSC proliferating"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("VSMCs")] <- "VSMCs"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Epithelial cells")] <- "Epithelial cells"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Erythroid")] <- "Erythroid"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Endothelial cells")] <- "Endothelial cells"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Monocytes")] <- "Monocytes"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Megakaryocytes")] <- "Megakaryocytes"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("HSCs")] <- "HSCs"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Macrophages")] <- "Macrophages"
RBMOs@meta.data$celltype[RBMOs@active.ident %in% c("Pericytes")] <- "Pericytes"
table(RBMOs@meta.data$celltype)

save(RBMOs,file = "/media/xiangmeida/meida/8-Organoids/RData/RBMOs.RData")
saveRDS(RBMOs, "/media/xiangmeida/meida/8-Organoids/RDS/RBMOs.rds")

# erythrocytes in RBMOs---------------------------------------------------------
RBMOs_erythrocytes = RBMOs[,RBMOs@meta.data$seurat_clusters %in% c(57,73,75)]
DimPlot(RBMOs_erythrocytes, reduction = "umap")

RBMOs_erythrocytes <- NormalizeData(RBMOs_erythrocytes)

RBMOs_erythrocytes <- FindVariableFeatures(RBMOs_erythrocytes, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(RBMOs_erythrocytes), 10)
plot1 <- VariableFeaturePlot(RBMOs_erythrocytes)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.RBMOs_erythrocytes <- rownames(RBMOs_erythrocytes)
RBMOs_erythrocytes <- ScaleData(RBMOs_erythrocytes, features = all.genes.RBMOs_erythrocytes)

RBMOs_erythrocytes <- RunPCA(RBMOs_erythrocytes, features = VariableFeatures(object = RBMOs_erythrocytes))
ElbowPlot(RBMOs_erythrocytes)

RBMOs_erythrocytes <- FindNeighbors(RBMOs_erythrocytes, dims = 1:10)
RBMOs_erythrocytes <- FindClusters(RBMOs_erythrocytes, resolution = 0.5)
head(Idents(RBMOs_erythrocytes), 5)

RBMOs_erythrocytes <- RunUMAP(RBMOs_erythrocytes, dims = 1:15)
DimPlot(RBMOs_erythrocytes, reduction = "umap",label = TRUE)


markers <- c("CD36","PTPRC",#Early erythroid
             "SOX6","CD44","SLC4A1",#Mid erythroid
             "GYPA","TFRC")#Late erythroid

Dotplot_anno(RBMOs_erythrocytes, features = markers, celltype_color = c("#E41A1C","#4DAF4A","#377EB8"),
             group = c( rep("Early\nerythroid",2),rep("Mid\nerythroid",3),
                        rep("Late\nerythroid",2)),
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             order = T)

RBMOs_erythrocytes <- RenameIdents(RBMOs_erythrocytes, `0` = "Early erythroid", 
                                   `1` = "Mid erythroid", `2` = "Late erythroid")
DimPlot(RBMOs_erythrocytes, label = TRUE)

