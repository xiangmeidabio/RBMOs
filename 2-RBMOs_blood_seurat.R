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
RBMOs_Blood.data <- Read10X(data.dir = "/media/xiangmeida/meida/8-Organoids/data/RBMOs-blood")
RBMOs_Blood <- CreateSeuratObject(counts = RBMOs_Blood.data, project = "RBMOs-Blood", min.cells = 3, min.features = 200)
RBMOs_Blood[["percent.mt"]] <- PercentageFeatureSet(RBMOs_Blood, pattern = "^MT-")

# Data quality control----------------------------------------------------------
VlnPlot(RBMOs_Blood, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RBMOs_Blood <- subset(RBMOs_Blood, subset =  nFeature_RNA < 7500 & nFeature_RNA > 200 & percent.mt < 15)

# data processing---------------------------------------------------------------
RBMOs_Blood <- NormalizeData(RBMOs_Blood)

RBMOs_Blood <- FindVariableFeatures(RBMOs_Blood, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(RBMOs_Blood), 10)
plot1 <- VariableFeaturePlot(RBMOs_Blood)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.RBMOs_Blood <- rownames(RBMOs_Blood)
RBMOs_Blood <- ScaleData(RBMOs_Blood, features = all.genes.RBMOs_Blood)

RBMOs_Blood <- RunPCA(RBMOs_Blood, features = VariableFeatures(object = RBMOs_Blood))

ElbowPlot(RBMOs_Blood)

RBMOs_Blood <- FindNeighbors(RBMOs_Blood, dims = 1:10)
RBMOs_Blood <- FindClusters(RBMOs_Blood, resolution = 2)
head(Idents(RBMOs_Blood), 5)
RBMOs_Blood <- RunUMAP(RBMOs_Blood, dims = 1:10)
DimPlot(RBMOs_Blood, reduction = "umap",label = TRUE)

# Cell annotation---------------------------------------------------------------
RBMOs_Blood <- RenameIdents(RBMOs_Blood, `0`="Erythroid cells",`2`="Erythroid cells",`10`="Erythroid cells",
                            `11` ="Erythroid cells", `20` ="Erythroid cells", `6`="Erythroid cells",
                            `7`="Erythroid cells", `12`="Erythroid cells", `13` ="Erythroid cells",
                            `4` ="Erythroid cells", `3`="Erythroid cells",`5`="Erythroid cells",
                            `8`="Erythroid cells", `14` ="Erythroid cells", `15` ="Erythroid cells",
                            `25`="Mast cells",`17`="Mast cells",`16`="Mast cells",`19`="Mast cells",
                            `1`="Eo/baso/neutrophil",`27`="Eo/baso/neutrophil",
                            `9` =  "Megakaryocytes", `24` ="Megakaryocytes",
                            `22` = "DCs",`18` = "DCs", `21` = "HSCs",
                            `23` = "Macrophages", `26` = "T cells")

DimPlot(RBMOs_Blood, label = TRUE)

DefaultAssay(RBMOs_Blood) <- "RNA"
source('/media/xiangmeida/meida/8-Organoids/Draw/Dotplot_anno.R')
markers <- c("CD34","THY1","ITGA6", #HSCs
             "TFRC","GYPA","SLC4A1", #Erythroid cells
             "CCR3","ITGAM","ADGRE1", #Eo/baso/neutrophil
             "CD33","KIT","ENPP3", #Mast cells
             "GP9","PPBP","ITGA2B", #Megakaryocytes
             "CD68","HLA-DRA", "CD163", "CD14", #Macrophages
             "ITGAX","CD83","LYZ", #DCs
             "CD3D","CD3E","CD3G","NCAM1") #T cells

Dotplot_anno(RBMOs_Blood, features = markers, celltype_color = c("#FF7F00","#4DAF4A","#E41A1C","#FFFF33","#A65628","#377EB8","#984EA3","#F781BF"),
             group = c(rep('HSCs',3),rep('Erythroid cells',3),rep('Eo/baso/neutrophil',3),
                       rep('Mast cells',3),rep('Megakaryocytes',3),rep('Macrophages',4),
                       rep('DCs',3),rep('T cells',4)),
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             order = T)

save(RBMOs_Blood,file = "/media/xiangmeida/meida/8-Organoids/RData/RBMOs_Blood.RData")
saveRDS(RBMOs_Blood, "/media/xiangmeida/meida/8-Organoids/RDS/RBMOs_Blood.rds")

# erythrocytes in RBMOs-blood---------------------------------------------------
RBMOs_Blood <- FindNeighbors(RBMOs_Blood, dims = 1:10)
RBMOs_Blood <- FindClusters(RBMOs_Blood, resolution = 0.5)
head(Idents(RBMOs_Blood), 5)
RBMOs_Blood <- RunUMAP(RBMOs_Blood, dims = 1:10)
DimPlot(RBMOs_Blood, reduction = "umap",label = TRUE)

RBMOs_Blood_RBC = RBMOs_Blood[,RBMOs_Blood@meta.data$seurat_clusters %in% c(0,1,2,3,5,9)]
DimPlot(RBMOs_Blood_RBC, reduction = "umap")

RBMOs_Blood_RBC <- NormalizeData(RBMOs_Blood_RBC)
RBMOs_Blood_RBC <- FindVariableFeatures(RBMOs_Blood_RBC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(RBMOs_Blood_RBC), 10)
plot1 <- VariableFeaturePlot(RBMOs_Blood_RBC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes.RBMOs_Blood_RBC <- rownames(RBMOs_Blood_RBC)
RBMOs_Blood_RBC <- ScaleData(RBMOs_Blood_RBC, features = all.genes.RBMOs_Blood_RBC)
RBMOs_Blood_RBC <- RunPCA(RBMOs_Blood_RBC, features = VariableFeatures(object = RBMOs_Blood_RBC))
ElbowPlot(RBMOs_Blood_RBC)
RBMOs_Blood_RBC <- FindNeighbors(RBMOs_Blood_RBC, dims = 1:15)
RBMOs_Blood_RBC <- FindClusters(RBMOs_Blood_RBC, resolution = 0.7)
head(Idents(RBMOs_Blood_RBC), 5)
RBMOs_Blood_RBC <- RunUMAP(RBMOs_Blood_RBC, dims = 1:15)
DimPlot(RBMOs_Blood_RBC, reduction = "umap",label = TRUE)

#Remove cell cycle
FeaturePlot(RBMOs_Blood_RBC, features = c("PCNA","MKI67"),  max.cutoff = 3, min.cutoff=0, cols = c("grey","red"), reduction = "umap")
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
RBMOs_Blood_RBC.seu <- CellCycleScoring(RBMOs_Blood_RBC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(RBMOs_Blood_RBC.seu@meta.data,20)

cycle <- table(RBMOs_Blood_RBC.seu@meta.data$Phase, RBMOs_Blood_RBC.seu@meta.data$old.ident)
RBMOs_Blood_RBC.seu@meta.data %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase)) + theme_minimal()
DimPlot(RBMOs_Blood_RBC.seu,reduction = "umap")
RBMOs_Blood_RBC.seu <- ScaleData(RBMOs_Blood_RBC.seu, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(RBMOs_Blood_RBC.seu))

RBMOs_Blood_RBC.seu <- RunPCA(RBMOs_Blood_RBC.seu, features = VariableFeatures(object = RBMOs_Blood_RBC.seu))
ElbowPlot(RBMOs_Blood_RBC.seu)
RBMOs_Blood_RBC.seu <- FindNeighbors(RBMOs_Blood_RBC.seu, dims = 1:10)
RBMOs_Blood_RBC.seu <- FindClusters(RBMOs_Blood_RBC.seu, resolution = 0.5)
head(Idents(RBMOs_Blood_RBC.seu), 5)

RBMOs_Blood_RBC.seu <- RunUMAP(RBMOs_Blood_RBC.seu, dims = 1:10)
DimPlot(RBMOs_Blood_RBC.seu, reduction = "umap",label = TRUE)

DimPlot(RBMOs_Blood_RBC.seu,reduction = "umap",group.by = "Phase")

# Cell annotation
markers <- c("PTPRC", #Early erythroid
             "SOX6", #Mid erythroid
             "SLC4A1") #Late erythroid

Dotplot_anno(RBMOs_Blood_RBC.seu, features = markers, celltype_color = dittoColors(),
             group = c( rep("Early erythroid",1),rep("Mid erythroid",1),
                        rep("Late erythroid",1)),
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             order = T)

# UMAP
DefaultAssay(RBMOs_Blood_RBC.seu) <- "RNA"
RBMOs_Blood_RBC.seu <- RenameIdents(RBMOs_Blood_RBC.seu, `0` = "Late erythroid", `1` = "Late erythroid", 
                                    `2` = "Mid erythroid",`3` = "Early erythroid",
                                    `4` = "Late erythroid", `5` = "Mid erythroid",
                                    `6` = "Late erythroid", `7` = "Late erythroid",`8` = "Late erythroid")

RBMOs_Blood_RBC.seu@meta.data$celltype <- ""
RBMOs_Blood_RBC.seu@meta.data$celltype[RBMOs_Blood_RBC.seu@active.ident %in% c("Early erythroid")] <- "Early erythroid"
RBMOs_Blood_RBC.seu@meta.data$celltype[RBMOs_Blood_RBC.seu@active.ident %in% c("Late erythroid")] <- "Late erythroid"
RBMOs_Blood_RBC.seu@meta.data$celltype[RBMOs_Blood_RBC.seu@active.ident %in% c("Mid erythroid")] <- "Mid erythroid"

DimPlot(RBMOs_Blood_RBC.seu, label = TRUE)

saveRDS(RBMOs_Blood_RBC.seu, "/media/xiangmeida/meida/8-Organoids/RDS/RBMOs_Blood_RBC.rds")




