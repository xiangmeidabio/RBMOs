library(CellChat)
library(Seurat)

# load data---------------------------------------------------------------------
load("/media/xiangmeida/meida/8-Organoids/RData/A_AC.RData")
A_AC@meta.data$group <- "BMO-2023"
DefaultAssay(A_AC) <- 'RNA'

load("/media/xiangmeida/meida/8-Organoids/RData/RBMOs.RData")
RBMOs@meta.data$group <- "RBMOs"
DefaultAssay(RBMOs) <- 'RNA'

# cell communication analysis (A_AC)--------------------------------------------
A_AC_input <- GetAssayData(A_AC, layer = 'data')
levels(A_AC)
A_AC_meta <- A_AC@meta.data[,c("group","celltype")]
colnames(A_AC_meta) <-  c("group","labels")

identical(colnames(A_AC_input),rownames(A_AC_meta)) 

A_AC.cellchat <- createCellChat(object = A_AC_input, meta = A_AC_meta, group.by = "labels")
levels(A_AC.cellchat@idents)
groupSize <- as.numeric(table(A_AC.cellchat@idents))
CellChatDB <- CellChatDB.human
A_AC.cellchat@DB <- CellChatDB

#Preprocess expression data for cell-cell communication analysis
A_AC.cellchat <- subsetData(A_AC.cellchat) 
future::plan("multisession", workers = 20)
A_AC.cellchat <- identifyOverExpressedGenes(A_AC.cellchat)
A_AC.cellchat <- identifyOverExpressedInteractions(A_AC.cellchat)

#Inferring cellular communication networks
A_AC.cellchat <- computeCommunProb(A_AC.cellchat, type = "triMean")
A_AC.cellchat <- filterCommunication(A_AC.cellchat, min.cells = 10)
A_AC.cellchat <- computeCommunProbPathway(A_AC.cellchat)
A_AC.cellchat <- aggregateNet(A_AC.cellchat)

#Data extraction
df.net <- subsetCommunication(A_AC.cellchat)

#Visualization
groupSize <- as.numeric(table(A_AC.cellchat@idents)) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(A_AC.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(A_AC.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

save(A_AC.cellchat,file = "/media/xiangmeida/meida/8-Organoids/RData/A_AC.cellchat.RData")

#Heatmap displays interactive data
pheatmap::pheatmap(A_AC.cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")

A_AC.cellchat@netP$pathways
pathways.show <- c("MK")

#Signal pathway
A_AC.cellchat <- netAnalysis_computeCentrality(A_AC.cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(A_AC.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

selected_pathways <- c("CD99","LAMININ",
                       "FN1","BMP","SPP1",
                       "EPHA","CADM","MPZ",
                       "JAM","PECAM1","NCAM",
                       "CDH5","THY1","GAP",
                       "ADGRG","ICAM","ADGRB",
                       "DESMOSOME","TENASCIN","CDH",
                       "CNTN","CLDN3","OCLN")

ht1 <- netAnalysis_signalingRole_heatmap(A_AC.cellchat, pattern = "outgoing", signaling = selected_pathways, width=5, height=12)
ht2 <- netAnalysis_signalingRole_heatmap(A_AC.cellchat, pattern = "incoming", signaling = selected_pathways, width=5, height=12)
ht1 + ht2

selected_pathways <- c("NOTCH","PDGF","TGFb","KIT",
                       "WNT","CSF","ACTIVIN")

ht1 <- netAnalysis_signalingRole_heatmap(A_AC.cellchat, pattern = "outgoing", signaling = selected_pathways, width=4, height=2)
ht2 <- netAnalysis_signalingRole_heatmap(A_AC.cellchat, pattern = "incoming", signaling = selected_pathways, width=4, height=2)
ht1 + ht2

# cell communication analysis (RBMOs)--------------------------------------------
RBMOs_input <- GetAssayData(RBMOs, layer = 'data')
levels(RBMOs)
RBMOs_meta <- RBMOs@meta.data[,c("group","celltype")]
colnames(RBMOs_meta) <-  c("group","labels")

identical(colnames(RBMOs_input),rownames(RBMOs_meta)) 

RBMOs.cellchat <- createCellChat(object = RBMOs_input, meta = RBMOs_meta, group.by = "labels")
levels(RBMOs.cellchat@idents)
groupSize <- as.numeric(table(RBMOs.cellchat@idents))
CellChatDB <- CellChatDB.human
RBMOs.cellchat@DB <- CellChatDB

#Preprocess expression data for cell-cell communication analysis
RBMOs.cellchat <- subsetData(RBMOs.cellchat) 
future::plan("multisession", workers = 20)
RBMOs.cellchat <- identifyOverExpressedGenes(RBMOs.cellchat)
RBMOs.cellchat <- identifyOverExpressedInteractions(RBMOs.cellchat)

#Inferring cellular communication networks
RBMOs.cellchat <- computeCommunProb(RBMOs.cellchat, type = "triMean")
RBMOs.cellchat <- filterCommunication(RBMOs.cellchat, min.cells = 10)
RBMOs.cellchat <- computeCommunProbPathway(RBMOs.cellchat)
RBMOs.cellchat <- aggregateNet(RBMOs.cellchat)

#Data extraction
df.net <- subsetCommunication(RBMOs.cellchat)

#Visualization
groupSize <- as.numeric(table(RBMOs.cellchat@idents)) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(RBMOs.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(RBMOs.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

save(RBMOs.cellchat,file = "/media/xiangmeida/meida/8-Organoids/RData/RBMOs.cellchat.RData")

#Heatmap displays interactive data
pheatmap::pheatmap(RBMOs.cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")

RBMOs.cellchat@netP$pathways
pathways.show <- c("MK")

#Signal pathway
RBMOs.cellchat <- netAnalysis_computeCentrality(RBMOs.cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(RBMOs.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

selected_pathways <- c("CD99","LAMININ",
                       "FN1","BMP","SPP1",
                       "EPHA","CADM","MPZ",
                       "JAM","PECAM1","NCAM",
                       "CDH5","THY1","GAP",
                       "ADGRG","ICAM","ADGRB",
                       "DESMOSOME","TENASCIN","CDH",
                       "CNTN","CLDN3","OCLN")

ht1 <- netAnalysis_signalingRole_heatmap(RBMOs.cellchat, pattern = "outgoing", signaling = selected_pathways, width=5, height=12)
ht2 <- netAnalysis_signalingRole_heatmap(RBMOs.cellchat, pattern = "incoming", signaling = selected_pathways, width=5, height=12)
ht1 + ht2

selected_pathways <- c("NOTCH","PDGF","TGFb","KIT",
                       "WNT","CSF","ACTIVIN")

ht1 <- netAnalysis_signalingRole_heatmap(RBMOs.cellchat, pattern = "outgoing", signaling = selected_pathways, width=4, height=2)
ht2 <- netAnalysis_signalingRole_heatmap(RBMOs.cellchat, pattern = "incoming", signaling = selected_pathways, width=4, height=2)
ht1 + ht2


