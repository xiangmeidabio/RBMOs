library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(SeuratWrappers)

# monocle3----------------------------------------------------------------------
cds2 <- as.cell_data_set(x = RBMOs_Blood_RBC.seu)

cds2@clusters$UMAP$clusters <- Idents(RBMOs_Blood_RBC.seu)[rownames(colData(cds2))]
cds2@clusters$UMAP$partitions <- factor(x = rep(1, length(rownames(colData(cds2)))), levels = 1)
names(cds2@clusters$UMAP$partitions) <- rownames(colData(cds2))
cds2 <- estimate_size_factors(cds2)
cds2@rowRanges@elementMetadata@listData$gene_short_name <- rownames(cds2)
rownames(cds2@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
colnames(cds2@int_colData@listData$reducedDims@listData$UMAP) <- NULL

cds2 <- learn_graph(cds2, use_partition = F)
cds2 <- order_cells(cds2)
ps_tim <- pseudotime(cds2)

plot_cells(cds2, color_cells_by = "pseudotime", label_cell_groups=F, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)+
  labs(x = "UMAP1", y = "UMAP2")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

color_celltype<- c("#33A02C","#FDBF6F","#A6CEE3")
plot_cells(cds2, label_groups_by_cluster = T, label_leaves = F, cell_size = 0.5,label_branch_points = F,
           show_trajectory_graph = TRUE,color_cells_by = "celltype",label_cell_groups = F)+
  scale_color_manual(values = color_celltype)+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  theme(plot.title = element_text(size = 16)
  )+labs(x = "UMAP1", y = "UMAP2")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(legend.text=element_text(size=12),
        legend.title=element_blank())+
  scale_color_manual(aesthetics ="white")

# Pseudo temporal differential genes--------------------------------------------
mycds2_res <- graph_test(cds2, neighbor_graph="principal_graph", cores=4)
genes <- row.names(subset(mycds2_res, q_value< 0.01 & morans_I > 0.2))
mycds2 <- cds2

plot_matrix <- exprs(mycds2)[match(genes,rownames(rowData(mycds2))),order(pseudotime(mycds2))]
dim(plot_matrix)

#Merge functions
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- genes;
dim(plot_matrix)

cutColumn_Means <- function(data_exp,cut
){
  plot_matrix_combin <- list()
  nums <- ncol(data_exp)/cut
  if (nums-round(nums, 0)==0){
    
    for (i in 1:length(seq(1, ncol(data_exp), cut))){
      num <- seq(1, ncol(data_exp), cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
      
    }
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
    
  }else{
    
    for (i in 1:length(seq(1, ncol(data_exp)-cut, cut))){
      num <- seq(1, ncol(data_exp)-cut, cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
    }
    
    plot_matrix_combin[[length(seq(1, ncol(data_exp)-cut, cut))+1]] <- as.data.frame(rowMeans(data_exp[,(max(num)+cut):ncol(data_exp)]))                       
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
  }
}

plot_test <- cutColumn_Means(plot_matrix,cut = 25)
dim(plot_test)

#Sort Settings
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

#Heatmap
p1 <- pheatmap::pheatmap(plot_test, 
                         useRaster = T,
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         cutree_rows=4,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                         clustering_callback = callback)
p1

#Decorate images
annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 3)))
row.names(annotation_row) <- rownames(plot_test)
rowcolor <- c("#85B22E","#E29827","#922927") 
names(rowcolor) <- c("1","2","3")
levels(annotation_row$Cluster)

ann_colors <- list(Cluster=rowcolor)

annotation_col = data.frame(celltype=factor(c(rep("Early erythroid",21),rep("Mid erythroid",26),rep("Late erythroid",117))))
row.names(annotation_col) = colnames(plot_test)

p2 <- pheatmap::pheatmap(plot_test, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(200),
                         annotation_colors=ann_colors,
                         annotation_row = annotation_row,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         #annotation_col = annotation_col,
                         main="Pseudotime")
p2





