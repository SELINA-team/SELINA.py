suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(dbscan))
suppressMessages(library(gtools))
suppressMessages(library(presto))

option_list <- list(
    make_option(c('--mode'), type = 'character', help = 'path of celltype'),
    make_option(c('--seurat'), type = 'character', help = 'path of seurat object'),
    make_option(c('--pred_label'), type = 'character', help = 'path of predicted cell type'),
    make_option(c('--pred_prob'), type = 'character', help = 'path of prediction probability'),
    make_option(c('--cell_cutoff'), type = 'integer', help = 'cutoff for cells with the same cell type in 10 nearest neighbor cells'),
    make_option(c('--prob_cutoff'), type = 'double', help = 'cutoff for prediction probability'),
    make_option(c('--path_out'), type = 'character', help = 'path of output'),
    make_option(c('--outprefix'), type = 'character', help = 'prefix of output files')
)

opt_parser <- OptionParser(option_list = option_list);
opts <- parse_args(opt_parser);

mode <- opts$mode
seurat <- opts$seurat
pred_label <- opts$pred_label
pred_prob <- opts$pred_prob
cell_cutoff <- opts$cell_cutoff
prob_cutoff <- opts$prob_cutoff
path_out <- opts$path_out
outprefix <- opts$outprefix


#find differentially expressed genes
FindMarkers <- function(object, celltypes, features = NULL, min.pct = 0.1, logfc.threshold = 0.25,
                                 only.pos = FALSE, return.thresh = 1e-2,
                                 slot = "data") {
  matrix = GetAssayData(object, slot = slot)
  features = rownames(matrix)
  y <- celltypes
  y <- factor(y)
  test.res = wilcoxauc(matrix, y)

  # Calculate logFC
  if (slot != "scale.data") {
    if (slot == "data") {
      X <- expm1(matrix)
    }
    group_sums <- sumGroups(X, y, 1)
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
      log(group_means[, g] + 1) - log(((cs - group_sums[g,]) / (length(y) - gs[g])) + 1)
    }))
  } else {
    group_sums = sumGroups(X, y, 1)
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
      group_means[, g] - ((cs - group_sums[g,]) / (length(y) - gs[g]))
    }))
  }

  test.res$avg_logFC <- as.vector(lfc)
  res <- test.res[, c("pval", "avg_logFC", "pct_in", "pct_out", "padj", "group", "feature")]
  res[, c("pct_in", "pct_out")] = round(res[, c("pct_in", "pct_out")] / 100, digits = 3)
  colnames(res) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "celltype", "gene")
  res <- res %>% dplyr::filter(.data$p_val < return.thresh &
                         abs(.data$avg_logFC) > logfc.threshold &
                         (.data$pct.1 > min.pct |
                         .data$pct.2 > min.pct) &
                         .data$gene %in% features)
  if (only.pos) {
    res <- res %>% dplyr::filter(.data$avg_logFC > 0)
  }
  res <- res %>% dplyr::arrange(.data$celltype, .data$p_val, desc(.data$avg_logFC))
  return(res)
}

#filter cells with low prediction probability
pred_filter <- function(SeuratObj,pred_label,pred_prob,cell_cutoff,prob_cutoff){
  neighbor_id <- kNN(SeuratObj@reductions[["umap"]]@cell.embeddings, k = 10)$id
  pred_prob <- apply(pred_prob,MARGIN = 1,max)
  filtered_label <- pred_label
  for (cell_index in 1:length(pred_label)){
    neighbor_cts <- pred_label[neighbor_id[cell_index,]]
    ct <- pred_label[cell_index] 
    prob <- pred_prob[cell_index]
    nsame <- length(which(neighbor_cts == ct))
    if (nsame<cell_cutoff | prob<prob_cutoff){
      filtered_label[cell_index] <- 'Unknown'
    }
  }
  return(filtered_label)
}

theme_box <- function(...){
  require(grid)
  theme(rect=element_rect(fill='transparent'),
  panel.background=element_rect(fill='transparent',color='black'),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  plot.title = element_text(family="Helvetica",size = 15, vjust = 0.5, hjust = 0.5, margin = margin(t=10,b=10)),
  axis.text.x = element_text(family="Helvetica",size=15, hjust = 0.5 ,vjust = 0.5, margin = margin(t=5,b=5),colour = "black"),
  axis.text.y = element_text(family="Helvetica",size=15, margin = margin(l=5,r=5),colour = "black"),
  axis.title.x = element_text(family="Helvetica",size = 14,colour = "black"),
  axis.title.y = element_text(family="Helvetica",size = 14,colour = "black"),
  axis.ticks.y = element_line(size=0.1),
  legend.title = element_text(family="Helvetica",size = 13,colour = "black"),
  legend.text  = element_text(family="Helvetica",size = 13,colour = "black"),
  legend.key.size = unit(20,'pt'),
  legend.position="none",
  legend.direction = "horizontal",
  legend.key = element_blank(),
  )
}

#generate plots and files indicating the prediction quality for each cluster
cluster_quality <- function(SeuratObj,filtered_label,pred_prob,path_out){
  pred_prob <- apply(pred_prob,MARGIN = 1,max)
  unknown_percent <- c()
  for (cluster in 1:length(unique(SeuratObj$seurat_clusters))){
    cluster_index <- which(SeuratObj$seurat_clusters == as.character(cluster-1))
    unknown_percent <- c(unknown_percent,length(which(filtered_label[cluster_index] == 'Unknown'))/length(cluster_index))
    names(unknown_percent)[cluster] = as.character(cluster-1)
  }
  pred_prob <- data.frame(prob = pred_prob, cluster = as.character(SeuratObj$seurat_clusters))
  unknown_percent <- data.frame(unknown_percent = unknown_percent, cluster = names(unknown_percent))
  pred_prob$cluster <- factor(pred_prob$cluster, levels <- mixedsort(unique(pred_prob$cluster)))
  unknown_percent$cluster <- factor(unknown_percent$cluster, levels <- mixedsort(unique(pred_prob$cluster)))
  png(file.path(path_out, paste0(outprefix, "_cluster_prob.png")), width = 300+120*length(unique(SeuratObj$seurat_clusters)), height = 1500, res = 300)
  print(
    ggplot(pred_prob) + 
    geom_boxplot(aes(x=cluster,y=prob), colour = '#4DBBD5CC', width = 0.6,outlier.shape = NA, lwd=0.2) +
    labs(title = '', y = 'Prob.', x = 'Cluster')+
    theme_box()
  )
  dev.off()
  png(file.path(path_out, paste0(outprefix, "_unknown_percent.png")), width = 300+120*length(unique(SeuratObj$seurat_clusters)), height = 1500, res = 300)
  print(
    ggplot(unknown_percent) + 
    geom_col(aes(x=cluster,y=unknown_percent), colour = '#91D1C2CC', fill = '#91D1C2CC', position=position_dodge(0.7), width=0.7) +
    labs(title = '', y = 'Unknown percent.', x = 'Cluster')+
    theme_classic() +
    theme(
        plot.title = element_text(family="Helvetica",size = 13, vjust = 0.5, hjust = 0.5, margin = margin(t=10,b=10)),
        axis.text.x = element_text(family="Helvetica",size=10, hjust = 0,vjust = 0.5, angle = -90, margin = margin(t=5,b=5)),
        axis.text.y = element_text(family="Helvetica",size=10, margin = margin(l=5,r=5)),
        axis.title.x = element_text(family="Helvetica",size = 12),
        axis.title.y = element_text(family="Helvetica",size = 12),
        legend.title = element_text(family="Helvetica",size = 11),
        legend.text  = element_text(family="Helvetica",size = 10),
        legend.key.size = unit(12,'pt')
        )
  )
  dev.off()
  write.table(pred_prob,file.path(path_out, paste0(outprefix, "_prob.txt")),col.names = TRUE,row.names = FALSE, quote = FALSE, sep = '\t')
  write.table(unknown_percent,file.path(path_out, paste0(outprefix, "_unknown_percent.txt")),col.names = TRUE,row.names = FALSE, quote = FALSE, sep = '\t')
  return()
}

#load data
SeuratObj <- readRDS(seurat)
pred_label <- read.table(pred_label, header = FALSE, sep = '\t')[, 1]
pred_prob <- read.table(pred_prob, header = TRUE, sep = '\t', row.names = 1)

#main step
if (mode == 'single') {
  #filter cells with low prediction score
  filtered_label <- pred_filter(SeuratObj,pred_label,pred_prob,cell_cutoff,prob_cutoff)
  #generate plots and files indicating the prediction quality for each cluster
  cluster_quality(SeuratObj,filtered_label,pred_prob,path_out)
  #find differentially expressed genes for each cell type
  if (length(unique(filtered_label[filtered_label!='Unknown']))>1){
    cluster.genes <- FindMarkers(object = SeuratObj[,filtered_label!='Unknown'], celltypes = filtered_label[filtered_label!='Unknown'])
    write.table(cluster.genes, file.path(path_out, paste0(outprefix, "_DiffGenes.tsv")), quote = F, sep = "\t", row.names = FALSE)
  }
  #output umap plot with predicted cell type labels
  SeuratObj$pred <- filtered_label
  p <- DimPlot(object = SeuratObj[,filtered_label!='Unknown'], label = TRUE, pt.size = 0.2, repel = TRUE, group.by = 'pred')
  ggsave(file.path(path_out, paste0(outprefix, "_pred.png")), p, width = 7, height = 5)
  write.table(data.frame(Cell=colnames(SeuratObj),Prediction=filtered_label,Cluster=SeuratObj$seurat_clusters),file.path(path_out, paste0(outprefix, "_predictions.txt")),col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

} else {
  SeuratObj$pred <- 'pred'
  prob <- apply(pred_prob,MARGIN = 1,max)
  pred_label[prob<prob_cutoff] <- 'Unknown'
  for (i in 1:length(pred_label)) {
    SeuratObj$pred[SeuratObj$seurat_clusters == as.character(i - 1)] = pred_label[i]
  }
  if (length(unique(SeuratObj$pred[SeuratObj$pred!='Unknown']))>1){
    cluster.genes <- FindMarkers(object = SeuratObj[,SeuratObj$pred!='Unknown'], celltypes = SeuratObj$pred[SeuratObj$pred!='Unknown'])
    write.table(cluster.genes, file.path(path_out, paste0(outprefix, "_DiffGenes.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)
  }
  p <- DimPlot(object = SeuratObj[,SeuratObj$pred!='Unknown'], label = TRUE, pt.size = 0.2, repel = TRUE, group.by = 'pred')
  ggsave(file.path(path_out, paste0(outprefix, "_pred.png")), p, width = 7, height = 5)
  write.table(data.frame(Cell=colnames(SeuratObj),Prediction=SeuratObj$pred,Cluster=SeuratObj$seurat_clusters),file.path(path_out, paste0(outprefix, "_predictions.txt")),col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
}