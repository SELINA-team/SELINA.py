suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(presto))

option_list <- list(
    make_option(c('--seurat'), type = 'character', help = 'path of seurat object'),
    make_option(c('--pred_label'), type = 'character', help = 'path of predicted cell type'),
    make_option(c('--pred_prob'), type = 'character', help = 'path of prediction probability'),
    make_option(c('--prob_cutoff'), type = 'double', help = 'cutoff for prediction probability'),
    make_option(c('--path_out'), type = 'character', help = 'path of output'),
    make_option(c('--outprefix'), type = 'character', help = 'prefix of output files')
)

opt_parser <- OptionParser(option_list = option_list);
opts <- parse_args(opt_parser);

seurat <- opts$seurat
pred_label <- opts$pred_label
pred_prob <- opts$pred_prob
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
pred_filter <- function(pred_label,pred_prob,prob_cutoff){
  pred_prob <- apply(pred_prob,MARGIN = 1,max)
  filtered_label <- pred_label
  for (cell_index in 1:length(pred_label)){
    prob <- pred_prob[cell_index]
    if (prob<prob_cutoff){
      filtered_label[cell_index] <- 'Unknown'
    }
  }
  return(filtered_label)
}

#load data
SeuratObj <- readRDS(seurat)
pred_label <- read.table(pred_label, header = FALSE, sep = '\t')[, 1]
pred_prob <- read.table(pred_prob, header = TRUE, sep = '\t', row.names = 1)

#main step
#filter cells with low prediction score
filtered_label <- pred_filter(pred_label,pred_prob,prob_cutoff)
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
