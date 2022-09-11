FindMarkers <- function(object, cluster, features = NULL, min.pct = 0.1, logfc.threshold = 0.25,
                                 only.pos = FALSE, return.thresh = 1e-2,
                                 slot = "data") {
  matrix = GetAssayData(object, slot = slot)
  features = rownames(matrix)
  y <- cluster
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
  colnames(res) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
  res <- res %>% dplyr::filter(.data$p_val < return.thresh &
                         abs(.data$avg_logFC) > logfc.threshold &
                         (.data$pct.1 > min.pct |
                         .data$pct.2 > min.pct) &
                         .data$gene %in% features)
  if (only.pos) {
    res <- res %>% dplyr::filter(.data$avg_logFC > 0)
  }
  res <- res %>% dplyr::arrange(.data$cluster, .data$p_val, desc(.data$avg_logFC))
  return(res)
}