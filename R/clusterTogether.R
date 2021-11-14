#' Cluster
#'
#' This is a function
#'
#' @param normalized_counts The count
#' @param DE_res The result from DESeq2 or edgeR
#' @param order the clustering order
#' @param tool tool used for DGE analysis, "DESeq2" or "edgeR"
#'
#' @examples
#' # empty
#'
#' @references
#'empty
#'
#' @export
#' @import stats
#' @import pheatmap
#' @import RColorBrewer
#' @import grDevices
#' @import cowplot

ClusterTogether <- function(normalized_counts, DE_res, order, tool) {
  if (!is.numeric(normalized_counts)) {
    stop("Normalized_counts should be a numeric matirx containing the normalized
         count of genes.")
  }

  if (!is.matrix(normalized_counts) & !is.data.frame(normalized_counts)) {
    stop("Normalized_counts should be a  matirx or data frame.")
  }

  # get the significant genes (adjusted p-value <= 0.05) that are differentially
  # expressed
  if (tool == "DESeq2"){
    diff_df <- as.data.frame(DE_res[which(DE_res$padj <= 0.05 & DE_res$log2FoldChange != 0),])
  } else if (tool == "edgeR") {
    diff_df <- as.data.frame(DE_res[which(DE_res$FDR <= 0.05 & DE_res$logFC != 0), ])
  }

  diff_counts <- log(normalized_counts[match(rownames(diff_df), rownames(normalized_counts)),] + 0.001, 3)

  if (is.null(order)) {
    SD <- apply(diff_counts[, c(1, 2)], 1, stats::sd)
    if (length(which(SD == 0)) != 0) {
      diff_counts <- diff_counts[-which(SD == 0), ]
    }

    p <- pheatmap::pheatmap(diff_counts[, c(1, 2)],
                  cluster_cols=FALSE,
                  cluster_rows=TRUE,
                  clustering_distance_rows = 'correlation',
                  clustering_method = 'complete',
                  cellwidth = 10,
                  fontsize = 6,
                  angle_col = c('315'),
                  border_color = "NA",
                  show_rownames = FALSE,
                  main = 'normalized count',
                  cellheight = 0.3
    )

    order <- p$tree_row$order
  }

  # get heatmap for counts
  # get the lower and upper limits used for color
  (limit <- stats::quantile(diff_counts, probs = c(0.05, 0.95)))
  my_colors <- c(limit[1] - 0.01 ,seq(limit[1], limit[2], by=0.01), limit[2] + 0.01)
  my_palette <- c(grDevices::colorRampPalette(colors = c("white", "yellow", "red"))
                  (n = length(my_colors) - 2), 'red')

  p1 <- pheatmap::pheatmap(diff_counts[order,],
                cluster_cols=FALSE,
                cluster_rows=FALSE,
                clustering_distance_rows = 'correlation',
                clustering_method = 'complete',
                cellwidth = 10,
                fontsize = 6,
                angle_col = c('315'),
                border_color = "NA",
                show_rownames = FALSE,
                main = 'normalized count',
                cellheight = 0.3,
                color = my_palette,
                breaks = my_colors
  )

  # re-order the DE result
  reordered_results <- diff_df[match(rownames(diff_counts[order,]), rownames(diff_df)),]

  if (tool == "DESeq2"){
    lfc <- reordered_results$log2FoldChange
    padj <- reordered_results$padj
  } else if (tool == "edgeR") {
    lfc <- reordered_results$logFC
    padj <- reordered_results$FDR
  }

  # get heatmap for log 2 fold change
  (limit <- stats::quantile(lfc, probs = c(0.05, 0.95)))
  my_colors <- c(limit[1] - 0.01 ,seq(limit[1], limit[2], by=0.01), limit[2] + 0.01)
  my_palette <- c('blue', grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(n = length(my_colors)-2), 'red')

  p2 <- pheatmap::pheatmap(lfc,
                 cluster_cols=FALSE,
                 cluster_rows=FALSE,
                 clustering_distance_rows = 'correlation',
                 clustering_method = 'complete',
                 cellwidth = 10,
                 fontsize = 6,
                 angle_col = c('315'),
                 border_color = "NA",
                 show_rownames = FALSE,
                 main = 'LFC',
                 cellheight = 0.3,
                 color = my_palette,
                 breaks = my_colors
  )

  # get heatmap for adjusted p-value
  my_colors <- c(seq(0,0.05,by=0.0001), 0.0501)
  my_palette <- c(grDevices::colorRampPalette(colors = c("yellow", "red"))
                  (n = length(my_colors)-2), 'red')

  p3 <- pheatmap::pheatmap(padj,
                 cluster_cols=FALSE,
                 cluster_rows=FALSE,
                 clustering_distance_rows = 'correlation',
                 clustering_method = 'complete',
                 cellwidth = 10,
                 fontsize = 6,
                 angle_col = c('315'),
                 border_color = "NA",
                 show_rownames = FALSE,
                 main = 'Padj',
                 cellheight = 0.3,
                 color = my_palette,
                 breaks = my_colors
  )

  cowplot::plot_grid(p1[[4]], p2[[4]], p3[[4]], nrow = 1, ncol = 3)
}
