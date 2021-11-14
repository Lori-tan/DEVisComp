#' Draw clustered heatmaps for DGE results.
#'
#' A function to draw the clustered heatmaps for normalized counts, log fold
#' changes, and adjusted p-values obtained from differential gene expression
#' analysis results with DESeq2 or edgeR.
#'
#' @param normalized_counts The normailized counts of the RNA-seq data.
#' @param DE_result The result from DESeq2 or edgeR, which should be a DESeqResults
#'    object or a data frame taken from the table component of the returned value
#'    of edgeR topTags function (edgeR::topTags(..)$table).
#' @param order The clustering order. If not provided, the function will cluster
#'    with correlation and complete method based on the normalized counts to obtain
#'    the order.  Default: NULL.
#' @param tool The tool used for DGE analysis, "DESeq2" or "edgeR".
#'
#' @return Plots three heatmaps side by side, in the order of normalized counts,
#'    log fold change, and padj.
#'
#' @examples
#' # Example 1: DGE analysis result from DESeq2
#' # Using airwayCounts and airwayMetadata available with package
#' # perform DGE analysis with DESeq2
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData=airwayCounts,
#'                                       colData=airwayMetadata,
#'                                       design=~dex,
#'                                       tidy=TRUE)
#' dds <- DESeq2::estimateSizeFactors(dds)
#'
#' normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
#'
#' dds@colData$dex <- relevel(dds@colData$dex, ref = "control")
#' dds <- DESeq2::DESeq(dds)
#'
#' res <- DESeq2::results(dds, tidy=TRUE)
#' shrinkResult <- DESeq2::lfcShrink(dds, coef = "dex_treated_vs_control",
#'                                   type = "apeglm")
#' # plot heatmaps
#' ClusterTogether(normalized_counts = normalized_counts,
#'                             DE_result = shrinkResult,
#'                             tool = "DESeq2")
#'
#' # Example 2: DGE analysis result from edgeR
#' # Using airwayCounts and airwayMetadata available with package
#' counts <- data.frame(airwayCounts[,-1], row.names = airwayCounts$ensgene)
#'
#' # perform DGE analysis with edgeR
#' diff_list <- edgeR::DGEList(counts, samples = airwayMetadata)
#'
#' normalized_counts <- edgeR::cpm(diff_list, normalized.lib.sizes = TRUE)
#'
#' dex <- factor(rep(c("control", "treated"), 4))
#' design <- model.matrix(~dex)
#' rownames(design) <- colnames(diff_list)
#'
#' diff_list <- edgeR::estimateDisp(diff_list, design, robust = TRUE)
#' fit <- edgeR::glmFit(diff_list, design)
#' lrt <- edgeR::glmLRT(fit)
#' res <- edgeR::topTags(lrt, n = dim(lrt)[1])$table
#'
#' # plot heatmaps
#' # ClusterTogether(normalized_counts = normalized_counts,
#' #                             DE_result = res,
#' #                            tool = "edgeR")
#'
#' @export
#' @import stats
#' @import pheatmap
#' @import RColorBrewer
#' @import grDevices
#' @import cowplot

ClusterTogether <- function(normalized_counts,
                            DE_result,
                            order = NULL,
                            tool) {

  # Performing checks of user input
  if (!is.numeric(normalized_counts)) {
    stop("Normalized_counts should be a numeric matirx containing the normalized
         count of genes.")
  }

  if (!is.matrix(normalized_counts) & !is.data.frame(normalized_counts)) {
    stop("Normalized_counts should be a  matirx or data frame.")
  }

  # get the significant genes (adjusted p-value <= 0.05) that are differentially
  # expressed
  if (tool == "DESeq2") {
    diff_df <- as.data.frame(DE_result[which(DE_result$padj <= 0.05
                                             & DE_result$log2FoldChange != 0),])
  } else if (tool == "edgeR") {
    diff_df <- as.data.frame(DE_result[which(DE_result$FDR <= 0.05
                                             & DE_result$logFC != 0), ])
  }

  diff_counts <- normalized_counts[match(rownames(diff_df),
                                         rownames(normalized_counts)),]
  # get the log 3 value of the counts in order to make the plots more readable.
  diff_counts <- log(diff_counts + 0.001, 3)

  # if the user doesn't provide clustering order, we cluster based on the
  # normalized counts of one replicate (control vs treated) and get the order
  if (is.null(order)) {
    # we only want the rows with a non-zero SD for clustering purpose
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

  # get heatmap for the counts
  # get the lower and upper limits used for color
  limit <- stats::quantile(diff_counts, probs = c(0.05, 0.95))
  my_colors <- c(limit[1] - 0.01,
                 seq(limit[1], limit[2], by=0.01),
                 limit[2] + 0.01)
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
  reordered_results <- diff_df[match(rownames(diff_counts[order,]),
                                     rownames(diff_df)),]

  if (tool == "DESeq2"){
    lfc <- reordered_results$log2FoldChange
    padj <- reordered_results$padj
  } else if (tool == "edgeR") {
    lfc <- reordered_results$logFC
    padj <- reordered_results$FDR
  }

  # get heatmap for log 2 fold change
  limit <- stats::quantile(lfc, probs = c(0.05, 0.95))
  my_colors <- c(limit[1] - 0.01,
                 seq(limit[1], limit[2], by=0.01),
                 limit[2] + 0.01)
  my_palette <- c('blue',
                  grDevices::colorRampPalette(
                    rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))
                  (n = length(my_colors)-2), 'red')

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

  # plot three heatmaps side by side
  cowplot::plot_grid(p1[[4]], p2[[4]], p3[[4]], nrow = 1, ncol = 3)
}

# [END]
