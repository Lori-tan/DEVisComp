#' Draw clustered heatmaps for DGE results.
#'
#' A function to draw the clustered heatmaps for normalized counts, log fold
#' changes, and adjusted p-values obtained from differential gene expression
#' analysis results with DESeq2 or edgeR. The function clusters with
#' correlation and complete method based on the normalized counts to obtain the
#' clustering order.
#'
#' @param normalizedCounts The normailized counts of the RNA-seq data.
#' @param DEResult The result from DESeq2 or edgeR, which should be a DESeqResults
#'    object or a data frame taken from the table component of the returned value
#'    of edgeR topTags function (edgeR::topTags(..)$table).
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
#' normalizedCounts <- DESeq2::counts(dds, normalized=TRUE)
#'
#' dds@colData$dex <- relevel(dds@colData$dex, ref = "control")
#' dds <- DESeq2::DESeq(dds)
#'
#' res <- DESeq2::results(dds, tidy=TRUE)
#' shrinkResult <- DESeq2::lfcShrink(dds, coef = "dex_treated_vs_control",
#'                                   type = "apeglm")
#' # plot heatmaps
#' ClusterTogether(normalizedCounts = normalizedCounts,
#'                             DEResult = shrinkResult,
#'                             tool = "DESeq2")
#'
#' # Example 2: DGE analysis result from edgeR
#' # Using airwayCounts and airwayMetadata available with package
#' counts <- data.frame(airwayCounts[,-1], row.names = airwayCounts$ensgene)
#'
#' # perform DGE analysis with edgeR
#' diffList <- edgeR::DGEList(counts, samples = airwayMetadata)
#'
#' normalizedCounts <- edgeR::cpm(diffList, normalized.lib.sizes = TRUE)
#'
#' dex <- factor(rep(c("control", "treated"), 4))
#' design <- model.matrix(~dex)
#' rownames(design) <- colnames(diffList)
#'
#' diffList <- edgeR::estimateDisp(diffList, design, robust = TRUE)
#' fit <- edgeR::glmFit(diffList, design)
#' lrt <- edgeR::glmLRT(fit)
#' res <- edgeR::topTags(lrt, n = dim(lrt)[1])$table
#'
#' # plot heatmaps
#' ClusterTogether(normalizedCounts = normalizedCounts,
#'                            DEResult = res,
#'                            tool = "edgeR")
#'
#' @references
#' Kolde, R. (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12.
#' https://CRAN.R-project.org/package=pheatmap
#'
#' Neuwirth, E. (2014). RColorBrewer: ColorBrewer Palettes. R package
#' version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer
#'
#' R Core Team (2021). R: A language and environment for statistical
#' computing. R Foundation for Statistical Computing, Vienna, Austria. URL
#' https://www.R-project.org/.
#'
#' Wilke, C. O. (2020). cowplot: Streamlined Plot Theme and Plot Annotations
#' for 'ggplot2'. R package version 1.1.1.
#' https://CRAN.R-project.org/package=cowplot
#'
#' @export
#' @import stats
#' @import pheatmap
#' @import RColorBrewer
#' @import grDevices
#' @import cowplot

ClusterTogether <- function(normalizedCounts,
                            DEResult,
                            tool) {

  # == Performing checks of user input ====
  if (! is.numeric(normalizedCounts)) {
    stop("normalizedCounts should be a numeric matirx containing the normalized
         count of genes.")
  }

  if (! is.matrix(normalizedCounts) & ! is.data.frame(normalizedCounts)) {
    stop("normalizedCounts should be a matirx or data frame.")
  }

  if (! setequal(rownames(normalizedCounts), rownames(DEResult))) {
    stop("Rownames of normalizedCounts should matches the rownames of DEResult,
         i.e. they should contain same genes.")
  }

  if (tool != "DESeq2" & tool != "edgeR") {
    stop("Tool should be either DESeq2 or edgeR.")
  }

  if (tool == "DESeq2" & ! "DESeqResults" %in% class(DEResult)) {
    stop("DEResult should be a DESeqResults object.")
  }

  if (tool == "edgeR" & (! is.data.frame(DEResult) |
       !all(attributes(DEResult)$names == c("logFC",
                                              "logCPM",
                                              "LR",
                                              "PValue",
                                              "FDR")
            )
       )
      ) {
    stop("DEResult should be a data frame taken from the table component of the
         returned value of edgeR topTags function.")
  }

  # == Body ====
  # get the significant genes (adjusted p-value <= 0.05) that are
  # differentially expressed
  if (tool == "DESeq2") {
    diffDf <- as.data.frame(DEResult[which(DEResult$padj <= 0.05
                                             & DEResult$log2FoldChange != 0),])
  } else if (tool == "edgeR") {
    diffDf <- as.data.frame(DEResult[which(DEResult$FDR <= 0.05
                                             & DEResult$logFC != 0), ])
  }

  diffCounts <- normalizedCounts[match(rownames(diffDf),
                                         rownames(normalizedCounts)),]
  # get the log 3 value of the counts in order to make the plots more readable.
  diffCounts <- log(diffCounts + 0.001, 3)

  # if the user doesn't provide clustering order, we cluster based on the
  # normalized counts of one replicate (control vs treated) and get the order

  # we only want the rows with a non-zero SD for clustering purpose
  SD <- apply(diffCounts[, c(1, 2)], 1, stats::sd)
  if (length(which(SD == 0)) != 0) {
    warning(sprintf("Removing %d rows when clustering due to zero SD in the first two column of normalized counts.",
                    length(which(SD == 0)) != 0))
    diffCounts <- diffCounts[-which(SD == 0), ]
  }

  p <- pheatmap::pheatmap(diffCounts[, c(1, 2)],
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


  # get heatmap for the counts
  # get the lower and upper limits used for color
  limit <- stats::quantile(diffCounts, probs = c(0.05, 0.95))
  myColors <- c(limit[1] - 0.01,
                 seq(limit[1], limit[2], by=0.01),
                 limit[2] + 0.01)
  myPalette <- c(grDevices::colorRampPalette(colors = c("white", "yellow", "red"))
                  (n = length(myColors) - 2), 'red')

  p1 <- pheatmap::pheatmap(diffCounts[order,],
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
                color = myPalette,
                breaks = myColors
  )

  # re-order the DE result
  reorderedResults <- diffDf[match(rownames(diffCounts[order,]),
                                     rownames(diffDf)),]

  if (tool == "DESeq2"){
    lfc <- reorderedResults$log2FoldChange
    padj <- reorderedResults$padj
  } else if (tool == "edgeR") {
    lfc <- reorderedResults$logFC
    padj <- reorderedResults$FDR
  }

  # get heatmap for log 2 fold change
  limit <- stats::quantile(lfc, probs = c(0.05, 0.95))
  myColors <- c(limit[1] - 0.01,
                 seq(limit[1], limit[2], by=0.01),
                 limit[2] + 0.01)
  myPalette <- c('blue',
                  grDevices::colorRampPalette(
                    rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))
                  (n = length(myColors)-2), 'red')

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
                 color = myPalette,
                 breaks = myColors
  )

  # get heatmap for adjusted p-value
  myColors <- c(seq(0,0.05,by=0.0001), 0.0501)
  myPalette <- c(grDevices::colorRampPalette(colors = c("yellow", "red"))
                  (n = length(myColors)-2), 'red')

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
                 color = myPalette,
                 breaks = myColors
  )

  # plot three heatmaps side by side
  cowplot::plot_grid(p1[[4]], p2[[4]], p3[[4]], nrow = 1, ncol = 3)
}

# [END]
