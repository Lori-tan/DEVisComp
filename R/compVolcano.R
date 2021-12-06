#' Make Volcano plots of the genes
#'
#' A function to make Volcano plots of the genes. The x-axis is the log2 fold
#' change and y-axis is -log10 padj. The log2 fold change and padj are taken from
#' DESeq2 and edgeR results respectively in the two Volcano plots shown as output.
#' Plots are colored by the tools which consider the genes as differentially
#' expressed. The results are meaningful to compare different tools only when the
#' parameters used for both tools when performing the DGE analysis are similar.
#'
#' @param deseq2Result The result from DESeq2, a DESeqResults object.
#' @param edgerResult A data frame taken from the table component of the returned
#'    value of edgeR topTags function (edgeR::topTags(..)$table).
#' @param cutoff padj cutoff. Default: 0.05
#'
#' @return Plots two Volcano plots side by side, which shows -log10 padj vs. log2
#'    fold change gotten from DESeq2 and edgeR results respectively.
#'
#' @examples
#' # Using airwayCounts and airwayMetadata available with package
#'
#' # perform DGE analysis with DESeq2
#' dds <- DESeq2::DESeqDataSetFromMatrix(countData=airwayCounts,
#'                                       colData=airwayMetadata,
#'                                       design=~dex,
#'                                       tidy=TRUE)
#' dds <- DESeq2::estimateSizeFactors(dds)
#'
#' dds@colData$dex <- relevel(dds@colData$dex, ref = "control")
#' dds <- DESeq2::DESeq(dds)
#'
#' res <- DESeq2::results(dds, tidy=TRUE)
#' deseq2Res <- DESeq2::lfcShrink(dds, coef = "dex_treated_vs_control",
#'                                   type = "apeglm")
#'
#' # perform DGE analysis with edgeR
#' counts <- data.frame(airwayCounts[,-1], row.names = airwayCounts$ensgene)
#'
#' diffList <- edgeR::DGEList(counts, samples = airwayMetadata)
#'
#' dex <- factor(rep(c("control", "treated"), 4))
#' design <- model.matrix(~dex)
#' rownames(design) <- colnames(diffList)
#'
#' diffList <- edgeR::estimateDisp(diffList, design, robust = TRUE)
#' fit <- edgeR::glmFit(diffList, design)
#' lrt <- edgeR::glmLRT(fit)
#' edgerRes <- edgeR::topTags(lrt, n = dim(lrt)[1])$table
#'
#' # make the Volcano plots
#' compVolcano(deseq2Result = deseq2Res,
#'        edgerResult = edgerRes)
#'
#' @export
#' @import ggplot2
#' @import gridExtra

compVolcano <- function(deseq2Result,
                        edgerResult,
                        cutoff = 0.05) {

  # Performing checks of user input
  if (!setequal(rownames(deseq2Result), rownames(edgerResult))) {
    stop("deseq2Result should contain same genes as edgerResult.")
  }

  if (cutoff <= 0 | cutoff > 1) {
    stop("Cutoff should be a positive number between 0 and 1.")
  }

  if (!"DESeqResults" %in% class(deseq2Result)) {
    stop("deseq2Result should be a DESeqResults object.")
  }

  if ((!is.data.frame(edgerResult) |
       !all(attributes(edgerResult)$names == c("logFC",
                                                "logCPM",
                                                "LR",
                                                "PValue",
                                                "FDR")))) {
    stop("edgerResult should be a data frame taken from the table component of the
         returned value of edgeR topTags function.")
  }

  # get all differentially expressed genes with a padj <= cutoff
  deseq2DiffGenes <- rownames(deseq2Result[which(
    deseq2Result$padj <= cutoff & deseq2Result$log2FoldChange != 0),])

  edgerDiffGenes <- rownames(edgerResult[which(
    edgerResult$FDR <= cutoff & edgerResult$logFC != 0),])

  allDiffGenes <- union(deseq2DiffGenes, edgerDiffGenes)

  # Mark genes by the tools which consider them as differentially expressed
  # first use the result data from DESeq2
  deseq2Data <- as.data.frame(cbind(deseq2Result,  tool=" "))

  sel <- match(deseq2DiffGenes, rownames(deseq2Data))
  deseq2Data$tool[sel] <- paste(deseq2Data$tool[sel], "deseq2", sep = " ")

  sel <- match(edgerDiffGenes, rownames(deseq2Data))
  deseq2Data$tool[sel] <- paste(deseq2Data$tool[sel], "edgeR", sep = " ")

  sel <- which(deseq2Data$tool == " ")
  deseq2Data$tool[sel] <- "  None"

  # use the result data from edgeR
  edgerData <- as.data.frame(cbind(edgerResult,  tool=" "))

  sel <- match(deseq2DiffGenes, rownames(edgerData))
  edgerData$tool[sel] <- paste(edgerData$tool[sel], "deseq2", sep = " ")

  sel <- match(edgerDiffGenes, rownames(edgerData))
  edgerData$tool[sel] <- paste(edgerData$tool[sel], "edgeR", sep = " ")

  sel <- which(edgerData$tool == " ")
  edgerData$tool[sel] <- "  None"

  # plot the Volcano plot, with log 2 fold change as x-axis and -log10(padj)
  # as y axis
  deseq2Data$padj <- -log(deseq2Data$padj, 10)

  log2FoldChange <- deseq2Data$log2FoldChange
  padj <- deseq2Data$padj
  tool <- deseq2Data$tool

  v1 <- ggplot2::ggplot(data=deseq2Data,
                  ggplot2::aes(x=log2FoldChange, y=padj, color=tool)) +
    ggplot2::geom_point(alpha=1, size=0.2) +
    ggplot2::xlim(c(min(deseq2Data$log2FoldChange),
                    max(deseq2Data$log2FoldChange))) +
    ggplot2::ylim(c(min(deseq2Data$padj), max(deseq2Data$padj))) +
    ggplot2::labs(title="adjusted p-value vs. log-fold change in DESeq2 data",
                  x="log2 fold change",
                  y="-log10 padj") +
    ggplot2::scale_color_manual(values = c("  None" = "grey",
                                           "  deseq2" = "#E69F00",
                                           "  deseq2 edgeR" = "#009E73",
                                           "  edgeR" = "red3"))


  edgerData$FDR <- -log(edgerData$FDR, 10)

  logFC <- edgerData$logFC
  FDR <- edgerData$FDR
  tool <- edgerData$tool

  v2 <- ggplot2::ggplot(data=edgerData,
                  ggplot2::aes(x=logFC, y=FDR, color=tool)) +
    ggplot2::geom_point(alpha=1, size=0.2) +
    ggplot2::xlim(c(min(edgerData$logFC), max(edgerData$logFC))) +
    ggplot2::ylim(c(min(edgerData$FDR), max(edgerData$FDR))) +
    ggplot2::labs(title="adjusted p-value vs. log-fold change in edgeR data",
                  x="log2 fold change",
                  y="-log10 padj") +
    ggplot2::scale_color_manual(values = c("  None" = "grey",
                                           "  deseq2" = "#E69F00",
                                           "  deseq2 edgeR" = "#009E73",
                                           "  edgeR" = "red3"))

  gridExtra::grid.arrange(v1, v2, nrow=1, ncol=2)
}

# [END]
