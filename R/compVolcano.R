#' Make Volcano plots of the genes
#'
#' A function to make Volcano plots of the genes. The x-axis is the log2 fold
#' change and y-axis is -log10 padj. The log2 fold change and padj are taken from
#' DESeq2 and edgeR results respectively in the two Volcano plots shown as output.
#' Plots are colored by the tools which consider the genes as differentially
#' expressed. The results are meaningful to compare different tools only when the
#' parameters used for both tools when performing the DGE analysis are similar.
#'
#' @param deseq2_result The result from DESeq2, a DESeqResults object.
#' @param edger_result A data frame taken from the table component of the returned
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
#' deseq2_res <- DESeq2::lfcShrink(dds, coef = "dex_treated_vs_control",
#'                                   type = "apeglm")
#'
#' # perform DGE analysis with edgeR
#' counts <- data.frame(airwayCounts[,-1], row.names = airwayCounts$ensgene)
#'
#' diff_list <- edgeR::DGEList(counts, samples = airwayMetadata)
#'
#' dex <- factor(rep(c("control", "treated"), 4))
#' design <- model.matrix(~dex)
#' rownames(design) <- colnames(diff_list)
#'
#' diff_list <- edgeR::estimateDisp(diff_list, design, robust = TRUE)
#' fit <- edgeR::glmFit(diff_list, design)
#' lrt <- edgeR::glmLRT(fit)
#' edger_res <- edgeR::topTags(lrt, n = dim(lrt)[1])$table
#'
#' # make the Volcano plots
#' compVolcano(deseq2_result = deseq2_res,
#'        edger_result = edger_res)
#'
#' @export
#' @import ggplot2
#' @import gridExtra

compVolcano <- function(deseq2_result,
                        edger_result,
                        cutoff = 0.05) {

  # get all differentially expressed genes with a padj <= cutoff
  deseq2_diff_genes <- rownames(deseq2_result[which(
    deseq2_result$padj <= cutoff & deseq2_result$log2FoldChange != 0),])

  edger_diff_genes <- rownames(edger_result[which(
    edger_result$FDR <= cutoff & edger_result$logFC != 0),])

  all_diff_genes <- union(deseq2_diff_genes, edger_diff_genes)

  # Mark genes by the tools which consider them as differentially expressed
  # first use the result data from DESeq2
  deseq2_data <- as.data.frame(cbind(deseq2_result,  tool=" "))

  sel <- match(deseq2_diff_genes, rownames(deseq2_data))
  deseq2_data$tool[sel] <- paste(deseq2_data$tool[sel], "deseq2", sep = " ")

  sel <- match(edger_diff_genes, rownames(deseq2_data))
  deseq2_data$tool[sel] <- paste(deseq2_data$tool[sel], "edgeR", sep = " ")

  sel <- which(deseq2_data$tool == " ")
  deseq2_data$tool[sel] <- "  None"

  # use the result data from edgeR
  edger_data <- as.data.frame(cbind(edger_result,  tool=" "))

  sel <- match(deseq2_diff_genes, rownames(edger_data))
  edger_data$tool[sel] <- paste(edger_data$tool[sel], "deseq2", sep = " ")

  sel <- match(edger_diff_genes, rownames(edger_data))
  edger_data$tool[sel] <- paste(edger_data$tool[sel], "edgeR", sep = " ")

  sel <- which(edger_data$tool == " ")
  edger_data$tool[sel] <- "  None"

  # plot the Volcano plot, with log 2 fold change as x-axis and -log10(padj)
  # as y axis
  deseq2_data$padj <- -log(deseq2_data$padj, 10)

  log2FoldChange <- deseq2_data$log2FoldChange
  padj <- deseq2_data$padj
  tool <- deseq2_data$tool

  v1 <- ggplot2::ggplot(data=deseq2_data,
                  ggplot2::aes(x=log2FoldChange, y=padj, color=tool)) +
    ggplot2::geom_point(alpha=1, size=0.2) +
    ggplot2::xlim(c(min(deseq2_data$log2FoldChange),
                    max(deseq2_data$log2FoldChange))) +
    ggplot2::ylim(c(min(deseq2_data$padj), max(deseq2_data$padj))) +
    ggplot2::labs(title="genes over mig-32 in deseq2",
                  x="log2 fold change",
                  y="-log10 padj") +
    ggplot2::scale_color_manual(values = c("  None" = "grey",
                                           "  deseq2" = "#E69F00",
                                           "  deseq2 edgeR" = "#009E73",
                                           "  edgeR" = "red3"))


  edger_data$FDR <- -log(edger_data$FDR, 10)

  logFC <- edger_data$logFC
  FDR <- edger_data$FDR
  tool <- edger_data$tool

  v2 <- ggplot2::ggplot(data=edger_data,
                  ggplot2::aes(x=logFC, y=FDR, color=tool)) +
    ggplot2::geom_point(alpha=1, size=0.2) +
    ggplot2::xlim(c(min(edger_data$logFC), max(edger_data$logFC))) +
    ggplot2::ylim(c(min(edger_data$FDR), max(edger_data$FDR))) +
    ggplot2::labs(title="genes over mig-32 in deseq2",
                  x="log2 fold change",
                  y="-log10 padj") +
    ggplot2::scale_color_manual(values = c("  None" = "grey",
                                           "  deseq2" = "#E69F00",
                                           "  deseq2 edgeR" = "#009E73",
                                           "  edgeR" = "red3"))

  gridExtra::grid.arrange(v1, v2, nrow=1, ncol=2)
}

# [END]
