#' Make MA plots of the genes
#'
#' A function to make MA plots of the genes. The x-axis is the gene counts and
#' y-axis is log fold change. The gene counts and log fold change are taken from
#' DESeq2 and edgeR results respectively in the two MA plots shown as output. Plots
#' are colored by the tools which consider the genes as differentially expressed.
#' The results are meaningful to compare different tools only when the parameters
#' used for both tools when performing the DGE analysis are similar.
#'
#' @param deseq2_result The result from DESeq2, a DESeqResults object.
#' @param edger_result A data frame taken from the table component of the returned
#'    value of edgeR topTags function (edgeR::topTags(..)$table).
#' @param cutoff padj cutoff. Default: 0.05
#'
#' @return Plots two MA plots side by side, which shows log fold change vs. counts
#' gotten from DESeq2 and edgeR results respectively.
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
#' # make the MA plots
#' compMA(deseq2_result = deseq2_res,
#'        edger_result = edger_res)
#'
#' @export
#' @import ggplot2
#' @import gridExtra

compMA <- function(deseq2_result,
                   edger_result,
                   cutoff = 0.05) {

  # Performing checks of user input
  if (!setequal(rownames(deseq2_result), rownames(edger_result))) {
    stop("deseq2_result should contain same genes as edger_result")
  }

  if (cutoff <= 0 | cutoff > 1) {
    stop("Cutoff should be a positive number between 0 and 1.")
  }

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

  # plot the MA plot, with the gene counts as x-axis and log 2 fold
  # change as y axis
  baseMean <- deseq2_data$baseMean
  log2FoldChange <- deseq2_data$log2FoldChange
  tool <- deseq2_data$tool

  m1 <- ggplot2::ggplot(data=deseq2_data,
                        ggplot2::aes(x=baseMean, y=log2FoldChange, color=tool)) +
    ggplot2::geom_point(alpha=1, size=0.2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
                        colour = "blue", size = 0.5) +
    ggplot2::ylim(c(min(deseq2_data$log2FoldChange),
                    max(deseq2_data$log2FoldChange))) +
    ggplot2::xlim(c(0.1, 50000)) +
    ggplot2::ggtitle("plot with DESeq2 data")  +
    ggplot2::scale_color_manual(values = c("  None" = "grey",
                                  "  deseq2" = "#E69F00",
                                  "  deseq2 edgeR" = "#009E73",
                                  "  edgeR" = "red3"))

  logCPM <- edger_data$logCPM
  logFC <- edger_data$logFC
  tool <- edger_data$tool

  m2 <- ggplot2::ggplot(data=edger_data,
                        ggplot2::aes(x=logCPM, y=logFC, color=tool)) +
    ggplot2::geom_point(alpha=1, size=0.2) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
                        colour = "blue", size = 0.5) +
    ggplot2::ylim(c(min(edger_data$logFC), max(edger_data$logFC))) +
    ggplot2::xlim(c(min(edger_data$logCPM), max(edger_data$logCPM))) +
    ggplot2::ggtitle("plot with edgeR data")  +
    ggplot2::scale_color_manual(values = c("  None" = "grey",
                                  "  deseq2" = "#E69F00",
                                  "  deseq2 edgeR" = "#009E73",
                                  "  edgeR" = "red3"))

  gridExtra::grid.arrange(m1, m2, nrow=1, ncol=2)
}

# [END]
