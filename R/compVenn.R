#' Make a Venn Diagram of the differentially expressed genes
#'
#' A function to make a Venn diagram of the differentially expressed genes got from
#' DESeq2 or edgeR. The results are meaningful to compare different tools only when
#' the parameters used for both tools when performing the DGE analysis are similar.
#'
#' @param deseq2_result The result from DESeq2, a DESeqResults object.
#' @param edger_result A data frame taken from the table component of the returned
#'    value of edgeR topTags function (edgeR::topTags(..)$table).
#' @param cutoff padj cutoff. Default: 0.05
#' @param filename Filename for image output. Default: "DESeq2 vs edgeR.png".
#'
#' @return Plots a figure to the file given by the filename argument.
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
#' # make the Venn diagram
#' # compVenn(deseq2_result = deseq2_res,
#' #          edger_result = edger_res)
#'
#' @export
#' @import VennDiagram

compVenn <- function(deseq2_result,
                     edger_result,
                     cutoff = 0.05,
                     filename = "DESeq2 vs edgeR.png") {

  # get all differentially expressed genes with a padj <= cutoff
  deseq2_diff_genes <- rownames(deseq2_result[which(
    deseq2_result$padj <= cutoff & deseq2_result$log2FoldChange != 0),])

  edger_diff_genes <- rownames(edger_result[which(
    edger_result$FDR <= cutoff & edger_result$logFC != 0),])

  # plot the Venn diagram
  VennDiagram::venn.diagram(list(DESeq2 = deseq2_diff_genes,
                                 edgeR = edger_diff_genes),
                            filename = filename,
                            imagetype = "png",
                            fill = c("red", "blue"),
                            cex = 1.75,
                            cat.cex = 2,
                            width = 5000,
                            height = 5000,
                            cat.pos = 0)


}

# [END]
