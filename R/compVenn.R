#' Make a Venn Diagram of the differentially expressed genes
#'
#' A function to make a Venn diagram of the differentially expressed genes got from
#' DESeq2 or edgeR. The results are meaningful to compare different tools only when
#' the parameters used for both tools when performing the DGE analysis are similar.
#'
#' @param deseq2Result The result from DESeq2, a DESeqResults object.
#' @param edgerResult A data frame taken from the table component of the returned
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
#' # make the Venn diagram
#' # compVenn(deseq2Result = deseq2Res,
#' #          edgerResult = edgerRes)
#'
#' @references
#'
#' Chen, H. (2021). VennDiagram: Generate High-Resolution Venn and Euler Plots.
#' R package version 1.7.0. https://CRAN.R-project.org/package=VennDiagram
#'
#' @export
#' @import VennDiagram

compVenn <- function(deseq2Result,
                     edgerResult,
                     cutoff = 0.05,
                     filename = "DESeq2 vs edgeR.png") {

  # == Performing checks of user input ====
  if (! setequal(rownames(deseq2Result), rownames(edgerResult))) {
    stop("deseq2Result should contain same genes as edgerResult.")
  }

  if (! is.numeric(cutoff) | cutoff <= 0 | cutoff > 1) {
    stop("Cutoff should be a positive number between 0 and 1.")
  }

  if (! is.character(filename)) {
    stop("Filename should be characters.")
  }

  if (! "DESeqResults" %in% class(deseq2Result)) {
    stop("deseq2Result should be a DESeqResults object.")
  }

  if ((! is.data.frame(edgerResult) |
       ! all(attributes(edgerResult)$names == c("logFC",
                                              "logCPM",
                                              "LR",
                                              "PValue",
                                              "FDR")))) {
    stop("edgerResult should be a data frame taken from the table component of the
         returned value of edgeR topTags function.")
  }

  # == Body ====
  # get all differentially expressed genes with a padj <= cutoff
  deseq2DiffGenes <- rownames(deseq2Result[which(
    deseq2Result$padj <= cutoff & deseq2Result$log2FoldChange != 0),])

  edgerDiffGenes <- rownames(edgerResult[which(
    edgerResult$FDR <= cutoff & edgerResult$logFC != 0),])

  # plot the Venn diagram
  VennDiagram::venn.diagram(list(DESeq2 = deseq2DiffGenes,
                                 edgeR = edgerDiffGenes),
                            filename = filename,
                            imagetype = "png",
                            fill = c("red", "blue"),
                            cex = 1.75,
                            cat.cex = 2,
                            width = 4000,
                            height = 4000,
                            cat.pos = 0,
                            main = "The number of differentally expressed genes detected\n by DESeq2 or/and edgeR",
                            main.cex = 2,
                            margin = 0.1)


}

# [END]
