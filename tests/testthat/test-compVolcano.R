library(DEVisComp)

test_that("invalid inputs", {

  # get DESeq2 results
  suppressWarnings(
    suppressMessages(
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=airwayCounts,
                                            colData=airwayMetadata,
                                            design=~dex,
                                            tidy=TRUE)))
  dds <- DESeq2::estimateSizeFactors(dds)

  deseq2NormalizedCounts <- DESeq2::counts(dds, normalized=TRUE)

  dds@colData$dex <- relevel(dds@colData$dex, ref = "control")
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, tidy=TRUE)
  deseq2Result <- DESeq2::lfcShrink(dds, coef = "dex_treated_vs_control",
                                     type = "apeglm", quiet = TRUE)

  # get edgeR results
  counts <- data.frame(airwayCounts[,-1], row.names = airwayCounts$ensgene)
  diffList <- edgeR::DGEList(counts, samples = airwayMetadata)
  edgerNormalizedCounts <- edgeR::cpm(diffList, normalized.lib.sizes = TRUE)
  dex <- factor(rep(c("control", "treated"), 4))
  design <- model.matrix(~dex)
  rownames(design) <- colnames(diffList)
  diffList <- edgeR::estimateDisp(diffList, design, robust = TRUE)
  fit <- edgeR::glmFit(diffList, design)
  lrt <- edgeR::glmLRT(fit)
  edgerResult <- edgeR::topTags(lrt, n = dim(lrt)[1])$table

  # Test
  expect_error(compVolcano(deseq2Result, edgerResult[-1, ]),
               "deseq2Result should contain same genes as edgerResult.")

  expect_error(compVolcano(deseq2Result, edgerResult, cutoff = "a"),
               "Cutoff should be a positive number between 0 and 1.")

  expect_error(compVolcano(deseq2Result, edgerResult, cutoff = -1),
               "Cutoff should be a positive number between 0 and 1.")

  expect_error(compVolcano(deseq2Result, edgerResult, cutoff = 2),
               "Cutoff should be a positive number between 0 and 1.")

  expect_error(compVolcano(edgerResult, edgerResult),
               "deseq2Result should be a DESeqResults object.")

  expect_error(compVolcano(deseq2Result, deseq2Result),
               "edgerResult should be a data frame")
})

# [END]
