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
  expect_error(ClusterTogether(1, deseq2Result, "DESeq2"),
               "normalizedCounts should be a matirx or data frame.")

  # mock one invalid matrix
  failCounts <- deseq2NormalizedCounts
  failCounts[1, 1] <- "a"
  expect_error(ClusterTogether(failCounts, deseq2Result, "DESeq2"),
               "normalizedCounts should be a numeric")

  # mock one counts that has different genes than the DE result
  failCounts <- deseq2NormalizedCounts[-1, ]
  expect_error(ClusterTogether(failCounts, deseq2Result, "DESeq2"),
               "Rownames of normalizedCounts")

  expect_error(ClusterTogether(deseq2NormalizedCounts, deseq2Result, "deseq2"),
               "Tool should be either DESeq2 or edgeR.")

  expect_error(ClusterTogether(deseq2NormalizedCounts, deseq2Result, "edgeR"),
               "DEResult should be a data frame")

  expect_error(ClusterTogether(edgerNormalizedCounts, edgerResult, "DESeq2"),
               "DEResult should be a DESeqResults object.")
})

# [END]
