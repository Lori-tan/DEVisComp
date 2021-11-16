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

  deseq2_normalized_counts <- DESeq2::counts(dds, normalized=TRUE)

  dds@colData$dex <- relevel(dds@colData$dex, ref = "control")
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, tidy=TRUE)
  deseq2_result <- DESeq2::lfcShrink(dds, coef = "dex_treated_vs_control",
                                     type = "apeglm", quiet = TRUE)

  # get edgeR results
  counts <- data.frame(airwayCounts[,-1], row.names = airwayCounts$ensgene)
  diff_list <- edgeR::DGEList(counts, samples = airwayMetadata)
  edger_normalized_counts <- edgeR::cpm(diff_list, normalized.lib.sizes = TRUE)
  dex <- factor(rep(c("control", "treated"), 4))
  design <- model.matrix(~dex)
  rownames(design) <- colnames(diff_list)
  diff_list <- edgeR::estimateDisp(diff_list, design, robust = TRUE)
  fit <- edgeR::glmFit(diff_list, design)
  lrt <- edgeR::glmLRT(fit)
  edger_result <- edgeR::topTags(lrt, n = dim(lrt)[1])$table

  # Test
  expect_error(compVolcano(deseq2_result, edger_result[-1, ]),
               "deseq2_result should contain same genes as edger_result.")

  expect_error(compVolcano(deseq2_result, edger_result, cutoff = "a"),
               "Cutoff should be a positive number between 0 and 1.")

  expect_error(compVolcano(deseq2_result, edger_result, cutoff = -1),
               "Cutoff should be a positive number between 0 and 1.")

  expect_error(compVolcano(deseq2_result, edger_result, cutoff = 2),
               "Cutoff should be a positive number between 0 and 1.")

  expect_error(compVolcano(edger_result, edger_result),
               "deseq2_result should be a DESeqResults object.")

  expect_error(compVolcano(deseq2_result, deseq2_result),
               "edger_result should be a data frame")
})

# [END]
