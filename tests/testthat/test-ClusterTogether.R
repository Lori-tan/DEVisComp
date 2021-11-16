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
  expect_error(ClusterTogether(1, deseq2_result, "DESeq2"),
               "Normalized_counts should be a matirx or data frame.")

  # mock one invalid matrix
  fail_counts <- deseq2_normalized_counts
  fail_counts[1, 1] <- "a"
  expect_error(ClusterTogether(fail_counts, deseq2_result, "DESeq2"),
               "Normalized_counts should be a numeric")

  # mock one counts that has different genes than the DE result
  fail_counts <- deseq2_normalized_counts[-1, ]
  expect_error(ClusterTogether(fail_counts, deseq2_result, "DESeq2"),
               "Rownames of normalized_counts")

  expect_error(ClusterTogether(deseq2_normalized_counts, deseq2_result, "deseq2"),
               "Tool should be either DESeq2 or edgeR.")

  expect_error(ClusterTogether(deseq2_normalized_counts, deseq2_result, "edgeR"),
               "DE_result should be a data frame")

  expect_error(ClusterTogether(edger_normalized_counts, edger_result, "DESeq2"),
               "DE_result should be a DESeqResults object.")
})

# [END]
