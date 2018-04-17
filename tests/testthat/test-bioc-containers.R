context("Bioconductor Containers")

if (!exists("a4")) {
  a4 <- Archs4Repository()
}

test_that("as.DGEList creates gene- and sample-level DGELists", {
  # This is a human dataset
  scovs <- c("Sample_title", "Sample_source_name_ch1")
  yg <- as.DGEList(a4, "GSE52564", feature_type = "gene",
                   sample_columns = scovs)
  expect_true(all(substr(rownames(yg), 1, 7) == "ENSMUSG"))
  for (cov in scovs) {
    expect_is(yg$samples[[cov]], "character", info = paste("yg:", cov))
  }


  yt <- as.DGEList(a4, "GSE52564", feature_type = "transcript",
                   sample_columns = scovs)
  expect_true(all(substr(rownames(yt), 1, 7) == "ENSMUST"))
  for (cov in scovs) {
    expect_is(yt$samples[[cov]], "character", info = paste("yt:", cov))
  }
})