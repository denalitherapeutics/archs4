context("Data Accessors")

if (!exists("a4")) {
  # This is loaded by the testthat/helper-all.R script when testthat is running
  # the unit tests, but included here for convenience when doing interactive
  # test development
  a4 <- Archs4Repository()
}

test_that("feature-level metadata retrieval works", {
  # retrieve gene-level feature information with unique symbols
  mg <- feature_info(a4, feature_type = "gene", source = "mouse",
                     distinct_symbol = TRUE)
  expect_true(all(substr(mg$ensembl_gene_id, 1, 7) == "ENSMUSG"))
  expect_is(mg$h5idx, "integer")
  expect_true(all(!is.na(mg$h5idx)))
  # no duplicated symbols
  expect_equal(sum(duplicated(mg$symbol) & !is.na(mg$symbol)), 0)

  hg <- feature_info(a4, feature_type = "gene", source = "human",
                     distinct_symbol = TRUE)
  expect_true(all(substr(hg$ensembl_gene_id, 1, 4) == "ENSG"))
  expect_is(hg$h5idx, "integer")
  expect_true(all(!is.na(hg$h5idx)))
  # no duplicated symbols
  expect_equal(sum(duplicated(hg$symbol) & !is.na(hg$symbol)), 0)
})