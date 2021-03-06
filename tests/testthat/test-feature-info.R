context("Feature Info")

if (!exists("a4")) {
  # This is loaded by the testthat/helper-all.R script when testthat is running
  # the unit tests, but included here for convenience when doing interactive
  # test development
  a4 <- Archs4Repository()
}

test_that("augmented feature info largely concordant", {
  a4 <- Archs4Repository("/Users/lianoglou/workspace/data/archs4/v4")
  a5 <- Archs4Repository("/Users/lianoglou/workspace/data/archs4/v5")

  v4.fi <- archs4_feature_info("gene", "mouse")
  v5.fi <- archs4_feature_info("gene", "mouse")
})

test_that("feature-level metadata retrieval works", {
  # retrieve gene-level feature information with unique symbols
  mg <- feature_info(a4, feature_type = "gene", source = "mouse",
                     distinct_symbol = TRUE)

  # all a4name entries should be non NA and nchar() >= 1
  expect_true(all(nchar(mg$a4name) >= 1))
  expect_true(all(!is.na(mg$a4name)))
  expect_is(mg$h5idx, "integer")
  expect_true(all(!is.na(mg$h5idx)))
  # no duplicated symbols
  # isna.symbol <- is.na(mg$symbol) |
  #   tolower(mg$symbol) == "na" |
  #   tolower(mg$symbol) == "null"
  # expect_equal(sum(duplicated(mg$symbol) & !isna.symbol), 0)
  # isna.ens <- is.na(mg$ensembl_id) |
  #   tolower(mg$ensembl_id) == "na" |
  #   tolower(mg$ensembl_id) == "null"
  expect_equal(sum(duplicated(mg$ensembl_id) & !is.na(mg$ensembl_id)), 0)

  # There are some entries that we couldn't get identifiers for, but the ones
  # we got should all be prefixed with the mouse prefix.
  mg.ens <- filter(mg, !is.na(ensembl_id))
  expect_true(all(substr(mg.ens$ensembl_id, 1, 7) == "ENSMUSG"))


  hg <- feature_info(a4, feature_type = "gene", source = "human",
                     distinct_symbol = TRUE)
  # all a4name entries should be non NA and nchar() >= 1
  expect_true(all(nchar(hg$a4name) >= 1))
  expect_true(all(!is.na(hg$a4name)))
  expect_is(hg$h5idx, "integer")
  expect_true(all(!is.na(hg$h5idx)))
  # no duplicated symbols
  expect_equal(sum(duplicated(hg$symbol) & !is.na(hg$symbol)), 0)

  hg.ens <- filter(hg, !is.na(ensembl_id))
  expect_true(all(substr(hg.ens$ensembl_id, 1, 4) == "ENSG"))
})
