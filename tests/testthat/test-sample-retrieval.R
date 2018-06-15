context("Sample and Covariate Retrieval")

if (!exists("a4")) {
  # This is loaded by the testthat/helper-all.R script when testthat is running
  # the unit tests, but included here for convenience when doing interactive
  # test development
  a4 <- Archs4Repository()
}


test_that("All expected samples come back when queried by series", {
  series <- "GSE52564"
  expected <- sample_table(a4) %>%
    filter(series_id == series)

  info <- sample_info(a4, series)
  full <- expected %>%
    full_join(info, by = c("series_id", "sample_id"))

  expect_equal(nrow(full), nrow(expected))
})

test_that("(archs4_)sample_status identifies samples missing from GEO series", {
  missing.none <- "GSE52564" # The Ben Barres dataset has all samples in ARCHS4
  # The blurton jones dataset (GSE89189) used to be missing some, but those have
  # been filled out
  missing.some <- "GSE43366" # The Chiu et al. SOD1 datasets is missing some

  ss.none <- series_status(a4, missing.none)
  expect_true(all(ss.none[["in_archs4"]]))

  ss.some <- series_status(a4, missing.some)
  expect_true(!all(ss.some[["in_archs4"]]))
})

test_that("sample_info call handlies missing IDs gracefully", {
  # GSE43366 should have 42 samples, but is missing 7 of them as of April 8, 2018
  # when we are using "v2" of the datasets.
  expected <- tibble(
    series_id = "GSE43366",
    sample_id = paste0("GSM10611", 43:84))
  ex.missing <- c("GSM1061148", "GSM1061149", "GSM1061150", "GSM1061151",
                  "GSM1061152", "GSM1061153", "GSM1061154")
  ex.present <- setdiff(expected$sample_id, ex.missing)

  # Tests a mix of existing and missing sample identifiers.
  # Missing samples should have NAs in a number of colums. One colume that is
  # always returned is `organism`.
  res <- expect_warning(sample_info(a4, expected$sample_id), "not found")

  # Expect that we have one row for each entry in `expected.all`
  # we don't join on series_id because the sample_id's that were queried for
  # and are missing come back with NA series_id values.
  xx <- expected %>%
    full_join(res, by = c("sample_id"))
  expect_equal(nrow(xx), nrow(expected))
  expect_setequal(res$sample_id, expected$sample_id)

  found <- filter(res, !is.na(organism))
  expect_setequal(found$sample_id, ex.present)

  missed <- filter(res, is.na(organism))
  expect_setequal(missed$sample_id, ex.missing)

  # When we query for all missing samples, we still return a tibble of the same
  # form as res.all, with all NAs where you expect them to be. We don't expect
  # to throw an error.
  amiss <- expect_warning({
    sample_info(a4, ex.missing)
  }, "not found")
  expect_setequal(ex.missing, amiss$sample_id)
  expect_true(all(is.na(amiss$organism)))
})

# Tests that are no longer necessary in v4+ matrices (they were made for v2).
# These tests:
# 1. Looked for missing samples in certain seriesexercised
# 2. Identified discordant sample covariates among the mouse and human datasets

# test_that("(archs4_)sample_info warns when querying series with missing samples", {
#   missing.none <- "GSE52564" # The Ben Barres dataset has all samples in ARCHS4
#   missing.some <- "GSE89189" # The blurton jones iPSC paper is missing some
#
#   # This series should have no missing samples
#   si.none <- expect_silent(sample_info(a4, missing.none))
#
#   # This series has a few missing samples
#   wregex <- sprintf("%s series .*missing samples", missing.some)
#   status.some <- expect_warning(sample_info(a4, missing.some), wregex)
#
#   # as.DGEList should also warn when we are missng samples
#   y <- expect_warning(as.DGEList(a4, missing.some), wregex)
# })

# The v4 data matrices have the same covariates in the data matrices among
# the mouse and human data.
# ------------------------------------------------------------------------------
# test_that("sample_info returns desired covariate columns", {
#   # the code here looks a bit convoluted because it should be cleaned/updated
#   # to support testing "universal" covariates, mouse- and human-only covariates
#   # as well.
#
#   ids.all <- tribble(
#     ~id,         ~type,            ~organism,        ~complete,
#     "GSE69354",  "series",         "mouse",          TRUE,
#     "GSE79525",  "series",         "mouse",          TRUE,
#     "GSE98041",  "series",         "mouse",          TRUE,
#     "GSE85702",  "series",         "mouse",          FALSE,
#     "GSE99095",  "series",         "human",          FALSE,
#     "GSE88681",  "series",         "human",          TRUE)
#
#   ids.query <- tribble(
#     ~id,          ~type,     ~organism,
#     "GSE69354",   "series",  "mouse",
#     "GSE88681",   "series",  "human",
#     "GSM1095128", "sample",  "mouse",
#     "GSM1095129", "sample",  "mouse",
#     "GSM1095130", "sample",  "mouse")
#
#   def.cols <- c("Sample_title", "Sample_source_name_ch1")
#   extra.cols <- c("Sample_molecule_ch1", "Sample_treatment_protocol_ch1",
#                   "Sample_description")
#
#   # The human data have these covariates that are not in mouse:
#   # * Sample_contact_laboratory
#   # * Sample_description
#   # * Sample_supplementary_file_2
#   h.only <- c("Sample_contact_laboratory", "Sample_description",
#               "Sample_supplementary_file_2")
#   #
#   # The mouse data have these covariates thata are not in human:
#   # * Sample_contact_state
#   # * Sample_growth_protocol_ch1
#   # * Sample_treatment_protocol_ch1
#   m.only <- c("Sample_contact_state", "Sample_growth_protocol_ch1",
#               "Sample_treatment_protocol_ch1")
#
#   all.cols <- c(def.cols, extra.cols)
#   info <- sample_info(a4, ids.query$id, all.cols)
#
#   for (col in all.cols) {
#     expect_is(info[[col]], "character", info = col)
#   }
#
#   info.m <- filter(info, organism == "mouse")
#   for (col in intersect(h.only, colnames(info))) {
#     expect_true(all(is.na(info.m[[col]])), info = col)
#   }
#
#   info.h <- filter(info, organism == "human")
#   for (col in intersect(m.only, colnames(info))) {
#     expect_true(all(is.na(info.h[[col]])), info = col)
#   }
# })

