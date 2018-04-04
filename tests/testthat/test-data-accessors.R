context("Data Accessors")

# if (!exists("a4")) a4 <- Archs4Repository()

test_that("sample_info returns desired covariate columns", {
  # the code here looks a bit convoluted because it should be cleaned/updated
  # to support testing "universal" covariates, mouse- and human-only covariates
  # as well.

  # these are all human for now
  ids <- tibble(
    id = c('GSE89189', 'GSE29943', "GSM1095128", "GSM1095129", "GSM1095130"),
    organism = "human")

  def.cols <- c("Sample_title", "Sample_source_name_ch1")
  extra.cols <- c("Sample_molecule_ch1", "Sample_treatment_protocol_ch1",
                  "Sample_description")

  # The human data have these covariates that are not in mouse:
  # * Sample_contact_laboratory
  # * Sample_description
  # * Sample_supplementary_file_2
  h.only <- c("Sample_contact_laboratory", "Sample_description",
              "Sample_supplementary_file_2")
  #
  # The mouse data have these covariates thata are not in human:
  # * Sample_contact_state
  # * Sample_growth_protocol_ch1
  # * Sample_treatment_protocol_ch1
  m.only <- c("Sample_contact_state", "Sample_growth_protocol_ch1",
              "Sample_treatment_protocol_ch1")

  all.cols <- c(def.cols, extra.cols)
  info <- sample_info(a4, ids$id, all.cols)

  for (col in all.cols) {
    expect_is(info[[col]], "character", info = col)
    if (col %in% m.only) {
      expect_true(all(is.na(info[[col]])), info = col)
    } else {
      expect_true(!any(is.na(info[[col]])), info = col)
    }
  }
})
