validate.data.dir <- function(datadir = getOption("archs4.datadir")) {
  assert_directory(datadir, 'r')
  out <- TRUE

  req.files <- archs4.files()
  req.files <- sapply(req.files, function(x) file.path(datadir, x))

  f.exists <- setNames(file.exists(req.files), names(req.files))
  mf <- names(req.files)[!f.exists]
  if (length(mf)) {
    if (any(c("human_gene_info", "mouse_gene_info") %in% mf)) {
      msg <- paste(
        "Missing augmented gene information in archs4 data directory:\n ",
        sprintf("[%s]", datadir),
        "\n\nRun the create_augmented_feature_info(datadir) function")
    } else {
      msg <- paste(
        "The following ARCHS4 files need to be downloaded into the datadir:\n ",
        sprintf("[%s]\n", datadir),
        paste(mf, collapse = ", "))
    }
    out <- msg
  }

  out
}

is_archs4_expression_file <- function(fn) {
  assert_character(fn)
  # last N characters must be '_matrix.h5'
  suffix <- '_matrix.h5$'
  isTRUE(length(grep(suffix, fn)[1L]) == 1L)
}

is_geo_series_id <- function(id) {
  assert_character(id)
  substr(id, 1, 3) == "GSE"
}

is_geo_sample_id <- function(id) {
  assert_character(id)
  substr(id, 1, 3) == "GSM"
}