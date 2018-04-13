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