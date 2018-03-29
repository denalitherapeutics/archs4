#' Classify a vector of sample or series GEO ID's as such
#'
#' GEO series identifiers all start with GSE and sample identifiers all
#' start with GSM. We use that to identify what types of identifiers are
#' passed into `id`
#'
#' @export
#'
#' @param id a character vector of `GSEnnnnn` or `GSMnnnnn` ids
#' @return a tibble of `unique(id)` indicating if the id is a `"series"`
#'   (GSEnnnnn), `"sample"` (GSMnnnnn), or `"unknown"`.
geo_id_type <- function(id) {
  id <- assert_character(id) %>% unique
  type <- case_when(
    is_geo_series_id(id) ~ "series",
    is_geo_sample_id(id) ~ "sample",
    TRUE                 ~ "unknown")
  tibble(id = id, type = type)
}

