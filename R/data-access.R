archs4.files <- function() {
  files <- c(
    human_gene = file.path('human_matrix.h5'),
    mouse_gene = file.path('mouse_matrix.h5'))
  files
}

#' Identify the file path on the system for specific ARCHS4 files.
#'
#' @export
#'
#' @param what what type of file? `"human_gene"` or `"mouse_gene"`
#' @return the path on the filesystem for the path asked for
archs4_file_path <- function(what, datadir = getOption("archs4.datadir")) {
  afiles <- archs4.files()
  assert_directory(datadir, 'r')
  assert_choice(what, names(afiles))

  fn <- file.path(datadir, afiles[what])
  if (!file.exists(fn)) {
    stop("Can not find archs4 file [", what, "]: ", fn)
  }

  fn
}

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


#' Retrieves information for samples by GSE series or sample IDs
#'
#' Fetch a tibble of series and sample information by querying the arcsh4
#' dataset by GEO sample (GSE) or sample (GSM) ids.
#'
#' @export
#'
#' @param id a character vector of GEO series or sample ids.
#' @param from the dataset to look up. If this is `"human"` or `"mouse"`, then
#'   the default gene expression file for that organism found in
#'   `getOption("archs4.datadir")` is used. Otherwise this can be a path to an
#'   hdf5 ARCHS4 expression file.
#' @param a tibble of sample information.
#' @return a tibble of series, sample_id, sample_title, and sample_name columns.
#'   If the query sample or series query can't be found, then there will be an
#'   `NA` value for these columns. The `query_type` column will indicat whether
#'   the row was returned from querying by series or by sample.
archs4_sample_info <- function(id, from = "human", ...) {
  input <- geo_id_type(id)
  bad.id <- filter(input, type == "unknown")
  if (nrow(bad.id)) {
    stop("Unknown identifiers: ", paste(bad.id$id, collapse = ", "))
  }

  assert_string(from)
  if (from %in% c("human", "mouse")) {
    h5.fn <- archs4_file_path(paste0(from, "_gene"))
  }
  stopifnot(is_archs4_expression_file(h5.fn))

  series.all <- h5read(h5.fn, "meta/Sample_series_id")
  sample.all <- h5read(h5.fn, "meta/Sample_geo_accession")

  out <- pmap_dfr(input, .fetch_sample_info,
                  series.all = series.all, sample.all = sample.all,
                  h5.fn = h5.fn)
  out
}

# helper function to `archs4_sample_info`
.fetch_sample_info <- function(id, type, series.all, sample.all, h5.fn,
                               with.description = FALSE) {
  assert_string(id)
  assert_choice(type, c("series", "sample"))
  assert_file_exists(h5.fn, "r")

  if (type == "series") {
    take <- which(series.all == id)
    sample_id <- sample.all[take]
    series_id <- id
  } else {
    take <- which(sample.all == id)
    sample_id <- id
    series_id <- series.all[take]
  }

  if (length(take)) {
    stitle <- h5read(h5.fn, "meta/Sample_title", list(take))
    sm <- h5read(h5.fn, "meta/Sample_source_name_ch1", list(take))
    if (with.description) {
      desc <- h5read(h5.fn, "meta/Sample_description", list(take))
    }
  } else {
    stitle <- sm <- desc <- NA_character_
    if (type == "series") sample_id <- sm else series_id <- sm
  }

  out <- tibble(series_id = id, sample_id = sample_id, sample_title = stitle,
                sample_name = sm)
  if (with.description) {
    out[['sample_description']] <- desc
  }
  out[['query_type']] <- type

  out
}

if (FALSE) {
  # TODO: create some tests out of these
  series_id <- 'GSE89189'
  id <- c('GSE89189', 'GSE29943', "GSM1095128", "GSM1095129", "GSM1095130")
  h5.fn <- archs4_file_path("human_gene")
}
