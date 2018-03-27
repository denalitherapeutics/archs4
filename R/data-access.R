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

#' Retrieves the feature (gene/transcript) information for the archs4 data
#'
#' Unfortunately only the gene symbols are stored in the gene-level count data
#' in `"meta/genes"`.
#'
#' @importFrom readr read_csv
#'
#' @param feature_type gene or transcript?
#' @param source human or mouse
#' @return a tibble of information
archs4_feature_info <- function(feature_type = c("gene", "transcript"),
                                source = c("human", "mouse"),
                                distinct_symbol = TRUE) {
  feature_type <- match.arg(feature_type)
  if (feature_type == "transcript") {
    stop("transcript level features not supported yet")
  }

  source <- match.arg(source)

  h5.fn <- archs4_file_path(paste(source, feature_type, sep = "_"))

  if (feature_type == "gene") {
    ainfo <- tibble(symbol = rhdf5::h5read(h5.fn, "meta/genes"),
                    h5idx = seq(symbol))

    meta <- local({
      fn <- sprintf("gencode.%s.basic.annotation.csv.gz",
                    if (source == "human") "v27" else "vM16")
      fn <- system.file("extdata", "feature-annotation", fn,
                        package = "archs4", mustWork = TRUE)
      coltypes <- "ccicc"

      out <- readr::read_csv(fn, col_types = coltypes)

      if (distinct_symbol) {
        # there are duplicate entries by symbol. here I pick one by arranging
        # by symbol and entrez_id, this puts NA entrez_ids last
        out <- arrange(out, symbol, entrez_id)
        out <- distinct(out, symbol, .keep_all = TRUE)
      }

      out
    })

    ainfo <- left_join(ainfo, meta, by = "symbol")
  }

  ainfo
}


#' Retrieves information for samples by GSE series or sample IDs
#'
#' Fetch a tibble of series and sample information by querying the arcsh4
#' dataset by GEO sample (GSE) or sample (GSM) ids.
#'
#' @export
#'
#' @param id a character vector of GEO series or sample ids.
#' @param source the dataset to look up. If this is `"human"` or `"mouse"`, then
#'   the default gene expression file for that organism found in
#'   `getOption("archs4.datadir")` is used. Otherwise this can be a path to an
#'   hdf5 ARCHS4 expression file.
#' @param a tibble of sample information.
#' @return a tibble of series, sample_id, sample_h5idx, sample_title, and
#'   sample_name columns. If the query sample or series query can't be found,
#'   then there will be an `NA` value for these columns. The `query_type` column
#'   will indicat whether the row was returned from querying by series or by
#'   sample.
archs4_sample_info <- function(id, source = "human", ...) {
  input <- geo_id_type(id)
  bad.id <- filter(input, type == "unknown")
  if (nrow(bad.id)) {
    stop("Unknown identifiers: ", paste(bad.id$id, collapse = ", "))
  }

  assert_string(source)
  if (source %in% c("human", "mouse")) {
    h5.fn <- archs4_file_path(paste0(source, "_gene"))
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
    stitle <- sm <- desc <- take <- NA_character_
    if (type == "series") sample_id <- sm else series_id <- sm
  }

  out <- tibble(series_id = id, sample_id = sample_id, sample_h5idx = take,
                sample_title = stitle, sample_name = sm)
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
