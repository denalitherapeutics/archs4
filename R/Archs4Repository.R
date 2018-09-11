#' An interface to a locally downloaded ARCHS4 dataset
#'
#' This instantiates an object that acts as a central broker to handle queries
#' against the ARCHS4 dataset. Please refer to the vignette for instructions
#' on how to setup a local directory to act as an Archs4Repository.
#'
#' @export
#'
#' @param datadir The directory that stores the ARCHS4 data.
#' @return an Arhcs4DataSet object
Archs4Repository <- function(datadir = getOption("archs4.datadir")) {
  kosher.dir <- archs4_local_data_dir_validate(echo = FALSE, datadir)
  if (!isTRUE(kosher.dir)) {
    stop(kosher.dir, call. = FALSE)
  }

  asi <- archs4_sample_table(datadir = datadir)
  gstats <- asi %>%
    filter(!is.na(sample_h5idx_gene)) %>%
    group_by(organism) %>%
    summarize(nseries = length(unique(series_id)),
              nsamples = length(unique(sample_id))) %>%
    mutate(feature_type = "gene")

  tstats <- asi %>%
    filter(!is.na(sample_h5idx_transcript)) %>%
    group_by(organism) %>%
    summarize(nseries = length(unique(series_id)),
              nsamples = length(unique(sample_id))) %>%
    mutate(feature_type = "transcript")

  out <- list(
    datadir = datadir,
    file_info = archs4_file_info(datadir),
    meta = archs4_meta(datadir),
    sample_covariates = archs4_sample_covariates(datadir),
    sample_stats = bind_rows(gstats, tstats),
    sample_table = asi)

  # You *know* we're going to make a "remote" or "service" version of an
  # Archs4Repository, in due time ...
  class(out) <- c("LocalArchs4Repository", "Archs4Repository")
  out
}

#' @export
#' @method print Archs4Repository
print.Archs4Repository <- function(x, ...) {
  cat(format(x, ...), "\n")
}

format.Archs4Repository <- function(x, ...) {
  asi <- sample_table(x)

  out <- paste0(
    "===========================================================\n",
    paste(class(x)[1L], "object\n"),
    "-----------------------------------------------------------\n",
    "datadir: ", x$datadir, "\n\n",
    "mouse data:\n",
    "----------\n",
    "       gene: ", nseries(x, "gene", "mouse"), " series; ",
                     nsamples(x, "gene", "mouse"), " samples\n",
    " transcript: ", nseries(x, "transcript", "mouse"), " series; ",
                     nsamples(x, "transcript", "mouse"), " samples\n\n",
    "human data:\n",
    "-----------\n",
    "       gene: ", nseries(x, "gene", "human"), " series; ",
                     nsamples(x, "gene", "human"), " samples\n",
    " transcript: ", nseries(x, "transcript", "human"), " series; ",
                     nsamples(x, "transcript", "human"), " samples\n",
    "===========================================================\n")
}

nsamples <- function(x, feature_type = c("gene", "transcript"),
                     source = archs4_sources()) {
  assert_class(x, "Archs4Repository")
  feature_type <- match.arg(feature_type)
  source <- match.arg(source)
  take <- x$sample_stats$organism == source &
    x$sample_stats$feature_type == feature_type
  x$sample_stats$nsamples[take]
}

nseries <- function(x, feature_type = c("gene", "transcript"),
                    source = archs4_sources()) {
  assert_class(x, "Archs4Repository")
  feature_type <- match.arg(feature_type)
  source <- match.arg(source)
  take <- x$sample_stats$organism == source &
    x$sample_stats$feature_type == feature_type
  x$sample_stats$nseries[take]
}

#' Retrieves the directory that contains the data for the Archs4Repository
#'
#' @export
#' @param x an `Archs4Repository`
datadir <- function(x, ...) {
  assert_class(x, "Archs4Repository")
  x$datadir
}

#' @export
#' @rdname archs4_feature_info
#' @param x an `Archs4Repository`
feature_info <- function(x, feature_type = "gene", source = "human",
                         augmented = TRUE, ...) {
  assert_class(x, "Archs4Repository")
  assert_choice(feature_type, c("gene", "transcript"))
  assert_choice(source, sources(x))
  archs4_feature_info(feature_type, source, augmented, datadir(x), ...)
}

#' Perform a loose/fuzzy lookup for a gene/transcript feature.
#'
#' This funciton facilitates exploratory data analyses by trying to find gene
#' or transcripts by different type of identifiers (symbol, ensembl_id, etc).
#'
#' @export
#' @param x An Archs4Repository
#' @param query a character string of feature names to look for
#' @param feature_type "gene" or "transcript"
#' @param source organism dataset to lookup
#' @return a tibble of features that match against the query. The first column
#'   is the value of the query itself. If no match is found for a query, its
#'   row is all `NA`.
#' @examples
#' a4 <- Archs4Repository()
#' features <- feature_lookup(a4, c("CFAP65", "PECR", "ENSG00000131408"),
#'                            feature_type = "gene", source ="human")
feature_lookup <- function(x, query, feature_type = "gene",
                           source = "human", ...) {
  assert_class(x, "Archs4Repository")
  assert_character(query, min.len = 1L, any.missing = FALSE)
  assert_choice(feature_type, c("gene", "transcript"))
  query <- unique(query)

  fi <- feature_info(x, feature_type, source, ...)
  if (feature_type == "gene") {
    search <- c("ensembl_id", "symbol", "entrez_id", "a4external")
  } else {
    search <- c("ensembl_id", "gene_id", "symbol")
  }
  search <- intersect(search, colnames(fi))

  idxs <- sapply(query, function(qry) .fuzzy_lookup(fi, search, qry))
  isna <- is.na(idxs)
  if (any(isna)) {
    warning("The following ", feature_type, " queries were not found",
            paste(query[isna], paste = ","))
    # idxs <- idxs[!isna]
    # query <- query[!isna]
  }
  bind_cols(
    tibble(query = query),
    fi[idxs,])
}

.fuzzy_lookup <- function(x, columns, query) {
  assert_data_frame(x)
  assert_subset(columns, colnames(x))
  assert_string(query)
  query <- tolower(query)

  idx <- NA_integer_
  for (cname in columns) {
    vals <- tolower(x[[cname]])
    i <- which(vals == query)
    if (length(i) > 1) {
      warning("More than one row for `", query, "` found -- taking first one",
              call. = TRUE)
      idx <- i[1L]
      break
    } else if (length(i) == 1) {
      idx <- i
    }
  }
  if (is.null(idx)) {
    warning("Could not find feature using provided query: `", oquery, "`",
            immediate. = TRUE)
  }
  idx
}


#' @export
#' @rdname archs4_file_info
#' @param x an `Archs4Repository`
file_info <- function(x) {
  assert_class(x, "Archs4Repository")
  x[["file_info"]]
}

#' @export
#' @rdname archs4_file_path
#' @param x an `Archs4Repository`
file_path <- function(x, key) {
  assert_class(x, "Archs4Repository")
  archs4_file_path(key, file_info = file_info(x), datadir = datadir(x))
}

#' Extract the read depth and normalization factors for the samples
#'
#' @export
#' @param x an `Archs4Repository`
#' @param with_a4libsize If `TRUE`, includes an `a4libsize` column, which
#'   was extracted from the `meta/reads_aligned` hdf5 file. Defaults to `FALSE`.
#' @return a tibble with sample_id, a4libsize, libsize, normfactor
libstats <- function(x, with_a4libsize = FALSE) {
  assert_class(x, "Archs4Repository")
  cols <- c("libsize", "normfactor")
  if (with_a4libsize) cols <- c("a4libsize", cols)
  sample_table(x) %>%
    select(sample_id, !!cols)
}

#' @export
#' @rdname archs4_sample_table
#' @param x an `Archs4Repository`
sample_table <- function(x, ...) {
  assert_class(x, "Archs4Repository")
  x$sample_table
}

#' @export
#' @rdname archs4_series_status
#' @param x an `Archs4Repository`
series_status <- function(x, id, ...) {
  assert_class(x, "Archs4Repository")
  archs4_series_status(id, sample_table = sample_table(x), datadir = datadir(x))
}

#' @export
#' @rdname archs4_meta
#' @param x an `Archs4Repository`
meta <- function(x) {
  assert_class(x, "Archs4Repository")
  x$meta
}

#' @export
#' @rdname archs4_sources
#' @param x an `Archs4Repository`
sources <- function(x) {
  assert_class(x, "Archs4Repository")
  meta(x)[["sources"]]
}

#' @export
#' @rdname archs4_sample_covariates
#' @param x an `Archs4Repository`
sample_covariates <- function(x, ...) {
  assert_class(x, "Archs4Repository")
  x$sample_covariates
}

#' @export
#' @rdname archs4_sample_info
#' @param x an `Archs4Repository`
sample_info <- function(x, id,
                        columns = c("Sample_title", "Sample_source_name_ch1"),
                        check_missing_samples = TRUE, ...) {
  assert_class(x, "Archs4Repository")
  archs4_sample_info(id, columns, sample_table(x), sample_covariates(x),
                     check_missing_samples, datadir(x), ...)
}

