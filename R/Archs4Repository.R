#' An interface to a locally downloaded ARCHS4 dataset
#'
#' This instantiates an object that acts as a central broker to handle queries
#' against the ARCHS4 dataset. You will need to download all of the files
#' enumerated in `archs4:::archs4.files()` and save them in `datadir`.
#'
#' @export
#'
#' @param datadir The directory that stores the ARCHS4 data.
#' @return an Arhcs4DataSet object
Archs4Repository <- function(datadir = getOption("archs4.datadir")) {
  kosher.dir <- validate.data.dir(datadir)
  if (!isTRUE(kosher.dir)) {
    stop(kosher.dir)
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
    sample_table = asi,
    sample_stats = bind_rows(gstats, tstats),
    sample_covariates = archs4_sample_covariates(datadir),
    datadir = datadir)
  class(out) <- "Archs4Repository"
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
    "Archs4Repository object\n",
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
                     source = c("human", "mouse")) {
  assert_class(x, "Archs4Repository")
  feature_type <- match.arg(feature_type)
  source <- match.arg(source)
  take <- x$sample_stats$organism == source &
    x$sample_stats$feature_type == feature_type
  x$sample_stats$nsamples[take]
}

nseries <- function(x, feature_type = c("gene", "transcript"),
                    source = c("human", "mouse")) {
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
feature_info <- function(x, feature_type = c("gene", "transcript"),
                         source = archs4_sources(),
                         distinct_symbol = TRUE, augmented = TRUE, ...) {
  assert_class(x, "Archs4Repository")
  feature_type <- match.arg(feature_type)
  source <- match.arg(source)
  archs4_feature_info(feature_type, source, distinct_symbol, augmented,
                      datadir(x), ...)
}

#' @export
#' @rdname archs4_sample_table
#' @param x an `Archs4Repository`
sample_table <- function(x, feature_type = c("gene", "transcript"), ...) {
  assert_class(x, "Archs4Repository")
  x$sample_table
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
                        ...) {
  assert_class(x, "Archs4Repository")

  # check sample metadata columns
  # scols <- archs4_sample_metadata_names(x$datadir)
  # columns <- unique(columns)
  # assert_character(columns, any.missing = FALSE, min.len = 1L)
  # assert_subset(columns, scols)

  archs4_sample_info(id, columns, sample_table(x), sample_covariates(x),
                     datadir(x), ...)
}

# Helper Functions -------------------------------------------------------------

#' Create meta information for the genes and transcripts in the ARCHS4 dataset.
#'
#' This function creates feature-level metadata for the features enumerated in
#' the `meta/genes` gene-level hdf5 file, and the `meta/transcript` transcript
#' identfiers in the transctipt-level hdf5 file for the mouse and human files
#' found in `datadir`.
#'
#' Gene symbols are the only piece of information provided for the row-level
#' identifieres for the gene count matrices. Furthermore, the gene symbol used
#' in mouse are in all uppercase, which is not how genes are referred to there.
#' In order to augment the gene symbol information with gene-level identifiers
#' and other information, we parse relatively recent GTFs provided by GENCODE.
#'
#' The fruits of the labor generated by this function are used by the
#' [archs4_feature_info()] function.
#'
#' @export
#' @seealso [archs4_feature_info()]
#'
#' @param datadir The directory that has the mouse and human expression
#'   hdf5 files. There will be `SPECIES_FEATURETYPE_augmented_info.csv.gz` files
#'   saved in this directory whe this function completes.
create_augmented_gene_info <- function(datadir = getOption("archs4.datadir")) {
  stop("TODO: Implement create_augmented_gene_info (Issue #1)")
}
