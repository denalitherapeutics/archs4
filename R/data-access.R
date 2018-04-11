#' Retrieves the feature (gene/transcript) information for the archs4 data
#'
#' Only the gene symbols (`meta/genes` in gene expression hd5 file) or entrez
#' transcript identifiers (`meta/transcript` for the transcript hdf5 file) are
#' stored in thse data. We use [create_augmented_feature_info()] function to
#' generate and store extra metadata for these features, which are then appended
#' to these identifiers with this function.
#'
#' @export
#' @importFrom readr read_csv
#' @seealso [create_augmented_feature_info()]
#'
#' @param feature_type gene or transcript?
#' @param source human or mouse
#' @param augmented include extra gene- or transcript-level features?
#'   Default: `TRUE`
#' @param datadir the directory that stores the ARCHS4 data files
#' @param ... pass through
#' @return a tibble of information
archs4_feature_info <- function(feature_type = c("gene", "transcript"),
                                source = archs4_sources(), augmented = TRUE,
                                datadir = getOption("archs4.datadir"), ...) {
  assert_flag(augmented)
  feature_type <- match.arg(feature_type)
  source <- match.arg(source)

  h5.fn <- archs4_file_path(paste(source, feature_type, sep = "_"))

  if (feature_type == "gene") {
    # I am using `a4name` instead of symbol, because the values stored here
    # aren't universally/technically symbols. In the mouse dataset, the "symbol"
    # names are all uppercase, which is a non-canonical something.
    ainfo <- tibble(a4name = rhdf5::h5read(h5.fn, "meta/genes"),
                    h5idx = seq(a4name))
  } else {
    ainfo <- tibble(ensembl_id_full = rhdf5::h5read(h5.fn, "meta/transcript"),
                    ensembl_id = sub("\\.\\d+$", "", ensembl_id_full),
                    length = rhdf5::h5read(h5.fn, "meta/transcriptlength"),
                    h5idx = seq(ensembl_id_full))
  }

  if (augmented) {
    aug.fn <- paste(source, feature_type, "info", sep = "_")
    aug.fn <- archs4_file_path(aug.fn)

    if (feature_type == "gene") {
      join <- "a4name"
      coltypes <- "cicccciic"
      meta <- readr::read_csv(aug.fn, col_types = coltypes)
    } else {
      join <- "ensembl_id_full"
      coltypes <- "cciicccccciic"
      meta <- readr::read_csv(aug.fn, col_types = coltypes)
    }
    # remove duplicate columns in meta table except for join column
    meta <- meta[, !colnames(meta) %in% setdiff(colnames(ainfo), join)]
    tmp <- left_join(ainfo, meta, by = join)
    stopifnot(
      nrow(tmp) == nrow(ainfo),
      all(tmp[[join]] == ainfo[[join]]))
    ainfo <- tmp
  }

  ainfo
}

archs4.files <- function() {
  # files <- c(
  #   human_gene = file.path('human_matrix.h5'),
  #   mouse_gene = file.path('mouse_matrix.h5'))

  files <- list(
    # human
    human_gene = list(fn = 'human_matrix.h5'),
    human_gene_info = list(fn = 'human_gene_augmented_info.csv.gz'),
    # human_transcript = list(fn='human_hiseq_eid_1.0.h5'),
    human_transcript = list(fn='human_hiseq_transcript_v2.h5'),
    human_transcript_info = list(fn='human_transcript_augmented_info.csv.gz'),
    # mouse
    mouse_gene = list(fn = 'mouse_matrix.h5'),
    mouse_gene_info = list(fn = 'mouse_gene_augmented_info.csv.gz'),
    # mouse_transcript = list(fn='mouse_hiseq_eid_1.0.h5'),
    mouse_transcript = list(fn='mouse_hiseq_transcript_v2.h5'),
    mouse_transcript_info = list(fn='mouse_transcript_augmented_info.csv.gz'))

  files
}

#' Identify the file path on the system for specific ARCHS4 files.
#'
#' @export
#'
#' @param what what type of file? `"human_gene"` or `"mouse_gene"`. This can
#'   be vectorized.
#' @return the path on the filesystem for the path asked for
archs4_file_path <- function(what, datadir = getOption("archs4.datadir")) {
  afiles <- archs4.files()
  assert_directory(datadir, 'r')
  bad.what <- setdiff(what, names(afiles))
  if (length(bad.what)) {
    stop("Unknown file queries: ", paste(bad.what, collapse = ", "))
  }

  afiles <- sapply(afiles, '[[', 'fn')
  fn <- file.path(datadir, afiles[what])
  fe <- file.exists(fn)
  if (!all(fe)) {
    bad.query <- unique(what[!fe])
    stop("Can not find archs4 file(s) on disk: ",
         paste(bad.query, collapse = ", "))
  }
  names(fn) <- what
  fn
}

#' Retrieves a tibble of sample-level covariates available in mouse and human data.
#'
#' Enumerate the sample covariates available in mouse and human data.
#' Note that the covariates available in human and mouse are the same
#' between the gene and transcript level files
#'
#' @export
#' @param datadir the directory that holds the archs4 data
archs4_sample_covariates <- function(datadir = getOption("archs4.datadir"),
                                     ...) {
  mouse.covs <- .sample_metadata_files("mouse_gene", datadir) %>%
    select(name) %>%
    mutate(mouse = TRUE)
  human.covs <- .sample_metadata_files("human_gene", datadir) %>%
    select(name) %>%
    mutate(human = TRUE)
  sample.covs <- mouse.covs %>%
    full_join(human.covs, by = "name") %>%
    mutate(mouse = ifelse(is.na(mouse), FALSE, TRUE),
           human = ifelse(is.na(human), FALSE, TRUE))
  sample.covs
}

#' Enumerates internal files from hdf5 binary that contain sample metadata
#'
#' This function is intentionally not exported
#'
#' @importFrom rhdf5 h5ls
#'
#' @param file hdf5 file in `datadar` to use to identify the internal "metal"
#'   files that correspond to sample-level metadata. You shouldn't need to
#'   change this as all hdf5 data files should carry the same metadata. This is
#'   here mainly for debugging purposes.
#' @param datadir the directory that holds the archs4 data
#' @param ... pass through
#' @return a vector of sample metadata names that are stored in archs4
.sample_metadata_files <- function(file = "mouse_gene",
                                   datadir = getOption("archs4.datadir"),
                                   ...) {
  h5.fn <- archs4_file_path(file)
  all.files <- as_tibble(h5ls(h5.fn))
  all.files[["idim"]] <- suppressWarnings(as.integer(all.files[["dim"]]))

  dat.dim <- local({
    dd <- filter(all.files, group == "/data", name == "expression")
    if (nrow(dd) != 1L) {
      stop("The hdf5 file is not an archs4 expression matrix ",
           sprintf("[%s]", h5.fn))
    }
    dd <- tolower(dd$dim)
    out <- sapply(strsplit(dd, " *x *")[[1L]], as.integer)
    setNames(out, c("rows", "columns"))
  })
  nsamples <- dat.dim["columns"]
  nfeats <- dat.dim["rows"]

  sample.meta <- filter(all.files, group == "/meta", idim == nsamples)
  mutate(sample.meta, h5name = paste(group, name, sep = "/"))
}

#' Retrieves information for samples by GSE series or sample IDs
#'
#' Fetch a tibble of series and sample information by querying the arcsh4
#' dataset by GEO sample (GSE) or sample (GSM) ids.
#'
#' @export
#'
#' @param id a character vector of GEO series or sample ids.
#' @param columns the names of the sample metadata columns desired. This
#'   defaults to `c("Sample_title", "Sample_source_name_ch1")`. The values
#'   in `columns` must be a subset of the values enumerated in
#'   [archs4_sample_metadata_names()].
#' @param sample_table the output from [archs4_sample_table()], which lists
#'   the series_id,sample_id combinations found in the ARCHS4 repository.
#' @param datadir the directory that holds the archs4 data
#' @return a tibble of series_id, sample_id, sample_h5idx, sample_title, and
#'   sample_name columns. If the query sample or series query can't be found,
#'   then there will be an `NA` value for these columns. The `query_type` column
#'   will indicat whether the row was returned from querying by series or by
#'   sample.
archs4_sample_info <- function(id,
                               columns = c("Sample_title", "Sample_source_name_ch1"),
                               sample_table = archs4_sample_table(datadir = datadir),
                               sample_covariates = archs4_sample_covariates(datadir),
                               datadir = getOption("archs4.datadir"), ...) {
  input <- geo_id_type(id)
  bad.id <- filter(input, type == "unknown")
  if (nrow(bad.id)) {
    stop("Malformed identifiers in query: ", paste(bad.id$id, collapse = ", "))
  }

  # check sample metadata columns
  columns <- unique(columns) %>%
    assert_character(any.missing = FALSE, min.len = 1L) %>%
    assert_subset(sample_covariates$name)

  # This hurts: I'm doing this to ensure that queries to this function where
  # all `id`s are not found in the ARCHS4 repository still return a tibble
  # of the right dimensions, but has NAs in most places
  dummy <- tibble(
    series_id = character(),
    sample_id = character(),
    query_type = character(),
    sample_h5idx_gene = integer(),
    sample_h5idx_transcript = integer(),
    organism = character())
  for (cname in columns) dummy[[cname]] <- character()

  series <- input %>%
    filter(type == "series") %>%
    rename(series_id = id) %>%
    left_join(sample_table, by = "series_id")
  bad.series <- series %>%
    filter(is.na(organism)) %>%
    distinct(series_id) %>%
    mutate(query_type = "series")
  if (nrow(bad.series)) {
    warning(nrow(bad.series), " series identifiers not found",
            immediate. = TRUE)
    series <- anti_join(series, bad.series, by = "series_id")
    bad.series <- bind_rows(bad.series, dummy)
  }

  samples <- input %>%
    filter(type == "sample") %>%
    rename(sample_id = id) %>%
    left_join(sample_table, by = "sample_id")
  bad.samples <- samples %>%
    filter(is.na(organism)) %>%
    distinct(sample_id) %>%
    mutate(query_type = "sample")
  if (nrow(bad.samples)) {
    warning(nrow(bad.samples), " sample identifiers not found",
            immediate. = TRUE)
    samples <- anti_join(samples, bad.samples, by = "sample_id")
    bad.samples <- bind_rows(bad.samples, dummy)
  }

  query <- series %>%
    bind_rows(samples) %>%
    select(series_id, sample_id, query_type = type,
           sample_h5idx_gene, sample_h5idx_transcript, organism)

  # Only perform this if there >= 1 series or sample identifiers were found
  if (nrow(query)) {
    res <- query %>%
      group_by(organism) %>%
      do({
        .with_sample_info(., columns, .$organism[1L], sample_covariates, datadir)
      }) %>%
      ungroup
  } else {
    res <- query
  }

  out <- res %>%
    bind_rows(bad.series) %>%
    bind_rows(bad.samples)

  out
}

#' Helper function that implements sample-level meatdata retrieval
#'
#' This function is only meant to be called within the `do({})` block in the
#' [archs4_sample_info()] function, and as such does no argument checking and
#' is intentionally not exported.
#'
#' @importFrom rhdf5 h5read
.with_sample_info <- function(x, columns, organism, sample_covariates, datadir) {
  organism <- organism[1L]
  h5g.fn <- archs4_file_path(paste0(organism, "_gene"), datadir)
  h5t.fn <- archs4_file_path(paste0(organism, "_transcript"), datadir)

  out <- sapply(columns, function(i) {
    rep(NA_character_, nrow(x))
  }, simplify = FALSE)

  # ensure that we only lookup covariates that are available for this organism
  scovs <- sample_covariates[["name"]][sample_covariates[[organism]]]
  columns <- intersect(columns, scovs)

  # we first try to get metadata from the gene matrix
  ginfo <- !is.na(x$sample_h5idx_gene)
  if (any(ginfo)) {
    i <- list(x$sample_h5idx_gene[ginfo])
    for (what in columns) {
      out[[what]][ginfo] <- h5read(h5g.fn, paste0("meta/", what), i)
    }
  }

  # anything left to fetch from the transcript hdf5?
  tinfo <- !ginfo
  if (any(tinfo)) {
    i <- list(x$sample_h5idx_transcript[tinfo])
    for (what in columns) {
      out[[what]][tinfo] <- h5read(h5t.fn, paste0("meta/", what), i)
    }
  }

  res <- bind_cols(out)
  bind_cols(x, res)
}

#' Lists the GEO series and samples available in the human and mouse datasets
#'
#' This function queries the human and mouse gene expression matrices from
#' the ARCHS4 data release and combines their GEO series and sample identifiers
#' into a long table, annotated with the organism the samples come from.
#'
#' This function executes very quickly (less thatn 0.10th of a second), so most
#' sample-level query functions in this package which you would think would
#' benefit from specifying human/mouse don't have to, as they will join into
#' this table to find out what species is being queried.
#'
#' @export
#' @param feature_type currently, the `"gene"` and `"transcript"` datasets
#'   are not the same.
#' @param unroll_series There are some malformed series identifiers, like
#'   `"GSE36025Xx-xXGSE49417Xx-xXGSE49847"` when the same sample_id appears
#'   in multiple series. When this is `TRUE` (default), the series_id's are
#'   unrolled and cleaned up.
#' @param datadir the directory that holds the archs4 data
#' @return a tibble of series_id, sample_id, species columns
archs4_sample_table <- function(feature_type = c("all", "gene", "transcript"),
                                unroll_series = TRUE,
                                datadir = getOption("archs4.datadir")) {
  feature_type <- match.arg(feature_type)

  if (feature_type == "all") {
    gst <- archs4_sample_table("gene", unroll_series, datadir)
    tst <- archs4_sample_table("transcript", unroll_series, datadir)
    res <- gst %>%
      full_join(tst, by = c("series_id", "sample_id"),
                suffix = c("_gene", "_transcript")) %>%
      mutate(organism = ifelse(is.na(organism_gene), organism_transcript,
                               organism_gene)) %>%
      mutate(organism_gene = NULL, organism_transcript = NULL)
  } else {
    res <- map_dfr(archs4_sources(), function(source) {
      h5.fn <- archs4_file_path(paste(source, feature_type, sep = "_"), datadir)
      dat <- tibble(series_id = trimws(h5read(h5.fn, "meta/Sample_series_id")),
                    sample_id = trimws(h5read(h5.fn, "meta/Sample_geo_accession")),
                    sample_h5idx = seq(sample_id),
                    organism = source)
      if (unroll_series) {
        regex <- "[^GSEM0-9]+"
        dat <- separate_rows(dat, series_id, sep = regex)
      }
    })
  }
  res
}

#' Lists the different sources ARCHS4 is built for
#'
#' We hardocde these values in a lot of places ... who knows if one day these
#' are updated?
#'
#' @export
#' @return a character vector listing the different sources (organisms) that
#'   the ARCHS4 repository has data for.
archs4_sources <- function() {
  c("mouse", "human")
}
