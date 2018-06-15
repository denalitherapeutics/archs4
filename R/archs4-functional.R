# These functions define a "high-level" interface to working with the data
# downloaded into an ARCHS4 Data Directory. They are the workhorses of this
# package.
#
# Although users can interact with the ARCHS4 data directory "directly" using
# these functions, most will find it more convenient to interact with ARCHS4
# using an `Archs4Repository`. Many of the functions here are also available
# against an `Archs4Repository` object, and are named by stripping the archs4_
# prefix here.
#
# For instance, the following are equivalent:
#
#     R> archs4_feature_info(datadir)
#
#   and
#
#     R> a4 <- Archs4Repository(datadir)
#     R> feature_info(a4)
#
# Interacting with these data via the `Archs4Repository` should also better
# future-proof your code.


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
archs4_feature_info <- function(feature_type = "gene", source = "human",
                                augmented = TRUE,
                                datadir = getOption("archs4.datadir"), ...) {
  assert_choice(feature_type, c("gene", "transcript"))
  assert_choice(source, archs4_sources(datadir))
  assert_flag(augmented)

  fkey <- paste(source, feature_type, sep = "_")
  h5.fn <- archs4_file_path(fkey, datadir = datadir)

  if (feature_type == "gene") {
    # I am using `a4name` instead of symbol, because the values stored here
    # aren't universally/technically symbols. In the mouse dataset, the "symbol"
    # names are all uppercase, which is a non-canonical something.
    ainfo <- tibble(a4name = rhdf5::h5read(h5.fn, "meta/genes"),
                    entrez_id = rhdf5::h5read(h5.fn, "meta/gene_entrezid"),
                    h5idx = seq(a4name))
    if (source == "mouse") {
      ainfo$ens_id <- as.character(rhdf5::h5read(h5.fn, "meta/gene_ensemblid"))
    }
  } else {
    ainfo <- tibble(ensembl_id_full = rhdf5::h5read(h5.fn, "meta/transcript"),
                    ensembl_id = sub("\\.\\d+$", "", ensembl_id_full),
                    length = rhdf5::h5read(h5.fn, "meta/transcriptlength"),
                    h5idx = seq(ensembl_id_full))
  }

  if (augmented) {
    aug.fn <- paste(source, feature_type, "info", sep = "_")
    aug.fn <- archs4_file_path(aug.fn, datadir = datadir)

    if (feature_type == "gene") {
      if (source == "human") {
        join <- "a4name"
        coltypes <- "cicccciic"
      } else {
        join <- "ens_id"
        coltypes <- "cciccccciic"
      }
      meta <- readr::read_csv(aug.fn, col_types = coltypes)
      meta <- rename(meta, ensembl_id = "gene_id")
      meta <- mutate(meta, feature_type = "gene")
    } else {
      join <- "ensembl_id_full"
      coltypes <- "cciicccccciic"
      meta <- readr::read_csv(aug.fn, col_types = coltypes)
      meta <- mutate(meta, feature_type = "transcript")
    }
    # remove duplicate columns in meta table except for join column
    meta <- meta[, !colnames(meta) %in% setdiff(colnames(ainfo), join)]
    tmp <- left_join(ainfo, meta, by = join)
    stopifnot(
      nrow(tmp) == nrow(ainfo),
      all(tmp[[join]] == ainfo[[join]]))
    ainfo <- tmp
  }

  ainfo$source <- source
  ainfo
}

#' Retrieves a table of files that back an Archs4Repository
#'
#' @description
#' A yaml iskept in the Archs4 data directory (`getOption("archs4.datadir")`)
#' that links keys, (ie. `mouse_gene`) to the name of the file in the directory.
#' This abstraction is introduced so that the version of these files can be
#' updated with new downloads, and the user simply has to modify the yaml file
#' so that they are used downstream
#'
#' Reference the "ARCHS4 Data Download" section in the vignette for more
#' information.
#'
#' @export
#' @param datadir the directory that stores the ARCHS4 data files
archs4_file_info <- function(datadir = getOption("archs4.datadir")) {
  yml <- archs4_meta(datadir = datadir)
  finfo <- yml[["files"]]

  take <- function(l, wut) {
    val <- l[[wut]]
    if (is.null(val)) NA_character_ else val
  }

  tibble(
    key = names(finfo),
    source = sapply(finfo, take, "source"),
    name = sapply(finfo, take, "name"),
    url = sapply(finfo, take, "url"),
    desription = sapply(finfo, take, "description"),
    file_path = file.path(datadir, name),
    file_exists = file.exists(file_path))
}

#' Identify the file path on the system for specific ARCHS4 files.
#'
#' By default, this function will throw an error if a file does not exist
#' upon lookup. To return the *expected* path to the, even if it does not
#' exist on the file system, set `stop_if_missing = FALSE`.
#'
#' @export
#' @importFrom stats setNames
#'
#' @param key the lookup key for the file, ie. `"human_gene"` or `"mouse_gene"`.
#'   The known keys are enumerated in `archs4_file_info()$key` column.
#' @param stop_if_missing defaults to `TRUE`, which causes this function to
#'   throw an error if the file does not exist at the expected `file_path`.
#'   Set this to `FALSE` to simply raise a warning
#' @param na_missing by default, we set paths to files that don't exist to
#'   `NA`. Set this to `FALSE` to retrieve the expected path of the missing
#'   file.
#' @param file_info the output from [archs4_file_info()], which enumerates the
#'   files used by the `Archs4Repository`.
#' @param datadir the directory that stores the ARCHS4 data files
#' @return a named (by `key`) character vector of paths to the filesystem that
#'   correspond to the entries in `key`.
archs4_file_path <- function(key, stop_if_missing = TRUE,  na_missing = TRUE,
                             file_info = archs4_file_info(datadir),
                             datadir = getOption("archs4.datadir")) {
  assert_character(key, min.len = 1L)
  assert_directory(datadir, 'r')
  assert_class(file_info, "data.frame")
  assert_names(c("key", "source", "file_path"),
               subset.of = colnames(file_info))
  query <- tibble(key = key) %>% left_join(file_info, by = "key")
  qbad <- filter(query, is.na(source))

  if (nrow(qbad)) {
    bad.key <- unique(qbad[["key"]])
    stop("Unknown file queries: ", paste(bad.key, collapse = ", "))
  }

  qmiss <- filter(query, !file_exists)
  if (nrow(qmiss)) {
    bad.key <- unique(qmiss[["key"]])
    msg <- paste("Can not find archs4 file(s) on disk: ",
                 paste(bad.key, collapse = ", "))
    if (stop_if_missing) stop(msg) else warning(msg)
    if (na_missing) {
      query <- mutate(query, file_path = ifelse(file_exists, file_path, NA))
    }
  }

  setNames(query[["file_path"]], query[["key"]])
}

#' Retrieves tibble of sample-level covariates available in mouse and human data.
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
#' @noRd
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
  h5.fn <- archs4_file_path(file, datadir = datadir)
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

#' Checks which samples from a single series are present/absent.
#'
#' Due to download or alignment issues, the ARCHS4 data processing pipeline may
#' not include all of the samples included in a particular GEO series. This
#' function will return a table with an `in_archs4` columns that indicates
#' whether or not a sample from a particular series is present in ARCHS4.
#'
#' @export
#'
#' @param id a single GEO series id, ie. `"GSEnnnnn"`
#' @param sample_table the output from [archs4_sample_table()], which lists
#'   the series_id,sample_id combinations found in the ARCHS4 repository.
#' @param datadir the directory that holds the archs4 data
#' @return a tibble of information for a series.
#' @examples
#' info <- archs4_series_status("GSE89189")
archs4_series_status <- function(id,
                                 sample_table = archs4_sample_table(datadir = datadir),
                                 datadir = getOption("archs4.datadir"),
                                 ...) {
  assert_string(id)
  if (geo_id_type(id)[["type"]] != "series") {
    stop("A GEO series id is required (ie. \"GSEnnnnn\")")
  }

  gq <- lookup_gse(id)
  geo <- tibble(series_id = id, sample_id = gq[["sample"]])

  archs <- filter(sample_table, series_id == id)
  archs <- transmute(archs, sample_id, in_archs4 = TRUE)

  out <- left_join(geo, archs, by = "sample_id")
  out[["in_archs4"]] <- !is.na(out[["in_archs4"]])
  out
}


#' Retrieves information for samples by GSE series or sample IDs
#'
#' @description
#' Fetch a tibble of series and sample information by querying the arcsh4
#' dataset by GEO sample (GSE) or sample (GSM) ids.
#'
#' For each unique GEO series identifier (`"GSEnnnn"`), we will check if the
#' ARCHS4 dataset is missing any of its samples when `check_missing_samples`
#' is set to `TRUE` (default).
#'
#' @export
#'
#' @param id a character vector of GEO series or sample ids.
#' @param columns the names of the sample metadata columns desired. This
#'   defaults to `c("Sample_title", "Sample_source_name_ch1")`. The values
#'   in `columns` must be a subset of the values enumerated in
#'   [archs4_sample_covariates()].
#' @param sample_table the output from [archs4_sample_table()], which lists
#'   the series_id,sample_id combinations found in the ARCHS4 repository.
#' @param sample_covariates the `data.frame`-definition of the sample covariates
#'   found in the ARCHS4 datasetes, which is constructed via a call to
#'   [archs4_sample_covariates()]. The parameter is included in here so that
#'   a cached version of this `data.frame` can be re-used.
#' @param check_missing_samples When `TRUE` (the default), this function will
#'   check every unique GEO series identifier (`"GSEnnnn"`) for missing samples
#'   by using an NCBI Rest service via a call to [archs4_series_status()],
#'   and [lookup_gse()].
#' @param datadir the directory that holds the archs4 data
#' @return a tibble of series_id, sample_id, sample_h5idx, sample_title, and
#'   sample_name columns. If the query sample or series query can't be found,
#'   then there will be an `NA` value for these columns. The `query_type` column
#'   will indicat whether the row was returned from querying by series or by
#'   sample.
#' @examples
#' si <- archs4_sample_info("GSE52564") # ben barres transcriptome db ...
archs4_sample_info <- function(id,
                               columns = c("Sample_title", "Sample_source_name_ch1"),
                               sample_table = archs4_sample_table(datadir = datadir),
                               sample_covariates = archs4_sample_covariates(datadir),
                               check_missing_samples = TRUE,
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

  with.missing.samples <- NULL
  if (check_missing_samples && nrow(series)) {
    for (s in unique(series[["series_id"]])) {
      check <- archs4_series_status(s, sample_table = sample_table,
                                    datadir = datadir)
      missing <- filter(check, !in_archs4)
      if (nrow(missing)) {
        msg <- paste(
          sprintf("The %s series has %d missing samples in ARCHS4:\n    *",
                  s, nrow(missing)),
          paste(missing[["sample_id"]], collapse ="\n    * "),
          "\n")
        warning(msg, immediate. = TRUE, call. = FALSE)
        # TODO: Add `with_missing_samples` parameter, or does it make sense
        #       to just automatically append these?
        # it hurts me to do this:
        with.missing.samples <- bind_rows(with.missing.samples, missing)
      }
    }
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

  col.order <- c("series_id", "sample_id", columns)
  col.order <- c(col.order, setdiff(colnames(out), col.order))
  out[, col.order]
}

#' Helper function that implements sample-level meatdata retrieval
#'
#' This function is only meant to be called within the `do({})` block in the
#' [archs4_sample_info()] function, and as such does no argument checking and
#' is intentionally not exported.
#'
#' @noRd
#'
#' @importFrom rhdf5 h5read
.with_sample_info <- function(x, columns, organism, sample_covariates, datadir) {
  organism <- organism[1L]
  h5g.fn <- archs4_file_path(paste0(organism, "_gene"), datadir = datadir)
  h5t.fn <- archs4_file_path(paste0(organism, "_transcript"), datadir = datadir)

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
#' @importFrom readr read_csv
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
      mutate(organism_gene = NULL, organism_transcript = NULL) %>%
      select(series_id, sample_id, organism, libsize, normfactor, everything())
  } else {
    res <- map_dfr(archs4_sources(), function(source) {
      fkey <- paste(source, feature_type, sep = "_")
      h5.fn <- archs4_file_path(fkey, datadir = datadir)
      dat <- tibble(series_id = trimws(h5read(h5.fn, "meta/Sample_series_id")),
                    sample_id = trimws(h5read(h5.fn, "meta/Sample_geo_accession")),
                    sample_h5idx = seq(sample_id),
                    organism = source)
      if (feature_type == "gene") {
        # Load library size and sample factors
        sfn <- file.path(dirname(h5.fn), paste0(source, "_gene-normfactors.csv"))
        # browser()
        if (file.exists(sfn)) {
          sf <- suppressWarnings(read_csv(sfn, col_types = "cciidd"))
          sf <- select(sf, series_id, sample_id, libsize, normfactor)
          dat <- left_join(dat, sf, by = c("series_id", "sample_id"))
        } else {
          dat$libsize <- h5read(h5.fn, "meta/reads_aligned")
          dat$normfactor <- 1
        }
      }
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
#' @param datadir the directory that holds the archs4 data
#' @return a character vector listing the different sources (organisms) that
#'   the ARCHS4 repository has data for.
archs4_sources <- function(datadir = getOption("archs4.datadir")) {
  archs4_meta(datadir)[["sources"]]
}

#' Retrieves the meta information associated with an ARCHS4 datadir
#'
#' @export
#' @importFrom yaml read_yaml
#'
#' @param datadir the directory that holds the archs4 data
#' @return a list-representation of the `meta.yaml` file
archs4_meta <- function(datadir = getOption("archs4.datadir")) {
  assert_directory(datadir, "r")
  fpath <- assert_file_exists(file.path(datadir, "meta.yaml"), "r", "yaml")
  yml <- yaml::read_yaml(fpath)

  # doing some validation on this file
  assert_character(yml[["sources"]])
  assert_list(yml[["files"]])
  yml
}
