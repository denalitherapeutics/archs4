#' Retrieves the feature (gene/transcript) information for the archs4 data
#'
#' Only the gene symbols (`meta/genes` in gene expression hd5 file) or entrez
#' transcript identifiers (`meta/ensemblid` for the transcript hdf5 file) are
#' stored in thse data. We use [create_augmented_gene_info()] function to
#' generate and store extra metadata for these features, which are then appended
#' to these identifiers with this function.
#'
#' @export
#' @importFrom readr read_csv
#' @seealso [create_augmented_gene_info()]
#'
#' @param feature_type gene or transcript?
#' @param source human or mouse
#' @return a tibble of information
archs4_feature_info <- function(feature_type = c("gene", "transcript"),
                                source = c("human", "mouse"),
                                distinct_symbol = TRUE,
                                datadir = getOption("archs4.datadir")) {
  feature_type <- match.arg(feature_type)
  source <- match.arg(source)

  if (source == "mouse" && feature_type == "transcript") {
    stop("The meta information for the mouse transcript data ",
         "(meta/ensemblid, meta/transriptlength) is incomplete. There are ",
         "only 98,492 entries in these meta entries, but 178136 tx estimates ",
         "in the data/expression file.")
  }
  if (feature_type == "transcript") {
    warning("Augmented transcript-level has not yet been generated (Issue #1)")
  }


  h5.fn <- archs4_file_path(paste(source, feature_type, sep = "_"))
  aug.fn <- paste(source, feature_type, "info", sep = "_")
  aug.fn <- archs4_file_path(aug.fn)

  if (feature_type == "gene") {
    coltypes <- "ccicc"
  } else {
    coltypes <- NULL
  }

  if (feature_type == "gene") {
    ainfo <- tibble(symbol = rhdf5::h5read(h5.fn, "meta/genes"),
                    h5idx = seq(symbol),
                    join = tolower(symbol))

    meta <- readr::read_csv(aug.fn, col_types = coltypes)
    meta[["join"]] <- tolower(meta[["symbol"]])
    if (distinct_symbol) {
      # there are duplicate entries by symbol. here I pick one by arranging
      # by symbol and entrez_id, this puts NA entrez_ids last
      meta <- arrange(meta, join, entrez_id)
      meta <- distinct(meta, join, .keep_all = TRUE)
    }

    ainfo <- left_join(select(ainfo, -symbol), meta, by = "join") %>%
      select(-join) %>%
      select(symbol, ensembl_gene_id, entrez_id, gene_biotype, everything())
  } else {
    ainfo <- tibble(ensembl_id_full = rhdf5::h5read(h5.fn, "meta/ensemblid"),
                    ensembl_id = sub("\\.\\d+$", "", ensembl_id_full),
                    h5idx = seq(ensembl_id_full))
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
    human_transcript = list(fn='human_hiseq_eid_1.0.h5'),
    human_transcript_info = list(fn='human_transcript_augmented_info.csv.gz'),
    # mouse
    mouse_gene = list(fn = 'mouse_matrix.h5'),
    mouse_gene_info = list(fn = 'mouse_gene_augmented_info.csv.gz'),
    mouse_transcript = list(fn='mouse_hiseq_eid_1.0.h5'),
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
    res <- map_dfr(c("human", "mouse"), function(source) {
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
#' @return a tibble of series_id, sample_id, sample_h5idx, sample_title, and
#'   sample_name columns. If the query sample or series query can't be found,
#'   then there will be an `NA` value for these columns. The `query_type` column
#'   will indicat whether the row was returned from querying by series or by
#'   sample.
archs4_sample_info <- function(id, sample_table = NULL,
                               with_description = FALSE,
                               datadir = getOption("archs4.datadir"), ...) {
  input <- geo_id_type(id)
  bad.id <- filter(input, type == "unknown")
  if (nrow(bad.id)) {
    stop("Malformed identifiers in query: ", paste(bad.id$id, collapse = ", "))
  }

  if (is.null(sample_table)) {
    sample_table <- archs4_sample_table(datadir = datadir)
  }

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
  }

  query <- bind_rows(series, samples) %>%
    select(series_id, sample_id, query_type = type,
           sample_h5idx_gene, sample_h5idx_transcript, organism)

  res <- query %>%
    group_by(organism) %>%
    do({
      h5g.fn <- archs4_file_path(paste0(.$organism[1L], "_gene"), datadir)
      h5t.fn <- archs4_file_path(paste0(.$organism[1L], "_transcript"), datadir)
      stitle <- sname <- desc <- character(nrow(.))

      # get data from gene hdf5 file where !is.na(sample_h5idx_gene)
      ginfo <- !is.na(.$sample_h5idx_gene)
      if (any(ginfo)) {
        i <- list(.$sample_h5idx_gene[ginfo])
        stitle[ginfo] <- h5read(h5g.fn, "meta/Sample_title", i)
        sname[ginfo] <- h5read(h5g.fn, "meta/Sample_title", i)
        if (with_description) {
          desc[ginfo] <- h5read(h5g.fn, "meta/Sample_description", i)
        }
      }
      # Anything left to grab from transcript hdf5 file?
      tinfo <- !ginfo
      if (any(tinfo)) {
        i <- list(.$sample_h5idx_transcript[tinfo])
        stitle[tinfo] <- h5read(h5t.fn, "meta/Sample_title", i)
        sname[tinfo] <- h5read(h5t.fn, "meta/Sample_title", i)
        if (with_description) {
          desc[tinfo] <- h5read(h5t.fn, "meta/Sample_description", i)
        }

      }
      .$sample_title <- stitle
      .$sample_name <- sname
      if (with_description) {
        .$sample_description <- desc
      }
      .
    }) %>%
    ungroup

  out <- res %>%
    bind_rows(bad.series) %>%
    bind_rows(bad.samples)

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
