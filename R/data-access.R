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

  fn <- file.path(datadir, afiles[what])
  fe <- file.exists(fn)
  if (!all(fe)) {
    bad.query <- unique(what[!fe])
    stop("Can not find archs4 file(s) on disk: ",
         paste(bad.query, collapse = ", "))
  }

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
#' @return a tibble of series_id, sample_id, species columns
archs4_sample_table <- function() {
  map_dfr(c("human", "mouse"), function(source) {
    h5.fn <- archs4_file_path(paste0(source, "_gene"))
    tibble(series_id = h5read(h5.fn, "meta/Sample_series_id"),
           sample_id = h5read(h5.fn, "meta/Sample_geo_accession"),
           sample_h5idx = seq(sample_id),
           organism = source)
  })
}

#' Retrieves the feature (gene/transcript) information for the archs4 data
#'
#' Unfortunately only the gene symbols are stored in the gene-level count data
#' in `"meta/genes"`, so we augment the gene information by creating tables
#' of gene information found in `extdata/feature-annotation`.
#'
#' @section Augmented Gene Information:
#' Unfortunately the gene symbol is the only piece of information provided for
#' the row-level in the gene count matrices. Furthermore, the gene symbol used
#' in mouse are in all uppercase, which is not how genes are referred to there.
#' In order to augment the gene symbol information with gene-level identifiers
#' and other information, we parse relatively recent GTFs provided by GENCODE.
#'
#' The `extdata/feature-annotation/create-gene-annotaiton.R` script parses
#' relatively recent versions of the GENCODE-basic annotations and combines
#' these with biomaRt queries to get gene-level meta informations for the
#' annotated genes in these assemblies.
#'
#' The gene symbols from both the ARCHS4 h5 matrix and our gene-level metadata
#' files are transformed to lowercase and joined to harmonize these results.
#' This may likely be wrong!
#'
#' @importFrom readr read_csv
#' @export
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
                    h5idx = seq(symbol),
                    join = tolower(symbol))

    meta <- local({
      fn <- sprintf("gencode.%s.basic.annotation.csv.gz",
                    if (source == "human") "v27" else "vM16")
      fn <- system.file("extdata", "feature-annotation", fn,
                        package = "archs4", mustWork = TRUE)
      coltypes <- "ccicc"

      out <- readr::read_csv(fn, col_types = coltypes)
      out[["join"]] <- tolower(out[["symbol"]])

      if (distinct_symbol) {
        # there are duplicate entries by symbol. here I pick one by arranging
        # by symbol and entrez_id, this puts NA entrez_ids last
        out <- arrange(out, join, entrez_id)
        out <- distinct(out, join, .keep_all = TRUE)
      }

      out
    })

    ainfo <- left_join(select(ainfo, -symbol), meta, by = "join") %>%
      select(-join) %>%
      select(symbol, ensembl_gene_id, entrez_id, gene_biotype, everything())
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
#' @return a tibble of series_id, sample_id, sample_h5idx, sample_title, and
#'   sample_name columns. If the query sample or series query can't be found,
#'   then there will be an `NA` value for these columns. The `query_type` column
#'   will indicat whether the row was returned from querying by series or by
#'   sample.
archs4_sample_info <- function(id, with_description = FALSE, ...) {
  input <- geo_id_type(id)
  bad.id <- filter(input, type == "unknown")
  if (nrow(bad.id)) {
    stop("Malformed identifiers in query: ", paste(bad.id$id, collapse = ", "))
  }

  asi <- archs4_sample_table()
  series <- input %>%
    filter(type == "series") %>%
    rename(series_id = id) %>%
    left_join(asi, by = "series_id")
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
    left_join(asi, by = "sample_id")
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
    select(series_id, sample_id, query_type = type, sample_h5idx, organism)

  # query$h5.fn <- archs4_file_path(paste0(query$organism, "_gene"))

  res <- query %>%
    group_by(organism) %>%
    do({
      h5.fn <- archs4_file_path(paste0(.$organism[1L], "_gene"))
      i <- list(.$sample_h5idx)
      .$sample_title <- h5read(h5.fn, "meta/Sample_title", i)
      .$sample_name <- h5read(h5.fn, "meta/Sample_source_name_ch1", i)
      if (with_description) {
        .$sample_description <- h5read(h5.fn, "meta/Sample_description", i)
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
