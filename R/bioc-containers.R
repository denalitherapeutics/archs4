#' Create a DGEList for the expression data of a series or set of samples.
#'
#'
#'
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom edgeR DGEList
#' @seealso [fetch_expression()]
#'
#' @param id a vector of series or sample id's.
#' @param features a feature-descriptor of the features you want to include
#'   counts for
#' @param sample_columns the names of the sample covariates that are stored
#'   in the ARCHS4 Dataset; a complete list of what covariates are available
#'   in the ARCHS4 dataset is found using the [archs4_sample_covariates()]
#'   function.
#' @param feature_type do you want `"gene"` or `"transcript"` level expression?
#' @param row_id either `"ensembl"` or `"symbol"`. If this is `"ensembl"` and
#'   `feature_type == "transcript"`, then we remove the rows from the count
#'   dataset that we couldn't map symbol -> ensembl_gene_id for.
#' @return a `DGEList` of results
#'
#' @examples
#' a4 <- Archs4Repository()
#' y <- as.DGEList(a4, "GSE89189", feature_type = "gene")
as.DGEList <- function(x, id, features = NULL,
                       sample_columns = c("Sample_title", "Sample_source_name_ch1"),
                       feature_type = c("gene", "transcript"),
                       row_id = c("ensembl", "symbol"),
                       check_missing_samples = TRUE) {
  assert_class(x, "Archs4Repository")
  feature_type <- match.arg(feature_type)
  row_id <- match.arg(row_id)

  if (!is.null(features)) {
    if (is.character(features)) {
      features <- feature_lookup(x, features, type = type)
      if (any(is.na(features$a4name))) {
        warning("Removing 'not found' features from query")
        features <- filter(features, !is.na(a4name))
      }
    }
    assert_data_frame(features, min.rows = 1L)
    assert_subset(c("ensembl_id", "a4name"), colnames(features))
    # features <- distinct(features, ensembl_id, .keep_all = TRUE)
    # ensembl_id is something of a "second class citizen". The a4name is
    # what has been there in the ARCHS4 data since "the beginning"
    features <- distinct(features, ensembl_id, .keep_all = TRUE)
  }

  # Identify the unique samples that are being queried -------------------------
  si <- sample_info(x, id, columns = sample_columns,
                    check_missing_samples = check_missing_samples)
  si <- as.data.frame(si, stringsAsFactors = FALSE)
  si <- distinct(si, sample_id, .keep_all = TRUE)
  rownames(si) <- si$sample_id

  # Check for identifiers that were not found
  not.found <- filter(si, is.na(organism))
  if (nrow(not.found)) {
    stop("The following samples could not be found: ",
         paste(sprintf("%s::%s", not.found$series_id, not.found$sample_id),
               collapse = "; "))
  }

  # Check that user is asking for samples that are all from the same species
  org <- unique(si$organism)
  if (length(org) != 1L) {
    stop("You are querying across species")
  }

  # Fetch feature meta information ---------------------------------------------
  finfo <- feature_info(x, feature_type, org, augmented = TRUE)
  finfo <- as.data.frame(finfo, stringsAsFactors = FALSE)
  if (feature_type == "gene") {
    if (row_id == "ensembl") {
      finfo <- filter(finfo, !is.na(ensembl_id))
      rownames(finfo) <- finfo[["ensembl_id"]]
    } else {
      rownames(finfo) <- finfo[["a4name"]]
    }
  } else {
    dup.ensid <- duplicated(finfo[["ensembl_id"]])
    if (any(dup.ensid)) {
      warning("Duplicated ensembl identifiers when version is removed, ",
              "rownames maintain their versioned id")
      rownames(finfo) <- finfo[["ensembl_id_full"]]
    } else {
      rownames(finfo) <- finfo[["ensembl_id"]]
    }
  }

  if (!is.null(features)) {
    # cf. note about ensembl_id being a second class citizen
    # finfo <- semi_join(finfo, features, by = "ensembl_id")
    # finfo <- semi_join(finfo, features, by = "a4name")
    # semi_join (and filter(?)) strips rownames!
    finfo <- subset(finfo, a4name %in% features[["a4name"]])
  }

  # Fetch the count data for the given samples ---------------------------------
  h5col <- paste0("sample_h5idx_", feature_type)
  isna <- is.na(si[[h5col]])
  if (any(isna)) {
    msg <- paste0(
      "The following samples do not have ", feature_type,
      "-level quantitation and will not be included in the expression ",
      "container:\n  ",
      paste(sprintf("%s::%s", si$series_id[isna], si$sample_id[isna]),
            collapse = "; "))
    warning(msg, immediate. = TRUE)
    si <- si[!isna,,drop = FALSE]
  }

  if (nrow(si) == 0L) {
    stop("No samples left to assemble expression data")
  }

  # Pull out pre-calcuated lib.size and norm.factor values. This has the
  # added benefit of pre-emptively identifying samples with weird (huge)
  # numbers that turn to NAs due to integer rollover, ie. we get random
  # errors of this variety when reading from the HDF5 file:
  #
  #   integer value -2^63 replaced NA. See the section
  #   'Large integer data types' in the 'rhdf5' vignette
  #   for more details.
  #
  # This happens when we fish data out from the hdf5 file, but also we
  # ran into this problem when running the `estimate_repository_norm_factors`
  # function when generating these lib.size and norm.factor values
  libinfo <- sample_table(a4) %>%
    select(series_id, sample_id, lib.size = libsize, norm.factors = normfactor,
           a4libsize) %>%
    distinct(sample_id, .keep_all = TRUE)
  libinfo <- select(si, series_id, sample_id) %>%
    left_join(libinfo, by = c("series_id", "sample_id"))
  na.overflow <- is.na(libinfo$lib.size) | is.na(libinfo$norm.factors)
  if (any(na.overflow)) {
    warning("Removing ", sum(na.overflow), " samples due to libsize NA overflow issues")
    libinfo <- filter(libinfo, !na.overflow)
    si <- subset(si, sample_id %in% libinfo$sample_id)
  }

  # avg.lsize <- mean(libinfo$lib.size, na.rm = TRUE)
  # avg.nfact <- mean(libinfo$norm.factors, na.rm = TRUE)
  # libinfo <- libinfo %>%
  #   transform(lib.size = ifelse(is.na(lib.size), a4libsize, lib.size),
  #             norm.factors = ifelse(is.na(norm.factors), avg.nfact, norm.factors))

  rownames(si) <- si[["sample_id"]]

  counts <- local({
    key <- paste(org, feature_type, sep = "_")
    h5.fn <- file_path(x, key)
    index <- list(finfo$h5idx, si[[h5col]])
    cnts <- rhdf5::h5read(h5.fn, "data/expression", index=index)
    colnames(cnts) <- rownames(si)
    # cnts <- cnts[finfo$h5idx,,drop = FALSE]
    rownames(cnts) <- rownames(finfo)
    cnts
  })

  out <- suppressWarnings(edgeR::DGEList(counts, genes = finfo, samples = si))
  xref <- match(colnames(out), libinfo$sample_id)
  if (any(is.na(xref))) {
    stop("Problem matching sample_id to libinfo data.frame")
  }
  if (!all(colnames(out) == libinfo$sample_id[xref])) {
    stop("Mismatch in outgoing DGEList to libinfo data.frame")
  }
  out$samples$lib.size <- libinfo$lib.size[xref]
  out$samples$norm.factors <- libinfo$norm.factors[xref]
  out
}

