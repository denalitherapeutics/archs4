#' Create a DGEList for the expression data of a series or set of samples.
#'
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom edgeR DGEList
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
    }
    assert_data_frame(features)
    assert_subset(c("ensembl_id"), colnames(features))
    features <- distinct(features, ensembl_id, .keep_all = TRUE)
  }

  # Identify the unique samples that are being queried -------------------------
  si <- sample_info(x, id, columns = sample_columns,
                    check_missing_samples = check_missing_samples)
  si <- as.data.frame(si, strinsAsFactors = FALSE)
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
    finfo <- semi_join(finfo, features, by = "ensembl_id")
  }

  rownames(finfo) <- finfo$ensembl_id

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

  # Tack on the pre-calculated lib.size and norm.factors from the
  # `estimate_repository_norm_factors` call.
  libinfo <- sample_table(a4) %>%
    select(series_id, sample_id, lib.size = libsize, norm.factors = normfactor) %>%
    distinct(sample_id, .keep_all = TRUE)
  libinfo <- left_join(select(si, series_id, sample_id),
                       libinfo, by = c("series_id", "sample_id"))
  libinfo <- libinfo %>%
    transform(lib.size = ifelse(is.na(lib.size), mean(lib.size, na.rm=TRUE),
                                lib.size),
              norm.factors = ifelse(is.na(norm.factors),
                                    mean(norm.factors, na.rm=TRUE),
                                    norm.factors))
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

