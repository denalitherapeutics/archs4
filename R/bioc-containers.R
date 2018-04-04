#' Create a DGEList for the expression data of a series or set of samples.
#'
#' @export
#'
#' @param id a vector of series or sample id's.
#' @param feature_type do you want `"gene"` or `"transcript"` level expression?
#' @param row_id either `"ensembl"` or `"symbol"`. If this is `"ensembl"` and
#'   `feature_type == "transcript"`, then we remove the rows from the count
#'   dataset that we couldn't map symbol -> ensembl_gene_id for.
#' @return a `DGEList` of results
#'
#' @examples
#' y <- as.DGEList("GSE89189", feature_type = "gene", source = "human")
as.DGEList <- function(x, id,
                       sample_columns = c("Sample_title", "Sample_source_name_ch1"),
                       feature_type = c("gene", "transcript"),
                       row_id = c("ensembl", "symbol")) {
  assert_class(x, "Archs4Repository")
  feature_type <- match.arg(feature_type)
  row_id <- match.arg(row_id)

  # Identify the unique samples that are being queried -------------------------
  si <- sample_info(x, id, columns = sample_columns)
  si <- as.data.frame(si, strinsAsFactors = FALSE)
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
  finfo <- archs4_feature_info(feature_type, org)
  finfo <- as.data.frame(finfo, stringsAsFactors = FALSE)
  finfo <- filter(finfo, !is.na(h5idx))
  if (feature_type == "gene") {
    if (row_id == "ensembl") {
      finfo <- filter(finfo, !is.na(ensembl_gene_id))
      rownames(finfo) <- finfo[["ensembl_gene_id"]]
    } else {
      rownames(finfo) <- finfo[["symbol"]]
    }
  } else {
    dup.ensid <- duplicated(finfo[["ensembl_id_full"]])
    if (any(dup.ensid)) {
      warning("Duplicated ensembl identifiers when version is removed, ",
              "rownames maintain their versioned id")
      rownames(finfo) <- finfo[["ensembl_id_full"]]
    } else {
      rownames(finfo) <- finfo[["ensembl_id"]]
    }
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

  counts <- local({
    h5.fn <- archs4_file_path(paste(org, feature_type, sep = "_"))
    index <- list(NULL, si[[h5col]])
    cnts <- rhdf5::h5read(h5.fn, "data/expression", index=index)
    colnames(cnts) <- rownames(si)
    cnts <- cnts[finfo$h5idx,,drop = FALSE]
    rownames(cnts) <- rownames(finfo)
    cnts
  })

  out <- edgeR::DGEList(counts, genes = finfo, samples = si)
  out
}

