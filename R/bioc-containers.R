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
as.DGEList <- function(id, feature_type = c("gene", "transcript"),
                       row_id = c("ensembl", "symbol")) {
  feature_type <- match.arg(feature_type)
  if (feature_type == "transcript") {
    stop("transcript level features not supported yet")
  }

  row_id <- match.arg(row_id)

  si <- archs4_sample_info(id)
  org <- unique(si$organism)
  if (length(org) != 1L) {
    stop("You are querying across species")
  }

  bad.id <- filter(si, is.na(sample_h5idx))
  if (nrow(bad.id)) {
    stop("The following samples could not be found: ",
         paste(sprintf("%s::%s", bad.id$series_id, bad.id$sample_id),
               collapse = "; "))
  }
  si <- as.data.frame(si, strinsAsFactors = FALSE)
  rownames(si) <- si$sample_id

  counts <- local({
    h5.fn <- archs4_file_path(paste(org, feature_type, sep = "_"))
    cnts <- rhdf5::h5read(h5.fn, "data/expression",
                          index=list(NULL, si$sample_h5idx))
    colnames(cnts) <- rownames(si)
    cnts
  })

  finfo <- archs4_feature_info(feature_type, org)
  finfo <- as.data.frame(finfo, stringsAsFactors = FALSE)
  if (feature_type == "gene") {
    if (row_id == "ensembl") {
      finfo <- filter(finfo, !is.na(ensembl_gene_id))
      counts <- counts[finfo$h5idx,,drop = FALSE]
      rownames(counts) <- finfo[["ensembl_gene_id"]]
      rownames(finfo) <- finfo[["ensembl_gene_id"]]
    } else {
      rownames(counts) <- finfo[["symbol"]]
      rownames(finfo) <- finfo[["symbol"]]
    }
  } else {
    stop("TODOO: what to do with transcript-level finfo")
  }

  out <- edgeR::DGEList(counts, genes = finfo, samples = si)
  out
}

