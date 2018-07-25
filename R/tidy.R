#' @noRd
#' @export
#' @importFrom edgeR cpm
#' @importFrom reshape2 melt
#' @importFrom broom tidy
#'
#' @method tidy DGEList
tidy.DGEList <- function(x, normalized.lib.sizes = TRUE, prior.count = 3, ...) {
  mats <- list(
    cpm = cpm(x, normalized.lib.sizes = normalized.lib.sizes,
              log = TRUE, prior.count = prior.count),
    count = x$counts)

  .tidy.core(mats, genes = x$genes, samples = x$samples)
}

#' @noRd
#' @export
#' @importFrom SummarizedExperiment assay assayNames rowData colData
#' @importFrom edgeR DGEList
#' @method tidy SummarizedExperiment
tidy.SummarizedExperiment <- function(x, assay_name = NULL, is_counts = TRUE,
                                      ...) {
  if (is.null(assay_name)) {
    assay_name <- assayNames(x)[1L]
  }
  assert_string(assay_name)
  assert_choice(assay_name, assayNames(x))

  dat <- assay(x, assay_name)
  ginfo <- as.data.frame(rowData(x))
  sinfo <- as.data.frame(colData(x))

  if (is_counts) {
    y <- calcNormFactors(DGEList(dat, samples = sinfo, genes = ginfo))
    out <- tidy(y, ...)
  } else {
    out <- .tidy.core(dat, genes = ginfo, samples = sinfo, ...)
  }

  out
}

#' @noRd
#' @export
#' @method tidy EList
tidy.EList <- function(x, ...)  {
  mats <- list(cpm = x$E)
  if (is.matrix(x$weights)) {
    mats$weight <- x$weights
    rownames(mats$weight) <- rownames(x)
    colnames(mats$weight) <- colnames(x)
  } else {
    names(mats)[1L] <- "value"
  }

  .tidy.core(mats, genes = x$genes, samples = x$targets)
}

.tidy.core <- function(mats, genes, samples, ...) {
  if (is.matrix(mats)) mats <- list(value = mats)
  stopifnot(is.list(mats))
  stopifnot(all(sapply(mats, is.matrix)))
  assert_named(mats, type = "unique")

  rnames <- rownames(mats[[1]])
  snames <- colnames(mats[[1]])
  genes$.gene_id <- rnames
  gid.col <- sapply(genes, function(xx) all(xx == rnames))
  gid.col <- colnames(genes)[which(gid.col)[1L]]
  if (gid.col != ".gene_id") genes$.gene_id <- NULL

  samples$.sample_id <- snames
  sid.col <- sapply(samples, function(xx) all(xx == snames))
  sid.col <- colnames(samples)[which(sid.col)[1L]]
  if (sid.col != ".sample_id") samples$.sample_id <- NULL

  adat.all <- lapply(names(mats), function(mname) {
    m <- mats[[mname]]
    stopifnot(all.equal(rownames(m), rnames))
    m <- melt(m)
    m <- transform(m, Var1 = as.character(Var1), Var2 = as.character(Var2))
    colnames(m) <- c(gid.col, sid.col, mname)
    m
  })
  adat <- do.call(cbind, adat.all)
  # if there were multiple matrices, there will be multiple sample_id columns
  # so we remove those
  adat <- adat[, !duplicated(colnames(adat))]
  out <- inner_join(adat, genes, by = gid.col)
  out <- inner_join(out, samples, by = sid.col)
  out
}
