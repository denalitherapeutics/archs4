.onLoad <- function(libname, pkgname) {
  ## Setup default option values
  opts <- options()

  pkg.opts <- list(
    archs4.datadir='~/.archs4data')
  toset <- !(names(pkg.opts) %in% names(opts))
  if (any(toset)) {
    options(pkg.opts[toset])
  }

  kosher <- validate.data.dir(getOption("archs4.datadir"))
  if (!isTRUE(kosher)) {
    message(kosher)
  }

  invisible()
}
