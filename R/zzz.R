.onLoad <- function(libname, pkgname) {
  ## Setup default option values
  opts <- options()

  pkg.opts <- list(
    archs4.datadir='~/.archs4data')
  toset <- !(names(pkg.opts) %in% names(opts))
  if (any(toset)) {
    options(pkg.opts[toset])
  }

  kosher <- archs4_local_data_dir_validate(getOption("archs4.datadir"))
  if (!isTRUE(kosher)) {
    message(kosher)
  }

  invisible()
}
