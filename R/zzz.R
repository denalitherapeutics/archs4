validate.data.dir <- function(datadir = getOption("archs4.datadir")) {
  assert_directory(datadir, 'r')
  out <- TRUE

  req.files <- archs4.files()
  req.files <- sapply(req.files, function(x) file.path(datadir, x))

  f.exists <- file.exists(req.files)
  if (!all(f.exists)) {
    out <- paste0("The following files are missing from datadir ",
                  "[", datadir, "]:\n  ",
         paste(names(req.files)[!f.exists], collapse = "\n  "))
  }

  out
}

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
