.onLoad <- function(libname, pkgname) {
  ## Setup default option values
  opts <- options()

  pkg.opts <- list(
    archs4.datadir='~/.archs4data')
  toset <- !(names(pkg.opts) %in% names(opts))
  if (any(toset)) {
    options(pkg.opts[toset])
  }

  ddir <- getOption("archs4.datadir")
  kosher <- archs4_local_data_dir_validate(echo = FALSE, datadir = ddir)
  if (!isTRUE(kosher)) {
    message(
      "Note that your default archs4 data directory is NOT setup correctly\n\n",
      "  * Run `archs4_local_data_dir_validate()` to diagnose\n",
      "  * Refer to the ARCHS4 Data Download section of the archs4 vignette ",
      "for more information\n\n",
      "Your default archs4 data directory (`getOption(\"archs4.datadir\")`) ",
      "is:\n\n  ", ddir, "\n\n")
  }

  invisible()
}
