#' Generates markdown enumeration of linked files to download for datadir
#'
#' @noRd
#' @importFrom yaml
#'
#' @param datadir the archs4 data dir
#' @return a character string
md_archs4_download_bullet_list <- function(datadir = getOption("archs4.datadir")) {
  mdat <- archs4_meta(datadir)
  files <- mdat$files
  sources <- sapply(files, "[[", "source")
  sources <- setdiff(sources, "computed")
  mdown <- lapply(sources, function(s) {
    header <- sprintf("* %s\n", s)
    items <- lapply(files, function(f) {
      if (f$source != s) return(NULL)
      sprintf("    - [%s](%s): %s", f$name, f$url, f$description)
    })
    items <- items[sapply(items, Negate(is.null))]
    items <- paste(items, collapse = "\n")
    paste0(header, items)
  })
  paste(mdown, collapse = "\n")
}
