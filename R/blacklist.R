#' Return a list of series/samples you may want to ignore for now
#'
#' The entries that appear here so far are because the data appear to come from
#' single-cell experiments
#'
#' @export
blacklist <- function() {
  fn <- system.file("extdata", "blacklist.csv", package = "archs4")
  read.csv(fn, stringsAsFactors = FALSE)
}
