# These are lower-level functions that support "the health" of the local
# datadir that is used to store the data required to drive a
# LocalArchs4Repository

#' Initialize a local datadir to act as an ARCHS4 data datadir
#'
#' @details
#' A local datadir needs to be created and initialized (wth a `meta.yaml`
#' file), to house ARCHS4 data for use in an Archs4Repository. This function
#' creates that datadir and copies an initial `meta.yaml` file.
#'
#' Please refer to the vignette section "ARCHS4 Data Download" for more
#' details.
#'
#' @export
#'
#' @param datadir the path to the datadir to create (or initialize) as an
#'   ARCHS4 data datadir.
#' @param stop_if_exists by default, this function will `stop` if `datadir`
#'   already exists. Set this to `FALSE` to continue. Setting it to `FALSE` is
#'   convenient to initialize the target `datadir` with a `meta.yaml` file.
#'   If a `meta.yaml` file already exists in `datadir`, then this function
#'   will stop unconditionally. Move the `datadir/meta.yaml` out of the way
#'   if you simply want to refresh it with the default version.
#' @return invisibly returns the path to the `meta.yaml` in the target
#'   `datadir`
archs4_local_data_dir_create <- function(datadir = getOption("archs4.datadir"),
                                         stop_if_exists = TRUE) {
  assert_character(datadir)
  d.exists <- file.exists(datadir)

  meta.in <- system.file("extdata", "meta.yaml", package = "archs4",
                         mustWork = TRUE)
  meta.to <- file.path(datadir, "meta.yaml")

  if (d.exists && !dir.exists(datadir)) {
    stop("Desired output datadir is already file(!): ", datadir)
  }

  if (d.exists) {
    if (stop_if_exists) {
      stop("Output datadir already exisits: ", datadir)
    } else {
      if (file.exists(meta.to)) {
        stop("meta.yaml file already exists in output datadir, ",
             "remove it if you want to replace it with the default meta.yaml")
      }
    }
  } else {
    parent.dir <- assert_directory(dirname(datadir), "w")
    dir.create(datadir)
  }
  file.copy(meta.in, meta.to)
  invisible(meta.to)
}

#' Check "the health" of a local ARCHS4 data datadir
#'
#' This function will notify the suer which files are missing from the
#' ARCHS4 data datadir, and what course of action they can use to
#' fix it.
#'
#' @export
#' @param echo echo validation diagnostic message to stdout via [base::cat()]
#' @param datadir the path to the datadir that stores local ARCHS4 data.
#'   Defaults to `getOption("archs4.datadir")`.
#' @return A string that indicates "what's wrong", or `TRUE` if validation
#'   succeeds.
archs4_local_data_dir_validate <- function(echo = TRUE,
                                           datadir = getOption("archs4.datadir")
                                           ) {
  msg <- character()
  if (!dir.exists(datadir)) {
    msg <- paste0(
      "datadir does not exists, run ",
      '`archs4_local_data_dir_create("', datadir, '")`\n')
    if (echo) cat(msg)
    return(invisible(msg))
  }

  meta.fn <- file.path(datadir, "meta.yaml")
  if (!file.exists(meta.fn)) {
    msg <- paste(
      "meta.yaml file is missing from the data datadir, run ",
      "`archs4_local_data_dir_create(datadir, stop_if_exists = FALSE)`\n")
    if (echo) cat(msg)
    return(invisible(msg))
  }

  finfo <- archs4_file_info(datadir)
  missing <- filter(finfo, source == "archs4" & !file_exists)
  if (nrow(missing)) {
    msg <- paste0(
      "The following ARCHS4 files are missing, please download them:\n",
      paste0(
        sprintf("  * %s: %s", missing[["name"]], missing[["url"]]),
        collapse = "\n"))
    if (echo) cat(msg)
    return(invisible(msg))
  }

  missing <- filter(finfo, source == "ensembl" & !file_exists)
  if (nrow(missing)) {
    msg <- paste0(
      "The following ensembl files are missing, please download them:\n",
      paste0(
        sprintf("  * %s: %s", missing[["name"]], missing[["url"]]),
        collapse = "\n"))
    if (echo) cat(msg)
    return(invisible(msg))
  }

  missing <- filter(finfo, source == "computed"  & !file_exists)
  if (nrow(missing)) {
    header <- "The following computed files are missing:"
    filez <- paste(sprintf("  * %s\n", missing[["name"]]), collapse = "")
    advice <- paste0(
      "You can create them by running:\n",
      "  `create_augmented_feature_info(\"", datadir, "\")`")
    msg <- sprintf("%s\n%s\n%s\n\n", header, filez, advice)
    if (echo) cat(msg)
    return(invisible(msg))
  }

  TRUE
}

