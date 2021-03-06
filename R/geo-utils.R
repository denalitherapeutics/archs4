#' Classify a vector of sample or series GEO ID's as such
#'
#' GEO series identifiers all start with GSE and sample identifiers all
#' start with GSM. We use that to identify what types of identifiers are
#' passed into `id`
#'
#' @export
#'
#' @param id a character vector of `GSEnnnnn` or `GSMnnnnn` ids
#' @return a tibble of `unique(id)` indicating if the id is a `"series"`
#'   (GSEnnnnn), `"sample"` (GSMnnnnn), or `"unknown"`.
geo_id_type <- function(id) {
  id <- assert_character(id) %>% unique
  type <- case_when(
    is_geo_series_id(id) ~ "series",
    is_geo_sample_id(id) ~ "sample",
    TRUE                 ~ "unknown")
  tibble(id = id, type = type)
}


#' Query NCBI GEO through its REST interface
#'
#' @export
#' @importFrom xml2 read_xml xml_validate xml_ns_strip xml_contents xml_find_all
#'   xml_text
#' @source https://www.ncbi.nlm.nih.gov/geo/info/download.html
#'
#' @param accession Scalar character, GEO identifier for a series (GSE), a
#'   sample (GSM) or a platform (GPL).
#' @param validate Scalar boolean, validate the retrieved xml file against
#'   NCBI's schema?
#' @return xml2::xml_document object
#' @examples
#' query_geo("GSE109171")
query_geo <- function(accession, target = c("self", "gsm", "gpl", "gse", "all"),
                      validate = FALSE, verbose = FALSE) {
  target = match.arg(target)
  geo_url <- sprintf(
    paste0("https://www.ncbi.nlm.nih.gov/geo/query/",
           "acc.cgi?acc=%s&targ=%s&view=%s&form=%s"),
    accession, target, "full", "xml")
  if (verbose) {
    message(sprintf("Retrieving %s", geo_url))
  }
  res <- xml2::read_xml(geo_url)

  if (validate) {
    schema_url <- strsplit(xml2::xml_attr(res, attr = "schemaLocation"),
                           split = " ", fixed = TRUE)[[1]][-1]
    valid_xml <- xml2::xml_validate(res, schema = xml2::read_xml(schema_url))
    if (!valid_xml) {
      stop(sprintf("Validation with schema %s failed", schema_url))
    }
  }
  return(res)
}

#' Retrieve information about a GEO series
#'
#' Queries NCBI GEO's REST interface to retrieve e.g. title, summary and the
#' list of samples for a GEO series.
#'
#' @export
#' @importFrom xml2 xml_contents xml_find_all xml_text
#' @importFrom stats setNames
#'
#' @param accession Scalar character, GEO series identifier e.g. GSE109171
#' @param fields Character vector specifying which fields to extract from the
#'   XML file returned by GEO
#' @param ... Additional arguments passed on to the `query_geo` function.
#' @return List the requested `fields`
#' @examples
#' if (interactive()) {
#'   lookup_gse("GSE109171")
#' }
lookup_gse <- function(accession,
                       fields = c("Accession", "Title", "Summary",
                                  "Overall-Design", "Type", "Pubmed-ID",
                                  "Sample"),
                       ...) {
  fields <- match.arg(fields, several.ok = TRUE)
  xml <- query_geo(accession, target = "gse", ...) %>%
    xml2::xml_ns_strip()
  series_fields <- setdiff(fields, "Sample")
  series <- purrr::map(
    setNames(series_fields, tolower(series_fields)),
    .f = function(field) {
      xml %>%
        xml2::xml_find_first(xpath = "Series") %>%
        xml2::xml_find_first(field) %>%
        xml2::xml_text(trim = TRUE)
    })
  sample_fields <- setdiff(fields, series_fields)
  samples <- purrr::map(
    setNames(sample_fields, tolower(sample_fields)),
    .f = function(field) {
      xml %>%
        xml2::xml_find_all(xpath = "Sample") %>%
        xml2::xml_text(trim = TRUE)
    })
  append(series, samples)
}

#' Retrieve sample annotations from NCBI's Biosample database
#'
#' This function uses the `rentrez` package to retrieve sample annotations
#' from NCBI's Biosample database.
#'
#' @export
#' @importFrom rentrez entrez_search entrez_fetch
#' @importFrom xml2 read_xml xml_find_all xml_text xml_attr xml_find_first
#'   xml_children
#' @importFrom tibble set_tidy_names
#' @importFrom readr type_convert cols
#' @source https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_
#'
#' @param x Character vector of sample identifiers to search the Biosample
#'   database for. Typically either `BioSample (SAMN)`, `SRA (SRS)` or
#'   `GEO (GSM)` accession numbers.
#' @param retmax Scalar integer, the maximum number of (total) matches to
#' retrieve from Entrez. See
#'   \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_}
#'   for details. The number of records that can be retrieved in one query
#'   must be < 100,000.
#' @return A tbl_df data.frame with sample annotations. Column names and numbers
#'   vary depending on the attributes available in the Biosample database.
#' @examples
#' if (interactive()) {
#'  # BioSample identifiers
#'  lookup_biosamples(c("GSM1947162", "GSM1947179"))
#'  # mixed SRS and GSM identifiers
#'  lookup_biosamples(c("SRS1171537", "SRS1171536", "GSM1947179"))
#'  # mixed samples from two different studies (with different attributes)
#'  lookup_biosamples(c("SRS1171537", "SRS1271536"))
#' }
lookup_biosamples <- function(x, retmax = 1e5 - 1L) {
  if (retmax >= 1e5) stop("retmax must be < 100,000")

  search_results <- rentrez::entrez_search(
    db = "biosample", retmax = retmax,
    term = sprintf("%s[ACCN]", paste(x, collapse = "[ACCN] OR ")))
  rentrez::entrez_fetch(db = "biosample", id = search_results$ids,
                        rettype = "xml", retmax = retmax,
                        parsed = FALSE) %>%
    xml2::read_xml() %>%
    xml2::xml_children() %>%
    purrr::map_df(.f = function(x) {
      tibble::data_frame(
        Title = xml2::xml_find_first(
          x,
          xpath = ".//Description/Title") %>%
          xml2::xml_text(trim = TRUE),
        OrganismName = xml2::xml_find_first(
          x,
          xpath = ".//Description/Organism/OrganismName") %>%
          xml2::xml_text(trim = TRUE),
        taxonomy_id = xml2::xml_find_all(
          x,
          xpath = ".//Description/Organism") %>%
          xml2::xml_attr("taxonomy_id"),
        BioSample = xml2::xml_find_first(
          x,
          xpath = './/Ids/Id[@db="BioSample"]') %>%
          xml2::xml_text(trim = TRUE),
        SRA = xml2::xml_find_first(
          x,
          xpath = './/Ids/Id[@db="SRA"]') %>%
          xml2::xml_text(trim = TRUE),
        GEO = xml2::xml_find_first(
          x,
          xpath = './/Ids/Id[@db="GEO"]') %>%
          xml2::xml_text(trim = TRUE),
        value = xml2::xml_find_all(
          x,
          xpath = ".//Attribute") %>%
          xml2::xml_text(trim = TRUE),
        key = xml2::xml_find_all(
          x,
          xpath = ".//Attribute") %>%
          xml2::xml_attr("attribute_name")
      )
    }) %>%
    dplyr::distinct() %>%
    tidyr::spread(key = key, value = value) %>%
    tibble::set_tidy_names(syntactic = TRUE, quiet = TRUE) %>%
    dplyr::rename_all(tolower) %>%
    readr::type_convert(col_types = readr::cols())
}

#' Retrieve metadata for an SRA accession
#'
#' This function uses the EBI's or the NCBI's REST APIs to retrieve information
#' about SRA data.
#' Study accessions (ERP, SRP, DRP, PRJ prefixes), experiment accessions
#' (ERX, SRX, DRX prefixes), sample accessions (ERS, SRS, DRS, SAM prefixes)
#' and run accessions (ERR, SRR, DRR prefixes) can be supplied.
#' For more information see \url{http://www.ebi.ac.uk/ena/browse/file-reports}
#' @param x SRA identifier
#' @param from Scalar character, specifying either \code{ncbi} or \code{ena} as
#' the source database
#' @return A tbl_df data.frame
#' @note The output data.frame will be different for the two source databases.
#' @importFrom readr read_tsv read_csv cols
#' @export
#' @examples
#' if (interactive()) {
#'    # retrieve study annotations
#'    retrieve_sra_metadata("SRP066489")
#'    # paired-end samples
#'    retrieve_sra_metadata("PRJEB2054", "ena") %>%
#'    dplyr::filter(sample_accession == "SAMEA728920")
#' }
retrieve_sra_metadata <- function(x, from = c("ena", "ncbi")) {
  from <- match.arg(from)
  runinfo <- switch(
    from,
    ncbi = {
      sra_url <- sprintf(
        paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?",
               "save=efetch&rettype=runinfo&db=sra&term=%s"), x)
      readr::read_csv(url(sra_url), col_types = readr::cols())
    },
    ena = {
      sra_url <- sprintf(paste0(
        "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s",
        "&result=read_run"), x)
      readr::read_tsv(url(sra_url), col_types = readr::cols())
    }
  )
  return(runinfo)
}
