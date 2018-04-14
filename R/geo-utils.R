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

#' Retrieve sample annotations from NCBI's Biosample database
#'
#' This function uses the \code{rentrez} package to retrieve sample annotations
#' from NCBI's Biosample database.
#' @param x Character vector of sample identifiers to search the Biosample
#' database for. Typically either \code{BioSample (SAMN)}, \code{SRA (SRS)} or
#' \code{GEO (GSM)} accession numbers.
#' @param retmax Scalar integer, the maximum number of (total) matches to
#' retrieve from Entrez. See
#' \url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_}
#' for details. The number of records that can be retrieved in one query
#' must be < 100,000.
#' @return A tbl_df data.frame with sample annotations. Column names and numbers
#' vary depending on the attributes available in the Biosample database.
#' @export
#' @importFrom magrittr "%>%" 
#' @importFrom rentrez entrez_search entrez_fetch
#' @importFrom xml2 read_xml xml_find_all xml_text xml_attr xml_find_first
#' xml_children
#' @importFrom purrr map
#' @importFrom tibble set_tidy_names
#' @importFrom tidyr unnest spread
#' @importFrom dplyr rename rename_all
#' @importFrom readr type_convert cols
#' @source https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_
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
    purrr::map(.f = function(x) {
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
    dplyr::bind_rows() %>%
    tidyr::spread(key = key, value = value) %>%
    tibble::set_tidy_names(syntactic = TRUE, quiet = TRUE) %>%
    dplyr::rename(cell_type = cell.type) %>%
    dplyr::rename_all(tolower) %>%
    readr::type_convert(col_types = readr::cols())
}