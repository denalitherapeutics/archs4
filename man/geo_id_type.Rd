% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo-utils.R
\name{geo_id_type}
\alias{geo_id_type}
\title{Classify a vector of sample or series GEO ID's as such}
\usage{
geo_id_type(id)
}
\arguments{
\item{id}{a character vector of \code{GSEnnnnn} or \code{GSMnnnnn} ids}
}
\value{
a tibble of \code{unique(id)} indicating if the id is a \code{"series"}
(GSEnnnnn), \code{"sample"} (GSMnnnnn), or \code{"unknown"}.
}
\description{
GEO series identifiers all start with GSE and sample identifiers all
start with GSM. We use that to identify what types of identifiers are
passed into \code{id}
}
