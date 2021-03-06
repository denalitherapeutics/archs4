% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Archs4Repository.R, R/archs4-functional.R
\name{sample_info}
\alias{sample_info}
\alias{archs4_sample_info}
\title{Retrieves information for samples by GSE series or sample IDs}
\usage{
sample_info(x, id, columns = c("Sample_title", "Sample_source_name_ch1"),
  check_missing_samples = TRUE, ...)

archs4_sample_info(id, columns = c("Sample_title",
  "Sample_source_name_ch1"), sample_table = archs4_sample_table(datadir =
  datadir), sample_covariates = archs4_sample_covariates(datadir),
  check_missing_samples = TRUE, datadir = getOption("archs4.datadir"),
  ...)
}
\arguments{
\item{x}{an \code{Archs4Repository}}

\item{id}{a character vector of GEO series or sample ids.}

\item{columns}{the names of the sample metadata columns desired. This
defaults to \code{c("Sample_title", "Sample_source_name_ch1")}. The values
in \code{columns} must be a subset of the values enumerated in
\code{\link[=archs4_sample_covariates]{archs4_sample_covariates()}}.}

\item{check_missing_samples}{When \code{TRUE} (the default), this function will
check every unique GEO series identifier (\code{"GSEnnnn"}) for missing samples
by using an NCBI Rest service via a call to \code{\link[=archs4_series_status]{archs4_series_status()}},
and \code{\link[=lookup_gse]{lookup_gse()}}.}

\item{sample_table}{the output from \code{\link[=archs4_sample_table]{archs4_sample_table()}}, which lists
the series_id,sample_id combinations found in the ARCHS4 repository.}

\item{sample_covariates}{the \code{data.frame}-definition of the sample covariates
found in the ARCHS4 datasetes, which is constructed via a call to
\code{\link[=archs4_sample_covariates]{archs4_sample_covariates()}}. The parameter is included in here so that
a cached version of this \code{data.frame} can be re-used.}

\item{datadir}{the directory that holds the archs4 data}
}
\value{
a tibble of series_id, sample_id, sample_h5idx, sample_title, and
sample_name columns. If the query sample or series query can't be found,
then there will be an \code{NA} value for these columns. The \code{query_type} column
will indicat whether the row was returned from querying by series or by
sample.
}
\description{
Fetch a tibble of series and sample information by querying the arcsh4
dataset by GEO sample (GSE) or sample (GSM) ids.

For each unique GEO series identifier (\code{"GSEnnnn"}), we will check if the
ARCHS4 dataset is missing any of its samples when \code{check_missing_samples}
is set to \code{TRUE} (default).
}
\examples{
si <- archs4_sample_info("GSE52564") # ben barres transcriptome db ...
}
