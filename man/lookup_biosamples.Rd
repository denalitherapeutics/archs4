% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geo-utils.R
\name{lookup_biosamples}
\alias{lookup_biosamples}
\title{Retrieve sample annotations from NCBI's Biosample database}
\source{
https://www.ncbi.nlm.nih.gov/books/NBK25499/#\emph{chapter4_ESearch}
}
\usage{
lookup_biosamples(x, retmax = 1e+05 - 1L)
}
\arguments{
\item{x}{Character vector of sample identifiers to search the Biosample
database for. Typically either \code{BioSample (SAMN)}, \code{SRA (SRS)} or
\code{GEO (GSM)} accession numbers.}

\item{retmax}{Scalar integer, the maximum number of (total) matches to
retrieve from Entrez. See
\url{https://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_}
for details. The number of records that can be retrieved in one query
must be < 100,000.}
}
\value{
A tbl_df data.frame with sample annotations. Column names and numbers
vary depending on the attributes available in the Biosample database.
}
\description{
This function uses the \code{rentrez} package to retrieve sample annotations
from NCBI's Biosample database.
}
\examples{
if (interactive()) {
 # BioSample identifiers
 lookup_biosamples(c("GSM1947162", "GSM1947179"))
 # mixed SRS and GSM identifiers
 lookup_biosamples(c("SRS1171537", "SRS1171536", "GSM1947179"))
 # mixed samples from two different studies (with different attributes)
 lookup_biosamples(c("SRS1171537", "SRS1271536"))
}
}
