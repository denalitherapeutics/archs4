% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression.R
\name{fetch_expression}
\alias{fetch_expression}
\title{Retrieve expression data for genes/transcripts across samples}
\usage{
fetch_expression(a4, features, samples = NULL, feature_type = "gene",
  source = NULL, feature_meta = "symbol",
  sample_meta = c("Sample_title", "Sample_source_name_ch1"),
  prior.count = 3, ...)
}
\arguments{
\item{a4}{An \code{Archs4Repository}}

\item{features}{a \code{tibble} of feature descriptors (as returned from
\code{\link[=feature_info]{feature_info()}}). Really this table simply needs the following columns:
\itemize{
\item \code{"a4name"} or \code{"ensembl_id" or};
\item \code{"feature_type"}: (gene or transcript)
\item \code{"source"} (human or mouse),;
}}

\item{samples}{a samples identifier: ie. a tibble with series_id and
sample_id columns.}

\item{source}{"mouse" or "human" pass this in explicitly if there is no
"source" column in your \code{features} or \code{samples} data.frame}

\item{feature_meta}{additional columns from
\code{feature_info(a4, feature_type = feature_type)}}

\item{sample_meta}{metadata information to include for the samples. These
values are extracted from the \code{meta/VALUE} files in the respective hdf5
data files. See \code{\link[=archs4_sample_info]{archs4_sample_info()}} for more details.}

\item{type}{"gene" or "transcript"-level data?}
}
\value{
a data.frame of expression data. Each row is an observation
(one gene one sample)
}
\description{
Returns the estimated counts for the genes/transcripts enumerated in the
\code{features} table for the samples enumerated in the \code{samples} table.
}
\details{
Note that the values returned are simply estimated counts. They are not
normalized by sequencing depth. For now, the only use for this function is
to compare how ratios of genes compare across samples.
}
\examples{
a4 <- Archs4Repository()
gnz <- feature_lookup(a4, c("CFAP65", "PECR"))
gexprs <- fetch_expression(a4, gnz)
}
\seealso{
\code{\link[=as.DGEList]{as.DGEList()}}
}
