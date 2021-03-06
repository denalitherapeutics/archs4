% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression.R
\name{estimate_norm_factors}
\alias{estimate_norm_factors}
\title{Calculates library size and norm factors for a specific dataset}
\usage{
estimate_norm_factors(a4, key, n = 500, logratioTrim = 0.3,
  sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10, p = 0.75,
  ...)
}
\description{
This function will serialize the results of the library size and
normalization factors into files inside \code{datadir(a4)}
}
