% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structure.R
\name{UN}
\alias{UN}
\title{Unstructured Covariance Matrix}
\usage{
UN(formula)
}
\arguments{
\item{formula}{Time and cluster variables. The left hand side of the formula is ignored.}
}
\description{
Unstructured covariance matrix
}
\details{
A typical formula would be \code{~time} or \code{~time|id}, indicating a time-specific variance parameter and a correlation parameter specific to each pair of times.
}
