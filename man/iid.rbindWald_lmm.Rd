% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iid.R
\name{iid.rbindWald_lmm}
\alias{iid.rbindWald_lmm}
\title{Extract the Influence Function From Wald Tests}
\usage{
\method{iid}{rbindWald_lmm}(x, effects = "Wald", ordering = NULL, transform.names = TRUE, ...)
}
\arguments{
\item{effects}{[character] should the influence function for the linear contrasts involved in the Wald test be output (\code{"Wald"}),
or for the linear mixed model parameters (\code{"all"})?}

\item{ordering}{[character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).}

\item{transform.names}{[logical] should the name of the coefficients be updated to reflect the transformation that has been used?
Only relevant when \code{effects="all"}.}

\item{...}{Not used. For compatibility with the generic method.}

\item{object}{a \code{rbindWald_lmm} object.}
}
\value{
A matrix with one row per observation and one column per parameter (\code{effects="Wald"} or \code{effects="all"}) or a logical value (\code{effects="test"}).
}
\description{
Extract the influence function of linear contrasts applied to an ML or REML estimator of parameters from a linear mixed model.
}
