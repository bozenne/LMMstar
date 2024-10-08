% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R
\name{coef.rbindWald_lmm}
\alias{coef.rbindWald_lmm}
\title{Extract Coefficients From Combined Wald Tests Applied to Linear Mixed Models}
\usage{
\method{coef}{rbindWald_lmm}(
  object,
  effects = "Wald",
  method = "none",
  ordering = NULL,
  transform.names = TRUE,
  backtransform = NULL,
  simplify = TRUE,
  ...
)
}
\arguments{
\item{object}{a \code{rbindWald_lmm} object.}

\item{effects}{[character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
or the value of the linear mixed model parameters (\code{"all"})?}

\item{method}{[character vector] should the estimated value for the linear contrasts be output (one of \code{"none"}, \code{"bonferroni"}, ..., \code{"single-step2"})
and/or pooled linear contrast estimate(s) (\code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}, \code{"p.rejection"})?
Only relevant when \code{effects = "Wald"}.}

\item{ordering}{[character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).}

\item{transform.names}{[logical] should the name of the coefficients be updated to reflect the transformation that has been used?
Only relevant when \code{effects="all"}.}

\item{backtransform}{[logical] should the estimates be back-transformed?}

\item{simplify}{[logical] should the output be a vector or a list with one element specific to each possible ordering (i.e. contrast or model).
Only relevant when argument \code{method} refers to multiple comparisons and not to a pooling method.}

\item{...}{Not used. For compatibility with the generic method.}
}
\description{
Combine estimated values across linear contrasts applied on parameters from different linear mixed models.
}
\details{
Argument \bold{effects}: when evaluating the proportion of rejected hypotheses (\code{effects="p.rejection"})
a \code{"single-step"} method will be used by default to evaluate the critical quantile.
This can be changed by adding adjustment method, e.g. \code{effects=c("bonferronin","p.rejection"}, in the argument.
}
