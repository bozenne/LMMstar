% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.tables.R
\name{model.tables.rbindWald_lmm}
\alias{model.tables.rbindWald_lmm}
\title{Statistical Inference From Combined Wald Tests}
\usage{
\method{model.tables}{rbindWald_lmm}(
  x,
  effects = "Wald",
  level = 0.95,
  df = NULL,
  method = NULL,
  columns = NULL,
  ordering = NULL,
  backtransform = NULL,
  transform.names = TRUE,
  simplify = TRUE,
  ...
)
}
\arguments{
\item{x}{a \code{mlmm} object.}

\item{effects}{[character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
the contrast matrix (\code{"contrast"}),
or the name/value/type of the underlying mixed model parameters (\code{"param"})?}

\item{level}{[numeric, 0-1] nominal coverage of the confidence intervals.}

\item{df}{[logical] Should a Student's t-distribution be used to model the distribution of the Wald statistic. Otherwise a normal distribution is used.}

\item{method}{[character] Should pointwise confidence intervals be output (\code{"none"}) or simultaneous confidence intervals (\code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"})?}

\item{columns}{[character vector] Columns to be output.
Can be any of \code{"estimate"}, \code{"se"}, \code{"df"}, \code{"quantile"}, \code{"lower"}, \code{"upper"}, \code{"statistic"}, \code{"null"}, \code{"p.value"}.}

\item{ordering}{[character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).}

\item{backtransform}{[logical] should the estimates, standard errors, and confidence intervals be backtransformed?}

\item{transform.names}{[logical] should the name of the coefficients be updated to reflect the transformation that has been used?
Only relevant when \code{effects="contrast"}.}

\item{simplify}{[logical] with argument \code{effects="Wald"}, omit from the output attributes containing additional information (e.g. approximation error made when adjusting p-values).
with argument \code{effects="contrast"} the output will be converted into a matrix (instead of a list of matrix) whenever possible.}

\item{...}{Not used. For compatibility with the generic method.}
}
\value{
A \code{data.frame} object.
}
\description{
Combine estimates, standard errors, degrees-of-freedom, confidence intervals (CIs) and p-values
relative to linear contrasts of parameters from different linear mixed models.
}
\details{
When \code{effects} is \code{"Wald"}, this function simply calls \code{\link{confint}} with a specific value for the argument \code{column}.
}
\keyword{methods}
