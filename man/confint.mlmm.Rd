% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint.R
\name{confint.mlmm}
\alias{confint.mlmm}
\title{Confidence Intervals for Multiple Linear Mixed Model.}
\usage{
\method{confint}{mlmm}(
  object,
  parm = NULL,
  level = 0.95,
  method = NULL,
  ordering = "parameter",
  ...
)
}
\arguments{
\item{object}{an \code{mlmm} object, output of \code{mlmm}.}

\item{parm}{Not used. For compatibility with the generic method.}

\item{level}{[numeric,0-1] the confidence level of the confidence intervals.}

\item{method}{[character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}, or \code{"pool"}.}

\item{ordering}{[character] should the output be ordered by type of parameter (\code{parameter}) or by model (\code{by}).
Only relevant for \code{mlmm} objects.}

\item{...}{other arguments are passed to \code{\link{confint.Wald_lmm}}.}
}
\description{
Compute confidence intervals for several linear mixed models.
}
\details{
Statistical inference following pooling is performed according to Rubin's rule whose validity requires the congeniality condition of Meng (1994).
}
\references{
Meng X. L.(1994). Multiple-imputation inferences with uncongenial sources of input. Statist. Sci.9, 538–58.
}
\keyword{methods}
