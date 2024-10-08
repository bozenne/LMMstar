% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logLik.R
\name{logLik.lmm}
\alias{logLik.lmm}
\title{Log-Likelihood From a Linear Mixed Model}
\usage{
\method{logLik}{lmm}(object, newdata = NULL, p = NULL, indiv = FALSE, ...)
}
\arguments{
\item{object}{a \code{lmm} object.}

\item{newdata}{[data.frame] dataset relative to which the log-likelihood should be computed. Only relevant if differs from the dataset used to fit the model.}

\item{p}{[numeric vector] value of the model coefficients at which to evaluate the log-likelihood. Only relevant if differs from the fitted values.}

\item{indiv}{[logical] Should the contribution of each cluster to the log-likelihood be output? Otherwise output the sum of all clusters of the derivatives.}

\item{...}{Not used. For compatibility with the generic method.}
}
\value{
A numeric value (total logLikelihood) or a vector of numeric values, one for each cluster (cluster specific logLikelihood).
}
\description{
Extract or compute the log-likelihood of a linear mixed model.
}
\details{
\bold{indiv}: only relevant when using maximum likelihood. Must be \code{FALSE} when using restricted maximum likelihood.
}
\keyword{methods}
