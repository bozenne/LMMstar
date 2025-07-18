% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate.R
\name{estimate.Wald_lmm}
\alias{estimate.Wald_lmm}
\title{Delta Method for Wald Tests}
\usage{
\method{estimate}{Wald_lmm}(x, f, df = x$args$df, level = 0.95, method.numDeriv = NULL, ...)
}
\arguments{
\item{x}{a \code{Wald_lmm} object.}

\item{f}{[function] function taking as input \code{object}, the \code{Wald_lmm} object, which outputs the parameter(s) of interest.
Can accept extra-arguments.}

\item{df}{[logical] Should degree-of-freedom, computed using Satterthwaite approximation, for the parameter of interest be output.
Can also be a numeric vector providing providing the degrees-of-freedom relative to each estimate.}

\item{level}{[numeric,0-1] the confidence level of the confidence intervals.}

\item{method.numDeriv}{[character] method used to approximate the gradient: either \code{"simple"} or \code{"Richardson"}.
Passed to \code{numDeriv::jacobian}.}

\item{...}{extra arguments passed to \code{f}.}
}
\description{
Estimate standard errors, confidence intervals, and p-values for a smooth transformation of parameters involved in Wald tests.
}
\details{
Based a first order delta method to evaluate the variance of the estimate.
The derivative of the transformation is evaluated using numerical differentiation (\code{numDeriv::jacobian}). \cr \cr

Compared to \code{estimate.lmm}, the \code{f} argument cannot take \code{p} as argument.
One should instead explicitely request the estimated contrasts in the function \code{f} using \code{coef(object)}.
}
\keyword{mhtest}
