% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R
\name{coef.Wald_lmm}
\alias{coef.Wald_lmm}
\title{Extract Coefficients From Wald Tests for Linear Mixed Model}
\usage{
\method{coef}{Wald_lmm}(object, effects = "Wald", backtransform = NULL, simplify = TRUE, ...)
}
\arguments{
\item{object}{a \code{Wald_lmm} object.}

\item{effects}{[character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
or the value of the linear mixed model parameters (\code{"all"})?}

\item{backtransform}{[logical] should the estimates be back-transformed?}

\item{simplify}{[logical] omit from the output the attribute containing the type of each parameter or contrast (mu/sigma/k/rho).}

\item{...}{Not used. For compatibility with the generic method.}
}
\description{
Extract coefficients from Wald tests applied to a linear mixed model.
}
