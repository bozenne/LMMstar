% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multcomp.R
\name{estfun}
\alias{estfun}
\alias{estfun.lmm}
\title{Extract the Score Function for Multcomp}
\usage{
\method{estfun}{lmm}(x, ...)
}
\arguments{
\item{x}{a \code{lmm} object.}

\item{...}{Not used. For compatibility with the generic method.}
}
\description{
Extract the Score Function for Multcomp. For internal use.
}
\examples{
## simulate data in the long format
set.seed(10)
dL <- sampleRem(100, n.times = 3, format = "long")

## fit Linear Mixed Model
eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)

## test multiple linear hypotheses
if(require(multcomp)){
multcomp::glht(eUN.lmm)
}
}
