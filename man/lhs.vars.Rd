% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{lhs.vars}
\alias{lhs.vars}
\title{Extract Variables in the Left Hand Side of a Formula}
\usage{
lhs.vars(formula)
}
\arguments{
\item{formula}{a formula}
}
\description{
Extract Variables in the Left Hand Side of a Formula
}
\examples{
lhs.vars(~X)
lhs.vars(Y~X)
lhs.vars(Y~X1+X2)
lhs.vars(Y~1)
lhs.vars(Y+Z~1)
lhs.vars(~1|id)
lhs.vars(Y~X|id)
lhs.vars(Y+Z~X|id)
}
