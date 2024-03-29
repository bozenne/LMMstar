% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.R
\name{add}
\alias{add}
\title{Add Columns to Output}
\usage{
add(...)
}
\arguments{
\item{...}{[character vector] name of the columns to be added to the default output.}
}
\value{
A character vector
}
\description{
Auxiliary function that can be used when specifying the argument \code{columns} (e.g. calling \code{confint.lmm}) to add columns.
}
\examples{
set.seed(10)
dW <- sampleRem(25, n.times = 1, format = "long")
e.lmm <- lmm(Y~X1, data = dW)

confint(e.lmm, columns = add("statistic"))
}
\keyword{utilities}
