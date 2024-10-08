% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.frame.R
\name{model.frame.lmm}
\alias{model.frame.lmm}
\title{Extracting the Model Frame from a Linear Mixed Model}
\usage{
\method{model.frame}{lmm}(
  formula,
  newdata = NULL,
  type = NULL,
  add.index = FALSE,
  na.rm = TRUE,
  ...
)
}
\arguments{
\item{formula}{a \code{lmm} object.}

\item{newdata}{[data.frame] dataset relative to which the model frame should be constructed.}

\item{type}{[character] by default returns the processed dataset used to fit the Linear Mixed Model (\code{NULL}).
Can be used to add rows relative to missing repetitions (\code{"add.NA"})
or obtain a dataset with unique sets of covariates (\code{"unique"}) with respect to the mean structure.}

\item{add.index}{[logical] should columns indexing the row number from the original dataset, time variable, cluster variable, strata variable
be added to the output?}

\item{na.rm}{[logical] should rows containing missing values for the variables used in the linear mixed model be removed?
Not relevant when argument type is \code{"unique"}.}

\item{...}{not used. For compatibility with the generic method.}
}
\description{
Contruct a data frame containing all the variables involved in a Linear Mixed Model.
}
\details{
Column \code{"XXindexXX"} refers to the row of the original dataset (i.e. passed to argument \code{data} when calling \code{\link{lmm}}).
When adding rows relative to missing repetitions, since there is no row in the original dataset, a negative sign is used.
}
\examples{
data("armd.wide", package = "nlmeU")
e.lmH <- lmm(visual52 ~ lesion, structure = IND(~treat.f), data = armd.wide)
model.frame(e.lmH)
model.frame(e.lmH, add.index = TRUE)
model.frame(e.lmH, type = "unique")

data("gastricbypassL", package = "LMMstar")
dfL.NNA <- na.omit(gastricbypassL)
e.lmm <- lmm(glucagonAUC ~ time, repetition = ~visit|id, data = dfL.NNA, df = FALSE)
model.frame(e.lmm, type = "unique")
model.frame(e.lmm, type = c("unique","correlation"))
model.frame(e.lmm, type = "add.NA", add.index = TRUE)
}
\keyword{methods}
