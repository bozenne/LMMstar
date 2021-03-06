% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.R
\name{residuals}
\alias{residuals}
\alias{residuals.lmm}
\title{Extract The Residuals From a Linear Mixed Model}
\usage{
\method{residuals}{lmm}(
  object,
  type.residual = "response",
  format = "long",
  data = NULL,
  p = NULL,
  keep.data = FALSE,
  plot = "none",
  engine.qqplot = "ggplot2",
  digit.cor = 2,
  size.text = 16,
  type.object = "lmm",
  ...
)
}
\arguments{
\item{object}{a \code{lmm} object.}

\item{type.residual}{[character] Should the raw residuals be output (\code{"response"}), the Pearson residuals (\code{"pearson"}), the studentized residuals (\code{"studentized"}), or normalized residuals (\code{"normalized"}).}

\item{format}{[character] Should the residuals be output relative as a vector (\code{"long"}), or as a matrix with in row the clusters and in columns the outcomes (\code{"wide"}).}

\item{data}{[data.frame] dataset relative to which the residuals should be computed. Only relevant if differs from the dataset used to fit the model.}

\item{p}{[numeric vector] value of the model coefficients at which to evaluate the residuals. Only relevant if differs from the fitted values.}

\item{keep.data}{[logical] Should the argument \code{data} be output along side the residuals? Only possible in the long format.}

\item{plot}{[character] Should a qqplot (\code{"qqplot"}), or a heatmap of the correlation between residuals  (\code{"correlation"}, require wide format), or a plot of residuals along the fitted values (\code{"scatterplot"}, require long format) be displayed?}

\item{engine.qqplot}{[character] Should ggplot2 or qqtest be used to display quantile-quantile plots? Only used when argument \code{plot} is \code{"qqplot"}.}

\item{digit.cor}{[integer, >0] Number of digit used to display the correlation coefficients? No correlation coefficient is displayed when set to 0. Only used when argument \code{plot} is \code{"correlation"}.}

\item{size.text}{[numeric, >0] Size of the font used to displayed text when using ggplot2.}

\item{type.object}{[character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.}

\item{...}{Not used. For compatibility with the generic method.}
}
\value{
When argument format is \code{"long"} and type.oobject is \code{"lmm"}, a vector containing the value of the residual realtive to each observation.
It is a matrix if the argument \code{type.residual} contains several values.
When argument format is \code{"wide"} and type.oobject is \code{"lmm"}, a data.frame with the value of the residual relative to each cluster (in rows) at each timepoint (in columns).
}
\description{
Extract or compute the residuals of a multivariate gaussian model.
}
\details{
The argument \code{type.residual} defines how the residuals are computed:
\itemize{
\item \code{"raw"}: observed outcome minus fitted value \eqn{\varepsilon = Y_{ij} - X_{ij} \hat{\beta}}.
\item \code{"pearson"}: each raw residual is divided by its modeled standard deviation \eqn{\varepsilon = \frac{Y_{ij} - X_{ij} \hat{\beta}}{\sqrt{\hat{\omega}_{ij}}}}.
\item \code{"studentized"}: same as \code{"pearson"} but excluding the contribution of the cluster in the modeled standard deviation  \eqn{\varepsilon = \frac{Y_{ij} - X_{ij} \hat{\beta}}{\sqrt{\hat{\omega}_{ij}-\hat{q}_{ij}}}}.
\item \code{"normalized"}: raw residuals are multiplied, within clusters, by the inverse of the (lower) Cholesky factor of the modeled residual variance covariance matrix \eqn{\varepsilon = ( Y_{i} - X_{i} \hat{\beta} )\hat{C}^{-1}}.
\item \code{"normalized2"}: same as \code{"normalized"} but excluding the contribution of the cluster in the modeled residual variance covariance matrix \eqn{\varepsilon = ( Y_{i} - X_{i} \hat{\beta} ) \hat{D}_i^{-1}}.
}
where
\itemize{
\item \eqn{X} the design matrix
\item \eqn{Y} the outcome
\item \eqn{\hat{\beta}} the estimated mean coefficients
\item \eqn{\hat{\Omega}} the modeled variance-covariance of the residuals and \eqn{\hat{\omega}} its diagonal elements
\item \eqn{\hat{C}} the lower Cholesky factor of \eqn{\hat{\Omega}}, i.e. \eqn{\hat{C} \hat{C}^{t} = \hat{\Omega}}
\item \eqn{\hat{Q}_i= X_i (X^{t}\hat{\Omega}X)^{-1}X_i^{t}} a cluster specific correction factor, approximating the contribution of cluster i to \eqn{\hat{Omega}}. Its diagonal elements are denoted \eqn{\hat{q}_i}.
\item \eqn{\hat{D}_i} the lower Cholesky factor of \eqn{\hat{\Omega}-\hat{Q}_i}
}
}
\examples{
## simulate data in the long format
set.seed(10)
dL <- sampleRem(100, n.times = 3, format = "long")

## fit Linear Mixed Model
eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)

## residuals
residuals(eUN.lmm, format = "long", type.residual = c("normalized","pearson"))
residuals(eUN.lmm, format = "long", type.residual = "all", keep.data = TRUE)
residuals(eUN.lmm, format = "wide", plot = "correlation")
residuals(eUN.lmm, format = "wide", type.residual = "normalized")
}
