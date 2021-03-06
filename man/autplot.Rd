% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoplot.R
\name{autoplot}
\alias{autoplot}
\alias{autoplot.lmm}
\title{Graphical Display For Linear Mixed Models}
\usage{
\method{autoplot}{lmm}(
  object,
  at = NULL,
  color = TRUE,
  ci = TRUE,
  alpha = NA,
  plot = TRUE,
  size.point = 3,
  size.line = 1,
  size.text = 16,
  position.errorbar = "identity",
  ...
)
}
\arguments{
\item{object}{a \code{lmm} object.}

\item{at}{[data.frame] values for the covariates at which to evaluate the fitted values.}

\item{color}{[character] name of the variable in the dataset used to color the curve.}

\item{ci}{[logical] should confidence intervals be displayed?}

\item{alpha}{[numeric, 0-1] When not NA, transparency parameter used to display the confidence intervals.}

\item{plot}{[logical] should the plot be displayed?}

\item{size.point}{[numeric, >0] the size of the point on the plot.}

\item{size.line}{[numeric, >0] the size of the line on the plot.}

\item{size.text}{[numeric, >0] Size of the font used to displayed text when using ggplot2.}

\item{position.errorbar}{[character] relative position of the errorbars.}

\item{...}{Not used. For compatibility with the generic method.}
}
\description{
Graphical Display For Linear Mixed Models
}
