#' @docType package
#' @title LMMstar package: repeated measurement models for discrete times
#' @name LMMstar-package
#' 
#' @description Companion R package for the course "Statistical analysis
#' of correlated and repeated measurements for health science
#' researchers" taught by the section of Biostatistics of the University
#' of Copenhagen. It implements linear mixed models where the model for the variance-covariance of the residuals
#' is specified via patterns (compound symmetry, unstructured, ...). Statistical inference for mean, variance, and correlation parameters
#' is performed based on the observed information and a Satterthwaite degrees of freedom.
#' Normalized residuals are provided to assess model misspecification.
#' Statistical inference can be performed for arbitrary linear or non-linear combination(s) of model coefficients.
#' Predictions can be computed conditional to covariates only or also to outcome values. 
#'
#' \strong{Notations}: the linear mixed model is denoted:
#''\deqn{ \mathbf{Y}_{i} = \mathbf{X}_{i}\beta+\boldsymbol{\varepsilon}_i }
##' where
##' \itemize{
##' \item \eqn{\mathbf{Y}}: vector of outcome.
##' \item \eqn{\mathbf{X}}: design matrix.
##' \item \eqn{\boldsymbol{\varepsilon}}: vector of residuals with 0-mean and variance \eqn{\Omega}.
##' \item \eqn{\beta}: estimated mean coefficients relative to \eqn{X}.
##' \item \eqn{\Omega}: the modeled variance-covariance of the residuals with diagonal elements \eqn{\omega}.
##' \item \eqn{i} indexes the cluster (level where replicates are assumed independent).
##' \item \eqn{j} indexes the repetitions, e.g. the variance of \eqn{\varepsilon_{ij}} is \eqn{\omega_{ij}}.
##' }
#'
#' \strong{Covariance patterns}: \eqn{\Omega} can be parametrized as: \itemize{
#' \item \code{"ID"}: identity (no correlation, constant variance).
#' \item \code{"IND"}: independent (no correlation, time-specific variance).
#' \item \code{"CS"}: compound symmetry (constant correlation and variance). Can also be used to specify a nested random effect structure or a block specific correlation and variance.
#' \item \code{"UN"}: unstructured (time-specific correlation, time-specific variable).
#' }
#' It possible to stratify the last two structure with respect to a categorical variable. \cr \cr
#'
#' \strong{Optimizer}: the default optimizer is \code{nlme::gls} which is restricted to certain covariance patterns.
#' To use the other covariance patterns switch to the optimizer \code{"FS"}.
#' This may fail for complex covariance patterns in small samples since \eqn{\Omega} is not constrained to be positive definite.
#' 
#' @importFrom ggplot2 autoplot
#' @importFrom nlme getData getVarCov intervals
#' @importFrom stats anova coef confint df dummy.coef fitted logLik model.matrix model.tables nobs residuals vcov
#' @importFrom lava estimate iid score information 
#' @importFrom sandwich estfun
#' @importFrom emmeans emm_basis recover_data 
NULL


