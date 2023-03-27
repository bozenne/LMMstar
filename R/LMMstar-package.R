##' @docType package
##' @title LMMstar package: repeated measurement models for discrete times
##' @name LMMstar-package
##' 
##' @description Companion R package for the course "Statistical analysis
##' of correlated and repeated measurements for health science
##' researchers" taught by the section of Biostatistics of the University
##' of Copenhagen. It implements linear mixed models where the model for the variance-covariance of the residuals
##' is specified via patterns (compound symmetry, toeplitz, unstructured, ...). Statistical inference for mean, variance, and correlation parameters
##' is performed based on the observed information and a Satterthwaite approximation of the degrees of freedom.
##' Normalized residuals are provided to assess model misspecification.
##' Statistical inference can be performed for arbitrary linear or non-linear combination(s) of model coefficients.
##' Predictions can be computed conditional to covariates only or also to outcome values. \cr \cr
##'
##' \strong{Notations}: the linear mixed model estimated by \code{\link{lmm}} is denoted:
##' \deqn{ \mathbf{Y}_{i} = \mathbf{X}_{i}\beta+\boldsymbol{\varepsilon}_i }
##' where
##' \itemize{
##' \item \eqn{\mathbf{Y}=(Y_1,\ldots,Y_m)}: vector of outcomes.
##' \item \eqn{\mathbf{X}=(X_1,\ldots,X_p)}: design matrix (extractor: \code{model.matrix.lmm}).
##' \item \eqn{\boldsymbol{\varepsilon}}: vector of residuals with 0-mean and variance \eqn{\Omega_i} (extractor: \code{\link{residuals.lmm}}).
##' \item \eqn{\beta}: estimated mean coefficients relative to \eqn{X} (extractor: \code{\link{coef.lmm}}).
##' \item \eqn{\Omega}: the modeled variance-covariance of the residuals with diagonal elements \eqn{\sigma^2_{j}}  (extractor: \code{\link{sigma.lmm}}).
##' \item \eqn{i} indexes the cluster (level where replicates are assumed independent).
##' \item \eqn{j} indexes the repetitions, e.g. the variance of \eqn{\varepsilon_{ij}} is \eqn{\sigma^2_{ij}}.
##' }
##'
##' \strong{Covariance patterns}: \eqn{\Omega} can be parametrized as: \itemize{
##' \item \code{"ID"}: identity (no correlation, constant variance).
##' \item \code{"IND"}: independent (no correlation, time-specific variance).
##' \item \code{"CS"}: compound symmetry (constant correlation and variance). Can also be used to specify a nested random effect structure or a block specific correlation and variance.
##' \item \code{"TOEPLITZ"}: toeplitz (lag-specific correlation, time-specific variance).
##' \item \code{"UN"}: unstructured (time-specific correlation, time-specific variance).
##' }
##' It possible to stratify each structure with respect to a categorical variable. \cr \cr
##'
##' \strong{Optimizer}: the default optimizer, \code{"FS"}, implements a fisher scoring algorithm descent with back-tracking in case of decreasing or undefined log-likelihood.
##' It does not constrain \eqn{\Omega} to be positive definite which may cause problem in small sample or complex models.
##' It is possible to use other optimizer: \code{nlme::gls} for certain covariance patterns or \code{stats::optim}.
##' 
##' @importFrom ggplot2 autoplot
##' @importFrom rlang .data
##' @importFrom nlme getData getVarCov intervals
##' @importFrom stats anova coef confint df dummy.coef fitted logLik model.matrix model.tables nobs profile residuals sigma vcov weights
##' @importFrom lava bootstrap estimate iid information manifest score
##' @importFrom sandwich estfun
##' @importFrom emmeans emm_basis recover_data 
NULL


