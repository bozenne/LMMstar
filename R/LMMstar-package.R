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
##' @details \strong{Notations}: the linear mixed model estimated by \code{\link{lmm}} is denoted:
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
##' \item \code{\link{ID}}: identity (no correlation, constant variance).
##' \item \code{\link{IND}}: independent (no correlation, time-specific variance).
##' \item \code{\link{CS}}: compound symmetry (constant correlation and variance). Can also be used to specify a nested random effect structure or a block specific correlation and variance.
##' \item \code{\link{RE}}: random effects.
##' \item \code{\link{TOEPLITZ}}: toeplitz (lag-specific correlation, time-specific variance).
##' \item \code{\link{UN}}: unstructured (time-specific correlation, time-specific variance).
##' }
##' It possible to stratify each structure with respect to a categorical variable. \cr \cr
##'
##' \strong{Optimizer}: the default optimizer, \code{"FS"}, implements a fisher scoring algorithm descent with back-tracking in case of decreasing or undefined log-likelihood.
##' It does not constrain \eqn{\Omega} to be positive definite which may cause problem in small sample or complex models.
##' It is possible to use other optimizer inferfaced by \code{optimx::optimx}.
##' 
##' \strong{Keywords}: documented methods/functions are classified according to the following keywords \itemize{
##' \item models: function fitting a statistical model based on a dataset (e.g. \code{\link{lmm}}, \code{\link{lmmCC}}, \code{\link{mlmm}}, \code{\link{mt.test}}, \code{\link{partialCor}})
##' \item htest: methods performing statistical inference based on an existing model (e.g. \code{\link{anova.lmm}}, \code{\link{estimate.lmm}}, \code{\link{effects.lmm}}, \code{\link{profile.lmm}}, \code{\link{proportion}}, \code{\link{rbind.Wald_lmm}}, \code{\link{resample.lmm}})
##' \item methods: extractors (e.g. \code{\link{coef.lmm}}, \code{\link{confint.lmm}}, \code{\link{df.residual.lmm}}, \code{\link{fitted.lmm}}, \code{\link{iid.lmm}}, \code{\link{information.lmm}}, \code{\link{levels.lmm}}, \code{\link{logLik.lmm}}, \code{\link{manifest.lmm}}, \code{\link{model.frame.lmm}}, \code{\link{model.matrix.lmm}}, \code{\link{model.tables.lmm}}, \code{\link{nobs.lmm}}, \code{\link{predict.lmm}}, \code{\link{ranef.lmm}}, \code{\link{residuals.lmm}}, \code{\link{score.lmm}}, \code{\link{sigma.lmm}}, \code{\link{summary.lmm}}, \code{\link{vcov.lmm}}, \code{\link{weights.Wald_lmm}})
##' \item utilities: function used to facilitate the user interface (e.g. \code{\link{add}}, \code{\link{baselineAdjustment}}, \code{\link{LMMstar.options}}, \code{\link{remove}}, \code{\link{scatterplot}}, \code{\link{summarize}}, \code{\link{summarizeNA}})
##' \item datasets: dataset stored in the package  (e.g. \code{\link{abetaW}})
##' \item hplot: graphical display (e.g. \code{\link{autoplot.lmm}} or \code{\link{plot.lmm}})
##' \item datagen: function for generating data sets (e.g. \code{\link{sampleRem}})
##' \item multivariate: covariance patterns (e.g. \code{\link{ID}}, \code{\link{IND}}, \code{\link{CS}}, \code{\link{RE}}, \code{\link{TOEPLITZ}}, \code{\link{UN}}, \code{\link{CUSTOM}})
##' }
##' 
##' @importFrom ggplot2 autoplot
##' @importFrom rlang .data
##' @importFrom nlme ranef
##' @importFrom stats anova coef confint df dummy.coef effects fitted logLik model.matrix model.tables nobs profile residuals sigma vcov weights
##' @importFrom lava bootstrap estimate iid information manifest score
##' @keywords internal 
"_PACKAGE"


