#' @docType package
#' @title LMMstar package: repeated measurement models for discrete times
#' @name LMMstar-package
#' 
#' @description Companion R package for the course "Statistical analysis
#' of correlated and repeated measurements for health science
#' researchers" taught by the section of Biostatistics of the University
#' of Copenhagen. It implements linear mixed models where the model for the variance-covariance of the residuals
#' is specified via patterns (compound symmetry, unstructured). Statistical inference for mean, variance, and correlation parameters
#' is performed based on the observed information and a Satterthwaite degrees of freedom.
#' Normalized residuals are provided to assess model misspecification.
#' Statistical inference can be performed for arbitrary linear combination(s) of model coefficients.
#' Predictions can be computed conditional to covariates only or also to outcome values. 
#'
#' Currently only four types of model for the residual variance-covariance matrix are available: \itemize{
#' \item \code{"ID"}: Identity (no correlation, constant variance)
#' \item \code{"IND"}: Independent (no correlation, time-specific variance)
#' \item \code{"CS"}: compound symmetry (constant correlation, constant variable)
#' \item \code{"UN"}: unstructured (time-specific correlation, time-specific variable)
#' }
#' It possible to stratify the last two structure with respect to a categorical variable.
#'
#' The package is based on the \code{nlme::gls} function and the PROC MIXED from the SAS software.
#' Adjustment for multiple comparisons is based on the multcomp package.
#'
#' @importFrom ggplot2 autoplot
#' @importFrom nlme getData getVarCov intervals
#' @importFrom stats anova coef confint df dummy.coef fitted logLik model.matrix model.tables nobs residuals vcov
#' @importFrom lava iid score information 
#' @importFrom sandwich estfun
#' @importFrom emmeans emm_basis recover_data 
NULL


