#' @docType package
#' @title LMMstar package: Helper functions for handling repeated measurements in R
#' @name LMMstar-package
#' 
#' @description Companion R package for the course "Statistical analysis
#' of correlated and repeated measurements for health science
#' researchers" taught by the section of Biostatistics of the University
#' of Copenhagen. It provides functions for computing summary statistics
#' and obtainign graphical displays of longitudinal data, as well as for
#' statistical modeling and statistical inference using multivariate gaussian model.
#'
#' Currently only two types of multivariate gaussian models are available: \itemize{
#' \item one with a compound symmetry structure for the residual variance-covariance matrix. This is equivalent to a random intercept model
#' \item one with a unstructured structure for the residual variance-covariance matrix
#' }
#' In addition, it possible to stratify the residual variance-covariance matrix (and the mean) with respect to a categorical variable.
#'
#' @importFrom ggplot2 autoplot
#' @importFrom nlme getData getVarCov intervals
#' @importFrom stats anova coef confint logLik model.matrix nobs residuals vcov
#' @importFrom lava iid score information 
#' @importFrom sandwich estfun
#' @importFrom emmeans emm_basis recover_data 
NULL


