### confint.mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2022 (14:42) 
## Version: 
## Last-Updated: jun 24 2022 (17:26) 
##           By: Brice Ozenne
##     Update #: 20
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.anova_lmm
##' @title Confidence Intervals for Multiple Linear Mixed Model.
##' @description Compute confidence intervals for several linear mixed models.
##' 
##' @param object an \code{mlmm} object, output of \code{mlmm}.
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}, or \code{"pool"}.
##' @param ... other arguments are passed to \code{\link{confint.anova_lmm}}.
##'
##' @details Statistical inference following pooling is performed according to Rubin's rule whose validity requires the congeniality condition of Meng (1994).
##'
##' @references
##' Meng X. L.(1994). Multiple-imputation inferences with uncongenial sources of input. Statist. Sci.9, 538–58.
##'
##' @export
confint.mlmm <- function(object, parm = NULL, level = 0.95, method = NULL, ...){

    if(is.null(method)){
        method <- "none"
    }

    return(confint.anova_lmm(object, parm = parm, level = level, method = method, ...))
    
}


##----------------------------------------------------------------------
### confint.mlmm.R ends here
