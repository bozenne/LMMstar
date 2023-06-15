### weights.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct 11 2022 (10:56) 
## Version: 
## Last-Updated: jun 15 2023 (16:58) 
##           By: Brice Ozenne
##     Update #: 19
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * weights.Wald_lmm (documentation)
##' @title Extract Weights Used to Pool Estimates
##' @description Extract weights used to pool estimates.
##' 
##' @param object a \code{Wald_lmm} object, output of \code{anova.lmm}, or \code{rbind.lmm}, or \code{mlmm}.
##' @param method [character] method for combining the estimates, one of \code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.rubin"}.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return a numeric vector whose elements sum to 1.
##'
##' @keywords methods
##'
##' @examples
##' set.seed(10)
##' dL <- sampleRem(250, n.times = 3, format = "long")
##'
##' e.mlmm <- mlmm(Y~X1+X2+X6, repetition = ~visit|id, data = dL,
##'                by = "X4", effects = "X1=0", structure = "CS")
##' weights(e.mlmm, method = "average")
##' weights(e.mlmm, method = "pool.fixse")
##' weights(e.mlmm, method = "pool.se")
##' weights(e.mlmm, method = "pool.gls")

## * weights.Wald_lmm (code)
##' @export
weights.Wald_lmm <- function(object, method, ...){

    valid.method <- c("average","pool.fixse","pool.se","pool.gls","pool.gls1","pool.rubin")
    if(missing(method)){
        stop("Argument \'method\' is missing.\n",
             "Should be one of \"",paste(valid.method, collapse = "\" \""),"\".\n")
    }
    if(length(method)!=1){
        stop("Argument \'method\' should have lenght 1.\n")
    }    
    if(method %in% valid.method == FALSE){
        stop("Argument \'method\' is missing.\n",
             "Should be one of \"",paste(valid.method, collapse = "\" \""),"\".\n")
    }    
    object.confint <- confint(object, method = method, ...)
    out <- stats::setNames(attr(object.confint,"contrast")[1,], names(coef(object)))
    return(out)

}

##----------------------------------------------------------------------
### weights.R ends here
