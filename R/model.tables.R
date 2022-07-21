### table.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (10:48) 
## Version: 
## Last-Updated: jul 21 2022 (16:59) 
##           By: Brice Ozenne
##     Update #: 21
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.tables.lmm (documentation)
##' @title Statistical Inference for Linear Mixed Model
##' @description Export estimates, standard errors, degrees of freedom, confidence intervals (CIs) and p-values for the mean coefficients of a linear mixed model. 
##' @name model.tables
##'
##' @param x a \code{lmm} object.
##' @param ... arguments to be passed to the \code{confint} method. Should not contain the argument \code{column}.
##' 
##' @details This function simply calls \code{\link{confint}} with a specific value for the argument \code{column}.
##' 
##' @export
model.tables.lmm <- function(x, ...){
    out <- confint(x, ..., columns = c("estimate","se","df","lower","upper","p.value"))
    attr(out, "backtransform") <- NULL
    class(out) <- "data.frame"
    return(out)
}

##' @export
model.tables.Wald_lmm <- function(x, ...){
    out <- confint(x, ..., columns = c("estimate","se","df","lower","upper","p.value"))
    attr(out, "backtransform") <- NULL
    attr(out, "error") <- NULL
    attr(out, "level") <- NULL
    attr(out, "method") <- NULL
    class(out) <- "data.frame"
    return(out)
}

##' @export
model.tables.mlmm <- function(x, ...){
    out <- confint(x, ..., columns = c("estimate","se","df","lower","upper","p.value"))
    attr(out, "backtransform") <- NULL
    class(out) <- "data.frame"
    return(out)
}

##----------------------------------------------------------------------
### model.tables.R ends here
