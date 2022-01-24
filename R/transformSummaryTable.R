### transformSummaryTable.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  3 2020 (18:29) 
## Version: 
## Last-Updated: jan 24 2022 (10:45) 
##           By: Brice Ozenne
##     Update #: 29
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @title Apply Transformation to Summary Table
#' @description Update summary table according to a transformation, e.g. log-transformtion.
#' P-values are left unchanged but estimates, standard errors, and confidence intervals are updated.
#'
#' @param object A data.frame with columns estimate, se, lower, upper.
#' @param transform the name of a transformation or a function.
#'
#' @return a data.frame
#' @export
transformSummaryTable <- function(object, transform = NULL){
    if(is.null(transform)){
        return(object)
    }else if(identical(transform,"atanh")){
        transform <- atanh
        dtransform <- function(x){1/(1-x^2)}
    }else if(identical(transform,"exp")){
        transform <- exp
        dtransform <- function(x){exp(x)}
    }else if(identical(transform,"log")){
        transform <- log
        dtransform <- function(x){1/x}
    }else if(identical(transform,"loglog")){
        transform <- function(x){log(log(x))}
        dtransform <- function(x){1/(-x*log(x))}
    }else if(identical(transform,"cloglog")){
        transform <- function(x){log(log(1-x))}
        dtransform <- function(x){1/(-(1-x)*log(1-x))}
    }else if(!is.null(attr(transform,"derivative"))){
        dtransform <- attr(transform,"derivative")
    }else{
        dtransform <- function(x){diag(numDeriv::jacobian(transform, x))}
    }
    object[,"se"] <- object[,"se"]*dtransform(object[,"estimate"])
    object[,"estimate"] <- transform(object[,"estimate"])
    if("lower" %in% names(object)){
        object[,"lower"] <- transform(object[,"lower"])
    }
    if("upper" %in% names(object)){
        object[,"upper"] <- transform(object[,"upper"])
    }
    return(object)
}

######################################################################
### transformSummaryTable.R ends here
