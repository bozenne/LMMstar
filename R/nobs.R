### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:41) 
## Version: 
## Last-Updated: mar 12 2024 (11:35) 
##           By: Brice Ozenne
##     Update #: 24
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * nobs.lmm
##' @title Number of Observations from a Linear Mixed Model
##' @description Extract the number of observations from a Linear Mixed Model
##'
##' @param object an lmm object.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A vector with 4 elements: \itemize{
##' \item \code{obs}: the number of repetitions with full data
##' \item \code{cluster}: the number of clusters with a least one repetition with full data
##' \item \code{missing.obs}: the number of repetitions with missing data
##' \item \code{missing.cluster}: the number of cluster with only missing data
##' }
##' @keywords methods
##' @export
nobs.lmm <- function(object, ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** extract
    return(c(obs = length(object$design$Y), ## number of observations without missing data
             cluster = length(object$design$index.cluster), ## number of cluster with a least one observed data
             missing.obs = length(object$index.na),
             missing.cluster = sum(is.na(object$cluster$index))
             )
           )
}

## * nobs.mlmm
##' @export
nobs.mlmm <- function(object, ...){

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** extract
    out <- do.call(rbind,lapply(object$model, stats::nobs))

    ## ** export
    return(out)
}
##----------------------------------------------------------------------
### nobs.R ends here
