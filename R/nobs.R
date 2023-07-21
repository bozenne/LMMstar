### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:41) 
## Version: 
## Last-Updated: jul 21 2023 (17:31) 
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

## * nobs.lmm
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
