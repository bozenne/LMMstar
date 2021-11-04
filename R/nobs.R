### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:41) 
## Version: 
## Last-Updated: nov  4 2021 (10:34) 
##           By: Brice Ozenne
##     Update #: 15
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
    return(c(obs = sum(object$design$cluster$nobs),
             cluster = object$design$cluster$n,
             missing = length(object$index.na)))
}


##----------------------------------------------------------------------
### nobs.R ends here
