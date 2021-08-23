### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:41) 
## Version: 
## Last-Updated: aug 23 2021 (17:18) 
##           By: Brice Ozenne
##     Update #: 14
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
nobs.lmm <- function(object, type.object = "lmm", ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    type.object <- match.arg(type.object, c("lmm","gls"))

    ## ** extract
    if(type.object == "lmm"){
        return(c(obs = sum(object$design$cluster$nobs),
                 cluster = object$design$cluster$n,
                 missing = length(object$index.na)))
    }else if(type.object == "gls"){
        if(object$strata$n==1){
            stats::nobs(object$gls[[1]])
        }else{
            lapply(object$gls,stats::nobs)
        }
    }
}


##----------------------------------------------------------------------
### nobs.R ends here
