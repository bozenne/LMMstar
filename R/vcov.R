### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: mar  5 2021 (21:38) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * vcov.lmm
vcov.lmm <- function(object, type = "lmm"){
    type <- match.arg(type, c("lmm","gls"))

    if(type=="lmm"){
        return(object$betavcov)
    }else if(type=="gls"){
        if(is.null(object$variable$strata)){
            return(vcov(object$gls[[1]]))
        }else{
            return(lapply(object$gls, vcov))
        }
    }
}

##----------------------------------------------------------------------
### vcov.R ends here
