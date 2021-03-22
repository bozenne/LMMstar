### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:41) 
## Version: 
## Last-Updated: mar 22 2021 (22:18) 
##           By: Brice Ozenne
##     Update #: 8
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
nobs.lmm <- function(object, type = "lmm"){
    
    ## ** normalize user imput
    type <- match.arg(type, c("lmm","lmm-obs","lmm-cluster","gls"))

    ## ** extract
    if(type == "lmm"){
        return(c(obs = sum(object$design$cluster$nobs),
                 cluster = object$design$cluster$n))
    }else if(type == "lmm-obs"){
        return(sum(object$design$cluster$nobs))
    }else if(type == "lmm-cluster"){
        return(object$design$cluster$n)
    }else if(type == "gls"){
        if(object$strata$n==1){
            stats::nobs(object$gls[[1]])
        }else{
            lapply(object$gls,stats::nobs)
        }
    }
}


##----------------------------------------------------------------------
### nobs.R ends here
