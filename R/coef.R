### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: mar 22 2021 (10:48) 
##           By: Brice Ozenne
##     Update #: 23
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * coef.lmm (code)
##' @export
coef.lmm <- function(object, effects = "mean", type = "lmm", strata = NULL){
    
    ## ** normalize user imput
    type <- match.arg(type, c("lmm","gls"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }

    ## ** extract
    if(type=="lmm"){

        out <- NULL
        if("mean" %in% effects){
            out <- c(out, object$param$mu)
        }
        if("variance" %in% effects){
            out <- c(out, object$param$sigma, object$param$k, object$param$cor)
        }
        if(!is.null(strata)){
            out <- out[object$param$strata[names(out)] %in% strata]
        }
        return(out)

    }else if(type=="gls"){

        if(is.null(strata) && is.null(object$variable$strata)){
            return(coef(object$gls[[1]]))
        }else{
            return(lapply(object$gls[which(object$strata$levels %in% strata)], coef))
        }

    }
}


##----------------------------------------------------------------------
### coef.R ends here
