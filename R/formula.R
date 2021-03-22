### formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:53) 
## Version: 
## Last-Updated: mar 22 2021 (22:07) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * formula.lmm (code)
##' @export
formula.lmm <- function(object, type = "lmm"){
## 
    type <- match.arg(type, c("lmm","lmm-mean","lmm-variance","gls"))

    if(type=="lmm"){
        return(list(mean = object$formula$mean,
                    variance = object$formula$var))
    }else if(type == "lmm-mean"){
        return(object$formula$mean)
    }else if(type == "lmm-variance"){
        return(object$formula$var)
    }else if(type=="gls"){
        if(object$strata$n==1){
            return(stats::formula(object$gls[[1]]))
        }else{
            return(lapply(object$gls, stats::formula))
        }
    }
}


##----------------------------------------------------------------------
### formula.R ends here
