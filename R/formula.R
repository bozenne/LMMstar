### formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:53) 
## Version: 
## Last-Updated: mar  5 2021 (21:55) 
##           By: Brice Ozenne
##     Update #: 1
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
    type <- match.arg(type, c("lmm","lmm-mean","lmm-variance","gls"))
    if(type=="lmm"){
        return(list(mean = object$mean.structure$formula,
                    variance = object$variance.structure$formula))
    }else if(type == "lmm-mean"){
        return(object$mean.structure$formula)
    }else if(type == "lmm-variance.structure"){
        return(object$variance.structure$formula)
    }else if(type=="gls"){
        if(is.null(object$variable$strata)){
            return(stats::formula(object$gls[[1]]))
        }else{
            return(lapply(object$gls, stats::formula))
        }
    }
}


##----------------------------------------------------------------------
### formula.R ends here
