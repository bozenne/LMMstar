### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: mar  5 2021 (21:36) 
##           By: Brice Ozenne
##     Update #: 9
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
coef.lmm <- function(object, effects = c("mean"), type = "lmm"){
    type <- match.arg(type, c("lmm","gls"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)

    if(type=="lmm"){
        browser()
        return(object$vcov)
    }else if(type=="gls"){
        if(is.null(object$variable$strata)){
            return(coef(object$gls[[1]]))
        }else{
            return(lapply(object$gls, coef))
        }
    }
}


##----------------------------------------------------------------------
### coef.R ends here
