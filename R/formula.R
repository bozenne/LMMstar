### formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:53) 
## Version: 
## Last-Updated: May 10 2021 (19:10) 
##           By: Brice Ozenne
##     Update #: 13
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
formula.lmm <- function(object, type.object = "lmm", effects = "mean"){
    ## 
    type.object <- match.arg(type.object, c("lmm","lmm-mean","lmm-variance","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)

    if(type.object == "lmm"){
        if("mean" %in% effects && "variance" %in% effects){
            return(list(mean = object$formula$mean,
                        variance = object$formula$var))
        }else if("mean" %in% effects){
            return(object$formula$mean)
        }else if("variance" %in% effects){
            return(object$formula$var)
        }
    }else if(type.object =="gls"){
        if(object$strata$n==1){
            return(stats::formula(object$gls[[1]]))
        }else{
            return(lapply(object$gls, stats::formula))
        }
    }
}


##----------------------------------------------------------------------
### formula.R ends here
