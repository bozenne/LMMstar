### formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:53) 
## Version: 
## Last-Updated: May 27 2021 (12:11) 
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

## * formula.lmm (code)
##' @export
formula.lmm <- function(x, type.object = "lmm", effects = "mean", ...){

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    type.x <- match.arg(type.object, c("lmm","lmm-mean","lmm-variance","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)

    ## ** extract formula
    if(type.object == "lmm"){
        if("mean" %in% effects && "variance" %in% effects){
            return(list(mean = x$formula$mean,
                        variance = x$formula$var))
        }else if("mean" %in% effects){
            return(x$formula$mean)
        }else if("variance" %in% effects){
            return(x$formula$var)
        }
    }else if(type.object =="gls"){
        if(x$strata$n==1){
            return(stats::formula(x$gls[[1]]))
        }else{
            return(lapply(x$gls, stats::formula))
        }
    }
}


##----------------------------------------------------------------------
### formula.R ends here
