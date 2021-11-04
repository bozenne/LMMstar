### formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:53) 
## Version: 
## Last-Updated: nov  4 2021 (10:32) 
##           By: Brice Ozenne
##     Update #: 16
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
formula.lmm <- function(x, effects = "mean", ...){

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)

    ## ** extract formula
    if("mean" %in% effects && "variance" %in% effects){
        return(list(mean = x$formula$mean,
                    variance = x$formula$var))
    }else if("mean" %in% effects){
        return(x$formula$mean)
    }else if("variance" %in% effects){
        return(x$formula$var)
    }
}


##----------------------------------------------------------------------
### formula.R ends here
