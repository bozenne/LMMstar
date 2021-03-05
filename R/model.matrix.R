### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: mar  5 2021 (21:53) 
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


## * model.matrix.lmm (code)
##' @export
model.matrix.lmm <- function(object, data = NULL, type = "lmm"){
    type <- match.arg(type, c("lmm","gls"))
    if(type=="lmm"){
        if(is.null(data)){
            return(object$X)
        }else{
            return(model.matrix(object$mean.structure$formula, data = data))
        }
    }else if(type=="gls"){
        if(is.null(object$variable$strata)){
            return(model.matrix(object$gls[[1]], data = data))
        }else{
            stop("Not implemented yet!")
        }
    }
}



##----------------------------------------------------------------------
### model.matrix.R ends here
