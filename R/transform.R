### transform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 15 2026 (09:53) 
## Version: 
## Last-Updated: okt 15 2026 (11:18) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @description Extract the transform.sigma, transform.k, transform.rho.
##' For internal use.
##' @noRd
##' @keywords internal

## * transform.lmm (code)
transform.lmm <- function(object, type = NULL){


    if(is.null(type)){
        available.transform <- paste("transform",c("sigma","k","rho"),sep=".")
        ## grep("transform",names(object$args$control),value  = TRUE)
        out <- object$args$control[available.transform]
        names(out) <- c("sigma","k","rho")
    }else if(type %in% c("sigma","k","rho")){
        out <- object$args$control[[paste("transform",type,sep=".")]]
    }else if(type == "name"){
        out <- names(object$opt$estimate)
    }else if(type == "value(originalName)"){
        out <- stats::setNames(object$opt$estimate, names(object$param))
    }else if(type == "value(transName)"){
        out <- object$opt$estimate
    }else if(type == "index"){
        out <- which(names(object$opt$estimate) != names(object$param))
    }
    
    return(out)
}

##----------------------------------------------------------------------
### transform.R ends here
