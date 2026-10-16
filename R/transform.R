### transform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 15 2026 (09:53) 
## Version: 
## Last-Updated: okt 16 2026 (15:28) 
##           By: Brice Ozenne
##     Update #: 19
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
transform.lmm <- function(`_data`, type = NULL, ...){

    object <- `_data`

    if(is.null(type)){
        available.transform <- paste("transform",c("sigma","k","rho"),sep=".")
        ## grep("transform",names(object$args$control),value  = TRUE)
        out <- object$args[available.transform]
        names(out) <- c("sigma","k","rho")
    }else if(type %in% c("sigma","k","rho")){
        out <- object$args[[paste("transform",type,sep=".")]]
    }else if(type == "name"){
        out <- names(object$opt$estimate)
    }else if(type == "value(originalName)"){
        out <- stats::setNames(object$opt$estimate, names(object$param))
    }else if(type == "value(transName)"){
        out <- object$opt$estimate
    }else if(type == "index"){
        out <- which(names(object$opt$estimate) != names(object$param))
    }else{
        stop("Argument \'type\' should be one of \"transform.sigma\", \"transform.k\", \"transform.rho\", \n",
             "\"name\", \"value(originalName)\", \"value(transName)\", \"index\". \n",
             "Proposed value: \"",type,". \n")
    }
    
    return(out)
}

##----------------------------------------------------------------------
### transform.R ends here
