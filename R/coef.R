### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: Apr 25 2021 (18:34) 
##           By: Brice Ozenne
##     Update #: 107
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * coef.lmm (documentation)
##' @title Extract Coefficients From a Linear Mixed Model
##' @description Extract coefficients from a linear mixed model.
##' @name coef
##'
##' @param object a \code{lmm} object.
##' @param effects [character] Should only coefficients relative to the mean (\code{"mean"})
##' or only coefficients relative to the variance-covariance structure (\code{"variance"}) be output, or both (\code{all}).
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param strata [character vector] When not \code{NULL}, only output coefficient relative to specific levels of the variable used to stratify the mean and covariance structure.
##' @param transform [0,1,2] Transformation used on the variance coefficient. See details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##' 
##'
##' @details \bold{transform}: \cr
##' \itemize{
##' \item 0 means no transformation i.e. ouput stanrdard error, ratio of standard errors, and correlations.
##' \item 1 means log/atanh transformation i.e. ouput log(stanrdard error), log(ratio of standard errors), and atanh(correlations).
##' \item 2 ouput variance coefficients and correlations.
##' }
##'
##' @return A vector with the value of the model coefficients.

## * coef.lmm (code)
##' @rdname coef
##' @export
coef.lmm <- function(object, effects = "all", type.object = "lmm", strata = NULL, transform = NULL, transform.names = TRUE, ...){
    options <- LMMstar.options()
        
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    
    ## ** extract
    if(type.object=="lmm"){

        out <- NULL
        if("mean" %in% effects){
            out <- c(out, object$param$mu)
        }
        if("variance" %in% effects){
            outVar <- do.call(reparametrize,
                               args = c(list(p = c(object$param$sigma, object$param$k, object$param$rho),
                                             param.type = object$param$type, param.strata = object$param$strata, time.levels = object$time$levels,
                                             Jacobian = FALSE, dJacobian = FALSE), transform))
            out <- c(out,outVar)
            rename <- attr(outVar, "rename")
        }else{
            rename <- NULL
        }

        ## post process
        if(!is.null(strata)){
            out <- out[object$param$strata[names(out)] %in% strata]
        }
        if(length(rename)>0){
            names(out)[match(names(rename),names(out))] <- as.character(rename)
        }

        return(out)

    }else if(type.object=="gls"){
        if(!is.null(transform)){
            stop("Cannot handle argument \'transform\' when argument \'type.object\' is \"gls\". \n")
        }
        if(length(effects)!=1 || "variance" %in% effects){
            stop("Cannot handle argument \'effects\' when argument \'type.object\' is \"gls\". \n")
        }

        if(is.null(strata) && is.null(object$variable$strata)){
            return(coef(object$gls[[1]]))
        }else{
            return(lapply(object$gls[which(object$strata$levels %in% strata)], coef))
        }

    }
}
##----------------------------------------------------------------------
### coef.R ends here
