### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: Apr 22 2021 (17:19) 
##           By: Brice Ozenne
##     Update #: 97
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
    if(is.null(transform)){
        transform <- options$transform
    }else if(transform %in% c(0,1,2) == FALSE){
        stop("Argument \'transform\' must be 0 (standard error parameters, correlation parameters), \n",
             "                               1 (log transformation of the standard error parameters, atanh transformation of the correlation parameters), \n",
             "                               2 (variance parameters, correlation parameters). \n")
    }

    ## ** extract
    if(type.object=="lmm"){

        out <- NULL
        if("mean" %in% effects){
            out <- c(out, object$param$mu)
        }
        if("variance" %in% effects){
            res.coefVar <- .coefVar(sigma = object$param$sigma, k = object$param$k, rho = object$param$rho,
                                    transform = transform, transform.names = transform.names,
                                    param.type = object$param$type, param.strata = object$param$strata,
                                    time.var = object$time$var, time.levels = object$time$levels)
            out <- c(out,res.coefVar)
            rename <- attr(res.coefVar, "rename")
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
        if(transform>0){
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

## * .coef
## apply transformation or inverse transformation to the variance coefficients
.coefVar <- function(sigma, k, rho, transform, transform.names,
                     param.type, param.strata, time.var, time.levels){


    rename <- NULL
    out <- NULL

    ## ** sigma
    if(transform == 0){
        out <- sigma
    }else  if(transform == 1){
        out <- log(sigma)
        if(transform.names){
            rename <- setNames(paste0("log(",names(sigma),")"),names(out))
        }
    }else if(transform == 2 && length(k)==0){
        out <- sigma^2
        if(transform.names){
            rename <- setNames(paste0(names(sigma),"^2"),names(out))
        }
    }else if(transform == -1){ ## back-transform (internal use only)
        out <- exp(sigma)
    }else if(transform == -2 && length(k)==0){ ## back-transform (internal use only)
        out <- sqrt(sigma)
    }

    ## ** k
    if(length(k)>0){
        if(transform == 0){
            out.k <-  k
        }else if(transform == 1){
            out.k <-  log(k)
            if(transform.names){
                rename <- c(rename,setNames(paste0("log(",names(k),")"),names(out.k)))
            }
        }else if(transform == 2){
            strata.type <- param.type[param.type %in% c("sigma","k")]
            strata.index <- param.strata[param.type %in% c("sigma","k")]
            n.strata <- unique(strata.index)
            out.k <- NULL
            for(iStrata in 1:n.strata){ ## iStrata <- 1 
                iType <- strata.type[strata.index==iStrata]
                iName.sigma <- names(iType[iType=="sigma"])
                iName.k <- names(iType[iType=="k"])
                iOut.k <- sigma[iName.sigma]^2*c(1,k[iName.k]^2)
                if(transform.names){
                    rename <- c(rename,setNames(paste0(iName.sigma,"^2:",time.var,time.levels),names(iOut.k)))
                }
                out.k <- c(out.k, iOut.k)
            }
        }else if(transform == -1){ ## back-transform (internal use only)
            out.k <-  exp(k)
        }else if(transform == -2){ ## back-transform (internal use only)
            strata.type <- param.type[param.type %in% c("sigma","k")]
            strata.index <- param.strata[param.type %in% c("sigma","k")]
            n.strata <- unique(strata.index)
            out.k <- NULL
            for(iStrata in 1:n.strata){ ## iStrata <- 1 
                iType <- strata.type[strata.index==iStrata]
                iName.sigma <- names(iType[iType=="sigma"])
                iName.k <- names(iType[iType=="k"])
                iOut.k <- sqrt(c(1,k[iName.k])/sigma[iName.sigma])
                out.k <- c(out.k, iOut.k)
            }
        }
        out <- c(out, out.k)
    }

    ## ** rho
    if(length(rho)>0){
        if(transform == 0){
            out.rho <- rho
        }else if(transform == 1){
            out.rho <- atanh(rho)
            if(transform.names){
                rename <- c(rename,setNames(paste0("atanh(",names(rho),")"),names(iOut.rho)))
            }
        }else if(transform == -1){ ## back-transform (internal use only)
            out.rho <- tanh(rho)
        }
        out <- c(out, out.rho)
    }

    ## ** export
    if(transform>0 && transform.names){
        attr(out,"rename") <- rename
    }
    return(out)
}
##----------------------------------------------------------------------
### coef.R ends here
