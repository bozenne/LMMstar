### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: May 10 2021 (15:42) 
##           By: Brice Ozenne
##     Update #: 165
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
##' @param effects [character] Should all coefficients be output (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param strata [character vector] When not \code{NULL}, only output coefficient relative to specific levels of the variable used to stratify the mean and covariance structure.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##' 
##'
##' @details \bold{transform.sigma}: \cr
##' \itemize{
##' \item \code{"none"} ouput residual standard error.
##' \item \code{"log"} ouput log-transformed residual standard error.
##' \item \code{"square"} ouput residual variance.
##' \item \code{"logsquare"} ouput log-transformed residual variance.
##' }
##'
##'  \bold{transform.k}: \cr
##' \itemize{
##' \item \code{"none"} ouput ratio between the residual standard error of the current level and the reference level.
##' \item \code{"log"} ouput log-transformed ratio between the residual standard errors.
##' \item \code{"square"} ouput ratio between the residual variances.
##' \item \code{"logsquare"} ouput log-transformed ratio between the residual variances.
##' \item \code{"sd"} ouput residual standard error of the current level.
##' \item \code{"logsd"} ouput residual log-transformed standard error of the current level.
##' \item \code{"var"} ouput residual variance of the current level.
##' \item \code{"logvar"} ouput residual log-transformed variance of the current level.
##' }
##' 
##'  \bold{transform.rho}: \cr
##' \itemize{
##' \item \code{"none"} ouput correlation coefficient.
##' \item \code{"atanh"} ouput correlation coefficient after tangent hyperbolic transformation.
##' \item \code{"cov"} ouput covariance coefficient.
##' }
##'
##' @return A vector with the value of the model coefficients.

## * coef.lmm (code)
##' @rdname coef
##' @export
coef.lmm <- function(object, effects = "all", type.object = "lmm", strata = NULL,
                     transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    options <- LMMstar.options()
    x.transform.sigma <- object$reparametrize$transform.sigma
    x.transform.k <- object$reparametrize$transform.k
    x.transform.rho <- object$reparametrize$transform.rho
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    
    if(is.null(transform.sigma)){
        transform.sigma <- options$transform.sigma
    }else if(!is.function(transform.sigma)){
        transform.sigma <- match.arg(transform.sigma, c("none","one","log","square","logsquare"))
    }

    if(is.null(transform.k)){
        transform.k <- options$transform.k
    }else if(!is.function(transform.k)){
        transform.k <- match.arg(transform.k, c("none","log","square","logsquare","sd","logsd","var","logvar"))
    }

    if(is.null(transform.rho)){
        transform.rho <- options$transform.rho
    }else if(!is.function(transform.rho)){
        transform.rho <- match.arg(transform.rho, c("none","atanh","cov"))
    }

    if(is.function(transform.sigma) || is.function(transform.k) || is.function(transform.rho)){
        test.notransform <- FALSE
    }else{
        test.notransform <- (transform.sigma==x.transform.sigma) && (transform.k==x.transform.k) && (transform.rho==x.transform.rho)
    }

    if(transform.rho == "cov" && ("variance" %in% effects == FALSE || "correlation" %in% effects == FALSE)){
        stop("Cannot use the argument \'transform.rho\' set to \"cov\" when \"variance\" or \"correlation\" is not in argument \'effect\'. \n")
    }
    
    ## ** extract
    if(type.object=="lmm"){

        out <- NULL
        if("mean" %in% effects){
            out <- c(out, object$param$mu)
        }

        if(any(c("variance","correlation") %in% effects)){
            pVar <- NULL
            if("variance" %in% effects){
                pVar <- c(pVar,object$param$sigma, object$param$k)
            }else if("correlation" %in% effects){
                pVar <- c(pVar,object$param$rho)
            }
            if(test.notransform){
                outVar <- object$reparametrize$p[names(pVar)]
                if(!is.null(object$reparametrize$newname)){
                    newname <- setNames(object$reparametrize$newname[match(names(pVar),names(object$reparametrize$p))], names(pVar))
                }else{
                    newname <- NULL
                }
            }else{
                ls.reparam <- .reparametrize(p = pVar,
                                             type = object$param$type[names(pVar)], strata = object$param$strata[names(pVar)], time.levels = object$time$levels,
                                             Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                             transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
                outVar <- ls.reparam$p
                if(ls.reparam$transform){
                    newname <- setNames(ls.reparam$newname,pVar)
                }else{
                    newname <- NULL
                }
            }
            out <- c(out,outVar)

        }else{
            newname <- NULL
        }

        ## post process
        if(!is.null(strata)){
            out <- out[object$param$strata[names(out)] %in% strata]
        }
        if(length(newname)>0){
            names(out)[match(names(outVar),names(out))] <- as.character(newname)
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
