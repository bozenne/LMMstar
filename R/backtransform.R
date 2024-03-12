### backtransform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul 15 2022 (10:04) 
## Version: 
## Last-Updated: mar 12 2024 (16:56) 
##           By: Brice Ozenne
##     Update #: 54
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .backtransform
##' @title Back-Transformation of Estimates
##' @description Back-transform estimates and confidence intervals (CIs).
##' @noRd
##' 
##' @param object [data.frame] table containing estimate (\code{estimate}), standard errors (\code{se}), lower and upper bounds of confidence intervals (\code{lower} and \code{upper}).
##' @param type.param [character vector] type of each parameter in the table: can be \code{mu}, \code{sigma}, \code{k}, or \code{rho}.
##' @param backtransform.names [character vector] name of each parameter in the table after backtransformation.
##' @param transform.mu,transform.sigma,transform.k,transform.rho [character or function] name of the transformation to be reverted or back-transformation to be apply to the mean, variance, or correlation parameters.
##' @param backtransform [logical vector] whether back-transformation should be apply to the mean, variance, or correlation parameters.
##'
##' @details If the option \code{transform.sigma} and/or  \code{transform.k} is one of \code{"log"}, \code{"logsd"}, \code{"logvar"}, \code{"logsqaure"},
##' the estimate and CIs are transformed back to the original scale by applying the exponential function. 
##' If the option \code{transform.rho} is \code{"atanh"}, the estimate and CIs are transformed back to the original scale by applying the tangent hyperbolic function.
##'
##' @return data.frame and an attribute \code{"message"} identifying that a backtransformation has been performed.
##' 
##' @keywords internal
.backtransform <- function(object, type.param, backtransform.names,
                           backtransform, transform.mu, transform.sigma, transform.k, transform.rho){

    ## ** prepare
    if(inherits(transform.k,"character") && transform.k %in% c("sd","var","logsd","logvar")){
        type.param[type.param=="sigma"] <- "k"
    }
    
    transform <- list(mu = transform.mu,
                      sigma = transform.sigma,
                      k = transform.k,
                      rho = transform.rho)[c("mu","sigma","k","rho") %in% type.param]

    if(length(backtransform)==4){
        transform <- transform[intersect(names(transform),c("mu","sigma","k","rho")[which(backtransform)])]
    }

    todo <- intersect(c("estimate","se", "lower", "upper"), names(object))
    message <- data.frame(matrix(nrow = NROW(transform), ncol = length(todo)+1, dimnames = list(names(transform), c(todo,"FUN"))))

    ## ** backtransform
    for(iType in names(transform)){ ## iType <- names(transform)[1]

        if(is.function(transform[[iType]])){

            iBacktransform <- transform[[iType]]
            if(!is.null(attr(transform[[iType]], "derivative"))) {
                iDbacktransform <- attr(transform[[iType]], "derivative")
            }else{
                iDbacktransform <- function(x){numDeriv::grad(func = transform[[iType]], x = x)}
            }
            message[iType,"FUN"] <- "user-defined"

        }else if(length(transform[[iType]])!=1 || any(!is.character(transform[[iType]]))){
            stop("Incorrect value for argument \'backtransform\'. \n",
                 "Should be a character string (of length 1) or a function. \n")
        }else if(transform[[iType]] %in% c("none","sd","var","square","cov")){

            message[iType,todo] <- FALSE
            next

        }else if(transform[[iType]] %in% c("log","logsd","logvar","logsquare")){

            iBacktransform <- exp
            iDbacktransform <- function(x){exp(x)}
            message[iType,"FUN"] <- "exp"

        }else if(transform[[iType]] == "exp"){

            iBacktransform <- log
            iDbacktransform <- function(x){1/x}
            message[iType,"FUN"] <- "log"

        }else if(transform[[iType]] == "tanh"){

            iBacktransform <- atanh
            iDbacktransform <- function(x){1/(1-x^2)}
            message[iType,"FUN"] <- "atanh"

        }else if(transform[[iType]] == "atanh"){

            iBacktransform <- tanh
            iDbacktransform <- function(x){4/(exp(x)+exp(-x))^2}
            message[iType,"FUN"] <- "tanh"
            
        }else{

            message[iType,"on"] <- list("NA")
            stop("Incorrect value for argument \'backtransform\' (unknown transformation). \n",
                 "Available transformations: \"exp\", \"log\", \"tanh\", \"atanh\". \n")

         }

        iIndex.type <- which(type.param==iType)

        if("se" %in% names(object)){ ## needs to be before estimate otherwise iDbacktransform is applied to the already transformed estimate
            object$se[iIndex.type] <- object$se[iIndex.type]*iDbacktransform(object$estimate[iIndex.type])
            message[iType,"se"] <- TRUE
        }
        if("estimate" %in% names(object)){
            object$estimate[iIndex.type] <- iBacktransform(object$estimate[iIndex.type])
            message[iType,"estimate"] <- TRUE
        }
        if("lower" %in% names(object)){
            object$lower[iIndex.type] <- iBacktransform(object$lower[iIndex.type])
            message[iType,"lower"] <- TRUE
        }
        if("upper" %in% names(object)){
            object$upper[iIndex.type] <- iBacktransform(object$upper[iIndex.type])
            message[iType,"upper"] <- TRUE
        }
        
    }

    ## ** rename
    if(!is.null(backtransform.names)){
        rownames(object) <- backtransform.names
    }
    
    ## ** export
    attr(object,"backtransform") <- message
    return(object)
}


##----------------------------------------------------------------------
### backtransform.R ends here
