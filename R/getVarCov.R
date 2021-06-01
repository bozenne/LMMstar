### getVarCov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: jun  1 2021 (16:24) 
##           By: Brice Ozenne
##     Update #: 75
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getVarCov.lmm (documentation)
##' @title Extract The Residuals Variance-Covariance Matrix From a Linear Mixed Model
##' @description Extract the unique set of residuals variance-covariance matrices or the one relative to specific clusters.
##' @name getVarCov
##' 
##' @param object a \code{lmm} object.
##' @param individual [character] identifier of the cluster for which to extract the residual variance-covariance matrix.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the residual variance-covariance matrix. Only relevant if differs from the fitted values.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param strata [character vector] When not \code{NULL} and argument \code{individual} is not specified, only output the residual variance-covariance matrix relative to specific levels of the variable used to stratify the mean and covariance structure.
##' @param simplifies [logical] When there is only one variance-covariance matrix, output a matrix instead of a list of matrices.
##' @param ... Not used. For compatibility with the generic method.
##'
##'
##' @return A list where each element contains a residual variance-covariance matrix.
##'
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit mixed model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' ## extract residuals variance covariance matrix
##' getVarCov(eUN.lmm)
##' getVarCov(eUN.lmm, individual = c("1","5"))

## * getVarCov.lmm
getVarCov.lmm <- function(object, individual = NULL, p = NULL, type.object = c("lmm","gls"), simplifies = TRUE, strata = NULL, ...){

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.null(p) && any(names(which(object$param$type %in% c("sigma","k","rho"))) %in% names(p) == FALSE)){
        stop("Incorrect argument \'p\' - it should be a vector with names containing all variance and correlation parameters. \n")
    }
    
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    if(!is.null(individual)){
        individual <- match.arg(individual, object$design$cluster$levels, several.ok = TRUE)
    }

    if(type.object == "lmm"){

        if(!is.null(p)){
            Omega <- .calc_Omega(object = object$design$X.var,
                                 param = p,
                                 keep.interim = FALSE)
        }else{
            Omega <- object$Omega
        }
        
        if(is.null(individual)){
            Omega.strata <- object$design$X.var$strata
            Omega.time <- object$design$X.var$index.time
            
            n.timePattern <- unlist(lapply(Omega.time, length))
            index.fulltime <- which(n.timePattern==max(n.timePattern))
            if(!is.null(strata)){
                index.strata <- which(Omega.strata %in% strata)
                out <- stats::setNames(Omega,Omega.strata)[intersect(index.fulltime,index.strata)]
            }else{
                out <- stats::setNames(Omega,Omega.strata)[index.fulltime]
            }
        }else{
            out <- Omega[stats::setNames(object$design$X.var$cluster,object$design$cluster$levels)[individual]]
        }
        for(iO in 1:length(out)){
            dimnames(out[[iO]]) <- list(attr(out[[iO]],"time"),attr(out[[iO]],"time"))
            attr(out[[iO]],"time") <- NULL
            attr(out[[iO]],"sd") <- NULL
            attr(out[[iO]],"cor") <- NULL
        }
        if(is.list(out) && length(out)==1 && simplifies){
            return(out[[1]])
        }else{
            return(out)
        }
        
    }else if(type.object == "gls"){
        if(object$strata$n==1){
            if(is.null(individual)){
                return(getVarCov(object$gls[[1]]))
            }else{
                return(getVarCov(object$gls[[1]], individual = individual))
            }
        }else{
            if(is.null(individual)){
                return(lapply(object$gls, getVarCov))
            }else{
                out <- stats::setNames(vector(mode = "list", length = length(individual)),individual)
                for(iStrata in 1:object$strata$n){ ## iStrata <- 1
                    iIndiv <- intersect(individual,names(object$design$index.strata[object$design$index.strata==iStrata]))
                    if(length(iIndiv)>0){
                        out[match(iIndiv,names(out))] <- getVarCov(object$gls[[iStrata]], individual = iIndiv)
                    }
                }
                return(out)
            }
        }

    }
        
}

##----------------------------------------------------------------------
### getVarCov.R ends here
