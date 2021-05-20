### getVarCov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: May 20 2021 (09:44) 
##           By: Brice Ozenne
##     Update #: 52
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

## * getVarCov.lmm
getVarCov.lmm <- function(object, individual = NULL, p = NULL, type.object = c("lmm","gls"), simplifies = TRUE, strata = NULL, ...){

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.null(p) && any(names(c(object$param$sigma,object$param$k,object$param$rho)) %in% names(p) == FALSE)){
        stop("Incorrect argument \'p\' - it should be a vector with names containing all variance and correlation parameters. \n")
    }
    
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }

    if(!is.null(individual)){
        individual <- match.arg(individual, object$cluster$levels, several.ok = TRUE)
    }

    if(type.object == "lmm"){

        if(!is.null(p)){
            Omega <- .calc_Omega(object = object$design$X.var,
                                 sigma = p[names(which(object$param$type=="sigma"))],
                                 k = p[names(which(object$param$type=="k"))],
                                 rho = p[names(which(object$param$type=="rho"))],
                                 keep.interim = FALSE)
        }else{
            Omega <- object$Omega
        }
        
        if(is.null(individual)){
            n.timePattern <- unlist(lapply(attr(object$design$X.var,"Upattern.time"),length))
            index.fulltime <- which(n.timePattern==max(n.timePattern))
            if(!is.null(strata)){
                index.strata <- which(attr(object$design$X.var,"UX.strata") %in% strata)
            }else{
                index.strata <- 1:attr(object$design$X.var,"nUpattern")
                strata <- attr(object$design$X.var,"UX.strata")
            }
            out <- setNames(Omega[intersect(index.fulltime,index.strata)],strata)
        }else{
            out <- Omega[object$design$index.vargroup[individual]]
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
                out <- setNames(vector(mode = "list", length = length(individual)),individual)
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
