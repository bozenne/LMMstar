### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (17:26) 
## Version: 
## Last-Updated: Jun  4 2021 (09:37) 
##           By: Brice Ozenne
##     Update #: 143
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * logLik.lmm (documentation)
##' @title Extract The Log-Likelihood From a Multivariate Gaussian Model
##' @description Extract or compute the log-likelihood of a multivariate gaussian model.
##' @name logLik
##' 
##' @param object a \code{lmm} object.
##' @param data [data.frame] dataset relative to which the log-likelihood should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the log-likelihood be output? Otherwise output the sum of all clusters of the derivatives. 
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related method.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the log-likelihood. Only relevant if differs from the fitted values.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \bold{transform}: \cr
##' \itemize{
##' \item 0 means no transformation i.e. ouput stanrdard error, ratio of standard errors, and correlations.
##' \item 1 means log/atanh transformation i.e. ouput log(stanrdard error), log(ratio of standard errors), and atanh(correlations).
##' \item 2 ouput variance coefficients and correlations.
##' }
##'
##' @details \bold{indiv}: only relevant when using maximum likelihood. Must be \code{FALSE} when using restricted maximum likelihood.
##' 
##' @return A numeric value
##' 

## * logLik
##' @rdname logLik
##' @export
logLik.lmm <- function(object, data = NULL, p = NULL, type.object = "lmm", indiv = FALSE, ...){

    ## ** normalize user input
    type.object <- match.arg(type.object, c("lmm","gls"))

    ## ** extract or recompute log-likelihood
    if(type.object=="lmm"){
        if(is.null(data) && is.null(p) && indiv == FALSE){
            out <- object$logLik
        }else{
            if(!is.null(data)){
                ff.allvars <- c(all.vars(object$formula$mean), all.vars(object$formula$var))
                if(any(ff.allvars %in% names(data) == FALSE)){
                    stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
                }

                design <- .model.matrix.lmm(formula.mean = object$formula$mean.design,
                                            formula.var = object$formula$var.design,
                                            data = data,
                                            var.outcome = object$outcome$var,
                                            var.strata = object$strata$var, U.strata = object$strata$levels,
                                            var.time = object$time$var, U.time = object$time$levels,
                                            var.cluster = object$cluster$var,
                                            structure = object$structure
                                            )
                Y <- design$Y
                X <- design$X.mean
                index.vargroup <- design$X.var$cluster
                index.cluster <- design$index.cluster
                index.time <- design$index.time
                X.var <- design$X.var
            }else{
                Y <- object$design$Y
                X <- object$design$X.mean
                index.vargroup <- object$design$X.var$cluster
                index.cluster <- object$design$index.cluster
                index.time <- object$design$index.time
                X.var <- object$design$X.var
            }
            if(!is.null(p)){
                if(any(duplicated(names(p)))){
                    stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
                }
                if(any(names(object$param$type) %in% names(p) == FALSE)){
                    stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(object$param$type)[names(object$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
                }
                beta <- p[object$param$type=="mu"]
                Omega <- .calc_Omega(object = X.var, param = p)
                precision <- lapply(Omega, solve)
            }else{
                beta <- object$param$value[object$param$type=="mu"]
                precision <- object$OmegaM1
            }

            out <- .logLik(X = X, residuals = Y - X %*% beta, precision = precision,
                           index.variance = index.vargroup, time.variance = index.time, index.cluster = index.cluster, 
                           indiv = indiv, REML = object$method.fit=="REML")
        
        } ## end if data, p
    }else if(type.object=="gls"){
        if(!is.null(data)){
            stop("Cannot handle argument \'data\' when argument \'type.object\' is \"gls\". \n")
        }
        if(!is.null(p)){
            stop("Cannot handle argument \'p\' when argument \'type.object\' is \"gls\". \n")
        }

        if(is.null(object$variable$strata)){
            out <- stats::logLik(object$gls[[1]])
        }else{
            out <- lapply(object$gls, stats::logLik)
        }
    }

    ## ** export
    return(out)
}

## * .logLik
.logLik <- function(X, residuals, precision,
                    index.variance, time.variance, index.cluster,
                    indiv, REML){

    if(indiv && REML){
        #  https://towardsdatascience.com/maximum-likelihood-ml-vs-reml-78cf79bef2cf
        stop("Cannot compute individual likelihood contribution with REML. \n")
    }
    
    ## ** prepare
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    n.mucoef <- NCOL(X)
    name.mucoef <- colnames(X)
    ll <- rep(NA, n.cluster)

    logidet.precision <- lapply(precision, function(iM){-log(det(iM))}) ## log(det(\Omega)) = - log(det(\Omega^-1))
    log2pi <- log(2*pi)
    REML.det <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
    
    ## ** compute score
    for(iId in 1:n.cluster){ ## iId <- 7
        iIndex <- attr(index.cluster,"sorted")[[iId]]
        ## iIndex <- which(index.cluster==iId)
        ## iIndex <- iIndex[order(time.variance[iIndex])] ## re-order observations according to the variance-covariance matrix

        iResidual <- residuals[iIndex,,drop=FALSE]
        iX <- X[iIndex,,drop=FALSE]
        iOmega <- precision[[index.variance[iId]]]
        ll[iId] <- - (NCOL(iOmega) * log2pi + logidet.precision[[index.variance[iId]]] + t(iResidual) %*% iOmega %*% iResidual)/2
        if(REML){
            REML.det <- REML.det + t(iX) %*% iOmega %*% iX
        }
    }

    ## ** export
    if(indiv){
        return(ll)
    }else{
        ll <- sum(ll)
        if(REML){
            ll <- ll - log(det(REML.det))/2 + n.mucoef * log2pi/2
        }
        return(ll)
    }
}


##----------------------------------------------------------------------
### logLik.R ends here
