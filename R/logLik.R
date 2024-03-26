### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (17:26) 
## Version: 
## Last-Updated: Mar 26 2024 (09:40) 
##           By: Brice Ozenne
##     Update #: 357
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * logLik.lmm (documentation)
##' @title Extract The Log-Likelihood From a Linear Mixed Model
##' @description Extract or compute the log-likelihood of a linear mixed model.
##' 
##' @param object a \code{lmm} object.
##' @param newdata [data.frame] dataset relative to which the log-likelihood should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the log-likelihood be output? Otherwise output the sum of all clusters of the derivatives. 
##' @param p [numeric vector] value of the model coefficients at which to evaluate the log-likelihood. Only relevant if differs from the fitted values.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \bold{indiv}: only relevant when using maximum likelihood. Must be \code{FALSE} when using restricted maximum likelihood.
##' 
##' @return A numeric value (total logLikelihood) or a vector of numeric values, one for each cluster (cluster specific logLikelihood).
##' 
##' @keywords methods

## * logLik.lmm (code)
##' @export
logLik.lmm <- function(object, newdata = NULL, p = NULL, indiv = FALSE, ...){

    ## ** extract or recompute log-likelihood
    if(is.null(newdata) && is.null(p) && indiv == FALSE){
        design <- object$design ## useful in case of NA
        out <- object$logLik
    }else{
        test.precompute <- !is.null(object$design$precompute.XX) && !indiv

        ## normalize newdata argument
        if(!is.null(newdata)){
            design <- stats::model.matrix(object, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- object$design
        }

        ## normalize p argument
        if(!is.null(p)){
            init <- .init_transform(p = p, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, 
                                    x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                                    table.param = object$design$param)
            theta <- init$p
        }else{
            theta <- object$param
        }

        ## evaluate log-lik
        out <- .moments.lmm(value = theta, design = design, time = object$time, method.fit = object$args$method.fit,
                            transform.sigma = "none", transform.k = "none", transform.rho = "none",
                            logLik = TRUE, score = FALSE, information = FALSE, vcov = FALSE, df = FALSE, indiv = indiv, 
                            trace = FALSE, precompute.moments = test.precompute)$logLik
    } 

    ## ** name and restaure NAs
    if(indiv){

        if(!is.numeric(object$cluster$levels)){
            names(out) <- object$cluster$levels[match(1:length(out),object$cluster$index)]
        }
        out <- restaureNA(out, index.na = object$index.na,
                          level = "cluster", cluster = object$cluster)        

    }

    ## ** export
    return(out)
}

## * logLik.mlmm (code)
##' @export
logLik.mlmm <- function(object, ...){

    return(lapply(object$model, logLik, ...))

}

## * .logLik
.logLik <- function(X, residuals, precision, Upattern.ncluster, weights, scale.Omega,
                    pattern, index.cluster, indiv, REML, precompute){

    ## ** extract information
    if(indiv && REML){##  https://towardsdatascience.com/maximum-likelihood-ml-vs-reml-78cf79bef2cf
        stop("Cannot compute individual likelihood contribution with REML. \n")
    }
    test.loopIndiv <- indiv || is.null(precompute)

    n.obs <- length(index.cluster)
    n.cluster <- length(pattern) ## number of clusters, may different from Upattern.ncluster which is the weight of each cluster
    n.mucoef <- NCOL(X)
    name.mucoef <- colnames(X)
    log2pi <- log(2*pi)
    if(REML){
        REML.det <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
    }
    logdet.precision <- attr(precision, "logdet")

    ## ** prepare output
    if(test.loopIndiv){
        ll <- rep(NA, n.cluster)
    }else{
        ll <- 0
    }
    if(any(is.na(logdet.precision))){ ## non positive definite residual variance covariance
        if(indiv){
            return(ll*NA)
        }else{
            return(NA)
        }
    }

    ## ** compute log-likelihood
    ## *** looping over individuals
    if(test.loopIndiv){
        ## loop
        for (iId in 1:n.cluster) { ## iId <- 1
            iIndex <- index.cluster[[iId]]
            iResidual <- residuals[iIndex, , drop = FALSE]
            iX <- X[iIndex, , drop = FALSE]
            iOmegaM1 <- precision[[pattern[iId]]] * scale.Omega[iId]
            iWeight <- weights[iId]
            ll[iId] <- - iWeight * (NCOL(iOmegaM1) * (log2pi-log(scale.Omega[iId])) - logdet.precision[[pattern[iId]]] + t(iResidual) %*% iOmegaM1 %*% iResidual)/2
            if (REML) {
                REML.det <- REML.det + iWeight * (t(iX) %*% iOmegaM1 %*% iX)
            }
            ## log(det(iOmegaM1)) - NCOL(iOmegaM1)*log(scale.Omega[iId])+logdet.precision[[pattern[iId]]]
        }
        if(!indiv){
            ll <- sum(ll)
        }
    }
    
    ## *** looping over covariance patterns
    if(!test.loopIndiv){
        ## precompute
        n.pattern <- length(Upattern.ncluster)

        ## loop
        for (iPattern in 1:n.pattern) { ## iPattern <- 1
            iOmegaM1 <- precision[[iPattern]]
            ll <- ll - 0.5 * unname(Upattern.ncluster[iPattern]) * (NCOL(iOmegaM1) * log2pi - logdet.precision[[iPattern]]) - 0.5 * sum(precompute$RR[[iPattern]] * iOmegaM1)
            if (REML) {
                ## compute (unique contribution, i.e. only lower part of the matrix)
                if(is.null(precompute$X.OmegaM1.X)){
                    iContribution <- as.double(iOmegaM1) %*% precompute$XX$pattern[[iPattern]]
                }else{
                    iContribution <- precompute$X.OmegaM1.X[[iPattern]]
                }
                ## fill the matrix
                REML.det <- REML.det + iContribution[as.vector(precompute$XX$key)]
            }
        }
    }

    ## ** export
    if(REML){
        ll <- ll - log(det(REML.det))/2 + n.mucoef * log2pi/2
    }
    return(ll)
}


##----------------------------------------------------------------------
### logLik.R ends here
