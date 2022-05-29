### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (17:26) 
## Version: 
## Last-Updated: May 29 2022 (23:37) 
##           By: Brice Ozenne
##     Update #: 287
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
##' @name logLik
##' 
##' @param object a \code{lmm} object.
##' @param data [data.frame] dataset relative to which the log-likelihood should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the log-likelihood be output? Otherwise output the sum of all clusters of the derivatives. 
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
##' @return A numeric value (total logLikelihood) or a vector of numeric values, one for each cluster (cluster specific logLikelihood).
##' 

## * logLik
##' @rdname logLik
##' @export
logLik.lmm <- function(object, data = NULL, p = NULL, indiv = FALSE, ...){

    ## ** extract or recompute log-likelihood
    if(is.null(data) && is.null(p) && indiv == FALSE){
        design <- object$design ## useful in case of NA
        out <- object$logLik
    }else{
        test.precompute <- !is.null(object$design$precompute.XX) && !indiv
            
        if(!is.null(data)){
            design <- stats::model.matrix(object, data = data, effects = "all", simplifies = FALSE)
        }else{
            design <- object$design
        }
          
        if(!is.null(p)){
            if(any(duplicated(names(p)))){
                stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
            }
            if(any(names(object$param) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(object$param)[names(object$param) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            p <- p[names(object$param)]
        }else{
            p <- object$param
        }
        out <- .moments.lmm(value = p, design = design, time = object$time, method.fit = object$method.fit,
                            transform.sigma = "none", transform.k = "none", transform.rho = "none",
                            logLik = TRUE, score = FALSE, information = FALSE, vcov = FALSE, df = FALSE, indiv = indiv, 
                            trace = FALSE, precompute.moments = test.precompute)$logLik
    } 

    ## ** restaure NAs and name
    if(indiv){
        if(is.null(data) && length(object$index.na)>0 && any(is.na(attr(object$index.na,"cluster.index")))){
            names(out) <- object$design$cluster$levels
            out.save <- out
            out <- stats::setNames(rep(NA, times = object$cluster$n), object$cluster$levels)
            out[rownames(out.save)] <- out.save

            if(is.numeric(design$cluster$levels.original)){
                names(out) <- NULL
            }
        }else if(!is.numeric(design$cluster$levels.original)){
            names(out) <- design$cluster$levels.original
        }
    }

    ## ** export
    return(out)
}

## * .logLik
.logLik <- function(X, residuals, precision, Upattern.ncluster, weights, scale.Omega,
                    index.variance, time.variance, index.cluster,
                    indiv, REML, precompute){

    ## ** extract information
    if(indiv && REML){##  https://towardsdatascience.com/maximum-likelihood-ml-vs-reml-78cf79bef2cf
        stop("Cannot compute individual likelihood contribution with REML. \n")
    }
    test.loopIndiv <- indiv || is.null(precompute)

    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    n.mucoef <- NCOL(X)
    name.mucoef <- colnames(X)
    log2pi <- log(2*pi)
    REML.det <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))

    ## ** prepare output
    if(test.loopIndiv){
        ll <- rep(NA, n.cluster)
    }else{
        ll <- 0
    }

    ## ** compute log-likelihood
    ## *** looping over individuals
    if(test.loopIndiv){
        ## precompute
        logidet.precision <- lapply(precision, function(iM) {-log(base::det(iM))})
        
        ## loop
        for (iId in 1:n.cluster) { ## iId <- 1
            iIndex <- index.cluster[[iId]]
            iResidual <- residuals[iIndex, , drop = FALSE]
            iX <- X[iIndex, , drop = FALSE]
            iOmegaM1 <- precision[[index.variance[iId]]] * scale.Omega[iId]
            iWeight <- weights[iId]

            ll[iId] <- - iWeight * (NCOL(iOmegaM1) * (log2pi-log(scale.Omega[iId])) + logidet.precision[[index.variance[iId]]] + t(iResidual) %*% iOmegaM1 %*% iResidual)/2
            if (REML) {
                REML.det <- REML.det + iWeight * (t(iX) %*% iOmegaM1 %*% iX)
            }
            ## log(det(iOmegaM1)) - NCOL(iOmegaM1)*log(scale.Omega[iId])+logidet.precision[[index.variance[iId]]]
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
            iLogDet.Omega <- log(base::det(iOmegaM1))
            ll <- ll - 0.5 * unname(Upattern.ncluster[iPattern]) * (NCOL(iOmegaM1) * log2pi - iLogDet.Omega) - 0.5 * sum(precompute$RR[[iPattern]] * iOmegaM1)
            if (REML) {
                ## compute (unique contribution, i.e. only lower part of the matrix)
                ## iContribution <- apply(precompute$XX$pattern[[iPattern]], MARGIN = 3, FUN = function(iM){sum(iM * iOmegaM1)})
                iContribution <- as.double(iOmegaM1) %*% matrix(precompute$XX$pattern[[iPattern]], nrow = length(iOmegaM1), ncol = dim(precompute$XX$pattern[[iPattern]])[3], byrow = FALSE)
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
