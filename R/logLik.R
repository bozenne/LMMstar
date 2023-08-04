### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (17:26) 
## Version: 
## Last-Updated: aug  3 2023 (16:16) 
##           By: Brice Ozenne
##     Update #: 390
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
##' @param data [data.frame] dataset relative to which the log-likelihood should be computed. Only relevant if differs from the dataset used to fit the model.
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
logLik.lmm <- function(object, data = NULL, p = NULL, indiv = FALSE, ...){

    ## ** extract or recompute log-likelihood
    if(is.null(data) && is.null(p) && indiv == FALSE){
        design <- object$design ## useful in case of NA
        out <- object$logLik
    }else{
        test.precompute <- !is.null(object$design$precompute.XX) && !indiv

        if(!is.null(data)){
            design <- stats::model.matrix(object, data = data, effects = "all", simplify = FALSE)
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
        out <- .moments.lmm(value = p, design = design, time = object$time, method.fit = object$args$method.fit,
                            transform.sigma = "none", transform.k = "none", transform.rho = "none",
                            logLik = TRUE, score = FALSE, information = FALSE, vcov = FALSE, df = FALSE, indiv = indiv, 
                            trace = FALSE, precompute.moments = test.precompute)$logLik
    } 

    ## ** name and restaure NAs
    if(indiv){

        if(!is.numeric(object$cluster$levels)){
            names(out) <- object$cluster$levels[match(1:length(out),object$cluster$index)]
        }
        out <- addNA(out, index.na = object$index.na,
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
.logLik <- function(X, residuals, precompute,
                    Upattern.ncluster, Upattern.ntime, weights, scale.Omega, 
                    pattern, index.cluster, indiv, REML){
    
    ## ** extract information
    if(indiv && REML){
        stop("Cannot compute cluster-specific likelihood contributions when using REML. \n")
    }
    test.loopIndiv <- indiv || !attr(precompute,"moments")
    log2pi <- log(2*pi)

    if(!indiv && any(is.na(precompute$Omega.logdet))){ ## non positive definite residual variance covariance
        return(NA)
    }

    ## ** compute log-likelihood
    ## *** looping over individuals
    if(test.loopIndiv){
        n.cluster <- length(pattern) ## number of clusters, may different from Upattern.ncluster which is the weight of each cluster
        ll <- rep(NA, n.cluster)
        
        for (iId in 1:n.cluster) { ## iId <- 1
            iResidual <- residuals[index.cluster[[iId]], , drop = FALSE]
            iOmegaM1 <- precompute$precision[[pattern[iId]]] * scale.Omega[iId]
            ll[iId] <- - weights[iId] * (NCOL(iOmegaM1) * (log2pi-log(scale.Omega[iId])) + precompute$Omega.logdet[[pattern[iId]]] + t(iResidual) %*% iOmegaM1 %*% iResidual)/2            
        }
        if(!indiv){
            ll <- sum(ll)
        }
    }

    ## *** looping over covariance patterns
    if(!test.loopIndiv){
        term1 <- Upattern.ncluster * (Upattern.ntime * log2pi + precompute$Omega.logdet)/2
        term2 <- sum(do.call(rbind,lapply(precompute$wRR,attr,"vectorwise")) * do.call(rbind,lapply(precompute$OmegaM1,attr,"vectorwise")))
        ll <- - term1 - term2

        ## microbenchmark(a = sum(do.call(rbind,lapply(precompute$wRR,attr,"vectorwise")) *do.call(rbind,lapply(precompute$OmegaM1,attr,"vectorwise"))),
        ##                b = mapply(x = precompute$wRR, y = precompute$OmegaM1, function(x,y){sum(attr(x,"vectorwise")*attr(y,"vectorwise"))}),
        ##                c = mapply(x = precompute$wRR, y = precompute$OmegaM1, function(x,y){sum(x*y)}))
    }

    ## ** export
    if(REML){
        ll <- ll - precompute$wXOmegaM1X.logdet/2 + NCOL(X) * log2pi/2
    }
    return(unname(ll))
}


##----------------------------------------------------------------------
### logLik.R ends here
