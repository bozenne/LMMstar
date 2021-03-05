### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (17:26) 
## Version: 
## Last-Updated: mar  5 2021 (23:07) 
##           By: Brice Ozenne
##     Update #: 34
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * logLik
logLik.lmm <- function(object, type = "lmm"){
    type <- match.arg(type, c("lmm","gls"))

    if(type=="lmm"){
        return(object$logLik)
    }else if(type=="gls"){
        if(is.null(object$variable$strata)){
            return(logLik(object$gls[[1]]))
        }else{
            return(lapply(object$gls, logLik))
        }
    }
}

## * .logLik
.logLik <- function(Y, X, beta, sigma, k, rho, precision,
                    index.variance, index.cluster, indiv, REML){

    if(indiv && REML){
        #  https://towardsdatascience.com/maximum-likelihood-ml-vs-reml-78cf79bef2cf
        stop("Cannot compute individual likelihood contribution with REML. \n")
    }
    
    ## ** prepare
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    n.allcoef <- NCOL(X)
    name.allcoef <- colnames(X)
    ll <- rep(NA, n.cluster)

    residuals <- Y - X %*% beta
    logidet.precision <- -log(sapply(precision, det))
    log2pi <- log(2*pi)
    REML.det <- 0
    
    ## ** compute score
    for(iId in 1:n.cluster){ ## iId <- 7
        iResidual <- residuals[index.cluster==iId,,drop=FALSE]
        iX <- X[index.cluster==iId,,drop=FALSE]
        iOmega <- precision[[index.variance[iId]]]        
        ll[iId] <- - (NCOL(iOmega) * log2pi + logidet.precision[index.variance[iId]] + t(iResidual) %*% iOmega %*% iResidual)/2
        if(REML){
            REML.det <- REML.det + t(iX) %*% iOmega %*% iX
        }
    }

    ## ** export
    if(indiv){
        return(ll)
    }else{
        if(REML){
            return(sum(ll) - log(det(REML.det))/2 + n.allcoef * log2pi/2)
        }else{
            return(sum(ll))
        }
    }
}


##----------------------------------------------------------------------
### logLik.R ends here
