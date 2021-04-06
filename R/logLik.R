### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (17:26) 
## Version: 
## Last-Updated: mar 23 2021 (11:51) 
##           By: Brice Ozenne
##     Update #: 59
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * logLik
logLik.lmm <- function(object, data = NULL, p = NULL, type = "lmm", indiv = FALSE, ...){

    ## ** normalize user input
    type <- match.arg(type, c("lmm","gls"))

    if(type=="lmm"){
        if(is.null(data) && is.null(p) && indiv == FALSE){
            out <- x$logLik
        }else{
            if(!is.null(data)){
                ff.allvars <- c(all.vars(x$formula$mean), all.vars(x$formula$var))
                if(any(ff.allvars %in% names(data) == FALSE)){
                    stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
                }

                design <- .model.matrix.lmm(formula.mean = x$formula$mean.design,
                                            formula.var = x$formula$var.design,
                                            data = data,
                                            var.outcome = x$outcome$var,
                                            var.strata = x$strata$var, U.strata = x$strata$levels,
                                            var.time = x$time$var, U.time = x$time$levels,
                                            var.cluster = x$cluster$var,
                                            structure = x$structure
                                            )
                X <- design$X.mean
                index.variance <- design$index.vargroup
                index.cluster <- design$index.cluster
                X.var <- design$X.var
            }else{
                X <- x$design$X.mean
                index.variance <- x$design$index.vargroup
                index.cluster <- x$design$index.cluster
                X.var <- x$design$X.var
            }
            if(!is.null(p)){
                if(any(names(x$param$type) %in% names(p) == FALSE)){
                    stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param$type)[names(x$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
                }
                beta <- p[names(x$param$mu)]
                Omega <- attr(X.var,"FUN.Omega")(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$cor)])
                precision <- lapply(Omega, solve)
            }else{
                beta <- x$param$mu
                precision <- x$OmegaM1
            }
            out <- .logLik(X = X, residuals = Y - X %*% beta, precision = precision,
                           index.variance = index.vargroup, index.cluster = index.cluster, 
                           indiv = indiv, REML = object$method.fit=="REML")
        
        } ## end if data, p
    }else if(type=="gls"){
        if(is.null(object$variable$strata)){
            out <- logLik(object$gls[[1]])
        }else{
            out <- lapply(object$gls, logLik)
        }
    }

    return(out)
}

## * .logLik
.logLik <- function(X, residuals, precision,
                    index.variance, index.cluster,
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

    logidet.precision <- -log(sapply(precision, det))
    log2pi <- log(2*pi)
    REML.det <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
    
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
        ll <- sum(ll)
        if(REML){
            ll <- ll - log(det(REML.det))/2 + n.mucoef * log2pi/2
        }
        return(ll)
    }
}


##----------------------------------------------------------------------
### logLik.R ends here
