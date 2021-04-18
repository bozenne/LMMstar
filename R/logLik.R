### logLik.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (17:26) 
## Version: 
## Last-Updated: Apr 16 2021 (16:40) 
##           By: Brice Ozenne
##     Update #: 93
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

    ## ** extract or recompute log-likelihood
    if(type=="lmm"){
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
                index.vargroup <- design$index.vargroup
                index.cluster <- design$index.cluster
                X.var <- design$X.var
            }else{
                Y <- object$design$Y
                X <- object$design$X.mean
                index.vargroup <- object$design$index.vargroup
                index.cluster <- object$design$index.cluster
                X.var <- object$design$X.var
            }
            if(!is.null(p)){
                if(any(names(object$param$type) %in% names(p) == FALSE)){
                    stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(object$param$type)[names(object$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
                }
                beta <- p[names(object$param$mu)]
                Omega <- attr(X.var,"FUN.Omega")(object = X.var, sigma = p[names(object$param$sigma)], k = p[names(object$param$k)], rho = p[names(object$param$cor)])
                precision <- lapply(Omega, solve)
            }else{
                beta <- object$param$mu
                precision <- object$OmegaM1
            }
            out <- .logLik(X = X, residuals = Y - X %*% beta, precision = precision,
                           index.variance = index.vargroup, index.cluster = index.cluster, 
                           indiv = indiv, REML = object$method.fit=="REML")
        
        } ## end if data, p
    }else if(type=="gls"){
        if(!is.null(data)){
            stop("Cannot handle argument \'data\' when argument \'type\' is \"gls\". \n")
        }
        if(!is.null(p)){
            stop("Cannot handle argument \'p\' when argument \'type\' is \"gls\". \n")
        }

        if(is.null(object$variable$strata)){
            out <- logLik(object$gls[[1]])
        }else{
            out <- lapply(object$gls, logLik)
        }
    }

    ## ** export
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

    logidet.precision <- lapply(precision, function(iM){-log(det(iM))}) ## log(det(\Omega)) = - log(det(\Omega^-1))
    log2pi <- log(2*pi)
    REML.det <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
    
    ## ** compute score
    for(iId in 1:n.cluster){ ## iId <- 7
        iResidual <- residuals[index.cluster==iId,,drop=FALSE]
        iX <- X[index.cluster==iId,,drop=FALSE]
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
