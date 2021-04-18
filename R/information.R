### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: Apr 16 2021 (17:58) 
##           By: Brice Ozenne
##     Update #: 154
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * information.lmm (code)
##' @export
information.lmm <- function(x, data = NULL, p = NULL, transform = NULL, indiv = FALSE, type = NULL, ...){
    options <- LMMstar.options()
    x.transform <- attr(e.lmm$design$X.var, "transform")

    ## ** normalize user input
    if(is.null(transform)){transform <- options$transform}
    if(is.null(type)){
        type <- options$type.information
    }else{
        type <- match.arg(type, c("expected","observed"))
    }

    ## ** extract or recompute information
    if(is.null(data) && is.null(p) && (indiv == FALSE) && (transform == x.transform)){
        out <- x$information
    }else{
        REML <- x$method.fit == "REML"

        if(!is.null(data) || (transform != x.transform)){
            ff.allvars <- c(all.vars(x$formula$mean), all.vars(x$formula$var))
            if(is.null(data)){
                data <- x$data
            }else if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            design <- .model.matrix.lmm(formula.mean = x$formula$mean.design,
                                        formula.var = x$formula$var.design,
                                        data = data,
                                        var.outcome = x$outcome$var,
                                        var.strata = x$strata$var, U.strata = x$strata$levels,
                                        var.time = x$time$var, U.time = x$time$levels,
                                        var.cluster = x$cluster$var,
                                        structure = x$structure,
                                        transform = transform
                                        )
            Y <- design$Y
            X <- design$X.mean
            index.vargroup <- design$index.vargroup
            index.cluster <- design$index.cluster
            X.var <- design$X.var
            pair.varcoef  <- design$param$pair.varcoef
        }else{
            Y <- x$design$Y
            X <- x$design$X.mean
            index.vargroup <- x$design$index.vargroup
            index.cluster <- x$design$index.cluster
            X.var <- x$design$X.var
            pair.varcoef  <- x$design$param$pair.varcoef
        }
        if(!is.null(p)){
            if(any(names(x$param$type) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param$type)[names(x$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            beta <- p[names(x$param$mu)]
            Omega <- attr(X.var,"FUN.Omega")(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$cor)], keep.interim = TRUE)
            dOmega <- attr(X.var,"FUN.dOmega")(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$cor)], Omega = Omega)
            if(REML || type == "observed"){
                d2Omega <- attr(X.var,"FUN.d2Omega")(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$cor)], Omega = Omega, dOmega = dOmega, pair = pair.varcoef)
            }
            precision <- lapply(Omega, solve)
        }else{
            beta <- x$param$mu
            precision <- x$OmegaM1
            dOmega <- x$dOmega
            d2Omega <- x$d2Omega
        }
        out <- .information(X = X, residuals = Y - X %*% beta, precision = precision, dOmega = dOmega, d2Omega = d2Omega,
                            index.variance = index.vargroup, index.cluster = index.cluster,
                            pair.varcoef = pair.varcoef, indiv = indiv, REML = REML, type = type)
        
    }
    return(out)
}

## * .information
## REML term
## d 0.5 tr[(X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 X)] = 0.5 tr[ (X \OmegaM1 d'\Omega \OmegaM1 X) (X \OmegaM1 X)^{-2} (X \OmegaM1 d\Omega \OmegaM1 X) ]
##                                                                 - 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d'\Omega \OmegaM1 d\Omega \OmegaM1 X) + (X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 d'\Omega \OmegaM1 X) ]
##                                                                 + 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d2\Omega \OmegaM1 X) ]
.information <- function(X, residuals, precision, dOmega, d2Omega,
                         index.variance, index.cluster,
                         pair.varcoef, indiv, REML, type){

    if(indiv && REML){
        stop("Not possible to compute individual information when using REML.\n")
    }

    ## ** prepare
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- names(dOmega[[1]])
    n.varcoef <- length(name.varcoef)
    name.allcoef <- c(name.mucoef,name.varcoef)
    n.allcoef <- length(name.allcoef)
    var.pattern <- names(dOmega)

    if(indiv){
        info <- array(NA, dim = c(n.cluster, n.allcoef, n.allcoef),
                     dimnames = list(NULL, name.allcoef, name.allcoef))
    }else{
        info <- matrix(0, nrow = n.allcoef, ncol = n.allcoef,
                      dimnames = list(name.allcoef, name.allcoef)
                      )
    }    

    ## precompute derivative for the variance
    npair.varcoef <- NCOL(pair.varcoef)

    dOmega.precomputed <- setNames(lapply(var.pattern, function(iPattern){
        iOut <- vector(mode = "list", length = npair.varcoef)

        iResiduals <- residuals[index.variance==iPattern,,drop=FALSE]
        obsVSexp <- diag(1, NCOL(iResiduals), NCOL(iResiduals)) -  precision[[iPattern]] %*% crossprod(iResiduals)/NROW(iResiduals)

        for(iPair in 1:npair.varcoef){
            iOut[[iPair]] <- tr(precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[1,iPair]]] %*% precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[2,iPair]]])
            if(type == "observed"){
                iOut[[iPair]] <- iOut[[iPair]] - tr(precision[[iPattern]] %*% d2Omega[[iPattern]][[iPair]] %*% obsVSexp)
            }
        }
        return(iOut)
    }), var.pattern)

    if(REML){
        dOmega.precomputed2 <- setNames(lapply(var.pattern, function(iPattern){
            iOut <- list(term1a = vector(mode = "list", length = npair.varcoef),
                         term1b = vector(mode = "list", length = npair.varcoef),
                         term2 = vector(mode = "list", length = npair.varcoef),
                         term3 = vector(mode = "list", length = npair.varcoef))                         

            for(iPair in 1:npair.varcoef){
                iOut$term1a[[iPair]] <- precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[1,iPair]]] %*% precision[[iPattern]]
                iOut$term1b[[iPair]] <- precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[2,iPair]]] %*% precision[[iPattern]]
                iOut$term2[[iPair]] <- precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[1,iPair]]] %*% precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[2,iPair]]] %*% precision[[iPattern]]
                iOut$term3[[iPair]] <- precision[[iPattern]] %*% d2Omega[[iPattern]][[iPair]] %*% precision[[iPattern]]
            }
            return(iOut)
        }), var.pattern)

        REML.num <- lapply(1:npair.varcoef, function(i){
            list(numerator1a = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)),
                 numerator1b = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)),
                 numerator2 = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)),
                 numerator3 = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)))
        })
        REML.denom <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
    }


    ## ** compute information
    for(iId in 1:n.cluster){ ## iId <- 7
        iPattern <- index.variance[iId]
        iIndex <- which(index.cluster==iId)

        iX <- X[iIndex,,drop=FALSE]
        iOmegaM1 <- precision[[iPattern]]
                              
        if(indiv){
            info[iId,name.mucoef,name.mucoef] <- t(iX) %*% iOmegaM1 %*% iX
        }else{
            info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + t(iX) %*% iOmegaM1 %*% iX
        }

        if(REML){
            REML.denom <- REML.denom + t(iX) %*% iOmegaM1 %*% iX
        }

        for(iPair in 1:npair.varcoef){ ## iPair <- 1
            if(indiv){
                info[iId,pair.varcoef[1,iPair],pair.varcoef[2,iPair]] <- 0.5 * dOmega.precomputed[[iPattern]][[iPair]]
                if(pair.varcoef[1,iPair] != pair.varcoef[2,iPair]){
                    info[iId,pair.varcoef[2,iPair],pair.varcoef[1,iPair]] <- 0.5 * dOmega.precomputed[[iPattern]][[iPair]]
                }
            }else{
                info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] <- info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] + 0.5 * dOmega.precomputed[[iPattern]][[iPair]]
                if(pair.varcoef[1,iPair] != pair.varcoef[2,iPair]){
                    info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] <- info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] + 0.5 * dOmega.precomputed[[iPattern]][[iPair]]
                }
                if(REML){
                    REML.num[[iPair]]$numerator1a <- REML.num[[iPair]]$numerator1a + t(iX) %*% (dOmega.precomputed2[[iPattern]]$term1a[[iPair]]) %*% iX
                    REML.num[[iPair]]$numerator1b <- REML.num[[iPair]]$numerator1b + t(iX) %*% (dOmega.precomputed2[[iPattern]]$term1b[[iPair]]) %*% iX
                    REML.num[[iPair]]$numerator2 <- REML.num[[iPair]]$numerator2 + t(iX) %*% (dOmega.precomputed2[[iPattern]]$term2[[iPair]]) %*% iX
                    REML.num[[iPair]]$numerator3 <- REML.num[[iPair]]$numerator3 + t(iX) %*% (dOmega.precomputed2[[iPattern]]$term3[[iPair]]) %*% iX
                }
            }
        }
    }

    ## ** export
    if(REML){
        REML.denom <- solve(REML.denom)
        for(iPair in 1:npair.varcoef){ ## iPair <- 1
            term1 <- tr(REML.denom %*% REML.num[[iPair]]$numerator1a %*% REML.denom %*% REML.num[[iPair]]$numerator1b)
            term2 <- tr(- 2 * REML.denom %*% REML.num[[iPair]]$numerator2)
            term3 <- tr(REML.denom %*% REML.num[[iPair]]$numerator3)
            print(c(info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]],term1,term2,term3))
            info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] <- info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] + 0.5 * (term1 + term2 + term3) 
            if(pair.varcoef[1,iPair] != pair.varcoef[2,iPair]){
                info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] <- info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] + 0.5 * (term1 + term2 + term3) 
            }
        }
    }
    return(info)

}


##----------------------------------------------------------------------
### information.R ends here
