### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: mar 26 2021 (09:37) 
##           By: Brice Ozenne
##     Update #: 57
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
information.lmm <- function(x, data = NULL, p = NULL, transform = TRUE, indiv = FALSE, ...){

    if(is.null(data) && is.null(p) && indiv == FALSE && transform == TRUE){
        out <- x$information
    }else{
        if(!is.null(data) || transform == FALSE){
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
                                        structure = x$structure,
                                        transform = transform
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
            Omega <- attr(X.var,"FUN.Omega")(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$cor)])
            precision <- lapply(Omega, solve)
        }else{
            precision <- x$OmegaM1
        }
        out <- .information(X = X, precision = precision,
                            index.variance = index.variance, index.cluster = index.cluster,
                            indiv = indiv, REML = object$method.fit == "REML")
        
    }
    return(out)
}

## * .information
## REML term
## d 0.5 tr[(X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 X)] = 0.5 tr[ (X \OmegaM1 d'\Omega \OmegaM1 X) (X \OmegaM1 X)^{-2} (X \OmegaM1 d\Omega \OmegaM1 X) ]
##                                                                 - 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d'\Omega \OmegaM1 d\Omega \OmegaM1 X) + (X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 d'\Omega \OmegaM1 X) ]
##                                                                 + 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d2\Omega \OmegaM1 X) ]
.information <- function(X, precision, dOmega,
                         index.variance, index.cluster,
                         indiv, REML){

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
    pair.varcoef <- .unorderedPairs(name.varcoef)
    npair.varcoef <- NCOL(pair.varcoef)
    dOmega.precomputed <- setNames(lapply(var.pattern, function(iPattern){
        iOut <- vector(mode = "list", length = npair.varcoef)
        for(iPair in 1:npair.varcoef){
            iOut[[iPair]] <- 0.5 * tr(precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[1,iPair]]] %*% precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[2,iPair]]])
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
                ## iOut$term3[[iPair]] <- TODO!!!
            }
            return(iOut)
        }), var.pattern)

        REML.det <- setNames(lapply(name.varcoef, function(i){
            list(numerator1 = matrix(0, nrow = n.varcoef, ncol = n.varcoef, dimnames = list(name.varcoef,name.varcoef)),
                 numerator2 = matrix(0, nrow = n.varcoef, ncol = n.varcoef, dimnames = list(name.varcoef,name.varcoef)),
                 numerator3 = matrix(0, nrow = n.varcoef, ncol = n.varcoef, dimnames = list(name.varcoef,name.varcoef)),
                 denominator = matrix(0, nrow = n.varcoef, ncol = n.varcoef, dimnames = list(name.varcoef,name.varcoef)))
        }), name.varcoef)
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

        for(iPair in 1:npair.varcoef){ ## iPair <- 1
            if(indiv){
                info[iId,pair.varcoef[1,iPair],pair.varcoef[2,iPair]] <- dOmega.precomputed[[iPattern]][[iPair]]
                if(pair.varcoef[1,iPair] != pair.varcoef[2,iPair]){
                    info[iId,pair.varcoef[2,iPair],pair.varcoef[1,iPair]] <- dOmega.precomputed[[iPattern]][[iPair]]
                }
            }else{
                info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] <- info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] + dOmega.precomputed[[iPattern]][[iPair]]
                if(pair.varcoef[1,iPair] != pair.varcoef[2,iPair]){
                    info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] <- info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] + dOmega.precomputed[[iPattern]][[iPair]]
                }
                browser()
                if(REML){
                    REML.det[[iPair]]$numerator1a <- REML.det[[iPair]]$numerator1a + t(iX) %*% (dOmega.precomputed2[[iPattern]]$term1a[[iPair]]) %*% iX
                    REML.det[[iPair]]$numerator1b <- REML.det[[iPair]]$numerator1b + t(iX) %*% (dOmega.precomputed2[[iPattern]]$term1b[[iPair]]) %*% iX
                    REML.det[[iPair]]$numerator2 <- REML.det[[iPair]]$numerator2 + t(iX) %*% (dOmega.precomputed2[[iPattern]]$term2[[iPair]]) %*% iX
                    ## REML.det[[iPair]]$numerator3 <- REML.det[[iPair]]$numerator3 + t(iX) %*% (dOmega.precomputed2[[iPattern]]$term3[[iVarcoef]]) %*% iX
                    REML.det[[iPair]]$denominator <- REML.det[[iPair]]$denominator + t(iX) %*% iOmegaM1 %*% iX
                }
            }
        }
    }

    ## ** export
    if(REML){
        for(iPair in 1:npair.varcoef){ ## iPair <- 1
            iDenomM1 <- solve(REML.det[[iPair]]$denominator)
            term1 <- 0.5 * tr(REML.det[[iPair]]$numerator1a %*% iDenomM1 %*% REML.det[[iPair]]$numerator1b)
            term2 <- - tr(iDenomM1 %*% REML.det[[iPair]]$numerator2)
            ## term3 <- 0.5 * tr(iDenomM1 %*% REML.det[[iPair]]$numerator3)
            info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] <- info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] + term1 + term2 + term3
            info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] <- info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] + term1 + term2 + term3
        }
    }
    return(info)

}

## * .unorderedPairs
## adapted from RecordLinkage package
.unorderedPairs <- function(x){
    n <- length(x)
    ls <- lapply(1:n, function(k){ rbind(x[k], x[k:n])})
    out <- array(unlist(ls), dim = c(2, n * (n + 1)/2))
    return(out)
}


##----------------------------------------------------------------------
### information.R ends here
