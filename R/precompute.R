### precompute.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 22 2021 (13:47) 
## Version: 
## Last-Updated: okt  3 2024 (11:19) 
##           By: Brice Ozenne
##     Update #: 211
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .precomputeXX
## Precompute square of the design matrix
.precomputeXX <- function(X, pattern, pattern.ntime, pattern.cluster, index.cluster){

    p <- NCOL(X)
    n.pattern <- length(pattern)
    out <- list(pattern = stats::setNames(lapply(pattern, function(iPattern){matrix(0, nrow = pattern.ntime[iPattern]*pattern.ntime[iPattern], ncol = p*(p+1)/2)}), pattern),
                key = matrix(as.numeric(NA),nrow=p,ncol=p,dimnames=list(colnames(X),colnames(X))))

    ## ** prepare key
    out$key[lower.tri(out$key,diag = TRUE)] <- 1:sum(lower.tri(out$key,diag = TRUE))
    out$key[upper.tri(out$key)] <- t(out$key)[upper.tri(out$key)]

    ## ** fill matrix
    for(iPattern in pattern){ ## iPattern <- pattern[1]
        iN.time <- pattern.ntime[iPattern]

        if(iN.time == 1){
            iX.summary <- crossprod(X[unlist(index.cluster[pattern.cluster[[iPattern]]]),,drop=FALSE])
            out$pattern[[iPattern]][1,] <- iX.summary[lower.tri(iX.summary, diag = TRUE)]

        }else{
            iX <- array(unlist(lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){X[iIndex,,drop=FALSE]})),
                        dim = c(iN.time,NCOL(X),length(index.cluster[pattern.cluster[[iPattern]]])),
                        dimnames = list(NULL,colnames(X),NULL))

            for(iCol1 in 1:p){ ## iCol1 <- 1
                for(iCol2 in 1:iCol1){ ## iCol2 <- 2
                    out$pattern[[iPattern]][,out$key[iCol1,iCol2]] <- tcrossprod(iX[,iCol1,],iX[,iCol2,])
                }
            }
        }
        ## Possible alternative: 
        ## iX <- do.call(rbind,lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){as.vector(X[iIndex,,drop=FALSE])}))
        ## iX.summary <- crossprod(iX)
        ## Issue how to properly store the elements (too many of them T^2 p^2 instead of T^2 p(p+1)/2)
    }
    return(out)
}

## * .precomputeXR
## Precompute design matrix times residuals
.precomputeXR <- function(X, residuals, pattern, pattern.ntime, pattern.cluster, index.cluster){

    p <- NCOL(X)
    name.mucoef <- colnames(X)
    n.pattern <- length(pattern)

    out <- stats::setNames(lapply(pattern, function(iPattern){
        matrix(0, nrow = pattern.ntime[iPattern]^2, ncol = p, dimnames = list(NULL,name.mucoef))
    }), pattern)

    for(iPattern in pattern){ ## iPattern <- pattern[1]
        iN.time <- pattern.ntime[iPattern]
        iIndex.cluster <- index.cluster[pattern.cluster[[iPattern]]]

        if(iN.time == 1){
            out[[iPattern]][] <- crossprod(X[unlist(iIndex.cluster),,drop=FALSE], residuals[unlist(iIndex.cluster)])[,1]
        }else{
            iResiduals <- do.call(cbind, lapply(iIndex.cluster, function(iIndex){residuals[iIndex,,drop=FALSE]}))
            iX <- array(unlist(lapply(iIndex.cluster, function(iIndex){X[iIndex,,drop=FALSE]})),
                        dim = c(iN.time,NCOL(X),length(index.cluster[pattern.cluster[[iPattern]]])),
                        dimnames = list(NULL,colnames(X),NULL))
            for(iCol in 1:p){ ## iCol <- 3
                out[[iPattern]][,iCol] <- as.vector(tcrossprod(iX[,iCol,], iResiduals))
            }    
        }
    }

    return(out)
}

## * .precomputeRR
## Precompute square of the residuals
.precomputeRR <- function(residuals, pattern.ntime, pattern, pattern.cluster, index.cluster){

    n.pattern <- length(pattern)
    out <- stats::setNames(vector(mode = "list", length = length(pattern)), pattern)

    for(iPattern in pattern){ ## iPattern <- pattern[1]
        out[[iPattern]] <- as.vector(tcrossprod(do.call(cbind,lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){residuals[iIndex,,drop=FALSE]}))))
    }
    
    return(out)
}

## * .precomputeOmega
## Precompute product between the residual variance-covariance matrix and its derivative
.precomputeOmega <- function(precision, dOmega, d2Omega, effects, pair.vcov,
                             REML, type.information, logLik, score, information, vcov, df){

    ## ** extract information
    pattern <- names(precision)
    n.pattern <- length(pattern)
    time.pattern <- lapply(precision, NROW)
    param.pattern <- lapply(dOmega,names)
    n.param.pattern <- sapply(param.pattern,length)
    param2.pattern <- lapply(pair.vcov,colnames)
    n.param2.pattern <- sapply(param2.pattern,length)

    ## ** special case
    if(((score==FALSE) && (information==FALSE) && (vcov==FALSE) && (df==FALSE)) || (("variance" %in% effects == FALSE) && ("correlation" %in% effects == FALSE))){
        return(list()) ## empty list
    }
    
    ## ** prepare output
    out <- list(OmegaM1.dOmega.OmegaM1 = lapply(pattern, function(iP){matrix(NA, nrow = time.pattern[[iP]]^2, ncol = n.param.pattern[iP], dimnames = list(NULL, param.pattern[[iP]]))}))
    names(out$OmegaM1.dOmega.OmegaM1) <- pattern
    if(score){
        out$tr.OmegaM1.dOmega <- lapply(pattern, function(iP){stats::setNames(rep(NA, n.param.pattern[iP]), param.pattern[[iP]])})
        names(out$tr.OmegaM1.dOmega) <- pattern
    }
    if((information || vcov)){

        out$tr.OmegaM1.dOmega.OmegaM1.dOmega <- lapply(pattern, function(iP){stats::setNames(rep(NA, n.param2.pattern[iP]), param2.pattern[[iP]])})
        names(out$tr.OmegaM1.dOmega.OmegaM1.dOmega) <- pattern
        if(type.information == "observed"){
            out$tr.OmegaM1.d2Omega <- lapply(pattern, function(iP){stats::setNames(rep(NA, n.param2.pattern[iP]), param2.pattern[[iP]])})
            names(out$tr.OmegaM1.d2Omega) <- pattern
        }

        if(REML || type.information=="observed"){
            out$OmegaM1.d2OmegaAndCo.OmegaM1 <- lapply(pattern, function(iP){matrix(NA, nrow = time.pattern[[iP]]^2, ncol = n.param2.pattern[iP], dimnames = list(NULL, param2.pattern[[iP]]))})
            names(out$OmegaM1.d2OmegaAndCo.OmegaM1) <- pattern
        }
    }

    ## ** pre-compute
    for(iPattern in pattern){ ## iPattern <- pattern[1]

        ## *** handle special case
        if (inherits(precision[[iPattern]], "try-error") || n.param.pattern[iPattern]==0) {
            next
        }

        ## *** evaluate matrix products
        iLS_dOmega.OmegaM1 <- lapply(dOmega[[iPattern]], function(iO){iO %*% precision[[iPattern]]})
        iLS_OmegaM1.dOmega.OmegaM1 <- lapply(iLS_dOmega.OmegaM1, function(iO){precision[[iPattern]] %*% iO})

        for(iP in param.pattern[[iPattern]]){ ## iP <- param.pattern[[iPattern]][1]
            if(score){
                out$tr.OmegaM1.dOmega[[iPattern]][iP] <- tr(iLS_dOmega.OmegaM1[[iP]])
            }
            out$OmegaM1.dOmega.OmegaM1[[iPattern]][,iP] <- as.vector(iLS_OmegaM1.dOmega.OmegaM1[[iP]])
        }

        if((information || vcov)){
            for(iP2 in 1:n.param2.pattern[iPattern]){ ## iP2 <- 2
                iParam2 <- param2.pattern[[iPattern]][iP2]
                iParam2.1 <- pair.vcov[[iPattern]][1,iParam2]
                iParam2.2 <- pair.vcov[[iPattern]][2,iParam2]

                out$tr.OmegaM1.dOmega.OmegaM1.dOmega[[iPattern]][iParam2] <- sum(iLS_OmegaM1.dOmega.OmegaM1[[iParam2.2]] * dOmega[[iPattern]][[iParam2.1]])
                if(type.information == "observed"){
                    out$tr.OmegaM1.d2Omega[[iPattern]][iParam2] <- sum(precision[[iPattern]] * d2Omega[[iPattern]][[iParam2]])
                }
                if(REML || type.information=="observed"){
                    idOmega.OmegaM1.dOmega <- iLS_dOmega.OmegaM1[[iParam2.2]] %*% dOmega[[iPattern]][[iParam2.1]]
                    out$OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]][,iParam2] <- as.double(precision[[iPattern]] %*% (d2Omega[[iPattern]][[iParam2]] - idOmega.OmegaM1.dOmega - t(idOmega.OmegaM1.dOmega)) %*% precision[[iPattern]])
                }
            }
        }
    }

    ## ** export
    return(out)
}

## * .precomputeREML
## Precompute REML terms
## Note: possible weights are already included in precompute$XX
.precomputeREML <- function(precision, dOmega, d2Omega, precompute, effects,
                            logLik, score, information, vcov, df){

    ## ** normalize user input
    if((score || information || vcov) && ("variance" %in% effects == FALSE) && ("correlation" %in% effects == FALSE)){
        score <- FALSE
        information <- FALSE
        vcov <- FALSE
    }

    ## ** prepare
    pattern <- names(precision)
    p <- NCOL(precompute$XX$key)
    name.meancoef <- colnames(precompute$XX$key)
    name.varcoef <- unique(unlist(lapply(dOmega,names)))
    XX.key <- as.vector(precompute$XX$key)

    out <- list(X.OmegaM1.X = rep(0, p*(p+1)/2))
    if(score || information || vcov){
        out$X.OmegaM1.dOmega.OmegaM1.X <- matrix(0, nrow = p*(p+1)/2, ncol = length(name.varcoef), dimnames = list(NULL,name.varcoef))
    }
    if(information || vcov){
        name.varcoef2 <- unique(unlist(lapply(d2Omega,names)))
        out$X.OmegaM1.d2OmegaAndCo.OmegaM1.X <- matrix(0, nrow = p*(p+1)/2, ncol = length(name.varcoef2), dimnames = list(NULL,name.varcoef2))
    }

    ## ** accumulate
    for (iPattern in pattern) { ## iPattern <- pattern[1]

        ## *** handle special case
        if (inherits(precision[[iPattern]], "try-error")) {
            next
        }
        
        ## *** evaluate matrix products
        iXX <- t(precompute$XX$pattern[[iPattern]])

        out$X.OmegaM1.X <- out$X.OmegaM1.X + (iXX %*% cbind(attr(precision[[iPattern]], "vectorize")))[,1]
        if(score || information || vcov){
            iName.varcoef <- colnames(precompute$Omega$OmegaM1.dOmega.OmegaM1[[iPattern]])
            out$X.OmegaM1.dOmega.OmegaM1.X[,iName.varcoef] <- out$X.OmegaM1.dOmega.OmegaM1.X[,iName.varcoef,drop=FALSE] + iXX %*% precompute$Omega$OmegaM1.dOmega.OmegaM1[[iPattern]]
        }

        if(information || vcov){
            iName.varcoef2 <- colnames(precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]])
            out$X.OmegaM1.d2OmegaAndCo.OmegaM1.X[,iName.varcoef2] <- out$X.OmegaM1.d2OmegaAndCo.OmegaM1.X[,iName.varcoef2,drop=FALSE] + iXX %*% precompute$Omega$OmegaM1.d2OmegaAndCo.OmegaM1[[iPattern]]
        }
    }

    ## ** global operations and reshape
    out$X.OmegaM1.X <- matrix(out$X.OmegaM1.X[XX.key], nrow = p, ncol = p, dimnames = list(name.meancoef,name.meancoef))
    if(logLik){
        X.OmegaM1.X_det <- det(out$X.OmegaM1.X)
        if(!is.na(X.OmegaM1.X_det) && X.OmegaM1.X_det>0){
            out$logdet_X.OmegaM1.X <- log(X.OmegaM1.X_det)
        }else{
            out$logdet_X.OmegaM1.X <- NA
        }
    }
    if(score || information || vcov){
        out$X.OmegaM1.X_M1 <- solve(out$X.OmegaM1.X)
    }
    out$X.OmegaM1.X <- NULL
    if(score || information || vcov){
        out$X.OmegaM1.dOmega.OmegaM1.X <- stats::setNames(apply(out$X.OmegaM1.dOmega.OmegaM1.X, MARGIN = 2, function(iRow){
            matrix(iRow[XX.key],nrow = p, ncol = p, dimnames = list(name.meancoef,name.meancoef))
        }, simplify = FALSE), name.varcoef)
    }
    if(information || vcov){
        out$X.OmegaM1.d2OmegaAndCo.OmegaM1.X <- stats::setNames(apply(out$X.OmegaM1.d2OmegaAndCo.OmegaM1.X, MARGIN = 2, function(iRow){
            matrix(iRow[XX.key],nrow = p, ncol = p, dimnames = list(name.meancoef,name.meancoef))
        }, simplify = FALSE), name.varcoef2)
    }

    
    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### precompute.R ends here
