### precompute.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 22 2021 (13:47) 
## Version: 
## Last-Updated: okt  8 2026 (15:32) 
##           By: Brice Ozenne
##     Update #: 367
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
    for(iP in 1:n.pattern){ ## iPattern <- pattern[1]
        iN.time <- pattern.ntime[iP]

        if(iN.time == 1){
            iX.summary <- crossprod(X[unlist(index.cluster[pattern.cluster[[iP]]]),,drop=FALSE])
            out$pattern[[iP]][1,] <- iX.summary[lower.tri(iX.summary, diag = TRUE)]

        }else{
            iX <- array(unlist(lapply(index.cluster[pattern.cluster[[iP]]], function(iIndex){X[iIndex,,drop=FALSE]})),
                        dim = c(iN.time,NCOL(X),length(index.cluster[pattern.cluster[[iP]]])),
                        dimnames = list(NULL,colnames(X),NULL))

            for(iCol1 in 1:p){ ## iCol1 <- 1
                for(iCol2 in 1:iCol1){ ## iCol2 <- 2
                    out$pattern[[iP]][,out$key[iCol1,iCol2]] <- tcrossprod(iX[,iCol1,],iX[,iCol2,])
                }
            }
        }
        ## Possible alternative: 
        ## iX <- do.call(rbind,lapply(index.cluster[pattern.cluster[[iP]]], function(iIndex){as.vector(X[iIndex,,drop=FALSE])}))
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

    out <- stats::setNames(lapply(1:n.pattern, function(iP){
        matrix(0, nrow = pattern.ntime[iP]^2, ncol = p, dimnames = list(NULL,name.mucoef))
    }), pattern)

    for(iP in 1:n.pattern){ ## iPattern <- pattern[1]
        iN.time <- pattern.ntime[iP]
        iIndex.cluster <- index.cluster[pattern.cluster[[iP]]]

        if(iN.time == 1){
            out[[iP]][] <- crossprod(X[unlist(iIndex.cluster),,drop=FALSE], residuals[unlist(iIndex.cluster)])[,1]
        }else{
            iResiduals <- do.call(cbind, lapply(iIndex.cluster, function(iIndex){residuals[iIndex,,drop=FALSE]}))
            iX <- array(unlist(lapply(iIndex.cluster, function(iIndex){X[iIndex,,drop=FALSE]})),
                        dim = c(iN.time,NCOL(X),length(index.cluster[pattern.cluster[[iP]]])),
                        dimnames = list(NULL,colnames(X),NULL))
            for(iCol in 1:p){ ## iCol <- 3
                out[[iP]][,iCol] <- as.vector(tcrossprod(iX[,iCol,], iResiduals))
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

    for(iP in 1:n.pattern){ ## iPattern <- pattern[1]
        out[[iP]] <- as.vector(tcrossprod(do.call(cbind,lapply(index.cluster[pattern.cluster[[iP]]], function(iIndex){residuals[iIndex,,drop=FALSE]}))))
    }
    
    return(out)
}

## * .precomputeOmega
## Precompute product between the residual variance-covariance matrix and its derivative
.precomputeOmega <- function(precision, dOmega, d2Omega, d3Omega,
                             effects, pair.vcov, triplet.vcov,
                             REML, type.information, logLik, score, information, vcov, df){

    ## ** extract information
    pattern <- names(precision)
    n.pattern <- length(pattern)
    time.pattern <- lapply(precision, NROW)
    param.pattern <- lapply(dOmega,names)
    nparam.pattern <- sapply(param.pattern,length)    
    npair.vcov <- colSums(pair.vcov[,pattern,drop=FALSE])
    namepair.vcov <- stats::setNames(lapply(pattern, function(iPattern){pair.vcov[pair.vcov[[iPattern]],"name"]}),pattern)
    ntriplet.vcov <- colSums(triplet.vcov[,pattern,drop=FALSE])
    nametriplet.vcov <- stats::setNames(lapply(pattern, function(iPattern){triplet.vcov[triplet.vcov[[iPattern]],"name"]}),pattern)

    ## ** special case
    if(((score==FALSE) && (information==FALSE) && (vcov==FALSE) && (df==FALSE)) || (("variance" %in% effects == FALSE) && ("correlation" %in% effects == FALSE))){
        return(list()) ## empty list
    }else if(df && identical(attr(df,"method"),"numeric")){
        df <- FALSE
    }
    
    ## ** prepare output
    out <- list(dOmegaM1 = lapply(pattern, function(iP){matrix(NA, nrow = time.pattern[[iP]]^2, ncol = nparam.pattern[iP], dimnames = list(NULL, param.pattern[[iP]]))}))
    names(out$dOmegaM1) <- pattern
    
    if(score){
        out$dTrace <- lapply(pattern, function(iP){stats::setNames(rep(NA, nparam.pattern[iP]), param.pattern[[iP]])})
        names(out$dTrace) <- pattern
    }

    if((information || vcov || df)){
        if(type.information=="observed" || REML){
            out$d2OmegaM1 <- lapply(pattern, function(iP){matrix(NA, nrow = time.pattern[[iP]]^2, ncol = npair.vcov[iP], dimnames = list(NULL,namepair.vcov[[iP]]))})
            names(out$d2OmegaM1) <- pattern
        }
        if(information || vcov){
            out$d2Trace <- lapply(pattern, function(iP){stats::setNames(rep(NA, npair.vcov[iP]),namepair.vcov[[iP]])})
            names(out$d2Trace) <- pattern
        }
    }
    
    if(df){
        if(type.information=="observed" || REML){
            out$d3OmegaM1 <- lapply(pattern, function(iP){matrix(NA, nrow = time.pattern[[iP]]^2, ncol = ntriplet.vcov[iP], dimnames = list(NULL,nametriplet.vcov[[iP]]))})
            names(out$d3OmegaM1) <- pattern
        }
        out$d3Trace <- lapply(pattern, function(iP){stats::setNames(rep(NA, ntriplet.vcov[iP]),nametriplet.vcov[[iP]])})
        names(out$d3Trace) <- pattern
        
    }
    
    ## ** pre-compute
    for(iPattern in pattern){ ## iPattern <- pattern[1]

        ## *** handle special case: inverse not defined or no covariance parameter
        if (inherits(precision[[iPattern]], "try-error") || nparam.pattern[iPattern]==0) {
            next
        }

        ## *** first derivative
        iLS_dOmega.OmegaM1 <- lapply(dOmega[[iPattern]], function(iO){iO %*% precision[[iPattern]]})
        iLS_OmegaM1.dOmega.OmegaM1 <- lapply(iLS_dOmega.OmegaM1, function(iO){precision[[iPattern]] %*% iO})

        for(iP in param.pattern[[iPattern]]){ ## iP <- param.pattern[[iPattern]][1]
            if(score){
                out$dTrace[[iPattern]][iP] <- tr(iLS_dOmega.OmegaM1[[iP]])
            }
            out$dOmegaM1[[iPattern]][,iP] <- as.vector(iLS_OmegaM1.dOmega.OmegaM1[[iP]])
        }

        ## *** second derivative
        if((information || vcov || df)){
            iPair.vcov <- pair.vcov[pair.vcov[,iPattern],c("param1","param2")]

            for(iP2 in 1:npair.vcov[iPattern]){ ## iP2 <- 29
                iParam2.1 <- iPair.vcov[iP2,"param1"]
                iParam2.2 <- iPair.vcov[iP2,"param2"]

                if(type.information=="observed" || REML){
                    idOmega.OmegaM1.dOmega <- iLS_dOmega.OmegaM1[[iParam2.2]] %*% dOmega[[iPattern]][[iParam2.1]]                    
                    out$d2OmegaM1[[iPattern]][,iP2] <- as.double(precision[[iPattern]] %*% (d2Omega[[iPattern]][[iP2]] - idOmega.OmegaM1.dOmega - t(idOmega.OmegaM1.dOmega)) %*% precision[[iPattern]])
                }
                if(information || vcov){                    
                    if(type.information == "expected"){
                        out$d2Trace[[iPattern]][iP2] <- -sum(iLS_OmegaM1.dOmega.OmegaM1[[iParam2.2]] * dOmega[[iPattern]][[iParam2.1]])
                    }else if(type.information == "observed"){
                        out$d2Trace[[iPattern]][iP2] <- sum(iLS_OmegaM1.dOmega.OmegaM1[[iParam2.2]] * dOmega[[iPattern]][[iParam2.1]]) - sum(precision[[iPattern]] * d2Omega[[iPattern]][[iP2]]) 
                    }
                }
            }
        }

        ## *** third derivative
        if(df && (is.null(attr(df,"method")) || attr(df,"method")=="analytic")){
            iTriplet.vcov <- triplet.vcov[triplet.vcov[,iPattern],c("name","param1","param2","param3")]

            for(iP3 in 1:ntriplet.vcov[iPattern]){ ## iP3 <- 1
                iName <- iTriplet.vcov[iP2,"name"]
                
                iParam3.1 <- iTriplet.vcov[iP2,"param1"]
                iParam3.2 <- iTriplet.vcov[iP2,"param2"]
                iParam3.3 <- iTriplet.vcov[iP2,"param3"]

                Upsilon <- list(iLS_OmegaM1.dOmega.OmegaM1[[iParam3.1]] %*% (iLS_dOmega.OmegaM1[[iParam3.2]] %*% dOmega[[iPattern]][[iParam3.3]] - d2Omega[[iPattern]][[paste(iParam3.2,iParam3.3,sep="_")]]),
                                iLS_OmegaM1.dOmega.OmegaM1[[iParam3.1]] %*% (iLS_dOmega.OmegaM1[[iParam3.3]] %*% dOmega[[iPattern]][[iParam3.2]] - d2Omega[[iPattern]][[paste(iParam3.2,iParam3.3,sep="_")]]),
                                iLS_OmegaM1.dOmega.OmegaM1[[iParam3.2]] %*% (iLS_dOmega.OmegaM1[[iParam3.1]] %*% dOmega[[iPattern]][[iParam3.3]] - d2Omega[[iPattern]][[paste(iParam3.1,iParam3.3,sep="_")]]),
                                iLS_OmegaM1.dOmega.OmegaM1[[iParam3.2]] %*% (iLS_dOmega.OmegaM1[[iParam3.3]] %*% dOmega[[iPattern]][[iParam3.1]] - d2Omega[[iPattern]][[paste(iParam3.1,iParam3.3,sep="_")]]),
                                iLS_OmegaM1.dOmega.OmegaM1[[iParam3.3]] %*% (iLS_dOmega.OmegaM1[[iParam3.1]] %*% dOmega[[iPattern]][[iParam3.2]] - d2Omega[[iPattern]][[paste(iParam3.1,iParam3.2,sep="_")]]),
                                iLS_OmegaM1.dOmega.OmegaM1[[iParam3.3]] %*% (iLS_dOmega.OmegaM1[[iParam3.2]] %*% dOmega[[iPattern]][[iParam3.1]] - d2Omega[[iPattern]][[paste(iParam3.1,iParam3.2,sep="_")]])
                                )

                if(type.information=="observed" || REML){
                    out$d3OmegaM1[[iPattern]][,iP3] <- as.double((precision[[iPattern]] %*% d3Omega[[iPattern]][[iP3]] + Reduce("+",Upsilon)) %*% precision[[iPattern]])
                }
                if(information || vcov){                    
                    if(type.information == "expected"){
                        out$d3Trace[[iPattern]][iP3] <- sum(diag(Upsilon[[1]])) + sum(diag(Upsilon[[2]]))
                    }else if(type.information == "observed"){
                        out$d3Trace[[iPattern]][iP3] <- sum(c(sum(iLS_OmegaM1.dOmega.OmegaM1[[iParam3.3]] * d2Omega[[iPattern]][[paste(iParam3.1,iParam3.2,sep="_")]]),
                                                              sum(precision[[iPattern]] * d3Omega[[iPattern]][[iP3]]),
                                                              sum(diag(Upsilon[[1]])) + sum(diag(Upsilon[[3]])))
                                                              ) 
                    }
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
.precomputeREML <- function(precision, dOmega, d2Omega, d3Omega, 
                            effects, param.vcov, pair.vcov, triplet.vcov, precompute,
                            logLik, score, information, vcov, df){

    ## ** normalize user input
    if((score || information || vcov || df) && ("variance" %in% effects == FALSE) && ("correlation" %in% effects == FALSE)){
        score <- FALSE
        information <- FALSE
        vcov <- FALSE
        df <- FALSE
    }else if(df){
        if(identical(attr(df,"method"),"numeric")){
            df <- FALSE
        }else{
            score <- TRUE
            information <- TRUE
        }
    }else if(information || vcov){
        score <- TRUE
        information <- TRUE
    }

    ## ** prepare
    pattern <- names(precision)
    p <- NCOL(precompute$XX$key)
    param.mean <- colnames(precompute$XX$key)
    XX.key <- as.vector(precompute$XX$key)

    out <- list()
    X.OmegaM1.X <- rep(0, p*(p+1)/2)
    if(score){
        X.dOmegaM1.X <- matrix(0, nrow = p*(p+1)/2, ncol = length(param.vcov), dimnames = list(NULL,param.vcov))
    }
    if(information){
        pair.varcoef <- unique(unlist(lapply(d2Omega,names)))
        X.d2OmegaM1.X <- matrix(0, nrow = p*(p+1)/2, ncol = NROW(pair.vcov), dimnames = list(NULL,pair.vcov$name))
    }
    if(df){
        triplet.varcoef <- unique(unlist(lapply(d3Omega,names)))
        X.d3OmegaM1.X <- matrix(0, nrow = p*(p+1)/2, ncol = NROW(triplet.vcov), dimnames = list(NULL,triplet.vcov$name))
    }

    ## ** accumulate
    for (iPattern in pattern) { ## iPattern <- pattern[1]

        ## *** handle special case
        if (inherits(precision[[iPattern]], "try-error")) {
            next
        }
        
        ## *** evaluate matrix products
        iXX <- t(precompute$XX$pattern[[iPattern]])
        browser()
        X.OmegaM1.X <- X.OmegaM1.X + (iXX %*% cbind(as.vector(precision[[iPattern]])))[,1]
        if(score){
            iParam.vcov <- names(dOmega[[iPattern]])
            X.dOmegaM1.X[,iParam.vcov] <- X.dOmegaM1.X[,iParam.vcov,drop=FALSE] + iXX %*% precompute$Omega$dOmegaM1[[iPattern]]
        }

        if(information){
            iPair.vcov <- names(d2Omega[[iPattern]])
            X.d2OmegaM1.X[,iPair.vcov] <- X.d2OmegaM1.X[,iPair.vcov,drop=FALSE] + iXX %*% precompute$Omega$d2OmegaM1[[iPattern]]
        }
        
        if(df){
            iTriplet.vcov <- names(d3Omega[[iPattern]])
            X.d3OmegaM1.X[,iTriplet.vcov] <- X.d3OmegaM1.X[,iTriplet.vcov,drop=FALSE] + iXX %*% precompute$Omega$d3OmegaM1[[iPattern]]
        }
        
    }

    ## ** global operations and reshape
    if(logLik){
        X.OmegaM1.X_det <- det(matrix(X.OmegaM1.X[XX.key], nrow = p, ncol = p, dimnames = list(param.mean,param.mean)))
        if(!is.na(X.OmegaM1.X_det) && X.OmegaM1.X_det>0){
            out$logdet_X.OmegaM1.X <- log(X.OmegaM1.X_det)
        }else{
            out$logdet_X.OmegaM1.X <- NA
        }
    }
    if(score){
        out$X.dOmegaM1.X_M1 <- solve(X.dOmegaM1.X)
    }
    if(information){        
        out$X.d2OmegaM1.X <- apply(X.d2OmegaM1.X, MARGIN = 2, function(iRow){
            matrix(iRow[XX.key],nrow = p, ncol = p, dimnames = list(param.mean,param.mean))
        }, simplify = FALSE)
    }
    if(df){
        out$X.d3OmegaM1.X <- apply(X.d3OmegaM1.X, MARGIN = 2, function(iRow){
            matrix(iRow[XX.key],nrow = p, ncol = p, dimnames = list(param.mean,param.mean))
        }, simplify = FALSE)
    }
    
    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### precompute.R ends here
