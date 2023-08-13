### precompute.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 22 2021 (13:47) 
## Version: 
## Last-Updated: aug  8 2023 (19:29) 
##           By: Brice Ozenne
##     Update #: 298
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
                key = matrix(as.numeric(NA),nrow=p,ncol=p,dimnames=list(colnames(X),colnames(X))),
                Xpattern = stats::setNames(vector(mode = "list", length = n.pattern),pattern))

    ## key
    out$key[lower.tri(out$key,diag = TRUE)] <- 1:sum(lower.tri(out$key,diag = TRUE))
    out$key[upper.tri(out$key)] <- t(out$key)[upper.tri(out$key)]

    ## fill matrix
    for(iPattern in pattern){ ## iPattern <- pattern[1]
        iTime <- pattern.ntime[iPattern]
        if(iTime==1){
            out$Xpattern[[iPattern]] <- do.call(rbind,lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){X[iIndex,,drop=FALSE]}))
            iX.summary <- crossprod(out$Xpattern[[iPattern]])
            ## out$key[lower.tri(out$key,diag = TRUE)]
            out$pattern[[iPattern]][1,] <- iX.summary[lower.tri(iX.summary, diag = TRUE)]
        }else{
            out$Xpattern[[iPattern]] <- array(unlist(lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){X[iIndex,,drop=FALSE]})),
                                              dim = c(iTime,NCOL(X),length(index.cluster[pattern.cluster[[iPattern]]])),
                                              dimnames = list(NULL,colnames(X),NULL))

            for(iCol1 in 1:p){ ## iCol1 <- 1
                for(iCol2 in 1:iCol1){ ## iCol2 <- 2
                    ## for(iId in pattern.cluster[[iPattern]]){
                    ##     out$pattern[[iPattern]][,,out$key[iCol1,iCol2]] <- out$pattern[[iPattern]][,,out$key[iCol1,iCol2]] + tcrossprod(X[index.cluster[[iId]],iCol1,drop=FALSE],X[index.cluster[[iId]],iCol2,drop=FALSE])
                    ## }
                    out$pattern[[iPattern]][,out$key[iCol1,iCol2]] <- tcrossprod(out$Xpattern[[iPattern]][,iCol1,],out$Xpattern[[iPattern]][,iCol2,])
                }
            }
        }
    }
    return(out)
}

## * .precomputeXR
## Precompute design matrix times residuals
.precomputeXR <- function(X, residuals, pattern, pattern.ntime, pattern.cluster, index.cluster){
    p <- NCOL(X[[1]])
    name.mucoef <- colnames(X[[1]])
    n.pattern <- length(pattern)

    out <- stats::setNames(lapply(pattern, function(iPattern){
        array(0, dim = c(pattern.ntime[iPattern], pattern.ntime[iPattern], ncol = p), dimnames = list(NULL,NULL,name.mucoef))
    }), pattern)

    for(iPattern in pattern){ ## iPattern <- pattern[1]

        iTime <- pattern.ntime[iPattern]
        iResiduals <- do.call(cbind, lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){residuals[iIndex,,drop=FALSE]}))
        iX <- X[[iPattern]]

        if(iTime == 1){
            out[[iPattern]][1,1,] <- iResiduals %*% iX
        }else{
            for(iCol in 1:p){ ## iCol1 <- 1
                ## for(iId in 1:length(pattern.cluster[[iPattern]])){ ## iId <- 1
                ##     out[[iPattern]][,,iCol] <- out[[iPattern]][,,iCol] + tcrossprod(iX[[iId]][,iCol,drop=FALSE],iResiduals[[iId]])
                ## }
                out[[iPattern]][,,iCol] <- tcrossprod(iX[,iCol,], iResiduals)
            }
        }

        attr(out[[iPattern]],"vectorwise") <- do.call(rbind,apply(out[[iPattern]], MARGIN = 3, FUN = as.vector, simplify = FALSE))
    }

    return(out)
}

## * .precomputeRR
## Precompute square of the residuals
.precomputeRR <- function(residuals, pattern.ntime, pattern, pattern.cluster, index.cluster){

    n.pattern <- length(pattern)
    out <- stats::setNames(lapply(pattern, function(iPattern){
        matrix(0, nrow = pattern.ntime[iPattern], ncol = pattern.ntime[iPattern])
    }), pattern)

    for(iPattern in pattern){ ## iPattern <- pattern[1]
        
        ## for(iId in pattern.cluster[[iPattern]]){
        ##     out[[iPattern]] <- out[[iPattern]] + tcrossprod(residuals[index.cluster[[iId]],,drop=FALSE])
        ## }
        out[[iPattern]] <- tcrossprod(do.call(cbind,lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){residuals[iIndex,,drop=FALSE]})))
        attr(out[[iPattern]],"vectorwise") <- as.vector(out[[iPattern]])
    }
    return(out)
}

## * .precomputeOmega
## Precompute functions of the residual variance-covariance matrix and its derivates
.precomputeOmega <- function(precompute, Omega, dOmega, d2Omega, 
                             method.fit, type.information, 
                             score, information, vcov, df){

    ## ** extract information
    n.pattern <- length(Omega)
    name.pattern <- names(Omega)
    if(method.fit == "REML"){
        vec.key <- as.vector(precompute$wXX$key)
        name.muparam <- colnames(precompute$wXX$key)
        n.muparam <- length(name.muparam)
        name.vcovparam <- attr(dOmega,"param")
        n.vcovparam <- length(name.vcovparam)
        name.pattern <- names(Omega)
        n.pattern <- length(name.pattern)
        allpair.vcov <- attr(d2Omega,"pair")
        nAllpair.vcov <- NCOL(allpair.vcov)
        precompute.moments <- attr(precompute,"moments")
    }

    ## ** prepare output
    out <- list(Omega.chol = stats::setNames(vector(mode = "list", length = n.pattern), name.pattern),
                OmegaM1 = stats::setNames(vector(mode = "list", length = n.pattern), name.pattern),
                Omega.logdet = stats::setNames(rep(NA, times = n.pattern), name.pattern))
    if(method.fit=="REML"){
        wXOmegaM1X <- matrix(0, nrow = n.muparam, ncol = n.muparam, dimnames = list(name.muparam, name.muparam))
        out$REML <- list(logdet = 0)
    }
    
    if(!is.null(dOmega) && (score || information || vcov || df)){
        ## !is.null(dOmega) for the case where only the mean parameters for the scores are asked

        ## trace
        out$dOmega_OmegaM1 <- matrix(NA, nrow = n.pattern, ncol = n.vcovparam, dimnames = list(name.pattern, name.vcovparam))
        ## full matrix
        out$OmegaM1_dOmega_OmegaM1 <- stats::setNames(vector(mode = "list", length = n.pattern), name.pattern)

        if(method.fit=="REML"){
            wXOmegaM1X.M1 <- matrix(0, nrow = n.muparam, ncol = n.muparam, dimnames = list(name.muparam, name.muparam))
            wXOmegaM1dOmegaOmegaM1X <- matrix(0, nrow = n.muparam, ncol = n.muparam, dimnames = list(name.muparam, name.muparam))
            out$REML$dOmega <- stats::setNames(rep(0, times = n.vcovparam), name.vcovparam)            
        }        
    }

    if(information || vcov || df){
        ## trace
        out$OmegaM1_dOmega_OmegaM1_dOmega <- array(NA, dim = c(n.pattern, n.vcovparam, n.vcovparam), dimnames = list(name.pattern, name.vcovparam, name.vcovparam))

        if(df || (method.fit == "REML" || type.information == "observed")){
            ## full matrix
            out$d2Omega_dOmega_OmegaM1_dOmega <- stats::setNames(vector(mode = "list", length = n.pattern), name.pattern)
            pair.vcov <- attr(d2Omega,"pair")
        }

        if(method.fit=="REML"){
            wXOmegaM1d2OmegaAndCoOmegaM1X <- matrix(0, nrow = n.muparam, ncol = n.muparam, dimnames = list(name.muparam, name.muparam))
            out$REML$d2OmegaAndCo <- stats::setNames(rep(0, times = n.vcovparam), name.vcovparam)            
        }        
    }

    if(df){
        browser()
    }

    ## ** Likelihood
    index.pattern <- 1:n.pattern
    for(iPattern in 1:n.pattern){ ## iPattern <- 1
        iOmega <- Omega[[iPattern]]
        attr(iOmega,"sd") <- NULL
        attr(iOmega,"cor") <- NULL

        ## *** Omega^{1/2}
        out$Omega.chol[[iPattern]] <- try(chol(iOmega), silent = TRUE)

        if(!inherits(out$Omega.chol[[iPattern]],"try-error")){
            ## *** det(Omega)
            out$Omega.logdet[iPattern] <- 2*sum(log(diag(out$Omega.chol[[iPattern]])))
            ## log(det(iOmega)) - out$Omega.logdet[[iPattern]]
        
            ## *** Omega^{-1}
            out$OmegaM1[[iPattern]] <- chol2inv(out$Omega.chol[[iPattern]])
            attr(out$OmegaM1[[iPattern]],"vectorwise") <- as.vector(out$OmegaM1[[iPattern]])

            ## *** REML
            if(method.fit=="REML"){
                if(precompute.moments){
                    iwXX <- precompute$wXX$pattern[[iPattern]]
                    wXOmegaM1X <- wXOmegaM1X + (attr(out$OmegaM1[[iPattern]],"vectorwise") %*% iwXX)[vec.key]
                }else{
                    iwX <- precompute$wX[[iPattern]]
                    wXOmegaM1X <- wXOmegaM1X + Reduce("+", lapply(iwX, function(iX){t(iX) %*% precompute$OmegaM1[[iPattern]] %*% iX}))
                }
            }
        }else{
            index.pattern <- setdiff(index.pattern, iPattern)
        }
        ## iOmega %*% out$OmegaM1[[iPattern]]
        ## out$OmegaM1[[iPattern]] %*% iOmega
    }

    if(method.fit=="REML"){
        out$REML$logdet <- log(det(wXOmegaM1X))
    }
    if(!information && !vcov && !df && (!score || is.null(dOmega))){
        return(out)
    }
    
    ## *** Score
    if(method.fit=="REML"){
        wXOmegaM1X.M1 <- solve(wXOmegaM1X)
        attr(wXOmegaM1X.M1,"vectorwise") <- as.vector(wXOmegaM1X.M1)
    }
    
    for(iParam1 in name.vcovparam){
    
        for(iPattern in index.pattern){
            dOmega_OmegaM1 <- dOmega[[iPattern]][[iParam1]] %*% out$OmegaM1[[iPattern]]
            out$dOmega_OmegaM1[iPattern,iParam1] <- tr(dOmega_OmegaM1[[iPattern]])

            out$OmegaM1_dOmega_OmegaM1[[iPattern]][[iParam1]] <- out$OmegaM1[[iPattern]] %*% dOmega_OmegaM1
        }

        if(method.fit == "REML"){
            if(precompute.moments){
                wXOmegaM1dOmegaOmegaM1X <- wXOmegaM1dOmegaOmegaM1X + (as.double(out$OmegaM1_dOmega_OmegaM1[[iPattern]][[iParam1]]) %*% iwXX)[vec.key]
            }else{
                wXOmegaM1dOmegaOmegaM1X <- wXOmegaM1dOmegaOmegaM1X + t(wX[[iPattern]]) %*% precompute$OmegaM1_dOmega_OmegaM1[[iPattern]][[iiParam]] %*% wX[[iPattern]]
            }
        }
    }

    if(!information && !vcov && !df){
        return(out)
    }

    ## *** Information
    for(iPair in xxx){
        browser()
    }

    if(!df){
        return(out)
    }

    ## *** Df (i.e. -dInformation)
    if(df){
        for(iTriplet in yyy){
            browser()
        }
    }
    ## for(iPattern in 1:n.pattern){ ## iPattern <- 1
    ##     iOmega <- Omega[[iPattern]]
    ##     attr(iOmega,"sd") <- NULL
    ##     attr(iOmega,"cor") <- NULL

        
        
    ##     ## *** Omega^{-1} d(Omega)/d(theta') Omega^{-1} d(Omega)/d(theta)
    ##     if(information || vcov || df){

    ##         iPair.vcov <- pair.vcov[[iPattern]]
    ##         iN.pair <- NCOL(iPair.vcov)

    ##         if(df || ((vcov || information) && (method.fit == "REML" || type.information == "observed"))){
    ##             out$OmegaM1_dOmega_OmegaM1_dOmega[[iPattern]] <- stats::setNames(lapply(1:iN.pair, function(iPair){ ## iPair <- 1
    ##                 return(out$OmegaM1_dOmega_OmegaM1[[iPattern]][[iPair.vcov[2,iPair]]] %*% dOmega[[iPattern]][[iPair.vcov[1,iPair]]])
    ##             }), colnames(iPair.vcov))
    ##             attr(out$OmegaM1_dOmega_OmegaM1_dOmega[[iPattern]],"tr") <- sapply(out$OmegaM1_dOmega_OmegaM1_dOmega[[iPattern]],tr)

    ##         }else{
    ##             out$OmegaM1_dOmega_OmegaM1_dOmega[[iPattern]] <- list()
    ##             attr(out$OmegaM1_dOmega_OmegaM1_dOmega[[iPattern]],"tr") <- stats::setNames(sapply(1:iN.pair, function(iPair){ ## fast trace calculation tr(AB) = sum(A*B)
    ##                 return(sum(out$OmegaM1_dOmega_OmegaM1[[iPattern]][[iPair.vcov[2,iPair]]] * dOmega[[iPattern]][[iPair.vcov[1,iPair]]]))
    ##             }), colnames(iPair.vcov))
    ##         }
    ##     }

    ##     ## *** Omega^{-1} (d2Omega/d(theta'theta) - d(Omega)/d(theta') Omega^{-1} d(Omega)/d(theta)) Omega^{-1}
    ##     if(df || ((vcov || information) && (method.fit == "REML" || type.information == "observed"))){
    ##         out$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]] <- stats::setNames(lapply(1:iN.pair, function(iPair){ ## iPair <- 5                
    ##             iOmega_M1_d2Omega <- out$OmegaM1[[iPattern]] %*% d2Omega[[iPattern]][[iPair]]
    ##             term1 <- out$OmegaM1_dOmega_OmegaM1_dOmega[[iPattern]][[iPair]] %*% out$OmegaM1[[iPattern]]
    ##             iOut <- iOmega_M1_d2Omega %*% out$OmegaM1[[iPattern]] - term1 - t(term1)
    ##             attr(iOut,"trModified") <- tr(out$OmegaM1_dOmega_OmegaM1_dOmega[[iPattern]][[iPair]] - iOmega_M1_d2Omega)
    ##             return(iOut)
    ##         }), colnames(iPair.vcov))
    ##         attr(out$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]], "trModified") <- sapply(out$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]], attr, "trModified")
    ##         attr(out$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]], "vectorwise") <- do.call(rbind,lapply(out$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]], function(iM){ as.vector(iM) }))
    ##     }
    ## }

    ## ** export
    return(out)

}

## * .precomputeREML
## Precompute functions of the residual variance-covariance matrix and the design matrix
.precomputeREML <- function(wX, precompute, 
                            Omega, dOmega, d2Omega,
                            score, information, vcov, df){

    ## ** extract information

    ## ** prepare output
    wXOmegaM1X <- matrix(0, nrow = n.muparam, ncol = n.muparam, dimnames = list(name.muparam, name.muparam))
    out <- list(wXOmegaM1X.logdet = 0)

    if(score || information || vcov || df){

        out$wXOmegaM1X.M1 <- matrix(0, nrow = n.muparam, ncol = n.muparam, dimnames = list(name.muparam, name.muparam))
        
        if(is.null(wX)){
            wXOmegaM1dOmegaOmegaM1X <- matrix(0, nrow = n.vcovparam, ncol = n.muparam*(n.muparam+1)/2, dimnames = list(name.vcovparam, NULL))
        }else{
            out$wXOmegaM1dOmegaOmegaM1X <- stats::setNames(lapply(name.vcovparam, FUN = function(iParam){
                matrix(0, n.muparam, n.muparam, dimnames = list(name.muparam, name.muparam))
            }), name.vcovparam)
        }

    }

    if(information || vcov || df){
        if(is.null(wX)){
            wXOmegaM1.d2OmegaAndCo.OmegaM1X <- matrix(0, nrow = nAllpair.vcov, ncol = n.muparam*(n.muparam+1)/2, dimnames = list(colnames(allpair.vcov), NULL))
        }else{
            out$wXOmegaM1.d2OmegaAndCo.OmegaM1X <- stats::setNames(lapply(1:nAllpair.vcov, FUN = function(iParam){
                matrix(0, n.muparam, n.muparam, dimnames = list(name.muparam, name.muparam))
            }), colnames(allpair.vcov))
        }
        
    }

    ## ** loop over patterns
    for(iPattern in 1:n.pattern){ ## iPattern <- 1
        if(inherits(precompute$Omega.chol[[iPattern]],"try-error")){next}
        iwXX <- precompute$wXX$pattern[[iPattern]]
        iParam.vcov <- precompute$OmegaM1_dOmega_O
megaM1[[iPattern]]

        if(is.null(wX)){
            wXOmegaM1X <- wXOmegaM1X + (attr(precompute$OmegaM1[[iPattern]],"vectorwise") %*% iwXX)[vec.key]
        }else{
            wXOmegaM1X <- wXOmegaM1X + Reduce("+", lapply(wX[[iPattern]], function(iX){t(iX) %*% precompute$OmegaM1[[iPattern]] %*% iX}))
        }        
        
        if(score || information || vcov || df){
            iParam <- names(precompute$OmegaM1_dOmega_OmegaM1[[iPattern]])
            if(is.null(wX)){
                wXOmegaM1dOmegaOmegaM1X[iParam,] <- wXOmegaM1dOmegaOmegaM1X[iParam,,drop=FALSE] + (attr(precompute$OmegaM1_dOmega_OmegaM1[[iPattern]],"vectorwise") %*% iwXX)
            }else{
                for(iiParam in iParam){
                    out$wXOmegaM1dOmegaOmegaM1X[[iiParam]] <- out$wXOmegaM1dOmegaOmegaM1X[[iiParam]] + Reduce("+", lapply(wX[[iPattern]], function(iX){
                        t(iX) %*% precompute$OmegaM1_dOmega_OmegaM1[[iPattern]][[iiParam]] %*% iX
                    }))
                }
            }
        }

        if(information || vcov || df){
            iPair.vcov <- names(precompute$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]])
            if(is.null(wX)){
                wXOmegaM1.d2OmegaAndCo.OmegaM1X[iPair.vcov,] <- wXOmegaM1.d2OmegaAndCo.OmegaM1X[iPair.vcov,,drop=FALSE] + (attr(precompute$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]], "vectorwise") %*% iwXX)
            }else{
                for(iiPair in iPair.vcov){
                    out$wXOmegaM1.d2OmegaAndCo.OmegaM1X[[iiPair]] <- out$wXOmegaM1.d2OmegaAndCo.OmegaM1X[[iiPair]] + Reduce("+", lapply(wX[[iPattern]], function(iX){
                        t(iX) %*% precompute$d2Omega_dOmega_OmegaM1_dOmega[[iPattern]][[iiPair]] %*% iX
                    }))
                }
            }
        }
    }

    ## ** reformat
    out$wXOmegaM1X.logdet <- log(det(wXOmegaM1X))
    if(score){
        out$wXOmegaM1X.M1 <- solve(wXOmegaM1X)
        attr(out$wXOmegaM1X.M1,"vectorwise") <- as.vector(out$wXOmegaM1X.M1)
            }

    if(is.null(wX)){
        if(score || information || vcov || df){
            out$wXOmegaM1dOmegaOmegaM1X <- apply(wXOmegaM1dOmegaOmegaM1X, MARGIN = 1, FUN = function(iRow){
                matrix(iRow[vec.key], n.muparam, n.muparam, dimnames = list(name.muparam, name.muparam))
            }, simplify = FALSE)
        }

        if(information || vcov || df){
            out$wXOmegaM1.d2OmegaAndCo.OmegaM1X <- apply(wXOmegaM1.d2OmegaAndCo.OmegaM1X, MARGIN = 1, FUN = function(iRow){
                matrix(iRow[vec.key], n.muparam, n.muparam, dimnames = list(name.muparam, name.muparam))
            }, simplify = FALSE)
        }
    }

    if(score || information || vcov || df){
        attr(out$wXOmegaM1dOmegaOmegaM1X,"vectorwise") <- do.call(rbind,lapply(out$wXOmegaM1dOmegaOmegaM1X, as.vector))
    }

    ## ** export
    return(out)
}

## * .precomputePairParam
.precomputePairParam <- function(structure, skeleton.param, type){
    
    param.constraint <- skeleton.param[!is.na(skeleton.param$constraint),"name"]
    if(type == "meanvcov"){
        skeleton.mu <- skeleton.param[skeleton.param$type=="mu",,drop=FALSE]
        freeparam.mu <- setdiff(skeleton.mu$name, param.constraint)
    }
    skeleton.vcov <- skeleton.param[skeleton.param$type %in% c("sigma","k","rho"),,drop=FALSE]
    freeparam.vcov <- setdiff(skeleton.vcov$name, param.constraint)
    
    skeleton.type <- c(mu = "mu", sigma = "sigmak", k = "sigmak", rho = "rho")[skeleton.param$type]
    param2type <- stats::setNames(skeleton.type, skeleton.param$name)
    
    ## ** form all pairs
    if(type == "meanvcov"){
        allpair <- unname(t(expand.grid(freeparam.mu, freeparam.vcov)))
        colnames(allpair) <- nlme::collapse(t(allpair), as.factor = FALSE)
        if(any(duplicated(colnames(allpair)))){
            stop("Something went wrong when identifying pairs between mean and variance-covariance parameters. \n",
                 "Duplicated names. \n")
        }
        attr(allpair, "key") <- matrix(NA, nrow = length(freeparam.mu), ncol = length(freeparam.vcov), dimnames = list(freeparam.mu,freeparam.vcov))
        for(iCol in 1:NCOL(allpair)){
            attr(allpair, "key")[allpair[1,iCol],allpair[2,iCol]] <- iCol
        }
            
    }else if(type == "vcov"){
        allpair <- unorderedPairs(freeparam.vcov)
        colnames(allpair) <- nlme::collapse(t(allpair), as.factor = FALSE)
        if(any(duplicated(colnames(allpair)))){
            stop("Something went wrong when identifying pairs between variance-covariance parameters. \n",
                 "Duplicated names. \n")
        }
        attr(allpair, "key") <- matrix(NA, nrow = length(freeparam.vcov), ncol = length(freeparam.vcov), dimnames = list(freeparam.vcov,freeparam.vcov))
        for(iCol in 1:NCOL(allpair)){
            attr(allpair, "key")[allpair[1,iCol],allpair[2,iCol]] <- iCol
            attr(allpair, "key")[allpair[2,iCol],allpair[1,iCol]] <- iCol
        }
        
    }

    ## ** form all pairs per pattern
    Upattern <- structure$Upattern
    if(type == "vcov"){

        if(is.null(structure$var) && is.null(structure$cor)){
            Mindicator <- NULL
        }else{
            if(is.null(structure$cor)){
                Mindicator <- stats::setNames(lapply(structure$var$Xpattern, function(iX){attr(iX,"ls.index.pair")}), Upattern$name)
            }else if(is.null(structure$var)){
                Mindicator <- stats::setNames(lapply(structure$cor$Xpattern, function(iX){attr(iX,"ls.index.pair")}), Upattern$name)
            }else if(!is.null(structure$var) && !is.null(structure$cor)){
                Mindicator <- stats::setNames(mapply(x = structure$var$Xpattern, y = structure$cor$Xpattern, function(x,y){
                    c(attr(x, "ls.index.pair"), attr(y, "ls.index.pair"))
                }, SIMPLIFY = FALSE), Upattern$name)
            }
        }
    }

    out <- stats::setNames(lapply(Upattern$name, function(iName){## iName <- Upattern$name[1]
        iPattern <- Upattern[Upattern$name == iName,,drop=FALSE]
        iParamVcov <- setdiff(iPattern$param[[1]], param.constraint)
        if(length(iParamVcov)==0){return(NULL)}

        if(type == "meanvcov"){
            iParamMu <- setdiff(c(skeleton.mu$name[skeleton.mu$index.strata == iPattern$index.strata], skeleton.mu$name[is.na(skeleton.mu$index.strata)]), param.constraint)
            if(length(iParamMu)==0){return(NULL)}
            iOut <- unname(t(expand.grid(iParamMu, iParamVcov)))
            colnames(iOut) <- nlme::collapse(t(iOut), as.factor = FALSE)
            attr(iOut, "key") <- matrix(NA, nrow = length(iParamMu), ncol = length(iParamVcov), dimnames = list(iParamMu,iParamVcov))
            for(iCol in 1:NCOL(iOut)){ ## iCol <- 1
                attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
            }
        }else if(type == "vcov"){
            iOut.all <- unorderedPairs(iParamVcov)
            colnames(iOut.all) <- nlme::collapse(t(iOut.all), as.factor = FALSE)
            if(is.null(Mindicator)){
                iOut <- iOut.all
            }else{
                iOut.Mindicator <- mapply(x = Mindicator[[iName]][iOut.all[1,]], y = Mindicator[[iName]][iOut.all[2,]], function(x,y){
                    xy.position <- intersect(x$position,y$position)
                    if(length(xy.position)==0){
                        return(NULL)
                    }else{
                        index.x <- match(xy.position, x$position)
                        index.y <- match(xy.position, y$position)
                        data.frame(position = xy.position,
                                   value1 = x[index.x,"value"],
                                   value2 = y[index.y,"value"],
                                   dvalue = x[index.x,"value"] * (y[index.y,"value"]-(x[index.x,"param"]==y[index.y,"param"])))
                    }
                }, SIMPLIFY = FALSE)                
                iTest.rm <- lengths(iOut.Mindicator)
                iOut <- iOut.all[,iTest.rm>0,drop=FALSE] ## remove combinaisons with no common times
                attr(iOut, "index.pair") <- stats::setNames(iOut.Mindicator[iTest.rm>0], colnames(iOut))
            }
            attr(iOut, "key") <- matrix(0, nrow = length(iParamVcov), ncol = length(iParamVcov), dimnames = list(iParamVcov,iParamVcov))
            for(iCol in 1:NCOL(iOut)){
                attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
                attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol]] <- iCol
            }
            attr(iOut,"type") <- matrix(param2type[as.vector(iOut)], nrow = 2, ncol = NCOL(iOut),
                                        dimnames = dimnames(iOut))
            collapse.type <- nlme::collapse(t(attr(iOut,"type")), as.factor = FALSE)
            attr(iOut,"typetype") <- tapply(1:length(collapse.type),collapse.type,base::identity,simplify = FALSE)
        } 
        attr(iOut,"index.global") <- match(colnames(iOut), colnames(allpair))        
        return(iOut)
    }),Upattern$name)

    
    ## ** export
    attr(out,"global") <- allpair
    return(out)
}

## * .precomputeTripletParam
.precomputeTripletParam <- function(structure){

    skeleton.vcov <- structure$param
    param.constraint <- skeleton.vcov[!is.na(skeleton.vcov$constraint),"name"]
    freeparam.vcov <- setdiff(skeleton.vcov$name, param.constraint)
    
    skeleton.type <- c(sigma = "sigmak", k = "sigmak", rho = "rho")[skeleton.vcov$type]
    param2type <- stats::setNames(skeleton.type, skeleton.vcov$name)
    
    ## ** form all triplets
    alltriplet <- unorderedTriplet(freeparam.vcov)
    colnames(alltriplet) <- nlme::collapse(t(alltriplet), as.factor = FALSE)
    if(any(duplicated(colnames(alltriplet)))){
        stop("Something went wrong when identifying triplet between variance-covariance parameters. \n",
             "Duplicated names. \n")
    }
    attr(alltriplet, "key") <- array(NA, dim = rep(length(freeparam.vcov),3), dimnames = list(freeparam.vcov,freeparam.vcov, freeparam.vcov))
    for(iElement in 1:NCOL(alltriplet)){ ## 
        attr(alltriplet, "key")[alltriplet[1,iElement],alltriplet[2,iElement],alltriplet[3,iElement]] <- iElement
        attr(alltriplet, "key")[alltriplet[1,iElement],alltriplet[3,iElement],alltriplet[2,iElement]] <- iElement
        attr(alltriplet, "key")[alltriplet[2,iElement],alltriplet[1,iElement],alltriplet[3,iElement]] <- iElement
        attr(alltriplet, "key")[alltriplet[2,iElement],alltriplet[3,iElement],alltriplet[1,iElement]] <- iElement
        attr(alltriplet, "key")[alltriplet[3,iElement],alltriplet[1,iElement],alltriplet[2,iElement]] <- iElement
        attr(alltriplet, "key")[alltriplet[3,iElement],alltriplet[2,iElement],alltriplet[1,iElement]] <- iElement
    }

    ## ** form all triplets per pattern
    Upattern <- structure$Upattern
    if(is.null(structure$var) && is.null(structure$cor)){
        Mindicator <- NULL
    }else{
        if(is.null(structure$cor)){
            Mindicator <- stats::setNames(lapply(structure$var$Xpattern, function(iX){attr(iX,"ls.index.pair")}), Upattern$name)
        }else if(is.null(structure$var)){
            Mindicator <- stats::setNames(lapply(structure$cor$Xpattern, function(iX){attr(iX,"ls.index.pair")}), Upattern$name)
        }else if(!is.null(structure$var) && !is.null(structure$cor)){
            Mindicator <- stats::setNames(mapply(x = structure$var$Xpattern, y = structure$cor$Xpattern, function(x,y){
                c(attr(x, "ls.index.pair"), attr(y, "ls.index.pair"))
            }, SIMPLIFY = FALSE), Upattern$name)
        }
        if(all(lengths(Mindicator)==0)){ ## CUSTOM structure without ls.index.triplet
            Mindicator <- NULL
        }
    }

    out <- stats::setNames(lapply(Upattern$name, function(iName){## iName <- Upattern$name[1]
        iPattern <- Upattern[Upattern$name == iName,,drop=FALSE]
        iParamVcov <- setdiff(iPattern$param[[1]], param.constraint)
        if(length(iParamVcov)==0){return(NULL)}

        iOut.all <- unorderedTriplet(iParamVcov)
        colnames(iOut.all) <- nlme::collapse(t(iOut.all), as.factor = FALSE)
        if(is.null(Mindicator)){
            iOut <- iOut.all
        }else{
            iOut.Mindicator <- mapply(x = Mindicator[[iName]][iOut.all[1,]], y = Mindicator[[iName]][iOut.all[2,]], z = Mindicator[[iName]][iOut.all[3,]], function(x,y,z){
                xyz.position <- intersect(intersect(x$position,y$position),z$position)
                if(length(xyz.position)==0){
                    return(NULL)
                }else{
                    index.x <- match(xyz.position, x$position)
                    index.y <- match(xyz.position, y$position)
                    index.z <- match(xyz.position, z$position)
                    data.frame(position = xyz.position,
                               value1 = x[index.x,"value"],
                               value2 = y[index.y,"value"],
                               value3 = z[index.z,"value"],
                               d2value = x[index.x,"value"] * (y[index.y,"value"] - (x[index.x,"param"]==y[index.y,"param"])) * (z[index.z,"value"] - (x[index.x,"param"]==z[index.z,"param"]) - (y[index.y,"param"]==z[index.z,"param"])))
                }
            }, SIMPLIFY = FALSE)
            iTest.rm <- lengths(iOut.Mindicator)
            iOut <- iOut.all[,iTest.rm>0,drop=FALSE] ## remove combinaisons with no common times
            attr(iOut, "index.triplet") <- stats::setNames(iOut.Mindicator[iTest.rm>0], colnames(iOut))
        }
        
        attr(iOut, "key") <- array(0, dim = rep(length(iParamVcov), 3), dimnames = list(iParamVcov,iParamVcov,iParamVcov))
        for(iCol in 1:NCOL(iOut)){
            attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol],iOut[3,iCol]] <- iCol
            attr(iOut, "key")[iOut[1,iCol],iOut[3,iCol],iOut[2,iCol]] <- iCol
            attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol],iOut[3,iCol]] <- iCol
            attr(iOut, "key")[iOut[2,iCol],iOut[3,iCol],iOut[1,iCol]] <- iCol
            attr(iOut, "key")[iOut[3,iCol],iOut[1,iCol],iOut[2,iCol]] <- iCol
            attr(iOut, "key")[iOut[3,iCol],iOut[2,iCol],iOut[1,iCol]] <- iCol
        }
        attr(iOut,"type") <- matrix(param2type[as.vector(iOut)], nrow = 3, ncol = NCOL(iOut),
                                    dimnames = dimnames(iOut))
        collapse.type <- nlme::collapse(t(attr(iOut,"type")), as.factor = FALSE)
        attr(iOut,"typetype") <- tapply(1:length(collapse.type),collapse.type,base::identity,simplify = FALSE)
        attr(iOut,"index.global") <- match(colnames(iOut), colnames(alltriplet))        
        return(iOut)
    }),Upattern$name)

    
    ## ** export
    attr(out,"global") <- alltriplet
    return(out)
}

##----------------------------------------------------------------------
### precompute.R ends here
