### df.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (10:34) 
## Version: 
## Last-Updated: okt 15 2026 (09:38) 
##           By: Brice Ozenne
##     Update #: 465
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * df.residual.lmm
##' @title Residual Degrees-of-Freedom From a Linear Mixed Model.
##' @description Estimate the residual degrees-of-freedom from a linear mixed model.
##' 
##' @param object a \code{lmm} object.
##' @param ... Passed to \code{residuals.lmm}.
##'
##' @details The residual degrees-of-freedom is estimated using the sum of squared normalized residuals.
##' 
##' @return A numeric value
##' 
##' @keywords methods
##' @export
df.residual.lmm <- function(object, ...){

    epsilonN <- stats::residuals(object, type = "normalized", ...)
    out <- as.double(crossprod(stats::na.omit(epsilonN)))
    return(out)
    
}

## * df.residual.mlmm
##' @title Residuals Degrees-of-Freedom From From Multiple Linear Mixed Model.
##' @description Combine the residuals degrees-of-freedom from group-specific linear mixed models.
##' 
##' @param object a \code{lmm} object.
##' @param ... Passed to \code{residuals.lmm}.
##'
##' @details The residual degrees-of-freedom is estimated separately for each model using the sum of squared normalized residuals.
##' 
##' @return A numeric vector with one element for each model.
##' 
##' @keywords methods
##' @export
df.residual.mlmm <- function(object, ...){

    out <- sapply(object$model,df.residual)
    return(out)
    
}

## * .df_analytic
.df_analytic <- function(X, residuals, weights, vcov, robust,
                         pattern, index.cluster, name.allcoef, pair.vcov, triplet.vcov,
                         REML, type.information, effects, precompute){

    ## ** extract information
    n.cluster <- length(pattern)
    n.mucoef <- NCOL(X)
    
    name.mucoef <- colnames(X)
    name.vcovcoef <- intersect(setdiff(name.allcoef,name.mucoef), unlist(lapply(precompute$Omega$dOmegaM1,function(iA){dimnames(iA)[[2]]})))
    n.vcovcoef <- length(name.vcovcoef)
    npair.vcov <- NROW(pair.vcov) ## number of pairs of variance parameters
    ntriplet.vcov <- NROW(triplet.vcov) ## number of triplet of variance parameters

    ## ** prepare output
    test.mean <- (n.mucoef>0) && ("mean" %in% effects)
    test.vcov <- (n.vcovcoef>0) && (("variance" %in% effects) || ("correlation" %in% effects))

    name.effects <- c(name.mucoef, name.vcovcoef)
    n.effects <- length(name.effects)
    dInformation <- array(0, dim = rep(n.mucoef + n.vcovcoef, 3), dimnames = list(name.allcoef, name.allcoef, name.allcoef))

    ## ** compute derivative of the information matrix
    if(REML && test.vcov){
        n.time <- sqrt(length(precompute$REML$X.OmegaM1.X_M1))
        M.REML.denom <- matrix(precompute$REML$X.OmegaM1.X_M1, ncol = n.time, nrow = n.time)
        REML.num1 <- do.call(cbind,apply(precompute$REML$X.dOmegaM1.X, MARGIN = 2, function(iVec){as.vector(M.REML.denom %*% matrix(iVec, ncol = n.time, nrow = n.time))}, simplify = FALSE))
        REML.num2 <- do.call(cbind,apply(precompute$REML$X.d2OmegaM1.X, MARGIN = 2, function(iVec){as.vector(M.REML.denom %*% matrix(iVec, ncol = n.time, nrow = n.time))}, simplify = FALSE))

        REML.num11 <- do.call(cbind,apply(REML.num1, MARGIN = 2, FUN = function(iCol){
            as.numeric(matrix(iCol, n.time, n.time) %*% M.REML.denom)
        }, simplify = FALSE))
        ## matrix(REML.num1[,5], n.time, n.time) - t(matrix(REML.num1[,5], n.time, n.time))
        ## matrix(REML.num11[,5], n.time, n.time) - t(matrix(REML.num11[,5], n.time, n.time))
        REML.num22 <- do.call(cbind,stats::setNames(lapply(1:NROW(pair.vcov), FUN = function(iP){
            as.numeric(matrix(precompute$REML$X.dOmegaM1.X[,pair.vcov$param1[iP]], n.time, n.time) %*% matrix(REML.num1[,pair.vcov$param2[iP]], n.time, n.time))
        }),pair.vcov$name))
        
    }
    ## ** derivative of the hessian (looping over covariance patterns)
    for (iG in 1:NROW(precompute$weights)) { ## iG <- 1

        iPattern <- precompute$weights[iG,"pattern"]
        browser()
        ## *** mean,mean,vcov
        if(test.mean & test.vcov){

            idOmegaM1 <- precompute$Omega$dOmegaM1[[iPattern]]
            iXX <- precompute$XX$pattern[[iPattern]]
            iX.dOmegaM1.X <- - t(iXX) %*% idOmegaM1 ## single covariance pattern
            ## precompute$REML$X.dOmegaM1.X ## over the whole dataset

            for(iParam in colnames(idOmegaM1)){ ## iParam <- colnames(idOmegaM1)[1]
                dInformation[name.mucoef,name.mucoef,iParam] <- dInformation[name.mucoef,name.mucoef,iParam] + iX.dOmegaM1.X[as.double(precompute$XX$key),iParam]
                dInformation[iParam,name.mucoef,name.mucoef] <- dInformation[iParam,name.mucoef,name.mucoef] + iX.dOmegaM1.X[as.double(precompute$XX$key),iParam]
                dInformation[name.mucoef,iParam,name.mucoef] <- dInformation[name.mucoef,iParam,name.mucoef] + iX.dOmegaM1.X[as.double(precompute$XX$key),iParam]
            }

        }
        ## *** mean,vcov,vcov
        if(test.mean & test.vcov & type.information=="observed"){

            id2OmegaM1 <- precompute$Omega$d2OmegaM1[[iPattern]]
            iXR <- precompute$XR[[iPattern]]
            iX.d2OmegaM1.R <- - t(iXR) %*% id2OmegaM1
            
            for(iPair in colnames(id2OmegaM1)){ ## iPair <- colnames(id2OmegaM1)[1]
                iParam <- pair.vcov[pair.vcov$name == iPair, c("param1","param2")]
                dInformation[iParam[1],iParam[2],name.mucoef] <- dInformation[iParam[1],iParam[2],name.mucoef] + iX.d2OmegaM1.R[,iPair]                
                dInformation[name.mucoef,iParam[1],iParam[2]] <- dInformation[name.mucoef,iParam[1],iParam[2]] + iX.d2OmegaM1.R[,iPair]                
                dInformation[iParam[2],name.mucoef,iParam[1]] <- dInformation[iParam[2],name.mucoef,iParam[1]] + iX.d2OmegaM1.R[,iPair]                
                if(iParam[1]!=iParam[2]){
                    dInformation[iParam[2],iParam[1],name.mucoef] <- dInformation[iParam[2],iParam[1],name.mucoef] + iX.d2OmegaM1.R[,iPair]
                    dInformation[name.mucoef,iParam[2],iParam[1]] <- dInformation[name.mucoef,iParam[2],iParam[1]] + iX.d2OmegaM1.R[,iPair]                
                    dInformation[iParam[1],name.mucoef,iParam[2]] <- dInformation[iParam[1],name.mucoef,iParam[2]] + iX.d2OmegaM1.R[,iPair]                
                }
            }
        }

        ## *** vcov,vcov,vcov (residual)
        if(test.vcov && type.information=="observed"){

            id3OmegaM1 <- precompute$Omega$d3OmegaM1[[iPattern]]
            iRR <- precompute$RR[[iPattern]]
            iR.d3OmegaM1.R <- - 0.5 * t(iRR) %*% id3OmegaM1
            
            for(iTriplet in colnames(id3OmegaM1)){ ## iTriplet <- colnames(id2OmegaM1)[1]
                iParam <- triplet.vcov[triplet.vcov$name == iTriplet, c("param1","param2","param3")]
                dInformation[iParam,iParam,iParam] <- dInformation[iParam,iParam,iParam] + iR.d3OmegaM1.R[iTriplet]                                
            }
        }

        ## *** vcov,vcov,vcov (trace)
        if(test.vcov){

            iKey <- c(stats::setNames(triplet.vcov$key1,triplet.vcov$name)[which(pair.vcov[[iPattern.vcov]])],
                      stats::setNames(triplet.vcov$key2,triplet.vcov$name)[intersect(which(pair.vcov[[iPattern.vcov]]),which(!is.na(pair.vcov$key2)))])
            
            iWeight <- precompute$weights[iG,"sum_w"]
            id3Trace <- precompute$Omega$d3Trace[[iPattern]]
            
            for(iTriplet in names(id3Trace)){ ## iTriplet <- names(id3Trace)[1]
                iParam <- triplet.vcov[triplet.vcov$name == iTriplet, c("param1","param2","param3")]
                dInformation[iParam,iParam,iParam] <- dInformation[iParam,iParam,iParam] - 0.5 * iWeight * id3Trace[iTriplet]                
            }
            
        }

        ## *** REML
        if(REML && test.vcov){
            
            REML.denom <- precompute$REML$X.OmegaM1.X_M1
            REML.num3 <- precompute$REML$X.d3OmegaM1.X[,triplet.vcov$name,drop=FALSE]
            
            REML.trace3 <- colSums(sweep(REML.num3, FUN = "*", MARGIN = 1, STATS = REML.denom))
            REML.trace2a <- colSums(REML.num1[,triplet.vcov$param3,drop=FALSE] * REML.num2[,paste(triplet.vcov$param1,triplet.vcov$param2,sep="_"),drop=FALSE])
            REML.trace2b <- colSums(REML.num2[,paste(triplet.vcov$param1,triplet.vcov$param3,sep="_"),drop=FALSE] * REML.num1[,triplet.vcov$param2,drop=FALSE])
            REML.trace2c <- colSums(REML.num1[,triplet.vcov$param1,drop=FALSE] * REML.num2[,paste(triplet.vcov$param2,triplet.vcov$param3,sep="_"),drop=FALSE])
            ## WARNING: tr(AB) = sum(A B^T) which may not be equal to sum(A B) when B is not symmetric like here B= M_1 M_2 which generally differs from M_2 M_1
            REML.trace1a <- colSums(REML.num11[,triplet.vcov$param3,drop=FALSE]  * REML.num22[,paste(triplet.vcov$param1,triplet.vcov$param2,sep="_"),drop=FALSE])
            ## GS3a <- apply(triplet.vcov[,c("param1","param2","param3")], MARGIN = 1, FUN = function(iParam){
            ##     tr(matrix(REML.num1[,iParam[3]], n.time, n.time) %*% matrix(REML.num1[,iParam[1]], n.time, n.time) %*% matrix(REML.num1[,iParam[2]], n.time, n.time))
            ## })
            ## range(REML.trace1a - GS3a)
            REML.trace1b <- colSums(REML.num11[,triplet.vcov$param2,drop=FALSE]  * REML.num22[,paste(triplet.vcov$param1,triplet.vcov$param3,sep="_"),drop=FALSE])
            ## GS3b <- apply(triplet.vcov[,c("param1","param2","param3")], MARGIN = 1, FUN = function(iParam){
            ##     tr(matrix(REML.num1[,iParam[1]], n.time, n.time) %*% matrix(REML.num1[,iParam[3]], n.time, n.time) %*% matrix(REML.num1[,iParam[2]], n.time, n.time))
            ## })
            ## GS3c <- apply(triplet.vcov[,c("param1","param2","param3")], MARGIN = 1, FUN = function(iParam){
            ##     tr(matrix(REML.num1[,iParam[2]], n.time, n.time) %*% matrix(REML.num1[,iParam[1]], n.time, n.time) %*% matrix(REML.num1[,iParam[3]], n.time, n.time))
            ## })
            ## range(REML.trace1b - GS3b)
            ## range(GS3c - GS3b)
            REML.trace <- stats::setNames(REML.trace1a + REML.trace1b + REML.trace2a + REML.trace2b + REML.trace2c + REML.trace3, triplet.vcov$name)
            
            for(iTriplet in triplet.vcov$name){ ## iTriplet <- triplet.vcov$name[1]
                iParam <- triplet.vcov[triplet.vcov$name == iTriplet, c("param1","param2","param3")]
                dInformation[iParam,iParam,iParam] <- dInformation[iParam,iParam,iParam] - 0.5 * REML.trace[iTriplet]
            }
        }
        
    }

    
    ## ** derivative of the variance covariance matrix
    n.effects <- length(name.effects)
    vcov.effects <- vcov[name.effects,name.effects,drop=FALSE]
    ## print(dInformation)
    
    A.dVcov <- array(0, dim = c(n.effects,n.effects,n.allcoef), dimnames = list(name.effects,name.effects,name.allcoef))
    for(iParam in 1:n.allcoef){ ## iParam <- 1
        A.dVcov[,,iParam] <- vcov.effects %*% dInformation[name.effects,name.effects,iParam] %*% vcov.effects
    }

    ## solve(crossprod(model.matrix(e.lmm, effects = "mean")))
    ## 4*coef(e.lmm)["sigma"]^2/stats::nobs(e.lmm)[1]
    ## ** degrees-of-freedom
    if(diag){
        df <- stats::setNames(sapply(1:n.effects, function(iP){
            2 * vcov.effects[iP,iP]^2 / (A.dVcov[iP,iP,] %*% vcov %*% A.dVcov[iP,iP,])
        }), name.effects)
    }else{
        df <- matrix(NA, nrow = n.effects, ncol = n.effects, dimnames = list(name.effects, name.effects))
        for(iParam in 1:n.effects){
            for(iiParam in 1:iParam){
                df[iParam,iiParam] <- 2 * vcov.effects[iParam,iiParam]^2 / (A.dVcov[iParam,iiParam,] %*% vcov %*% A.dVcov[iiParam,iParam,])
                if(iParam != iiParam){
                    df[iiParam,iParam] <- df[iParam,iiParam]
                }
            }
        }
    }
    
    ## ** export
    attr(df,"dVcov") <- A.dVcov
    return(df)    
}
## * .df_numDeriv
.df_numDeriv <- function(param, design, time, method.fit, type.information,
                         transform.sigma, transform.k, transform.rho,
                         effects, robust, vcov, method.numDeriv){  ## if robust is 1 then robust s.e. are computed but df are relative to model-based s.e.

    ## ** prepare vector of parameters
    name.allcoef <- design$param$name
    n.allcoef <- length(name.allcoef)
    name.effects <- names(param)
    n.effects <- length(name.effects)

    ## ** warper for computing information
    get.element <- ifelse(robust,"vcov","information")
    FUN_element <- function(p, as.double){
        
        iMoment <- .moments.lmm(value = p, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                logLik = FALSE, score = FALSE, information = (robust!=2), vcov = (robust==2), df = FALSE, indiv = FALSE, effects = effects, robust = robust,
                                trace = FALSE, method.numDeriv = method.numDeriv, transform.names = FALSE)

        if(as.double){
            return(as.double(iMoment[[get.element]]))
        }else{
            return(iMoment[[get.element]])
        }
    }

    ## ** variance-covariance matrix
    if(is.null(vcov) || (robust == 1)){ ## should be recomputed when robust equals 1 as df are relative to model-based s.e.
        if(robust){ 
            vcov <- FUN_element(param, as.double = FALSE)
        }else{
            vcov <- solve(FUN_element(param, as.double = FALSE))
        }
    }
    

    ## ** derivative of the information using numerical derivative
    M.dElement <- numDeriv::jacobian(func = function(p){FUN_element(p, as.double = TRUE)}, x = param, method = method.numDeriv)
    colnames(M.dElement) <- name.allcoef
    
    ## ** reshape into an array (model-based vcov: chain rule with inversion)
    A.dVcov <- array(NA, dim = c(n.effects,n.effects,n.allcoef), dimnames = list(name.effects,name.effects,name.allcoef))
    for(iParam in name.allcoef){ ## iParam <- name.allcoef[1]
        if(robust){
            A.dVcov[,,iParam] <- matrix(M.dElement[,iParam], nrow = n.allcoef, ncol = n.allcoef, dimnames = list(name.allcoef,name.allcoef))[name.effects,name.effects,drop=FALSE]
        }else{
            A.dVcov[,,iParam] <- - (vcov %*% matrix(M.dElement[,iParam], nrow = n.allcoef, ncol = n.allcoef) %*% vcov)[name.effects,name.effects,drop=FALSE]
        }
    }

    ## ** degrees-of-freedom
    df <- .df_contrast(contrast = NULL, vcov.param = vcov, dVcov.param = A.dVcov)
    
    ## ** export
    attr(df,"dVcov") <- A.dVcov
    return(df)
}

## * helpers
## ** .df_contrast
##' @description Evaluate degrees-of-freedom for a linear combination of parameters
##' @param contrast [matrix] contrast matrix (n,p)
##' @param vcov.param [matrix] matrix (p,p)
##' @param dVcov.param [array] array (p,p,p) or (n,n,p)
##' @param return.vcov [logical] 
##' @noRd
.df_contrast <- function(contrast, vcov.param, dVcov.param, return.vcov = FALSE){

    ## ** normalize user input
    ## *** p: name and number of parameters
    n.param <- NCOL(vcov.param)
    name.param <- colnames(vcov.param)

    ## *** n: name and number of linear combinations
    if(is.null(contrast)){        
        n.contrast <- dim(dVcov.param)[1]
        name.contrast <- dimnames(dVcov.param)[[1]]
    }else if(!is.matrix(contrast) && is.vector(contrast)){
        contrast <- rbind(contrast)
        n.contrast <- 1
        name.contrast <- NULL
    }else{
        n.contrast <- NROW(contrast)
        name.contrast <- rownames(contrast)
    }
    

    ## ** variance of the contrast (Sigma_beta)
    if(is.null(contrast)){
        sigma2.contrast <- diag(vcov.param)[name.contrast]
    }else{
        sigma2.contrast <- rowSums((contrast %*% vcov.param[colnames(contrast),colnames(contrast),drop=FALSE]) * contrast)
        ## diag(contrast %*% vcov.param[colnames(contrast),colnames(contrast),drop=FALSE] %*% t(contrast))
    }

    ## ** Gradient of the variance of the contrast w.r.t. the original parameters (dSigma_\beta/d\theta)
    if(is.null(contrast)){
        Mpair_dVcov.beta <- do.call(rbind,lapply(1:n.contrast, function(iContrast){dVcov.param[iContrast,iContrast,]}))
    }else{
        index.red <- which(colSums(contrast!=0)>0)
        contrast_red <- contrast[,index.red,drop=FALSE]
        Mpair_dVcov.beta <- do.call(cbind,lapply(1:n.param, function(iParam){rowSums((contrast_red %*% dVcov.param[index.red,index.red,iParam]) * contrast_red)}))
    }

    ## ** Satterthwaite approximation of the degrees-of-freedom
    ## delta method on \beta = C\theta:  Var[Sigma_\beta] = [dSigma_\beta/d\theta] [Sigma_\theta] [dSigma_\beta/dtheta]'
    denum <- rowSums((Mpair_dVcov.beta %*% vcov.param) * Mpair_dVcov.beta)
    ## same as diag(Mpair_dVcov.beta %*% vcov.param %*% t(Mpair_dVcov.beta))
    out <- 2*sigma2.contrast^2 / denum
    out[denum==0] <- Inf

    ## ** export
    if(is.null(name.contrast)){
        names(out) <- name.contrast
    }
    if(return.vcov){
        attr(out,"vcov") <- Mpair_dVcov.beta %*% vcov.param %*% t(Mpair_dVcov.beta)
    }
    return(out)
}

## ** .unorderedTriplet
## form all combinations
## .unorderedTriplet(1:5)
.unorderedTriplet <- function(x, distinct = FALSE){
    n <- length(x)
    ls <- lapply(1:n, function(k){do.call(cbind,lapply(k:n, function(l){rbind(x[k], x[l], x[l:n])}))})
    out <- do.call(cbind,ls)
    if(distinct){
        out <- out[,apply(out,2,function(iCol){sum(duplicated(iCol))<2}),drop=FALSE]
    }
    return(out)
}

## ** .meanCovTriplet
.meanCovTriplet <- function(structure, X, index.cluster, index.clusterTime){
    
    out <- list()
    n.pattern <- NROW(structure$Upattern)
    activeparam <- c(colnames(X),structure$param[is.na(structure$param$constraint),"name"])
    
    ## *** pairs of covariance parameters
    ## collect for each pattern the vector of parameters for a given cell in the residual variance-covariance matrix
    ls.Mparam.vcovvcov <- pairParamVcovVcov(structure)
    
    ## unique combinations of variance-covariance parameters across patterns
    Mparam.vcovvcov <- unique(do.call(rbind,ls.Mparam.vcovvcov))
    
    ## table with pairs for each pattern
    structure$pair.vcovvcov <- structure$pair.vcovvcov

    ## *** pairs of mean-covariance parameters
    ls.Mparam.meanvcov <- lapply(structure$Upattern$name, function(iName){
        iCluster <- structure$Upattern[structure$Upattern$name==iName,"index.cluster"][[1]]
        ## collect, for each pattern, the vector of mean parameters
        iX.mean <- X[unlist(index.cluster[iCluster]),,drop=FALSE] != 0
        ## collect, for each pattern, the vector of variance parameters
        iX.var <- structure$var$X[unlist(index.cluster[iCluster]),,drop=FALSE] != 0
        if(is.null(structure$cor)){
            iUX.meanvar <- unique(cbind(iX.mean,iX.var))
            iOut <- do.call(rbind,apply(iUX.meanvar, MARGIN = 1, function(iRow){c(names(iRow)[which(iRow)],rep("",sum(iRow==FALSE)))}, simplify = FALSE))
        }else{
            iNtime <- structure$Upattern[structure$Upattern$name == iName,"n.time"]
            iName.var <- structure$Upattern[structure$Upattern$name == iName,"var"]
            iName.cor <- structure$Upattern[structure$Upattern$name == iName,"cor"]
            iUX.mean <- unique(cbind(XXindexXX = unlist(index.clusterTime[iCluster]),iX.mean))
            
            iUX.varcor <- data.frame(XXindex1XX=rep(1:iNtime, times = iNtime),
                                     XXindex2XX=rep(1:iNtime, each = iNtime),
                                     apply(structure$var$Xpattern[[iName.var]], MARGIN = 2, FUN = rep, times = iNtime), 
                                     apply(structure$var$Xpattern[[iName.var]], MARGIN = 2, FUN = rep, each = iNtime), 
                                     do.call(cbind, apply(structure$cor$Xpattern[[iName.cor]], MARGIN = 3, FUN = as.vector, simplify = FALSE))
                                     )
            iOut <- cbind(iUX.varcor, matrix("", nrow = NROW(iUX.varcor), ncol = NCOL(iX.mean), dimnames = list(NULL, colnames(iX.mean))))
            for(iR in 1:NROW(iOut)){ ## iR <- 1
                iTest <- colSums(iUX.mean[iUX.mean[,"XXindexXX"] %in% iUX.varcor[iR,c("XXindex1XX","XXindex1XX")],-1,drop=FALSE])
                iOut[iR,names(iTest[iTest>0])] <- names(iTest[iTest>0])
            }
            iOut$XXindex1XX <- NULL
            iOut$XXindex2XX <- NULL
        }
        return(as.matrix(iOut))
    })
    names(ls.Mparam.meanvcov) <- structure$Upattern$name
    ## unique combinations of mean and covariance parameters across patterns
    Mparam.meanvcov <- unique(do.call(rbind,ls.Mparam.meanvcov))

    ## *** triplets of covariance parameters
    ## form all triplets in each pattern
    ls.triplet.vcovvcov <- apply(Mparam.vcovvcov, MARGIN = 1, FUN = function(iParam){
        unorderedTriplet(intersect(activeparam,iParam)) ## intersect will remove duplicates (e.g. sigma, one, sigma --> sigma, one) and keep a constant ordering
    }, simplify = FALSE)

    ## unique triplets across all patterns
    Utriplet.vcovvcov <- unique(t(do.call(cbind,ls.triplet.vcovvcov)))
    out$vcov3 <- data.frame(paste(Utriplet.vcovvcov[,1],Utriplet.vcovvcov[,2],Utriplet.vcovvcov[,3],sep = "_"),
                                    Utriplet.vcovvcov,
                                    matrix(NA, nrow = NROW(Utriplet.vcovvcov), ncol = n.pattern))
    names(out$vcov3) <- c("name","param1","param2","param3",structure$Upattern$name)
    ## find which triplet appears in which pattern
    for(iPattern in structure$Upattern$name){ ## iP <- structure$Upattern$name[1]
        ## put parameter from the pattern first
        iParam.pattern <- structure$Upattern$param[[which(structure$Upattern$name == iPattern)]]
        iParam.active <- c(iParam.pattern,setdiff(activeparam,iParam.pattern))
        ## convert to numeric
        iUtriplet.vcov3.num <- matrix(as.numeric(factor(Utriplet.vcovvcov, levels = iParam.active)), nrow = NROW(Utriplet.vcovvcov), ncol = NCOL(Utriplet.vcovvcov))
        ## check all parameters belong to pattern
        out$vcov3[,iPattern] <- rowSums(iUtriplet.vcov3.num <= length(iParam.pattern))==3
    }

    ## *** triplets of mean-covariance parameters
    ## form all triplets in each pattern
    ls.triplet.meanvcov <- apply(Mparam.meanvcov, MARGIN = 1, FUN = function(iParam){
        unorderedTriplet(intersect(activeparam,iParam)) ## intersect will remove duplicates (e.g. sigma, one, sigma --> sigma, one) and keep a constant ordering
    }, simplify = FALSE)
    ## unique triplets across all patterns
    Utriplet.meanvcov <- unique(t(do.call(cbind,ls.triplet.meanvcov)))
    Utriplet.meanvcov.num <- matrix(as.numeric(factor(Utriplet.meanvcov, levels = activeparam)), nrow = NROW(Utriplet.meanvcov), ncol = NCOL(Utriplet.meanvcov))
    Utriplet.meanvcov.nMu <- rowSums(Utriplet.meanvcov.num <= sum(activeparam %in% colnames(X)))
    ## mean covariance covariance
    Utriplet.mean1vcov2 <- Utriplet.meanvcov[Utriplet.meanvcov.nMu==1,,drop=FALSE]
    out$mean1vcov2 <- data.frame(paste(Utriplet.mean1vcov2[,1],Utriplet.mean1vcov2[,2],Utriplet.mean1vcov2[,3],sep = "_"),
                                        Utriplet.mean1vcov2,
                                        matrix(NA, nrow = NROW(Utriplet.mean1vcov2), ncol = n.pattern))
    names(out$mean1vcov2) <- c("name","param1","param2","param3",structure$Upattern$name)
    ## mean mean covariance
    Utriplet.mean2vcov1 <- Utriplet.meanvcov[Utriplet.meanvcov.nMu==2,,drop=FALSE]
    out$mean2vcov1 <- data.frame(paste(Utriplet.mean2vcov1[,1],Utriplet.mean2vcov1[,2],Utriplet.mean2vcov1[,3],sep = "_"),
                                         Utriplet.mean2vcov1,
                                         matrix(NA, nrow = NROW(Utriplet.mean2vcov1), ncol = n.pattern))
    names(out$mean2vcov1) <- c("name","param1","param2","param3",structure$Upattern$name)

    for(iPattern in structure$Upattern$name){ ## iP <- 1
        ## put parameter from the pattern first
        iParam.pattern <- intersect(activeparam, as.vector(ls.Mparam.meanvcov[[iPattern]]))
        iParam.active <- c(iParam.pattern,setdiff(activeparam,iParam.pattern))
        ## convert to numeric
        iUtriplet.mean1vcov2.num <- matrix(as.numeric(factor(Utriplet.mean1vcov2, levels = iParam.active)),
                                          nrow = NROW(Utriplet.mean1vcov2), ncol = NCOL(Utriplet.meanvcov))
        iUtriplet.mean2vcov1.num <- matrix(as.numeric(factor(Utriplet.mean2vcov1, levels = iParam.active)),
                                          nrow = NROW(Utriplet.mean2vcov1), ncol = NCOL(Utriplet.meanvcov))
        ## check all parameters belong to pattern
        out$mean1vcov2[,iPattern] <- rowSums(iUtriplet.mean1vcov2.num <= length(iParam.pattern))==3
        out$mean2vcov1[,iPattern] <- rowSums(iUtriplet.mean2vcov1.num <= length(iParam.pattern))==3
    }

    ## *** export
    return(out)
}

    
##----------------------------------------------------------------------
### df.R ends here

