### df.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (10:34) 
## Version: 
## Last-Updated: nov 12 2021 (13:49) 
##           By: Brice Ozenne
##     Update #: 132
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .df_analytic
.df_analytic <- function(residuals, precision, dOmega, d2Omega, vcov, Upattern.ncluster,
                         index.variance, time.variance, index.cluster, name.varcoef, name.allcoef,
                         pair.meanvarcoef, pair.varcoef, indiv, REML, type.information, name.effects, robust, diag,
                         precompute){
    
    if(is.null(precompute)){
        stop("Cannot compute degrees of freedom analytically when \'precompute.moments\' is set to FALSE. \n",
             "Use the function LMMstar.options to update \'precompute.moments\'. \n")
    }
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    n.varcoef <- lapply(name.varcoef, length)
    n.allcoef <- length(name.allcoef)
    name.allvarcoef <- name.allcoef[name.allcoef %in% unique(unlist(name.varcoef))] ## make sure the ordering is correct
    n.allvarcoef <- length(name.allvarcoef)
    name.mucoef <- setdiff(name.allcoef, name.allvarcoef)
    n.mucoef <- length(name.mucoef)
    U.pattern <- names(dOmega)
    n.pattern <- length(U.pattern)

    if(type.information == "expected"){
        triplet.meanvarcoef <- stats::setNames(lapply(U.pattern, function(iPattern){ ## iPattern <- "1.1"
            iGrid <- expand.grid(1:NCOL(pair.varcoef[[iPattern]]),name.allcoef)
            iTriplet <- rbind(pair.varcoef[[iPattern]][1,iGrid[,1]],pair.varcoef[[iPattern]][2,iGrid[,1]],as.character(iGrid[,2]))
            return(iTriplet[,apply(iTriplet, 2, function(iT){sum(iT %in% name.mucoef)<= 1}),drop=FALSE])
        }), U.pattern)
        triplet.varcoef <- stats::setNames(lapply(U.pattern, function(iPattern){ ## iPattern <- "1.1"
            iGrid <- expand.grid(1:NCOL(pair.varcoef[[iPattern]]),name.varcoef[[iPattern]])
            iTriplet <- rbind(pair.varcoef[[iPattern]][1,iGrid[,1]],pair.varcoef[[iPattern]][2,iGrid[,1]],as.character(iGrid[,2]))
            return(iTriplet[,apply(iTriplet, 2, function(iT){sum(iT %in% name.mucoef)<= 1}),drop=FALSE])
        }), U.pattern)
    }else if(type.information == "observed"){
        triplet.meanvarcoef <- stats::setNames(lapply(U.pattern, function(iPattern){ ## iPattern <- "1.1"
            iParam <- c(name.mucoef,name.varcoef[[iPattern]])
            iTriplet <- .unorderedTriplet(iParam)
            return(iTriplet[,apply(iTriplet, 2, function(iT){sum(iT %in% name.mucoef)<= 1}),drop=FALSE])
        }), U.pattern)
        triplet.varcoef <- lapply(name.varcoef, .unorderedTriplet)
    }
    
    

    dhessian <- array(0, dim = rep(n.allcoef,3), dimnames = list(name.allcoef, name.allcoef, name.allcoef))
    
    for (iPattern in U.pattern) { ## iPattern <- name.pattern[1]
        iOmega <- precision[[iPattern]]
        iTime <- NCOL(iOmega)
        iTime2 <- length(iOmega)
        iName.varcoef <- name.varcoef[[iPattern]]
        iN.varcoef <- length(iName.varcoef)

        iX <- matrix(unlist(precompute$XX$pattern[[iPattern]]), nrow = iTime2, ncol = dim(precompute$XX$pattern[[iPattern]])[3], byrow = FALSE)
                    
        ## *** mean,mean,vcov
        for(iParam in name.varcoef[[iPattern]]){ ## iParam <- name.varcoef[[iPattern]][1]
            iValue <- (as.double(iOmega %*% dOmega[[iPattern]][[iParam]] %*% iOmega)  %*% iX)[as.double(precompute$XX$key)]
            dhessian[name.mucoef,name.mucoef,iParam] <- dhessian[name.mucoef,name.mucoef,iParam] + iValue
        }
        ## *** mean,vcov,vcov
        ## for(iParam in name.varcoef[[iPattern]]){ ## iParam <- name.varcoef[[iPattern]][1]
        ## }
        
        ## *** vcov,vcov,vcov
        for(iTriplet in 1:NCOL(triplet.varcoef[[iPattern]])){ ## iTriplet <- 1
            iParam <- triplet.varcoef[[iPattern]][,iTriplet]
            termA <- precision[[iPattern]] %*% dOmega[[iPattern]][[iParam[3]]] %*% precision[[iPattern]] %*% dOmega[[iPattern]][[iParam[2]]] %*% precision[[iPattern]] %*% dOmega[[iPattern]][[iParam[1]]]
            index.pairB1 <- unique(c(which((pair.varcoef[[iPattern]][1,]==iParam[[2]])*(pair.varcoef[[iPattern]][2,]==iParam[[3]])>0),
                                     which((pair.varcoef[[iPattern]][1,]==iParam[[3]])*(pair.varcoef[[iPattern]][2,]==iParam[[2]])>0)))
            index.pairB2 <- unique(c(which((pair.varcoef[[iPattern]][1,]==iParam[[1]])*(pair.varcoef[[iPattern]][2,]==iParam[[3]])>0),
                                     which((pair.varcoef[[iPattern]][1,]==iParam[[3]])*(pair.varcoef[[iPattern]][2,]==iParam[[1]])>0)))
            termB <- precision[[iPattern]] %*% d2Omega[[iPattern]][[index.pairB1]] %*% precision[[iPattern]] %*% dOmega[[iPattern]][[iParam[1]]]
            termC <- precision[[iPattern]] %*% dOmega[[iPattern]][[iParam[2]]] %*% precision[[iPattern]] %*% d2Omega[[iPattern]][[index.pairB2]]
            iValue <- -0.5 * Upattern.ncluster[iPattern] * tr(-2*termA + termB + termC)
            ## if(all(iParam==c("rho","rho","rho"))){browser()}
            dhessian[iParam[1],iParam[2],iParam[3]] <- dhessian[iParam[1],iParam[2],iParam[3]] + iValue ## a b c
            if(iParam[1]!=iParam[2]){ 
                dhessian[iParam[2],iParam[1],iParam[3]] <- dhessian[iParam[2],iParam[1],iParam[3]] + iValue ## b a c
            }
            if(iParam[1]!=iParam[3] && type.information == "observed"){
                dhessian[iParam[3],iParam[2],iParam[1]] <- dhessian[iParam[3],iParam[2],iParam[1]] + iValue ## c b a
            }
            if(iParam[2]!=iParam[3] && type.information == "observed"){
                dhessian[iParam[1],iParam[3],iParam[2]] <- dhessian[iParam[1],iParam[3],iParam[2]] + iValue ## a c b
            }
            if(iParam[1]!=iParam[2] && iParam[2]!=iParam[3] && type.information == "observed"){
                dhessian[iParam[2],iParam[3],iParam[1]] <- dhessian[iParam[2],iParam[3],iParam[1]] + iValue ## b c a
                dhessian[iParam[3],iParam[1],iParam[2]] <- dhessian[iParam[3],iParam[1],iParam[2]] + iValue ## c a b
            }
        }
    }
    ## ** derivative of the variance covariance matrix
    n.effects <- length(name.effects)
    vcov.effects <- vcov[name.effects,name.effects,drop=FALSE]
    ## print(dhessian)
    
    A.dVcov <- array(0, dim = c(n.effects,n.effects,n.allcoef), dimnames = list(name.effects,name.effects,name.allcoef))
    for(iParam in 1:n.allcoef){ ## iParam <- 1
        A.dVcov[,,iParam] <- vcov.effects %*% dhessian[name.effects,name.effects,iParam] %*% vcov.effects
    }

    ## solve(crossprod(model.matrix(e.lmm, effects = "mean")))
    ## 4*coef(e.lmm)["sigma"]^2/stats::nobs(e.lmm)[1]
    ## ** degrees of freedom
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
.df_numDeriv <- function(value, reparametrize,
                         design, time, method.fit, type.information,
                         transform.sigma, transform.k, transform.rho, effects, robust, diag,
                         precompute.moments = precompute.moments, method.numDeriv = method.numDeriv){

    ## ** prepare vector of parameters
    param.value <- value
    param.type <- design$param$type
    param.strata <- design$param$strata
    name.allcoef <- names(param.type)
    n.allcoef <- length(param.type)

    param.nameVar <- name.allcoef[param.type %in% c("sigma","k","rho")]
    param.nameMean <- name.allcoef[param.type %in% c("mu")]

    test.transform <- (transform.sigma != "none") || (transform.k != "none") || (transform.rho != "none")

    param.trans.value <- c(param.value[param.nameMean],reparametrize$p)[name.allcoef]

    ## ** warper for computing information
    FUN_information <- function(p, as.double){

        if(test.transform){ ## back-transform
            backp <- .reparametrize(p = p[param.nameVar],
                                    type = param.type[param.nameVar], strata = param.strata[param.nameVar], 
                                    Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                    transform.sigma = transform.sigma,
                                    transform.k = transform.k,
                                    transform.rho = transform.rho,
                                    transform.names = FALSE)
            p[param.nameVar] <- backp$p
        }

        iMoment <- .moments.lmm(value = p, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("mean","variance","correlation"), robust = robust,
                                trace = FALSE, precompute.moments = precompute.moments, method.numDeriv = method.numDeriv, transform.names = FALSE)

        if(as.double){
            return(as.double(iMoment$information))
        }else{
            return(iMoment$information)
        }
    }

    ## ** variance-covariance matrix
    ## full variance covariance matrix
    info <- FUN_information(param.trans.value, as.double = FALSE)
    vcov <- solve(info)

    ## ** derivative of the information using numerical derivative
    if(type.information == "observed"){
        M.dInfo <- numDeriv::jacobian(func = function(p){FUN_information(p, as.double = TRUE)}, x = param.trans.value, method = method.numDeriv)
        colnames(M.dInfo) <- name.allcoef
    }else{
        M.dInfo <- numDeriv::jacobian(func = function(p){FUN_information(c(param.value[param.nameMean],p)[name.allcoef], as.double = TRUE)}, x = param.trans.value[param.nameVar], method = method.numDeriv)
        colnames(M.dInfo) <- param.nameVar
    }
    ## A.print <- array(0, dim = c(n.allcoef,n.allcoef,n.allcoef), dimnames = list(name.allcoef,name.allcoef,name.allcoef))
    ## for(iParam in 1:NCOL(M.dInfo)){ ## iParam <- 1
    ##     iName <- colnames(M.dInfo)[iParam]
    ##     A.print[,,iName] <- - matrix(M.dInfo[,iName], nrow = n.allcoef, ncol = n.allcoef, dimnames = list(name.allcoef,name.allcoef))
    ## }

    ## ** derivative of the variance covariance matrix
    name.effects <- attr(effects,"original.names")
    n.effects <- length(name.effects)
    A.dVcov <- array(0, dim = c(n.effects,n.effects,n.allcoef), dimnames = list(name.effects,name.effects,name.allcoef))
    for(iParam in 1:NCOL(M.dInfo)){ ## iParam <- 1
        iName <- colnames(M.dInfo)[iParam]
        A.dVcov[,,iName] <- - (vcov %*% matrix(M.dInfo[,iName], nrow = n.allcoef, ncol = n.allcoef) %*% vcov)[name.effects,name.effects,drop=FALSE]
    }

    ## solve(crossprod(model.matrix(e.lmm, effects = "mean")))
    ## 4*coef(e.lmm)["sigma"]^2/stats::nobs(e.lmm)[1]
    ## ** degrees of freedom
    if(diag){
        df <- stats::setNames(sapply(name.effects, function(iP){
            2 * vcov[iP,iP]^2 / (A.dVcov[iP,iP,] %*% vcov %*% A.dVcov[iP,iP,])
        }), name.effects)
    }else{
        df <- matrix(NA, nrow = n.effects, ncol = n.effects, dimnames = list(name.effects, name.effects))
        for(iParam in name.effects){
            for(iiParam in name.effects[1:which(iParam==name.effects)]){
                df[iParam,iiParam] <- 2 * vcov[iParam,iiParam]^2 / (A.dVcov[iParam,iiParam,] %*% vcov %*% A.dVcov[iiParam,iParam,])
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

## * .unorderedTriplet
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

##----------------------------------------------------------------------
### df.R ends here
