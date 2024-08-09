### df.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (10:34) 
## Version: 
## Last-Updated: aug  7 2024 (18:42) 
##           By: Brice Ozenne
##     Update #: 259
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
.df_analytic <- function(residuals, precision, dOmega, d2Omega, vcov, Upattern.ncluster,
                         pattern, index.cluster, name.varcoef, name.allcoef,
                         pair.meanvarcoef, pair.varcoef, indiv, REML, type.information, name.effects, robust, diag,
                         precompute){
    
    ## ** extract information
    if(is.null(precompute)){
        stop("Cannot compute degrees-of-freedom analytically when \'precompute.moments\' is set to FALSE. \n",
             "Use the function LMMstar.options to update \'precompute.moments\'. \n")
    }
    n.obs <- length(index.cluster)
    n.cluster <- length(pattern)
    n.allcoef <- length(name.allcoef)
    name.allvarcoef <- name.allcoef[name.allcoef %in% unique(unlist(name.varcoef))] ## make sure the ordering is correct
    n.allvarcoef <- length(name.allvarcoef)
    name.mucoef <- setdiff(name.allcoef, name.allvarcoef)
    n.mucoef <- length(name.mucoef)
    U.pattern <- names(dOmega)
    n.pattern <- length(U.pattern)

    ## ** prepare
    ## *** third derivative of Omega
    if(type.information == "observed" || REML){
        browser()
    }

    ## *** sets of parameters
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

        iX <- precompute$XX$pattern[[iPattern]]
                    
        ## *** mean,mean,vcov
        for(iParam in name.varcoef[[iPattern]]){ ## iParam <- name.varcoef[[iPattern]][1]
            iValue <- (as.double(iOmega %*% dOmega[[iPattern]][[iParam]] %*% iOmega) %*% iX)[as.double(precompute$XX$key)]
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
.df_numDeriv <- function(reparametrize,
                         value,  design, time, method.fit, type.information,
                         transform.sigma, transform.k, transform.rho,
                         effects, robust, 
                         precompute.moments, method.numDeriv){

    ## ** prepare vector of parameters
    param.value <- value
    param.type <- stats::setNames(design$param$type, design$param$name)
    sigma <- stats::setNames(design$param$sigma, design$param$name)
    k.x <- stats::setNames(design$param$k.x, design$param$name)
    k.y <- stats::setNames(design$param$k.y, design$param$name)
    name.allcoef <- design$param$name
    n.allcoef <- length(param.type)

    param.nameVar <- name.allcoef[param.type %in% c("sigma","k","rho")]
    param.nameMean <- name.allcoef[param.type %in% c("mu")]

    test.transform <- (transform.sigma != "none") || (transform.k != "none") || (transform.rho != "none")

    param.trans.value <- c(param.value[param.nameMean],reparametrize$p)[name.allcoef]
    name.effects <- attr(effects,"original.output")
    n.effects <- length(name.effects)

    ## ** warper for computing information
    get.element <- ifelse(robust,"vcov","information")
    FUN_element <- function(p, as.double){

        if(test.transform){ ## back-transform
            backp <- .reparametrize(p = p[param.nameVar],
                                    type = param.type[param.nameVar],
                                    sigma = sigma[param.nameVar], 
                                    k.x = k.x[param.nameVar], 
                                    k.y = k.y[param.nameVar], 
                                    Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                    transform.sigma = transform.sigma,
                                    transform.k = transform.k,
                                    transform.rho = transform.rho,
                                    transform.names = FALSE)
            p[param.nameVar] <- backp$p
        }

        iMoment <- .moments.lmm(value = p, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                logLik = FALSE, score = FALSE, information = !robust, vcov = robust, df = FALSE, indiv = FALSE, effects = effects, robust = robust,
                                trace = FALSE, precompute.moments = precompute.moments, method.numDeriv = method.numDeriv, transform.names = FALSE)

        if(as.double){
            return(as.double(iMoment[[get.element]]))
        }else{
            return(iMoment[[get.element]])
        }
    }

    ## ** variance-covariance matrix
    if(robust){
        vcov.all <- FUN_element(param.trans.value, as.double = FALSE)
    }else{
        vcov.all <- solve(FUN_element(param.trans.value, as.double = FALSE))
    }
    vcov <- vcov.all[name.effects,name.effects,drop=TRUE]

    ## ** derivative of the information using numerical derivative
    M.dElement <- numDeriv::jacobian(func = function(p){FUN_element(p, as.double = TRUE)}, x = param.trans.value, method = method.numDeriv)
    colnames(M.dElement) <- name.allcoef
    
    ## ** reshape into an array (model-based vcov: chain rule with inversion)
    A.dVcov <- array(NA, dim = c(n.effects,n.effects,n.allcoef), dimnames = list(name.effects,name.effects,name.allcoef))
    for(iParam in name.allcoef){ ## iParam <- name.allcoef[1]
        if(robust){
            A.dVcov[,,iParam] <- matrix(M.dElement[,iParam], nrow = n.allcoef, ncol = n.allcoef, dimnames = list(name.allcoef,name.allcoef))[name.effects,name.effects,drop=FALSE]
        }else{
            A.dVcov[,,iParam] <- - (vcov.all %*% matrix(M.dElement[,iParam], nrow = n.allcoef, ncol = n.allcoef) %*% vcov.all)[name.effects,name.effects,drop=FALSE]
        }
    }

    ## ** degrees-of-freedom
    df <- .df_contrast(contrast = NULL, vcov.param = vcov.all, dVcov.param = A.dVcov)
    
    ## ** export
    attr(df,"dVcov") <- A.dVcov
    return(df)
}

## * .df_contrast
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

