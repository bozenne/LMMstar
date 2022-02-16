### structure-initialization.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:20) 
## Version: 
## Last-Updated: feb 16 2022 (17:28) 
##           By: Brice Ozenne
##     Update #: 147
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * initialization
##' @title Initialize Variance-Covariance Structure
##' @description Initialize the parameters of the variance-covariance structure using residual variance and correlations.
##' @noRd
##'
##' @param structure [structure]
##' @param residuals [vector] vector of residuals.
##' @param Xmean [matrix] design matrix for the mean effects used to estimate the residual degrees of freedom for variance calculation.
##'
##' @keywords internal
##' 
##' @examples
##' data(gastricbypassW, package = "LMMstar")
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$gender <- c("M","F")[as.numeric(gastricbypassL$id) %% 2+1]
##' dd <- gastricbypassL[!duplicated(gastricbypassL[,c("time","gender")]),]
##'
##' eDD.lm <- lm(weight ~ visit, data = dd)
##' eGas.lm <- lm(weight ~ visit*gender, data = gastricbypassL)
##' 
##' ## independence
##' Sid1 <- .skeleton(IND(~1, var.time = "time"), data = dd)
##' Sid4 <- .skeleton(IND(~1|id, var.time = "time"), data = dd)
##' Sdiag1 <- .skeleton(IND(~visit), data = dd)
##' Sdiag4 <- .skeleton(IND(~visit|id), data = dd)
##' Sdiag24 <- .skeleton(IND(~visit+gender|id, var.time = "time"), data = gastricbypassL)
##'
##' .initialize(Sid1, residuals = residuals(eDD.lm))
##' ## sd(residuals(eDD.lm))
##' .initialize(Sdiag4, residuals = residuals(eDD.lm))
##' ## tapply(residuals(eDD.lm),dd$visit,sd)
##' .initialize(Sdiag24, residuals = residuals(eGas.lm))
##' ## tapply(residuals(eGas.lm),interaction(gastricbypassL[,c("visit","gender")]),sd)
##' 
##' ## compound symmetry
##' Scs4 <- .skeleton(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' Scs24 <- .skeleton(CS(gender~time|id), data = gastricbypassL)
##' 
##' .initialize(Scs4, residuals = residuals(eGas.lm))
##' ## cor(gastricbypassW[,c("weight1","weight2","weight3","weight4")])
##' .initialize(Scs24, residuals = residuals(eGas.lm))
##' 
##' ## unstructured
##' Sun4 <- .skeleton(UN(~visit|id), data = gastricbypassL)
##' Sun24 <- .skeleton(UN(gender~visit|id), data = gastricbypassL)
##' 
##' .initialize(Sun4, residuals = residuals(eGas.lm))
##' .initialize(Sun24, residuals = residuals(eGas.lm))
`.initialize` <-
    function(object, residuals, Xmean) UseMethod(".initialize")

## * initialization.IND
.initialize.IND <- function(object, residuals, Xmean){

    param.type <- stats::setNames(object$param$type,object$param$name)
    param.strata <- stats::setNames(object$param$strata,object$param$name)
    Upattern.name <- object$X$Upattern$name
    cluster.pattern <- object$X$cluster.pattern

    ## combine all residuals and all design matrices
    M.res <- do.call(rbind,lapply(object$X$var, function(iPattern){ ## iPattern <- object$X$var[[1]]
        X.iPattern <- iPattern
        cluster.iPattern <- attr(iPattern,"index.cluster")
        obs.iPattern <- unlist(attr(iPattern,"index.obs"))
        attr(X.iPattern,"indicator.param") <- NULL
        attr(X.iPattern,"index.cluster") <- NULL
        attr(X.iPattern,"index.obs") <- NULL
        iOut <- cbind(index = obs.iPattern, residuals = residuals[obs.iPattern], do.call(rbind,rep(list(X.iPattern),length(cluster.iPattern))))
        return(iOut)
    }))

    ## extract information
    epsilon2 <- M.res[,"residuals"]^2
    X <- M.res[,-(1:2),drop=FALSE]
    paramVar.type <- param.type[colnames(X)]
    paramVar.strata <- param.strata[colnames(X)]
    n.strata <- length(unique(paramVar.strata))
    n.obs <- NROW(X)

    ## small sample correction (inflate residuals)
    if(!is.null(Xmean)){
        ## n - df
        M.indexsigma <- do.call(rbind,lapply(object$X$var, function(iPattern){ ## iPattern <- object$X$var[[2]]
            iUX <- interaction(as.data.frame(iPattern), drop = TRUE)
            iDW <- data.frame(UX = 1:NROW(iUX), X = iUX, do.call(cbind,attr(iPattern, "index.obs")))
            iDL <- stats::reshape(iDW, direction = "long", idvar = c("UX","X"), varying = setdiff(names(iDW),c("UX","X")), v.names = "index")
            return(iDL[,c("X","index")])
        }))
        ## vec.hat <- diag(Xmean %*% solve(t(Xmean) %*% Xmean) %*% t(Xmean))
        vec.hat <- rowSums(Xmean %*% solve(t(Xmean) %*% Xmean) * Xmean)
        p <- tapply(M.indexsigma$index,M.indexsigma$X,function(iIndex){sum(vec.hat[iIndex])})
        n.UX <- table(M.indexsigma$X)
        ## range(M.res[,"index"] - M.indexsigma[,"index"])
        epsilon2.ssc <- epsilon2 * (n.UX/(n.UX-p))[M.indexsigma$X]
    }else{
        epsilon2.ssc <- epsilon2
    }

    ## fit
    e.res <- stats::lm.fit(y=epsilon2.ssc,x=X) 
    if(all(paramVar.type=="sigma")){
        out <- sqrt(e.res$coef)
    }else{
        ## try to move from additive to full interaction model
        ls.Z <- lapply(1:n.strata, function(iStrata){ ## iStrata <- 1
            iParamVar.type <- paramVar.type[paramVar.strata == iStrata]
            iX <- X[,paramVar.strata==iStrata,drop=FALSE]
            if(any("k" %in% iParamVar.type)){
                iX[,iParamVar.type=="sigma"] <- iX[,iParamVar.type=="sigma"] - rowSums(iX[,iParamVar.type!="sigma",drop=FALSE])
            }
            return(iX)
        })
        
        Z <- do.call(cbind,ls.Z)[,colnames(X)]
        eTest.res <- stats::lm.fit(y=epsilon2.ssc,x=Z)

        if(all(abs(e.res$fitted.value-eTest.res$fitted.value)<1e-6)){
            ls.out <- lapply(1:n.strata, function(iStrata){
                iParamVar.type <- paramVar.type[paramVar.strata == iStrata]
                iOut <- sqrt(eTest.res$coef[names(iParamVar.type)])
                if(any("k" %in% iParamVar.type)){
                    iOut[iParamVar.type=="k"] <- iOut[iParamVar.type=="k"]/iOut[iParamVar.type=="sigma"]
                }
                return(iOut)
            })
            out <- unlist(ls.out)[names(paramVar.type)]
        }else{ ## failure of the full interaction model. Use a log transform
            e.res <- stats::lm.fit(y=log(epsilon2.ssc),x=X)
            out <- exp(0.5*e.res$coef)
        }
    }

    ##  standardize residuals
    if(identical(attr(residuals,"studentized"),TRUE)){
        attr(residuals,"studentized") <- NULL
        attr(out,"studentized") <- rep(NA,n.obs)
        attr(out,"studentized")[M.res[,"index"]] <- M.res[,"residuals"]/exp(X %*% log(out))
    }

    ## export
    return(out)
}

## * initialize.CS
.initialize.CS <- function(object, residuals, Xmean){

    out <- stats::setNames(rep(NA, NROW(object$param)), object$param$name)

    ## extract information
    param.type <- stats::setNames(object$param$type,object$param$name)
    param.strata <- stats::setNames(object$param$strata,object$param$name)
    Upattern.name <- object$X$Upattern$name
    cluster.pattern <- object$X$cluster.pattern

    ## estimate variance and standardize residuals
    attr(residuals,"studentized") <- TRUE ## to return studentized residuals
    sigma <- .initialize.IND(object = object, residuals = residuals, Xmean = Xmean)
    residuals.studentized <- attr(sigma, "studentized")
    attr(sigma, "studentized") <- NULL
    out[names(sigma)] <- sigma

    if(is.null(object$X$cor)){return(out)}
    ## combine all residuals and all design matrices
    M.prodres <- do.call(rbind,lapply(1:length(object$X$cor), function(iPattern){ ## iPattern <- 1
        X.iPattern <- object$X$cor[[iPattern]]
        if(is.null(X.iPattern)){return(NULL)}

        ## number of timepoints
        iN.time <- length(attr(object$X$cor[[iPattern]],"index.Utime"))
        ## position in the matrix associated to each correlation parameter
        indicator.param <- attr(X.iPattern,"indicator.param")
        ## index of the residuals belonging to each individual
        obs.iPattern <- do.call(rbind,attr(X.iPattern,"index.obs"))
        ## identify non-duplicated pairs of observation (here restrict matrix to its  upper part)
        iIndexUpper <- which(attr(X.iPattern,"index.pairtime")[1,]<attr(X.iPattern,"index.pairtime")[2,])
        iPair <- attr(X.iPattern,"index.pairtime")[,iIndexUpper,drop=FALSE]
        ## convert index among all time to index among all time available for this pair
        iIndex.Utime <- attr(X.iPattern, "index.Utime")
        iPair2 <- matrix(NA, nrow = 2, ncol = NCOL(iPair))
        iPair2[] <- as.numeric(factor(iPair, levels = iIndex.Utime, labels = 1:length(iIndex.Utime)))
        iIndex.pair <- iPair2[1,] + (iPair2[2,]-1)*iN.time
                
        iLs.out <- lapply(names(indicator.param), function(iParam){
            if(length(indicator.param[[iParam]])>0){
                iIndexParam.pair <- iPair2[,iIndex.pair %in% indicator.param[[iParam]],drop=FALSE]
                
                iOut <- data.frame(prod = sum(residuals.studentized[obs.iPattern[,iIndexParam.pair[1,]]]*residuals.studentized[obs.iPattern[,iIndexParam.pair[2,]]]),
                                   sum1 = sum(residuals.studentized[obs.iPattern[,iIndexParam.pair[1,]]]),
                                   sum2 = sum(residuals.studentized[obs.iPattern[,iIndexParam.pair[2,]]]),
                                   sums1 = sum(residuals.studentized[obs.iPattern[,iIndexParam.pair[1,]]]^2),
                                   sums2 = sum(residuals.studentized[obs.iPattern[,iIndexParam.pair[2,]]]^2),
                                   n = NROW(obs.iPattern)*NCOL(iIndexParam.pair),
                                   param = iParam)
            }else{
                iOut <- data.frame(prod = 0, sum1 = 0, sum2 = 0, sums1 = 0, sums2 = 0,
                                   n = 0,
                                   param = iParam)
            }
        })


        return(do.call(rbind,iLs.out))
    }))
    ## estimate correlation
    param.rho <- names(param.type)[param.type=="rho"]

    for(iRho in param.rho){ ## iRho <- param.rho[3]
        iIndex <- which(M.prodres$param==iRho)
        iN <- sum(M.prodres[iIndex,"n"])
        iNum <- sum(M.prodres[iIndex,"prod"])/iN-(sum(M.prodres[iIndex,"sum1"])/iN)*(sum(M.prodres[iIndex,"sum2"])/iN)
        iDenom1 <- sum(M.prodres[iIndex,"sums1"])/iN-(sum(M.prodres[iIndex,"sum1"])/iN)^2
        iDenom2 <- sum(M.prodres[iIndex,"sums2"])/iN-(sum(M.prodres[iIndex,"sum2"])/iN)^2
        out[iRho] <- iNum/sqrt(iDenom1*iDenom2)
    }

    ## export
    return(out)
}

## * initialize.UN
.initialize.UN <- .initialize.CS



##----------------------------------------------------------------------
### structure-initialization.R ends here
