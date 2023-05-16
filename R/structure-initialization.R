### structure-initialization.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:20) 
## Version: 
## Last-Updated: May 14 2023 (11:19) 
##           By: Brice Ozenne
##     Update #: 321
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
##' @param index.cluster [list of numeric vectors] position of the observations of each cluster.
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
    function(object, residuals, Xmean, index.cluster) UseMethod(".initialize")
`.initialize2` <-
    function(object, Omega) UseMethod(".initialize2")

## * initialize.ID
.initialize.ID <- function(object, residuals, Xmean, index.cluster){

    param.type <- stats::setNames(object$param$type,object$param$name)
    param.strata <- stats::setNames(object$param$index.strata,object$param$name)
    Upattern.name <- object$X$Upattern$name

    ## combine all residuals and all design matrices
    M.res <- do.call(rbind,lapply(1:length(object$X$Xpattern.var), function(iPattern){ ## iPattern <- 1
        X.iPattern <- object$X$Xpattern.var[[iPattern]]
        cluster.iPattern <- attr(X.iPattern,"index.cluster")
        obs.iPattern <- unlist(index.cluster[cluster.iPattern])
        iOut <- cbind(index.pattern = iPattern, index.obs = obs.iPattern,
                      residuals = residuals[obs.iPattern], do.call(rbind,rep(list(X.iPattern),length(cluster.iPattern))))
        return(iOut)
    }))

    ## extract information
    epsilon2 <- M.res[,"residuals"]^2
    X <- M.res[,-(1:3),drop=FALSE]
    paramVar.type <- param.type[colnames(X)]
    paramVar.strata <- param.strata[colnames(X)]
    n.strata <- length(unique(paramVar.strata))
    n.obs <- NROW(X)

    ## small sample correction (inflate residuals)
    if(!is.null(Xmean) && NCOL(Xmean)>0){
        ## n - df
        ## vec.hat <- diag(Xmean %*% solve(t(Xmean) %*% Xmean) %*% t(Xmean))
        vec.hat <- rowSums(Xmean %*% solve(t(Xmean) %*% Xmean) * Xmean)
        
        M.indexsigma <- do.call(rbind,lapply(object$X$Xpattern.var, function(iPattern){ ## iPattern <- object$X$Xpattern.var[[1]]
            
            iUX <- interaction(as.data.frame(iPattern), drop = TRUE, sep = ":")
            ## NEW
            iCluster <- attr(iPattern, "index.cluster")
            iNCluster <- length(iCluster)
            iDL <- data.frame(UX = rep(1:NROW(iUX), iNCluster), X = rep(iUX, iNCluster), index = unlist(index.cluster[attr(iPattern, "index.cluster")]))
            ## OLD 
            ## iDW <- data.frame(UX = 1:NROW(iUX), X = iUX, do.call(cbind,index.cluster[attr(iPattern, "index.cluster")]))
            ## iDL2 <- stats::reshape(iDW, direction = "long", idvar = "UX", varying = setdiff(names(iDW),c("UX","X")), v.names = "index")
            ## all(iDL2$X==iDL$X);all(iDL2$index==iDL$index)
            return(iDL[,c("X","index")])
        }))
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

        if(all(abs(e.res$fitted.value-eTest.res$fitted.value)<1e-6) || any(epsilon2.ssc<=0)){
            ls.out <- lapply(1:n.strata, function(iStrata){
                iParamVar.type <- paramVar.type[paramVar.strata == iStrata]
                iOut <- sqrt(eTest.res$coef[names(iParamVar.type)])                
                if(any("k" %in% iParamVar.type)){
                    if(abs(iOut[iParamVar.type=="sigma"])<1e-6){
                        stop("Cannot initialize covariance structure: no residual variability for the reference level. \n",
                             "Consider using a simplified covariance structure, e.g. homoschedastic. \n")
                    }
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
        attr(out,"studentized")[M.res[,"index.obs"]] <- M.res[,"residuals"]/exp(X %*% log(out))
    }

    ## check values
    if(any(abs(out)<1e-10)){
        warning("Some of the variance parameter are initialized to a nearly null value. \n",
                "Parameters: \"",paste(names(out[which(abs(out)<1e-10)]), collapse = "\", \""),"\". \n")
    }
    if(any(out< -1e-10)){
        warning("Some of the variance parameter are initialized to a negative value. \n",
                "Parameters: \"",paste(names(out[which(out < -1e-10)]), collapse = "\", \""),"\". \n")
    }

    ## export
    return(out)
}

## * initialize2.ID
.initialize2.ID <- function(object, Omega){

    param.type <- stats::setNames(object$param$type,object$param$name)
    param.strata <- stats::setNames(object$param$strata,object$param$name)
    Upattern <- object$X$Upattern
    Upattern.name <- Upattern$name
    Omega.diag <- diag(Omega)

    ## ** combine all design matrices
    ls.XY <- stats::setNames(lapply(Upattern.name, function(iPattern){ ## iPattern <- Upattern.name[1]
        iX <- object$X$Xpattern.var[[Upattern[Upattern$name==iPattern,"var"]]]
        attr(iX,c("index.cluster")) <- NULL
        attr(iX,c("index.strata")) <- NULL
        attr(iX,c("param")) <- NULL
        attr(iX,c("indicator.param")) <- NULL
        attr(iX,c("Mindicator.param")) <- NULL
        iOut <- list(X = iX,
                     Y = Omega.diag[object$X$Upattern$time[[iPattern]]],
                     n = rep(Upattern[Upattern$name==iPattern,"n.cluster"], Upattern[Upattern$name==iPattern,"n.time"])
                     )
    }), Upattern.name)
    
    X.Omega <- do.call(rbind,lapply(ls.XY,"[[","X"))
    logY.Omega <- log(do.call(c,lapply(ls.XY,"[[","Y")))
    n.Omega <- do.call(c,lapply(ls.XY,"[[","n"))

    ## ** log linear regression
    df.data <- data.frame(Y = logY.Omega, X.Omega)
    form.txt <- paste0("Y~0+",paste(names(df.data)[-1], collapse = "+")) ## NOTE: normalize name in presence of interactions (sigma:1 -> sigma.1)
    e.lm <- stats::lm(stats::as.formula(form.txt),
                      data = data.frame(Y = logY.Omega, X.Omega),
                      weights = n.Omega)
    out <- stats::setNames(sqrt(exp(coef(e.lm))), colnames(X.Omega))
    
    ## ** check values
    param.sigma <- names(param.type)[param.type=="sigma"]
    if(any(abs(out[param.sigma])<1e-10)){
        warning("Some of the variance parameter are initialized to a nearly null value. \n",
                "Parameters: \"",paste(param.sigma[which(abs(out[param.sigma])<1e-10)], collapse = "\", \""),"\". \n")
    }

    ## ** export
    return(out)
}


## * initialize.IND, initialize2.IND
.initialize.IND <- .initialize.ID
.initialize2.IND <- .initialize2.ID

## * initialize.CS
.initialize.CS <- function(object, residuals, Xmean, index.cluster){
    out <- stats::setNames(rep(NA, NROW(object$param)), object$param$name)

    ## extract information
    param.type <- stats::setNames(object$param$type,object$param$name)
    param.strata <- stats::setNames(object$param$index.strata,object$param$name)
    Upattern.name <- object$X$Upattern$name

    ## estimate variance and standardize residuals
    attr(residuals,"studentized") <- TRUE ## to return studentized residuals
    if("sigma" %in% param.type){
        sigma <- .initialize.IND(object = object, residuals = residuals, Xmean = Xmean, index.cluster = index.cluster)
        residuals.studentized <- attr(sigma, "studentized")
        attr(sigma, "studentized") <- NULL
        out[names(sigma)] <- sigma
    }else{
        residuals.studentized <- residuals
    }

    if(is.null(object$X$Xpattern.cor)){return(out)}
    ## combine all residuals and all design matrices
    M.prodres <- do.call(rbind,lapply(1:length(object$X$Xpattern.cor), function(iPattern){ ## iPattern <- 1
        X.iPattern <- object$X$Xpattern.cor[[iPattern]]
        if(is.null(X.iPattern)){return(NULL)}
        ## index of the residuals belonging to each individual
        obs.iPattern <- do.call(rbind,index.cluster[attr(X.iPattern,"index.cluster")])
        ## identify non-duplicated pairs of observation (here restrict matrix to its  upper part)
        iAllPair <- attr(X.iPattern,"index.pair")
        iPair <- iAllPair[iAllPair[,"col"]<iAllPair[,"row"],,drop=FALSE]
        iParam <- unique(iPair$param)
        iPair$param <- as.numeric(factor(iPair$param, levels = iParam))

        if(NROW(iPair)<=NROW(obs.iPattern)){ ## more individuals than pairs
            iLs.out <- apply(iPair, 1, function(iRow){
                iOut <- data.frame(prod = sum(residuals.studentized[obs.iPattern[,iRow[1]]]*residuals.studentized[obs.iPattern[,iRow[2]]]),
                                   sum1 = sum(residuals.studentized[obs.iPattern[,iRow[1]]]),
                                   sum2 = sum(residuals.studentized[obs.iPattern[,iRow[2]]]),
                                   sums1 = sum(residuals.studentized[obs.iPattern[,iRow[1]]]^2),
                                   sums2 = sum(residuals.studentized[obs.iPattern[,iRow[2]]]^2),
                                   n = NROW(obs.iPattern),
                                   param = iRow[3])
                return(iOut)
            }, simplify = FALSE)
        }else{ ## more pairs than individuals
            iLs.out <- apply(obs.iPattern, 1, function(iRow){ ## iRow <- obs.iPattern[1,]
                iLSDF <- split(data.frame(row = residuals.studentized[iRow[iPair[,"row"]]],
                                          col = residuals.studentized[iRow[iPair[,"col"]]],
                                          param = iPair[,"param"]),
                               iPair[,"param"])
                iOut <- lapply(iLSDF, function(iiDF){
                    data.frame(prod = sum(iiDF[,1]*iiDF[,2]),
                               sum1 = sum(iiDF[,1]),
                               sum2 = sum(iiDF[,2]),
                               sums1 = sum(iiDF[,1]^2),
                               sums2 = sum(iiDF[,2]^2),
                               n=NROW(iiDF),
                               param = iiDF[1,3])})
                return(do.call(rbind,iOut))
            }, simplify = FALSE)
        }
        iDf.out <- do.call(rbind,iLs.out)
        iDf.out$param <- iParam[iDf.out$param]
        return(iDf.out)
    }))
    ## estimate correlation
    param.rho <- names(param.type)[param.type=="rho"]

    e.rho <- unlist(lapply(split(M.prodres, M.prodres$param), function(iDF){ ## iDF <- split(M.prodres, M.prodres$param)[[3]]
        iN <- sum(iDF$n)
        iNum <- sum(iDF$prod)/iN-(sum(iDF$sum1)/iN)*(sum(iDF$sum2)/iN)
        iDenom1 <- sum(iDF$sums1)/iN-(sum(iDF$sum1)/iN)^2
        iDenom2 <- sum(iDF$sums2)/iN-(sum(iDF$sum2)/iN)^2
        return(iNum/sqrt(iDenom1*iDenom2))
    }))

    ## take care of extreme cases, e.g. 0 variability
    if(any(is.na(e.rho))){
        e.rho[is.na(e.rho)] <- 0
    }
    if(any(is.infinite(e.rho))){
        e.rho[is.infinite(e.rho)] <- 0
    }
    out[names(e.rho)] <- e.rho

    ## export
    return(out)
}

## * initialize2.CS
.initialize2.CS <- function(object, Omega){
    out <- stats::setNames(rep(NA, NROW(object$param)), object$param$name)

    ## ** extract information
    param.type <- stats::setNames(object$param$type,object$param$name)
    param.strata <- stats::setNames(object$param$strata,object$param$name)
    Upattern <- object$X$Upattern
    Upattern.name <- Upattern$name

    ## ** variance
    sigma <- .initialize2.IND(object = object, Omega = Omega)
    out[names(sigma)] <- sigma

    ## ** correlation
    Rho <- stats::cov2cor(Omega)

    ls.XY <- stats::setNames(lapply(Upattern.name, function(iPattern){ ## iPattern <- Upattern.name[1]
        iX <- object$X$Xpattern.cor[[Upattern[Upattern$name==iPattern,"cor"]]]
        if(NROW(iX)==0){return(NULL)}
        index.vec2matrix <- attr(iX, "index.vec2matrix")
        index.time <- Upattern$time[[iPattern]] ## NOTE: use Upattern instead of attr(iX, "index.time") as the later is not define when there is ambiguity
                                                ##       e.g. CS with one missing data where the same pattern holds for time 1,2,4 and 1,3,4
        index.pair <- attr(iX, "index.pair")
        iOut <- list(X = cbind(param=index.pair$param),
                     Y = Rho[index.time,index.time][index.vec2matrix],
                     n = rep(Upattern[Upattern$name==iPattern,"n.cluster"], length(index.vec2matrix))
                     )
    }), Upattern.name)
    X.Omega <- do.call(rbind,lapply(ls.XY,"[[","X"))
    atanhY.Omega <- atanh(do.call(c,lapply(ls.XY,"[[","Y")))
    n.Omega <- do.call(c,lapply(ls.XY,"[[","n"))

    ## ** log linear regression
    df.data <- data.frame(Y = atanhY.Omega, X.Omega)
    df.data$param <- factor(df.data$param)
    if(length(levels(df.data$param))==1){
        out[levels(df.data$param)] <- stats::weighted.mean(tanh(df.data$Y), w = n.Omega)
    }else{
        e.lm <- stats::lm(Y~0+param,
                          data = df.data,
                          weights = n.Omega)
        out[levels(df.data$param)] <- as.numeric(tanh(coef(e.lm)))
    }
    
    ## ** export    
    return(out)
}

## * initialize.RE, initialize2.RE
.initialize.RE <- .initialize.CS
.initialize2.RE <- .initialize2.CS

## * initialize.TOEPLITZ, initialize2.TOEPLITZ
.initialize.TOEPLITZ <- .initialize.CS
.initialize2.TOEPLITZ <- .initialize2.CS

## * initialize.UN, initialize2.UN
.initialize.UN <- .initialize.CS
.initialize2.UN <- .initialize2.CS

## * initialize.EXP
.initialize.EXP <- function(object, residuals, Xmean, index.cluster){
    out <- stats::setNames(rep(NA, NROW(object$param)), object$param$name)

    ## extract information
    param.type <- stats::setNames(object$param$type,object$param$name)
    param.strata <- stats::setNames(object$param$strata,object$param$name)
    Upattern.name <- object$X$Upattern$name
    regressor <- stats::setNames(object$param[object$param$type=="rho","code"],object$param[object$param$type=="rho","name"])
    
    ## estimate variance and standardize residuals
    attr(residuals,"studentized") <- TRUE ## to return studentized residuals
    if("sigma" %in% param.type){
        sigma <- .initialize.IND(object = object, residuals = residuals, Xmean = Xmean, index.cluster = index.cluster)
        residuals.studentized <- attr(sigma, "studentized")
        attr(sigma, "studentized") <- NULL
        out[names(sigma)] <- sigma
    }else{
        residuals.studentized <- residuals
    }

    if(is.null(object$X$Xpattern.cor)){return(out)}
    ## combine all residuals and all design matrices
    M.prodres <- do.call(rbind,lapply(1:length(object$X$Xpattern.cor), function(iPattern){ ## iPattern <- 1
        X.iPattern <- object$X$Xpattern.cor[[iPattern]]
        if(is.null(X.iPattern)){return(NULL)}
        ## index of the residuals belonging to each individual
        obs.iPattern <- do.call(rbind,index.cluster[attr(X.iPattern,"index.cluster")])
        ## identify non-duplicated pairs of observation (here restrict matrix to its  upper part)
        iAllPair <- attr(X.iPattern,"index.pair")
        iPair <- iAllPair[iAllPair[,"col"]<iAllPair[,"row"],,drop=FALSE]
        iParam <- unique(iPair$param)
        iPair$param <- as.numeric(factor(iPair$param, levels = iParam))
        iPair$time <- X.iPattern[,regressor[iParam]]

        ## if(NROW(iPair)<=NROW(obs.iPattern)){ ## more individuals than pairs
            iLs.out <- apply(iPair, 1, function(iRow){
                iOut <- data.frame(prod = sum(residuals.studentized[obs.iPattern[,iRow[1]]]*residuals.studentized[obs.iPattern[,iRow[2]]]),
                                   sum1 = sum(residuals.studentized[obs.iPattern[,iRow[1]]]),
                                   sum2 = sum(residuals.studentized[obs.iPattern[,iRow[2]]]),
                                   sums1 = sum(residuals.studentized[obs.iPattern[,iRow[1]]]^2),
                                   sums2 = sum(residuals.studentized[obs.iPattern[,iRow[2]]]^2),
                                   n = NROW(obs.iPattern),
                                   param = iRow[3],
                                   time = iRow[4])
                return(iOut)
            }, simplify = FALSE)
        ## }else{ ## more pairs than individuals
        ##     iLs.out <- apply(obs.iPattern, 1, function(iRow){ ## iRow <- obs.iPattern[1,]
        ##         iLSDF <- split(data.frame(row = residuals.studentized[iRow[iPair[,"row"]]],
        ##                                   col = residuals.studentized[iRow[iPair[,"col"]]],
        ##                                   param = iPair[,"param"],
        ##                                   time = iPair[,"time"]),
        ##                        iPair[,"time"])
        ##         iOut <- lapply(iLSDF, function(iiDF){
        ##             data.frame(prod = sum(iiDF[,1]*iiDF[,2]),
        ##                        sum1 = sum(iiDF[,1]),
        ##                        sum2 = sum(iiDF[,2]),
        ##                        sums1 = sum(iiDF[,1]^2),
        ##                        sums2 = sum(iiDF[,2]^2),
        ##                        n=NROW(iiDF),
        ##                        param = iiDF[1,3],
        ##                        time = iiDF[1,4])})
        ##         return(do.call(rbind,iOut))
        ##     }, simplify = FALSE)
        ## }
        iDf.out <- do.call(rbind,iLs.out)
        iDf.out$param <- iParam[iDf.out$param]
        return(iDf.out)
    }))
    
    ## estimate correlation
    param.rho <- names(param.type)[param.type=="rho"]

    e.rho <- unlist(lapply(split(M.prodres, M.prodres$param), function(iDF){ ## iDF <- split(M.prodres, M.prodres$param)[[3]]

        iNum <- iDF$prod/iDF$n-(iDF$sum1/iDF$n)*(iDF$sum2/iDF$n)
        iDenom1 <- iDF$sums1/iDF$n-(iDF$sum1/iDF$n)^2
        iDenom2 <- iDF$sums2/iDF$n-(iDF$sum2/iDF$n)^2
        iRho <- iNum/sqrt(iDenom1*iDenom2)
        
        ## rougth approximation
        iRho.initMin <- -log(max(iRho))/mean(iDF$time) 
        iRho.initMean <- -log(mean(iRho))/mean(iDF$time) 
        iRho.initMax <- -log(min(iRho))/mean(iDF$time) 
        if(iRho.initMax<=0){return(0)}

        errorFun <- function(x){sum(iRho - exp(-x*iDF$time))}
        error.initMin <- errorFun(iRho.initMin)
        error.initMean <- errorFun(iRho.initMean)
        error.initMax <- errorFun(iRho.initMax)
        if(error.initMean<0){
            lower <- iRho.initMean
            if(error.initMax>0){
                upper <- iRho.initMax
            }else{
                return(0)
            }
        }else if(error.initMin<0){
            lower <- iRho.initMin
            if(error.initMean>0){
                upper <- iRho.initMean
            }else if(error.initMax>0){
                upper <- iRho.initMax
            }else{
                return(0)
            }
        }else{
            return(0)
        }

        return(stats::uniroot(f = errorFun, lower = lower, upper = upper)$root)
    }))

    ## take care of extreme cases, e.g. 0 variability
    if(any(is.na(e.rho))){
        e.rho[is.na(e.rho)] <- 0
    }
    if(any(is.infinite(e.rho))){
        e.rho[is.infinite(e.rho)] <- 0
    }
    out[names(e.rho)] <- e.rho

    ## export
    return(out)
}


## * initialize.CUSTOM
.initialize.CUSTOM <- function(object, residuals, Xmeans, index.cluster){

    out <- stats::setNames(rep(NA, NROW(object$param)), object$param$name)
    if(!is.null(object$init.sigma) && any(!is.na(object$init.sigma))){
        out[names(object$init.sigma[!is.na(object$init.sigma)])] <- object$init.sigma[!is.na(object$init.sigma)]
    }
    if(!is.null(object$init.rho) && any(!is.na(object$init.rho))){
        out[names(object$init.rho[!is.na(object$init.rho)])] <- object$init.rho[!is.na(object$init.rho)]
    }

    return(out)
}


##----------------------------------------------------------------------
### structure-initialization.R ends here
