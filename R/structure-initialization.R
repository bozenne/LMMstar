### structure-initialization.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:20) 
## Version: 
## Last-Updated: maj 10 2024 (11:47) 
##           By: Brice Ozenne
##     Update #: 427
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
    function(object, method.fit, residuals, Xmean, index.cluster) UseMethod(".initialize")
`.initialize2` <-
    function(object, index.clusterTime, Omega) UseMethod(".initialize2")

## * initialize.ID
.initialize.ID <- function(object, method.fit, residuals, Xmean, index.cluster){

    structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    param.type <- stats::setNames(structure.param$type,structure.param$name)
    param.strata <- stats::setNames(structure.param$index.strata,structure.param$name)
    Upattern.name <- object$Upattern$name

    ## combine all residuals and all design matrices
    M.res <- do.call(rbind,lapply(1:length(object$var$Xpattern), function(iPattern){ ## iPattern <- 1
        X.iPattern <- object$var$Xpattern[[iPattern]]
        cluster.iPattern <- attr(X.iPattern,"index.cluster")
        obs.iPattern <- unlist(index.cluster[cluster.iPattern])
        iOut <- cbind(index.lp = object$var$lp[obs.iPattern],
                      index.obs = obs.iPattern,
                      index.strata = attr(X.iPattern,"index.strata"),
                      residuals = residuals[obs.iPattern],
                      do.call(rbind,rep(list(X.iPattern),length(cluster.iPattern))))
        return(iOut)
    }))

    ## extract information
    epsilon2 <- M.res[,"residuals"]^2
    X <- M.res[,-(1:4),drop=FALSE]
    paramVar.type <- param.type[colnames(X)]
    paramVar.strata <- param.strata[colnames(X)]
    n.strata <- length(unique(paramVar.strata))
    n.obs <- NROW(X)

    ## small sample correction (inflate residuals)
    if(method.fit == "REML" && !is.null(Xmean) && NCOL(Xmean)>0){
        ## n - df
        ## vec.hat <- diag(Xmean %*% solve(t(Xmean) %*% Xmean) %*% t(Xmean))
        vec.hat <- rowSums(Xmean %*% solve(t(Xmean) %*% Xmean) * Xmean)
        strata.mu <- attr(Xmean,"strata")
        if(any(is.na(strata.mu))){ ## partially or non-stratified mean structure: considered as non-stratified
            index.clusterStrata <-  sapply(object$var$Xpattern,attr,"index.strata")[object$var$pattern]
            index.strata <- index.clusterStrata[attr(index.cluster,"vectorwise")]
            p <- tapply(vec.hat,index.strata,sum)
            n.UX <- table(index.strata)
            epsilon2.ssc <- epsilon2 * (n.UX/(n.UX-p))[M.res[,"index.strata"]]
        }else{ ## fully stratified mean and variance structure
            p <- tapply(1:n.obs,object$var$lp,function(iIndex){sum(vec.hat[iIndex])})
            n.UX <- table(object$var$lp)
            epsilon2.ssc <- epsilon2 * (n.UX/(n.UX-p))[M.res[,"index.lp"]]
        }
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
.initialize2.ID <- function(object, index.clusterTime, Omega){

    structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    param.type <- stats::setNames(structure.param$type,structure.param$name)
    param.strata <- stats::setNames(structure.param$index.strata,structure.param$name)
    Upattern <- object$Upattern
    Upattern.name <- Upattern$name
    Omega.diag <- diag(Omega)

    ## ** combine all design matrices
    ls.XY <- stats::setNames(lapply(Upattern.name, function(iPattern){ ## iPattern <- Upattern.name[1]
        iX <- object$var$Xpattern[[Upattern[Upattern$name==iPattern,"var"]]]

        ## NOTE: this handles the case where the pattern is the same for time 1,2,4 and 1,2,3 but maybe not the initialization matrix
        iIndex.cluster <- attr(iX,"index.cluster")
        iIndex.clusterTime <- index.clusterTime[iIndex.cluster]
        iPattern.clusterTime <- nlme::collapse(do.call(rbind,iIndex.clusterTime), as.factor = TRUE)
        iPatternTime.n <- table(iPattern.clusterTime)
        
        attr(iX,c("index.cluster")) <- NULL
        attr(iX,c("index.strata")) <- NULL
        attr(iX,c("param")) <- NULL
        attr(iX,c("indicator.param")) <- NULL
        attr(iX,c("Mindicator.param")) <- NULL

        iOut <- list(X = NULL, Y = NULL, n = NULL)
        for(iPattern.time in levels(iPattern.clusterTime)){ ## iPattern.time <- levels(iPattern.clusterTime)[1]
            iOut$X <- rbind(iOut$X,iX)
            iOut$Y <- c(iOut$Y,Omega.diag[iIndex.clusterTime[[which(iPattern.clusterTime==iPattern.time)[1]]]])
            iOut$n <- c(iOut$n,rep(iPatternTime.n[[iPattern.time]], Upattern[Upattern$name==iPattern,"n.time"]))
        }
        return(iOut)
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
.initialize.CS <- function(object, method.fit, residuals, Xmean, index.cluster){

    structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    out <- stats::setNames(rep(NA, NROW(structure.param)), structure.param$name)

    ## ** extract information
    param.type <- stats::setNames(structure.param$type,structure.param$name)
    param.strata <- stats::setNames(structure.param$index.strata,structure.param$name)
    Upattern.name <- object$Upattern$name

    ## ** estimate variance and standardize residuals
    attr(residuals,"studentized") <- TRUE ## to return studentized residuals
    if("sigma" %in% param.type || "k" %in% param.type){
        sigma <- .initialize.IND(object = object, method.fit = method.fit, residuals = residuals, Xmean = Xmean, index.cluster = index.cluster)
        residuals.studentized <- attr(sigma, "studentized")
        attr(sigma, "studentized") <- NULL
        out[names(sigma)] <- sigma
    }else{
        residuals.studentized <- residuals
    }

    ## ** combine all residuals and all design matrices
    M.prodres <- do.call(rbind,lapply(1:length(object$cor$Xpattern), function(iPattern){ ## iPattern <- 1
        X.iPattern <- object$cor$Xpattern[[iPattern]]
        if(is.null(X.iPattern)){return(NULL)}
        ## index of the residuals belonging to each individual
        obs.iPattern <- do.call(rbind,index.cluster[attr(X.iPattern,"index.cluster")])
        ## identify non-duplicated pairs of observation (here restrict matrix to its  upper part)
        iAllPair <- attr(X.iPattern,"index.pair")
        iPair <- iAllPair[iAllPair[,"col"]<iAllPair[,"row"] & iAllPair$param %in% structure.param$name,,drop=FALSE]
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

    ## ** estimate correlation
    param.rho <- names(param.type)[param.type=="rho"]
    if(length(param.rho)==0){return(out)}

    e.rho <- unlist(lapply(split(M.prodres, M.prodres$param), function(iDF){ ## iDF <- split(M.prodres, M.prodres$param)[[3]]
        iN <- sum(iDF$n)
        iNum <- sum(iDF$prod)/iN-(sum(iDF$sum1)/iN)*(sum(iDF$sum2)/iN)
        iDenom1 <- sum(iDF$sums1)/iN-(sum(iDF$sum1)/iN)^2
        iDenom2 <- sum(iDF$sums2)/iN-(sum(iDF$sum2)/iN)^2
        return(iNum/sqrt(iDenom1*iDenom2))
    }))

    ## ** take care of extreme cases, e.g. 0 variability
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
.initialize2.CS <- function(object, index.clusterTime, Omega){

    structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    out <- stats::setNames(rep(NA, NROW(structure.param)), structure.param$name)

    ## ** extract information
    param.type <- stats::setNames(structure.param$type,structure.param$name)
    param.strata <- stats::setNames(structure.param$index.strata,structure.param$name)
    Upattern <- object$Upattern
    Upattern.name <- Upattern$name

    ## ** variance
    sigma <- .initialize2.IND(object = object, index.clusterTime = index.clusterTime, Omega = Omega)
    out[names(sigma)] <- sigma

    ## ** correlation
    Rho <- stats::cov2cor(Omega)

    ls.XY <- stats::setNames(lapply(Upattern.name, function(iPattern){ ## iPattern <- Upattern.name[2]
        iX <- object$cor$Xpattern[[Upattern[Upattern$name==iPattern,"cor"]]]
        if(NROW(iX)==0){return(NULL)}
        index.vec2matrix <- attr(iX, "index.vec2matrix")
        index.pair <- attr(iX, "index.pair")

        ## NOTE: this handles the case where the pattern is the same for time 1,2,4 and 1,2,3 but maybe not the initialization matrix
        iIndex.cluster <- attr(iX,"index.cluster")
        iIndex.clusterTime <- index.clusterTime[iIndex.cluster]
        iPattern.clusterTime <- nlme::collapse(do.call(rbind,iIndex.clusterTime), as.factor = TRUE)
        iPatternTime.n <- table(iPattern.clusterTime)

        iOut <- list(X = NULL, Y = NULL, n = NULL)
        for(iPattern.time in levels(iPattern.clusterTime)){ ## iPattern.time <- levels(iPattern.clusterTime)[1]
            iIndex.time <- iIndex.clusterTime[[which(iPattern.clusterTime==iPattern.time)[1]]]

            iOut$X <- rbind(iOut$X,cbind(param = index.pair$param))
            iOut$Y <- c(iOut$Y, Rho[iIndex.time,iIndex.time,drop=FALSE][index.vec2matrix])
            iOut$n <- c(iOut$n,rep(iPatternTime.n[[iPattern.time]], NROW(index.pair)))
        }
        return(iOut)
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
## .initialize.EXP <- function(object, residuals, Xmean, index.cluster){
##     structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
##     out <- stats::setNames(rep(NA, NROW(structure.param)), structure.param$name)

##     ## ** extract information
##     param.type <- stats::setNames(structure.param$type,structure.param$name)
##     param.strata <- stats::setNames(structure.param$index.strata,structure.param$name)
##     Upattern.name <- object$X$Upattern$name
##     regressor <- stats::setNames(object$param[object$param$type=="rho","code"],object$param[object$param$type=="rho","name"])
    
##     ## estimate variance and standardize residuals
##     attr(residuals,"studentized") <- TRUE ## to return studentized residuals
##     if("sigma" %in% param.type){
##         sigma <- .initialize.IND(object = object, residuals = residuals, Xmean = Xmean, index.cluster = index.cluster)
##         residuals.studentized <- attr(sigma, "studentized")
##         attr(sigma, "studentized") <- NULL
##         out[names(sigma)] <- sigma
##     }else{
##         residuals.studentized <- residuals
##     }

##     if(is.null(object$X$Xpattern.cor)){return(out)}
##     ## combine all residuals and all design matrices
##     M.prodres <- do.call(rbind,lapply(1:length(object$X$Xpattern.cor), function(iPattern){ ## iPattern <- 1
##         X.iPattern <- object$X$Xpattern.cor[[iPattern]]
##         if(is.null(X.iPattern)){return(NULL)}
##         ## index of the residuals belonging to each individual
##         obs.iPattern <- do.call(rbind,index.cluster[attr(X.iPattern,"index.cluster")])
##         ## identify non-duplicated pairs of observation (here restrict matrix to its  upper part)
##         iAllPair <- attr(X.iPattern,"index.pair")
##         iPair <- iAllPair[iAllPair[,"col"]<iAllPair[,"row"],,drop=FALSE]
##         iParam <- unique(iPair$param)
##         iPair$param <- as.numeric(factor(iPair$param, levels = iParam))
##         iPair$time <- X.iPattern[,regressor[iParam]]

##         ## if(NROW(iPair)<=NROW(obs.iPattern)){ ## more individuals than pairs
##             iLs.out <- apply(iPair, 1, function(iRow){
##                 iOut <- data.frame(prod = sum(residuals.studentized[obs.iPattern[,iRow[1]]]*residuals.studentized[obs.iPattern[,iRow[2]]]),
##                                    sum1 = sum(residuals.studentized[obs.iPattern[,iRow[1]]]),
##                                    sum2 = sum(residuals.studentized[obs.iPattern[,iRow[2]]]),
##                                    sums1 = sum(residuals.studentized[obs.iPattern[,iRow[1]]]^2),
##                                    sums2 = sum(residuals.studentized[obs.iPattern[,iRow[2]]]^2),
##                                    n = NROW(obs.iPattern),
##                                    param = iRow[3],
##                                    time = iRow[4])
##                 return(iOut)
##             }, simplify = FALSE)
##         ## }else{ ## more pairs than individuals
##         ##     iLs.out <- apply(obs.iPattern, 1, function(iRow){ ## iRow <- obs.iPattern[1,]
##         ##         iLSDF <- split(data.frame(row = residuals.studentized[iRow[iPair[,"row"]]],
##         ##                                   col = residuals.studentized[iRow[iPair[,"col"]]],
##         ##                                   param = iPair[,"param"],
##         ##                                   time = iPair[,"time"]),
##         ##                        iPair[,"time"])
##         ##         iOut <- lapply(iLSDF, function(iiDF){
##         ##             data.frame(prod = sum(iiDF[,1]*iiDF[,2]),
##         ##                        sum1 = sum(iiDF[,1]),
##         ##                        sum2 = sum(iiDF[,2]),
##         ##                        sums1 = sum(iiDF[,1]^2),
##         ##                        sums2 = sum(iiDF[,2]^2),
##         ##                        n=NROW(iiDF),
##         ##                        param = iiDF[1,3],
##         ##                        time = iiDF[1,4])})
##         ##         return(do.call(rbind,iOut))
##         ##     }, simplify = FALSE)
##         ## }
##         iDf.out <- do.call(rbind,iLs.out)
##         iDf.out$param <- iParam[iDf.out$param]
##         return(iDf.out)
##     }))
    
##     ## estimate correlation
##     param.rho <- names(param.type)[param.type=="rho"]

##     e.rho <- unlist(lapply(split(M.prodres, M.prodres$param), function(iDF){ ## iDF <- split(M.prodres, M.prodres$param)[[3]]

##         iNum <- iDF$prod/iDF$n-(iDF$sum1/iDF$n)*(iDF$sum2/iDF$n)
##         iDenom1 <- iDF$sums1/iDF$n-(iDF$sum1/iDF$n)^2
##         iDenom2 <- iDF$sums2/iDF$n-(iDF$sum2/iDF$n)^2
##         iRho <- iNum/sqrt(iDenom1*iDenom2)
        
##         ## rougth approximation
##         iRho.initMin <- -log(max(iRho))/mean(iDF$time) 
##         iRho.initMean <- -log(mean(iRho))/mean(iDF$time) 
##         iRho.initMax <- -log(min(iRho))/mean(iDF$time) 
##         if(iRho.initMax<=0){return(0)}

##         errorFun <- function(x){sum(iRho - exp(-x*iDF$time))}
##         error.initMin <- errorFun(iRho.initMin)
##         error.initMean <- errorFun(iRho.initMean)
##         error.initMax <- errorFun(iRho.initMax)
##         if(error.initMean<0){
##             lower <- iRho.initMean
##             if(error.initMax>0){
##                 upper <- iRho.initMax
##             }else{
##                 return(0)
##             }
##         }else if(error.initMin<0){
##             lower <- iRho.initMin
##             if(error.initMean>0){
##                 upper <- iRho.initMean
##             }else if(error.initMax>0){
##                 upper <- iRho.initMax
##             }else{
##                 return(0)
##             }
##         }else{
##             return(0)
##         }

##         return(stats::uniroot(f = errorFun, lower = lower, upper = upper)$root)
##     }))

##     ## take care of extreme cases, e.g. 0 variability
##     if(any(is.na(e.rho))){
##         e.rho[is.na(e.rho)] <- 0
##     }
##     if(any(is.infinite(e.rho))){
##         e.rho[is.infinite(e.rho)] <- 0
##     }
##     out[names(e.rho)] <- e.rho

##     ## export
##     return(out)
## }


## * initialize2.CUSTOM
.initialize2.CUSTOM <- function(object, index.clusterTime, Omega){

    out <- stats::setNames(rep(NA, NROW(object$param)), object$param$name)
    if(!is.null(object$init.sigma) && any(!is.na(object$init.sigma))){
        out[names(object$init.sigma[!is.na(object$init.sigma)])] <- object$init.sigma[!is.na(object$init.sigma)]
    }
    if(!is.null(object$init.rho) && any(!is.na(object$init.rho))){
        out[names(object$init.rho[!is.na(object$init.rho)])] <- object$init.rho[!is.na(object$init.rho)]
    }

    return(out)
}


## * initializeLMER
.initializeLMER <- function(formula, structure, data,
                            param, method.fit, weights, scale.Omega){

    ## ** check feasibility
    requireNamespace("lme4")
    if(!inherits(structure,"RE")){
        stop("Initializer \"lmer\" only available for random effect structures.")
    }
    if(!is.na(structure$name$strata)){
        stop("Initializer \"lmer\" cannot handle multiple strata.")
    }

    if(!is.null(scale.Omega)){
        stop("Initializer \"lmer\" cannot handle weighted residual variance-covariance matrix (argument \'scale.Omega\'). \n")
    }

    ## ** estimation via lmer
    formula.lmer <- updateFormula(formula, add.x = structure$ranef$terms)
    if(is.null(weights)){
        e.lmer <- lme4::lmer(formula.lmer, data = data, REML = method.fit=="REML")
    }else{
        e.lmer <- lme4::lmer(formula.lmer, data = data, REML = method.fit=="REML", weights = weights)
    }

    ## ** extract coefficients
    start <- stats::setNames(rep(as.numeric(NA),length(param$name)), param$name)

    ## *** mean
    mu.lmer <- nlme::fixef(e.lmer)
    if(any(sort(names(mu.lmer)) != sort(names(start)[param$type=="mu"]))){
        stop("Cannot use lmer for initialization: something went wrong when retrieving the mean parameters. \n",
             "Names with lmer: \"",paste(names(mu.lmer),collapse="\", \""),"\". \n",
             "Names with lmm: \"",paste(names(start)[param$type=="mu"],collapse="\", \""),"\". \n", sep = "")
    }
    start[names(mu.lmer)] <- mu.lmer

    ## *** variance
    tau.lmer <- as.data.frame(nlme::VarCorr(e.lmer))
    if(any(sum(param$type=="sigma")!=1)){
        stop("Cannot use lmer for initialization: something went wrong when retrieving the variance parameter. \n",
             "Number of variance parameters with lmer: 1. \n",
             "Number of variance parameters with lmm: ",sum(param$type=="sigma"),". \n", sep = "")
    }
    start[param$type=="sigma"] <-  sqrt(sum(tau.lmer$vcov))

    ## *** correlation
    tau.lmer$name <- sapply(strsplit(tau.lmer$grp, split = ":", fixed = TRUE),"[",1) ## lmer collapse names within the hierarchy session:(day:patient)
    if(any(sum(param$type=="rho")!=sum(tau.lmer$name!="Residual"))){
        stop("Cannot use lmer for initialization: something went wrong when retrieving the correlation parameter. \n",
             "Number of random effect parameters with lmer: ",sum(tau.lmer$name!="Residual"),". \n",
             "Number of correlation parameters with lmm: ",sum(param$type=="rho"),". \n", sep = "")
    }
    if(any(sort(unlist(structure$ranef$hierarchy, use.names = FALSE))!=sort(tau.lmer[tau.lmer$name!="Residual","name"]))){
        stop("Cannot use lmer for initialization: something went wrong when retrieving the correlation parameter. \n",
             "Name of the random effect parameters with lmer: \"",paste(tau.lmer[tau.lmer$name!="Residual","name"], collapse = "\", \""),"\". \n",
             "Name of the correlation parameters with lmm: \"",paste(unlist(structure$ranef$hierarchy, use.names = FALSE), collapse = "\", \""),"\". \n", sep = "")
    }
    rho.lmer <- stats::setNames(tau.lmer[tau.lmer$name!="Residual","vcov"] / sum(tau.lmer$vcov), tau.lmer[tau.lmer$name!="Residual","name"])


    n.hierarchy <- length(structure$ranef$param)
    for(iH in 1:n.hierarchy){ ## iH <- 3
        start[structure$ranef$param[[iH]][,1]] <- cumsum(rho.lmer[rownames(structure$ranef$param[[iH]])])
    }

    ## ** export
    return(start)
    
}
##----------------------------------------------------------------------
### structure-initialization.R ends here
