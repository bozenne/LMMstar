### structure-initialize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:20) 
## Version: 
## Last-Updated: apr 15 2026 (16:18) 
##           By: Brice Ozenne
##     Update #: 738
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * initialize
##' @title Initialize Variance-Covariance Structure
##' @description Initialize the parameters of the variance-covariance structure using residual variance and correlations.
##' @noRd
##'
##' @param structure [structure]
##' @param residuals [vector] vector of residuals.
##' @param Xmean [matrix] design matrix for the mean effects used to estimate the residual degrees-of-freedom for variance calculation.
##' @param index.cluster [list of numeric vectors] position of the observations of each cluster.
##' @param init.cor [1,2] method to initialize the correlation parameters.
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
##' Sid1 <- .skeleton(IND(~1, var.ordering = "time"), data = dd)
##' Sid4 <- .skeleton(IND(~1|id, var.ordering = "time"), data = dd)
##' Sdiag1 <- .skeleton(IND(~visit), data = dd)
##' Sdiag4 <- .skeleton(IND(~visit|id), data = dd)
##' Sdiag24 <- .skeleton(IND(~visit+gender|id, var.ordering = "time"), data = gastricbypassL)
##'
##' .initialize(Sid1, residuals = residuals(eDD.lm))
##' ## sd(residuals(eDD.lm))
##' .initialize(Sdiag4, residuals = residuals(eDD.lm))
##' ## tapply(residuals(eDD.lm),dd$visit,sd)
##' .initialize(Sdiag24, residuals = residuals(eGas.lm))
##' ## tapply(residuals(eGas.lm),interaction(gastricbypassL[,c("visit","gender")]),sd)
##' 
##' ## compound symmetry
##' Scs4 <- .skeleton(CS(~1|id, var.ordering = "time"), data = gastricbypassL)
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
    function(object, init.cor, method.fit, residuals, Xmean, index.cluster) UseMethod(".initialize")
`.initialize2` <-
    function(object, index.clusterTime, Omega) UseMethod(".initialize2")

## * initialize.ID
.initialize.ID <- function(object, init.cor, method.fit, residuals, Xmean, index.cluster){
    
    ## ** extract information
    structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    param.type <- stats::setNames(structure.param$type,structure.param$name)
    param.strata <- stats::setNames(structure.param$index.strata,structure.param$name)
    Upattern.var <- getGroups(object, form = "variance", data = "pattern")
    
    ## ** combine all residuals and all design matrices
    M.res <- do.call(rbind,lapply(Upattern.var, function(iPattern){ ## iPattern <- 1
        cluster.iPattern <- getGroups(object, form = "variance", level = iPattern, data = "cluster")
        strata.iPattern <- getGroups(object, form = "variance", level = iPattern, data = "strata")
        lp.iPattern <- getGroups(object, form = "variance", level = iPattern, data = "lp")
        obs.iPattern <- unlist(index.cluster[cluster.iPattern])
        
        iOut <- data.frame(index.lp = lp.iPattern, ## recycled to match cluster length
                           index.obs = obs.iPattern,
                           index.strata = strata.iPattern,  ## recycled to match cluster length
                           residuals = residuals[obs.iPattern])
        return(iOut)
    }))

    ## ** extract information
    epsilon2 <- M.res$residuals^2
    X <- object$var$lp2X[M.res$index.lp,,drop=FALSE]
    paramVar.type <- param.type[colnames(X)]
    paramVar.strata <- param.strata[colnames(X)]
    n.strata <- length(unique(paramVar.strata))
    n.obs <- NROW(X)

    ## ** small sample correction (inflate residuals, n-df)
    vec.hat <- rowSums(Xmean %*% solve(t(Xmean) %*% Xmean) * Xmean)[M.res$index.obs]
    if(method.fit == "REML" && !is.null(Xmean) && NCOL(Xmean)>0){
        M.res$index.lpstrata <- paste(M.res$index.lp,M.res$index.strata,sep=".")
        p <- tapply(vec.hat, M.res$index.lpstrata,sum)
        n.UX <- table(M.res$index.lpstrata)
        epsilon2.ssc <- epsilon2 * (n.UX/(n.UX-p))[M.res$index.lpstrata]        
    }else{
        epsilon2.ssc <- epsilon2
    }

    ## ** fit
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

    ## ** standardize residuals
    if(identical(attr(residuals,"studentized"),TRUE)){
        attr(residuals,"studentized") <- NULL
        attr(out,"studentized") <- rep(NA,n.obs)
        attr(out,"studentized")[M.res[,"index.obs"]] <- M.res[,"residuals"]/exp(X %*% log(out))
    }

    ## ** check values
    if(any(abs(out)<1e-10)){
        warning("Some of the variance parameter are initialized to a nearly null value. \n",
                "Parameters: \"",paste(names(out[which(abs(out)<1e-10)]), collapse = "\", \""),"\". \n")
    }
    if(any(out< -1e-10)){
        warning("Some of the variance parameter are initialized to a negative value. \n",
                "Parameters: \"",paste(names(out[which(out < -1e-10)]), collapse = "\", \""),"\". \n")
    }

    ## ** export
    attr(out,"df") <- vec.hat
    return(out)
}

## * initialize2.ID
.initialize2.ID <- function(object, index.clusterTime, Omega){

    ## ** extract information
    structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    param.type <- stats::setNames(structure.param$type,structure.param$name)
    Upattern.var <- getGroups(object, form = "variance", data = "pattern")
        
    ## ** combine all design matrices
    ls.XY <- stats::setNames(lapply(Upattern.var, function(iPattern){ ## iPattern <- Upattern.var[2]
        ## index of the observations belonging to each cluster
        cluster.iPattern <- getGroups(object, form = "variance", level = iPattern, data = "cluster")
        ncluster.iPattern <- length(cluster.iPattern)
        ## repetitions corresponding to each cluster
        ## may not be identical despite same Omega (e.g. CS structure) as the first cluster maybe be 1,3,4 while the second is 1,2,3
        time.iPattern <- sapply(index.clusterTime[cluster.iPattern],paste,collapse="")
        tableTime.iPattern <- table(time.iPattern)
        Y.iPattern <- unlist(lapply(names(tableTime.iPattern), function(iTime){ ## iTime <- names(tableTime.iPattern)[1]
            ## use the first cluster of the pattern with a given vector of times
            iIndex <- index.clusterTime[[cluster.iPattern[which(time.iPattern==iTime)][1]]]
            diag(Omega[iIndex,iIndex,drop=FALSE])
        }))
        ## design matrix with parameters for the pattern
        X.iPattern <- getGroups(object, form = "variance", level = iPattern, data = "X")
        iOut <- list(X = do.call(rbind,replicate(X.iPattern, n = length(tableTime.iPattern), simplify = FALSE)),
                     Y = Y.iPattern,
                     n = do.call(c,lapply(tableTime.iPattern, rep, times = NROW(X.iPattern))))
        
        return(iOut)
    }), Upattern.var)

    X.Omega <- do.call(rbind,lapply(ls.XY,"[[","X"))
    logY.Omega <- log(do.call(c,lapply(ls.XY,"[[","Y")))
    n.Omega <- do.call(c,lapply(ls.XY,"[[","n"))

    ## ** log linear regression
    df.data <- data.frame(Y = logY.Omega, X.Omega)
    form.txt <- paste0("Y~0+",paste(names(df.data)[-1], collapse = "+")) ## NOTE: normalize name in presence of interactions (sigma:1 -> sigma.1)
    e.lm <- stats::lm(stats::as.formula(form.txt),
                      data = data.frame(Y = logY.Omega, X.Omega),
                      weights = n.Omega)
    out <- stats::setNames(sqrt(exp(stats::coef(e.lm))), colnames(X.Omega))
    
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
.initialize.CS <- function(object, init.cor, method.fit, residuals, Xmean, index.cluster){

    structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    out <- stats::setNames(rep(NA, NROW(structure.param)), structure.param$name)

    ## ** extract information
    param.type <- stats::setNames(structure.param$type,structure.param$name)
    param.rho <- names(param.type)[param.type=="rho"]
    Upattern.cor <- getGroups(object, form = "correlation", data = "pattern")

    ## ** estimate variance and standardize residuals
    attr(residuals,"studentized") <- TRUE ## to return studentized residuals
    if("sigma" %in% param.type || "k" %in% param.type){
        sigma <- .initialize.IND(object = object, method.fit = method.fit, residuals = residuals, Xmean = Xmean, index.cluster = index.cluster)
        residuals.studentized <- attr(sigma, "studentized")
        residuals.df <- attr(sigma, "df")
        attr(sigma, "studentized") <- NULL
        attr(sigma, "df") <- NULL
        out[names(sigma)] <- sigma
    }else{
        residuals.studentized <- residuals
    }

    ## ** combine all residuals and all design matrices
    M.prodres <- do.call(rbind,lapply(Upattern.cor, function(iPattern){ ## iPattern <- 1
        ## parametrisation of the correlation structure
        X.iPattern <- getGroups(object, form = "correlation", level = iPattern, data = "Xpattern")[,,"rho"] ## from array to matrix
        if(length(X.iPattern) %in% 0:1){return(NULL)} ## handle pattern with single timepoint
        ## identify non-duplicated pairs of observations (here restrict matrix to its upper part)
        X.iPattern[lower.tri(X.iPattern)] <- "one"
        iPair <- data.frame(which(X.iPattern!="one", arr.ind = TRUE), param = X.iPattern[which(X.iPattern!="one")])
        iPair$param.num <- as.numeric(factor(iPair$param, levels = param.rho))
        ## index of the observations belonging to each cluster
        cluster.iPattern <- getGroups(object, form = "variance", level = iPattern, data = "cluster")
        ncluster.iPattern <- length(cluster.iPattern)            
        
        if(NROW(iPair)<=ncluster.iPattern){ ## more individuals than pairs

            ls.orderingObs <- tapply(unlist(index.cluster[cluster.iPattern]),rep(1:NCOL(X.iPattern), ncluster.iPattern), FUN = identity, simplify = FALSE)
            
            iLs.out <- lapply(1:NROW(iPair), function(iP){ ## iP <- 1
                iRow <- iPair[iP,"row"]
                iCol <- iPair[iP,"col"]
                iOut <- data.frame(param = iPair[iP,"param"],
                                   prod = sum(residuals.studentized[ls.orderingObs[[iRow]]]*residuals.studentized[ls.orderingObs[[iCol]]]),
                                   sum1 = sum(residuals.studentized[ls.orderingObs[[iRow]]]),
                                   sum2 = sum(residuals.studentized[ls.orderingObs[[iCol]]]),
                                   sums1 = sum(residuals.studentized[ls.orderingObs[[iRow]]]^2),
                                   sums2 = sum(residuals.studentized[ls.orderingObs[[iCol]]]^2),
                                   df1 = sum(residuals.df[ls.orderingObs[[iRow]]]),
                                   df2 = sum(residuals.df[ls.orderingObs[[iCol]]])
                                   )
                return(iOut)
            })
            
        }else{ ## more pairs than individuals
            iLs.out <- lapply(cluster.iPattern, function(iC){ ## iC <- cluster.iPattern[1]
                ## residual and df for each element of all possible pairs
                iResRow <- residuals.studentized[index.cluster[[iC]][iPair$row]]
                iResCol <- residuals.studentized[index.cluster[[iC]][iPair$col]]
                iDfRow <- residuals.df[index.cluster[[iC]][iPair$row]]
                iDfCol <- residuals.df[index.cluster[[iC]][iPair$col]]

                iOut <- data.frame(param = tapply(iPair$param,iPair$param,unique),
                                   prod = tapply(iResRow * iResCol,iPair$param,sum),
                                   sum1 = tapply(iResRow,iPair$param,sum),
                                   sum2 = tapply(iResCol,iPair$param,sum),
                                   sums1 = tapply(iResRow^2,iPair$param,sum),
                                   sums2 = tapply(iResCol^2,iPair$param,sum),
                                   df1 = tapply(iDfRow,iPair$param,sum),
                                   df2 = tapply(iDfCol,iPair$param,sum))
                return(iOut)
            })
        }
        iDf.out <- do.call(rbind,iLs.out)
        iLs.out <- by(iDf.out[-1], iDf.out$param, colSums, simplify = FALSE)
        return(data.frame(pattern = iPattern, param = names(iLs.out), n = ncluster.iPattern, do.call(rbind,iLs.out)))
    }))
    
    ## ** estimate correlation
    param.rho <- names(param.type)[param.type=="rho"]
    if(length(param.rho)==0){return(out)}

    e.rho <- unlist(lapply(split(M.prodres, M.prodres$param), function(iDF){

        if(init.cor==1){ 
            ## *** method 1: average ordering-specific correlations (exact formula for ML when no missing values)
            iRho <- iDF$prod/sqrt((iDF$n-iDF$df1)*(iDF$n-iDF$df2))
            
            iMeanRho.pattern <- tapply(iRho, iDF$pattern, mean)
            iNobs.pattern <- tapply((iDF$n-iDF$df1), iDF$pattern, sum)+tapply((iDF$n-iDF$df2), iDF$pattern, sum)
            iHeterochedastic.pattern <- (tapply(iDF$sums1, iDF$pattern, sum)+tapply(iDF$sums2, iDF$pattern, sum))/iNobs.pattern
            iP.pattern <- tapply(iRho, iDF$pattern, length) ## p(p-1)/2
            iN.pattern <- tapply(iDF$n, iDF$pattern, unique)

            iOut <- sum(iN.pattern*iP.pattern*iMeanRho.pattern/sum(iN.pattern*iP.pattern*iHeterochedastic.pattern))
            
        }else if(init.cor==2){
            ## *** method 2: compute overall correlation
            iNobs <- sum(iDF$n)
            iNum <- sum(iDF$prod)/iNobs-(sum(iDF$sum1)/iNobs)*(sum(iDF$sum2)/iNobs)
            iDenom1 <- sum(iDF$sums1)/iNobs-(sum(iDF$sum1)/iNobs)^2
            iDenom2 <- sum(iDF$sums2)/iNobs-(sum(iDF$sum2)/iNobs)^2

            iOut <- iNum/sqrt(iDenom1*iDenom2)
        }
        return(iOut)
    }))

    if(any(is.na(e.rho)) || any(is.infinite(e.rho))){ ## take care of extreme cases, e.g. 0 variability
        e.rho[is.na(e.rho) | is.infinite(e.rho)] <- 0
    }
    out[names(e.rho)] <- e.rho
    ## export
    return(out)
}

## * initialize2.CS
.initialize2.CS <- function(object, index.clusterTime, Omega){
    
    ## ** variance
    structure.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    out <- stats::setNames(rep(NA, NROW(structure.param)), structure.param$name)

    sigma <- .initialize2.IND(object = object, index.clusterTime = index.clusterTime, Omega = Omega)
    out[names(sigma)] <- sigma

    ## ** correlation
    Rho <- stats::cov2cor(Omega)

    param.type <- stats::setNames(structure.param$type,structure.param$name)
    param.strata <- stats::setNames(structure.param$index.strata,structure.param$name)
    param.rho <- names(param.type)[param.type=="rho"]
    Upattern.cor <- getGroups(object, form = "correlation", data = "pattern")

    ls.XY <- stats::setNames(lapply(Upattern.cor, function(iPattern){ ## iPattern <- Upattern.cor[1]

        ## design matrix with parameters for the pattern
        X.iPattern <- getGroups(object, form = "correlation", level = iPattern, data = "Xpattern")[,,"rho"]
        if(length(X.iPattern) %in% 0:1){return(NULL)} ## handle pattern with single timepoint
        ## identify non-duplicated pairs of observations (here restrict matrix to its upper part)
        X.iPattern[lower.tri(X.iPattern)] <- "one"
        iPair <- data.frame(which(X.iPattern!="one", arr.ind = TRUE), param = X.iPattern[which(X.iPattern!="one")])
        iPair$param.num <- as.numeric(factor(iPair$param, levels = param.rho))
        ## index of the observations belonging to each cluster
        cluster.iPattern <- getGroups(object, form = "correlation", level = iPattern, data = "cluster")
        ncluster.iPattern <- length(cluster.iPattern)
        ## repetitions corresponding to each cluster
        ## may not be identical despite same Omega (e.g. CS structure) as the first cluster maybe be 1,3,4 while the second is 1,2,3
        time.iPattern <- sapply(index.clusterTime[cluster.iPattern],paste,collapse="")
        tableTime.iPattern <- table(time.iPattern)
        Y.iPattern <- unlist(lapply(names(tableTime.iPattern), function(iTime){ ## iTime <- names(tableTime.iPattern)[1]
            ## use the first cluster of the pattern with a given vector of times
            iIndex <- index.clusterTime[[cluster.iPattern[which(time.iPattern==iTime)][1]]]
            Rho[iIndex,iIndex,drop=FALSE][which(X.iPattern!="one")]            
        }))
        iOut <- list(X = rep(iPair$param, times = length(tableTime.iPattern)),
                     Y = Y.iPattern,
                     n = do.call(c,lapply(tableTime.iPattern, rep, times = NROW(X.iPattern)*(NROW(X.iPattern)-1)/2)))
        return(iOut)
    }), Upattern.cor)

    X.Omega <- do.call(c,lapply(ls.XY,"[[","X"))
    atanhY.Omega <- atanh(do.call(c,lapply(ls.XY,"[[","Y")))
    n.Omega <- do.call(c,lapply(ls.XY,"[[","n"))

    ## ** log linear regression
    df.data <- data.frame(Y = atanhY.Omega, param = X.Omega)
    df.data$param <- factor(df.data$param)
    if(length(levels(df.data$param))==1){
        out[levels(df.data$param)] <- stats::weighted.mean(tanh(df.data$Y), w = n.Omega)
    }else{
        e.lm <- stats::lm(Y~0+param,
                          data = df.data,
                          weights = n.Omega)
        out[levels(df.data$param)] <- as.numeric(tanh(stats::coef(e.lm)))
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
                            param, method.fit, weights){

    ## ** check feasibility
    requireNamespace("lme4")
    if(!inherits(structure,"RE")){
        stop("Initializer \"lmer\" only available for random effect structures.")
    }
    if(!is.na(structure$name$strata)){
        stop("Initializer \"lmer\" cannot handle multiple strata.")
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
### structure-initialize.R ends here
