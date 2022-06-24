### structure-skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: jun 24 2022 (15:08) 
##           By: Brice Ozenne
##     Update #: 2191
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * skeleton
##' @title Parametrization of Covariance Structure
##' @description Parametrization of Covariance Structure
##' @noRd
##'
##' @param structure [structure]
##' @param data [data.frame] dataset
##'
##' @keywords internal
##' 
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$gender <- c("M","F")[as.numeric(gastricbypassL$id) %% 2+1]
##' dd <- gastricbypassL[!duplicated(gastricbypassL[,c("time","gender")]),]
##' 
##' ## independence
##' .skeleton(IND(~1, var.time = "time"), data = dd)
##' .skeleton(IND(~time), data = dd)
##' .skeleton(IND(~time|id), data = dd)
##' .skeleton(IND(~time+gender|id, var.time = "time"), data = dd)
##' .skeleton(IND(gender~time|id, var.time = "time"), data = dd)
##' 
##' ## compound symmetry
##' .skeleton(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' .skeleton(CS(~time|id), data = gastricbypassL)
##' .skeleton(CS(gender~time|id), data = gastricbypassL)
##' 
##' ## unstructured
##' .skeleton(UN(~visit|id), data = gastricbypassL)
##' .skeleton(UN(~time|id), data = gastricbypassL)
##' .skeleton(UN(~visit+gender|id, var.time = "time"), data = gastricbypassL)
##' .skeleton(UN(gender~visit|id), data = gastricbypassL)
`.skeleton` <-
    function(structure, data, indexData) UseMethod(".skeleton")

## * skeleton.ID
.skeleton.ID <- function(structure, data, indexData = NULL){

    ## ** prepare
    if(is.null(indexData)){
        indexData <- .extractIndexData(data = data, structure = structure)
    }
    strata.var <- structure$name$strata
    U.cluster <- indexData$U.cluster
    U.time <- indexData$U.time
    U.strata <- indexData$U.strata
    n.strata <- length(U.strata)
    
    index.clusterTime <- indexData$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    index.clusterStrata <- indexData$index.clusterStrata ## vector of index relative to which strata each cluster belong to
    index.cluster <- indexData$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    X.var <- structure$X$var
    X.cor <- structure$X$cor

    ## ** param
    ## *** sigma
    outSigma <- .initSigma(X.var = X.var, strata.var = strata.var, U.strata = U.strata, n.strata = n.strata)
    X.var <- outSigma$X    
    param.sigma <- outSigma$param    
    strata.sigma <- outSigma$strata    
    code.sigma <- outSigma$code
    level.sigma <- stats::setNames(as.list(outSigma$level),param.sigma)

    ## *** k
    outK <- .initK(X.var = X.var, strata.var = strata.var, U.strata = U.strata, n.strata = n.strata,
                   param.sigma = param.sigma, strata.sigma = strata.sigma)
    X.var <- outK$X    
    param.k <- outK$param    
    strata.k <- outK$strata    
    code.k <- outK$code
    level.k <- stats::setNames(as.list(outK$level),param.k)

    ## *** param rho
    if(NROW(X.cor)>0){
        outRho <- .initRho(data = data, X.cor = X.cor, X.var = c(outK, list(strata.sigma = strata.sigma)), heterogeneous = structure$heterogeneous, 
                           U.cluster = U.cluster, index.cluster = index.cluster,
                           U.time = U.time, index.clusterTime = index.clusterTime, 
                           strata.var = strata.var, U.strata = U.strata, index.clusterStrata = index.clusterStrata, n.strata = n.strata)
        param.rho <- outRho$param[!duplicated(outRho$code)]
        strata.rho <- outRho$strata[!duplicated(outRho$code)]
        code.rho <- outRho$code[!duplicated(outRho$code)]
        level.rho <- tapply(outRho$level,outRho$code,function(iX){iX}, simplify = FALSE)[code.rho]
    }else{
        param.rho <- NULL
        strata.rho <- NULL
        code.rho <- NULL
        level.rho <- NULL
    }

    ## ** gather parameters
    structure$param <- data.frame(name = c(param.sigma,param.k,param.rho),
                                  strata = as.numeric(c(strata.sigma,strata.k,strata.rho)),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("k",length=length(param.k)),rep("rho",length=length(param.rho))),
                                  level = NA,
                                  code = c(code.sigma,code.k,code.rho),
                                  code.x = as.factor(NA),
                                  code.y = as.factor(NA),
                                  sigma = as.character(NA),
                                  k.x = as.character(NA),
                                  k.y = as.character(NA),                                  
                                  stringsAsFactors = FALSE)
    
    structure$param$level <- c(level.sigma,level.k,level.rho)
    attr(structure$param,"level.var") <- outSigma$code
    if(any(structure$param$type=="k")){
        structure$param[structure$param$type=="k","sigma"] <- names(strata.sigma)[structure$param[structure$param$type=="k","strata"]]
        attr(structure$param,"level.var") <- c(attr(structure$param,"level.var"), outK$code)
    }
    if(any(structure$param$type=="rho")){
        structure$param$code.x <- stats::setNames(vector(mode = "list", length = NROW(structure$param)), structure$param$code)
        ls.tempo <- tapply(factor(outRho$lp.x, levels = 1:length(attr(outRho,"levels")), labels = attr(outRho,"levels")), outRho$code, function(iX){iX}, simplify = FALSE)
        structure$param$code.x[names(ls.tempo)] <- ls.tempo
        names(structure$param$code.x) <- structure$param$name

        structure$param$code.y <- stats::setNames(vector(mode = "list", length = NROW(structure$param)), structure$param$code)
        ls.tempo <- tapply(factor(outRho$lp.y, levels = 1:length(attr(outRho,"levels")), labels = attr(outRho,"levels")), outRho$code, function(iX){iX}, simplify = FALSE)
        structure$param$code.y[names(ls.tempo)] <- ls.tempo
        names(structure$param$code.y) <- structure$param$name
        
        structure$param[structure$param$type=="rho","sigma"] <- outRho$sigma[!duplicated(outRho$code)]
        structure$param[structure$param$type=="rho","k.x"] <- outRho$k.x[!duplicated(outRho$code)]
        structure$param[structure$param$type=="rho","k.y"] <- outRho$k.y[!duplicated(outRho$code)]
        attr(structure$param,"level.cor") <- attr(outRho,"levels")
    }
    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL

    ## ** export
    structure$X$var <- X.var
    structure$X$cor <- X.cor
    return(structure)
}


## * skeleton.IND
.skeleton.IND <- .skeleton.ID

## * skeleton.CS
.skeleton.CS <- .skeleton.IND

## * skeleton.UN
.skeleton.UN <- .skeleton.CS

## * skeleton.CUSTOM
.skeleton.CUSTOM <- function(structure, data, indexData = NULL){

    ## ** gather parameters
    structure$param <- NULL
    if(!is.null(structure$FCT.sigma)){
        structure$param <- rbind(structure$param,
                                 data.frame(name = names(structure$init.sigma),
                                            strata = 1,
                                            type = "sigma",
                                            level = as.character(NA),
                                            code = as.character(NA),
                                            code.x = NA,
                                            code.y = NA,
                                            sigma = as.character(NA),
                                            k.x = as.character(NA),
                                            k.y = as.character(NA)))
    }

    if(!is.null(structure$FCT.rho)){
        structure$param <- rbind(structure$param,
                                 data.frame(name = names(structure$init.rho),
                                            strata = 1,
                                            type = "rho",
                                            level = as.character(NA),
                                            code = as.character(NA),
                                            code.x = NA,
                                            code.y = NA,
                                            sigma = as.character(NA),
                                            k.x = as.character(NA),
                                            k.y = as.character(NA)))
    }
    structure$param$code.x <- rep(list(NULL), NROW(structure$param))
    structure$param$code.y <- rep(list(NULL), NROW(structure$param))

    ## ** pattern
    return(structure)
}

## * helpers
## ** .initSigma
.initSigma <- function(X.var, strata.var, U.strata, n.strata, sep = ":"){


    if(n.strata==1){
        param.sigma <- "sigma"
        index.sigma <- which(attr(X.var,"assign")==0)
        if(!identical(colnames(X.var)[index.sigma],"(Intercept)")){
            stop("Could not find the intercept in the design matrix for the variance.\n")
        }
        colnames(X.var)[index.sigma] <- param.sigma
        strata.sigma <- stats::setNames(1,param.sigma)        
    }else{
        param.sigma <- paste0("sigma:",U.strata)
        index.sigma <- which(attr(X.var,"assign")==1)
        if(!identical(colnames(X.var)[index.sigma],paste0(strata.var,U.strata))){
            stop("Could not find the strata-specific intercepts in the design matrix for the variance.\n")
        }
        colnames(X.var)[index.sigma] <- param.sigma
        strata.sigma <- stats::setNames(1:n.strata,param.sigma)
        attr(X.var,"order")[index.sigma] <- 0
        
    }

    ## find code
    X.varsigma <- unique(X.var[rowSums(X.var[,param.sigma,drop=FALSE])==rowSums(X.var),,drop=FALSE])
    code.sigma <- stats::setNames(as.character(interaction(as.data.frame(X.varsigma), drop = TRUE, sep = sep)),
                                  apply(X.varsigma, 1, function(iRow){names(which(iRow==1))}))
    code.sigma <- as.character(code.sigma[param.sigma])

    ## find level
    M.level <- attr(X.var,"M.level")
    if(NCOL(M.level)==0){ ## nothing
        level.sigma <- ""
    }else if(n.strata==1 || NCOL(M.level)>1){ ## no strata or strata with covariate
        level.sigma <- paste0(".",as.character(interaction(M.level[index.sigma,,drop=FALSE],sep=sep, drop = TRUE)))
    }else{ ## only strata
        level.sigma <- paste0(sep,as.character(interaction(M.level[index.sigma,,drop=FALSE],sep=sep, drop = TRUE)))
    }

    return(list(X = X.var,
                param = param.sigma,
                strata = strata.sigma,
                code = code.sigma,
                level = level.sigma))

}

## ** .initK
.initK <- function(X.var, strata.var, U.strata, n.strata,
                   param.sigma, strata.sigma, sep = c(".",":")){

    index.k <- which(attr(X.var,"assign")>(n.strata>1))

    if(length(index.k)>0){
        M.level <- attr(X.var,"M.level")
        oldnames <- as.character(interaction(as.data.frame(lapply(1:NCOL(M.level), function(iVar){paste0(names(M.level)[iVar],M.level[index.k,iVar])})), sep = sep[2], drop = TRUE) )

        if(n.strata==1 || NCOL(M.level)>1){ ## no strata or strata with covariate
            level.k <- paste0(sep[1],as.character(interaction(M.level[index.k,,drop=FALSE], sep = sep[2], drop = TRUE)))
        }else{ ## only strata
            level.k <- as.character(interaction(M.level[index.k,,drop=FALSE], sep = sep[2], drop = TRUE))
        }
        
        if(!identical(colnames(X.var)[index.k],oldnames)){
            stop("Could not find the k parameters in the design matrix for the variance.\n")
        }
        colnames(X.var)[index.k] <- paste0("k",level.k)
        param.k <- colnames(X.var)[index.k]

        ## find code
        X.vark <- unique(X.var[rowSums(X.var[,param.sigma,drop=FALSE])!=rowSums(X.var),,drop=FALSE])
        code.k <- stats::setNames(as.character(interaction(as.data.frame(X.vark), drop = TRUE, sep = sep[2])),
                                  apply(X.vark[,param.k,drop=FALSE], 1, function(iRow){names(which(iRow==1))}))[param.k]
        code.k <- as.character(code.k[param.k])

        if(n.strata == 1){
            strata.k <- stats::setNames(rep(1,length(param.k)),param.k)
        }else{
            strata.k <- stats::setNames(match(M.level[index.k,strata.var],U.strata),param.k)
        }
    }else{
        param.k <- NULL
        strata.k <- NULL
        code.k <- NULL
        level.k <- NULL
    }

    return(list(X = X.var,
                param = param.k,
                strata = strata.k,
                code = code.k,
                level = level.k))
}

## ** .initRho
## for each cluster compute all pairwise difference in covariates to find the parameters
.initRho <- function(data, X.cor, X.var, heterogeneous, 
                     U.cluster, index.cluster,
                     U.time, index.clusterTime, 
                     strata.var, U.strata, index.clusterStrata, n.strata, sep = c(":")){

    fct.envir <- environment()
    n.time <- length(U.time)

    ## *** linear predictor for each observation and then grouped by cluster
    if(length(all.vars(attr(X.cor,"formula")))>0){
        ## reorder by cluster and time and then find unique rows
        nobs.cluster <- table(attr(index.cluster,"vectorwise"))
        UX.cor <- unique(X.cor[unlist(index.cluster[names(sort(nobs.cluster, decreasing = TRUE))]),,drop=TRUE])
        Ulp.cor <- sort(as.character(interaction(as.data.frame(UX.cor), sep = sep, drop=TRUE)), decreasing = TRUE)
    }else{
        Ulp.cor <- "1"
    }
    lpObs.cor <- interaction(as.data.frame(X.cor), sep = sep, drop=TRUE)
    lpnObs.cor <- as.numeric(factor(lpObs.cor, levels = Ulp.cor))
    lpnCluster.cor <- stats::setNames(lapply(U.cluster, function(iC){
        lpnObs.cor[index.cluster[[iC]]]
    }), U.cluster)
    attr(lpnCluster.cor, "levels") <- Ulp.cor

    ## *** identify the correlation parameters
    out <- vector(mode = "list", length = n.strata)

    for(iStrata in 1:n.strata){

        ## **** subset relative to strata
        if(n.strata==1){
            iCluster <- U.cluster
            iIndex.cluster <- index.cluster
            iIndex.clusterTime <- index.clusterTime
            iLpnCluster.cor <- lpnCluster.cor         
        }else{
            iCluster <- names(index.clusterStrata)[which(index.clusterStrata==iStrata)]
            iIndex.cluster <- index.cluster[iCluster]
            iIndex.clusterTime <- index.clusterTime[iCluster]
            iLpnCluster.cor <- lpnCluster.cor[iCluster]            
        }
        
        if(!missing(X.var)){
            iX.sigma <- X.var$X[,names(which(X.var$strata.sigma==iStrata)),drop=FALSE]
            iX.k <- X.var$X[,names(which(X.var$strata==iStrata)),drop=FALSE]
        }
        
        
        ## no replicates
        if(all(sapply(iLpnCluster.cor,length)<=1)){
            next
        }

        ## **** find unique linear predictor values within cluster
        iULpIndex.cor <- stats::setNames(vector(mode = "list", length = length(iCluster)), iCluster)
        iULpCluster.cor <- stats::setNames(lapply(iCluster, function(iId){ ## iId <- iCluster[2]
            ## unique levels
            iTest <- which(!duplicated(iLpnCluster.cor[[iId]]))
            iCode <- iLpnCluster.cor[[iId]][iTest]
            ## add duplicates
            if(length(iTest)>1 && length(iLpnCluster.cor[[iId]])>length(iTest)){
                iTest2 <- (1:length(iLpnCluster.cor[[iId]]))[-iTest][which(!duplicated(iLpnCluster.cor[[iId]][-iTest]))]
                iTest <- sort(c(iTest,iTest2))
                iCode <- iLpnCluster.cor[[iId]][iTest]
            }
            fct.envir$iULpIndex.cor[[iId]] <- sort(iTest)
            return(iCode)
        }), iCluster)

        ## **** remove cluster with single value or identical to other cluster in term of linear predictors
        iIndex.unique <- intersect(which(!duplicated(iULpCluster.cor)),
                                   which(sapply(iIndex.cluster,length)>1))
        iCluster2 <- iCluster[iIndex.unique]

        iN.pair <- unique(sapply(iULpIndex.cor[iIndex.unique], length))
        ls.pair <- vector(mode = "list", length = max(iN.pair))
        ls.pair[iN.pair] <- lapply(iN.pair, function(iN){.unorderedPairs(1:iN, distinct = TRUE)})
        
        ## **** contrast all pairs
        for(iC in iCluster2){ ## iC <- iCluster2[1]
            iCindex  <- index.cluster[[iC]][iULpIndex.cor[[iC]]]
            iCX.cor <- X.cor[iCindex,,drop=FALSE]
            iData <- data[iCindex,,drop=FALSE]

            if(NROW(iCX.cor)==1){ ## SAME LINEAR PREDICTOR FOR ALL OBSERVATIONS WITHIN CLUSTER
                if(heterogeneous){ ## make sure it is compatible with the next case (in case some cluster have only a single pair of observations)
                    iDF.diff <- as.data.frame(cbind("R",iCX.cor), drop = TRUE)
                }else{
                    iDF.diff <- as.data.frame(matrix(c("R",rep(iStrata,NCOL(iCX.cor))), nrow = 1, ncol = 1+NCOL(iCX.cor)))
                }
                iCode <- as.character(interaction(iDF.diff, drop=TRUE))
                ## name difference according to the covariate values
                if(length(attr(X.cor,"M.level"))==0 || identical(names(attr(X.cor,"M.level")),strata.var)){
                    if(n.strata>1){
                        iLevel <- paste0(":",U.strata[iStrata])
                    }else{
                        iLevel <- ""
                    }
                }else{
                    iLevel <- paste0("(",as.character(interaction(iData[,names(attr(X.cor,"M.level")),drop=FALSE],drop=TRUE)),")")
                    if(n.strata>1){
                        iLevel <- paste0(iLevel,":",U.strata[iStrata])
                    }
                }
                iName <- paste0("rho",iLevel)

                iOut <- data.frame(lp.x = lpnCluster.cor[[iC]][1],
                                   lp.y = lpnCluster.cor[[iC]][2],
                                   strata = iStrata,
                                   code = iCode,
                                   level = iLevel,
                                   param = iName,
                                   sigma = NA,
                                   k.x = NA,
                                   k.y = NA)

                if(!missing(X.var)){ ## add corresponding variance parameters
                    iCX.sigma <- iX.sigma[index.cluster[[iC]][1:2],,drop=FALSE]
                    iOut$sigma <- colnames(iCX.sigma)[colSums(iCX.sigma)!=0]

                    if(length(iX.k)>0){
                        iCX.k <- iX.k[index.cluster[[iC]][1:2],,drop=FALSE]
                        if(sum(iCX.k[1,,drop=FALSE]!=0)){
                            iOut$k.x <- colnames(iCX.k)[iCX.k[1,,drop=FALSE]==1]
                        }
                        if(sum(iCX.k[2,,drop=FALSE]!=0)){
                            iOut$k.y <- colnames(iCX.k)[iCX.k[2,,drop=FALSE]==1]
                        }
                    }
                }
                out[[iStrata]] <- rbind(out[[iStrata]], iOut)                   
                    
            }else{ ## DIFFERENT LINEAR PREDICTORS WITHIN CLUSTER
                iPair.time <- ls.pair[[NROW(iCX.cor)]]
                iM <- matrix(lpnCluster.cor[[iC]][iULpIndex.cor[[iC]]][iPair.time], ncol = 2, byrow = TRUE, dimnames = list(NULL, c("x","y")))
                
                iDF.diff <- as.data.frame(do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){ ## iCol <- 1

                    if(iM[iCol,"x"] < iM[iCol,"y"]){
                        ## make sure same correlation coefficient for (X=1,X=0) and (X=0,X=1)
                        iCX.cor1 <- iCX.cor[iPair.time[1,iCol],,drop=FALSE]
                        iCX.cor2 <- iCX.cor[iPair.time[2,iCol],,drop=FALSE]
                    }else if(iM[iCol,"x"] > iM[iCol,"y"]){
                        ## make sure same correlation coefficient for (X=1,X=0) and (X=0,X=1)
                        iCX.cor1 <- iCX.cor[iPair.time[2,iCol],,drop=FALSE]
                        iCX.cor2 <- iCX.cor[iPair.time[1,iCol],,drop=FALSE]
                    }else{
                        iCX.cor1 <- iCX.cor[min(iPair.time[,iCol]),,drop=FALSE]
                        iCX.cor2 <- iCX.cor[max(iPair.time[,iCol]),,drop=FALSE]
                    }
                    if(heterogeneous){
                        if(all(iCX.cor1==iCX.cor2)){return(cbind("R",iCX.cor1))}else{return(cbind(paste0("D",paste(iCX.cor1,collapse="")),iCX.cor2-iCX.cor1))}
                    }else{
                        if(all(iCX.cor1==iCX.cor2)){return(matrix(c("R",rep(iStrata,NCOL(iCX.cor1))), nrow = 1, ncol = 1+NCOL(iCX.cor1)))}else{return(matrix(c("D",as.numeric(iCX.cor2!=iCX.cor1)), nrow = 1))}
                    }
                })))
                iCode <- as.character(interaction(iDF.diff, drop=TRUE))
                ## name difference according to the covariate values
                iName.covcor <- setdiff(names(attr(X.cor,"M.level")),strata.var)
                iCov <- as.character(interaction(iData[,iName.covcor,drop=FALSE],drop=TRUE))

                index.iUCode <- which(!duplicated(iCode))
                iULevel <- stats::setNames(sapply(index.iUCode,function(iCol){
                    return(paste0("(",paste(unique(c(iCov[min(iPair.time[,iCol])],iCov[max(iPair.time[,iCol])])),collapse=","),")"))
                }), iCode[index.iUCode])
                iLevel <- iULevel[iCode]
                                
                if(n.strata>1){
                    iLevel <- paste0(iLevel,":",U.strata[iStrata])
                }
                iName <- paste0("rho",iLevel)

                index.unique <- which(!duplicated(cbind(iCode, iM)))
                iOut <- data.frame(lp.x = iM[index.unique,"x"],
                                   lp.y = iM[index.unique,"y"],
                                   strata = iStrata,
                                   code = iCode[index.unique],
                                   level = iLevel[index.unique],
                                   param = iName[index.unique],
                                   sigma = NA,
                                   k.x = NA,
                                   k.y = NA)

                if(!missing(X.var)){ ## add corresponding variance parameters
                    iCX.sigma <- iX.sigma[iCindex,,drop=FALSE]
                    iOut$sigma <- colnames(iCX.sigma)[colSums(iCX.sigma)!=0]

                    if(length(iX.k)>0){
                        iCX.k <- iX.k[iCindex,,drop=FALSE]
                        iNames.k <- colnames(iCX.k)
                        iName2 <- sapply(index.unique,function(iCol){ ## iCol <- 4
                            iName.x <- iNames.k[which(iCX.k[min(iPair.time[,iCol]),]==1)]
                            if(length(iName.x)==0){iName.x <- NA}
                            iName.y <- iNames.k[which(iCX.k[max(iPair.time[,iCol]),]==1)]
                            if(length(iName.y)==0){iName.y <- NA}
                            return(c(iName.x,iName.y))
                        })
                        iOut$k.x <- iName2[1,]
                        iOut$k.y <- iName2[2,]
                    }
                }
                out[[iStrata]] <- rbind(out[[iStrata]], iOut)
            
            }            
        }
        out[[iStrata]] <- out[[iStrata]][!duplicated(out[[iStrata]]),,drop=FALSE]
    }

    ## ** export
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    out <- out[order(out$lp.x,out$lp.y,out$strata),]
    attr(out, "levels") <- Ulp.cor
    return(out)
}

## ** .colnameOrder
## reorder the variable in the column name
.colnameOrder <- function(X, strata.var, n.strata){

    attr(X,"original.colnames") <- colnames(X)
    if(n.strata>1){
        attr(X,"ls.level") <- lapply(attr(X,"ls.level"), function(iL){
            iL[,c(setdiff(names(iL),strata.var),strata.var),drop=FALSE]
        })
        attr(X,"M.level") <- attr(X,"M.level")[,c(setdiff(names(attr(X,"M.level")),strata.var),strata.var),drop=FALSE]
        attr(X,"reference.level") <- attr(X,"reference.level")[,c(setdiff(names(attr(X,"M.level")),strata.var),strata.var),drop=FALSE]
        attr(X,"term.labels") <- unname(unlist(lapply(attr(X,"ls.level"), function(iL){paste(names(iL),collapse = ":")})))
        X.newname <- unname(sapply(attr(X,"ls.level"), function(iL){ ## iL <- attr(X,"ls.level")[[3]]
            iL[is.na(iL)] <- ""
            return(paste(paste0(names(iL),as.character(iL)),collapse = ":"))
        }))
        colnames(X) <- X.newname
    }
    
    return(X)    
}


##----------------------------------------------------------------------
### structure-skeleton.R ends here
