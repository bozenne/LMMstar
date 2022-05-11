### structure-skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: maj 11 2022 (13:47) 
##           By: Brice Ozenne
##     Update #: 1694
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

## * skeleton.IND
.skeleton.IND <- function(structure, data, indexData = NULL){

    ## ** prepare
    if(is.null(indexData)){
        indexData <- .extractIndexData(data = data, structure = structure)
    }
    time.var <- indexData$time.var
    strata.var <- indexData$strata.var
    U.time <- structure$U.time
    U.cluster <- structure$U.cluster
    U.strata <- structure$U.strata
    n.strata <- length(U.strata)
    index.clusterTime <- indexData$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    order.clusterTime <- indexData$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- indexData$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    X.var <- structure$X$var
    
    ## ** define parameters
    ## *** sigma
    outSigma <- .initSigma(X.var = X.var, strata.var = strata.var, U.strata = U.strata, n.strata = n.strata)
    X.var <- outSigma$X    
    param.sigma <- outSigma$param    
    strata.sigma <- outSigma$strata    

    ## *** k
    outK <- .initK(X.var = X.var, strata.var = strata.var, U.strata = U.strata, n.strata = n.strata,
                   param.sigma = param.sigma, strata.sigma = strata.sigma)
    X.var <- outK$X    
    param.k <- outK$param    
    strata.k <- outK$strata    
    name2 <- outK$name2

    ## ** gather parameters
    strata.param <- c(strata.sigma,strata.k)[colnames(X.var)]

    structure$param <- data.frame(name = c(param.sigma,param.k),
                                  strata = c(strata.sigma,strata.k),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("k",length=length(param.k))),
                                  name2 = as.character(NA),
                                  sigma = as.character(NA),
                                  k.x = as.character(NA),
                                  k.y = as.character(NA),
                                  stringsAsFactors = FALSE)
    structure$param$name2 <- if(is.null(name2)){structure$param$name}else{name2} ## sigma only if no k parameters
    if(any(structure$param$type=="k")){
        structure$param[structure$param$type=="k","sigma"] <- names(strata.sigma)[structure$param[structure$param$type=="k","strata"]]
    }

    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL

    ## ** export
    structure$X$var <- X.var
    return(structure)
}


## * skeleton.CS
.skeleton.CS <- function(structure, data, indexData = NULL){

    ## ** prepare
    if(is.null(indexData)){
        indexData <- .extractIndexData(data = data, structure = structure)
    }
    time.var <- indexData$time.var
    strata.var <- indexData$strata.var
    U.time <- structure$U.time
    U.cluster <- structure$U.cluster
    U.strata <- structure$U.strata
    n.strata <- length(U.strata)
    index.clusterTime <- indexData$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    index.clusterStrata <- indexData$index.clusterStrata ## vector of index relative to which strata each cluster belong to
    order.clusterTime <- indexData$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- indexData$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    X.var <- structure$X$var
    X.cor <- structure$X$cor

    ## ** param
    ## *** sigma
    outSigma <- .initSigma(X.var = X.var, strata.var = strata.var, U.strata = U.strata, n.strata = n.strata)
    X.var <- outSigma$X    
    param.sigma <- outSigma$param    
    strata.sigma <- outSigma$strata    

    ## *** k
    outK <- .initK(X.var = X.var, strata.var = strata.var, U.strata = U.strata, n.strata = n.strata,
                   param.sigma = param.sigma, strata.sigma = strata.sigma)
    X.var <- outK$X    
    param.k <- outK$param    
    strata.k <- outK$strata    
    name2 <- outK$name2

    ## *** param rho
    outRho <- .initRho(data = data, X.cor = X.cor, X.var = c(outK, list(strata.sigma = strata.sigma)), heterogeneous = structure$heterogeneous, 
                       U.cluster = U.cluster, index.cluster = index.cluster,
                       U.time = U.time, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime,
                       U.strata = U.strata, index.clusterStrata = index.clusterStrata, n.strata = n.strata)

    if(n.strata == 1 && length(outRho$param)==1){
        outRho$param[] <- "rho"
    }else if(n.strata == 1 && all(!duplicated(outRho$strata))){
        outRho$param[] <- paste0("rho",":",U.strata[outRho$strata])
    }
    param.rho <- unname(outRho$param)
    strata.rho <- stats::setNames(outRho$strata, param.rho)
    lpdiff.rho <- stats::setNames(outRho$code, param.rho)
    
    ## ** gather parameters
    strata.param <- c(c(strata.sigma,strata.k)[colnames(X.var)],strata.rho)
    structure$param <- data.frame(name = c(param.sigma,param.k,param.rho),
                                  strata = c(strata.sigma,strata.k,strata.rho),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("k",length=length(param.k)),rep("rho",length=length(param.rho))),
                                  name2 = as.character(NA),
                                  sigma = as.character(NA),
                                  k.x = as.character(NA),
                                  k.y = as.character(NA),                                  
                                  stringsAsFactors = FALSE)
    if(!is.null(name2)){
        structure$param$name2[structure$param$type %in% c("sigma","k")] <- name2
    }
    if(any(structure$param$type=="k")){
        structure$param[structure$param$type=="k","sigma"] <- names(strata.sigma)[structure$param[structure$param$type=="k","strata"]]
    }
    structure$param[structure$param$type=="rho","sigma"] <- outRho$sigma
    structure$param[structure$param$type=="rho","k.x"] <- outRho$k.x
    structure$param[structure$param$type=="rho","k.y"] <- outRho$k.y
    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL

    attr(structure$param, "table.rho") <- outRho

    ## ** export
    structure$X$var <- X.var
    structure$X$cor <- X.cor
    return(structure)
}

## * skeleton.UN
.skeleton.UN <- .skeleton.CS

## * skeleton.CUSTOM
.skeleton.CUSTOM <- function(structure, data, indexData = NULL){

    ## ** prepare
    if(is.null(indexData)){
        indexData <- .extractIndexData(data = data, structure = structure)
    }
    time.var <- indexData$time.var
    strata.var <- indexData$strata.var
    U.time <- structure$U.time
    U.cluster <- structure$U.cluster
    U.strata <- structure$U.strata
    n.strata <- length(U.strata)
    
    index.clusterTime <- indexData$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    index.clusterStrata <- indexData$index.clusterStrata ## vector of index relative to which strata each cluster belong to
    order.clusterTime <- indexData$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- indexData$index.cluster ## list of positions of the observation belonging to each cluster in the dataset
    cor.var <- structure$name$cor[[1]]

    X.var <- structure$X$var
    X.cor <- structure$X$cor

    ## ** param
    structure$param <- NULL
    if(!is.null(structure$FCT.sigma)){
        structure$param <- rbind(structure$param,
                                 data.frame(name = names(structure$init.sigma),
                                            strata = 1,
                                            type = "sigma",
                                            time = NA,
                                            name2 = as.character(NA)))
    }
    if(!is.null(structure$FCT.rho)){
        structure$param <- rbind(structure$param,
                                 data.frame(name = names(structure$init.rho),
                                            strata = 1,
                                            type = "rho",
                                            time = NA,
                                            name2 = as.character(NA)))
    }


    ## ** prepare for  patterns
    if(!is.null(X.cor)){
        outRho <- .initRho(data = data, X.cor = X.cor, heterogeneous = structure$heterogeneous, 
                           U.cluster = U.cluster, index.cluster = index.cluster,
                           U.time = U.time, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime,
                           U.strata = U.strata, index.clusterStrata = index.clusterStrata, n.strata = n.strata)
        browser()
        attr(structure$param,"lpdiff.rho") <- outRho$code        
    }

    ## ** pattern
    return(structure)
}

## * helpers
## ** .initSigma
.initSigma <- function(X.var, strata.var, U.strata, n.strata){

    if(n.strata==1){
        param.sigma <- "sigma"
        if(!identical(colnames(X.var)[attr(X.var,"assign")==0],"(Intercept)")){
            stop("Could not find the intercept in the design matrix for the variance.\n")
        }
        colnames(X.var)[attr(X.var,"assign")==0] <- param.sigma
        strata.sigma <- stats::setNames(1,param.sigma)
    }else{
        param.sigma <- paste0("sigma:",U.strata)
        if(!identical(colnames(X.var)[attr(X.var,"assign")==1],paste0(strata.var,U.strata))){
            stop("Could not find the strata-specific intercepts in the design matrix for the variance.\n")
        }
        colnames(X.var)[attr(X.var,"assign")==1] <- param.sigma
        strata.sigma <- stats::setNames(1:n.strata,param.sigma)
        attr(X.var,"order")[attr(X.var,"assign")==1] <- 0
    }


    return(list(X = X.var,
                param = param.sigma,
                strata = strata.sigma))

}

## ** .initK
.initK <- function(X.var, strata.var, U.strata, n.strata,
                   param.sigma, strata.sigma, sep = "_X_XX_X_"){

    index.k <- which(attr(X.var,"assign")>(n.strata>1))

    if(length(index.k)>0){
        M.level <- attr(X.var,"M.level")
        oldnames <- as.character(interaction(as.data.frame(lapply(1:NCOL(M.level), function(iVar){paste0(names(M.level)[iVar],M.level[index.k,iVar])})),sep=sep, drop = TRUE) )
        newnames <- as.character(interaction(M.level[index.k,,drop=FALSE],sep=sep, drop = TRUE))

        index.sigma <- which(attr(X.var,"assign")<=(n.strata>1))
        name2 <- rep(NA, length(param.sigma)+length(index.k))
        name2[index.sigma] <- as.character(interaction(M.level[index.sigma,,drop=FALSE],sep=sep,drop=TRUE))
        name2[index.k] <- newnames
        
        if(!identical(colnames(X.var)[index.k],oldnames)){
            stop("Could not find the k parameters in the design matrix for the variance.\n")
        }
        colnames(X.var)[index.k] <- paste("k",newnames,sep=".")
        param.k <- colnames(X.var)[index.k]

        if(n.strata == 1){
            strata.k <- stats::setNames(rep(1,length(param.k)),param.k)
        }else{
            strata.k <- stats::setNames(match(M.level[index.k,strata.var],U.strata),param.k)

            ## reorder by strata
            strata.sigmak <- stats::setNames(rep(NA,NCOL(X.var)),colnames(X.var))
            strata.sigmak[param.sigma] <- strata.sigma
            strata.sigmak[param.k] <- strata.k
            save.attr <- attributes(X.var)
            X.var <- X.var[,order(strata.sigmak),drop=FALSE]
            attr(X.var,"assign") <- save.attr$assign[order(strata.sigmak)]
            attr(X.var,"order") <- save.attr$order[order(strata.sigmak)]-1
            attr(X.var,"contrasts") <- save.attr$contrasts
            attr(X.var,"ls.level") <- save.attr$ls.level
            attr(X.var,"original.colnames") <- save.attr$original.colnames[order(strata.sigmak)]
            attr(X.var,"reference.level") <- save.attr$reference.level
            attr(X.var,"M.level") <- save.attr$M.level

            param.k <- param.k[order(strata.k)]
            name2[index.k] <- name2[index.k][order(strata.k)]
            strata.k <- strata.k[order(strata.k)]
        }
    }else{
        param.k <- NULL
        strata.k <- NULL
        name2 <- NULL
    }
    return(list(X = X.var,
                param = param.k,
                strata = strata.k,
                name2 = name2))
}

## ** .initRho
## for each cluster compute all pairwise difference in covariates to find the parameters
.initRho <- function(data, X.cor, X.var, heterogeneous, 
                     U.cluster, index.cluster,
                     U.time, index.clusterTime, order.clusterTime,
                     U.strata, index.clusterStrata, n.strata, sep = "_X_XX_X_"){

    fct.envir <- environment()
    n.time <- length(U.time)
        
    ## *** linear predictor for each observation and then grouped by cluster
    lpObs.cor <- interaction(as.data.frame(X.cor),sep=sep,drop=TRUE)
    Ulp.cor <- sort(levels(lpObs.cor)) 
    lpnObs.cor <- as.numeric(factor(lpObs.cor, levels = Ulp.cor))
    lpnCluster.cor <- stats::setNames(lapply(U.cluster, function(iC){
        lpnObs.cor[index.cluster[[iC]][order.clusterTime[[iC]]]]
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
            iOrder.clusterTime <- order.clusterTime
            iLpnCluster.cor <- lpnCluster.cor         
        }else{
            iCluster <- names(index.clusterStrata)[which(index.clusterStrata==iStrata)]
            iIndex.cluster <- index.cluster[iCluster]
            iIndex.clusterTime <- index.clusterTime[iCluster]
            iOrder.clusterTime <- order.clusterTime[iCluster]
            iLpnCluster.cor <- lpnCluster.cor[iCluster]            
        }
        if(!missing(X.var)){
            iX.sigma <- X.var$X[,names(which(X.var$strata.sigma==iStrata)),drop=FALSE]
            iX.k <- X.var$X[,names(which(X.var$strata==iStrata)),drop=FALSE]
        }
        ## (re-order by the number of observed timepoints)
        iCluster <- iCluster[order(sapply(iIndex.clusterTime,length) - sapply(iIndex.clusterTime,min)/(n.time+1), decreasing = TRUE)]
        
        ## no replicates
        if(all(sapply(iLpnCluster.cor,length)<=1)){
            next
        }

        ## **** find unique linear predictor values within cluster
        iULpIndex.cor <- setNames(vector(mode = "list", length = length(iLpnCluster.cor)), iCluster)
        
        iULpCluster.cor <- stats::setNames(lapply(names(iLpnCluster.cor), function(iId){ ## iId <- names(iLpnCluster.cor)[1]
            iTest <- which(!duplicated(iLpnCluster.cor[[iId]]))
            iCode <- iLpnCluster.cor[[iId]][iTest]
            fct.envir$iULpIndex.cor[[iId]] <- iTest
            return(iCode)
        }), iCluster)

        ## **** remove cluster identical to other cluster in term of linear predictors
        iIndex.unique <- which(!duplicated(iULpCluster.cor))
        iCluster2 <- iCluster[iIndex.unique]
        iULpIndex.cor2 <- iULpIndex.cor[iIndex.unique]
        iULpCluster.cor2 <- iULpCluster.cor[iIndex.unique]

        iN.pair <- unique(sapply(iULpIndex.cor2, length))
        ls.pair <- vector(mode = "list", length = max(iN.pair))
        ls.pair[iN.pair] <- lapply(iN.pair, function(iN){.unorderedPairs(1:iN, distinct = TRUE)})
        
        ## **** contrast all pairs
        for(iC in iCluster2){ ## iC <- iCluster2[1]
            iCindex  <- index.cluster[[iC]][order.clusterTime[[iC]]][iULpIndex.cor2[[iC]]]
            iCX.cor <- X.cor[iCindex,,drop=FALSE]
            iData <- data[iCindex,,drop=FALSE]
                
            if(NROW(iCX.cor)==1){ ## SAME LINEAR PREDICTOR FOR ALL OBSERVATIONS WITHIN CLUSTER
                iCode <- paste0("R",as.character(interaction(as.data.frame(iCX.cor), drop = TRUE)))
                       
                if(is.null(out[[iStrata]]) || iCode %in% out[[iStrata]]$code == FALSE){
                    
                    ## name difference according to the covariate values
                    if(length(attr(X.cor,"M.level"))==0){
                        iName <- "rho"
                    }else{
                        iName <- paste0("rho(",as.character(interaction(iData[,names(attr(X.cor,"M.level")),drop=FALSE],drop=TRUE)),")")
                    }
                    if(n.strata>1){
                        iName <- paste0(iName,":",U.strata[iStrata])
                    }
                    iOut <- data.frame(lp.x = lpnCluster.cor[[iC]][1],
                                       lp.y = lpnCluster.cor[[iC]][2],
                                       strata = iStrata,
                                       code = iCode,
                                       param = iName,
                                       sigma = NA,
                                       k.x = NA,
                                       k.y = NA)

                    if(!missing(X.var)){ ## add corresponding variance parameters
                        iCX.sigma <- iX.sigma[index.cluster[[iC]][order.clusterTime[[iC]]][1:2],,drop=FALSE]
                        iOut$sigma <- colnames(iCX.sigma)[colSums(iCX.sigma)!=0]

                        if(length(iX.k)>0){
                            iCX.k <- iX.k[index.cluster[[iC]][order.clusterTime[[iC]]][1:2],,drop=FALSE]
                            if(sum(iCX.k[1,,drop=FALSE]!=0)){
                                iOut$k.x <- colnames(iCX.k)[iCX.k[1,,drop=FALSE]==1]
                            }
                            if(sum(iCX.k[2,,drop=FALSE]!=0)){
                                iOut$k.y <- colnames(iCX.k)[iCX.k[2,,drop=FALSE]==1]
                            }
                        }
                    }
                    
                    out[[iStrata]] <- rbind(out[[iStrata]], iOut)                   
                    
                }
            }else{ ## DIFFERENT LINEAR PREDICTORS WITHIN CLUSTER
                iPair.time <- ls.pair[[NROW(iCX.cor)]]
                iDF.diff <- as.data.frame(do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){ ## iCol <- 1
                    iCX.cor1 <- iCX.cor[min(iPair.time[,iCol]),,drop=FALSE]
                    iCX.cor2 <- iCX.cor[max(iPair.time[,iCol]),,drop=FALSE]
                    if(heterogeneous){
                        if(all(iCX.cor1==iCX.cor2)){return(cbind("R",iCX.cor1))}else{return(cbind(paste0("D",paste(iCX.cor1,collapse="")),iCX.cor2-iCX.cor1))}
                    }else{
                        if(all(iCX.cor1==iCX.cor2)){return(matrix(c("R",""), nrow = 1, ncol = 2))}else{return(cbind("D",sum(iCX.cor2!=iCX.cor1)))}
                    }
                })))

                iCode <- as.character(interaction(iDF.diff, drop=TRUE))
                iIndex.store <- which(iCode %in% out[[iStrata]]$code == FALSE)
                if(length(iIndex.store)>0){

                    ## name difference according to the covariate values
                    iCov <- as.character(interaction(iData[,names(attr(X.cor,"M.level")),drop=FALSE],drop=TRUE))
                    iName <- sapply(1:NCOL(iPair.time),function(iCol){
                            return(paste0("rho(",paste(unique(c(iCov[min(iPair.time[,iCol])],iCov[max(iPair.time[,iCol])])),collapse=","),")"))
                    })
                    if(n.strata>1){
                        iName <- paste0(iName,":",U.strata[iStrata])
                    }
                    iM <- matrix(lpnCluster.cor[[iC]][iPair.time], ncol = 2, byrow = FALSE, dimnames = list(NULL, c("x","y")))[iIndex.store,]
                    iOut <- data.frame(lp.x = iM[,"x"],
                                       lp.y = iM[,"y"],
                                       strata = rep(iStrata, length(iIndex.store)),
                                       code = iCode[iIndex.store],
                                       param = iName[iIndex.store],
                                       sigma = NA,
                                       k.x = NA,
                                       k.y = NA)

                    if(!missing(X.var)){ ## add corresponding variance parameters
                        iCX.sigma <- iX.sigma[iCindex,,drop=FALSE]
                        iOut$sigma <- colnames(iCX.sigma)[colSums(iCX.sigma)!=0]

                        if(length(iX.k)>0){
                            iCX.k <- iX.k[iCindex,,drop=FALSE]

                            iName2 <- sapply(1:NCOL(iPair.time),function(iCol){ ## iCol <- 1
                                iName.x <- names(which(iCX.k[min(iPair.time[,iCol]),]==1))
                                if(length(iName.x)==0){iName.x <- NA}
                                iName.y <- names(which(iCX.k[max(iPair.time[,iCol]),]==1))
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
            
       }
        
    }

    ## ** export
    out <- do.call(rbind, out)
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
        colnames(X) <- unname(sapply(attr(X,"ls.level"), function(iL){paste(paste0(names(iL),as.character(iL)),collapse = ":")}))
    }
    
    return(X)    
}


##----------------------------------------------------------------------
### structure-skeleton.R ends here
