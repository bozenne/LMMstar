### structure-skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: apr 13 2022 (18:40) 
##           By: Brice Ozenne
##     Update #: 1476
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
                                  time = NA,
                                  stringsAsFactors = FALSE)
    structure$param$name2 <- if(is.null(name2)){structure$param$name}else{name2} ## sigma only if no k parameters

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
    if(is.null(X.cor)){ ## no repetition
        param.rho <- NULL
        strata.rho <- NULL
        lpdiff.rho <- NULL
    }else if(n.strata==1){ ## 1 strata with repetition
        if(is.na(structure$name$cor)){ ## no covariate
            param.rho <- "rho"
            strata.rho <- stats::setNames(1,param.rho)
            colnames(X.cor) <- param.rho
            lpdiff.rho <- stats::setNames(paste0("R.",interaction(as.data.frame(unique(X.cor)),drop=FALSE)),param.rho)
        }else{ ## covariates
            outBuild <- .buildRho(data = data, X.cor = X.cor, heterogeneous = structure$heterogeneous, 
                                  U.cluster = U.cluster, index.cluster = index.cluster, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime)
            param.rho <- outBuild$param 
            strata.rho <- outBuild$strata 
            lpdiff.rho <- outBuild$code 
        }        
    }else{ ## some strata with repetition (only relevant when using gls optimizer otherwise handle as covariate without strata)
        param.rho <- NULL
        strata.rho <- NULL
        lpdiff.rho <- NULL
        index.rho <- NULL
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            iIndex <- which(data$XXstrata.indexXX == iStrata)
            if(any(duplicated(data[iIndex,"XXclusterXX"]))){ ## check if strata has repetition
                iParam <- paste0("rho:",U.strata[iStrata])

                index.rho <- c(index.rho,iIndex[1])
                param.rho <- c(param.rho,iParam)
                strata.rho <- c(strata.rho,stats::setNames(iStrata,iParam))
                colnames(X.cor)[iStrata] <- iParam
            }
        }
        lpdiff.rho <- stats::setNames(paste0("R.",interaction(as.data.frame(X.cor[index.rho,,drop=FALSE]),drop=FALSE)),param.rho)
    }

    ## ** gather parameters
    strata.param <- c(c(strata.sigma,strata.k)[colnames(X.var)],strata.rho)
    structure$param <- data.frame(name = c(param.sigma,param.k,param.rho),
                                  strata = c(strata.sigma,strata.k,strata.rho),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("k",length=length(param.k)),rep("rho",length=length(param.rho))),
                                  time = NA,
                                  stringsAsFactors = FALSE)

    structure$param$name2 <- as.character(NA)
    if(!is.null(name2)){
        structure$param$name2[structure$param$type %in% c("sigma","k")] <- name2
    }
    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL
    attr(structure$param,"lpdiff.rho") <- lpdiff.rho

    ## ** export
    structure$X$var <- X.var
    structure$X$cor <- X.cor
    return(structure)
}

## * skeleton.UN
.skeleton.UN <- function(structure, data, indexData = NULL){

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
    cor.var <- structure$name$cor[[1]]

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
    if(is.null(X.cor)){
        lpdiff.rho <- NULL
        param.rho <- NULL
        strata.rho <- NULL
    }else if(n.strata==1){

        outBuild <- .buildRho(data = data, X.cor = X.cor, heterogeneous = TRUE,
                              U.cluster = U.cluster, index.cluster = index.cluster, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime)
        param.rho <- outBuild$param 
        strata.rho <- outBuild$strata 
        lpdiff.rho <- outBuild$code 

    }else{
        pair.time <- .unorderedPairs(1:length(U.time), distinct = TRUE)
        pair.name <- apply(pair.time,2,function(iT){paste0("(",U.time[iT[1]],",",U.time[iT[2]],")")})
        param.rho <- unlist(lapply(U.strata, function(iStrata){paste0(paste0("rho",pair.name),":",iStrata)}))
        strata.rho <- stats::setNames(unlist(lapply(1:n.strata, function(iStrata){rep(iStrata, NCOL(pair.time))})), param.rho)
        time.rho <- do.call(cbind,lapply(1:n.strata, function(iStrata){pair.time}))
        colnames(time.rho) <- param.rho

        Data.rho <- do.call(rbind,lapply(U.strata, function(iS){ ## iS <- U.strata[1]
            iIndex <- which(attr(X.cor,"M.level")[,strata.var]==iS)
            iIndex.order <- iIndex[match(attr(X.cor,"M.level")[iIndex,cor.var], U.time)]
            iData <- attr(X.cor,"M.level")[iIndex.order,,drop=FALSE]
            iData[[strata.var]] <- factor(iData[[strata.var]], levels = U.strata)
            iData[[cor.var]] <- factor(iData[[cor.var]], levels = U.time)
            return(iData)
        }))

        UX.cor <- .colnameOrder(.model.matrix_regularize(structure$formula$cor, data = Data.rho, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)[,colnames(X.cor),drop=FALSE]

        lpdiff.rho <- stats::setNames(unlist(lapply(U.strata, function(iS){ ## iS <- U.strata[1]
            iUX.cor <- UX.cor[Data.rho[[strata.var]]==iS,,drop=FALSE]
            if(NROW(iUX.cor)==1){return(NULL)}
            iLPdiff.rho <- apply(pair.time,2,function(iCol){paste0("D.",paste(iUX.cor[max(iCol),]-iUX.cor[min(iCol),],collapse="."))})
            return(iLPdiff.rho)
        })),param.rho)
    }

    ## ** gather parameters
    strata.param <- c(c(strata.sigma,strata.k)[colnames(X.var)],strata.rho)
    structure$param <- data.frame(name = c(param.sigma,param.k,param.rho),
                                  strata = c(strata.sigma,strata.k,strata.rho),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("k",length=length(param.k)),rep("rho",length=length(param.rho))),
                                  time = NA,
                                  stringsAsFactors = FALSE)
    structure$param$name2 <- as.character(NA)
    if(!is.null(name2)){
        structure$param$name2[structure$param$type %in% c("sigma","k")] <- name2
    }

    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL
    attr(structure$param,"lpdiff.rho") <- lpdiff.rho

    ## ** export
    structure$X$var <- X.var
    structure$X$cor <- X.cor
    return(structure)
}

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
        if(n.strata==1){
            outBuild <- .buildRho(data = data, X.cor = X.cor, heterogeneous = TRUE,
                                  U.cluster = U.cluster, index.cluster = index.cluster, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime)
            attr(structure$param,"lpdiff.rho") <- outBuild$code
        }else{
            stop("CUSTOM structure does not currently handle stratification. \n")
        }
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
                   param.sigma, strata.sigma){

    index.k <- which(attr(X.var,"assign")>(n.strata>1))

    if(length(index.k)>0){
        M.level <- attr(X.var,"M.level")
        oldnames <- as.character(interaction(as.data.frame(lapply(1:NCOL(M.level), function(iVar){paste0(names(M.level)[iVar],M.level[index.k,iVar])})),sep=":", drop = TRUE) )
        newnames <- as.character(interaction(M.level[index.k,,drop=FALSE],sep=":", drop = TRUE))

        index.sigma <- which(attr(X.var,"assign")<=(n.strata>1))
        name2 <- rep(NA, length(param.sigma)+length(index.k))
        name2[index.sigma] <- as.character(interaction(M.level[index.sigma,,drop=FALSE],sep=":",drop=TRUE))
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



## ** .buildRho
## identify all possible pairwise combinations 
.buildRho <- function(data, X.cor, heterogeneous,
                      U.cluster, index.cluster, index.clusterTime, order.clusterTime){

    lp.cor <- as.character(interaction(as.data.frame(X.cor),drop=TRUE))
    LPcluster.cor <- lapply(U.cluster, function(iC){
        lp.cor[index.cluster[[iC]][order.clusterTime[[iC]]]]
    })
    indexCluster.cor <- which(!duplicated(LPcluster.cor))

    seq.length <- sort(unique(sapply(index.cluster,length)))
    ls.pair.time <- lapply(seq.length, function(iN){.unorderedPairs(1:iN, distinct = TRUE)})

    ## for each cluster compute all pairwise difference in covariates
    all.lpdiff.rho <- unlist(lapply(U.cluster[indexCluster.cor], function(iC){ ## iC <- U.cluster[1]
        iX <- X.cor[index.cluster[[iC]][order.clusterTime[[iC]]],,drop=FALSE]
        if(NROW(iX)==1){return(NULL)}
        iPair.time <- ls.pair.time[[which(NROW(iX)==seq.length)]]
        iDF.diff <- as.data.frame(do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){ ## iCol <- 1
            iX1 <- iX[min(iPair.time[,iCol]),,drop=FALSE]
            iX2 <- iX[max(iPair.time[,iCol]),,drop=FALSE]
            if(heterogeneous){
                if(all(iX1==iX2)){return(cbind("R",iX1))}else{return(cbind("D",iX2-iX1))}
            }else{
                if(all(iX1==iX2)){return(matrix(c("R",""), nrow = 1, ncol = 2))}else{return(cbind("D",sum(iX2!=iX1)))}
            }
        })))
        iVec.diff <- as.character(interaction(iDF.diff, drop=TRUE))
        if(length(names(attr(X.cor,"M.level")))>0){ ## name difference according to the covariate values
            iCov <- as.character(interaction(data[index.cluster[[iC]][order.clusterTime[[iC]]],names(attr(X.cor,"M.level")),drop=FALSE],drop=TRUE))
            names(iVec.diff) <- sapply(1:NCOL(iPair.time),function(iCol){paste0("(",iCov[min(iPair.time[,iCol])],",",iCov[max(iPair.time[,iCol])],")")})
        }
        return(iVec.diff)
    }))

    test.duplicated <- duplicated(all.lpdiff.rho)
    lpdiff.rho <- all.lpdiff.rho[test.duplicated==FALSE]
    names(lpdiff.rho) <- paste0("rho",names(lpdiff.rho))
    param.rho <- names(lpdiff.rho)
    strata.rho <- stats::setNames(rep(1,length(param.rho)),param.rho)

    time.rho <- do.call(cbind,lapply(U.cluster[indexCluster.cor], function(iC){ ## iC <- "105"
        .unorderedPairs(index.clusterTime[[iC]][order.clusterTime[[iC]]], distinct = TRUE)
    }))[,test.duplicated==FALSE,drop=FALSE]

    reorder.rho <- order(time.rho[1,],time.rho[2,])

    out <- list(time = time.rho[,reorder.rho,drop=FALSE],
                code = lpdiff.rho[reorder.rho],
                param = param.rho[reorder.rho],
                strata = strata.rho[reorder.rho])

    return(out)
}
##----------------------------------------------------------------------
### structure-skeleton.R ends here
