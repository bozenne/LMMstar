### skeletonStructure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: sep 15 2021 (18:40) 
##           By: Brice Ozenne
##     Update #: 530
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * skeletonStructure
##' @title Parametrization of Covariance Structure
##' @description Parametrization of Covariance Structure
##'
##' @param structure structure
##' @param data dataset
##'
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$gender <- c("M","F")[as.numeric(gastricbypassL$id) %% 2+1]
##' dd <- gastricbypassL[!duplicated(gastricbypassL[,c("time","gender")]),]
##' 
##' ## independence
##' skeletonStructure(IND(~1, var.time = "time"), data = dd)
##' skeletonStructure(IND(~time), data = dd)
##' skeletonStructure(IND(~time|id), data = dd)
##' skeletonStructure(IND(~time+gender|id, var.time = "time"), data = dd)
##' 
##' ## compound symmetry
##' skeletonStructure(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' skeletonStructure(CS(~time|id), data = gastricbypassL)
##' skeletonStructure(CS(gender~time|id), data = gastricbypassL)
##' 
##' ## unstructured
##' skeletonStructure(UN(~visit|id), data = gastricbypassL)
##' skeletonStructure(UN(~time|id), data = gastricbypassL)
##' skeletonStructure(UN(~visit+gender|id, var.time = "time"), data = gastricbypassL)
##' skeletonStructure(UN(gender~visit|id), data = gastricbypassL)
##' 
##' @export
`skeletonStructure` <-
    function(structure, data) UseMethod("skeletonStructure")

## * skeletonStructure.IND
##' @export
skeletonStructure.IND <- function(structure, data){

    ## ** prepare
    outInit <- .initSkeleton(data = data, structure = structure)
    data <- outInit$data ## convert variables to factor 
    time.var <- outInit$time.var
    structure$U.time <- outInit$U.time
    n.time <- length(structure$U.time)
    cluster.var <- outInit$cluster.var
    structure$U.cluster <- outInit$U.cluster
    strata.var <- outInit$strata.var
    structure$U.strata <- outInit$U.strata
    index.clusterTime <- outInit$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    order.clusterTime <- outInit$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- outInit$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    formula.var <- structure$formula$var
    X.var <- model.matrix_regularize(formula.var, data = data)

    ## ** param
    ## *** sigma
    param.sigma <- "sigma"
    if(!identical(colnames(X.var)[attr(X.var,"assign")==0],"(Intercept)")){
        stop("Could not find the intercept in the design matrix for the variance.\n")
    }
    colnames(X.var)[attr(X.var,"assign")==0] <- param.sigma
    strata.sigma <- stats::setNames(1,param.sigma)
        
    ## *** k
    index.k <- which(attr(X.var,"assign")>0)
    if(length(index.k)>0){
        colnames(X.var)[index.k] <- paste("k",colnames(X.var)[index.k],sep=".")
        param.k <- colnames(X.var)[index.k]
        strata.k <- stats::setNames(rep(1,length(param.k)),param.k)
    }else{
        param.k <- NULL
        strata.k <- NULL
    }

    strata.param <- c(strata.sigma,strata.k)[colnames(X.var)]

    structure$param <- data.frame(name = c(param.sigma,param.k),
                                  strata = c(strata.sigma,strata.k),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("k",length=length(param.k))),
                                  time = NA)
    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL

    ## ** pattern
    outPattern <- .patternStructure(X.var = X.var, X.cor = NULL, lpdiff.rho = NULL, data = data,
                                    time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                    index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                    strata.var = strata.var, strata.param = setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)

    structure$param$time <- outPattern$time.param
    structure$X <- list(Upattern = outPattern$Upattern,
                        pattern.cluster = outPattern$pattern.cluster,
                        var = outPattern$X.var,
                        cor = outPattern$X.cor)

    ## ** export
    return(structure)
}


## * skeletonStructure.CS
##' @export
skeletonStructure.CS <- function(structure, data){

    ## ** prepare
    outInit <- .initSkeleton(data = data, structure = structure)
    data <- outInit$data ## convert variables to factor 
    time.var <- outInit$time.var
    structure$U.time <- outInit$U.time
    n.time <- length(structure$U.time)
    cluster.var <- outInit$cluster.var
    structure$U.cluster <- outInit$U.cluster
    strata.var <- outInit$strata.var
    structure$U.strata <- outInit$U.strata
    n.strata <- length(structure$U.strata)
    index.clusterTime <- outInit$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    order.clusterTime <- outInit$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- outInit$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    formula.var <- structure$formula$var
    formula.cor <- structure$formula$cor
    
    X.var <- model.matrix_regularize(formula.var, data = data)
    X.cor <- model.matrix_regularize(formula.cor, data = data)
    
    ## ** param
    ## *** sigma
    if(n.strata==1){
        param.sigma <- "sigma"
        if(!identical(colnames(X.var),"(Intercept)")){
            stop("Could not find the intercept in the design matrix for the variance.\n")
        }
        colnames(X.var) <- param.sigma
        strata.sigma <- stats::setNames(1,param.sigma)
    }else{
        param.sigma <- paste0("sigma:",structure$U.strata)
        if(!identical(colnames(X.var),paste0(strata.var,structure$U.strata))){
            stop("Could not find the strata-specific intercepts in the design matrix for the variance.\n")
        }
        colnames(X.var) <- param.sigma
        strata.sigma <- stats::setNames(1:n.strata,param.sigma)
    }

    ## *** param rho
    if(n.time==1){
        param.rho <- NULL
        strata.rho <- NULL
        X.cor <- NULL
    }else if(n.strata==1){
        param.rho <- "rho"
        strata.rho <- stats::setNames(1,param.rho)
        colnames(X.cor) <- param.rho
        lpdiff.rho <- setNames(as.character(interaction(as.data.frame(unique(X.cor)),drop=FALSE)),param.rho)
    }else{
        param.rho <- paste0("rho:",structure$U.strata)
        colnames(X.cor) <- param.rho
        strata.rho <- stats::setNames(1:n.strata,param.rho)
        lpdiff.rho <- setNames(as.character(sort(interaction(as.data.frame(unique(X.cor)),drop=FALSE))),param.rho)
    }
    strata.param <- c(strata.sigma,strata.rho)[c(colnames(X.var),colnames(X.cor))]

    structure$param <- data.frame(name = c(param.sigma,param.rho),
                                  strata = c(strata.sigma,strata.rho),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("rho",length=length(param.rho))),
                                  time = NA)
    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL

    ## ** pattern
    outPattern <- .patternStructure(X.var = X.var, X.cor = X.cor, lpdiff.rho = lpdiff.rho, data = data,
                                    time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                    index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                    strata.var = strata.var, strata.param = setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)
    structure$param$time <- outPattern$time.param
    structure$X <- list(Upattern = outPattern$Upattern,
                        pattern.cluster = outPattern$pattern.cluster,
                        var = outPattern$X.var,
                        cor = outPattern$X.cor)

    ## ** export
    return(structure)
}

## * skeletonStructure.UN
##' @export
skeletonStructure.UN <- function(structure, data){

    ## ** prepare
    outInit <- .initSkeleton(data = data, structure = structure)
    data <- outInit$data ## convert variables to factor 
    time.var <- outInit$time.var
    structure$U.time <- outInit$U.time
    n.time <- length(structure$U.time)
    cluster.var <- outInit$cluster.var
    structure$U.cluster <- outInit$U.cluster
    strata.var <- outInit$strata.var
    structure$U.strata <- outInit$U.strata
    n.strata <- length(structure$U.strata)
    index.clusterTime <- outInit$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    order.clusterTime <- outInit$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- outInit$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    formula.var <- structure$formula$var
    formula.cor <- structure$formula$cor

    X.var <- model.matrix_regularize(formula.var, data = data, augmodel = TRUE)
    X.cor <- model.matrix_regularize(formula.cor, data = data, augmodel = TRUE)

    ## ** param
    ## *** sigma
    if(n.strata==1){
        param.sigma <- "sigma"
        if(!identical(colnames(X.var)[attr(X.var,"assign")==0],"(Intercept)")){
            stop("Could not find the intercept in the design matrix for the variance.\n")
        }
        colnames(X.var)[attr(X.var,"assign")==0] <- param.sigma
        strata.sigma <- stats::setNames(1,param.sigma)
    }else{
        param.sigma <- paste0("sigma:",structure$U.strata)
        if(!identical(colnames(X.var)[attr(X.var,"assign")==1],paste0(strata.var,structure$U.strata))){
            stop("Could not find the strata-specific intercepts in the design matrix for the variance.\n")
        }
        colnames(X.var)[attr(X.var,"assign")==1] <- param.sigma
        strata.sigma <- stats::setNames(1:n.strata,param.sigma)
    }

    ## *** k
    index.k <- which(attr(X.var,"assign")>(n.strata>1))
    if(length(index.k)>0){
        colnames(X.var)[index.k] <- paste("k",colnames(X.var)[index.k],sep=".")
        param.k <- colnames(X.var)[index.k]
        if(n.strata == 1){
            strata.k <- stats::setNames(rep(1,length(param.k)),param.k)
        }else{
            strata.k <- stats::setNames(match(attr(X.var,"M.level")[index.k,strata.var],structure$U.strata),param.k)

            ## reorder by strata
            strata.sigmak <- setNames(rep(NA,NCOL(X.var)),colnames(X.var))
            strata.sigmak[param.sigma] <- strata.sigma
            strata.sigmak[param.k] <- strata.k
            X.var <- X.var[,order(strata.sigmak),drop=FALSE]
            param.k <- param.k[order(strata.k)]
            strata.k <- strata.k[param.k]
        }
    }else{
        param.k <- NULL
        strata.k <- NULL
    }

    ## *** param rho
    if(n.time==1){
        pairn.rho <- NULL
        lpdiff.rho <- NULL
        param.rho <- NULL
        strata.rho <- NULL
        X.cor <- NULL
    }else if(n.strata==1){
        lp.cor <- as.character(interaction(as.data.frame(X.cor),drop=TRUE))
        LPcluster.cor <- lapply(structure$U.cluster, function(iC){
            lp.cor[index.cluster[[iC]][order.clusterTime[[iC]]]]
        })
        indexCluster.cor <- which(!duplicated(LPcluster.cor))

        ## for each cluster compute all pairwise difference in covariates
        all.lpdiff.rho <- unlist(lapply(structure$U.cluster[indexCluster.cor], function(iC){ ## iC <- 1
            iX <- X.cor[index.cluster[[iC]][order.clusterTime[[iC]]],,drop=FALSE]
            iPair.time <- .unorderedPairs(1:NROW(iX), distinct = TRUE)
            iDF.diff <- as.data.frame(do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){
                iX1 <- iX[min(iPair.time[,iCol]),,drop=FALSE]
                iX2 <- iX[max(iPair.time[,iCol]),,drop=FALSE]
                if(all(iX1==iX2)){return(iX1)}else{return(iX2-iX1)}
            })))
            iVec.diff <- as.character(interaction(iDF.diff, drop=TRUE))
            iCov <- as.character(interaction(data[index.cluster[[iC]][order.clusterTime[[iC]]],names(attr(X.cor,"M.level")),drop=FALSE],drop=TRUE))
            names(iVec.diff) <- sapply(1:NCOL(iPair.time),function(iCol){paste0("(",iCov[min(iPair.time[,iCol])],",",iCov[max(iPair.time[,iCol])],")")})
            return(iVec.diff)
        }))
        test.duplicated <- duplicated(all.lpdiff.rho)
        lpdiff.rho <- all.lpdiff.rho[test.duplicated==FALSE]
        names(lpdiff.rho) <- paste0("rho",names(lpdiff.rho))

        param.rho <- names(lpdiff.rho)
        strata.rho <- stats::setNames(rep(1,length(param.rho)),param.rho)
        time.rho <- do.call(cbind,lapply(structure$U.cluster[indexCluster.cor], function(iC){
            .unorderedPairs(index.clusterTime[[iC]], distinct = TRUE)
        }))[,test.duplicated==FALSE]
    }else{
        pair.time <- .unorderedPairs(1:length(structure$U.time), distinct = TRUE)
        pair.name <- apply(pair.time,2,function(iT){paste0("(",structure$U.time[iT[1]],",",structure$U.time[iT[2]],")")})
        param.rho <- unlist(lapply(structure$U.strata, function(iStrata){paste0(paste0("rho",pair.name),":",iStrata)}))
        strata.rho <- setNames(unlist(lapply(1:n.strata, function(iStrata){rep(iStrata, NCOL(pair.time))})), param.rho)
        time.rho <- do.call(cbind,lapply(1:n.strata, function(iStrata){pair.time}))
        colnames(time.rho) <- param.rho

        lpdiff.rho <- stats::setNames(unlist(lapply(structure$U.strata, function(iS){ ## iS <- structure$U.strata[1]
            iIndex <- which(attr(X.cor,"M.level")[,strata.var]==iS)
            iIndex.order <- iIndex[match(attr(X.cor,"M.level")[iIndex,time.var], structure$U.time)]
            iData <- attr(X.cor,"M.level")[iIndex.order,,drop=FALSE]
            iData[[strata.var]] <- factor(iData[[strata.var]], levels = structure$U.strata)
            iData[[time.var]] <- factor(iData[[time.var]], levels = structure$U.time)
            iUX.cor <- model.matrix(structure$formula$cor, data = iData)[,colnames(X.cor),drop=FALSE]
            iLPdiff.rho <- apply(pair.time,2,function(iCol){paste(iUX.cor[max(iCol),]-iUX.cor[min(iCol),],collapse=".")})
            return(iLPdiff.rho)
        })),param.rho)
    }

    strata.param <- c(c(strata.sigma,strata.k)[colnames(X.var)],strata.rho)

    structure$param <- data.frame(name = c(param.sigma,param.k,param.rho),
                                  strata = c(strata.sigma,strata.k,strata.rho),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("k",length=length(param.k)),rep("rho",length=length(param.rho))),
                                  time = NA)
    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL

    ## ** pattern
    outPattern <- .patternStructure(X.var = X.var, X.cor = X.cor, lpdiff.rho = lpdiff.rho, data = data,
                                    time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                    index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                    strata.var = strata.var, strata.param = setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)
    structure$param$time <- outPattern$time.param
    structure$X <- list(Upattern = outPattern$Upattern,
                        pattern.cluster = outPattern$pattern.cluster,
                        var = outPattern$X.var,
                        cor = outPattern$X.cor)

    ## ** export
    return(structure)
}

## * helpers
## ** .initSkeleton
.initSkeleton <- function(data, structure){
    
    ## *** find variable names
    time.var <- structure$name$time
    cluster.var <- structure$name$cluster
    strata.var <- structure$name$strata
    all.var <- na.omit(c(time.var,cluster.var,strata.var,unlist(lapply(structure$formula,all.vars))))
    
    if(length(all.var)>0){
        for(iVar in all.var){
            if(is.factor(data[[iVar]])){
                data[[iVar]] <- droplevels(data[[iVar]])
            }else{
                data[[iVar]] <- as.factor(data[[iVar]])
            }
        }
    }

    ## *** find unique levels of each variable
    if(!is.na(cluster.var)){
        if(is.factor(data[[cluster.var]])){
            U.cluster <- levels(data[[cluster.var]])
        }else{
            U.cluster <- as.character(unique(data[[cluster.var]]))
        }
    }else{
        U.cluster <- paste0("id",1:NROW(data))
    }
    
    if(!is.na(time.var)){
        if(is.factor(data[[time.var]])){
            U.time <- levels(data[[time.var]])
        }else{
            U.time <- as.character(unique(data[[time.var]]))
        }
    }else if(!is.na(cluster.var)){
        U.time <- paste0("T",1:max(table(data[[cluster.var]])))
    }else{
        U.time <- "T1"
    }

    if(!is.na(strata.var)){
        if(is.factor(data[[strata.var]])){
            U.strata <- levels(data[[strata.var]])
        }else{
            U.strata <- as.character(unique(data[[strata.var]]))
        }
    }else{
        U.strata <- "S1"
    }

    ## *** find position of each cluster
    if(!is.na(cluster.var)){
        index.cluster <- tapply(1:NROW(data),data[[cluster.var]],function(iI){iI})
    }else{
        index.cluster <- setNames(as.list(1:NROW(data)),U.cluster)
    }

    ## *** find time corresponding to each cluster
    if(is.na(cluster.var) && is.na(time.var)){
        index.clusterTime <- setNames(as.list(rep(1,NROW(data))),U.cluster)
    }else if(!is.na(cluster.var) && !is.na(time.var)){
        index.clusterTime <- tapply(data[[time.var]],data[[cluster.var]],function(iT){as.numeric(factor(iT,levels = U.time))})
    }else if(!is.na(time.var)){
        index.clusterTime <- setNames(as.list(as.numeric(factor(data[[time.var]], levels = U.time))),U.cluster)
    }else if(!is.na(cluster.var)){
        index.clusterTime <- tapply(data[[cluster.var]],data[[cluster.var]],function(iT){1:length(iT)})
    }
    order.clusterTime <- lapply(index.clusterTime, order)
    
    ## *** export
    out <- list(data = data,
                time.var = time.var,
                U.time = U.time,
                cluster.var = cluster.var,
                U.cluster = U.cluster,
                strata.var = strata.var,
                U.strata = U.strata,
                index.cluster = index.cluster,
                index.clusterTime = index.clusterTime,
                order.clusterTime = order.clusterTime
                )

    return(out)
}

## ** .patternStructure
.patternStructure <- function(X.var, X.cor, lpdiff.rho, data,
                              time.var, index.clusterTime, order.clusterTime, U.time,
                              index.cluster, U.cluster,
                              strata.var, strata.param, U.strata){

    ## re-order by strata
    strata.param <- sort(strata.param)
    name.param <- names(strata.param)
    time.param <- setNames(vector(mode="list", length = length(strata.param)), name.param)
    if(is.na(time.var)){
        stop("Unknown time variable. \n")
    }
    ## *** find all variance patterns
    level.var <- as.numeric(droplevels(interaction(as.data.frame(X.var))))
    ls.patternVar.cluster <- lapply(U.cluster, function(iC){ ## iC <- U.cluster[[1]]
        iIndex.obs <- index.cluster[[iC]][order.clusterTime[[iC]]]
        iLevel <- paste(level.var[iIndex.obs], collapse=".")
        attr(iLevel,"index.alltime") <- index.clusterTime[[iC]][order.clusterTime[[iC]]]
        if(!is.na(strata.var)){
            attr(iLevel,"index.strata") <- which(unique(data[index.cluster[[iC]],strata.var])==U.strata)
        }else{
            attr(iLevel,"index.strata") <- 1
        }
        return(iLevel)
    })
    vec.patternVar.cluster <- unlist(ls.patternVar.cluster)
    Uindex.patternVar.cluster <- which(!duplicated(vec.patternVar.cluster))[order(vec.patternVar.cluster[!duplicated(vec.patternVar.cluster)])]
    UpatternVar.cluster <- vec.patternVar.cluster[Uindex.patternVar.cluster]
    patternVar.cluster <- as.numeric(factor(vec.patternVar.cluster, levels = UpatternVar.cluster))

    X.var2 <- setNames(lapply(Uindex.patternVar.cluster,function(iIndex){
        X.var[index.cluster[[iIndex]],,drop=FALSE]
    }),UpatternVar.cluster)
    
    for(iVar in colnames(X.var)){
        time.param[[iVar]] <- sort(match(unique(as.character(data[X.var[,iVar]==1,time.var])),U.time))[1]
    }

    ## *** find all correlation patterns
    if(is.null(X.cor)){
        UpatternCor.cluster <- "1"
        patternCor.cluster <- setNames(rep(1, length = length(U.cluster)), U.cluster)
        X.cor2 <- NULL
    }else{
        lp.cor <- as.character(interaction(as.data.frame(X.cor),drop=TRUE))
        patternCor.cluster <- sapply(U.cluster, function(iC){
            paste(lp.cor[index.cluster[[iC]][order.clusterTime[[iC]]]],collapse="|")
        })
        indexCluster.cor <- which(!duplicated(patternCor.cluster))
        UpatternCor.cluster <- sort(unique(patternCor.cluster))
        patternCor.cluster <- as.numeric(factor(patternCor.cluster, levels = UpatternCor.cluster))

        X.cor2 <- lapply(U.cluster[indexCluster.cor], function(iC){ ## iC <- 1
            ## extract design
            iX <- X.cor[index.cluster[[iC]][order.clusterTime[[iC]]],,drop=FALSE]
            ## generate all pairs for the matrix (including symmetric terms)
            iPair.time <- cbind(.unorderedPairs(1:NROW(iX), distinct = TRUE),.unorderedPairs(NROW(iX):1, distinct = TRUE))
            ## compute difference in covariates to identify coefficient
            iDiff <- as.character(interaction(as.data.frame(
                do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){ ## iCol <- 4
                    iX1 <- iX[min(iPair.time[,iCol]),,drop=FALSE]
                    iX2 <- iX[max(iPair.time[,iCol]),,drop=FALSE]
                    if(all(iX1==iX2)){return(iX1)}else{return(iX2-iX1)}                        
                }))
            ), drop = TRUE))
            ## generate design matrix relative to the coefficient
            if(length(lpdiff.rho)==1 && all(iDiff==lpdiff.rho)){
                iZ <- matrix(1, nrow = length(iDiff), ncol = 1, dimnames = list(NULL,names(lpdiff.rho)))
            }else{
                iZ <- model.matrix(~0+X,data.frame(X = factor(names(lpdiff.rho)[match(iDiff,lpdiff.rho)], levels = names(lpdiff.rho))))
                colnames(iZ) <- names(lpdiff.rho)
            }
            rownames(iZ) <- paste("(",iPair.time[1,],",",iPair.time[2,],")", sep = "")

            ## add index to build the matrix format and time
            attr(iZ,"index.vec2matrix") <- iPair.time[1,] + NROW(iX)*(iPair.time[2,]-1)
            attr(iZ,"index.pairtime") <- iPair.time
            attr(iZ,"index.Utime") <- index.clusterTime[[iC]][order.clusterTime[[iC]]]
            if(!is.na(strata.var)){
                attr(iZ,"index.strata") <- which(unique(data[index.cluster[[iC]],strata.var])==U.strata)
            }else{
                attr(iZ,"index.strata") <- 1
            }
            ## iTest <- matrix(as.character(NA),4,4); iTest[attr(iZ,"index")] <- names(lpdiff.rho)[iZ %*% 1:NCOL(iZ)]
            return(iZ)            
        })
        names(X.cor2) <- UpatternCor.cluster

        for(iPattern in 1:length(X.cor2)){ ## iPattern <- 1
            for(iParam in colnames(X.cor2[[iPattern]])){ ## iParam <- "rho"
                ## find all pair of times corresponding to the parameter
                iPairs <- attr(X.cor2[[iPattern]],"index.pairtime")[,which(as.logical(X.cor2[[iPattern]][,iParam])),drop=FALSE]
                if(length(iPairs)>0){ ## do not continue if the parameter belongs to another strata
                    ## find the pairs with the earlier time corresponding to the parameter
                    iPairs.min <- apply(iPairs,2,min)
                    iPair.min <- sort(iPairs[,intersect(which(iPairs.min==min(iPairs.min)),which.min(colSums(iPairs)))])
                    if(length(time.param[[iParam]])==0){
                        time.param[[iParam]] <- c(iPair.min)
                    }else if(min(time.param[[iParam]])>iPair.min[1] || min(time.param[[iParam]])==iPair.min[1] && sum(time.param[[iParam]])>sum(iPair.min)){
                        time.param[[iParam]] <- c(iPair.min)
                    }
                }
            }
        }
    }
    
    ## ***  assemble variance and correlation patterns to identify unique patterns
    pattern.cluster <- paste(patternVar.cluster,patternCor.cluster,sep=".")
    test.duplicated <- duplicated(pattern.cluster)
    Upattern <- data.frame(name = pattern.cluster[test.duplicated==FALSE],
                           var = patternVar.cluster[test.duplicated==FALSE],
                           cor = patternCor.cluster[test.duplicated==FALSE],
                           time = NA,
                           strata = NA,
                           cluster = NA,
                           param = NA)
    attr(Upattern,"name.UpatternVar") <- UpatternVar.cluster
    attr(Upattern,"name.UpatternCor") <- UpatternCor.cluster
    pattern.cluster <- as.numeric(factor(pattern.cluster, levels = Upattern$name))

    ## *** strata and time associated to each unique pattern
    Upattern$time <- lapply(ls.patternVar.cluster[which(!duplicated(pattern.cluster))],attr,"index.alltime")
    Upattern$strata <- lapply(ls.patternVar.cluster[which(!duplicated(pattern.cluster))],attr,"index.strata")
    Upattern$cluster <- lapply(1:length(Upattern$name),function(iN){which(iN==pattern.cluster)})
    Upattern$param <- lapply(1:length(Upattern$name),function(iN){
        iParam <- names(which(colSums(X.var2[[iN]])>0))
        if(!is.null(X.cor2)){
            iParam <- c(iParam,names(which(colSums(X.cor2[[iN]])>0)))
        }
        return(iParam)
    })

    names(Upattern$time) <- Upattern$name
    names(Upattern$strata) <- Upattern$name
    names(Upattern$cluster) <- Upattern$name
    names(Upattern$param) <- Upattern$name

    ## *** export
    out <- list(pattern.cluster = pattern.cluster,
                Upattern = Upattern,
                X.cor = X.cor2,
                X.var = X.var2,
                time.param = time.param)
    return(out)
}


##----------------------------------------------------------------------
### skeletonStructure.R ends here
