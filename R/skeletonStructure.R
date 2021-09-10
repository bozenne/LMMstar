### skeletonStructure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: sep  9 2021 (18:18) 
##           By: Brice Ozenne
##     Update #: 294
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
##'
##' ## independence
##' skeletonStructure(IND(~1), data = gastricbypassL)
##' skeletonStructure(IND(~time), data = gastricbypassL)
##' skeletonStructure(IND(~time|id), data = gastricbypassL)
##' skeletonStructure(IND(~time:gender|id), data = gastricbypassL)
##' 
##' ## compound symmetry
##' skeletonStructure(CS(~1|id), data = gastricbypassL)
##' skeletonStructure(CS(~time|id), data = gastricbypassL)
##' skeletonStructure(CS(gender~time|id), data = gastricbypassL)
##' 
##' ## unstructured
##' skeletonStructure(UN(~time|id), data = gastricbypassL)
##' skeletonStructure(UN(~time*gender|id), data = gastricbypassL)
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
    time.var <- outInit$time.var
    U.time <- outInit$U.time
    n.time <- length(U.time)
    cluster.var <- outInit$cluster.var
    U.cluster <- outInit$U.cluster
    strata.var <- outInit$strata.var
    U.strata <- outInit$U.strata
    index.clusterTime <- outInit$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    order.clusterTime <- outInit$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- outInit$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    for(iVar in all.vars(structure$formula.var)){
        data[[iVar]] <- as.factor(data[[iVar]])
    }
    X.var <- .model.matrix_noVarName(structure$formula.var, data = data, lastVar = strata.var)

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

    ## ** patterns
    structure$pattern <- .patternStructure(X.var = X.var, X.cor = NULL, data = data,
                                           index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = U.time,
                                           index.cluster = index.cluster, U.cluster = U.cluster,
                                           strata.var = strata.var, strata.param = strata.param, U.strata = U.strata)
    ## ** export
    structure$param.name <- c(param.sigma,param.k)
    structure$param.strata <- c(strata.sigma,strata.k)
    structure$param.time <- structure$pattern$param.time
    structure$pattern$param.time <- NULL
    structure$X.var <-  X.var
    structure$X.cor <-  NULL
    structure$U.time <- U.time
    return(structure)
}


## * skeletonStructure.CS
##' @export
skeletonStructure.CS <- function(structure, data){

    ## ** prepare
    outInit <- .initSkeleton(data = data, structure = structure)
    time.var <- outInit$time.var
    U.time <- outInit$U.time
    n.time <- length(U.time)
    cluster.var <- outInit$cluster.var
    U.cluster <- outInit$U.cluster
    strata.var <- outInit$strata.var
    U.strata <- outInit$U.strata
    n.strata <- length(U.strata)
    index.clusterTime <- outInit$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    order.clusterTime <- outInit$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- outInit$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    X.var <- .model.matrix_noVarName(structure$formula.var, data = data, lastVar = strata.var)
    X.cor <- .model.matrix_noVarName(structure$formula.cor, data = data, lastVar = strata.var)
    
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
        param.sigma <- paste0("sigma:",U.strata)
        if(!identical(colnames(X.var)[attr(X.var,"assign")==1],U.strata)){
            stop("Could not find the strata-specific intercepts in the design matrix for the variance.\n")
        }
        colnames(X.var)[attr(X.var,"assign")==1] <- param.sigma
        strata.sigma <- stats::setNames(1:n.strata,param.sigma)
    }

    ## *** k
    if(NCOL(X.var)>1){
        ## find name of the coeffient and update design matrix
        index.k <- 2:NCOL(X.var)
        time.k <- rep(as.character(NA), length(index.k))
        for(iTime in 1:n.time){ ## iTime <- 1
            iPattern <- paste0("^",time.var,U.time[iTime],"$")
            iTest <- grepl(pattern = iPattern, colnames(X.var)[index.k])
            if(any(iTest)){
                time.k[iTest] <- U.time[iTime]
            }
        }
        colnames(X.var)[index.k] <- paste0("k.",gsub(paste0("^",time.var),"",colnames(X.var)[index.k]))
        param.k <- colnames(X.var)[index.k]

        ## find strata of the coeffient
        strata.k <- stats::setNames(rep(1,length(param.k)),param.k)
    }else{
        param.k <- NULL
        strata.k <- NULL
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
    }else{
        param.rho <- paste0("rho:",U.strata)
        colnames(X.cor) <- param.rho
        strata.rho <- stats::setNames(1:n.strata,param.rho)
    }
    strata.param <- c(strata.sigma,strata.rho)[c(colnames(X.var),colnames(X.cor))]

    ## ** patterns
    structure$pattern <- .patternStructure(X.var = X.var, X.cor = X.cor, data = data,
                                           index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = U.time,
                                           index.cluster = index.cluster, U.cluster = U.cluster,
                                           strata.var = strata.var, strata.param = strata.param, U.strata = U.strata)
    
    ## ** export
    structure$param.name <- c(param.sigma,param.rho)
    structure$param.strata <- c(strata.sigma,strata.rho)
    structure$param.time <- structure$pattern$param.time
    structure$pattern$param.time <- NULL
    structure$X.var <-  X.var
    structure$X.cor <-  NULL
    structure$U.time <- U.time
    return(structure)
}

## * skeletonStructure.UN
##' @export
skeletonStructure.UN <- function(structure, data){

    ## ** prepare
    outInit <- .initSkeleton(data = data, structure = structure)
    time.var <- outInit$time.var
    U.time <- outInit$U.time
    n.time <- length(U.time)
    cluster.var <- outInit$cluster.var
    U.cluster <- outInit$U.cluster
    strata.var <- outInit$strata.var
    U.strata <- outInit$U.strata
    n.strata <- length(U.strata)
    index.clusterTime <- outInit$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    order.clusterTime <- outInit$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- outInit$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    X.var <- .model.matrix_noVarName(structure$formula.var, data = data, lastVar = strata.var)
    X.cor <- .model.matrix_noVarName(structure$formula.cor, data = data, lastVar = strata.var)

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
        param.sigma <- paste0("sigma:",U.strata)
        if(!identical(colnames(X.var)[attr(X.var,"assign")==1],U.strata)){
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
            strata.k <- stats::setNames(match(attr(X.var,"M.level")[index.k,strata.var],U.strata),param.k)

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
        LPcluster.cor <- lapply(U.cluster, function(iC){
            lp.cor[index.cluster[[iC]][order.clusterTime[[iC]]]]
        })
        indexCluster.cor <- which(!duplicated(LPcluster.cor))

        ## for each cluster compute all pairwise difference in covariates
        all.lpdiff.rho <- unlist(lapply(U.cluster[indexCluster.cor], function(iC){ ## iC <- 1
            iX <- X.cor[index.cluster[[iC]][order.clusterTime[[iC]]],,drop=FALSE]
            iPair.time <- .unorderedPairs(1:NROW(iX), distinct = TRUE)
            iDF.diff <- as.data.frame(do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){iX[iPair.time[2,iCol],,drop=FALSE]-iX[iPair.time[1,iCol],,drop=FALSE]})))
            iVec.diff <- as.character(interaction(iDF.diff, drop=TRUE))

            iCov <- as.character(interaction(data[index.cluster[[iC]][order.clusterTime[[iC]]],names(attr(X.cor,"M.level")),drop=FALSE],drop=TRUE))
            names(iVec.diff) <- sapply(1:NCOL(iPair.time),function(iCol){paste0("(",iCov[iPair.time[2,iCol]],",",iCov[iPair.time[1,iCol]],")")})
            return(iVec.diff)
        }))
        test.duplicated <- duplicated(all.lpdiff.rho)
        lpdiff.rho <- all.lpdiff.rho[test.duplicated==FALSE]

        param.rho <- paste0("rho",names(all.lpdiff.rho)[test.duplicated==FALSE])
        strata.rho <- stats::setNames(rep(1,length(param.rho)),param.rho)
        time.rho <- do.call(cbind,lapply(U.cluster[indexCluster.cor], function(iC){
            .unorderedPairs(index.clusterTime[[iC]], distinct = TRUE)
        }))[,test.duplicated==FALSE]
        
    }else{
        pair.time <- .unorderedPairs(1:length(U.time), distinct = TRUE)
        pair.name <- apply(pair.time,2,function(iT){paste0("(",U.time[iT[1]],",",U.time[iT[2]],")")})
        param.rho <- unlist(lapply(U.strata, function(iStrata){paste0(paste0("rho",pair.name),":",iStrata)}))
        strata.rho <- setNames(unlist(lapply(1:n.strata, function(iStrata){rep(iStrata, NCOL(pair.time))})), param.rho)
        time.rho <- do.call(cbind,lapply(1:n.strata, function(iStrata){pair.time}))
        colnames(time.rho) <- param.rho

        lpdiff.rho <- stats::setNames(unlist(lapply(U.strata, function(iS){ ## iS <- U.strata[1]
            iIndex <- which(attr(X.cor,"M.level")[,strata.var]==iS)
            iIndex.order <- iIndex[match(attr(X.cor,"M.level")[iIndex,time.var], U.time)]
            iData <- attr(X.cor,"M.level")[iIndex.order,,drop=FALSE]
            iData[[strata.var]] <- factor(iData[[strata.var]], levels = U.strata)
            iData[[time.var]] <- factor(iData[[time.var]], levels = U.time)
            iUX.cor <- .model.matrix_noVarName(structure$formula.cor, data = iData, lastVar = strata.var)
            iLPdiff.rho <- apply(pair.time,2,function(iCol){paste(iUX.cor[iCol[2],iIndex]-iUX.cor[iCol[1],iIndex],collapse=".")})
            return(iLPdiff.rho)
        })),param.rho)
    }

    strata.param <- c(c(strata.sigma,strata.k)[colnames(X.var)],strata.rho)
    attr(strata.param,"lpdiff.rho") <- lpdiff.rho

    ## ** patterns
    structure$pattern <- .patternStructure(X.var = X.var, X.cor = X.cor, data = data,
                                           index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = U.time,
                                           index.cluster = index.cluster, U.cluster = U.cluster,
                                           strata.var = strata.var, strata.param = strata.param, U.strata = U.strata)
    browser()
    

    
    ## ** export
    structure$param.name <- c(param.sigma,param.rho)
    structure$param.strata <- c(strata.sigma,strata.rho)
    structure$param.time <- structure$pattern$param.time
    structure$pattern$param.time <- NULL
    structure$X.var <-  X.var
    structure$X.cor <-  NULL
    structure$U.time <- U.time
    return(structure)
}

## * helpers
## ** .initSkeleton
.initSkeleton <- function(data, structure){
    
    ## *** find variable names
    time.var <- structure$time.var
    cluster.var <- structure$cluster
    strata.var <- structure$strata

    ## *** find unique levels of each variable
    if(!is.null(cluster.var)){
        U.cluster <- as.character(unique(data[[cluster.var]]))
    }else{
        U.cluster <- paste0("id",1:NROW(data))
    }
    
    if(!is.null(time.var)){
        U.time <- as.character(unique(data[[time.var]]))
    }else if(!is.null(cluster.var)){
        U.time <- paste0("T",1:max(table(data[[cluster.var]])))
    }else{
        U.time <- "T1"
    }

    if(!is.null(strata.var)){
        U.strata <- as.character(unique(data[[strata.var]]))
    }else{
        U.strata <- "S1"
    }

    ## *** find position of each cluster
    if(!is.null(cluster.var)){
        index.cluster <- tapply(1:NROW(data),data[[cluster.var]],function(iI){iI})
    }else{
        index.cluster <- setNames(as.list(1:NROW(data)),U.cluster)
    }

    ## *** find time corresponding to each cluster
    if(is.null(cluster.var) && is.null(time.var)){
        index.clusterTime <- setNames(as.list(rep(1,NROW(data))),U.cluster)
    }else if(!is.null(cluster.var) && !is.null(time.var)){
        index.clusterTime <- tapply(data[[time.var]],data[[cluster.var]],function(iT){as.numeric(factor(iT,levels = U.time))})
    }else if(!is.null(time.var)){
        index.clusterTime <- setNames(as.list(as.numeric(factor(data[[time.var]], levels = U.time))),U.cluster)
    }else if(!is.null(cluster.var)){
        index.clusterTime <- tapply(data[[cluster.var]],data[[cluster.var]],function(iT){1:length(iT)})
    }
    order.clusterTime <- lapply(index.clusterTime, order)
    
    ## *** export
    out <- list(time.var = time.var,
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
.patternStructure <- function(X.var, X.cor, data,
                              index.clusterTime, order.clusterTime, U.time,
                              index.cluster, U.cluster,
                              strata.var, strata.param, U.strata){

    ## for each individual
    ## should be split: first variance patterns
    ##            then: correlation pattern
    ## then assemble
    X.vcov <- cbind(X.var,X.cor)
    level.vcov <- as.numeric(droplevels(interaction(as.data.frame(X.vcov))))

    ## *** pattern of each cluster 
    pattern.cluster <- sapply(U.cluster, function(iC){ ## iC <- U.cluster[[1]]
        iIndex.obs <- index.cluster[[iC]][order.clusterTime[[iC]]]
        return(paste(level.vcov[iIndex.obs], collapse="."))
    })

    ## *** unique set of patterns
    indexEx.Upattern <- which(duplicated(pattern.cluster)==FALSE)
    Upattern <- pattern.cluster[indexEx.Upattern]

    ## *** representative element of each pattern (first cluster corresponding to each unique pattern)
    clusterEx.Upattern <- setNames(names(pattern.cluster[indexEx.Upattern]),Upattern)
    obsEx.Upattern <- lapply(clusterEx.Upattern, function(iC){index.cluster[[iC]][order.clusterTime[[iC]]]})

    ## *** strata associated to each unique pattern
    if(!is.null(strata.var)){
        strata.Upattern <- stats::setNames(match(sapply(obsEx.Upattern, function(iObs){unique(data[iObs,strata.var])}), U.strata), Upattern) 
    }else{
        strata.Upattern <- stats::setNames(rep(1,length(Upattern)), Upattern)
    }

    ## *** timepoints associated to each pattern
    indexTime.Upattern <- stats::setNames(lapply(clusterEx.Upattern, function(iC){index.clusterTime[[iC]][order.clusterTime[[iC]]]}), Upattern)
    ## time.Upattern <- stats::setNames(lapply(clusterEx.Upattern, function(iC){U.time[order.clusterTime[[iC]]]}), Upattern)
    
    ## *** design matrix associated to each pattern
    name.param <- names(strata.param)
    index.varparam <- which(name.param %in% colnames(X.var))
    X.var.Upattern <- stats::setNames(lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]
        iIndex.obs <- obsEx.Upattern[[iPattern]] ## select observations of the representative (re-ordered in increasing time)
        iIndex.param <- intersect(which(strata.param %in% strata.Upattern[iPattern]), index.varparam) ## select variance parameters from the same strata
        return(X.var[iIndex.obs,name.param[iIndex.param],drop=FALSE])
    }), Upattern)

    if(!is.null(X.cor)){
        name.corparam <- setdiff(name.param, colnames(X.var))
        browser()
        if(length(U.strata)==1){
            X.cor.Upattern <- stats::setNames(lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]
                iX.tempo <- X.cor[obsEx.Upattern[[iPattern]],,drop=FALSE] ## design matrix for time
                iPairs <- .unorderedPairs(1:NROW(iX.tempo), distinct = TRUE) ## generate all pairs of times
                iOut <- matrix(0, nrow = NCOL(iPairs), ncol = length(name.corparam),
                               dimnames = list(paste0("(",iPairs[1,],",",iPairs[2,],")"),name.corparam)) ## design matrix for distinct pairs of time
                iLPdiff <- as.character(apply(iPairs,2,function(iCol){paste(iX.tempo[iCol[2],]-iX.tempo[iCol[1],],collapse=".")}))
                attr(strata.param,"lpdiff")
                browser()
            
                
                iOut
                return(iOut)
            }), Upattern)
        }else{
            browser()
            X.cor.Upattern <- stats::setNames(lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]
                iIndex.obs <- obsEx.Upattern[[iPattern]] ## select observations of the representative (re-ordered in increasing time)
                iIndex.param <- intersect(name.param[strata.param %in% strata.Upattern[iPattern]], index.corparam) ## select correlation parameters from the same strata
                iX.tempo <- X.cor[iIndex.obs,name.param[iIndex.param],drop=FALSE] ## design matrix for time
                iPairs <- .unorderedPairs(1:NROW(iX.tempo), distinct = TRUE) ## generate all pairs of times
                iOut <- matrix(0, nrow = NCOL(iPairs), ncol = NCOL(iX.tempo),
                               dimnames = list(paste0("(",iPairs[1,],",",iPairs[2,],")"),colnames(iX.tempo))) ## design matrix for distinct pairs of time
            
            iUPairs.indexparam <- as.numeric(as.factor(apply(iPairs,2,function(iCol){paste(iX.tempo[iCol[2],]-iX.tempo[iCol[1],],collapse=".")})))
            iOut
            return(iOut)
            }), Upattern)
        }
        
        
    }else{
        X.cor.Upattern <- NULL
    }
    
    ## *** time associated to the coefficients
    name.param <- names(strata.param)
    param.time <- setNames(rep(NA,length(name.param)),name.param)
    
    ## variance
    name.varparam <- intersect(colnames(X.var),name.param)
    for(iParam in name.varparam){ ## iParam <- name.varparam[2]
        iLS.time <- lapply(Upattern, function(iPattern){indexTime.Upattern[[iPattern]][which(X.var.Upattern[[iPattern]][,iParam]>0)]})
        param.time[iParam] <- min(unlist(iLS.time))
    }


    ## correlation
    if(!is.null(X.cor)){
        X.cor.Upattern
        browser()
    }

    ## *** export
    out <- list(pattern.cluster = pattern.cluster,
                Upattern = Upattern,
                strata.Upattern = strata.Upattern,
                indexTime.Upattern = indexTime.Upattern,
                X.var.Upattern = X.var.Upattern,
                X.cor.Upattern = X.cor.Upattern,
                param.time = param.time)
    return(out)
}


##----------------------------------------------------------------------
### skeletonStructure.R ends here
