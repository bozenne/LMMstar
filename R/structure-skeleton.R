### structure-skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: sep 22 2021 (17:29) 
##           By: Brice Ozenne
##     Update #: 929
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
##'
##' @param structure [structure]
##' @param data [data.frame] dataset
##'
##' @keywords internal
##' 
##' @examples
##' \dontrun{
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
##' }
##' @export
`.skeleton` <-
    function(structure, data) UseMethod(".skeleton")

## * skeleton.IND
##' @export
.skeleton.IND <- function(structure, data){

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
    structure$U.strata <- outInit$U.strata
    index.clusterTime <- outInit$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    order.clusterTime <- outInit$order.clusterTime ## list of index for re-ordering the observations within cluster 
    index.cluster <- outInit$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    formula.var <- structure$formula$var
    X.var <- .colnameOrder(model.matrix_regularize(formula.var, data = data, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)

    if(!is.null(structure$param)){
        outPattern <- .findUpatterns(X.var = X.var, X.cor = X.cor, data = data, Upattern = structure$X$Upattern,
                                     time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                     index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                     strata.var = strata.var, strata.param = stats::setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)
        return(outPattern)
    }

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
        attr(X.var,"order")[attr(X.var,"assign")==1] <- 0
    }
        
    ## *** k
    index.k <- which(attr(X.var,"assign")>(n.strata>1))
    if(length(index.k)>0){
        M.level <- attr(X.var,"M.level")
        
        oldnames <- as.character(interaction(as.data.frame(lapply(1:NCOL(M.level), function(iVar){paste0(names(M.level)[iVar],M.level[index.k,iVar])})),sep=":", drop = TRUE) )
        newnames <- as.character(interaction(attr(X.var,"M.level")[index.k,,drop=FALSE],sep=":", drop = TRUE))

        if(!identical(colnames(X.var)[index.k],oldnames)){
            stop("Could not find the k parameters in the design matrix for the variance.\n")
        }
        colnames(X.var)[index.k] <- paste("k",newnames,sep=".")
        param.k <- colnames(X.var)[index.k]
        if(n.strata == 1){
            strata.k <- stats::setNames(rep(1,length(param.k)),param.k)
        }else{
            strata.k <- stats::setNames(match(attr(X.var,"M.level")[index.k,strata.var],structure$U.strata),param.k)

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
            attr(X.var,"reorder") <- order(strata.sigmak)
            attr(X.var,"reference.level") <- save.attr$reference.level
            attr(X.var,"M.level") <- save.attr$M.level
            param.k <- param.k[order(strata.k)]
            strata.k <- strata.k[param.k]
        }
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
    outPattern <- .findUpatterns(X.var = X.var, X.cor = X.cor, lpdiff.rho = lpdiff.rho, data = data,
                                  time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                  index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                  strata.var = strata.var, strata.param = stats::setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)
    structure$param$time <- outPattern$time.param
    structure$X <- list(Upattern = outPattern$Upattern,
                        pattern.cluster = outPattern$pattern.cluster,
                        cluster.pattern = outPattern$cluster.pattern,
                        var = outPattern$X.var,
                        cor = outPattern$X.cor)
    structure$pair.varcoef <- outPattern$pair.varcoef

    ## ** export
    return(structure)
}


## * skeleton.CS
##' @export
.skeleton.CS <- function(structure, data){

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
    
    X.var <- .colnameOrder(model.matrix_regularize(formula.var, data = data, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)
    if(n.time==1 || all(sapply(order.clusterTime,length)==1)){ ## only one timepoint per cluster so remove the correlation formula
        X.cor <- NULL
    }else{
        X.cor <- .colnameOrder(model.matrix_regularize(formula.cor, data = data, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)
    }    

    if(!is.null(structure$param)){
        outPattern <- .findUpatterns(X.var = X.var, X.cor = X.cor, lpdiff.rho = attr(structure$param,"lpdiff.rho"), data = data, Upattern = structure$X$Upattern,
                                     time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                     index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                     strata.var = strata.var, strata.param = stats::setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)
        return(outPattern)
    }

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
        attr(X.var,"order") <- rep(0,NCOL(X.var))
    }

    ## *** param rho
    if(is.null(X.cor)){
        param.rho <- NULL
        strata.rho <- NULL
    }else if(n.strata==1){
        param.rho <- "rho"
        strata.rho <- stats::setNames(1,param.rho)
        colnames(X.cor) <- param.rho
        lpdiff.rho <- stats::setNames(paste0("R.",interaction(as.data.frame(unique(X.cor)),drop=FALSE)),param.rho)
    }else{
        param.rho <- paste0("rho:",structure$U.strata)
        colnames(X.cor) <- param.rho
        strata.rho <- stats::setNames(1:n.strata,param.rho)
        lpdiff.rho <- stats::setNames(paste0("R.",interaction(as.data.frame(unique(X.cor)),drop=FALSE)),param.rho)
    }
    strata.param <- c(strata.sigma,strata.rho)[c(colnames(X.var),colnames(X.cor))]

    structure$param <- data.frame(name = c(param.sigma,param.rho),
                                  strata = c(strata.sigma,strata.rho),
                                  type = c(rep("sigma",length=length(param.sigma)),rep("rho",length=length(param.rho))),
                                  time = NA)
    structure$param <- structure$param[order(structure$param$strata),,drop=FALSE]
    rownames(structure$param) <- NULL

    ## ** pattern
    outPattern <- .findUpatterns(X.var = X.var, X.cor = X.cor, lpdiff.rho = lpdiff.rho, data = data,
                                 time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                 index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                 strata.var = strata.var, strata.param = stats::setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)
    structure$param$time <- outPattern$time.param
    structure$X <- list(Upattern = outPattern$Upattern,
                        pattern.cluster = outPattern$pattern.cluster,
                        cluster.pattern = outPattern$cluster.pattern,
                        var = outPattern$X.var,
                        cor = outPattern$X.cor)
    structure$pair.varcoef <- outPattern$pair.varcoef

    ## ** export
    return(structure)
}

## * skeleton.UN
##' @export
.skeleton.UN <- function(structure, data){
    
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

    ## ** design matrix and cluster specific patterns
    X.var <- .colnameOrder(model.matrix_regularize(formula.var, data = data, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)

    if(n.time==1 || all(sapply(order.clusterTime,length)==1)){  ## only one timepoint so remove time variable from the formula (and cor)
        X.cor <- NULL
    }else{
        X.cor <- .colnameOrder(model.matrix_regularize(formula.cor, data = data, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)
    }

    if(!is.null(structure$param)){
        outPattern <- .findUpatterns(X.var = X.var[,attr(structure$X$var,"reorder"),drop=FALSE], X.cor = X.cor, lpdiff.rho = attr(structure$param,"lpdiff.rho"),
                                     data = data, Upattern = structure$X$Upattern,
                                     time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                     index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                     strata.var = strata.var, strata.param = stats::setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)
        return(outPattern)
    }

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
        M.level <- attr(X.var,"M.level")
        oldnames <- as.character(interaction(as.data.frame(lapply(1:NCOL(M.level), function(iVar){paste0(names(M.level)[iVar],M.level[index.k,iVar])})),sep=":", drop = TRUE) )
        newnames <- as.character(interaction(attr(X.var,"M.level")[index.k,,drop=FALSE],sep=":", drop = TRUE))

        if(!identical(colnames(X.var)[index.k],oldnames)){
            stop("Could not find the k parameters in the design matrix for the variance.\n")
        }
        colnames(X.var)[index.k] <- paste("k",newnames,sep=".")
        param.k <- colnames(X.var)[index.k]
        if(n.strata == 1){
            strata.k <- stats::setNames(rep(1,length(param.k)),param.k)
        }else{
            strata.k <- stats::setNames(match(attr(X.var,"M.level")[index.k,strata.var],structure$U.strata),param.k)

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
            attr(X.var,"reference.level") <- save.attr$reference.level
            attr(X.var,"reorder") <- order(strata.sigmak)
            attr(X.var,"M.level") <- save.attr$M.level
            param.k <- param.k[order(strata.k)]
            strata.k <- strata.k[param.k]
        }
    }else{
        param.k <- NULL
        strata.k <- NULL
    }

    ## *** param rho
    if(is.null(X.cor)){
        pairn.rho <- NULL
        lpdiff.rho <- NULL
        param.rho <- NULL
        strata.rho <- NULL
    }else if(n.strata==1){
        lp.cor <- as.character(interaction(as.data.frame(X.cor),drop=TRUE))
        LPcluster.cor <- lapply(structure$U.cluster, function(iC){
            lp.cor[index.cluster[[iC]][order.clusterTime[[iC]]]]
        })
        indexCluster.cor <- which(!duplicated(LPcluster.cor))

        ## for each cluster compute all pairwise difference in covariates
        all.lpdiff.rho <- unlist(lapply(structure$U.cluster[indexCluster.cor], function(iC){ ## iC <- 1
            iX <- X.cor[index.cluster[[iC]][order.clusterTime[[iC]]],,drop=FALSE]
            if(NROW(iX)==1){return(NULL)}
            iPair.time <- .unorderedPairs(1:NROW(iX), distinct = TRUE)
            iDF.diff <- as.data.frame(do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){
                iX1 <- iX[min(iPair.time[,iCol]),,drop=FALSE]
                iX2 <- iX[max(iPair.time[,iCol]),,drop=FALSE]
                if(all(iX1==iX2)){return(cbind("R",iX1))}else{return(cbind("D",iX2-iX1))}
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

        reorder.rho <- order(time.rho[1,],time.rho[2,])
        time.rho <- time.rho[,reorder.rho,drop=FALSE]
        lpdiff.rho <- lpdiff.rho[reorder.rho]
        param.rho <- param.rho[reorder.rho]
        strata.rho <- strata.rho[reorder.rho]
        time.rho <- time.rho[reorder.rho]
    }else{
        pair.time <- .unorderedPairs(1:length(structure$U.time), distinct = TRUE)
        pair.name <- apply(pair.time,2,function(iT){paste0("(",structure$U.time[iT[1]],",",structure$U.time[iT[2]],")")})
        param.rho <- unlist(lapply(structure$U.strata, function(iStrata){paste0(paste0("rho",pair.name),":",iStrata)}))
        strata.rho <- stats::setNames(unlist(lapply(1:n.strata, function(iStrata){rep(iStrata, NCOL(pair.time))})), param.rho)
        time.rho <- do.call(cbind,lapply(1:n.strata, function(iStrata){pair.time}))
        colnames(time.rho) <- param.rho

        Data.rho <- do.call(rbind,lapply(structure$U.strata, function(iS){ ## iS <- structure$U.strata[1]
            iIndex <- which(attr(X.cor,"M.level")[,strata.var]==iS)
            iIndex.order <- iIndex[match(attr(X.cor,"M.level")[iIndex,time.var], structure$U.time)]
            iData <- attr(X.cor,"M.level")[iIndex.order,,drop=FALSE]
            iData[[strata.var]] <- factor(iData[[strata.var]], levels = structure$U.strata)
            iData[[time.var]] <- factor(iData[[time.var]], levels = structure$U.time)
            return(iData)
        }))

        UX.cor <- .colnameOrder(model.matrix_regularize(structure$formula$cor, data = Data.rho, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)[,colnames(X.cor),drop=FALSE]

        lpdiff.rho <- stats::setNames(unlist(lapply(structure$U.strata, function(iS){ ## iS <- structure$U.strata[1]
            iUX.cor <- UX.cor[Data.rho[[strata.var]]==iS,,drop=FALSE]
            if(NROW(iUX.cor)==1){return(NULL)}
            iLPdiff.rho <- apply(pair.time,2,function(iCol){paste0("D.",paste(iUX.cor[max(iCol),]-iUX.cor[min(iCol),],collapse="."))})
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
    attr(structure$param,"lpdiff.rho") <- lpdiff.rho

    ## ** pattern
    outPattern <- .findUpatterns(X.var = X.var, X.cor = X.cor, lpdiff.rho = lpdiff.rho, data = data,
                                 time.var = time.var, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime, U.time = structure$U.time,
                                 index.cluster = index.cluster, U.cluster = structure$U.cluster,
                                 strata.var = strata.var, strata.param = stats::setNames(structure$param$strata,structure$param$name), U.strata = structure$U.strata)
    structure$param$time <- outPattern$time.param
    structure$X <- list(Upattern = outPattern$Upattern,
                        pattern.cluster = outPattern$pattern.cluster,
                        cluster.pattern = outPattern$cluster.pattern,
                        var = outPattern$X.var,
                        cor = outPattern$X.cor)
    structure$pair.varcoef <- outPattern$pair.varcoef

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
    all.var <- stats::na.omit(c(time.var,cluster.var,strata.var,unlist(lapply(structure$formula,all.vars))))

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
        index.cluster <- stats::setNames(as.list(1:NROW(data)),U.cluster)
    }

    ## *** find time corresponding to each cluster
    if(is.na(cluster.var) && is.na(time.var)){
        index.clusterTime <- stats::setNames(as.list(rep(1,NROW(data))),U.cluster)
    }else if(!is.na(cluster.var) && !is.na(time.var)){
        index.clusterTime <- tapply(data[[time.var]],data[[cluster.var]],function(iT){as.numeric(factor(iT,levels = U.time))})
    }else if(!is.na(time.var)){
        index.clusterTime <- stats::setNames(as.list(as.numeric(factor(data[[time.var]], levels = U.time))),U.cluster)
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

## ** .buildPatterns
## convert design matrix into character string with associated strata and time
.buildVarPatterns <- function(X.var, data, Upattern,
                              U.cluster, index.cluster, index.clusterTime, order.clusterTime,
                              U.strata, strata.var){

    level.var <- interaction(as.data.frame(X.var),drop=TRUE)
    if(!is.null(Upattern)){
        level.var <- factor(level.var, levels = attr(Upattern,"levels.UpatternVar"))
    }
    Ulevel.var <- levels(level.var)
    numlevel.var <- as.numeric(level.var)

    out <- lapply(U.cluster, function(iC){ ## iC <- U.cluster[[1]]
        iIndex.obs <- index.cluster[[iC]][order.clusterTime[[iC]]]
        iLevel <- paste(numlevel.var[iIndex.obs], collapse=".")
        attr(iLevel,"index.alltime") <- index.clusterTime[[iC]][order.clusterTime[[iC]]]
        if(!is.na(strata.var)){
            attr(iLevel,"index.strata") <- which(unique(data[index.cluster[[iC]],strata.var])==U.strata)
        }else{
            attr(iLevel,"index.strata") <- 1
        }
        return(iLevel)
    })
    attr(out,"levels") <- Ulevel.var
    attr(out,"reorder") <- attr(X.var,"reorder")
    return(out)
    
}

.buildCorPatterns <- function(X.cor, Upattern,
                              U.cluster, index.cluster, order.clusterTime){

    level.cor <- interaction(as.data.frame(X.cor),drop=TRUE)
    if(!is.null(Upattern)){
        level.cor <- factor(level.cor, levels = attr(Upattern,"levels.UpatternCor"))
    }
    Ulevel.cor <- levels(level.cor)
    numlevel.cor <- as.numeric(level.cor)
    
    out <- sapply(U.cluster, function(iC){
        paste(numlevel.cor[index.cluster[[iC]][order.clusterTime[[iC]]]],collapse="|")
    })
    attr(out,"levels") <- Ulevel.cor
    return(out)
}


## ** .findUpatterns
.findUpatterns <- function(X.var, X.cor, Upattern = NULL, lpdiff.rho, data,
                           time.var, index.clusterTime, order.clusterTime, U.time,
                           index.cluster, U.cluster,
                           strata.var, strata.param, U.strata){

    ## re-order by strata
    if(is.null(Upattern)){
        strata.param <- sort(strata.param)
        name.param <- names(strata.param)
        time.param <- stats::setNames(vector(mode="list", length = length(strata.param)), name.param)
        if(is.na(time.var)){
            stop("Unknown time variable. \n")
        }
    }

    ## *** find all variance patterns
    ## variance pattern for all clusters
    ls.patternVar.cluster <- .buildVarPatterns(X.var = X.var, data = data, Upattern = Upattern,
                                               U.cluster = U.cluster, index.cluster = index.cluster, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime,
                                               U.strata = U.strata, strata.var = strata.var)
    patternVar.cluster <- unlist(ls.patternVar.cluster) 
    ## find unique patterns and order them
    if(is.null(Upattern)){
        out.orderUPatternVar <- .orderUpattern(pattern.cluster = patternVar.cluster, time.cluster = index.clusterTime, order.cluster = order.clusterTime, U.time = U.time) 
        UpatternVar.cluster <- out.orderUPatternVar$name
    }else{
        UpatternVar.cluster <- attr(Upattern,"name.UpatternVar")
    }
    ## update vector giving the index of the pattern corresponding to each cluster
    patternVar.cluster <- as.numeric(factor(patternVar.cluster, levels = UpatternVar.cluster))
    ## generate design matrix associated to each pattern
    if(is.null(Upattern)){
        X.var2 <- stats::setNames(lapply(out.orderUPatternVar$cluster,function(iIndex){ ## iIndex <- Uindex.patternVar.cluster[1]
            iX <- X.var[index.cluster[[iIndex]],,drop=FALSE]
            iIndex.param <- which(colSums(iX)!=0)
            attr(iX,"indicator.param") <- stats::setNames(lapply(iIndex.param,function(iCol){which(tcrossprod(iX[,iCol],rep(1,NROW(iX))) + t(tcrossprod(iX[,iCol],rep(1,NROW(iX)))) > 0)}),
                                                          names(iIndex.param))
            return(iX)
        }),UpatternVar.cluster)
        attr(X.var2,"assign") <- attr(X.var,"assign")
    }
    ## first timepoint where each variance parameter appear
    if(is.null(Upattern)){
        for(iVar in colnames(X.var)){
            time.param[[iVar]] <- sort(match(unique(as.character(data[X.var[,iVar]==1,time.var])),U.time))[1]
        }
    }
    ## *** find all correlation patterns
    if(is.null(X.cor)){
        UpatternCor.cluster <- "1"
        patternCor.cluster <- stats::setNames(rep(1, length = length(U.cluster)), U.cluster)
        X.cor2 <- NULL
    }else{
        ## correlation pattern for all clusters
        ls.patternCor.cluster <- .buildCorPatterns(X.cor = X.cor, Upattern = Upattern,
                                                   U.cluster = U.cluster, index.cluster = index.cluster, order.clusterTime = order.clusterTime) 
        ## find unique patterns and order them
        if(is.null(Upattern)){
            out.orderUPatternCor <- .orderUpattern(pattern.cluster = ls.patternCor.cluster, time.cluster = index.clusterTime, order.cluster = order.clusterTime, U.time = U.time) 
            UpatternCor.cluster <- out.orderUPatternCor$name
        }else{
            UpatternCor.cluster <- attr(Upattern,"name.UpatternCor")
        }
        ## update vector giving the index of the pattern corresponding to each cluster
        patternCor.cluster <- as.numeric(factor(ls.patternCor.cluster, levels = UpatternCor.cluster))
        ## generate design matrix associated to each pattern
        if(is.null(Upattern)){
            X.cor2 <- lapply(out.orderUPatternCor$cluster, function(iC){ ## iC <- 1
                ## extract design
                iX <- X.cor[index.cluster[[iC]][order.clusterTime[[iC]]],,drop=FALSE]
                iTime <- index.clusterTime[[iC]][order.clusterTime[[iC]]]
                if(NROW(iX)==1){return(NULL)}
                ## generate all pairs for the matrix (including symmetric terms)
                iPair.time <- cbind(.unorderedPairs(iTime, distinct = TRUE),.unorderedPairs(rev(iTime), distinct = TRUE))
                ## compute difference in covariates to identify coefficient
                iDiff <- as.character(interaction(as.data.frame(
                    do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){ ## iCol <- 2
                        iPair <- iPair.time[,iCol]
                        iX1 <- iX[iTime == min(iPair),,drop=FALSE]
                        iX2 <- iX[iTime == max(iPair),,drop=FALSE]
                        if(all(iX1==iX2)){return(cbind(ref="R",iX1))}else{return(cbind(ref = "D",iX2-iX1))}                        
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
            iPair.index <- matrix(NA, nrow = 2, ncol = NCOL(iPair.time))
            iPair.index[] <- as.numeric(factor(iPair.time, levels = iTime, labels = 1:length(iTime)))
            iIndex.param <- which(colSums(iZ)!=0)
            attr(iZ,"index.vec2matrix") <- iPair.index[1,] + NROW(iX)*(iPair.index[2,]-1)
            attr(iZ,"indicator.param") <- stats::setNames(lapply(iIndex.param,function(iI){attr(iZ,"index.vec2matrix")[which(iZ[,iI]>0)]}),names(iIndex.param))
            attr(iZ,"index.pairtime") <- iPair.time
            attr(iZ,"index.Utime") <- iTime
            if(!is.na(strata.var)){
                attr(iZ,"index.strata") <- which(unique(data[index.cluster[[iC]],strata.var])==U.strata)
            }else{
                attr(iZ,"index.strata") <- 1
            }
            ## iTest <- matrix(as.character(NA),4,4); iTest[attr(iZ,"index")] <- names(lpdiff.rho)[iZ %*% 1:NCOL(iZ)]
            return(iZ)            
        })
        names(X.cor2) <- UpatternCor.cluster
        }
        ## first pair of timepoints where each correlation parameter appear
        if(is.null(Upattern)){
            for(iPattern in 1:length(X.cor2)){ ## iPattern <- 7
                if(is.null(X.cor2[[iPattern]])){next}
                for(iParam in colnames(X.cor2[[iPattern]])){ ## iParam <- "rho(0,52)"
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
    }
    
    ## ***  assemble variance and correlation patterns to identify unique patterns
    pattern.cluster <- paste(patternVar.cluster,patternCor.cluster,sep=".")
    if(is.null(Upattern)){
        index.Nduplicated <- which(duplicated(pattern.cluster)==FALSE)
        index.Nduplicated2 <- index.Nduplicated[order(pattern.cluster[index.Nduplicated])]
        Upattern <- data.frame(name = pattern.cluster[index.Nduplicated2],
                               var = patternVar.cluster[index.Nduplicated2],
                               cor = patternCor.cluster[index.Nduplicated2],
                               time = NA,
                               strata = NA,
                               ncluster = NA,
                               param = NA)
        attr(Upattern,"name.UpatternVar") <- UpatternVar.cluster
        attr(Upattern,"name.UpatternCor") <- UpatternCor.cluster
        attr(Upattern,"levels.UpatternVar") <- attr(ls.patternVar.cluster,"levels")
        attr(X.var2,"reorder") <- attr(ls.patternVar.cluster,"reorder")
        attr(Upattern,"levels.UpatternCor") <- attr(ls.patternCor.cluster,"levels")
    }else{
        pattern.cluster <- as.character(factor(pattern.cluster, levels = Upattern$name))
        cluster.pattern <- stats::setNames(lapply(Upattern$name,function(iN){which(iN==pattern.cluster)}),Upattern$name)
        out <- list(Upattern = Upattern,
                    pattern.cluster = stats::setNames(pattern.cluster,U.cluster),
                    cluster.pattern = cluster.pattern,
                    cor = NULL,
                    var = NULL)
        return(out)
    }
    
    ## *** strata and time associated to each unique pattern
    Upattern$time <- lapply(ls.patternVar.cluster[index.Nduplicated2],attr,"index.alltime")
    Upattern$strata <- sapply(ls.patternVar.cluster[index.Nduplicated2],attr,"index.strata")
    Upattern$ncluster <- sapply(1:length(Upattern$name),function(iN){sum(iN==pattern.cluster)})
    cluster.pattern <- lapply(Upattern$name,function(iN){which(iN==pattern.cluster)})
    Upattern$param <- lapply(1:length(Upattern$name),function(iN){ ## iN <- 1
        iParam <- names(which(colSums(X.var2[[Upattern[iN,"var"]]])>0))
        if(!is.null(X.cor2) && !is.null(X.cor2[[Upattern[iN,"cor"]]])){
            iParam <- c(iParam,names(which(colSums(X.cor2[[Upattern[iN,"cor"]]])>0)))
        }
        return(iParam)
    })

    names(Upattern$time) <- Upattern$name
    names(cluster.pattern) <- Upattern$name
    names(Upattern$param) <- Upattern$name

    for(iPatternVar in 1:length(X.var2)){## iPatternVar <- 1
        attr(X.var2[[iPatternVar]],"index.cluster") <- unname(sort(unlist(cluster.pattern[which(Upattern$var==iPatternVar)])))        
        attr(X.var2[[iPatternVar]],"index.obs") <- stats::setNames(lapply(attr(X.var2[[iPatternVar]],"index.cluster"), function(iC){
            index.cluster[[iC]][order.clusterTime[[iC]]]
        }),U.cluster[attr(X.var2[[iPatternVar]],"index.cluster")])
    }
    if(!is.null(X.cor2)){
        for(iPatternCor in 1:length(X.cor2)){## iPatternVar <- 1
            if(!is.null(X.cor2[[iPatternCor]])){
                attr(X.cor2[[iPatternCor]],"index.cluster") <- unname(sort(unlist(cluster.pattern[which(Upattern$cor==iPatternCor)])))        
                attr(X.cor2[[iPatternCor]],"index.obs") <- stats::setNames(lapply(attr(X.cor2[[iPatternCor]],"index.cluster"), function(iC){
                    index.cluster[[iC]][order.clusterTime[[iC]]]
                }),U.cluster[attr(X.cor2[[iPatternCor]],"index.cluster")])
            }
        }
        attr(X.cor2,"assign") <- rep(1,length(unique(unlist(lapply(X.cor2,colnames)))))
    }

    ## *** pair of variance-covariance coefficients
    pair.varcoef <- stats::setNames(lapply(Upattern$name, function(iPattern){## ## iPattern <- structure$X$Upattern$name[1]
        iParam <- Upattern$param[[iPattern]]

        iOut <- .unorderedPairs(iParam)
        attr(iOut, "key") <- matrix(NA, nrow = length(iParam), ncol = length(iParam), dimnames = list(iParam,iParam))
        for(iCol in 1:NCOL(iOut)){
            attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
            attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol]] <- iCol
        }
        return(iOut)
    }),Upattern$name)

    ## *** export
    out <- list(pattern.cluster = stats::setNames(pattern.cluster,U.cluster),
                cluster.pattern = cluster.pattern,
                Upattern = Upattern,
                X.cor = X.cor2,
                X.var = X.var2,
                time.param = time.param,
                pair.varcoef = pair.varcoef)
    return(out)
}


## ** colnameOrder
## reorder the variable in the column name
.colnameOrder <- function(X, strata.var, n.strata){

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

## ** .orderUpattern
## .orderMatTime(time.cluster = list(1:5,1,1:3,2:4), order.cluster = NULL, U.time = 1:5)
.orderUpattern <- function(pattern.cluster, time.cluster, order.cluster, U.time){

    indexU.cluster <- which(!duplicated(pattern.cluster))
    score.pattern <- pattern.cluster[indexU.cluster]

    ## score the time vector corresponding to each cluster
    n.time <- length(U.time)
    n.Ucluster <- length(indexU.cluster)
    score.time <- do.call(rbind,lapply(indexU.cluster, function(iC){ ## iC <- 1
        if(!is.null(order.cluster)){
            iTime <- time.cluster[[iC]][order.cluster[[iC]]]
        }else{
            iTime <- time.cluster[[iC]]
        }
        if(n.time == length(iTime)){
            c(-1,0) ## no missing, no missing
        }else{
            c("firstmis" = n.time-min(setdiff(1:n.time, iTime)),"nmis"=n.time-length(iTime)) ## n-position of the first missing value, number of missing values
        }
    }))

    out <- list(cluster = indexU.cluster[order(score.time[,1],score.time[,2],score.pattern)])
    out$name <- pattern.cluster[out$cluster]
    return(out)
}



##----------------------------------------------------------------------
### structure-skeleton.R ends here
