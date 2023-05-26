### structure-skeletonRho.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 11 2023 (13:27) 
## Version: 
## Last-Updated: maj 26 2023 (17:07) 
##           By: Brice Ozenne
##     Update #: 208
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .skeletonRho
##' @title Parametrization of Correlation Structure.
##' @description Parametrization of the correlation structure.
##' @noRd
##'
##' @param structure [structure]
##' @param data [data.frame] dataset.
##' @param U.cluster [character vector] cluster levels.
##' @param index.cluster [list of integer vector] position of each cluster in the dataset.
##' @param U.time [character vector] time levels.
##' @param index.clusterTime [list of integer vector] time index for each cluster.
##' @param U.strata [character vector] time levels.
##' @param index.clusterStrata [list of integer vector] strata index for each cluster.
##' @param sep [character vector of length 2] characters used to name the variance parameters,
##' the first is used to specify the covariate level whereas
##' the second is used to specify the strata level (when more than 1 strata).
##'
##' @keywords internal

`.skeletonRho` <-
    function(structure, data, 
             U.cluster, index.cluster,
             U.time, index.clusterTime, 
             U.strata, index.clusterStrata, sep) UseMethod(".skeletonRho")

## * .skeletonRho.ID
.skeletonRho.ID <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata, sep){

    ## no correlation parameter
    return(structure)

}

## * .skeletonRho.IND
.skeletonRho.IND <- .skeletonRho.ID

## * .skeletonRho.CS
.skeletonRho.CS <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata, sep = c(".",":")){

    ## ** extract information
    ## parameters
    index.sigma <- structure$param[structure$param$type=="sigma","index.level"]
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strata.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    index.k <- structure$param[structure$param$type=="k","index.level"]
    param.k <- structure$param[structure$param$type=="k","name"]
    strata.k <- structure$param[structure$param$type=="k","index.strata"]

    ## design matrix (reduce sample size to unique replicates)
    XpairPattern <- .pairPatternX(structure$X$cor, lp.obs = structure$X$lp.cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    if(length(index.k)>0){
        X.var.k <- rowSums(sweep(structure$X$var[,param.k,drop=FALSE], FUN = "*", STATS = 1:length(param.k), MARGIN = 2))
        if(any(X.var.k>length(param.k))){
            stop("Something went wrong when identifying the variance multiplier parameters. \n")
        }
        param.k.obs <- c(NA, param.k)[X.var.k+1]
    }
    ## variables
    strata.var <- structure$name$strata
    n.strata <- length(U.strata)
    xfactor <- structure$xfactor$cor

    ## levels
    M.level <- attr(structure$X$cor, "M.level")
    if(length(xfactor)>0){
        for(iVar in names(xfactor)){
            M.level[[iVar]] <- factor(M.level[[iVar]], levels = xfactor[[iVar]])
        }
    }
    regressor.var <- setdiff(colnames(M.level),strata.var) ## non-strata variable
    n.regressor <- length(regressor.var)
    level.cor <- rownames(XpairPattern$lp2X)
    
    ## others
    formula <- attr(structure$X$cor, "formula")
    X.level <- model.matrix(formula, M.level)
    heterogeneous <- structure$heterogeneous

    ## ** special case without regressors
    if(n.regressor==0){
        browser()
        if(n.strata==1){
            strata.rho <- 1
            param.rho <- "rho"
        }else{
            strata.rho <- unname(which(sapply(XpairPattern$LpU.strata,NROW)>0))
            param.rho <- paste0("rho",sep[2],U.strata[strata.rho])
        }
        structure.rho <- data.frame(name = param.rho,
                                    index.strata = strata.rho,
                                    type = rep("rho",length=length(active.strata)),
                                    index.level = as.numeric(NA),
                                    level = level.rho,
                                    code = code.rho,
                                    code.x = code.x.rho,
                                    code.y = code.y.rho,
                                    sigma = param.sigma[match(strata.rho,strata.sigma)],
                                    k.x = k.x,
                                    k.y = k.y,                                  
                                    stringsAsFactors = FALSE)
        structure$param <- rbind(structure$param, structure.rho)
        return(structure)
    }
    
    ## ** identify and name parameters
    level.rho <- NULL
    param.rho <- NULL
    strata.rho <- NULL
    code.rho <- NULL
    code.x.rho <- NULL
    code.y.rho <- NULL
    k.x <- NULL
    k.y <- NULL

        for(iS in 1:n.strata){ ## iS <- 1
            iStrata <- U.strata[iS]
        
            if(n.strata==1){
                iM.level <- M.level
                iX.level <- X.level
            }else{
                iIndex.strata <- M.level[[strata.var]]==iStrata
                iM.level <- M.level[iIndex.strata,regressor.var,drop=FALSE]
                iX.level <- X.level[iIndex.strata,,drop=FALSE]
            }
            ## form all unique pairs of covariate values
            iVec.level <- interaction(iM.level, sep = sep[1])
            iM.pairs <- .unorderedPairs(as.character(iVec.level), distinct = TRUE)
            ## deduce parameters        
            iLevel.rho <- paste0("(",iM.pairs[1,],",",iM.pairs[2,],")")
            if(n.strata>1){
                iLevel.rho <- paste0(iLevel.rho,sep[2],iStrata)
            }
            iParam.rho <- paste0("rho",iLevel.rho)
            iStrata.rho <- rep(iS, length(iParam.rho))
            ## generate the design matrix associated to each time and reduce it to a vector, e.g. time1=1,time2=0,time3=0 --> 100
            iUlp.level <- factor(interaction(as.data.frame(iX.level), sep = "", drop=TRUE), level.cor)
            ## deduce the pairs in term of design matrix 
            iM.pairs.num <- .unorderedPairs(1:length(iVec.level), distinct = TRUE)
            iUlp.diff <- as.character(interaction(as.data.frame(iX.level[iM.pairs.num[2,],]-iX.level[iM.pairs.num[1,],]), sep = sep[1], drop= TRUE))
            ## make sure same correlation coefficient for (X=1,X=0) and (X=0,X=1)

            ## generate code
            iCode.rho <- paste0("D",iUlp.level[iM.pairs.num[1,]],sep[1],iUlp.diff)
            iCode.x.rho <- match(as.character(iUlp.level[iM.pairs.num[1,]]), level.cor)
            iCode.y.rho <- match(as.character(iUlp.level[iM.pairs.num[2,]]), level.cor)
            ## compare to available observations
            iDiff.obs <- as.character(interaction(as.data.frame(XpairPattern$diffU.strata[[iS]]), sep = sep[1], drop= TRUE))
            iMatch.obs <- match(iUlp.diff,iDiff.obs)

            if(any(is.na(iMatch.obs))){
                stop("Something went wrong when generating the correlation parameters for the UN pattern. \n",
                     "Missing parameter(s) compared to what is observed. \n")
            }
            iIndex.keep <- iUlp.diff %in% iDiff.obs
            iMindex.obs <- attr(XpairPattern$LpU.strata[[iS]],"index")[iMatch.obs,,drop=FALSE]

            ## export
            level.rho <- c(level.rho,iLevel.rho[iIndex.keep])
            param.rho <- c(param.rho,iParam.rho[iIndex.keep])
            strata.rho <- c(strata.rho,iStrata.rho[iIndex.keep])
            code.rho <- c(code.rho,iCode.rho[iIndex.keep])
            code.x.rho <- c(code.x.rho,iCode.x.rho[iIndex.keep])
            code.y.rho <- c(code.y.rho,iCode.y.rho[iIndex.keep])
            if(length(param.k)>0){
                k.x <- c(k.x, param.k.obs[iMindex.obs[,1]])
                k.y <- c(k.y, param.k.obs[iMindex.obs[,2]])
            }
        
        
    }

    ## ** update
    structure.rho <- data.frame(name = param.rho,
                                index.strata = strata.rho,
                                type = rep("rho",length=length(param.rho)),
                                index.level = as.numeric(NA),
                                level = level.rho,
                                code = code.rho,
                                code.x = code.x.rho,
                                code.y = code.y.rho,
                                sigma = param.sigma[match(strata.rho,strata.sigma)],
                                k.x = k.x,
                                k.y = k.y,                                  
                                stringsAsFactors = FALSE)
    structure$param <- rbind(structure$param, structure.rho)

    ## ** export
    attr(structure$param, "level.cor") <- level.cor
    return(structure)

        ## ## **** contrast all pairs
        ## for(iC in iCluster2){ ## iC <- iCluster2[1]
        ##     iCindex  <- index.cluster[[iC]][iULpIndex.cor[[iC]]]
        ##     iCX.cor <- X.cor[iCindex,,drop=FALSE]
        ##     iData <- data[iCindex,,drop=FALSE]
        ##     if(NROW(iCX.cor)==1){
        ##         iPair.time <- matrix(1, nrow = 2, ncol = 1)
        ##     }else{
        ##         iPair.time <- ls.pair[[NROW(iCX.cor)]]
        ##     }
        ##     iM <- matrix(lpnCluster.cor[[iC]][iULpIndex.cor[[iC]]][iPair.time], ncol = 2, byrow = TRUE, dimnames = list(NULL, c("x","y")))

        ##     iDF.diff <- as.data.frame(do.call(rbind,lapply(1:NCOL(iPair.time),function(iCol){ ## iCol <- 1

        ##         if(iM[iCol,"x"] < iM[iCol,"y"]){
        ##             ## make sure same correlation coefficient for (X=1,X=0) and (X=0,X=1)
        ##             iCX.cor1 <- iCX.cor[iPair.time[1,iCol],,drop=FALSE]
        ##             iCX.cor2 <- iCX.cor[iPair.time[2,iCol],,drop=FALSE]
        ##         }else if(iM[iCol,"x"] > iM[iCol,"y"]){
        ##             ## make sure same correlation coefficient for (X=1,X=0) and (X=0,X=1)
        ##             iCX.cor1 <- iCX.cor[iPair.time[2,iCol],,drop=FALSE]
        ##             iCX.cor2 <- iCX.cor[iPair.time[1,iCol],,drop=FALSE]
        ##         }else{
        ##             iCX.cor1 <- iCX.cor[min(iPair.time[,iCol]),,drop=FALSE]
        ##             iCX.cor2 <- iCX.cor[max(iPair.time[,iCol]),,drop=FALSE]
        ##         }

        ##             if(heterogeneous>=1){
        ##                 if(all(iCX.cor1==iCX.cor2)){
        ##                     return(cbind("R",iCX.cor1))
        ##                 }else{
        ##                     return(cbind(paste0("D",paste(iCX.cor1,collapse="")),iCX.cor2-iCX.cor1))
        ##                 }
        ##             }else if(heterogeneous>=0){
        ##                 if(all(iCX.cor1==iCX.cor2)){
        ##                     return(matrix(c("R",rep(iStrata,NCOL(iCX.cor1))), nrow = 1, ncol = 1+NCOL(iCX.cor1)))
        ##                 }else{
        ##                     return(matrix(c("D",as.numeric(iCX.cor2!=iCX.cor1)), nrow = 1))
        ##                 }
        ##             }else if(heterogeneous<0){
        ##                 if(all(iCX.cor1==iCX.cor2)){ ## all equal
        ##                     return(matrix(c("R",rep(iStrata,NCOL(iCX.cor1))), nrow = 1, ncol = 1+NCOL(iCX.cor1)))
        ##                 }else{ ## at least one equal 
        ##                     test.common <- setdiff(colnames(iCX.cor1)[which(iCX.cor1==iCX.cor2)], "(Intercept)")
        ##                     if(length(test.common)>0){
        ##                         return(matrix(c("D",as.numeric(iCX.cor2!=iCX.cor1)), nrow = 1))
        ##                     }else{
        ##                         return(matrix(NA,nrow = 1, ncol = 1+NCOL(iCX.cor1)))
        ##                     }
        ##                 }
        ##             }

                
        ##     })))

        ##     iCode <- as.character(interaction(iDF.diff, drop=TRUE))
        ##     ## name difference according to the covariate values
        ##     iName.covcor <- setdiff(names(attr(X.cor,"M.level")),strata.var)
        ##     if(length(iName.covcor)>0){
        ##         iCov <- as.character(interaction(iData[,iName.covcor,drop=FALSE],drop=TRUE))
        ##         index.iUCode <- intersect(which(!duplicated(iCode)),which(!is.na(iCode)))
        ##         iULevel <- stats::setNames(sapply(index.iUCode,function(iCol){
        ##             return(paste0("(",paste(unique(c(iCov[min(iPair.time[,iCol])],iCov[max(iPair.time[,iCol])])),collapse=","),")"))
        ##         }), iCode[index.iUCode])
        ##         iLevel <- iULevel[iCode]
        ##         if(n.strata>1){
        ##             iLevel <- paste0(iLevel,":",U.strata[iStrata])
        ##         }
        ##     }else if(n.strata==1){
        ##         iLevel <- ""
        ##     }else if(n.strata>1){
        ##         iLevel <- paste0(":",U.strata[iStrata])
        ##     }

}

## * skeletonRho.RE
.skeletonRho.RE <- .skeletonRho.CS

## * skeletonRho.TOEPLITZ
.skeletonRho.TOEPLITZ <- function(structure, data, 
                                U.cluster, index.cluster,
                                U.time, index.clusterTime, 
                                U.strata, index.clusterStrata, sep = c(".",":")){

    browser()


    if(toeplitz){
        index.XcolTime <- which(attr(X.cor,"term.labels")==attr(X.cor,"variable")[1])
        if(block){
            index.XcolBlock <- which(attr(X.cor,"term.labels")==attr(X.cor,"variable")[2])
        }else{
            index.XcolBlock <- NULL
        }
    }
if(toeplitz){
                    if(block){
                        if(iCX.cor2[,index.XcolBlock] == iCX.cor1[,index.XcolBlock]){ ## same block
                            if(heterogeneous=="UN"){
                                return(cbind("R",iCX.cor2[,index.XcolBlock],paste(iCX.cor1,collapse=""),iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime]))
                            }else if(heterogeneous=="LAG"){
                                return(cbind("R",iCX.cor2[,index.XcolBlock],abs(iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime])))
                            }else if(heterogeneous=="CS"){
                                return(cbind("R",iCX.cor2[,index.XcolBlock],0))
                            }
                        }else{ ## different block
                            if(heterogeneous=="UN"){
                                if(iCX.cor2[,index.XcolTime]==iCX.cor1[,index.XcolTime]){
                                    return(cbind("D",abs(iCX.cor2[,index.XcolBlock]-iCX.cor1[,index.XcolBlock]),"0",iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime]))
                                }else{
                                    return(cbind("D",abs(iCX.cor2[,index.XcolBlock]-iCX.cor1[,index.XcolBlock]),paste(iCX.cor1,collapse=""),iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime]))
                                }                                
                            }else if(heterogeneous=="LAG"){
                                return(cbind("D",abs(iCX.cor2[,index.XcolBlock]-iCX.cor1[,index.XcolBlock]),abs(iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime])))
                            }else if(heterogeneous=="CS"){
                                return(cbind("D",abs(iCX.cor2[,index.XcolBlock]-iCX.cor1[,index.XcolBlock]),as.numeric(iCX.cor2[,index.XcolTime]!=iCX.cor1[,index.XcolTime])))
                            }
                        }
                    }else{
                        return(cbind("D",abs(iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime])))
                    }

                }
    ## no correlation parameter
    return(structure)

}

## * skeletonRho.UN
.skeletonRho.UN <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata, sep = c(".",":")){

    ## ** extract information
    ## parameters
    index.sigma <- structure$param[structure$param$type=="sigma","index.level"]
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strata.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    index.k <- structure$param[structure$param$type=="k","index.level"]
    param.k <- structure$param[structure$param$type=="k","name"]
    strata.k <- structure$param[structure$param$type=="k","index.strata"]

    ## design matrix (reduce sample size to unique replicates)
    cor.column <- colnames(structure$X$cor)
    XpairPattern <- .pairPatternX(structure$X$cor, lp.obs = structure$X$lp.cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)
    if(length(index.k)>0){
        X.var.k <- rowSums(sweep(structure$X$var[,param.k,drop=FALSE], FUN = "*", STATS = 1:length(param.k), MARGIN = 2))
        if(any(X.var.k>length(param.k))){
            stop("Something went wrong when identifying the variance multiplier parameters. \n")
        }
        param.k.obs <- c(NA, param.k)[X.var.k+1]
    }
    ## variables
    strata.var <- structure$name$strata
    n.strata <- length(U.strata)
    xfactor <- structure$xfactor$cor

    ## levels
    M.level <- attr(structure$X$cor, "M.level")
    if(length(xfactor)>0){
        for(iVar in names(xfactor)){
            M.level[[iVar]] <- factor(M.level[[iVar]], levels = xfactor[[iVar]])
        }
    }
    time.var <- setdiff(colnames(M.level),strata.var) ## non-strata variable
    level.cor <- rownames(XpairPattern$lp2X)
    
    ## others
    formula <- attr(structure$X$cor, "formula")
    X.level <- model.matrix(formula, M.level)[,cor.column,drop=FALSE] ## subset in case some columns have been removed (e.g. strata to be ignored as no pair)
        
    ## ** identify and name parameters
    level.rho <- NULL
    param.rho <- NULL
    strata.rho <- NULL
    code.rho <- NULL
    code.x.rho <- NULL
    code.y.rho <- NULL
    k.x <- NULL
    k.y <- NULL
    
    for(iS in 1:n.strata){ ## iS <- 2
        iStrata <- U.strata[iS]

        if(n.strata==1){
            iM.level <- M.level
            iX.level <- X.level
        }else{
            iIndex.strata <- M.level[[strata.var]]==iStrata
            iM.level <- M.level[iIndex.strata,time.var,drop=FALSE]
            iX.level <- X.level[iIndex.strata,,drop=FALSE]
            if(NROW(iM.level)==1){next} ## case where one strata has no repetition
        }
        ## form all unique pairs of covariate values
        iVec.level <- interaction(iM.level, sep = sep[1], drop = TRUE)
        iM.pairs <- .unorderedPairs(as.character(iVec.level), distinct = TRUE)
        ## deduce parameters        
        iLevel.rho <- paste0("(",iM.pairs[1,],",",iM.pairs[2,],")")
        if(n.strata>1){
            iLevel.rho <- paste0(iLevel.rho,sep[2],iStrata)
        }
        iParam.rho <- paste0("rho",iLevel.rho)
        iStrata.rho <- rep(iS, length(iParam.rho))
        ## generate the design matrix associated to each time and reduce it to a vector, e.g. time1=1,time2=0,time3=0 --> 100
        iUlp.level <- factor(interaction(as.data.frame(iX.level), sep = "", drop=TRUE), level.cor)
        ## deduce the pairs in term of design matrix 
        iM.pairs.num <- .unorderedPairs(1:length(iVec.level), distinct = TRUE)
        iUlp.diff <- as.character(interaction(as.data.frame(iX.level[iM.pairs.num[2,],]-iX.level[iM.pairs.num[1,],]), sep = sep[1], drop= TRUE))
        ## generate code
        iCode.rho <- paste0("D",iUlp.level[iM.pairs.num[1,]],sep[1],iUlp.diff)
        iCode.x.rho <- match(as.character(iUlp.level[iM.pairs.num[1,]]), level.cor)
        iCode.y.rho <- match(as.character(iUlp.level[iM.pairs.num[2,]]), level.cor)
        ## compare to available observations
        iDiff.obs <- as.character(interaction(as.data.frame(XpairPattern$diffU.strata[[iS]]), sep = sep[1], drop= TRUE))
        iMatch.obs <- match(iUlp.diff,iDiff.obs)

        ## iDiff.obs %in% iUlp.diff
        ## if(any(is.na(iMatch.obs))){
        ##     stop("Something went wrong when generating the correlation parameters for the UN pattern. \n",
        ##          "Missing parameter(s) compared to what is observed. \n")
        ## }

        iIndex.keep <- iUlp.diff %in% iDiff.obs
        iMindex.obs <- attr(XpairPattern$LpU.strata[[iS]],"index")[iMatch.obs,,drop=FALSE]

        ## export
        level.rho <- c(level.rho,iLevel.rho[iIndex.keep])
        param.rho <- c(param.rho,iParam.rho[iIndex.keep])
        strata.rho <- c(strata.rho,iStrata.rho[iIndex.keep])
        code.rho <- c(code.rho,iCode.rho[iIndex.keep])
        code.x.rho <- c(code.x.rho,iCode.x.rho[iIndex.keep])
        code.y.rho <- c(code.y.rho,iCode.y.rho[iIndex.keep])
        if(length(param.k)>0){
            k.x <- c(k.x, param.k.obs[iMindex.obs[iIndex.keep,1]])
            k.y <- c(k.y, param.k.obs[iMindex.obs[iIndex.keep,2]])
        }
        
        
    }
    ## ** update
    structure.rho <- data.frame(name = param.rho,
                                index.strata = strata.rho,
                                type = rep("rho",length=length(param.rho)),
                                index.level = as.numeric(NA),
                                level = level.rho,
                                code = code.rho,
                                code.x = code.x.rho,
                                code.y = code.y.rho,
                                sigma = param.sigma[match(strata.rho,strata.sigma)],
                                k.x = k.x,
                                k.y = k.y,                                  
                                stringsAsFactors = FALSE)
    structure$param <- rbind(structure$param, structure.rho)
    rownames(structure$param) <- NULL

    ## ** export
    attr(structure$param, "level.cor") <- level.cor
    return(structure)
}

## * skeletonRho.EXP
.skeletonRho.EXP <- function(structure, data, 
                           U.cluster, index.cluster,
                           U.time, index.clusterTime, 
                           U.strata, index.clusterStrata, sep = c(".",":")){

    browser()
    ## *** param rho
    regressor <- colnames(X.cor)[which(attr(X.cor, "assign") == max(attr(X.cor, "assign")))]

    if(n.strata==1){
        param.rho <- "lambda"
        strata.rho <- 1
        code.rho <- regressor
        level.rho <- ""
        if(structure$heterogeneous){
            param.rho <- c(param.rho,"nugget")
            strata.rho <- c(strata.rho,1)
            code.rho <- c(code.rho,NA)
            level.rho <- c(code.rho,"")
        }
    }else{
        param.rho <- paste0("lambda",U.strata)
        strata.rho <- 1:n.strata
        code.rho <- regressor
        level.rho <- U.strata
        if(structure$heterogeneous){
            param.rho <- c(param.rho,paste0("nugget",U.strata))
            strata.rho <- c(strata.rho,1:n.strata)
            code.rho <- c(code.rho,rep(NA, n.strata))
            level.rho <- c(level.rho,U.strata)
        }
    }

}

## * helper
## ** .pairPatternX
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##'
##' X.cor <- model.matrix(~1, data = gastricbypassL)
##' .pairPatternX(X.cor,
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))))
##' 
##' X.cor <- model.matrix(~visit, data = gastricbypassL)
##' .pairPatternX(X.cor,
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))))
.pairPatternX <- function(object, lp.obs = NULL,
                          U.cluster, index.cluster,
                          U.strata, index.clusterStrata,
                          sep = c("")){

    ## ** normalize user input
    n.strata <- length(U.strata)

    ## code for linear predictor for each observation
    if(is.null(lp.obs)){
        lp.obs <- interaction(as.data.frame(object), sep = "", drop=TRUE)
    }
    ## same but grouped by cluster
    if(!is.null(attr(lp.obs, "cluster"))){
        lp.cluster <- attr(lp.obs, "cluster")
    }else{
        lp.cluster <- lapply(U.cluster, function(iC){ lp.obs[index.cluster[[iC]]] })
    }

    ## ** extract pattern
    ## unique code
    Ulp.obs <- sort(as.character(unique(object)), decreasing = TRUE)
    ## list containing for each strata the index of the clusters with unique linear predictor pattern
    indexLpU.clusterStrata <- lapply(1:n.strata, function(iS){ ## iS <- 2
        iIndex.cluster <- which(index.clusterStrata[U.cluster]==iS)
        return(as.double(iIndex.cluster[which(!duplicated(lp.cluster[iIndex.cluster]))]))
    })
    ## list containing for each strata the unique pairs of linear predictors
    ## (and the position in the dataset of these linear predictors as an attribute)
    LpU.strata <- lapply(1:n.strata, function(iS){ ## iS <- 1
        iM <- base::t(do.call(cbind,lapply(lp.cluster[indexLpU.clusterStrata[[iS]]], .unorderedPairs, distinct = TRUE)))
        iTest <- duplicated(iM)
        iOut <- iM[!iTest,,drop=FALSE]
        attr(iOut,"index") <- base::t(do.call(cbind,lapply(index.cluster[indexLpU.clusterStrata[[iS]]], .unorderedPairs, distinct = TRUE)))[!iTest,,drop=FALSE]
        return(iOut)
    })
    ## design matrix associated to each code
    index.duplicated <- duplicated(lp.obs)
    UX.obs <- object[!index.duplicated,,drop=FALSE]
    rownames(UX.obs) <- lp.obs[!index.duplicated]

    ## list containing for each strata the unique difference in design matrix between pairs
    diffU.strata <- lapply(LpU.strata, function(iS){
        iM <- object[attr(iS,"index")[,2],] - object[attr(iS,"index")[,1],]
        rownames(iM) <- NULL
        return(iM)
    })

    ## ** export
    out <- list(LpU.strata = stats::setNames(LpU.strata, U.strata),
                diffU.strata = stats::setNames(diffU.strata, U.strata),
                indexCluster.LpU.strata = stats::setNames(indexLpU.clusterStrata, U.strata),
                lp2X = UX.obs[levels(lp.obs),,drop=FALSE]
                )
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
### structure-skeletonRho.R ends here
