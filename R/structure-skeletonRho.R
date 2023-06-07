### structure-skeletonRho.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 11 2023 (13:27) 
## Version: 
## Last-Updated: jun  7 2023 (09:46) 
##           By: Brice Ozenne
##     Update #: 435
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
    cor.column <- colnames(structure$X$cor)
    XpairPattern <- .pairPatternX(structure$X$cor, data = data, lp.obs = structure$X$lp.cor,
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
    reg.var <- setdiff(structure$name$cor[[1]], NA)
    n.reg <- length(reg.var)
    n.strata <- length(U.strata)

    ## levels
    X.level <- XpairPattern$lp2X
    M.level <- XpairPattern$lp2data
    strata.level <- XpairPattern$index.strata.lp2X
    level.cor <- rownames(X.level)

    ## ** special case (no strata, no covariate)
    if(length(reg.var)==0){
        vec.strataRho <- which(sapply(XpairPattern$LpU.strata,length)>0)
        if(n.strata==1){
            level.strataRho <- ""
        }else{
            level.strataRho <- paste0(":",U.strata[vec.strataRho])
        }
        structure.rho <- data.frame(name = paste0("rho",level.strataRho),
                                    index.strata = vec.strataRho,
                                    type = "rho",
                                    constraint = as.numeric(NA),
                                    level = level.strataRho,
                                    code = paste0("R.",vec.strataRho,".",vec.strataRho),
                                    code.x = NA,
                                    code.y = NA,
                                    sigma = param.sigma[match(vec.strataRho, strata.sigma)],
                                    k.x = NA,
                                    k.y = NA,                                  
                                    stringsAsFactors = FALSE)        
        structure.rho$code.x <- as.list(vec.strataRho)
        structure.rho$code.y <- as.list(vec.strataRho)
        structure$param <- rbind(structure$param, structure.rho)
        rownames(structure$param) <- NULL
        attr(structure$param, "Xcode.xy") <- setNames(1:NROW(XpairPattern$lp2X),rownames(XpairPattern$lp2X))
        return(structure)
    }
    
    
    ## ** identify and name parameters
    
    ## *** combine linear predictor accross strata
    diffLp <- do.call(rbind,XpairPattern$diffU.strata) ## difference in linear predictor index for pair of observation
    indexLp <- do.call(rbind,XpairPattern$LpU.strata) ## linear predictor index for each observation of each pair
    lp.x <- XpairPattern$lp2X[indexLp[,1],,drop=FALSE] ## design matrox for one observation of each pair of observations
    data.x <- XpairPattern$lp2data[indexLp[,1],,drop=FALSE] ## data for one observation of each pair of observations
    data.y <- XpairPattern$lp2data[indexLp[,2],,drop=FALSE] ## data for the other observation of each pair of observations
    strataLp <- XpairPattern$index.strata.lp2X[indexLp[,1]] ## strata for pair of observation (would be the same with indexLp[,2])
    n.Lp <- NROW(lp.x)

    ## *** identify pairs with equal or non-equal linear predictors
    index.equal <- which(rowSums(diffLp!=0)==0)
    index.unequal <- setdiff(1:n.Lp, index.equal)
    index.0 <- NULL ## cells in the design matrix constrained to be 0
    
    ## *** generate code
    code.rho <- rep(NA, length = n.Lp)
    level.rho <- rep(NA, length = n.Lp)
    if(length(index.equal)>0){
        if(structure$type == "homogeneous"){            
            code.rho[index.equal] <- paste("R",strataLp[index.equal],sep=sep[2])
            level.rho[index.equal] <- ""
        }else if(structure$type == "heterogeneous"){
            code.rho[index.equal] <- paste("R",lp.x[index.equal],sep=sep[2])
            level.rho[index.equal] <- paste0("(",interaction(as.data.frame(data.x[index.equal,,drop=FALSE]), sep = sep[1], drop = TRUE),")")
        } 
    }
    if(length(index.unequal)>0){
        if(structure$type == "homogeneous"){
            ## contrast at the design matrix level
            test.n0 <- 1*(data.x[index.unequal,reg.var,drop=FALSE]!=data.y[index.unequal,reg.var,drop=FALSE])
            code.rho[index.unequal] <- paste("D",strataLp[index.unequal],interaction(as.data.frame(test.n0),sep=sep[1]),sep=sep[2])
            
            index.unequal.red <- which(!duplicated(code.rho[index.unequal]))
            code2level <- stats::setNames(sapply(index.unequal.red, function(iUnequal){ ## iUnequal <- index.unequal[1]
                iData.x <- data.x[iUnequal,,drop=FALSE]
                iData.y <- data.y[iUnequal,,drop=FALSE]
                paste0("(",paste(iData.x[iData.x!=iData.y], collapse = sep[1]),",",paste(iData.y[iData.x!=iData.y], collapse = sep[1]),")")                
            }),code.rho[index.unequal.red])
            level.rho[index.unequal] <- code2level[code.rho[index.unequal]]
            
        }else if(structure$type == "heterogeneous"){
            code.rho[index.unequal] <- paste("D",lp.x[index.unequal],interaction(as.data.frame(diffLp[index.unequal,,drop=FALSE]), sep = sep[1], drop = TRUE),sep=sep[2])
            level.rho[index.unequal] <- paste0("(",interaction(as.data.frame(data.x[index.unequal,,drop=FALSE]), sep = sep[1], drop = TRUE),
                                               ",",interaction(as.data.frame(data.y[index.unequal,,drop=FALSE]), sep = sep[1], drop = TRUE),
                                               ")")
        }
        if(length(unique(structure$group.type))>1){
            sumGroup.0 <- do.call(cbind,tapply(names(structure$group.type), structure$group.type, function(iCol){
                rowSums(diffLp[index.unequal,iCol,drop=FALSE]==0)
            }, simplify = FALSE))
            
            if(n.strata==1){
                index.0 <- index.unequal[which(rowSums(sumGroup.0!=0)==0)]
            }else{
                if(length(strata.var)>1){
                    stop("Cannot handle multiple strata variables with this covariance structure. \n")
                }
                col.strata <- attr(structure$X$cor,"M.level")[strata.var]
                index.0 <- unlist(lapply(1:n.strata, FUN = function(iS){ ## iS <- 1
                    iCol <- rownames(col.strata)[col.strata==U.strata[iS]] ## columns in the design matrix corresponding 
                    iIndex <- index.unequal[strataLp[index.unequal]==iS]
                    iIndex[which(rowSums(diffLp[iIndex,iCol,drop=FALSE]==0)==0)]
                }))
            }            
        }
    }
    test.rho <- !duplicated(code.rho)
    code.Urho <- code.rho[test.rho]
    
    ## ***  name parameters
    level.Urho <-  level.rho[test.rho]
    strata.Urho <- XpairPattern$index.strata.lp2X[indexLp[test.rho,1]]
    if(n.strata>1){
        level.Urho <- paste0(level.Urho,sep[2],strata.Urho)
    }

    ## ***  retrive k
    if(length(index.k)>0){
        indexObs <- do.call(rbind,lapply(XpairPattern$LpU.strata,attr,"index"))
        k.x <- param.k.obs[indexObs[,1]]
        k.y <- param.k.obs[indexObs[,2]]
    }else{
        k.x <- NA
        k.y <- NA
    }

    ## ***  collect    
    structure.rho <- data.frame(name = paste0("rho",level.Urho),
                                index.strata = strata.Urho,
                                type = rep("rho",length=length(level.Urho)),
                                constraint = as.numeric(NA),
                                level = level.Urho,
                                code = code.Urho,
                                code.x = NA,
                                code.y = NA,
                                sigma = param.sigma[match(strata.Urho,strata.sigma)],
                                k.x = k.x[test.rho],
                                k.y = k.y[test.rho],                                  
                                stringsAsFactors = FALSE)

    ## save information about correlation parameters fixed to 0 (e.g. crossed random effects)
    if(length(index.0)>0){
        Ucode.rho0 <- unique(code.rho[index.0])
        if(length(intersect(Ucode.rho0,code.rho[-index.0]))){
            stop("Something went wrong with the constraint when identifying the correlation parameters. \n")
        }
        structure.rho$constraint[structure.rho$code %in% Ucode.rho0] <- 0
    }
    ## ensure proper ordering
    code.x <- tapply(indexLp[,1], code.rho, base::identity, simplify = FALSE)
    code.y <- tapply(indexLp[,2], code.rho, base::identity, simplify = FALSE)
    structure.rho$code.x[match(names(code.x),structure.rho$code)] <- code.x
    structure.rho$code.y[match(names(code.y),structure.rho$code)] <- code.y

    ## ** export
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)
    return(structure)
}

## * skeletonRho.RE
.skeletonRho.RE <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata, sep = c(".",":")){

    ## ** identify variance and correlation parameters
    structure <- .skeletonRho.CS(structure = structure, data = data, 
                                 U.cluster = U.cluster, index.cluster = index.cluster,
                                 U.time = U.time, index.clusterTime = index.clusterTime, 
                                 U.strata = U.strata, index.clusterStrata = index.clusterStrata, sep = sep)

    ## ** relate correlation parameters to random effects
    if(!is.null(structure$ranef$hierarchy)){
        browser()
    }
    
    ## ** export
    return(structure)
}

## * skeletonRho.TOEPLITZ
.skeletonRho.TOEPLITZ <- function(structure, data, 
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
    XpairPattern <- .pairPatternX(structure$X$cor, data = data, lp.obs = structure$X$lp.cor,
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
    reg.var <- setdiff(structure$name$cor[[1]], NA)
    n.strata <- length(U.strata)

    ## levels
    X.level <- XpairPattern$lp2X
    M.level <- XpairPattern$lp2data
    strata.level <- XpairPattern$index.strata.lp2X
    level.cor <- rownames(X.level)

    ## ** special case (no strata, no covariate)
    if(length(reg.var)==0){
        vec.strataRho <- which(sapply(XpairPattern$LpU.strata,length)>0)
        if(n.strata==1){
            level.strataRho <- ""
        }else{
            level.strataRho <- paste0(":",U.strata[vec.strataRho])
        }
        structure.rho <- data.frame(name = paste0("rho",level.strataRho),
                                    index.strata = vec.strataRho,
                                    type = "rho",
                                    index.level = as.numeric(NA),
                                    level = level.strataRho,
                                    code = paste0("R.",vec.strataRho,".",vec.strataRho),
                                    code.x = NA,
                                    code.y = NA,
                                    sigma = param.sigma[match(vec.strataRho, strata.sigma)],
                                    k.x = NA,
                                    k.y = NA,                                  
                                    stringsAsFactors = FALSE)        
        structure.rho$code.x <- as.list(vec.strataRho)
        structure.rho$code.y <- as.list(vec.strataRho)
        structure$param <- rbind(structure$param, structure.rho)
        rownames(structure$param) <- NULL
        attr(structure$param, "Xcode.xy") <- setNames(1:NROW(XpairPattern$lp2X),rownames(XpairPattern$lp2X))
        return(structure)
    }
    
    ## ** identify and name parameters
    
    ## *** extract linear predictor
    indexLp <- do.call(rbind,XpairPattern$LpU.strata) ## linear predictor index for each observation of each pair
    lp.x <- XpairPattern$lp2X[indexLp[,1],,drop=FALSE] ## design matrox for one observation of each pair of observations
    data.x <- XpairPattern$lp2data[indexLp[,1],,drop=FALSE] ## data for one observation of each pair of observations
    data.y <- XpairPattern$lp2data[indexLp[,2],,drop=FALSE] ## data for the other observation of each pair of observations
    diffLp <- do.call(rbind,XpairPattern$diffU.strata) ## difference in linear predictor index for pair of observation
    strataLp <- XpairPattern$index.strata.lp2X[indexLp[,1]] ## strata for pair of observation (would be the same with indexLp[,2])
    n.Lp <- NROW(indexLp)
    n.reg <- length(reg.var)
        
    ## *** generate code
    code.rho <- rep(NA, length = n.Lp)
    level.rho <- rep(NA, length = n.Lp)
    index.equal <- which(rowSums(diffLp!=0)==0)
    if(length(index.equal)>0){        
        if(structure$heterogeneous){
            code.rho[index.equal] <- paste("R",lp.x[index.equal],sep=sep[2])
            level.rho[index.equal] <- paste0("(",interaction(as.data.frame(data.x[index.equal,,drop=FALSE]), sep = sep[1], drop = TRUE),")")
        }else{
            code.rho[index.equal] <- paste("R",strataLp[index.equal],sep=sep[2])
            level.rho[index.equal] <- ""
        }
    }
        
    index.unequal <- setdiff(1:n.Lp, index.equal)
    if(length(index.unequal)>0){
        if(structure$heterogeneous){
            code.rho[index.unequal] <- paste("D",lp.x[index.unequal],diffLp[index.unequal],sep=sep[2])
            level.rho[index.unequal] <- paste0("(",interaction(as.data.frame(data.x[index.unequal,,drop=FALSE]), sep = sep[1], drop = TRUE),
                                              ",",interaction(as.data.frame(data.y[index.unequal,,drop=FALSE]), sep = sep[1], drop = TRUE),
                                              ")")
        }else{
            code.rho[index.unequal] <- paste("D",strataLp[index.unequal],as.numeric(diffLp[index.unequal]!=0),sep=sep[2])

            iData.x <- data.x[index.unequal,,drop=FALSE]
            iData.x[,diffLp[index.unequal]==0] <- ""
            iData.y <- data.y[index.unequal,,drop=FALSE]
            iData.y[,diffLp[index.unequal]==0] <- ""
            level.rho[index.unequal] <- paste0("(",interaction(as.data.frame(iData.x), sep = sep[1], drop = TRUE),
                                              ",",interaction(as.data.frame(iData.y), sep = sep[1], drop = TRUE),
                                              ")")
        }            
        
    }
    test.rho <- !duplicated(code.rho)
    code.Urho <- code.rho[test.rho]

    
    ## ***  name parameters
    level.Urho <-  level.rho[test.rho]
    strata.Urho <- XpairPattern$index.strata.lp2X[indexLp[test.rho,1]]
    if(n.strata>1){
        level.Urho <- paste0(level.Urho,sep[2],strata.rho)
    }

    ## ***  retrive k
    if(length(index.k)>0){
        indexObs <- do.call(rbind,lapply(XpairPattern$LpU.strata,attr,"index"))
        k.x <- param.k.obs[indexObs[,1]]
        k.y <- param.k.obs[indexObs[,2]]
    }else{
        k.x <- NA
        k.y <- NA
    }
    
    ## ***  collect
    structure.rho <- data.frame(name = paste0("rho",level.Urho),
                                index.strata = strata.Urho,
                                type = rep("rho",length=length(level.Urho)),
                                index.level = as.numeric(NA),
                                level = level.Urho,
                                code = code.Urho,
                                code.x = tapply(indexLp[,1],code.rho,base::identity, simplify = FALSE),
                                code.y = tapply(indexLp[,2],code.rho,base::identity, simplify = FALSE),
                                sigma = param.sigma[match(strata.Urho,strata.sigma)],
                                k.x = k.x,
                                k.y = k.y,                                  
                                stringsAsFactors = FALSE)
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)

    ## ** export
    attr(structure$param, "Xcode.xy") <- setNames(1:NROW(XpairPattern$lp2X),rownames(XpairPattern$lp2X))
    return(structure)


##     if(toeplitz){
##         index.XcolTime <- which(attr(X.cor,"term.labels")==attr(X.cor,"variable")[1])
##         if(block){
##             index.XcolBlock <- which(attr(X.cor,"term.labels")==attr(X.cor,"variable")[2])
##         }else{
##             index.XcolBlock <- NULL
##         }
##     }
## if(toeplitz){
##                     if(block){
##                         if(iCX.cor2[,index.XcolBlock] == iCX.cor1[,index.XcolBlock]){ ## same block
##                             if(heterogeneous=="UN"){
##                                 return(cbind("R",iCX.cor2[,index.XcolBlock],paste(iCX.cor1,collapse=""),iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime]))
##                             }else if(heterogeneous=="LAG"){
##                                 return(cbind("R",iCX.cor2[,index.XcolBlock],abs(iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime])))
##                             }else if(heterogeneous=="CS"){
##                                 return(cbind("R",iCX.cor2[,index.XcolBlock],0))
##                             }
##                         }else{ ## different block
##                             if(heterogeneous=="UN"){
##                                 if(iCX.cor2[,index.XcolTime]==iCX.cor1[,index.XcolTime]){
##                                     return(cbind("D",abs(iCX.cor2[,index.XcolBlock]-iCX.cor1[,index.XcolBlock]),"0",iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime]))
##                                 }else{
##                                     return(cbind("D",abs(iCX.cor2[,index.XcolBlock]-iCX.cor1[,index.XcolBlock]),paste(iCX.cor1,collapse=""),iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime]))
##                                 }                                
##                             }else if(heterogeneous=="LAG"){
##                                 return(cbind("D",abs(iCX.cor2[,index.XcolBlock]-iCX.cor1[,index.XcolBlock]),abs(iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime])))
##                             }else if(heterogeneous=="CS"){
##                                 return(cbind("D",abs(iCX.cor2[,index.XcolBlock]-iCX.cor1[,index.XcolBlock]),as.numeric(iCX.cor2[,index.XcolTime]!=iCX.cor1[,index.XcolTime])))
##                             }
##                         }
##                     }else{
##                         return(cbind("D",abs(iCX.cor2[,index.XcolTime]-iCX.cor1[,index.XcolTime])))
##                     }

##                 }
##     ## no correlation parameter
##     return(structure)

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
    XpairPattern <- .pairPatternX(structure$X$cor, data = data, lp.obs = structure$X$lp.cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)
    if(length(index.k)>0){
        param.k.obs <- structure$param$name[match(structure$X$lp.var, structure$param$code)]
        if(any(is.na(param.k.obs))){
            stop("Something went wrong when identifying the variance parameters. \n")
        }
        param.k.obs[param.k.obs %in% param.sigma] <- NA
    }
    ## variables
    time.var <- structure$name$cor[[1]]
    strata.var <- structure$name$strata
    n.strata <- length(U.strata)

    ## levels
    X.level <- XpairPattern$lp2X
    M.level <- XpairPattern$lp2data
    strata.level <- XpairPattern$index.strata.lp2X
    level.cor <- rownames(X.level)
    
    ## ** identify and name parameters
    ## name parameters
    indexLp <- do.call(rbind,XpairPattern$LpU.strata)
    level.rho <- paste0("(",XpairPattern$lp2data[indexLp[,1],time.var],",",XpairPattern$lp2data[indexLp[,2],time.var],")")
    strata.rho <- XpairPattern$index.strata.lp2X[indexLp[,1]]
    if(n.strata>1){
        level.rho <- paste0(level.rho,sep[2],strata.rho)
    }

    ## generate code
    diffLp <- do.call(rbind,XpairPattern$diffU.strata)
    code.rho <- paste("D",strata.rho,interaction(as.data.frame(diffLp),drop=TRUE,sep=sep[1]), sep = sep[2])

    ## retrive k
    if(length(index.k)>0){
        indexObs <- do.call(rbind,lapply(XpairPattern$LpU.strata,attr,"index"))
        k.x <- param.k.obs[indexObs[,1]]
        k.y <- param.k.obs[indexObs[,2]]
    }else{
        k.x <- NA
        k.y <- NA
    }

    ## collect
    structure.rho <- data.frame(name = paste0("rho",level.rho),
                                index.strata = strata.rho,
                                type = rep("rho",length=length(level.rho)),
                                constraint = as.numeric(NA),
                                level = level.rho,
                                code = code.rho,
                                code.x = NA,
                                code.y = NA,
                                sigma = param.sigma[match(strata.rho,strata.sigma)],
                                k.x = k.x,
                                k.y = k.y,                                  
                                stringsAsFactors = FALSE)
    structure.rho$code.x <- as.list(indexLp[,1])
    structure.rho$code.y <- as.list(indexLp[,2])
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)

    ## ** export    
    attr(structure$param, "Xcode.xy") <- setNames(1:NROW(XpairPattern$lp2X),rownames(XpairPattern$lp2X))
    return(structure)
}

## * skeletonRho.EXP
.skeletonRho.EXP <- function(structure, data, 
                           U.cluster, index.cluster,
                           U.time, index.clusterTime, 
                           U.strata, index.clusterStrata, sep = c(".",":")){

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
##' @return A list with the following elements \itemize{
##' \item LpU.strata [list of matrices] linear predictor (and position of the observations in the attribute index) used for the pairwise differences.
##' \item diffU.strata [list of matrices] pairwise difference between the linear predictor index
##' \item indexCluster.LpU.strata [list of vectors] index of the clusters involved in the linear predictor, by strata.
##' \item lp2data [data.frame] dataset (columns) associated to each linear predictor (rows)
##' \item lp2X [matrix] design matrix (columns) associated to each linear predictor (rows)
##' \item index.strata.lp2X [vector] strata associated with each lp level.
##' }
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
.pairPatternX <- function(object, data, lp.obs = NULL,
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
    ## list containing for each strata the index of the clusters with unique linear predictor pattern
    indexLpU.clusterStrata <- lapply(1:n.strata, function(iS){ ## iS <- 1
        iIndex.cluster <- which(index.clusterStrata[U.cluster]==iS)
        return(as.double(iIndex.cluster[which(!duplicated(lp.cluster[iIndex.cluster]))]))
    })
    ## list containing for each strata the unique pairs of linear predictors
    ## (and the position in the dataset of these linear predictors as an attribute)
    LpU.strata <- lapply(1:n.strata, function(iS){ ## iS <- 1
        iN.cluster <- length(indexLpU.clusterStrata[[iS]])
        iLp.cluster <- vector(mode = "list", length = iN.cluster)
        iIndex.cluster <- vector(mode = "list", length = iN.cluster)
        iTable.cluster <- table(factor(levels = levels(lp.obs)))
        ## select linear predictor, remove triplicated levels, and only keep levels not already seen
        for(iC in 1:iN.cluster){ ## iC <- 2
            iCluster <- indexLpU.clusterStrata[[iS]][iC]
            iIndex <- which(!triplicated(lp.cluster[[iCluster]]))
            iTable <- table(lp.cluster[[iCluster]][iIndex])
            if(any(iTable>iTable.cluster)){
                iLp.cluster[[iC]] <- lp.cluster[[iCluster]][iIndex]
                iIndex.cluster[[iC]] <- index.cluster[[iCluster]][iIndex]
                iTable.cluster <- pmax(iTable,iTable.cluster)
            }            
        }
        iIndex.cluster <- iIndex.cluster[sapply(iLp.cluster,length)>0]
        iLp.cluster <- iLp.cluster[sapply(iLp.cluster,length)>0]

        ## form all pairs
        iM <- base::t(do.call(cbind,lapply(iLp.cluster, .unorderedPairs, distinct = TRUE)))
        iTest <- duplicated(iM)
        iOut <- iM[!iTest,,drop=FALSE]
        attr(iOut,"index") <- base::t(do.call(cbind,lapply(iIndex.cluster, .unorderedPairs, distinct = TRUE)))[!iTest,,drop=FALSE]
        return(iOut)
    })
    ## design matrix associated to each code
    index.duplicated <- duplicated(lp.obs)
    UX.obs <- object[!index.duplicated,,drop=FALSE]
    rownames(UX.obs) <- lp.obs[!index.duplicated]

    ## list containing for each strata the unique difference in design matrix between pairs
    diffU.strata <- lapply(LpU.strata, function(iS){
        iM <- object[attr(iS,"index")[,2],,drop=FALSE] - object[attr(iS,"index")[,1],,drop=FALSE]
        rownames(iM) <- NULL
        return(iM)
    })

    ## data corresponding to each linear predictor
    M.lp2index <- do.call(rbind,lapply(LpU.strata, function(iM){
        rbind(data.frame(index.lp = iM[,1,drop=FALSE], index = attr(iM,"index")[,1,drop=FALSE]),
              data.frame(index.lp = iM[,2,drop=FALSE], index = attr(iM,"index")[,2,drop=FALSE]))
    }))
    M.Ulp2index <- M.lp2index[!duplicated(M.lp2index$index.lp),,drop=FALSE]
    M.Ulp2index$lp <- factor(levels(lp.obs)[M.Ulp2index$index.lp], levels = levels(lp.obs))
    M.Ulp2index$index.strata <- index.clusterStrata[attr(index.cluster,"vectorwise")[M.Ulp2index$index]]

    ## ** export    
    out <- list(LpU.strata = stats::setNames(LpU.strata, U.strata),
                diffU.strata = stats::setNames(diffU.strata, U.strata),
                indexCluster.LpU.strata = stats::setNames(indexLpU.clusterStrata, U.strata),
                lp2data = data[M.Ulp2index[match(levels(lp.obs),M.Ulp2index$lp),"index"],attr(object,"variable"),drop=FALSE],
                lp2X = UX.obs[levels(lp.obs),,drop=FALSE],
                index.strata.lp2X = M.Ulp2index[match(levels(lp.obs),M.Ulp2index$lp),"index.strata"]
                )
    rownames(out$lp2data) <- levels(lp.obs)
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
