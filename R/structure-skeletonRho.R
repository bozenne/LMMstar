### structure-skeletonRho.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 11 2023 (13:27) 
## Version: 
## Last-Updated: jul 13 2023 (11:01) 
##           By: Brice Ozenne
##     Update #: 610
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

    param.k <- structure$param[structure$param$type=="k","name"]
    strata.k <- structure$param[structure$param$type=="k","index.strata"]

    ## design matrix (reduce sample size to unique replicates)
    XpairPattern <- .pairPatternX(structure$cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    ## k parameters
    param.k.augmented <- structure$param$name[match(rownames(structure$var$lp2X)[structure$var$lp],structure$param$code)]
    param.k.augmented[param.k.augmented %in% param.sigma] <- NA

    ## variables
    strata.var <- structure$name$strata
    reg.var <- setdiff(structure$name$cor[[1]], NA)
    n.reg <- length(reg.var)
    n.strata <- length(U.strata)

    ## levels
    M.level.cor <- structure$cor$lp2data
    lp2X.cor <- structure$cor$lp2X
    lp2data.cor <- structure$cor$lp2data
    
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
                                    lp.x = NA,
                                    lp.y = NA,
                                    sigma = param.sigma[match(vec.strataRho, strata.sigma)],
                                    k.x = NA,
                                    k.y = NA,                                  
                                    stringsAsFactors = FALSE)        
        structure.rho$lp.x <- as.list(vec.strataRho)
        structure.rho$lp.y <- as.list(vec.strataRho)
        structure$param <- rbind(structure$param, structure.rho)
        rownames(structure$param) <- NULL
        return(structure)
    }
    
    
    ## ** identify and name parameters

    ## *** combine linear predictor accross strata
    diffLp <- do.call(rbind,XpairPattern$diffU.strata) ## difference in linear predictor index for pair of observation
    indexLp <- do.call(rbind,XpairPattern$LpU.strata) ## linear predictor index for each observation of each pair
    lp.x <- lp2X.cor[indexLp[,1],,drop=FALSE] ## design matrox for one observation of each pair of observations
    data.x <- lp2data.cor[indexLp[,1],,drop=FALSE] ## data for one observation of each pair of observations
    data.y <- lp2data.cor[indexLp[,2],,drop=FALSE] ## data for the other observation of each pair of observations
    strataLp <- index.clusterStrata[attr(index.cluster,"vectorwise")[indexLp[,1]]] ## strata for pair of observation (would be the same with indexLp[,2])
    n.Lp <- NROW(lp.x)

    ## *** identify pairs with equal or non-equal linear predictors
    index.equal <- which(rowSums(diffLp!=0)==0)
    index.unequal <- setdiff(1:n.Lp, index.equal)
    index.0 <- NULL ## cells in the design matrix constrained to be 0

    ## *** generate code
    code.rho <- rep(NA, length = n.Lp)
    level.rho <- rep(NA, length = n.Lp)

    ## pairs containing the same linear predictor
    if(length(index.equal)>0){
        if(structure$type == "homogeneous"){            
            code.rho[index.equal] <- paste("R",strataLp[index.equal],sep=sep[2])
            level.rho[index.equal] <- ""
        }else if(structure$type == "heterogeneous"){
            code.rho[index.equal] <- paste("R",lp.x[index.equal],sep=sep[2])
            level.rho[index.equal] <- paste0("(",nlme::collapse(data.x[index.equal,,drop=FALSE], sep = sep[1], as.factor = FALSE),")")
        } 
    }

    ## pairs with distinct linear predictors
    if(length(index.unequal)>0){
        if(structure$type == "homogeneous"){
            ## contrast at the design matrix level
            test.n0 <- 1*(data.x[index.unequal,reg.var,drop=FALSE]!=data.y[index.unequal,reg.var,drop=FALSE])
            code.rho[index.unequal] <- paste("D",strataLp[index.unequal],nlme::collapse(test.n0,sep=sep[1],as.factor=FALSE),sep=sep[2])
            
            index.unequal.red <- index.unequal[!duplicated(code.rho[index.unequal])]
            code2level <- stats::setNames(sapply(index.unequal.red, function(iUnequal){ ## iUnequal <- index.unequal.red[1]
                iData.x <- data.x[iUnequal,,drop=FALSE]
                iData.y <- data.y[iUnequal,,drop=FALSE]
                iOut <- paste0("(",paste(iData.x[iData.x!=iData.y], collapse = sep[1]),",",paste(iData.y[iData.x!=iData.y], collapse = sep[1]),")")
                return(iOut)
            }),code.rho[index.unequal.red])
            level.rho[index.unequal] <- code2level[code.rho[index.unequal]]
            
        }else if(structure$type == "heterogeneous"){
            code.rho[index.unequal] <- paste("D",lp.x[index.unequal],nlme::collapse(diffLp[index.unequal,,drop=FALSE], sep = sep[1], as.factor = FALSE),sep=sep[2])
            level.rho[index.unequal] <- paste0("(",nlme::collapse(data.x[index.unequal,,drop=FALSE], sep = sep[1], as.factor = FALSE),
                                               ",",nlme::collapse(data.y[index.unequal,,drop=FALSE], sep = sep[1], as.factor = FALSE),
                                               ")")
        }

        ## constrain independence in the case of different groups of correlation structure
        ## i.e. remove some of the correlation parameters
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
                col.strata <- M.level.cor[strata.var]
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
    strata.Urho <- strataLp[test.rho]
    if(n.strata>1){
        level.Urho <- paste0(level.Urho,sep[2],strata.Urho)
    }

    ## ***  retrive k
    if(length(param.k)>0){
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
                                lp.x = NA,
                                lp.y = NA,
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
    structure.rho$lp.x[match(names(code.x),structure.rho$code)] <- code.x
    structure.rho$lp.y[match(names(code.y),structure.rho$code)] <- code.y

    ## ** export
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)
    return(structure)
}

## * .skeletonRho.RE
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
    n.strata <- length(U.strata)

    ## *** extract hierarchy
    hierarchy <- structure$ranef$hierarchy
    name.hierarchy <- unname(unlist(hierarchy))
    n.hierarchy <- length(name.hierarchy) 
    if(any(duplicated(name.hierarchy))){
        stop("Duplicated name(s) for the random effects: \"",paste(name.hierarchy[duplicated(name.hierarchy)], collapse = "\", \""),"\". \n",
             "Cannot associate correlation parameters and random effects. \n")
    }

    if(n.hierarchy==1){
        param.rho <- structure$param[structure$param$type=="rho" & is.na(structure$param$constraint),,drop=FALSE]
        structure$ranef$param <- list(matrix(param.rho$name[param.rho$index.strata], nrow = n.hierarchy, ncol = n.strata,
                                             dimnames = list(name.hierarchy,U.strata)))
    }else{
        ## prepare output
        structure$ranef$param <- lapply(hierarchy, function(iH){
            matrix(as.character(NA), nrow = length(iH), ncol = n.strata,
                   dimnames = list(iH,U.strata))
        })
        
        ## *** convert from index of the linear predictor to design matrix
        M.level.cor <- attr(structure$cor$X,"M.level")
        lp2X.cor <- structure$cor$lp2X
        var.cluster <- attr(structure$name$cluster,"original")
        if(!is.null(var.cluster) && var.cluster %in% name.hierarchy && var.cluster %in% colnames(lp2X.cor) == FALSE){
            ## add cluster variable when missing
            lp2X.cor <- cbind(1,lp2X.cor)
            colnames(lp2X.cor)[1] <- var.cluster
            M.level.cor <- cbind(TRUE,M.level.cor)
            colnames(M.level.cor)[1] <- var.cluster
        }
        ## split hierarchy by strata (since the variable are name time:strata0 and time:strata1 instead of time)
        if(n.strata==1){
            name.hierarchy.strata <- list(stats::setNames(name.hierarchy,name.hierarchy))
        }else{
            name.hierarchy.strata <- lapply(U.strata, function(iStrata){ ## iStrata <- U.strata[1]
                iM.level <- M.level.cor[M.level.cor[[structure$name$strata]]==iStrata,setdiff(name.hierarchy, var.cluster),drop=FALSE]
                if(var.cluster %in% name.hierarchy){
                    iOut <- stats::setNames(c(var.cluster,rownames(iM.level)),
                                            c(var.cluster,colnames(iM.level)[which(iM.level==TRUE, arr.ind = TRUE)[,"col"]]))
                }else{
                    iOut <- stats::setNames(rownames(iM.level),
                                            colnames(iM.level)[which(iM.level==TRUE, arr.ind = TRUE)[,"col"]])
                }
                return(iOut)
            })
        }

        ## *** find for each strata and hierarchy of variable the correlation parameter with
        ## - constant value within the variables of the hierarchy
        ## - variable values for the variables outside of the hierarchy
        index.hierarchy <- unname(unlist(lapply(1:length(hierarchy), function(iH){rep(iH, length(hierarchy[[iH]]))})))

        for(iS in 1:n.strata){ ## iS <- 1
            iParam.rho <- structure$param[structure$param$type=="rho" & is.na(structure$param$constraint) & structure$param$index.strata==iS,,drop=FALSE]
            iName.rho <- iParam.rho$name
            iLevel.rho <- iParam.rho$level
            iCodeX.rho <- stats::setNames(iParam.rho$lp.x, iName.rho)
            iCodeY.rho <- stats::setNames(iParam.rho$lp.y, iName.rho)
            
            iEqual <- do.call(rbind,lapply(iName.rho, function(iRho){ ## iRho <- iName.rho[1]
                colSums(lp2X.cor[iCodeX.rho[[iRho]],name.hierarchy.strata[[iS]],drop=FALSE] - lp2X.cor[iCodeY.rho[[iRho]],name.hierarchy.strata[[iS]],drop=FALSE] != 0)==0
            }))
            rownames(iEqual) <- iName.rho

            for(iH in 1:n.hierarchy){ ## iH <- 1
                iHierarchy <- hierarchy[[index.hierarchy[iH]]]
                iPosition <- which(iHierarchy==name.hierarchy[iH])
                iVarIdentical <- iHierarchy[1:iPosition]
                iVarDifferent <- setdiff(structure$ranef$vars, iVarIdentical)
                iVarIdentical2 <- name.hierarchy.strata[[iS]][iVarIdentical]
                iVarDifferent2 <- name.hierarchy.strata[[iS]][iVarDifferent]
                if(length(iVarDifferent)==0){
                    iTest <- rowSums(iEqual[,iVarIdentical2,drop=FALSE]==FALSE)
                }else{
                    iTest <- rowSums(iEqual[,iVarIdentical2,drop=FALSE]==FALSE) + rowSums(iEqual[,iVarDifferent2,drop=FALSE]==TRUE)
                }
                iRho <- names(which(iTest==0))
                if(length(iRho)!=1){
                    stop("Something when wrong when associating the correlation parameters and random effects. \n",
                         "hierarchy: ",paste(cumhierarchy[[iH]], collapse="\", \""),"\n",
                         "correlation parameter(s): \"",paste(iRho, collapse="\", \""),"\"\n")
                }
                structure$ranef$param[[index.hierarchy[iH]]][iPosition,iS] <- iRho
            }
        }
    }

    ## ** export
    return(structure)
}

## * .skeletonRho.TOEPLITZ
.skeletonRho.TOEPLITZ <- function(structure, data, 
                                  U.cluster, index.cluster,
                                  U.time, index.clusterTime, 
                                  U.strata, index.clusterStrata, sep = c(".",":")){

    ## ** extract information
    ## parameters
    index.sigma <- structure$param[structure$param$type=="sigma","index.level"]
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strata.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    param.k <- structure$param[structure$param$type=="k","name"]
    strata.k <- structure$param[structure$param$type=="k","index.strata"]

    ## design matrix (reduce sample size to unique replicates)
    XpairPattern <- .pairPatternX(structure$cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    if(length(param.k)>0){
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
                                    lp.x = NA,
                                    lp.y = NA,
                                    sigma = param.sigma[match(vec.strataRho, strata.sigma)],
                                    k.x = NA,
                                    k.y = NA,                                  
                                    stringsAsFactors = FALSE)        
        structure.rho$lp.x <- as.list(vec.strataRho)
        structure.rho$lp.y <- as.list(vec.strataRho)
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
            level.rho[index.equal] <- paste0("(",nlme::collapse(data.x[index.equal,,drop=FALSE], sep = sep[1], as.factor = FALSE),")")
        }else{
            code.rho[index.equal] <- paste("R",strataLp[index.equal],sep=sep[2])
            level.rho[index.equal] <- ""
        }
    }
        
    index.unequal <- setdiff(1:n.Lp, index.equal)
    if(length(index.unequal)>0){
        if(structure$heterogeneous){
            code.rho[index.unequal] <- paste("D",lp.x[index.unequal],diffLp[index.unequal],sep=sep[2])
            level.rho[index.unequal] <- paste0("(",nlme::collapse(data.x[index.unequal,,drop=FALSE], sep = sep[1], as.factor = FALSE),
                                              ",",nlme::collapse(data.y[index.unequal,,drop=FALSE], sep = sep[1], as.factor = FALSE),
                                              ")")
        }else{
            code.rho[index.unequal] <- paste("D",strataLp[index.unequal],as.numeric(diffLp[index.unequal]!=0),sep=sep[2])

            iData.x <- data.x[index.unequal,,drop=FALSE]
            iData.x[,diffLp[index.unequal]==0] <- ""
            iData.y <- data.y[index.unequal,,drop=FALSE]
            iData.y[,diffLp[index.unequal]==0] <- ""
            level.rho[index.unequal] <- paste0("(",nlme::collapse(iData.x, sep = sep[1], as.factor = FALSE),
                                              ",",nlme::collapse(iData.y, sep = sep[1], as.factor = FALSE),
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
    if(length(param.k)>0){
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
                                lp.x = tapply(indexLp[,1],code.rho,base::identity, simplify = FALSE),
                                lp.y = tapply(indexLp[,2],code.rho,base::identity, simplify = FALSE),
                                sigma = param.sigma[match(strata.Urho,strata.sigma)],
                                k.x = k.x,
                                k.y = k.y,                                  
                                stringsAsFactors = FALSE)
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)

    ## ** export
    return(structure)
}

## * .skeletonRho.UN
.skeletonRho.UN <- function(structure, data, 
                            U.cluster, index.cluster,
                            U.time, index.clusterTime, 
                            U.strata, index.clusterStrata, sep = c(".",":")){

    ## ** extract information
    ## parameters
    index.sigma <- structure$param[structure$param$type=="sigma","index.level"]
    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strata.sigma <- structure$param[structure$param$type=="sigma","index.strata"]

    param.k <- structure$param[structure$param$type=="k","index.level"]
    param.k <- structure$param[structure$param$type=="k","name"]
    strata.k <- structure$param[structure$param$type=="k","index.strata"]

    ## design matrix (reduce sample size to unique replicates)
    XpairPattern <- .pairPatternX(structure$cor,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    ## k parameters
    param.k.augmented <- structure$param$name[match(rownames(structure$var$lp2X)[structure$var$lp],structure$param$code)]
    param.k.augmented[param.k.augmented %in% param.sigma] <- NA

    ## variables
    time.var <- structure$name$cor[[1]]
    strata.var <- structure$name$strata
    n.strata <- length(U.strata)

    ## levels
    M.level <- structure$cor$lp2data

    ## ** identify and name parameters
    ## name parameters
    indexLp <- do.call(rbind,XpairPattern$LpU.strata)
    level.rho <- paste0("(",M.level[indexLp[,1],time.var],",",M.level[indexLp[,2],time.var],")")
    if(n.strata==1){
        strata.rho <- U.strata
    }else if(n.strata>1){
        strata.rho <- M.level[indexLp[,1],strata.var]
        level.rho <- paste0(level.rho,sep[2],strata.rho)        
    }
    ## generate code
    diffLp <- do.call(rbind,XpairPattern$diffU.strata)
    code.rho <- paste("D",strata.rho,nlme::collapse(diffLp, as.factor=FALSE, sep=sep[1]), sep = sep[2])

    ## retrive k
    if(length(param.k)>0){
        indexObs <- do.call(rbind,lapply(XpairPattern$LpU.strata,attr,"index"))
        k.x <- param.k.augmented[indexObs[,1]]
        k.y <- param.k.augmented[indexObs[,2]]
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
                                lp.x = NA,
                                lp.y = NA,
                                sigma = param.sigma[match(strata.rho,strata.sigma)],
                                k.x = k.x,
                                k.y = k.y,                                  
                                stringsAsFactors = FALSE)

    structure.rho$lp.x <- as.list(indexLp[,1])
    structure.rho$lp.y <- as.list(indexLp[,2])
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)

    ## ** export    
    return(structure)
}

## * .skeletonRho.EXP
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
##' @noRd
##' @return A list with the following elements \itemize{
##' \item LpU.strata [list of matrices] linear predictor (and position of the observations in the attribute index) used for the pairwise differences.
##' \item diffU.strata [list of matrices] pairwise difference between the linear predictor index
##' \item indexCluster.LpU.strata [list of vectors] index of the clusters involved in the linear predictor, by strata.
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
##' .pairPatternX(list(X=X.cor),
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))))
.pairPatternX <- function(object, 
                          U.cluster, index.cluster,
                          U.strata, index.clusterStrata,
                          sep = c("")){

    ## ** normalize user input
    object.X <- object$X
    object.lp <- object$lp
    object.pattern <- object$pattern
    object.pattern2lp <- object$pattern2lp
    n.strata <- length(U.strata)

    ## code for linear predictor for each observation
    if(is.null(object.lp)){
        lp.obs <- nlme::collapse(object.X, sep = "", as.factor = TRUE)
    }else{
        lp.obs <- object.lp
    }

    ## same but grouped by cluster
    if(is.null(object.pattern)){
        pattern <- lapply(U.cluster, function(iC){ lp.obs[index.cluster[[iC]]] })
    }else{
        pattern <- object.pattern
    }

    ## ** extract pattern
    ## list containing for each strata the index of the clusters with unique linear predictor pattern
    indexLpU.clusterStrata <- lapply(1:n.strata, function(iS){ ## iS <- 1
        iIndex.cluster <- which(index.clusterStrata[U.cluster]==iS)
        return(as.double(iIndex.cluster[which(!duplicated(pattern[iIndex.cluster]))]))
    })
    ## list containing for each strata the unique pairs of linear predictors
    ## (and the position in the dataset of these linear predictors as an attribute)
    LpU.strata <- lapply(1:n.strata, function(iS){ ## iS <- 1
        iN.cluster <- length(indexLpU.clusterStrata[[iS]])
        iPattern <- vector(mode = "list", length = iN.cluster)
        iIndex.cluster <- vector(mode = "list", length = iN.cluster)
        iTable.cluster <- table(factor(levels = sort(unique(lp.obs))))
        ## select linear predictor, remove triplicated levels, and only keep levels not already seen
        for(iC in 1:iN.cluster){ ## iC <- 1
            iCluster <- indexLpU.clusterStrata[[iS]][iC]
            iLp <- object.pattern2lp[[pattern[[iCluster]]]]
            iIndex <- which(!triplicated(iLp))
            iTable <- table(iLp[iIndex])
            if(any(iTable>iTable.cluster[names(iTable)])){
                iPattern[[iC]] <- iLp[iIndex]
                iIndex.cluster[[iC]] <- index.cluster[[iCluster]][iIndex]
                iTable.cluster[names(iTable)] <- pmax(iTable,iTable.cluster[names(iTable)])
            }            
        }
        iIndex.cluster <- iIndex.cluster[sapply(iPattern,length)>0]
        iPattern <- iPattern[sapply(iPattern,length)>0]

        ## form all pairs
        iM <- base::t(do.call(cbind,lapply(iPattern, unorderedPairs, distinct = TRUE)))
        iTest <- duplicated(iM)
        iOut <- iM[!iTest,,drop=FALSE]
        attr(iOut,"index") <- base::t(do.call(cbind,lapply(iIndex.cluster, unorderedPairs, distinct = TRUE)))[!iTest,,drop=FALSE]
        return(iOut)
    })

    ## list containing for each strata the unique difference in design matrix between pairs
    diffU.strata <- lapply(LpU.strata, function(iS){
        iM <- object.X[attr(iS,"index")[,2],,drop=FALSE] - object.X[attr(iS,"index")[,1],,drop=FALSE]
        rownames(iM) <- NULL
        return(iM)
    })

    ## ** export    
    out <- list(LpU.strata = stats::setNames(LpU.strata, U.strata),
                diffU.strata = stats::setNames(diffU.strata, U.strata),
                indexCluster.LpU.strata = stats::setNames(indexLpU.clusterStrata, U.strata)
                )
    return(out)
}

## ** .colnameOrder
##' @description Reorder the column names such that the strata variable is at the end
##' @noRd
##' 
##' @examples
##' df <- data.frame(day = c(1, 1, 2, 2, 1, 1, 2, 2),
##'                  gender = c("1", "1", "1", "1", "0", "0", "0", "0"),
##'                  session = c(1, 2, 3, 4, 1, 2, 3, 4))
##' X <- stats::model.matrix(~0 + session:gender + day:gender, df)
##' colnames(X) ## "session:gender0" "session:gender1" "gender0:day" "gender1:day"
##'
##' X2 <- .model.matrix_regularize(~0 + session:gender + day:gender, df, augmodel = TRUE)
##' colnames(X2) ## "session:gender0" "session:gender1" "gender0:day" "gender1:day"
##' 
##' X3 <- .colnameOrder(X2, strata.var = "gender", n.strata = 2)
##' colnames(X3) ## "session:gender0" "session:gender1" "day:gender0" "day:gender1"    
##'
##' stats::model.matrix(~ 0+(day+session):gender, df[order(df$gender),])

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
        rownames(attr(X,"M.level")) <- X.newname
        names(attr(X,"ls.level")) <- X.newname
    }
    
    return(X)    
}


##----------------------------------------------------------------------
### structure-skeletonRho.R ends here
