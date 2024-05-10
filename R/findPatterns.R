### findPatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:06) 
## Version: 
## Last-Updated: maj 10 2024 (11:45) 
##           By: Brice Ozenne
##     Update #: 916
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

`.findUpatterns` <-function(structure, 
                            index.clusterTime, U.time,
                            index.cluster, U.cluster,
                            index.clusterStrata, U.strata){
    UseMethod(".findUpatterns")
}

## * .findUpatterns.ID
.findUpatterns.ID <- function(structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata){

    ## ** extract from object
    X.var <- structure$var$X
    structure.param <- structure$param[structure$param$type %in% c("sigma","k"),]
    param.sigma <- structure.param[structure.param$type=="sigma","name"]
    param.k <- structure.param[structure.param$type=="k","name"]
    n.cluster <- length(U.cluster)

    ## ** summarize variance patterns
    ## variance pattern per cluster
    pattern.var <- structure$var$pattern
    ## linear predictors associated to each variance pattern
    pattern2lp.var <- structure$var$pattern2lp
    ## correspondance between linear predictors and design matrix
    lp2X.var <- structure$var$lp2X
    ## number of variance patterns
    n.pattern.var <- length(pattern2lp.var)
    ## name of the variance pattern
    Upattern.var <- 1:n.pattern.var
    
    Upattern <- data.frame(name = as.character(Upattern.var),
                           var = Upattern.var,
                           cor = NA,
                           index.strata = NA,
                           n.time = unname(lengths(pattern2lp.var)),
                           n.cluster = NA,
                           param  = NA)

    Upattern$index.strata <- unname(tapply(index.clusterStrata, pattern.var, unique)[Upattern$var])
    Upattern$n.cluster <- unname(table(pattern.var)[Upattern$var])

    if(NCOL(X.var)>0){
        ## find param with non-0 value in the design matrix
        Upattern$param <- lapply(1:n.pattern.var, function(iP){ ## iP <- 1
            names(which(colSums(abs(lp2X.var[pattern2lp.var[[iP]],,drop=FALSE]))>0))
        })       
    }else{
        Upattern$param <- vector(mode = "list", length = n.pattern.var)
    }
    
    ## ** characterize each pattern
    if(NCOL(X.var)>0){
        Xpattern.var <- lapply(Upattern$var,function(iP){ ## iP <- Upattern$name[1]
            iC.all <- which(pattern.var==iP)
            iC <- iC.all[1]
            iX <- X.var[index.cluster[[iC]],,drop=FALSE]
            iIndex.param <- which(colSums(abs(iX))!=0)
            iParam.k <- intersect(names(iIndex.param),param.k)
            iParam.sigma <- intersect(names(iIndex.param),param.sigma)
        
            attr(iX,"index.cluster") <- iC.all
            attr(iX,"index.strata") <- unname(index.clusterStrata[iC])
            attr(iX,"param") <- names(iIndex.param)
            attr(iX,"Mindicator.param") <- stats::setNames(lapply(iIndex.param,function(iCol){tcrossprod(iX[,iCol],rep(1,NROW(iX))) + t(tcrossprod(iX[,iCol],rep(1,NROW(iX))))}),
                                                           names(iIndex.param))
            ## attr(iX,"indicator.param") <- lapply(attr(iX,"Mindicator.param"),function(iM){which(iM>0)})
        
            return(iX)
        })
        
    }else{
        Xpattern.var <- lapply(Upattern$var, function(iP){
            iC.all <- which(pattern.var==iP)
            iC <- iC.all[1]
            iM <- matrix(nrow = length(index.cluster[[iC]]), ncol = 0)
            attr(iM,"index.cluster") <- iC.all
            return(iM)
        })
    }

    ## ** export
    structure$pattern <- as.character(pattern.var)
    attr(structure$pattern,"list") <- tapply(1:n.cluster, structure$pattern, base::identity)
    structure$Upattern <- Upattern
    structure$var$Xpattern <- Xpattern.var
    return(structure)
}

## * .findUpatterns.IND
.findUpatterns.IND <- .findUpatterns.ID


## * .findUpatterns.UN
.findUpatterns.UN <- function(structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata){

    sep <- LMMstar.options()$sep["pattern"]

    ## ** identify unique var patterns
    Upatterns.init <- .findUpatterns.IND(structure = structure, 
                                         index.clusterTime = index.clusterTime, U.time = U.time,
                                         index.cluster = index.cluster, U.cluster = U.cluster,
                                         index.clusterStrata = index.clusterStrata, U.strata = U.strata)
    UpatternVar <- Upatterns.init$Upattern
    Xpattern.var <- Upatterns.init$var$Xpattern
    
    ## ** identify unique cor patterns
    X.cor <- structure$cor$X
    structure.param <- structure$param[structure$param$type=="rho",]
    param.rho <- structure.param$name[is.na(structure.param$constraint)] ## index.level=0 indicates correlation fixed to 0, i.e. not a real parameter (e.g. crossed random effects)
    n.cluster <- length(U.cluster)
    if(length(param.rho)==0){return(Upatterns.init)}

    ## *** group clusters by linear predictor patterns
    ## *** (to later work on a single representative cluster per pattern)

    ## linear predictor per observation
    lp.cor <- structure$cor$lp
    n.lp.cor <- NROW(structure$cor$lp2X)
    ## correlation pattern per cluster
    pattern.cor <- structure$cor$pattern
    ## linear predictors  associated to each correlation pattern
    pattern2lp.cor <- structure$cor$pattern2lp
    ## list of cluster per correlation pattern
    ls.indexPatternCluster.cor <- tapply(1:n.cluster,pattern.cor,identity)
    ## representative cluster for each correlation pattern
    indexPatternCluster1.cor <- sapply(ls.indexPatternCluster.cor,"[",1)
    ## name of the correlation pattern
    Upattern.cor <- 1:length(pattern2lp.cor)
    ## number of observations per correlation pattern
    nobs.Upattern.cor <- lengths(index.cluster[indexPatternCluster1.cor]) 

    ## *** form all pairs of observations within representative cluster
    ## *** and deduce correlation parameter
    
    ## identify all pairs
    pattern.pairwise <- do.call(rbind,lapply(indexPatternCluster1.cor[nobs.Upattern.cor>1], function(iC){ ## iC <- 1
        iPair <- base::t(unorderedPairs(1:length(index.cluster[[iC]]), distinct = TRUE))
        iDf <- data.frame(cluster = iC,
                          strata = unname(index.clusterStrata[iC]),
                          pattern = pattern.cor[iC],
                          index.x = iPair[,1],
                          index.y = iPair[,2],
                          obs.x = index.cluster[[iC]][iPair[,1]],
                          obs.y = index.cluster[[iC]][iPair[,2]])
        iDf$lp.x <- lp.cor[iDf$obs.x]
        iDf$lp.y <- lp.cor[iDf$obs.y]
        return(iDf)
    }))
    rownames(pattern.pairwise) <- NULL
        
    ## matrix converting the pair of linear predictors into correlation parameters
    M.lp2rho <- matrix(as.character(NA), nrow = n.lp.cor, ncol = n.lp.cor)
    for(iR in 1:NROW(structure.param)){
        M.lp2rho[structure.param$lp.x[[iR]] + n.lp.cor*(structure.param$lp.y[[iR]]-1)] <- structure.param$name[[iR]]
        M.lp2rho[structure.param$lp.y[[iR]] + n.lp.cor*(structure.param$lp.x[[iR]]-1)] <- structure.param$name[[iR]]
    }
    ## convert linear predictor into correlation parameter
    pattern.pairwise$name <- M.lp2rho[pattern.pairwise$lp.x + n.lp.cor*(pattern.pairwise$lp.y-1)]

    ## remove correlation parameters fixed to 0 (e.g. crossed random effects)
    constraint.pairwise <- pattern.pairwise[pattern.pairwise$name %in% setdiff(structure.param$name,param.rho),,drop=FALSE]
    pattern.pairwise <- pattern.pairwise[pattern.pairwise$name %in% param.rho,,drop=FALSE]
    ## Upattern.pairwise <- pattern.pairwise[!duplicated(pattern.pairwise[,c("pattern","name")]),]

    ## *** summarize correlation patterns
    UpatternCor <- data.frame(name = NA,
                              var = NA,
                              cor = Upattern.cor,
                              index.strata = unname(index.clusterStrata[indexPatternCluster1.cor]),
                              n.time = unname(nobs.Upattern.cor),
                              n.cluster = unname(lengths(ls.indexPatternCluster.cor)),
                              param = NA)
    UpatternCor$param <- tapply(pattern.pairwise$name,factor(pattern.pairwise$pattern,Upattern.cor),unique,simplify=FALSE)
    
    ## *** characterize each correlation pattern
    Xpattern.cor <- lapply(Upattern.cor, function(iPattern){ ## iPattern <- Upattern.cor[1]

        iCluster <- indexPatternCluster1.cor[[iPattern]]
        iIndex.cluster <- index.cluster[[iCluster]]
        iRep <- length(iIndex.cluster)        
        if(iRep==1){return(NULL)}
        iOut <- X.cor[iIndex.cluster,,drop=FALSE]
        iPattern.pairwise <- pattern.pairwise[pattern.pairwise$pattern == iPattern,]
        iConstraint.pairwise <- constraint.pairwise[constraint.pairwise$pattern == iPattern,]
        iUparam <- unique(iPattern.pairwise$name)
        attr(iOut, "index.cluster") <- ls.indexPatternCluster.cor[[iPattern]]
        attr(iOut, "index.time") <- index.clusterTime[[iCluster]]
        attr(iOut, "index.strata") <- index.clusterStrata[[iCluster]]
        attr(iOut, "index.pair") <- rbind(data.frame(row =  iPattern.pairwise$index.x,
                                                     col = iPattern.pairwise$index.y,
                                                     param = iPattern.pairwise$name),
                                          data.frame(row =  iPattern.pairwise$index.y,
                                                     col = iPattern.pairwise$index.x,
                                                     param = iPattern.pairwise$name))
        attr(iOut, "index.vec2matrix") <- attr(iOut, "index.pair")$row + (attr(iOut, "index.pair")$col-1)*iRep
        attr(iOut, "param") <- intersect(param.rho,iUparam) ## ignore parameter constrained to a specific value
        attr(iOut, "indicator.param") <- tapply(attr(iOut, "index.vec2matrix"),attr(iOut, "index.pair")$param,identity,simplify=FALSE)
        ## attr(iOut, "Mindicator.param") <- tapply(attr(iOut, "index.vec2matrix"),attr(iOut, "index.pair")$param,function(iIndex){
        ##     iM <- matrix(0, nrow = iRep, ncol = iRep)
        ##     iM[iIndex] <- 1
        ##     return(iM)
        ## }, simplify=FALSE)
        attr(iOut, "Omega.cor") <- matrix(NA, nrow = iRep, ncol = iRep)
        diag(attr(iOut, "Omega.cor")) <- 1
        ## case with NA, i.e. correlation parameter constrained to be 0
        if(NROW(iConstraint.pairwise)>0){
            iConstraint.param <- stats::setNames(structure.param[!is.na(structure.param$constraint),"constraint"], structure.param[!is.na(structure.param$constraint),"name"])
            attr(iOut, "Omega.cor")[iConstraint.pairwise$index.x + iRep * (iConstraint.pairwise$index.y - 1)] <- unname(iConstraint.param[iConstraint.pairwise$name])
            attr(iOut, "Omega.cor")[iConstraint.pairwise$index.y + iRep * (iConstraint.pairwise$index.x - 1)] <- unname(iConstraint.param[iConstraint.pairwise$name])
        }

        ## table(attr(iOut, "Omega.cor")[unlist(attr(iOut, "indicator.param"))], useNA = "always")
        ## table(attr(iOut, "Omega.cor")[-unlist(attr(iOut, "indicator.param"))], useNA = "always")

        return(iOut)
    })

    ## ** joint variance and correlation patterns
    vec.pattern <- paste0(structure$var$pattern,sep,structure$cor$pattern)
    test.Upattern <- !duplicated(vec.pattern)
    index.Upattern <- which(test.Upattern)[order(vec.pattern[test.Upattern])] ## re-order index so patterns appears in alphabetic order (instead of order driven by the dataset)

    Upattern <- data.frame(name = vec.pattern[index.Upattern],
                           var = structure$var$pattern[index.Upattern],
                           cor = structure$cor$pattern[index.Upattern],
                           index.strata = NA,
                           n.time = NA,
                           n.cluster = NA,
                           param = NA)

    matchVar <- match(Upattern$var,UpatternVar$var)
    matchCor <- match(Upattern$cor,UpatternCor$cor)
    Upattern$index.strata <- UpatternCor[matchCor,"index.strata"]
    Upattern$n.cluster <- UpatternCor[matchCor,"n.cluster"]
    Upattern$n.time <- UpatternCor[matchCor,"n.time"]
    Upattern$param <- unname(mapply(x = UpatternVar$param[matchVar], y = UpatternCor$param[matchCor], FUN = function(x,y){unname(c(x,y))}, SIMPLIFY = FALSE))
    
    ## ** export
    structure$pattern <- vec.pattern
    attr(structure$pattern,"list") <- tapply(1:n.cluster, vec.pattern, base::identity)
    structure$Upattern <- Upattern
    attr(structure$Upattern,"sep") <- sep
    structure$var$Xpattern <- Xpattern.var
    structure$cor$Xpattern <- Xpattern.cor
    return(structure)
}

## * .findUpatterns.CS
.findUpatterns.CS <- .findUpatterns.UN

## * .findUpatterns.RE
.findUpatterns.RE <- .findUpatterns.UN

## * .findUpatterns.TOEPLITZ
.findUpatterns.TOEPLITZ <- .findUpatterns.UN

## * .findUpatterns.EXP
.findUpatterns.EXP <- .findUpatterns.UN


## * .findUpatterns_CUSTOM
.findUpatterns.CUSTOM <- function(structure, 
                                  index.clusterTime, U.time,
                                  index.cluster, U.cluster,
                                  index.clusterStrata, U.strata){

    X.var <- structure$var$X
    X.cor <- structure$cor$X
    param <- structure$param
    sep <- LMMstar.options()$sep["pattern"]

    ## ** identify unique patterns
    cluster.pattern.var <- structure$var$pattern
    cluster.pattern.cor <- structure$cor$pattern
    
    if(!is.null(cluster.pattern.cor)){
        cluster.pattern <- paste(cluster.pattern.var,
                                 cluster.pattern.cor,
                                 sep = sep)       
    }else{
        cluster.pattern <- cluster.pattern.var
        cluster.pattern.cor <- rep(NA, length(cluster.pattern.var))
    }
    index.Upattern <- which(!duplicated(cluster.pattern))
    
    Upattern <- data.frame(name = cluster.pattern[index.Upattern],
                           var = cluster.pattern.var[index.Upattern],
                           cor = cluster.pattern.cor[index.Upattern],
                           index.strata = index.clusterStrata[index.Upattern],
                           n.time = lengths(structure$var$pattern2lp),
                           n.cluster = NA,
                           param = NA)

    Upattern$n.cluster <- unname(table(cluster.pattern)[Upattern$name])
    Upattern$param <- unname(tapply(param$name, param$index.strata, base::identity, simplify = FALSE))
    Upattern <- Upattern[order(Upattern$index.strata,Upattern$n.time),,drop=FALSE]
    rownames(Upattern) <- NULL

    ## ** X pattern
    cluster.pattern.num <- as.numeric(factor(cluster.pattern, Upattern$name))
    attr(cluster.pattern.num,"list") <- tapply(1:length(cluster.pattern.num), cluster.pattern, base::identity)[Upattern$name]

    cluster.Upattern.var <- which(!duplicated(cluster.pattern.var))
    cluster.Upattern.var <- cluster.Upattern.var[order(cluster.pattern.var[cluster.Upattern.var])]

    Xpattern.var <- lapply(cluster.Upattern.var, function(iC){ ## iC <- cluster.Upattern.var[1]
        iP <- cluster.pattern.var[iC]
        iOut <- structure$var$lp2data[structure$var$pattern2lp[[iP]],,drop=FALSE]
        attr(iOut,"index.cluster") <- which(cluster.pattern.var==iP)
        attr(iOut,"index.strata") <- index.clusterStrata[iC]
        attr(iOut,"param") <- param[param$index.strata  == index.clusterStrata[iC] & param$type == "sigma","name"]
        return(iOut)
    })

    if(any(!is.na(cluster.pattern.cor))){
    cluster.Upattern.cor <- which(!duplicated(cluster.pattern.cor))
    cluster.Upattern.cor <- cluster.Upattern.cor[order(cluster.pattern.cor[cluster.Upattern.cor])]

    Xpattern.cor <- lapply(cluster.Upattern.cor, function(iC){ ## iC <- cluster.Upattern.cor[1]
        iP <- cluster.pattern.cor[iC]
        iOut <- structure$cor$lp2data[structure$cor$pattern2lp[[iP]],,drop=FALSE]
        attr(iOut,"index.cluster") <- which(cluster.pattern.cor==iP)
        attr(iOut,"index.strata") <- index.clusterStrata[iC]
        attr(iOut,"param") <- param[param$index.strata  == index.clusterStrata[iC] & param$type == "rho","name"]
        return(iOut)
    })
    }else{
        Xpattern.cor <- NULL
    }

    ## ** export
    structure$pattern <- cluster.pattern.num
    structure$Upattern <- Upattern
    attr(structure$Upattern,"sep") <- sep
    structure$var$Xpattern <- Xpattern.var
    structure$cor$Xpattern <- Xpattern.cor
    return(structure)
}

##----------------------------------------------------------------------
### findPatterns.R ends here
