### findPatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:06) 
## Version: 
## Last-Updated: jul 21 2023 (15:02) 
##           By: Brice Ozenne
##     Update #: 741
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
                            index.clusterStrata, U.strata,
                            sep) UseMethod(".findUpatterns")

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
                           n.time = unname(sapply(pattern2lp.var,length)),
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
            iM <- matrix(nrow = length(index.cluster[[iC]]), ncol = 0)
            attr(iM,"index.cluster") <- iC.all
            return(iM)
        })
    }

    ## ** rename patterns
    if(!is.na(structure$name$strata)){ ## WARNING may not be unique, e.g. in presence of missing values
        Upattern$group <- as.character(U.strata[Upattern$index.strata])
    }

    ## ** export
    structure$pattern <- as.character(pattern.var)
    attr(structure$pattern,"list") <- tapply(1:n.cluster, structure$pattern, base::identity)
    structure$Upattern <- Upattern
    structure$var$Xpattern <- Xpattern.var
    return(structure)
}

## * .findUpatterns.IND
.findUpatterns.IND <- function(structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata){

    sep <- LMMstar.options()$sep["pattern"]

    structure <- .findUpatterns.ID(structure = structure, 
                                   index.clusterTime = index.clusterTime, U.time = U.time,
                                   index.cluster = index.cluster, U.cluster = U.cluster,
                                   index.clusterStrata = index.clusterStrata, U.strata = U.strata)

    ## rename patterns (to add possible extra covariates)
    structure$Upattern$group <- .namePatternCov(Upattern = structure$Upattern, structure = structure, sep = sep)

    ## export
    return(structure)

}


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
    nobs.Upattern.cor <- sapply(index.cluster[indexPatternCluster1.cor],length) 

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
    }
    ## convert linear predictor into correlation parameter
    pattern.pairwise$name <- M.lp2rho[pattern.pairwise$lp.x + n.lp.cor*(pattern.pairwise$lp.y-1)]

    ## remove correlation parameters fixed to 0 (e.g. crossed random effects)
    pattern.pairwise <- pattern.pairwise[pattern.pairwise$name %in% param.rho,,drop=FALSE]
    ## Upattern.pairwise <- pattern.pairwise[!duplicated(pattern.pairwise[,c("pattern","name")]),]

    ## *** summarize correlation patterns
    UpatternCor <- data.frame(name = NA,
                              var = NA,
                              cor = Upattern.cor,
                              index.strata = unname(index.clusterStrata[indexPatternCluster1.cor]),
                              n.time = unname(nobs.Upattern.cor),
                              n.cluster = unname(sapply(ls.indexPatternCluster.cor,length)),
                              param = NA)
    UpatternCor$param <- tapply(pattern.pairwise$name,pattern.pairwise$pattern,unique)

    ## *** characterize each correlation pattern
    Xpattern.cor <- lapply(Upattern.cor, function(iPattern){ ## iPattern <- Upattern.cor[1]

        iCluster <- indexPatternCluster1.cor[[iPattern]]
        iIndex.cluster <- index.cluster[[iCluster]]
        iRep <- length(iIndex.cluster)        
        if(iRep==1){return(NULL)}
        iOut <- X.cor[iIndex.cluster,,drop=FALSE]
        iPattern.pairwise <- pattern.pairwise[pattern.pairwise$pattern == iPattern,]
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
        if(any(iUparam  %in% attr(iOut, "param") == FALSE)){
            for(iParam in setdiff(iUparam,attr(iOut, "param"))){ ## iParam <- iUparam[1]
                attr(iOut, "Omega.cor")[attr(iOut, "indicator.param")[[iParam]]] <- structure.param[structure.param$name==iParam,"constraint"]
            }
        }        
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

    ## ** rename patterns
    if(!is.na(structure$name$strata)){ ## WARNING may not be unique, e.g. in presence of missing values
        Upattern$group <- as.character(U.strata[Upattern$index.strata])
    }

    ## ** export
    structure$pattern <- vec.pattern
    attr(structure$pattern,"list") <- tapply(1:n.cluster, vec.pattern, base::identity)
    structure$Upattern <- Upattern
    attr(structure$Upattern,"sep") <- sep
    structure$var$Xpattern <- Xpattern.var
    structure$cor$Xpattern <- Xpattern.cor
    return(structure)
}

## * .findUpatterns.RE
.findUpatterns.RE <- .findUpatterns.UN

## * .findUpatterns.TOEPLITZ
.findUpatterns.TOEPLITZ <- .findUpatterns.UN

## * .findUpatterns.EXP
.findUpatterns.EXP <- .findUpatterns.UN

## * .findUpatterns.CS
.findUpatterns.CS <- function(structure = structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata){

    structure <- .findUpatterns.UN(structure = structure, 
                                   index.clusterTime = index.clusterTime, U.time = U.time,
                                   index.cluster = index.cluster, U.cluster = U.cluster,
                                   index.clusterStrata = index.clusterStrata, U.strata = U.strata)

    ## rename patterns (to add possible extra covariates)
    structure$Upattern$nameCov <- .namePatternCov(Upattern = structure$Upattern, structure = structure, sep = attr(structure$Upattern,"sep"))

    ## export
    return(structure)

}

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

    ## *** var
    lpObs.var <- nlme::collapse(X.var, sep = sep, as.factor = TRUE)
    lpnObs.var <- as.numeric(lpObs.var)
    lpnCluster.var <- stats::setNames(sapply(U.cluster, function(iC){
        paste(lpnObs.var[index.cluster[[iC]]], collapse=".")
    }), U.cluster)

    name.pattern.var <- sort(unique(lpnCluster.var))
    cluster.pattern.var <- which(!duplicated(lpnCluster.var))
    cluster.pattern.var <- cluster.pattern.var[match(lpnCluster.var[cluster.pattern.var], name.pattern.var)]
    n.pattern.var <- length(name.pattern.var)

    ## *** cor
    if(!is.null(X.cor)){
        lpObs.cor <- nlme::collapse(X.cor, sep = sep, as.factor = TRUE)
        lpnObs.cor <- as.numeric(lpObs.cor)
        lpnCluster.cor0 <- stats::setNames(lapply(U.cluster, function(iC){
            if(length(index.cluster[[iC]])>1){
                return(lpnObs.cor[index.cluster[[iC]]])
            }else{
                return(NULL)
            }
        }), U.cluster)
        lpnCluster.cor <- sapply(lpnCluster.cor0, paste, collapse = ".")
        name.pattern.cor <- unique(lpnCluster.cor)
        name.pattern.cor <- sort(name.pattern.cor)
        
        cluster.pattern.cor <- which(!duplicated(lpnCluster.cor))
        cluster.pattern.cor <- cluster.pattern.cor[match(lpnCluster.cor[cluster.pattern.cor], name.pattern.cor)]
        n.pattern.cor <- length(name.pattern.cor)
    }else{
        name.pattern.cor <- NULL
    }

    ## *** var and cor
    if(!is.null(X.cor)){
        lpnCluster <- paste(as.numeric(factor(lpnCluster.var, name.pattern.var)),
                            as.numeric(factor(lpnCluster.cor, name.pattern.cor)),
                            sep = sep)       
    }else{
        lpnCluster <- as.character(as.numeric(factor(lpnCluster.var, name.pattern.var)))
    }
    test.duplicated <- duplicated(lpnCluster)
    name.pattern <- lpnCluster[!test.duplicated]
    order.pattern <- order(name.pattern)
    name.pattern <- name.pattern[order.pattern]
    n.pattern <- length(name.pattern)

    pattern.cluster <- data.frame(index.cluster = U.cluster,
                                  pattern = factor(name.pattern[match(lpnCluster, name.pattern)], levels = name.pattern))

    ## *** var/cor pattern corresponding to pattern
    pattern.indexVar <- lpnCluster.var[!test.duplicated][order.pattern]
    if(!is.null(X.cor)){
        pattern.indexCor <- lpnCluster.cor[!test.duplicated][order.pattern]
    }else{
        pattern.indexCor <- NA
    }
    Upattern <- data.frame(name = name.pattern,
                           var = pattern.indexVar,
                           cor = pattern.indexCor,
                           index.strata = NA,
                           n.time = NA,
                           time = NA,
                           param = NA,
                           index.cluster = NA,
                           n.cluster = NA)

    Upattern$index.cluster <- stats::setNames(lapply(name.pattern, function(iPattern){which(lpnCluster == iPattern)}), name.pattern)
    Upattern$n.cluster <- sapply(Upattern$index.cluster, length)
    Upattern$n.time <- sapply(Upattern$index.cluster, function(iId){length(index.clusterTime[[iId[1]]])})
    Upattern$time <- lapply(Upattern$index.cluster,function(iC){index.clusterTime[[iC[1]]]})
    Upattern$index.strata <- sapply(Upattern$index.cluster, function(iId){unname(unique(index.clusterStrata[iId]))}, simplify = FALSE)
    attr(Upattern, "level.var") <- name.pattern.var
    attr(Upattern, "level.cor") <- name.pattern.cor
    pattern.cluster$var <- stats::setNames(stats::setNames(Upattern$var, Upattern$name)[pattern.cluster$pattern], pattern.cluster$cluster)
    if(!is.null(X.cor)){
        pattern.cluster$cor <- stats::setNames(stats::setNames(Upattern$cor, Upattern$name)[pattern.cluster$pattern], pattern.cluster$cluster)
    }
    Upattern$param <- stats::setNames(rep(list(param$name), NROW(Upattern)), Upattern$name)

    ## ** characterize each pattern

    ## *** var
    Xpattern.var <- stats::setNames(lapply(name.pattern.var,function(iP){ ## iP <- name.pattern.var[1]
        iC.all <- pattern.cluster$index.cluster[which(pattern.cluster$var==iP)]
        iC <- iC.all[1]
        iX <- X.var[index.cluster[[iC]],,drop=FALSE]
        attr(iX,"index.cluster") <- iC.all
        attr(iX,"index.strata") <- unname(index.clusterStrata[iC])
        return(iX)
    }),name.pattern.var)

    ## *** cor
    if(!is.null(X.cor)){
        UlpnCluster.cor0 <- lpnCluster.cor0[cluster.pattern.cor]
        
        iN.pair <- unique(sapply(UlpnCluster.cor0, length))
        ls.pair <- vector(mode = "list", length = max(iN.pair))

        ls.pair[iN.pair[iN.pair>0]] <- lapply(iN.pair[iN.pair>0], function(iN){unorderedPairs(1:iN, distinct = TRUE)})
        Xpattern.cor <- stats::setNames(lapply(setdiff(name.pattern.cor,""),function(iP){ ## iP <- name.pattern.cor[1]
            ## ignore pattern "", i.e. patterns with no pair
            iC.all <- pattern.cluster$index.cluster[which(pattern.cluster$cor==iP)]
            iC <- iC.all[1]
            iX <- X.cor[index.cluster[[iC]],,drop=FALSE]            
            attr(iX, "index.cluster") <- iC.all
            if(sum(!duplicated(index.clusterTime[iC.all]))==1){
                attr(iX,"index.time") <- index.clusterTime[[iC]]
            }
            attr(iX, "index.strata") <- unname(index.clusterStrata[iC])

            return(iX)
            ## iM <- matrix(0,length(iC.code),length(iC.code)); iM[attr(iX.pairwise, "index.vec2matrix")] <- iX.pairwise %*% 1:NROW(param.rho)
        }),setdiff(name.pattern.cor,""))

    }else{
        Xpattern.cor <- NULL
    }

    ## ** export
    structure$Upattern <- Upattern
    structure$var$Xpattern <- Xpattern.var
    structure$cor$Xpattern <- Xpattern.cor
    structure$pattern <- pattern.cluster
    
    return(structure)
}

## * .namePatternCov
.namePatternCov <- function(Upattern, structure, sep){

    all.cov <- structure$name$strata
    if(all(!is.na(Upattern$var))){
        all.cov <- c(all.cov,structure$name$var[[1]])
    }
    if(all(!is.na(Upattern$cor))){
        all.cov <- c(all.cov,structure$name$cor[[1]])
    }
    all.Ucov <- setdiff(unique(stats::na.omit(all.cov)),c(structure$name$time,attr(structure$name$time,"original")))
    if(length(all.cov)>0){
        out <- sapply(structure$var$pattern2lp, function(iLp){ ## iLp <- 1:4
            iX <- structure$var$lp2data[iLp,all.Ucov,drop=FALSE]
            iName <- sapply(iX, function(iVec){if(length(unique(iVec))>1){NA}else{as.character(iVec[1])}})
            return(paste(na.omit(iName), collapse = sep))
        })        
    }else{
        out <- NULL
    }
    return(out)
}
##----------------------------------------------------------------------
### findPatterns.R ends here
