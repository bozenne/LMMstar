### findPatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:06) 
## Version: 
## Last-Updated: apr 10 2026 (17:59) 
##           By: Brice Ozenne
##     Update #: 998
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
                            sep){
    UseMethod(".findUpatterns")
}

## * .findUpatterns.ID
.findUpatterns.ID <- function(structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata,
                              sep){

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
                           n.time = NA,
                           index.cluster = NA,
                           param  = NA)

    Upattern$index.strata <- unname(tapply(index.clusterStrata, pattern.var, unique)[Upattern$var])
    Upattern$index.cluster <- unname(tapply(U.cluster, pattern.var, unique)[Upattern$var])
    Upattern$n.time <- lengths(pattern2lp.var) ## some patterns may be the same but for different combinations of times (e.g. CS for 3 timepoints but could be 1:3 or 2:4)

    if(NCOL(X.var)>0){
        ## find param with non-0 value in the design matrix
        Upattern$param <- lapply(1:n.pattern.var, function(iP){ ## iP <- 1
            names(which(colSums(abs(lp2X.var[pattern2lp.var[[iP]],,drop=FALSE]))>0))
        })       
    }else{
        Upattern$param <- vector(mode = "list", length = n.pattern.var)
    }

    ## ** characterize each variance pattern
    Xpattern.var <- lapply(1:NROW(Upattern), function(iP){ ## iP <- 1
        iLp <- structure$var$pattern2lp[[iP]]
        iSigma <- apply(structure$var$lp2X[iLp,param.sigma,drop=FALSE], MARGIN = 1, function(iRow){param.sigma[iRow>0]})
        if(length(unique(iSigma))!=1){
            stop("Something went wrong when identifying the sigma parameters in the residual covariance patterns. \n",
                 "Consider contacting the package manager with a reproducible example. \n")
        }
        if(length(param.k)>0){
            iK <- apply(structure$var$lp2X[iLp,param.k,drop=FALSE], MARGIN = 1, function(iRow){c(param.k[iRow>0],"one")[1]})
            iOut <- array(NA_character_, dim = c(rep(length(iLp),2),3), dimnames = list(NULL,NULL,c("sigma","k1","k2")))
            iOut[,,"sigma"] <- unique(iSigma)
            iOut[,,"k1"] <- matrix(iK, byrow = TRUE, nrow = length(iLp), ncol = length(iLp))
            iOut[,,"k2"] <- matrix(iK, byrow = FALSE, nrow = length(iLp), ncol = length(iLp))
        }else{
            iOut <- array(unique(iSigma), dim = c(rep(length(iLp),2),1), dimnames = list(NULL,NULL,"sigma"))
        }
        return(iOut)
        
    })

    ## ** export
    structure$pattern <- as.character(pattern.var)
    attr(structure$pattern,"list") <- tapply(1:n.cluster, structure$pattern, base::identity)
    structure$Upattern <- Upattern
    structure$var$Xpattern <- Xpattern.var
    return(structure)
}

## * .findUpatterns.IND
.findUpatterns.IND <- .findUpatterns.ID


## * .findUpatterns.CS
.findUpatterns.CS <- function(structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata,
                              sep){

    ## ** identify unique var patterns
    Upatterns.init <- .findUpatterns.IND(structure = structure, 
                                         index.clusterTime = index.clusterTime, U.time = U.time,
                                         index.cluster = index.cluster, U.cluster = U.cluster,
                                         index.clusterStrata = index.clusterStrata, U.strata = U.strata)
    
    UpatternVar <- Upatterns.init$Upattern
    Xpattern.var <- Upatterns.init$var$Xpattern
    if(sum("rho" %in% structure$param$type & is.na(structure$param$constraint)) == 0){return(Upatterns.init)} ## no 'free' correlation parameter

    ## ** identify unique cor patterns
    n.cluster <- length(U.cluster)
    X.cor <- structure$cor$X
    structure.param <- structure$param[structure$param$type=="rho",]
    ## from linear predictor to correlation parameter
    lp2rho <- do.call(c,lapply(1:NROW(structure.param$name), function(iRow){ ## iRow <- 1
        stats::setNames(rep(structure.param$name[iRow], length(structure.param$code[[iRow]])), as.character(structure.param$code[[iRow]]))
    }))

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
        iDf$lp <- paste(lp.cor[iDf$obs.x], lp.cor[iDf$obs.y], sep=sep["rho.name"])
        return(iDf)
    }))
    rownames(pattern.pairwise) <- NULL
    pattern.pairwise$param.rho <- lp2rho[pattern.pairwise$lp]

    ## *** summarize correlation patterns
    UpatternCor <- data.frame(name = NA,
                              var = NA,
                              cor = Upattern.cor,
                              index.strata = NA,
                              n.time = NA,
                              index.cluster = NA,
                              param = NA)
    UpatternCor$index.strata <- unname(tapply(index.clusterStrata, pattern.cor, unique)[UpatternCor$cor])
    UpatternCor$index.cluster <- unname(tapply(U.cluster, pattern.cor, unique)[UpatternCor$cor])
    UpatternCor$n.time <- lengths(pattern2lp.cor)

    ## use intersect instead of unique to keep common ordering of the model parameters, i.e. avoid to have rho1, rho2 in one pattern and rho2, rho1 in the other
    UpatternCor$param <- tapply(pattern.pairwise$param.rho, INDEX = factor(pattern.pairwise$pattern,Upattern.cor),
                                FUN = function(iVec){intersect(structure.param$name,iVec)},simplify=FALSE) 

    ## *** characterize each correlation pattern
    Xpattern.cor <- lapply(1:NROW(UpatternCor), function(iP){ ## iP <- 1
        iN.time <- UpatternCor$n.time[iP]
        iPattern.pairwise <- pattern.pairwise[pattern.pairwise$pattern == UpatternCor$cor[iP],,drop=FALSE]
        iM <- matrix(NA_character_, nrow = iN.time, ncol = iN.time)
        iM[iPattern.pairwise[,"index.x"] + iN.time*(iPattern.pairwise[,"index.y"]-1)] <- iPattern.pairwise[,"param.rho"]
        iM[iPattern.pairwise[,"index.y"] + iN.time*(iPattern.pairwise[,"index.x"]-1)] <- iPattern.pairwise[,"param.rho"]
        diag(iM) <- "one"

        iOut <- array(NA_character_, dim = c(iN.time, iN.time, 1))
        iOut[,,1] <- iM
        return(iOut)
    })

    ## ** joint variance and correlation patterns
    vec.pattern <- paste0(structure$var$pattern,sep["pattern"],structure$cor$pattern)
    test.Upattern <- !duplicated(vec.pattern)
    index.Upattern <- which(test.Upattern)[order(vec.pattern[test.Upattern])] ## re-order index so patterns appears in alphabetic order (instead of order driven by the dataset)

    Upattern <- data.frame(name = vec.pattern[index.Upattern],
                           var = structure$var$pattern[index.Upattern],
                           cor = structure$cor$pattern[index.Upattern],
                           index.strata = NA,
                           n.time = NA,
                           index.cluster = NA,
                           param = NA)

    matchVar <- match(Upattern$var,UpatternVar$var)
    matchCor <- match(Upattern$cor,UpatternCor$cor)
    Upattern$index.strata <- UpatternCor[matchCor,"index.strata"]
    Upattern$n.time <- UpatternCor[matchCor,"n.time"]
    Upattern$index.cluster <- tapply(U.cluster, factor(vec.pattern, levels = Upattern$name), identity, simplify = FALSE)
    Upattern$param <- unname(mapply(x = UpatternVar$param[matchVar], y = UpatternCor$param[matchCor], FUN = function(x,y){unname(c(x,y))}, SIMPLIFY = FALSE))

    ## ** export
    structure$pattern <- vec.pattern
    structure$Upattern <- Upattern
    attr(structure$Upattern,"sep") <- sep["pattern"]
    structure$var$Xpattern <- Xpattern.var
    structure$cor$Xpattern <- Xpattern.cor
    return(structure)
}

## * .findUpatterns.RE
.findUpatterns.RE <- .findUpatterns.CS

## * .findUpatterns.UN
.findUpatterns.UN <- .findUpatterns.CS

## * .findUpatterns.DUN
.findUpatterns.DUN <- .findUpatterns.CS

## * .findUpatterns.TOEPLITZ
.findUpatterns.TOEPLITZ <- .findUpatterns.UN

## * .findUpatterns.EXP
.findUpatterns.EXP <- .findUpatterns.UN


## * .findUpatterns_CUSTOM
.findUpatterns.CUSTOM <- function(structure, 
                                  index.clusterTime, U.time,
                                  index.cluster, U.cluster,
                                  index.clusterStrata, U.strata,
                                  sep){

    X.var <- structure$var$X
    X.cor <- structure$cor$X
    param <- structure$param

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
