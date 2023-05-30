### findPatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:06) 
## Version: 
## Last-Updated: maj 30 2023 (19:38) 
##           By: Brice Ozenne
##     Update #: 552
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
                              index.clusterStrata, U.strata,
                              sep = "."){


    ## ** extract from object
    X.var <- structure$X$var
    structure.param <- structure$param[structure$param$type %in% c("sigma","k"),]
    param.sigma <- structure.param[structure.param$type=="sigma","name"]
    param.k <- structure.param[structure.param$type=="k","name"]
    Ulp.var <- structure.param$code[structure.param$index.level]
    n.cluster <- length(U.cluster)

    ## ** associate a pattern to each cluster
    if(NCOL(X.var)>0){
        ## linear predictor for the variance associated to each observation [full]
        lpObs.var <- structure$X$lp.var
        if(!all(levels(lpObs.var) %in% Ulp.var)){
            warning("Something went wrong when extracting the variance patterns. \n")
        }
        ## linear predictor for the variance associated to each observation [integer]
        lpnObs.var <- as.numeric(factor(lpObs.var, levels = Ulp.var))
        ## combine linear predictors for the variance within cluster
        lpnCluster.var <- stats::setNames(sapply(U.cluster, function(iC){
            paste(lpnObs.var[index.cluster[[iC]]], collapse = sep)
        }), U.cluster)

        ## find unique values, i.e. variance pattern
        name.pattern.var <- sort(unique(lpnCluster.var))
        ## variance pattern of each cluster
        lpnCluster.factor <- factor(lpnCluster.var, name.pattern.var)
        lpnCluster.char <- as.character(lpnCluster.factor)
        lpnCluster.num <- as.numeric(lpnCluster.factor)
        
        pattern.cluster <- data.frame(index.cluster = 1:n.cluster,
                                      pattern = as.factor(lpnCluster.num),
                                      var = lpnCluster.char)

    }else{ ## no variance
        ## variance pattern as a function of the cluster size
        pattern.cluster <- data.frame(index.cluster = names(index.cluster),
                                      pattern = sapply(index.cluster,length),
                                      var = paste0("length",sapply(index.cluster,length)))
       
    }

    ## ** summarize patterns
    index.unique <- which(!duplicated(pattern.cluster$pattern)) ## index of individuals with distinct patterns
    index.unique.order <- index.unique[order(pattern.cluster$pattern[index.unique])] ## re-order index so patterns appears in alphabetic order (instead of order driven by the dataset)
    n.pattern.var <- length(name.pattern.var)

    Upattern <- data.frame(name = as.character(pattern.cluster$pattern[index.unique.order]),
                           var = as.character(pattern.cluster$var[index.unique.order]),
                           cor = NA,
                           index.strata = unname(index.clusterStrata[index.unique.order]),
                           n.time = unname(sapply(index.clusterTime[index.unique.order],length)),
                           index.time = NA,
                           param  = NA,
                           n.cluster = NA,
                           index.cluster = NA)

    Upattern$index.cluster <- tapply(pattern.cluster$index.cluster, pattern.cluster$pattern, identity, simplify = FALSE)[Upattern$name]
    Upattern$index.time <- lapply(Upattern$index.cluster,function(iC){index.clusterTime[[iC[1]]]})
    Upattern$n.cluster <- table(pattern.cluster$var)[Upattern$var]

    if(NCOL(X.var)>0){
        ## find param with non-0 value in the design matrix
        Upattern$param <- lapply(Upattern$index.cluster, function(iC){
            names(which(colSums(abs(X.var[index.cluster[[iC[1]]],,drop=FALSE]))>0))
        })
        attr(Upattern, "level.var") <- name.pattern.var
    }else{
        Upattern$param <- vector(mode = "list", length = n.pattern.var)
    }
    rownames(structure$X$Upattern) <- NULL

    ## ** characterize each pattern
    if(NCOL(X.var)>0){
        Xpattern.var <- stats::setNames(lapply(Upattern$name,function(iP){ ## iP <- Upattern$name[1]
            iC.all <- pattern.cluster$index.cluster[which(pattern.cluster$pattern==iP)]
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
            attr(iX,"indicator.param") <- lapply(attr(iX,"Mindicator.param"),function(iM){which(iM>0)})
        
            return(iX)
        }),Upattern$var)
        
    }else{
        Xpattern.var <- stats::setNames(lapply(Upattern$name, function(iName){
            iM <- matrix(nrow = Upattern[Upattern$name==iName,"n.cluster"], ncol = 0)
            attr(iM,"index.cluster") <- Upattern[Upattern$name==iName,"index.cluster"]
            attr(iM,"index.time") <- Upattern[Upattern$name==iName,"index.time"]
            return(iM)
        }), Upattern$name)
    }

    ## ** export
    structure$X$Upattern <- Upattern
    structure$X$Xpattern.var <- Xpattern.var
    structure$X$pattern.cluster <- pattern.cluster
    return(structure)
}

## * .findUpatterns.IND
.findUpatterns.IND <- .findUpatterns.ID

## * .findUpatterns.CS
.findUpatterns.CS <- function(structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata,
                              sep = "."){

    ## ** identify unique var patterns
    Upatterns.init <- .findUpatterns.IND(structure = structure, 
                                         index.clusterTime = index.clusterTime, U.time = U.time,
                                         index.cluster = index.cluster, U.cluster = U.cluster,
                                         index.clusterStrata = index.clusterStrata, U.strata = U.strata,
                                         sep = sep)$X
    UpatternVar <- Upatterns.init$Upattern
    Xpattern.var <- Upatterns.init$Xpattern.var
    patternVar.cluster <- Upatterns.init$pattern.cluster

    ## ** identify unique cor patterns
    X.cor <- structure$X$cor
    structure.param <- structure$param[structure$param$type=="rho",]
    param.rho <- structure.param$name
    level.cor <- attr(structure.param,"level.cor")
    n.cluster <- length(U.cluster)

    corLP.cluster <- attr(structure$X$lp.cor,"cluster")
    
    ## *** group clusters by linear predictor patterns
    ## *** (to later work on a single representative cluster per pattern)
    ## vector summarizing the (combined) linear predictor of each cluster 
    patternCorLP.cluster <- sapply(corLP.cluster, function(iLp){ 
        paste0(as.numeric(iLp), collapse = sep)
    })
    ## list of cluster per (combined) linear predictor
    patternCorLP.index <- tapply(1:n.cluster,patternCorLP.cluster,identity)
    ## unique (combined) linear predictors
    patternCorLP.index1 <- sapply(patternCorLP.index,"[",1)
    UpatternCorLP <- patternCorLP.cluster[patternCorLP.index1]
    nObs.index1 <- sapply(index.cluster[patternCorLP.index1],length) ## number of observations per pattern

    ## *** form all pairs of observations within representative cluster
    ## *** and deduce correlation parameter
    ## identify all pairs
    pattern.pairwise <- do.call(rbind,lapply(patternCorLP.index1[nObs.index1>0], function(iC){ ## iC <- 1
        iPair <- base::t(.unorderedPairs(1:length(index.cluster[[iC]]), distinct = TRUE))
        data.frame(cluster = iC,
                   index.x = iPair[,1],
                   index.y = iPair[,2],
                   obs.x = index.cluster[[iC]][iPair[,1]],
                   obs.y = index.cluster[[iC]][iPair[,2]])
    }))
    pattern.pairwise$pattern <- patternCorLP.cluster[pattern.pairwise$cluster]
    rownames(pattern.pairwise) <- NULL

    ## convert design matrix into linear predictor
    pattern.pairwise$code.x <- as.character(interaction(as.data.frame(X.cor[pattern.pairwise$obs.x,,drop=FALSE]), sep = ""))
    pattern.pairwise$code.y <- as.character(interaction(as.data.frame(X.cor[pattern.pairwise$obs.y,,drop=FALSE]), sep = ""))
    ## convert linear predictor into correlation parameter
    indexX.lp <- attr(structure.param,"Xcode.xy")
    n.lp <- max(indexX.lp) ## number of distinct values for the linear predictor

    code2param <- do.call(rbind,apply(structure.param, 1, function(iDf){data.frame(code.x = iDf$code.x, code.y = iDf$code.y, name = iDf$name)}))
    M.param <- matrix(as.character(NA), nrow = n.lp, ncol = n.lp)
    M.param[code2param$code.x + (n.lp-1)*code2param$code.y] <- code2param$name
    
    pattern.pairwise$name <- M.param[indexX.lp[pattern.pairwise$code.x] + (n.lp-1)*indexX.lp[pattern.pairwise$code.y]]

    ## *** summarize correlation patterns
    UpatternCor <- data.frame(name = NA,
                              var = NA,
                              cor = UpatternCorLP,
                              index.strata = index.clusterStrata[patternCorLP.index1],
                              n.time = unname(sapply(index.clusterTime[patternCorLP.index1],length)),
                              index.time = NA,
                              param = NA,
                              index.cluster = NA,
                              n.cluster = unname(sapply(patternCorLP.index,length)))
    UpatternCor$index.cluster <- patternCorLP.index[UpatternCor$cor]
    UpatternCor$index.time <- lapply(UpatternCor$index.cluster,function(iC){index.clusterTime[[iC[1]]]})
    UpatternCor$param <- tapply(pattern.pairwise$name,pattern.pairwise$pattern,unique, simplify = FALSE)

    patternCor.cluster <- do.call(rbind,lapply(UpatternCorLP, function(iPattern){
        data.frame(index.cluster = patternCorLP.index[[iPattern]], cor = iPattern)
    }))
    patternCor.cluster <- patternCor.cluster[order(patternCor.cluster$index.cluster),,drop=FALSE]

    ## *** characterize each correlation pattern
    Xpattern.cor <- stats::setNames(lapply(UpatternCorLP, function(iPattern){ ## iPattern <- UpatternCorLP[1]

        iCluster <- patternCorLP.index1[[iPattern]]
        iIndex.cluster <- index.cluster[[iCluster]]
        iRep <- length(iIndex.cluster)        
        if(iRep==1){return(NULL)}
        iOut <- X.cor[iIndex.cluster,,drop=FALSE]
        iPattern.pairwise <- pattern.pairwise[pattern.pairwise$pattern == iPattern,]
        attr(iOut, "index.cluster") <- patternCorLP.index[[iPattern]]
        attr(iOut, "index.time") <- index.clusterTime[[iCluster]]
        attr(iOut, "index.strata") <- index.clusterStrata[[iCluster]]
        attr(iOut, "index.pair") <- rbind(data.frame(row =  iPattern.pairwise$index.x,
                                                     col = iPattern.pairwise$index.y,
                                                     param = iPattern.pairwise$name),
                                          data.frame(row =  iPattern.pairwise$index.y,
                                                     col = iPattern.pairwise$index.x,
                                                     param = iPattern.pairwise$name))
        attr(iOut, "index.vec2matrix") <- attr(iOut, "index.pair")$row + (attr(iOut, "index.pair")$col-1)*iRep
        attr(iOut, "param") <- unique(iPattern.pairwise$name)
        attr(iOut, "indicator.param") <- tapply(attr(iOut, "index.vec2matrix"),attr(iOut, "index.pair")$param,identity,simplify=FALSE)
        if(any(is.na(attr(iOut, "index.pair")$param))){
            attr(iOut, "indicator.param")[[as.character(NA)]] <- attr(iOut, "index.vec2matrix")[is.na(attr(iOut, "index.pair")$param)]
        }
        attr(iOut, "Mindicator.param") <- tapply(attr(iOut, "index.vec2matrix"),attr(iOut, "index.pair")$param,function(iIndex){
            iM <- matrix(0, nrow = iRep, ncol = iRep)
            iM[iIndex] <- 1
            return(iM)
        }, simplify=FALSE)

        return(iOut)
    }), UpatternCorLP)

    ## ** joint variance and correlation patterns
    pattern.name <- paste0(as.numeric(factor(patternVar.cluster$var,attr(UpatternVar,"level.var"))),
                           ":",
                           as.numeric(factor(patternCor.cluster$cor,UpatternCorLP)))
    pattern.cluster <- data.frame(index.cluster = 1:n.cluster,
                                  pattern = pattern.name,
                                  var = patternVar.cluster$var,
                                  cor = patternCor.cluster$cor)
    index.Upattern <- which(!duplicated(pattern.name))
    index.Upattern.order <- index.Upattern[order(pattern.cluster$pattern[index.Upattern])] ## re-order index so patterns appears in alphabetic order (instead of order driven by the dataset)

    Upattern <- data.frame(name = pattern.cluster$pattern[index.Upattern.order],
                           var = pattern.cluster$var[index.Upattern.order],
                           cor = pattern.cluster$cor[index.Upattern.order],
                           index.strata = NA,
                           n.time = NA,
                           index.time = NA,
                           param = NA,
                           index.cluster = NA,
                           n.cluster = NA)
    Upattern$index.cluster <- tapply(pattern.cluster$index.cluster,pattern.cluster$pattern,identity,simplify=FALSE)
    Upattern$n.cluster <- unname(sapply(Upattern$index.cluster,length))
    Upattern.index.cluster1 <- unname(sapply(Upattern$index.cluster,"[",1))
    Upattern$index.strata <- unname(index.clusterStrata[Upattern.index.cluster1])
    Upattern$index.time <- unname(index.clusterTime[Upattern.index.cluster1])
    Upattern$n.time <- unname(sapply(Upattern$index.time,length))
    Upattern$param <- unname(mapply(x = UpatternVar$param[match(Upattern$var,UpatternVar$var)], y = UpatternCor$param[match(Upattern$cor,UpatternCor$cor)], FUN = function(x,y){unname(c(x,y))}, SIMPLIFY = FALSE))
    attr(Upattern,"level.var") <- attr(UpatternVar,"level.var")
    attr(Upattern,"level.cor") <- UpatternCorLP

    ## ** export
    structure$X$Upattern <- Upattern
    structure$X$Xpattern.var <- Xpattern.var
    structure$X$Xpattern.cor <- Xpattern.cor
    structure$X$pattern.cluster <- pattern.cluster
    return(structure)
}

## * .findUpatterns.RE
.findUpatterns.RE <- .findUpatterns.CS

## * .findUpatterns.TOEPLITZ
.findUpatterns.TOEPLITZ <- .findUpatterns.CS

## * .findUpatterns.UN
.findUpatterns.UN <- .findUpatterns.CS

## * .findUpatterns.EXP
.findUpatterns.EXP <- .findUpatterns.CS

## * .findUpatterns_CUSTOM
.findUpatterns.CUSTOM <- function(structure, 
                                  index.clusterTime, U.time,
                                  index.cluster, U.cluster,
                                  index.clusterStrata, U.strata,
                                  sep = c(":","::X::")){

    X.var <- structure$X$var
    X.cor <- structure$X$cor
    param <- structure$param

    ## ** identify unique patterns

    ## *** var
    lpObs.var <- interaction(as.data.frame(X.var), drop = TRUE, sep = sep[1])
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
        lpObs.cor <- interaction(as.data.frame(X.cor), drop = TRUE, sep = sep[1])
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
                            sep = sep[1])       
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

        ls.pair[iN.pair[iN.pair>0]] <- lapply(iN.pair[iN.pair>0], function(iN){.unorderedPairs(1:iN, distinct = TRUE)})
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
    structure$X$Upattern <- Upattern
    structure$X$Xpattern.var <- Xpattern.var
    structure$X$Xpattern.cor <- Xpattern.cor
    structure$X$pattern.cluster <- pattern.cluster
    
    return(structure)
}
##----------------------------------------------------------------------
### findPatterns.R ends here
