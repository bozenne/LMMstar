### findPatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:06) 
## Version: 
## Last-Updated: apr 13 2022 (10:31) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .findUpatterns
.findUpatterns <- function(X.var, X.cor, Upattern = NULL, lpdiff.rho, data, heterogeneous,
                           time.var, index.clusterTime, order.clusterTime, U.time,
                           index.cluster, U.cluster,
                           strata.var, strata.param, U.strata){

    ## ** re-order by strata
    strata.param <- sort(strata.param)
    name.param <- names(strata.param)
    time.param <- stats::setNames(vector(mode="list", length = length(strata.param)), name.param)

    ## ** find all variance patterns
    ## variance pattern for all clusters
    ls.patternVar.cluster <- .buildPatterns(X = X.var, data = data, Upattern = Upattern,
                                            U.cluster = U.cluster, index.cluster = index.cluster, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime,
                                            U.strata = U.strata, strata.var = strata.var)
    
    ## find unique patterns and order them
    out.orderUPatternVar <- .orderUpattern(ls.pattern.cluster = ls.patternVar.cluster,
                                           time.cluster = index.clusterTime,
                                           order.cluster = order.clusterTime,
                                           U.time = U.time) 
    UpatternVar.cluster <- out.orderUPatternVar$name
    ## update vector giving the index of the pattern corresponding to each cluster
    patternVar.cluster <- as.numeric(factor(unlist(ls.patternVar.cluster), levels = UpatternVar.cluster))
    ## generate design matrix associated to each pattern
    X.var2 <- stats::setNames(lapply(out.orderUPatternVar$cluster,function(iIndex){ ## iIndex <- out.orderUPatternVar$cluster[1]
        iX <- X.var[index.cluster[[iIndex]][order.clusterTime[[iIndex]]],,drop=FALSE]
        iIndex.param <- which(colSums(iX)!=0)
        attr(iX,"indicator.param") <- stats::setNames(lapply(iIndex.param,function(iCol){which(tcrossprod(iX[,iCol],rep(1,NROW(iX))) + t(tcrossprod(iX[,iCol],rep(1,NROW(iX)))) > 0)}),
                                                      names(iIndex.param))
        return(iX)
    }),UpatternVar.cluster)
    attr(X.var2,"assign") <- attr(X.var,"assign")

    ## first timepoint where each variance parameter appear
    for(iVar in colnames(X.var)){
        time.param[[iVar]] <- sort(match(unique(as.character(data[X.var[,iVar]==1,time.var])),U.time))[1]
    }

    ## ** find all correlation patterns
    if(is.null(X.cor)){
        UpatternCor.cluster <- "1"
        patternCor.cluster <- stats::setNames(rep(1, length = length(U.cluster)), U.cluster)
        X.cor2 <- NULL
    }else{
        ## correlation pattern for all clusters
        ls.patternCor.cluster <- .buildPatterns(X = X.cor, data = data, Upattern = Upattern,
                                                U.cluster = U.cluster, index.cluster = index.cluster, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime,
                                                U.strata = U.strata, strata.var = strata.var)
        ## find unique patterns and order them
        out.orderUPatternCor <- .orderUpattern(ls.pattern.cluster = ls.patternCor.cluster, time.cluster = index.clusterTime, order.cluster = order.clusterTime, U.time = U.time) 
        UpatternCor.cluster <- out.orderUPatternCor$name
        ## update vector giving the index of the pattern corresponding to each cluster
        patternCor.cluster <- as.numeric(factor(unlist(ls.patternCor.cluster), levels = UpatternCor.cluster))
        ## generate design matrix associated to each pattern
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
                    if(heterogeneous){
                        if(all(iX1==iX2)){return(cbind("R",iX1))}else{return(cbind("D",iX2-iX1))}
                    }else{
                        if(all(iX1==iX2)){return(matrix(c("R",""),nrow=1,ncol=2))}else{return(cbind("D",sum(iX2!=iX1)))}
                    }

                }))
            ), drop = TRUE))
            ## generate design matrix relative to the coefficient
            if(length(lpdiff.rho)==1 && all(iDiff==lpdiff.rho)){
                iZ <- matrix(1, nrow = length(iDiff), ncol = 1, dimnames = list(NULL,names(lpdiff.rho)))
            }else{
                iZ <- stats::model.matrix(~0+X,data.frame(X = factor(names(lpdiff.rho)[match(iDiff,lpdiff.rho)], levels = names(lpdiff.rho)),stringsAsFactors = FALSE))
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
        ## first pair of timepoints where each correlation parameter appear
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

    ## **  assemble variance and correlation patterns to identify unique patterns
    pattern.cluster <- paste(patternVar.cluster,patternCor.cluster,sep=".")
    index.Nduplicated <- which(duplicated(pattern.cluster)==FALSE)
    index.Nduplicated2 <- index.Nduplicated[order(patternVar.cluster[index.Nduplicated],patternCor.cluster[index.Nduplicated])]
    Upattern <- data.frame(name = pattern.cluster[index.Nduplicated2],
                           var = patternVar.cluster[index.Nduplicated2],
                           cor = patternCor.cluster[index.Nduplicated2],
                           time = NA,
                           strata = NA,
                           ncluster = NA,
                           param = NA,
                           example = index.Nduplicated2,
                           stringsAsFactors = FALSE)
    attr(Upattern,"name.UpatternVar") <- UpatternVar.cluster
    attr(Upattern,"levels.UpatternVar") <- attr(ls.patternVar.cluster,"levels")
    attr(X.var2,"M.level") <- attr(X.var,"M.level")
    attr(X.var2,"original.colnames") <- stats::setNames(attr(X.var,"original.colnames"), colnames(X.var))
    if(!is.null(X.cor)){
        attr(Upattern,"name.UpatternCor") <- UpatternCor.cluster
        attr(Upattern,"levels.UpatternCor") <- attr(ls.patternCor.cluster,"levels")
        attr(X.cor2,"M.level") <- attr(X.cor,"M.level")
        attr(X.cor2,"original.colnames") <- attr(X.cor,"original.colnames")
    }
    
    ## ** strata and time associated to each unique pattern
    Upattern$time <- lapply(ls.patternVar.cluster[index.Nduplicated2],attr,"index.alltime")
    Upattern$strata <- sapply(ls.patternVar.cluster[index.Nduplicated2],attr,"index.strata")
    Upattern$ncluster <- sapply(Upattern$name,function(iN){sum(iN==pattern.cluster)})
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

    ## ** pair of variance-covariance coefficients
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

    ## ** export
    out <- list(pattern.cluster = stats::setNames(pattern.cluster,U.cluster),
                cluster.pattern = cluster.pattern,
                Upattern = Upattern,
                X.cor.pairwise = X.cor2,
                X.var = X.var2,
                time.param = time.param,
                pair.varcoef = pair.varcoef)
    return(out)
}


## * .buildPatterns
## convert design matrix into character string with associated strata and time
.buildPatterns <- function(X, data, Upattern,
                           U.cluster, index.cluster, index.clusterTime, order.clusterTime,
                           U.strata, strata.var){
    
    ## *** combine all covariate values at the observation levels and convert to positive integer
    if(NROW(X)==1){
        level.var <- interaction(as.list(X),drop=TRUE)
    }else{
        level.var <- interaction(as.data.frame(X),drop=TRUE)
    }

    Ulevel.var <- levels(level.var)
    numlevel.var <- as.numeric(level.var)

    ## *** structure the covariate pattern by cluster and add strata and time indexes
    patternVar.cluster <- stats::setNames(lapply(U.cluster, function(iC){ ## iC <- U.cluster[[1]]
        iIndex.obs <- index.cluster[[iC]][order.clusterTime[[iC]]]
        iLevel <- paste(numlevel.var[iIndex.obs], collapse=".")
        attr(iLevel,"index.level") <- numlevel.var[iIndex.obs]
        attr(iLevel,"index.alltime") <- index.clusterTime[[iC]][order.clusterTime[[iC]]]
        if(!is.na(strata.var)){
            attr(iLevel,"index.strata") <- which(unique(data[index.cluster[[iC]],strata.var])==U.strata)
        }else{
            attr(iLevel,"index.strata") <- 1
        }
        return(iLevel)
    }), U.cluster)
    attr(patternVar.cluster,"levels") <- Ulevel.var

    ## *** export
    return(patternVar.cluster)
}

## * .orderUpattern
## .orderMatTime(time.cluster = list(1:5,1,1:3,2:4), order.cluster = NULL, U.time = 1:5)
.orderUpattern <- function(ls.pattern.cluster, time.cluster, order.cluster, U.time){

    pattern.cluster <- unlist(ls.pattern.cluster)
    indexU.cluster <- which(!duplicated(pattern.cluster))
    score.pattern <- pattern.cluster[indexU.cluster]
    if(length(indexU.cluster)==1){
        out <- list(cluster = indexU.cluster, name = pattern.cluster[indexU.cluster])
    }
    
    ## score the time vector corresponding to each cluster
    n.time <- length(U.time)
    n.Ucluster <- length(indexU.cluster)

    ls.score.time <- lapply(indexU.cluster, function(iC){ ## iC <- 1
        if(!is.null(order.cluster)){
            iTime <- time.cluster[[iC]][order.cluster[[iC]]]
        }else{
            iTime <- time.cluster[[iC]]
        }
            
        if(n.time == length(iTime)){
            iOut <- c(-1,0) ## no missing, no missing
        }else{
            ## n-position of the first missing value, number of missing values
            iOut <- c("firstmis" = n.time-min(setdiff(1:n.time, iTime)),"nmis"=n.time-length(iTime)) 
        }
        iLevel <- matrix(NA, nrow = 1, ncol = length(U.time), dimnames = list(NULL, U.time))
        iLevel[1,iTime] <-  attr(ls.pattern.cluster[[iC]],"index.level")
        attr(iOut,"level") <- iLevel
        return(iOut)
    })
    score.time <- do.call(rbind, ls.score.time)
    score.pattern2 <- do.call(base::order, args = c(as.list(as.data.frame(do.call(rbind,lapply(ls.score.time,attr, "level"))))))
    ## better than just score.pattern because it properly handle numeric i.e. 9<10 instead of "9" > "10"
    out <- list(cluster = indexU.cluster[order(score.time[,1],score.time[,2],score.pattern2)])
    out$name <- pattern.cluster[out$cluster]
    return(out)
}

##----------------------------------------------------------------------
### findPatterns.R ends here
