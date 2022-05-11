### findPatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:06) 
## Version: 
## Last-Updated: maj 11 2022 (15:34) 
##           By: Brice Ozenne
##     Update #: 108
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .findUpatterns
.findUpatterns <- function(X.var, X.cor, table.rho,
                           Upattern = NULL, data, heterogeneous,
                           time.var, index.clusterTime, order.clusterTime, U.time,
                           index.cluster, U.cluster,
                           strata.var, index.clusterStrata, strata.param, U.strata,
                           sep = "_X_XX_X_"){


    out <- list()
    
    ## ** identify unique patterns

    ## *** var
    lpObs.var <- interaction(as.data.frame(X.var), drop = TRUE, sep = sep)
    Ulp.var <- sort(levels(lpObs.var))
    lpnObs.var <- as.numeric(factor(lpObs.var, levels = Ulp.var)) 
    lpnCluster.var <- stats::setNames(sapply(U.cluster, function(iC){
        paste(lpnObs.var[index.cluster[[iC]][order.clusterTime[[iC]]]], collapse=".")
    }), U.cluster)
    attr(lpnCluster.var, "levels") <- Ulp.var

    name.pattern.var <- sort(unique(lpnCluster.var))
    cluster.pattern.var <- which(!duplicated(lpnCluster.var))
    cluster.pattern.var <- cluster.pattern.var[match(lpnCluster.var[cluster.pattern.var], name.pattern.var)]
    n.pattern.var <- length(name.pattern.var)

    ## *** cor
    if(!is.null(X.cor)){
        lpObs.cor <- interaction(as.data.frame(X.cor), drop = TRUE, sep = sep)
        Ulp.cor <- sort(levels(lpObs.cor)) 
        lpnObs.cor <- as.numeric(factor(lpObs.cor, levels = Ulp.cor)) 
        lpnCluster.cor0 <- stats::setNames(lapply(U.cluster, function(iC){
            lpnObs.cor[index.cluster[[iC]][order.clusterTime[[iC]]]]
        }), U.cluster)
        lpnCluster.cor <- sapply(lpnCluster.cor0, paste, collapse = ".")
        attr(lpnCluster.cor, "levels") <- Ulp.cor
        if(!identical(attr(table.rho,"levels"), Ulp.cor)){
            warning("Something went wrong when extracting the correlation patterns. \n")
        }
        
        name.pattern.cor <- sort(unique(lpnCluster.cor))
        cluster.pattern.cor <- which(!duplicated(lpnCluster.cor))
        cluster.pattern.cor <- cluster.pattern.cor[match(lpnCluster.cor[cluster.pattern.cor], name.pattern.cor)]
        n.pattern.cor <- length(name.pattern.cor)
    }

    ## *** var and cor
    if(!is.null(X.cor)){
        lpnCluster <- paste(as.numeric(factor(lpnCluster.var, name.pattern.var)),
                            as.numeric(factor(lpnCluster.cor, name.pattern.cor)),
                            sep = ".")       
    }else{
        lpnCluster <- as.character(as.numeric(factor(lpnCluster.var, name.pattern.var)))
    }
    test.duplicated <- duplicated(lpnCluster)
    name.pattern <- lpnCluster[!test.duplicated]
    order.pattern <- order(name.pattern)
    name.pattern <- name.pattern[order.pattern]
    n.pattern <- length(name.pattern)

    ## *** individual belonging to each pattern
    out$pattern.cluster <- data.frame(cluster = U.cluster,
                                      index.cluster = as.numeric(as.factor(U.cluster)),
                                      pattern = factor(name.pattern[match(lpnCluster, name.pattern)], levels = name.pattern))
browser()
    ## *** var/cor pattern corresponding to pattern
    pattern.indexVar <- lpnCluster.var[!test.duplicated][order.pattern]
    if(!is.null(X.cor)){
        pattern.indexCor <- lpnCluster.cor[!test.duplicated][order.pattern]
    }else{
        pattern.indexCor <- NULL
    }
    out$Upattern <- data.frame(name = factor(name.pattern, name.pattern),
                               var = factor(pattern.indexVar, name.pattern.cor),
                               cor = factor(pattern.indexCor, name.pattern.var),
                               strata = NA,
                               n.time = NA,
                               param = NA,
                               cluster = NA,
                               n.cluster = NA)
    out$Upattern$cluster <- stats::setNames(lapply(name.pattern, function(iPattern){which(lpnCluster == iPattern)}),
                                            name.pattern)
    out$Upattern$n.cluster <- sapply(out$Upattern$cluster, length)
    out$Upattern$n.time <- lapply(out$Upattern$cluster, function(iId){length(index.clusterTime[[iId[1]]])})
    out$Upattern$strata <- lapply(out$Upattern$cluster, function(iId){index.clusterStrata[iId[1]]})
    attr(out$Upattern, "lp.UpatternVar") <- Ulp.var
    attr(out$Upattern, "lp.UpatternCor") <- Ulp.cor

    out$pattern.cluster$var <- stats::setNames(stats::setNames(out$Upattern$var, out$Upattern$name)[out$pattern.cluster$pattern], out$pattern.cluster$cluster)
    if(!is.null(X.cor)){
        out$pattern.cluster$cor <- stats::setNames(stats::setNames(out$Upattern$cor, out$Upattern$name)[out$pattern.cluster$pattern], out$pattern.cluster$cluster)
    }
    ## ** characterize each pattern

    ## *** var
    out$Xpattern.var <- stats::setNames(lapply(name.pattern.var,function(iP){ ## iP <- name.pattern.var[1]
        iC.all <- out$pattern.cluster$cluster[which(out$pattern.cluster$var==iP)]
        iC <- iC.all[1]
        iX <- X.var[index.cluster[[iC]][order.clusterTime[[iC]]],,drop=FALSE]
        iIndex.param <- which(colSums(iX)!=0)
        attr(iX,"cluster") <- iC.all
        attr(iX,"strata") <- index.clusterStrata[iC]
        attr(iX,"cluster") <- index.clusterStrata[iC]
        attr(iX,"param") <- names(iIndex.param)
        attr(iX,"indicator.param") <- stats::setNames(lapply(iIndex.param,function(iCol){which(tcrossprod(iX[,iCol],rep(1,NROW(iX))) + t(tcrossprod(iX[,iCol],rep(1,NROW(iX)))) > 0)}),
                                                      names(iIndex.param))
        return(iX)
    }),name.pattern.var)

    out$Upattern$param <- stats::setNames(lapply(out$Upattern$name, function(iP){ ## iP <- out$Upattern$name[2]
        iP.var <- out$Upattern[["var"]][[which(out$Upattern$name==iP)]]
        return(attr(out$Xpattern.var[[iP.var]],"param"))
    }), out$Upattern$name)

    ## *** cor
    if(!is.null(X.cor)){
        lp.xy <- paste0(table.rho[,"lp.x"],".",table.rho[,"lp.y"])
        UlpnCluster.cor0 <- lpnCluster.cor0[cluster.pattern.cor]
        
        iN.pair <- unique(sapply(UlpnCluster.cor0, length))
        ls.pair <- vector(mode = "list", length = max(iN.pair))
        ls.pair[iN.pair] <- lapply(iN.pair, function(iN){.unorderedPairs(1:iN, distinct = TRUE)})

        out$Xpattern.cor <- stats::setNames(lapply(name.pattern.cor,function(iP){ ## iP <- name.pattern.cor[1]
            iC.all <- out$pattern.cluster$cluster[which(out$pattern.cluster$cor==iP)]
            iC <- iC.all[1]
            iX <- X.cor[index.cluster[[iC]][order.clusterTime[[iC]]],,drop=FALSE]

            iC.code <- lpnCluster.cor0[[iC]]
            iC.pair <- ls.pair[[length(iC.code)]]
            iC.codepair <- matrix(iC.code[iC.pair], ncol = 2, byrow = FALSE, dimnames = list(NULL,c("lp.x","lp.y")))
            iC.param <- table.rho$param[match(paste0(iC.codepair[,"lp.x"],".",iC.codepair[,"lp.y"]),lp.xy)]

            iC.table <- rbind(data.frame(row = iC.pair[1,], col = iC.pair[2,], param = iC.param),
                              data.frame(row = iC.pair[2,], col = iC.pair[1,], param = iC.param))
            attr(iX, "cluster") <- iC.all
            attr(iX, "strata") <- index.clusterStrata[iC]
            attr(iX, "index.pair") <- iC.table
            attr(iX, "index.vec2matrix") <- c(iC.table[,"row"] + NROW(iX) * (iC.table[,"col"] - 1))
            attr(iX, "param") <- sort(unique(iC.table$param))
            attr(iX, "indicator.param") <- stats::setNames(lapply(attr(iX, "param"),function(iP){attr(iX, "index.vec2matrix")[which(iC.table$param==iP)]}),
                                                           attr(iX, "param"))
            return(iX)
        }),name.pattern.var)

        out$Upattern$param <- stats::setNames(lapply(out$Upattern$name, function(iP){ ## iP <- out$Upattern$name[1]
            c(out$Upattern[["param"]][[which(out$Upattern$name==iP)]],
              attr(out$Xpattern.cor[[out$Upattern[["cor"]][[which(out$Upattern$name==iP)]]]],"param"))
        }), out$Upattern$name)
    }

    ## ** pair of variance-covariance coefficients
    out$pair.varcoef <- stats::setNames(lapply(out$Upattern$name, function(iPattern){## ## iPattern <- structure$X$Upattern$name[1]
        iParam <- out$Upattern$param[[iPattern]]

        iOut <- .unorderedPairs(iParam)
        attr(iOut, "key") <- matrix(NA, nrow = length(iParam), ncol = length(iParam), dimnames = list(iParam,iParam))
        for(iCol in 1:NCOL(iOut)){
            attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
            attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol]] <- iCol
        }
        return(iOut)
    }),out$Upattern$name)
browser()
    ## ** export
    return(out)
}



##----------------------------------------------------------------------
### findPatterns.R ends here
