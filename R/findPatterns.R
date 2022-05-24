### findPatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:06) 
## Version: 
## Last-Updated: maj 24 2022 (19:20) 
##           By: Brice Ozenne
##     Update #: 241
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .findUpatterns
.findUpatterns <- function(X.var, X.cor, param,
                           Upattern = NULL, data, heterogeneous,
                           time.var, index.clusterTime, U.time,
                           index.cluster, U.cluster,
                           strata.var, index.clusterStrata, U.strata,
                           sep = c(":","::X::")){

    out <- list()
    Ulp.var <- attr(param,"level.var")
    Ulp.cor <- attr(param,"level.cor")

    ## ** identify unique patterns

    ## *** var
    lpObs.var <- interaction(as.data.frame(X.var), drop = TRUE, sep = sep[1])
    if(!identical(sort(levels(lpObs.var)), sort(Ulp.var))){
        warning("Something went wrong when extracting the variance patterns. \n")
    }
    lpnObs.var <- as.numeric(factor(lpObs.var, levels = Ulp.var)) 
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
        if(!identical(sort(levels(lpObs.cor)), sort(Ulp.cor))){
            warning("Something went wrong when extracting the correlation patterns. \n")
        }
        lpnObs.cor <- as.numeric(factor(lpObs.cor, levels = Ulp.cor)) 
        lpnCluster.cor0 <- stats::setNames(lapply(U.cluster, function(iC){
            if(length(index.cluster[[iC]])>1){
                return(lpnObs.cor[index.cluster[[iC]]])
            }else{
                return(NULL)
            }
        }), U.cluster)
        lpnCluster.cor <- sapply(lpnCluster.cor0, paste, collapse = ".")
        name.pattern.cor <- unique(lpnCluster.cor)
        if(!heterogeneous){
            name.pattern.cor <- name.pattern.cor[order(nchar(as.character(name.pattern.cor)))]
        }else{
            name.pattern.cor <- sort(name.pattern.cor)
        }
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

    out$pattern.cluster <- data.frame(index.cluster = U.cluster,
                                      pattern = factor(name.pattern[match(lpnCluster, name.pattern)], levels = name.pattern))

    ## *** var/cor pattern corresponding to pattern
    pattern.indexVar <- lpnCluster.var[!test.duplicated][order.pattern]
    if(!is.null(X.cor)){
        pattern.indexCor <- lpnCluster.cor[!test.duplicated][order.pattern]
    }else{
        pattern.indexCor <- NA
    }
    out$Upattern <- data.frame(name = name.pattern,
                               var = pattern.indexVar,
                               cor = pattern.indexCor,
                               index.strata = NA,
                               n.time = NA,
                               param = NA,
                               index.cluster = NA,
                               n.cluster = NA)

    out$Upattern$index.cluster <- stats::setNames(lapply(name.pattern, function(iPattern){which(lpnCluster == iPattern)}), name.pattern)
    out$Upattern$n.cluster <- sapply(out$Upattern$index.cluster, length)
    out$Upattern$n.time <- sapply(out$Upattern$index.cluster, function(iId){length(index.clusterTime[[iId[1]]])})
    out$Upattern$index.strata <- sapply(out$Upattern$index.cluster, function(iId){unname(unique(index.clusterStrata[iId]))}, simplify = FALSE)
    attr(out$Upattern, "level.var") <- name.pattern.var
    attr(out$Upattern, "level.cor") <- name.pattern.cor
    attr(out$Upattern, "lp.UpatternVar") <- Ulp.var
    out$pattern.cluster$var <- stats::setNames(stats::setNames(out$Upattern$var, out$Upattern$name)[out$pattern.cluster$pattern], out$pattern.cluster$cluster)
    if(!is.null(X.cor)){
        attr(out$Upattern, "lp.UpatternCor") <- Ulp.cor
        out$pattern.cluster$cor <- stats::setNames(stats::setNames(out$Upattern$cor, out$Upattern$name)[out$pattern.cluster$pattern], out$pattern.cluster$cluster)
    }

    ## ** characterize each pattern

    ## *** var
    param.sigma <- param[param$type=="sigma","name"]
    param.k <- param[param$type=="k","name"]
    out$Xpattern.var <- stats::setNames(lapply(name.pattern.var,function(iP){ ## iP <- name.pattern.var[1]
        iC.all <- out$pattern.cluster$index.cluster[which(out$pattern.cluster$var==iP)]
        iC <- iC.all[1]
        iX <- X.var[index.cluster[[iC]],,drop=FALSE]
        iIndex.param <- which(colSums(iX)!=0)
        iParam.k <- intersect(names(iIndex.param),param.k)
        iParam.sigma <- intersect(names(iIndex.param),param.sigma)

        attr(iX,"index.cluster") <- iC.all
        attr(iX,"index.strata") <- unname(index.clusterStrata[iC])
        attr(iX,"param") <- names(iIndex.param)
        attr(iX,"indicator.param") <- stats::setNames(lapply(iIndex.param,function(iCol){which(tcrossprod(iX[,iCol],rep(1,NROW(iX))) + t(tcrossprod(iX[,iCol],rep(1,NROW(iX)))) > 0)}),
                                                      names(iIndex.param))
        attr(iX,"Mindicator.param") <- c(stats::setNames(lapply(iParam.sigma, function(iParam){
            return(matrix(2, nrow = NROW(iX), ncol = NROW(iX)))
        }),iParam.sigma),
        stats::setNames(lapply(iParam.k, function(iParam){
            iM <- matrix(0, nrow = NROW(iX), ncol = NROW(iX))
            iM[attr(iX,"indicator.param")[[iParam]]] <- 1
            diag(iM) <- 2*diag(iM)
            return(iM)
        }),iParam.k))

        return(iX)
    }),name.pattern.var)

    out$Upattern$param <- stats::setNames(lapply(out$Upattern$name, function(iP){ ## iP <- out$Upattern$name[2]
        iP.var <- out$Upattern[["var"]][[which(out$Upattern$name==iP)]]
        return(attr(out$Xpattern.var[[iP.var]],"param"))
    }), out$Upattern$name)

    ## *** cor
    if(!is.null(X.cor)){
        param.rho <- param[param$type=="rho",,drop=FALSE]
        ls.lp.xy <- mapply(x = param.rho$code.x, y = param.rho$code.y, FUN = paste, sep = "::X::", SIMPLIFY = FALSE)
        lp.xy <- unlist(lapply(names(ls.lp.xy), function(iN){stats::setNames(ls.lp.xy[[iN]],rep(iN, length(ls.lp.xy[[iN]])))}))
        UlpnCluster.cor0 <- lpnCluster.cor0[cluster.pattern.cor]
        
        iN.pair <- unique(sapply(UlpnCluster.cor0, length))
        ls.pair <- vector(mode = "list", length = max(iN.pair))

        ls.pair[iN.pair[iN.pair>0]] <- lapply(iN.pair[iN.pair>0], function(iN){.unorderedPairs(1:iN, distinct = TRUE)})
        out$Xpattern.cor <- stats::setNames(lapply(setdiff(name.pattern.cor,""),function(iP){ ## iP <- name.pattern.cor[1]
            ## ignore pattern "", i.e. patterns with no pair
            iC.all <- out$pattern.cluster$index.cluster[which(out$pattern.cluster$cor==iP)]
            iC <- iC.all[1]
            iX <- X.cor[index.cluster[[iC]],,drop=FALSE]
            iC.code <- lpnCluster.cor0[[iC]]
            iC.pair <- ls.pair[[length(iC.code)]]
            iC.codepair <- matrix(iC.code[iC.pair], ncol = 2, byrow = TRUE, dimnames = list(NULL,c("lp.x","lp.y")))
            iC.lp.xy <- paste0(attr(param,"level.cor")[iC.codepair[,"lp.x"]],sep[2],attr(param,"level.cor")[iC.codepair[,"lp.y"]])
            iC.param <- names(lp.xy)[match(iC.lp.xy,lp.xy)]
            ## unique(iC.lp.xy)

            iC.table <- rbind(data.frame(row = iC.pair[1,], col = iC.pair[2,], param = iC.param),
                              data.frame(row = iC.pair[2,], col = iC.pair[1,], param = iC.param))
            ## iX.pairwise <- matrix(0, nrow = NROW(iC.table), ncol = length(param.rho$name),
            ##                       dimnames = list(paste0("(",iC.table$row,",",iC.table$col,")"), param.rho$name))
            ## iX.pairwise[1:NROW(iC.table) + (match(iC.table$param , param.rho$name)-1) * NROW(iC.table)] <- 1
            attr(iX, "index.cluster") <- iC.all
            attr(iX, "index.strata") <- unname(index.clusterStrata[iC])
            attr(iX, "index.pair") <- iC.table
            attr(iX, "index.vec2matrix") <- c(iC.table[,"row"] + NROW(iX) * (iC.table[,"col"] - 1))
            attr(iX, "param") <- sort(unique(iC.table$param))
            attr(iX, "indicator.param") <- stats::setNames(lapply(attr(iX, "param"),function(iP){attr(iX, "index.vec2matrix")[which(iC.table$param==iP)]}),
                                                           attr(iX, "param"))
            attr(iX,"Mindicator.param") <- stats::setNames(lapply(attr(iX, "param"), function(iParam){
                iM <- matrix(0, nrow = NROW(iX), ncol = NROW(iX))
                iM[attr(iX,"indicator.param")[[iParam]]] <- 1
                return(iM)
            }),attr(iX, "param"))

            return(iX)
            ## iM <- matrix(0,length(iC.code),length(iC.code)); iM[attr(iX.pairwise, "index.vec2matrix")] <- iX.pairwise %*% 1:NROW(param.rho)
        }),setdiff(name.pattern.cor,""))

        out$Upattern$param <- stats::setNames(lapply(out$Upattern$name, function(iP){ ## iP <- out$Upattern$name[1]
            iPattern.cor <- out$Upattern$cor[which(out$Upattern$name==iP)]
            iParam <- out$Upattern$param[[which(out$Upattern$name==iP)]]
            if(iPattern.cor==""){
                ## pattern "", i.e. patterns with no pair
                return(iParam)
            }else{
                return(c(iParam, attr(out$Xpattern.cor[[iPattern.cor]],"param")))
            }
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

    ## ** export
    return(out)
}



##----------------------------------------------------------------------
### findPatterns.R ends here
