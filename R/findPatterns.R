### findPatterns.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:06) 
## Version: 
## Last-Updated: maj 12 2023 (11:35) 
##           By: Brice Ozenne
##     Update #: 482
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
                              sep = c(".","::X::")){


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
            paste(lpnObs.var[index.cluster[[iC]]], collapse=".")
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
    index.unique <- which(!duplicated(pattern.cluster$pattern))
    n.pattern.var <- length(name.pattern.var)

    Upattern <- data.frame(name = as.character(pattern.cluster$pattern[index.unique]),
                           var = as.character(pattern.cluster$var[index.unique]),
                           cor = NA,
                           index.strata = unname(index.clusterStrata[index.unique]),
                           n.time = unname(sapply(index.clusterTime[index.unique],length)),
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
        }),name.pattern.var)
    
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
                              sep = c(".","::X::")){

    browser()
    ## ** extract from object
    X.var <- structure$X$var
    X.cor <- structure$X$cor
    param <- structure$param
    Ulp.var <- attr(param,"level.var")
    Ulp.cor <- attr(param,"level.cor")

    heterogeneous <- structure$heterogeneous

    ## ** identify unique patterns

    ## *** var
    if(NCOL(X.var)>0){
        lpObs.var <- interaction(as.data.frame(X.var), drop = TRUE, sep = sep[1])
        if(!all(levels(lpObs.var) %in% Ulp.var)){
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
    }

    ## *** cor
    if(!is.null(X.cor)){
        lpObs.cor <- interaction(as.data.frame(X.cor), drop = TRUE, sep = sep[1])

        if(!is.null(Ulp.cor)){
            if(!all(levels(lpObs.cor) %in% sort(Ulp.cor))){
                warning("Something went wrong when extracting the correlation patterns. \n")
            }
            lpnObs.cor <- as.numeric(factor(lpObs.cor, levels = Ulp.cor))
        }else{
            lpnObs.cor <- as.numeric(lpObs.cor)
        }

        lpnCluster.cor0 <- stats::setNames(lapply(U.cluster, function(iC){
            if(length(index.cluster[[iC]])>1){
                return(lpnObs.cor[index.cluster[[iC]]])
            }else{
                return(NULL)
            }
        }), U.cluster)
        lpnCluster.cor <- sapply(lpnCluster.cor0, paste, collapse = ".")
        name.pattern.cor <- unique(lpnCluster.cor)
        if(!is.character(heterogeneous) && !heterogeneous){
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
        if(NCOL(X.var)>0){
            lpnCluster <- paste(as.numeric(factor(lpnCluster.var, name.pattern.var)),
                                as.numeric(factor(lpnCluster.cor, name.pattern.cor)),
                                sep = sep[1])
        }else{
            lpnCluster <- as.numeric(factor(lpnCluster.cor, name.pattern.cor))
            name.pattern.var <- NULL
        }
    }else if(NCOL(X.var)>0){
        lpnCluster <- as.character(as.numeric(factor(lpnCluster.var, name.pattern.var)))
    }else{ ## no variance nor correlation structure
        pattern.cluster <- data.frame(index.cluster = names(index.cluster),
                                      pattern = sapply(index.cluster,length),
                                      var = paste0("length",sapply(index.cluster,length)))
        Upattern.nRep <- unique(pattern.cluster$pattern)        
        Upattern.name <- paste0("length",Upattern.nRep)

        structure$X$Upattern <- data.frame(name = 1:length(Upattern.nRep),
                                           var = Upattern.name,
                                           cor = NA,
                                           index.strata = 1,
                                           n.time = Upattern.nRep,
                                           time = NA,
                                           param  = NA,
                                           index.cluster = NA,
                                           n.cluster = table(pattern.cluster$var)[Upattern.name])
        structure$X$Upattern$param <- list(NULL)
        structure$X$Upattern$index.cluster <- split(pattern.cluster$index.cluster,pattern.cluster$var)[Upattern.name]
        structure$X$Upattern$time <- lapply(structure$X$Upattern$index.cluster,function(iC){index.clusterTime[[iC[1]]]})
        rownames(structure$X$Upattern) <- NULL

        structure$X$Xpattern.var <- stats::setNames(lapply(Upattern.name, function(iPattern){ ## iPattern <- Upattern.name[1]
            iM <- matrix(nrow = Upattern.nRep[iPattern==Upattern.name], ncol = 0)
            attr(iM,"index.cluster") <- pattern.cluster$index.cluster[pattern.cluster$var==iPattern]
            attr(iM,"index.time") <- index.clusterTime[[attr(iM,"index.cluster")[1]]]
            return(iM)
        }),Upattern.name)
        structure$X$Xpattern.cor <- NULL
        structure$X$pattern.cluster <- pattern.cluster
        return(structure)
    }
    test.duplicated <- duplicated(lpnCluster)
    name.pattern <- lpnCluster[!test.duplicated]
    order.pattern <- order(name.pattern)
    name.pattern <- name.pattern[order.pattern]
    n.pattern <- length(name.pattern)

    pattern.cluster <- data.frame(index.cluster = U.cluster,
                                  pattern = factor(name.pattern[match(lpnCluster, name.pattern)], levels = name.pattern))

    ## *** var/cor pattern corresponding to pattern
    if(NCOL(X.var)>0){
        pattern.indexVar <- lpnCluster.var[!test.duplicated][order.pattern]
    }else{
        pattern.indexVar <- NA
    }
    if(!is.null(X.cor)){
        pattern.indexCor <- lpnCluster.cor[!test.duplicated][order.pattern]
    }else{
        pattern.indexCor <- NA
    }
    Upattern <- data.frame(name = name.pattern,
                           var = pattern.indexVar,
                           cor = pattern.indexCor,
                           index.strata = NA,
                           time = NA,
                           n.time = NA,
                           param = NA,
                           index.cluster = NA,
                           n.cluster = NA)

    Upattern$index.cluster <- stats::setNames(lapply(name.pattern, function(iPattern){which(lpnCluster == iPattern)}), name.pattern)
    Upattern$n.cluster <- sapply(Upattern$index.cluster, length)
    Upattern$time <- lapply(Upattern$index.cluster, function(iId){
        index.clusterTime[[iId[1]]]
    })
    Upattern$n.time <- sapply(Upattern$time,length)
    Upattern$index.strata <- sapply(Upattern$index.cluster, function(iId){unname(unique(index.clusterStrata[iId]))}, simplify = TRUE)
    attr(Upattern, "level.var") <- name.pattern.var
    attr(Upattern, "level.cor") <- name.pattern.cor
    attr(Upattern, "lp.UpatternVar") <- Ulp.var
    pattern.cluster$var <- stats::setNames(stats::setNames(Upattern$var, Upattern$name)[pattern.cluster$pattern], pattern.cluster$cluster)
    if(!is.null(X.cor)){
        attr(Upattern, "lp.UpatternCor") <- Ulp.cor
        pattern.cluster$cor <- stats::setNames(stats::setNames(Upattern$cor, Upattern$name)[pattern.cluster$pattern], pattern.cluster$cluster)
    }

    ## ** characterize each pattern

    ## *** var
    param.sigma <- param[param$type=="sigma","name"]
    param.k <- param[param$type=="k","name"]
    Xpattern.var <- stats::setNames(lapply(name.pattern.var,function(iP){ ## iP <- name.pattern.var[1]
        iC.all <- pattern.cluster$index.cluster[which(pattern.cluster$var==iP)]
        iC <- iC.all[1]
        iX <- X.var[index.cluster[[iC]],,drop=FALSE]
        iIndex.param <- which(colSums(iX)!=0)
        iParam.k <- intersect(names(iIndex.param),param.k)
        iParam.sigma <- intersect(names(iIndex.param),param.sigma)
        
        attr(iX,"index.cluster") <- iC.all
        attr(iX,"index.strata") <- unname(index.clusterStrata[iC])
        attr(iX,"param") <- names(iIndex.param)
        attr(iX,"Mindicator.param") <- stats::setNames(lapply(iIndex.param,function(iCol){tcrossprod(iX[,iCol],rep(1,NROW(iX))) + t(tcrossprod(iX[,iCol],rep(1,NROW(iX))))}),
                                                       names(iIndex.param))
        attr(iX,"indicator.param") <- lapply(attr(iX,"Mindicator.param"),function(iM){which(iM>0)})
        
        return(iX)
    }),name.pattern.var)

    Upattern$param <- stats::setNames(lapply(Upattern$name, function(iP){ ## iP <- Upattern$name[2]
        iP.var <- Upattern[["var"]][[which(Upattern$name==iP)]]
        return(attr(Xpattern.var[[iP.var]],"param"))
    }), Upattern$name)

    ## *** cor
    if(!is.null(X.cor)){
        param.rho <- param[param$type=="rho",,drop=FALSE]
        ls.lp.xy <- mapply(x = param.rho$code.x, y = param.rho$code.y, FUN = paste, sep = "::X::", SIMPLIFY = FALSE)
        lp.xy <- unlist(lapply(names(ls.lp.xy), function(iN){stats::setNames(ls.lp.xy[[iN]],rep(iN, length(ls.lp.xy[[iN]])))}))
        UlpnCluster.cor0 <- lpnCluster.cor0[cluster.pattern.cor]
        
        iN.pair <- unique(sapply(UlpnCluster.cor0, length))
        ls.pair <- vector(mode = "list", length = max(iN.pair))

        ls.pair[iN.pair[iN.pair>0]] <- lapply(iN.pair[iN.pair>0], function(iN){.unorderedPairs(1:iN, distinct = TRUE)})
        Xpattern.cor <- stats::setNames(lapply(setdiff(name.pattern.cor,""),function(iP){ ## iP <- name.pattern.cor[2]
            ## ignore pattern "", i.e. patterns with no pair
            iC.all <- pattern.cluster$index.cluster[which(pattern.cluster$cor==iP)]
            iC <- iC.all[1]
            iX <- X.cor[index.cluster[[iC]],,drop=FALSE]
            iC.code <- lpnCluster.cor0[[iC]]
            iC.pair <- ls.pair[[length(iC.code)]]
            iC.codepair <- matrix(iC.code[iC.pair], ncol = 2, byrow = TRUE, dimnames = list(NULL,c("lp.x","lp.y")))
            if(!is.null(attr(param,"level.cor"))){
                iC.lp.xy <- paste0(attr(param,"level.cor")[iC.codepair[,"lp.x"]],sep[2],attr(param,"level.cor")[iC.codepair[,"lp.y"]])
                iC.param <- names(lp.xy)[match(iC.lp.xy,lp.xy)] ## can contain NA when asking a new structure without having seen a particular pair
            }else{
                iC.param <- paste(param.rho[param.rho$index.strata == index.clusterStrata[iC],"name"], collapse = "|")
            }
            ## unique(iC.lp.xy)
            iC.table <- rbind(data.frame(row = iC.pair[1,], col = iC.pair[2,], param = iC.param),
                              data.frame(row = iC.pair[2,], col = iC.pair[1,], param = iC.param))
            ## iX.pairwise <- matrix(0, nrow = NROW(iC.table), ncol = length(param.rho$name),
            ##                       dimnames = list(paste0("(",iC.table$row,",",iC.table$col,")"), param.rho$name))
            ## iX.pairwise[1:NROW(iC.table) + (match(iC.table$param , param.rho$name)-1) * NROW(iC.table)] <- 1
            attr(iX, "index.cluster") <- iC.all
            if(sum(!duplicated(index.clusterTime[iC.all]))==1){
                attr(iX,"index.time") <- index.clusterTime[[iC]]
            }
            attr(iX, "index.strata") <- unname(index.clusterStrata[iC])
            attr(iX, "index.pair") <- iC.table
            attr(iX, "index.vec2matrix") <- c(iC.table[,"row"] + NROW(iX) * (iC.table[,"col"] - 1))
            if(structure$type=="CS" && structure$heterogeneous<0){
                attr(iX, "param") <- stats::na.omit(unique(iC.param))
            }else{
                attr(iX, "param") <- unique(iC.param)
            }                
            attr(iX, "indicator.param") <- stats::setNames(lapply(attr(iX, "param"),function(iP){
                if(is.na(iP)){ ## deal with NA for the case prediction where a pair of time may not have been observed
                    return(attr(iX, "index.vec2matrix")[which(is.na(iC.table$param))])
                }else{
                    return(attr(iX, "index.vec2matrix")[which(iC.table$param==iP)])
                }
            }), attr(iX, "param"))
            attr(iX,"Mindicator.param") <- stats::setNames(lapply(1:length(attr(iX, "param")), function(iParam){
                iM <- matrix(0, nrow = NROW(iX), ncol = NROW(iX))
                iM[attr(iX,"indicator.param")[[iParam]]] <- 1
                return(iM)
            }),attr(iX, "param"))
            return(iX)
            ## iM <- matrix(0,length(iC.code),length(iC.code)); iM[attr(iX.pairwise, "index.vec2matrix")] <- iX.pairwise %*% 1:NROW(param.rho)
        }),setdiff(name.pattern.cor,""))

        Upattern$param <- stats::setNames(lapply(Upattern$name, function(iP){ ## iP <- Upattern$name[1]
            iPattern.cor <- Upattern$cor[which(Upattern$name==iP)]
            iParam <- Upattern$param[[which(Upattern$name==iP)]]
            if(iPattern.cor==""){
                ## pattern "", i.e. patterns with no pair
                return(iParam)
            }else{
                return(c(iParam, attr(Xpattern.cor[[iPattern.cor]],"param")))
            }
        }), Upattern$name)
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

## * .findUpatterns.RE
.findUpatterns.RE <- function(structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata,
                              sep = c(".","::X::")){


    browser()
    ## ** random effects
    if(inherits(structure,"RE")){
        ## nested
        ## decode: R.1.1.1 innermost element (R  =same school/class/student, 1.1.1=strata 1)
        ##       : D.0.0.1                   (0.0=same school/class,             1=different student)
        ##       : D.0.1.1                   (0  =same school,                 1.1=different student/class)

        ## crossed
        ## decode: D.0.1.0                   (same region)
        ##       : D.0.0.1                   (same timepoint)

        Xvar.ranef <- colnames(X.cor)
        structure.paramRho <- structure$param[structure$param$type=="rho",,drop=FALSE]
        
        structure$ranef$param <- do.call(rbind,lapply(1:length(U.strata), function(iS){ ## iS <- 1
            iPattern.cor <- Upattern[Upattern$index.strata==iS,"cor"]
            iRho <- structure.paramRho[structure.paramRho$index.strata==iS,"name"]
            iRho.code <- structure.paramRho[structure.paramRho$index.strata==iS,"code"]

            iRho.splitcode <- do.call(rbind,strsplit(iRho.code, split = ".", fixed = TRUE))
            iRho.splitcode2 <- iRho.splitcode[,-1,drop=FALSE]
            iRho.splitcode2[iRho.splitcode[,1]=="R"] <- "0" ## specific coding within block
            dimnames(iRho.splitcode2) <- list(iRho,Xvar.ranef)
            
            Xvar.ranef[rowSums(iRho.splitcode2=="1")+1]
            apply(iRho.splitcode2)
            if(length(iRho.code)==1){
                iRho.ordered <- iRho
                iRho.vars <- structure$ranef$vars
                iRho.column <- names(which(colSums(do.call(rbind,Xpattern.cor[iPattern.cor]))>0))
            }else if(attr(structure$ranef,"nested")){
                
                ##structure$ranef$hierarchy[[1]]
                iRho.columns <- rep(as.character(NA), length = length(iRho))
                
                index.firstRho <- which(grepl("R",iRho.code, fixed = TRUE))
                if(length(index.firstRho)!=1){
                    stop("Cannot not link random effect structure to correlation parameters. \n",
                         "Issue with the reference level. \n")
                }
                iRho.columns[-index.firstRho] <- sapply(iRho.splitcode, function(iCode){rev(Xvar.ranef)[iCode[-1]=="1"][1]})
                iRho.columns[index.firstRho] <- setdiff(Xvar.ranef, iRho.columns[-index.firstRho])
                
                    
                iRho.ordered <- iRho[index.firstRho]
                if(length(index.firstRho)!=0){
                    stop("Cannot not link random effect structure to correlation parameters. \n",
                         "Multiple reference levels. \n")
                }else{
                    iIndicator <- sapply(iRho.splitcode, function(iCode){Xvar.ranef[iCode[-1]=="1"][1]})
                    if(any(sort(iIndicator)!=sort(iRho.vars[-1]))){
                        stop("Cannot not link random effect structure to correlation parameters. \n",
                             "Mismatch between planned and reconstructed variables. \n")
                    }
                    match(iRho.vars,iIndicator,)
                    iRho.ordered <- c(iRho.ordered, iRho[-index.firstRho])                    
                }
            }else if(attr(structure$ranef,"crossed")){

                iRho.splitcode <- strsplit(iRho.code, split = ".", fixed = TRUE)
                iIndicator <- lapply(iRho.splitcode, function(iCode){Xvar.ranef[iCode[-1]=="1"]})
                if(any(sapply(iIndicator,length)==0)){
                    stop("Cannot not link random effect structure to correlation parameters. \n",
                         "Correlation parameter no linked to any random effect.")
                }
                if(any(sapply(iIndicator,length)>1)){
                    stop("Cannot not link random effect structure to correlation parameters. \n",
                         "Correlation parameter linked to multiple random effect.")
                }
                if(any(duplicated(unlist(iIndicator)))){
                    stop("Cannot not link random effect structure to correlation parameters. \n",
                         "Correlation parameters linked to the same random effect.")
                }
                unlist(iIndicator)
                Xvar.ranef
                
            }
            iDF <- data.frame(strata = iS, param = iRho.ordered, level = iRho.vars, column = iRho.column)
            return(iDF)
        }))
    }

}


## * .findUpatterns.TOEPLITZ
.findUpatterns.TOEPLITZ <- .findUpatterns.CS

## * .findUpatterns.UN
.findUpatterns.UN <- function(structure, 
                              index.clusterTime, U.time,
                              index.cluster, U.cluster,
                              index.clusterStrata, U.strata,
                              sep = c(".","::X::")){

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
    Ulp.cor <- structure.param$code
    n.cluster <- length(U.cluster)

    corLP.cluster <- attr(structure$X$lp.cor,"cluster")
    
    ## *** associate a correlation pattern to each cluster
    ## vector summarizing the (combined) linear predictor of each cluster 
    patternCorLP.cluster <- sapply(corLP.cluster, function(iLp){ 
        paste0(as.numeric(iLp), collapse = sep[1])
    })
    ## list of cluster per (combined) linear predictor
    patternCorLP.index <- tapply(1:n.cluster,patternCorLP.cluster,identity)
    ## unique (combined) linear predictors
    patternCorLP.index1 <- sapply(patternCorLP.index,"[",1)
    UpatternCorLP <- patternCorLP.cluster[patternCorLP.index1]
    ## unique design matrix for pairwise comparisons
    UpatternCorDiff <- lapply(patternCorLP.index1, function(iC){ ## iC <- 1
        iIndex <- .unorderedPairs(1:length(index.cluster[[iC]]), distinct = TRUE)
        iIndex.cluster <- rbind(index.cluster[[iC]][iIndex[1,]], index.cluster[[iC]][iIndex[2,]])
        iX.cor1 <- X.cor[iIndex.cluster[1,],,drop=FALSE]
        iDiff <- X.cor[iIndex.cluster[2,],,drop=FALSE] - iX.cor1
        iCode <- paste0("D",interaction(as.data.frame(iX.cor1),sep=""),sep[1],interaction(as.data.frame(iDiff),sep=sep[1]))
        iParam <- structure.param
        iMatch <- match(iCode,structure.param$code)
        if(any(is.na(iMatch))){
            warning("Something went wrong when extracting the correlation patterns. \n")
        }
        return(data.frame(index1 = iIndex[1,],
                          index2 = iIndex[2,],
                          ## index.time1 = attr(index.clusterTime,"vectorwise")[iIndex.cluster[1,]],
                          ## index.time2 = attr(index.clusterTime,"vectorwise")[iIndex.cluster[2,]],
                          code = iCode,
                          param  = structure.param[iMatch,"name"]))
    })
    
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
    UpatternCor$param <- lapply(UpatternCorDiff,function(iPattern){unique(iPattern$param)})

    patternCor.cluster <- do.call(rbind,lapply(UpatternCorLP, function(iPattern){
        data.frame(index.cluster = patternCorLP.index[[iPattern]], cor = iPattern)
    }))
    patternCor.cluster <- patternCor.cluster[order(patternCor.cluster$index.cluster),,drop=FALSE]

    ## *** characterize each correlation pattern
    Xpattern.cor <- lapply(UpatternCorLP, function(iPattern){ ## iPattern <- UpatternCorLP[1]

        iCluster <- patternCorLP.index1[[iPattern]]
        iIndex.cluster <- index.cluster[[iCluster]]
        iRep <- length(iIndex.cluster)
        iOut <- X.cor[iIndex.cluster,,drop=FALSE]
        attr(iOut, "index.cluster") <- patternCorLP.index[[iPattern]]
        attr(iOut, "index.time") <- index.clusterTime[[iCluster]]
        attr(iOut, "index.strata") <- index.clusterStrata[[iCluster]]
        attr(iOut, "index.pair") <- rbind(data.frame(row =  UpatternCorDiff[[iPattern]][,c("index1")],
                                                     col = UpatternCorDiff[[iPattern]][,c("index2")],
                                                     param = UpatternCorDiff[[iPattern]][,c("param")]),
                                          data.frame(row =  UpatternCorDiff[[iPattern]][,c("index2")],
                                                     col = UpatternCorDiff[[iPattern]][,c("index1")],
                                                     param = UpatternCorDiff[[iPattern]][,c("param")]))
        attr(iOut, "index.vec2matrix") <- attr(iOut, "index.pair")$row + (attr(iOut, "index.pair")$col-1)*iRep
        attr(iOut, "param") <- unique(UpatternCorDiff[[iPattern]]$param)
        attr(iOut, "indicator.param") <- tapply(attr(iOut, "index.vec2matrix"),attr(iOut, "index.pair")$param,identity,simplify=FALSE)
        attr(iOut, "Mindicator.param") <- tapply(attr(iOut, "index.vec2matrix"),attr(iOut, "index.pair")$param,function(iIndex){
            iM <- matrix(0, nrow = iRep, ncol = iRep)
            iM[iIndex] <- 1
            return(iM)
        }, simplify=FALSE)

        return(iOut)
    })

    ## ** joint variance and correlation patterns
    
    pattern.name <- paste0(as.numeric(factor(patternVar.cluster$var,attr(UpatternVar,"level.var"))),
                           ":",
                           as.numeric(factor(patternCor.cluster$cor,UpatternCorLP)))
    pattern.cluster <- data.frame(index.cluster = 1:n.cluster,
                                  pattern = pattern.name,
                                  var = patternVar.cluster$var,
                                  cor = patternCor.cluster$cor)
    index.Upattern <- which(!duplicated(pattern.name))
    
    Upattern <- data.frame(name = pattern.cluster$pattern[index.Upattern],
                           var = pattern.cluster$var[index.Upattern],
                           cor = pattern.cluster$cor[index.Upattern],
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
    Upattern$param <- 
    UpatternVar$param[UpatternVar$var==Upattern$var]
    UpatternCor$param[UpatternCor$cor==Upattern$cor]
    attr(Upattern,"level.var") <- attr(UpatternVar,"level.var")
    attr(Upattern,"level.cor") <- UpatternCorLP
        
browser()
    ## ** export
    structure$X$Upattern <- Upattern
    structure$X$Xpattern.var <- Xpattern.var
    structure$X$Xpattern.cor <- Xpattern.cor
    structure$X$pattern.cluster <- pattern.cluster
    return(structure)
}

## * .findUpatterns.EXP
.findUpatterns.EXP <- .findUpatterns.ID

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
