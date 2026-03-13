### structure-skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: mar 13 2026 (14:37) 
##           By: Brice Ozenne
##     Update #: 2602
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * skeleton
##' @title Parametrization of Covariance Structure
##' @description Parametrization of Covariance Structure
##' @noRd
##'
##' @param structure [structure]
##' @param data [data.frame] dataset
##' @param indexData pre-computed quantity
##'
##' @keywords internal
##' 
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$gender <- c("M","F")[as.numeric(gastricbypassL$id) %% 2+1]
##' dd <- gastricbypassL[!duplicated(gastricbypassL[,c("time","gender")]),]
##' 
##' ## independence
##' .skeleton(IND(~1, var.cluster = "id", var.time = "time"), data = dd)
##' .skeleton(IND(gender~1, var.cluster = "id", var.time = "time"), data = dd)
##' .skeleton(IND(~time, var.time = "time", var.cluster = "id"), data = dd)
##' .skeleton(IND(~time+gender, var.time = "time", var.cluster = "id"), data = dd)
##' .skeleton(IND(gender~time, var.time = "time", var.cluster = "id"), data = dd)
##' 
##' ## compound symmetry
##' .skeleton(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' .skeleton(CS(~time|id), data = gastricbypassL)
##' .skeleton(CS(gender~time|id), data = gastricbypassL)
##' 
##' ## unstructured
##' .skeleton(UN(NULL, var.cluster = "id", var.time = "visit", TRUE), data = dd)
##' .skeleton(UN(~gender, var.cluster = "id", var.time = "visit", TRUE), data = dd)
`.skeleton` <-
    function(structure, data, indexData, options) UseMethod(".skeleton")

## * skeleton.ID
.skeleton.ID <- function(structure, data, indexData = NULL, options = NULL){

    ## ** initialize (not use by LMMstar, only if direct call from the user like in the examples of the documentation)
    if(is.null(indexData)){
        indexData <- .extractIndexData(data = data, structure = structure)
    }
    if(is.null(structure$var) && is.null(structure$cor)){
        outDesign <- .vcov.matrix.lmm(structure = structure, data = data, index.cluster = index.cluster, drop.X = options$drop.X, sep = options$sep["lp"])
        structure$xfactor <- outDesign$xfactor
        structure$var <- outDesign$var
    }

    ## ** prepare
    if(is.null(options)){
        options <- LMMstar.options()
    }
    U.strata <- indexData$U.strata
    
    ## **  variance parameters
    structure <- .skeletonSigma(structure, U.strata = U.strata)
    structure <- .skeletonK(structure, U.strata = U.strata, sep = options$sep[c("k.cov","k.strata")])

    ## ** export
    structure$param <- structure$param[order(structure$param$index.strata),,drop=FALSE]
    rownames(structure$param) <- NULL
    return(structure)
}


## * skeleton.IND
.skeleton.IND <- .skeleton.ID

## * skeleton.CS
.skeleton.CS <- function(structure, data, indexData = NULL, options = NULL){


    ## ** initialize (not use by LMMstar, only if direct call from the user like in the examples of the documentation)
    if(is.null(indexData)){
        indexData <- .extractIndexData(data = data, structure = structure)
    }
    if(is.null(structure$var) && is.null(structure$cor)){
        outDesign <- .vcov.matrix.lmm(structure = structure, data = data, index.cluster = index.cluster, drop.X = options$drop.X, sep = options$sep["lp"])
        structure$xfactor <- outDesign$xfactor
        structure$var <- outDesign$var
        structure$cor <- outDesign$cor
    }


    ## ** prepare
    if(is.null(options)){
        options <- LMMstar.options()
    }
    U.cluster <- indexData$U.cluster
    U.time <- indexData$U.time
    U.strata <- indexData$U.strata
    n.strata <- length(U.strata)

    index.clusterTime <- indexData$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    index.clusterStrata <- indexData$index.clusterStrata ## vector of index relative to which strata each cluster belong to
    index.cluster <- indexData$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    ## **  variance parameters
    structure <- .skeletonSigma(structure, U.strata = U.strata)
    structure <- .skeletonK(structure, U.strata = U.strata, sep = options$sep[c("k.cov","k.strata")])

    param.sigma <- structure$param[structure$param$type=="sigma","name"]
    strataIndex.sigma <- structure$param[structure$param$type=="sigma","index.strata"]
    param.k <- structure$param[structure$param$type=="k","name"]

    ## **  correlation parameters
    ## correlation covariates that are not time nor strata
    name.cov <- na.omit(setdiff(structure$name$cor[[1]], c(structure$name$strata,structure$name$time)))

    ## *** identify unique pairs of linear predictors
    XpairPattern <- .pairPatternX(structure$cor,
                                  name.cov = name.cov, name.time = structure$name$time,
                                  U.cluster = U.cluster, index.cluster = index.cluster,
                                  U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    ## *** between-block (pairs with distinct linear predictors)
    if(NROW(XpairPattern$lpB) > 0){ 
        outCrossBlock <- switch(structure$class["correlation.cross"],
                                ID = .skeletonRho.ID(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "B", sep = options$sep),
                                CS = .skeletonRho.CS(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "B", sep = options$sep),
                                AR1 = .skeletonRho.AR1(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "B", sep = options$sep),
                                EXP = .skeletonRho.EXP(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "B", sep = options$sep),
                                TOEPLITZ = .skeletonRho.TOEPLITZ(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "B", sep = options$sep),
                                DUN = .skeletonRho.DUN(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "B", sep = options$sep),
                                UN = .skeletonRho.UN(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "B", sep = options$sep)
                                )
        if(!is.null(outCrossBlock$nesting)){
            structure$nesting <- outCrossBlock$nesting
        }
    }else{
        outCrossBlock <- list(code = NULL,
                              level = NULL,
                              strata = NULL,
                              nesting = NULL,
                              constraint = NULL,
                              indexObs = NULL)
    }

    ## *** within-block (pairs with same linear predictors)
    if(NROW(XpairPattern$lpW) > 0){ 
        outDiagBlock <- switch(structure$class["correlation"],
                               ID = .skeletonRho.ID(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "W", sep = options$sep),
                               CS = .skeletonRho.CS(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "W", sep = options$sep),
                               AR1 = .skeletonRho.AR1(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "W", sep = options$sep),
                               EXP = .skeletonRho.EXP(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "W", sep = options$sep),
                               TOEPLITZ = .skeletonRho.TOEPLITZ(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "W", sep = options$sep),
                               DUN = .skeletonRho.DUN(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "W", sep = options$sep),
                               UN = .skeletonRho.UN(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = "W", sep = options$sep)
                               )
    }else{
        outDiagBlock <- list(code = NULL,
                             level = NULL,
                             strata = NULL,
                             nesting = NULL,
                             constraint = NULL,
                             indexObs = NULL)
    }

    ## ***  retrieve k
    if(length(param.k)>0){
        ## from linear predictor (variance) to k parameter
        lp2param.k <- structure$param$name[match(rownames(structure$var$lp2X),structure$param$code)]
        lp2param.k[lp2param.k %in% param.sigma] <- NA
        
        ## from observation used to create the correlation parameter to k parameter
        k.xW <- lapply(outDiagBlock$indexObs, function(iM){lp2param.k[structure$var$lp[iM[,1]]]})
        k.yW <- lapply(outDiagBlock$indexObs, function(iM){lp2param.k[structure$var$lp[iM[,2]]]})
        
        k.xB <- lapply(outCrossBlock$indexObs, function(iM){lp2param.k[structure$var$lp[iM[,1]]]})
        k.yB <- lapply(outCrossBlock$indexObs, function(iM){lp2param.k[structure$var$lp[iM[,2]]]})
    }
    
    ## ***  collect    
    code.rho <- unname(c(outDiagBlock$code, outCrossBlock$code))
    level.rho <- unname(c(outDiagBlock$level, outCrossBlock$level))
    strata.rho <- unname(c(outDiagBlock$strata, outCrossBlock$strata))

    structure.rho <- data.frame(name = paste0("rho",level.rho),
                                index.strata = strata.rho,
                                type = rep("rho",length=length(level.rho)),
                                constraint = NA_real_,
                                level = level.rho,
                                code = NA,
                                sigma = param.sigma[match(strata.rho,strataIndex.sigma)],
                                k.x = NA,
                                k.y = NA,                                  
                                stringsAsFactors = FALSE)
    structure.rho$code <- code.rho
    if("k" %in% structure$param$type){
        structure.rho$k.x <- c(k.xW, k.xB)
        structure.rho$k.y <- c(k.yW, k.yB)
    }

    ## *** export
    rownames(structure.rho) <- NULL
    structure$param <- rbind(structure$param, structure.rho)
    return(structure)
}

## * skeleton.RE
.skeleton.RE <- .skeleton.CS

## * skeleton.TOEPLITZ
.skeleton.TOEPLITZ <- .skeleton.CS

## * skeleton.UN
.skeleton.UN <- .skeleton.CS

## * skeleton.DUN
.skeleton.DUN <- .skeleton.CS

## * skeleton.EXP
.skeleton.EXP <- .skeleton.CS

## * skeleton.AR1
.skeleton.AR1 <- .skeleton.CS


## * skeleton.CUSTOM
.skeleton.CUSTOM <- function(structure, data, indexData = NULL, options = NULL){

    ## ** prepare
    var.strata <- structure$name$strata
    if(is.null(options)){
        options <- LMMstar.options()
    }
    if(is.null(indexData)){
        indexData <- .extractIndexData(data = data, structure = structure)
    }
    U.cluster <- indexData$U.cluster
    U.time <- indexData$U.time
    U.strata <- indexData$U.strata
    n.strata <- length(U.strata)
    
    sep.strata <- c(sigma = unname(options$sep["k.strata"]),
                    rho = unname(options$sep["rho.strata"]))
    
    ## ** gather parameters
    structure$param <- NULL
    if(!is.null(structure$FCT.sigma)){
        if(is.na(var.strata)){
            name.sigma <- names(structure$init.sigma)
            strata.sigma <- 1
        }else{
            init.name.sigma <- names(structure$init.sigma)
            init.n.sigma <- length(init.name.sigma)

            ls.name.sigma <- lapply(U.strata, function(iStrata){stats::setNames(paste(init.name.sigma,iStrata,sep=sep.strata["sigma"]), rep(iStrata, init.n.sigma))})
            name.sigma <- unname(unlist(ls.name.sigma, use.names = FALSE))
            strata.sigma <- match(unlist(lapply(ls.name.sigma, names)), U.strata)

            structure$init.sigma <- unlist(lapply(U.strata, function(iStrata){stats::setNames(structure$init.sigma[init.name.sigma],paste(init.name.sigma,iStrata,sep=sep.strata["sigma"]))}))
        }
        structure$param <- rbind(structure$param,
                                 data.frame(name = name.sigma,
                                            index.strata = strata.sigma,
                                            type = "sigma",
                                            constraint = as.numeric(NA),
                                            level = as.character(NA),
                                            code = as.character(NA),
                                            lp.x = NA,
                                            lp.y = NA,
                                            sigma = as.character(NA),
                                            k.x = as.character(NA),
                                            k.y = as.character(NA)))
        
    }

    if(!is.null(structure$FCT.rho)){
        if(is.na(var.strata)){
            name.rho <- names(structure$init.rho)
            strata.rho <- 1
        }else{
            init.name.rho <- names(structure$init.rho)
            init.n.rho <- length(init.name.rho)

            ls.name.rho <- lapply(U.strata, function(iStrata){stats::setNames(paste(init.name.rho,iStrata,sep=sep.strata["rho"]), rep(iStrata, init.n.rho))})
            name.rho <- unname(unlist(ls.name.rho, use.names = FALSE))
            strata.rho <- match(unlist(lapply(ls.name.rho, names)), U.strata)

            structure$init.rho <- unlist(lapply(U.strata, function(iStrata){stats::setNames(structure$init.rho[init.name.rho],paste(init.name.rho,iStrata,sep=sep.strata["rho"]))}))
         }
        structure$param <- rbind(structure$param,
                                 data.frame(name = name.rho,
                                            index.strata = strata.rho,
                                            type = "rho",
                                            constraint = as.numeric(NA),
                                            level = as.character(NA),
                                            code = as.character(NA),
                                            lp.x = NA,
                                            lp.y = NA,
                                            sigma = as.character(NA),
                                            k.x = as.character(NA),
                                            k.y = as.character(NA)))
    }

    ## ** export
    return(structure)
}





##----------------------------------------------------------------------
### structure-skeleton.R ends here
