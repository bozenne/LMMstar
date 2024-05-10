### structure-skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: maj 10 2024 (11:42) 
##           By: Brice Ozenne
##     Update #: 2502
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

    ## ** prepare
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
    
    index.clusterTime <- indexData$index.clusterTime ## list of index relative to the time at which the observations are made within cluster
    index.clusterStrata <- indexData$index.clusterStrata ## vector of index relative to which strata each cluster belong to
    index.cluster <- indexData$index.cluster ## list of positions of the observation belonging to each cluster in the dataset

    if(is.null(structure$var) && is.null(structure$cor)){
        outDesign <- .vcov.matrix.lmm(structure = structure, data = data, index.cluster = index.cluster, drop.X = options$drop.X, sep = options$sep["lp"])
        structure$xfactor <- outDesign$xfactor
        structure$var <- outDesign$var
        structure$cor <- outDesign$cor
    }

    ## **  param
    structure <- .skeletonSigma(structure, U.strata = U.strata)
    structure <- .skeletonK(structure, U.strata = U.strata)
    structure <- .skeletonRho(structure, data = data, 
                              U.cluster = U.cluster, index.cluster = index.cluster,
                              U.time = U.time, index.clusterTime = index.clusterTime, 
                              U.strata = U.strata, index.clusterStrata = index.clusterStrata)

    ## ** export
    structure$param <- structure$param[order(structure$param$index.strata),,drop=FALSE]
    rownames(structure$param) <- NULL
    return(structure)
}


## * skeleton.IND
.skeleton.IND <- .skeleton.ID

## * skeleton.CS
.skeleton.CS <- .skeleton.ID

## * skeleton.RE
.skeleton.RE <- .skeleton.ID

## * skeleton.TOEPLITZ
.skeleton.TOEPLITZ <- .skeleton.ID

## * skeleton.UN
.skeleton.UN <- .skeleton.ID

## * skeleton.EXP
.skeleton.EXP <- .skeleton.ID

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
