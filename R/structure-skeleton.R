### structure-skeleton.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  8 2021 (17:56) 
## Version: 
## Last-Updated: jul 11 2023 (12:14) 
##           By: Brice Ozenne
##     Update #: 2486
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
    function(structure, data, indexData) UseMethod(".skeleton")

## * skeleton.ID
.skeleton.ID <- function(structure, data, indexData = NULL){

    ## ** prepare
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
        outDesign <- .vcov.matrix.lmm(structure = structure, data = data, index.cluster = index.cluster, drop.X = LMMstar.options()$drop.X)
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
.skeleton.CUSTOM <- function(structure, data, indexData = NULL){

    ## ** gather parameters
    structure$param <- NULL
    if(!is.null(structure$FCT.sigma)){
        structure$param <- rbind(structure$param,
                                 data.frame(name = names(structure$init.sigma),
                                            strata = 1,
                                            type = "sigma",
                                            constraint = as.numeric(NA),
                                            level = as.character(NA),
                                            code = as.character(NA),
                                            index.lp.x = NA,
                                            index.lp.y = NA,
                                            sigma = as.character(NA),
                                            k.x = as.character(NA),
                                            k.y = as.character(NA)))
    }

    if(!is.null(structure$FCT.rho)){
        structure$param <- rbind(structure$param,
                                 data.frame(name = names(structure$init.rho),
                                            strata = 1,
                                            type = "rho",
                                            constraint = as.numeric(NA),
                                            level = as.character(NA),
                                            code = as.character(NA),
                                            index.lp.x = NA,
                                            index.lp.y = NA,
                                            sigma = as.character(NA),
                                            k.x = as.character(NA),
                                            k.y = as.character(NA)))
    }
    structure$param$index.lp.x <- rep(list(NULL), NROW(structure$param))
    structure$param$index.lp.y <- rep(list(NULL), NROW(structure$param))

    ## ** pattern
    return(structure)
}





##----------------------------------------------------------------------
### structure-skeleton.R ends here
