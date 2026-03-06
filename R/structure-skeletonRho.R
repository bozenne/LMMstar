### structure-skeletonRho.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 11 2023 (13:27) 
## Version: 
## Last-Updated: mar  6 2026 (15:52) 
##           By: Brice Ozenne
##     Update #: 1650
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .skeletonRho
##' @title Parametrization of Correlation Structure.
##' @description Parametrization of the correlation structure.
##' @noRd
##'
##' @param structure [structure]
##' @param XpairPattern [list] linear predictor, strata, index of observations for the unique combinations of pairs.
##' @param U.strata [character vector] time levels.
##' @param name.cov [character vector] covariate influencing the correlation, other than strata and time.
##' @param block [character] Diagonal block \code{"W"} (pair of observations with same linear predictor, up to time)
##' or off diagonal block \code{"B"} (pair of observations with different linear predictor, up to time).
##' @param sep [character vector of length 2] characters used to name the variance parameters,
##' the first is used to specify the covariate level whereas
##' the second is used to specify the strata level (when more than 1 strata).
##'
##' @keywords internal

`.skeletonRho` <-
    function(structure, XpairPattern, U.strata, name.cov, block, sep) UseMethod(".skeletonRho")

## * .skeletonRho.ID (cross random effects)
.skeletonRho.ID <- function(structure, XpairPattern, U.strata, name.cov, block, sep){

    ## ** run 'normal' CS
    out <- .crossBlockCS(structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = block, sep = sep)

    ## ** add constraint
    lp2data <- structure$cor$lp2data
    lpB <- XpairPattern$lpB

    diffS.nesting <- rowSums(lp2data[lpB[,1],name.cov,drop=FALSE]==lp2data[lpB[,2],name.cov,drop=FALSE], as.factor = FALSE)
    out$constraint <- unname(ifelse(diffS.nesting==0,0,NA))

    ## ** export
    return(out)

}

## * .skeletonRho.CS
.skeletonRho.CS <- function(structure, XpairPattern, U.strata, name.cov, block, sep){

    ## ** extract data
    lp <- XpairPattern[[paste0("lp",block)]]
    strata <- XpairPattern[[paste0("strata",block)]]
    index <- XpairPattern[[paste0("index",block)]]
    lp2data <- structure$cor$lp2data
    name.cor <- unlist(structure$name$cor)
    name.strata <-  structure$name$strata
    name.cluster <- attr(structure$name$cluster,"original")
    name.nesting <- structure$nesting
    n.strata <- length(U.strata)
    
    ## ** find name for correlation coefficient, e.g. rho(id)
    if(block=="W"){
        cluster.var <- stats::na.omit(c(attr(structure$name$cluster,"original"),structure$name$cluster))[1]
        if(cluster.var=="XXcluster.indexXX"){
            label.clusterW <- NULL            
        }else if(!is.null(structure$nesting)){ ## instead of rho(id) which should be for the cross block element, uses rho(id/...) for the diagonal block
            label.clusterW <- paste(rev(name.nesting),collapse = "/")
        }else{
            label.clusterW <- cluster.var
        }
    }
    
    ## ** identify parameters
    if(block == "W"){
        if(structure$twin){ ## identical blocks on the diagonal or single block (no covariates)
            ## internal indentifier for the parameter (one parameter per strata)
            code.rho <- tapply(1:NROW(lp), strata, FUN = function(iVec){
                paste("W",lp[iVec,1,drop=FALSE],lp[iVec,2,drop=FALSE],sep=sep["rho.name"])
            }, simplify = FALSE)
            ## strata relative to the parameter
            strata.rho <- tapply(strata, strata, FUN = unique)
            ## provide a name for the parameter
            if(n.strata==1){
                level.rho <- paste(paste0("(",label.clusterW,")"), sep = sep[3])
            }else{
                level.rho <- paste(paste0("(",label.clusterW,")"),U.strata[strata.rho], sep = sep[3])
            }
            ## index of the corresponding observations
            indexObs.rho <- tapply(1:NROW(index), strata, FUN = function(iVec){index[iVec,,drop=FALSE]})
        }else{ ## distinct blocks on the diagonal
            ## internal indentifier for the parameter (one parameter per strata and covariate, i.e. per linear predictor)
            code.rho <- as.list(paste("W",lp[,1],lp[,2],sep=sep["rho.name"]))
            ## strata relative to the parameter
            strata.rho <- strata
            ## provide a name for the parameter
            level.rho <- paste(paste0("(",nlme::collapse(lp2data[lp[,1],name.cov,drop=FALSE], sep = sep[1], as.factor = FALSE),")"),U.strata[strata], sep = sep[3])
            ## index of the corresponding observations
            indexObs.rho <- apply(index, MARGIN = 1, FUN = identity, simplify = FALSE)
        }        

    }else if(block == "B"){

        if(structure$twin){ ## same correlation for all blocks at the same nesting level
            diff.nesting <- (lp2data[lp[,1],name.cov,drop=FALSE]==lp2data[lp[,2],name.cov,drop=FALSE])
            diffS.nesting <- collapse(diff.nesting, sep = sep["lp"], as.factor = FALSE)
            index.rho <- paste(diffS.nesting, strata, sep = sep["rho.strata"])
        }else{ ## correlation specific to each pair of linear predictors
            index.rho <- paste(collapse(lp2data[lp[,1],,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                               collapse(lp2data[lp[,2],,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                               strata, sep = sep["rho.strata"])
        }
        
        ## internal identifier for the parameter        
        code.rho <- tapply(paste(block,lp[,1],lp[,2], sep=sep["rho.name"]), index.rho, FUN = identity, simplify = FALSE)
        ## strata
        strata.rho <- tapply(strata, index.rho, unique)
        ## index of the corresponding observations
        indexObs.rho <- tapply(1:NROW(index), index.rho, function(iVec){index[iVec,,drop=FALSE]})

        ## provide a name for the parameter
        if(structure$twin){
            ## order the covariates to find the ordering of the nesting
            diffE.nesting <- cbind(diff.nesting, matrix(TRUE, ncol = 1, nrow = NROW(diff.nesting), dimnames = list(NULL, name.cluster)))
            name.nesting <- names(sort(colSums(diffE.nesting), decreasing = FALSE)) ## outer nesting is when fewest pairs have equal levels
            ## sanity check
            lsDiff.nesting <- apply(cbind(FALSE, diffE.nesting[,name.nesting,drop=FALSE]),1,diff, simplify = FALSE)
            if(any(unlist(lsDiff.nesting)<0)){
                ## if one specifies id/baseline/hour and that some observations from the same subject have same baseline but different hours
                ## while others have same hours and different baseline
                ## hour should be changed to be relative to baseline, e.g. (no,0h) (no,1h) (yes,0h) (yes,1h) ----> (no,no.0h) (no,no.1h) (yes,no.0h) (yes,no.1h)
                warning("Variables \"",paste(name.nesting,collapse="\", \""),"\" do not define a proper nesting. \n")
            }
            ## naming
            level.rho <- paste0("(",sapply(lsDiff.nesting, function(iVec){paste(rev(names(iVec[iVec>0])),collapse="/")}),")")
            if(n.strata>1){
                level.rho <- paste(level.rho, U.strata[strata.rho], sep = sep["rho.strata"])
            }
        }else{ 
            level.rho <- paste0("(",
                                nlme::collapse(lp2data[lp[,1],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                                ",",
                                nlme::collapse(lp2data[lp[,2],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                                ")")
            name.nesting <- NULL
        }        
    }   

    ## ** constraint
    constraint.rho <- rep(NA, length(code.rho))

    ## ** export
    return(list(code = code.rho,
                level = level.rho,
                strata = unname(strata.rho),
                nesting = name.nesting,
                constraint = constraint.rho,
                indexObs = indexObs.rho))
}

## * .skeletonRho.RE
.skeletonRho.RE <- function(structure, XpairPattern, U.strata, name.cov, block, sep){

    ## ** handle special case (no repetition)
    if(all(lengths(index.cluster)==1)){return(structure)}

    ## ** identify variance and correlation parameters
    structure <- .skeletonRho.CS(structure = structure, XpairPattern = XpairPattern, U.strata = U.strata, name.cov = name.cov, block = block, sep = sep)

    ## ** relate correlation parameters to random effects
    n.strata <- length(U.strata)

    ## *** extract hierarchy
    hierarchy <- structure$ranef$hierarchy
    n.hierarchy <- length(hierarchy) 
    name.ranef <- unname(unlist(hierarchy))
    n.ranef <- length(name.ranef) 
    if(any(duplicated(name.ranef))){
        stop("Duplicated name(s) for the random effects: \"",paste(name.ranef[duplicated(name.ranef)], collapse = "\", \""),"\". \n",
             "Cannot associate correlation parameters and random effects. \n")
    }

    if(n.ranef==1){
        param.rho <- structure$param[structure$param$type=="rho" & is.na(structure$param$constraint),,drop=FALSE]
        structure$ranef$param <- list(matrix(param.rho$name[param.rho$index.strata], nrow = n.ranef, ncol = n.strata,
                                             dimnames = list(name.ranef,U.strata)))
    }else{

        ## *** prepare output
        structure$ranef$param <- lapply(hierarchy, function(iH){
            matrix(as.character(NA), nrow = length(iH), ncol = n.strata,
                   dimnames = list(iH,U.strata))
        })

        
        ## *** identify the variable constant outside the hierarchy and varying inside the hierarchy
        structure.activeRho <- structure$param[which(structure$param$type == "rho" & is.na(structure$param$constrain)),c("name","index.strata","code")]
        activeRho2var <- structure$rho$param2var[structure.activeRho$code,,drop=FALSE]
        code2name <- stats::setNames(structure.activeRho$name, structure.activeRho$code)
        var2col <- attr(structure$rho$param2var, "var2col")
        
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            iActiveRho2var <- activeRho2var[structure.activeRho$index.strata==iStrata,,drop=FALSE]

            for(iH in 1:n.hierarchy){ ## iH <- 1                
                iHierarchy <- hierarchy[[iH]]
                iIndex.equal <- which(rowSums(iActiveRho2var[,unlist(var2col[iHierarchy]),drop=FALSE])>0)
                if(length(iIndex.equal)!=length(iHierarchy)){
                    stop("Something when wrong when associating the correlation parameters and random effects. \n",
                         "hierarchy: ",paste(iHierarchy, collapse="\", \""),"\n",
                         "correlation parameter(s): \"",paste(rownames(iActiveRho2var)[iIndex.equal], collapse="\", \""),"\"\n")
                }

                ## with multiple variable inside a hierarchy, random effects are ordered according to how many variables vary within the hierarchy
                iOrder.ranef <- sort(rowSums(iActiveRho2var[iIndex.equal,,drop=FALSE]))
                structure$ranef$param[[iH]][,iStrata] <-  code2name[names(iOrder.ranef)]
            }
        }
    }

    ## ** export
    return(structure)
}

## * .skeletonRho.TOEPLITZ
.skeletonRho.TOEPLITZ <- function(structure, XpairPattern, U.strata, name.cov, block, sep){

    ## ** extract data
    lp <- XpairPattern[[paste0("lp",block)]]
    strata <- XpairPattern[[paste0("strata",block)]]
    index <- XpairPattern[[paste0("index",block)]]
    if(block=="B" && "lp2data.cross" %in% names(structure$cor)){
        lp2data <- structure$cor$lp2data.cross
    }else{ ## block W or same design matrix for both blocks
        lp2data <- structure$cor$lp2data
    }
    
    name.cor <- unlist(structure$name$cor)
    name.time <- structure$name$time
    name.strata <-  structure$name$strata
    name.cluster <- attr(structure$name$cluster,"original")
    name.nesting <- structure$nesting
    n.strata <- length(U.strata)
    
    ## ** identify parameters
    diff.time <- abs(lp2data[lp[,1],name.time]-lp2data[lp[,2],name.time])

    if(structure$twin){ ## shared correlation parameters for all blocks
        index.rho <- paste(diff.time, strata, sep = sep["rho.strata"])
    }else{ ## specific correlation parameters for each block
        index.rho <- paste(diff.time,
                           collapse(lp2data[lp[,1],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                           collapse(lp2data[lp[,2],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                           strata, sep = sep["rho.strata"])
    }
    ## internal identifier for the parameter
    code.rho <- tapply(paste(block,lp[,1],lp[,2], sep=sep["rho.name"]), index.rho, FUN = identity, simplify = FALSE)
    ## strata
    strata.rho <- tapply(strata, index.rho, unique)
    ## index of the corresponding observations
    indexObs.rho <- tapply(1:NROW(index), index.rho, function(iVec){index[iVec,,drop=FALSE]})

    ## provide a name for the parameter
    if(structure$twin){
        if(is.na(structure$class["correlation.cross"])){
            level.rho <- paste0("(",tapply(diff.time, index.rho, unique),")")
        }else{
            level.rho <- paste0(block,"(",tapply(diff.time, index.rho, unique),")")
        }
    }else{
        level.rho <- tapply(paste0("(",diff.time,"|",collapse(lp2data[lp[,1],name.cov,drop=FALSE], sep = sep["lp"]),",",collapse(lp2data[lp[,2],name.cov,drop=FALSE], sep = sep["lp"]),")"),
                            index.rho, unique)
    }
    if(n.strata>1){
        level.rho <- paste(level.rho, U.strata[strata.rho], sep = sep["rho.strata"])
    }

    ## ** constraint
    constraint.rho <- rep(NA, length(code.rho))

    ## ** export
    return(list(code = code.rho,
                level = level.rho,
                strata = unname(strata.rho),
                nesting = NULL,
                constraint = constraint.rho,
                indexObs = indexObs.rho))
}

## * .skeletonRho.UN
.skeletonRho.UN <- function(structure, XpairPattern, U.strata, name.cov, block, sep){

    ## ** extract data
    lp <- XpairPattern[[paste0("lp",block)]]
    strata <- XpairPattern[[paste0("strata",block)]]
    index <- XpairPattern[[paste0("index",block)]]
    if(block=="B" && "lp2data.cross" %in% names(structure$cor)){
        lp2data <- structure$cor$lp2data.cross
    }else{ ## block W or same design matrix for both blocks
        lp2data <- structure$cor$lp2data
    }
    
    name.cor <- unlist(structure$name$cor)
    name.time <- structure$name$time
    name.strata <-  structure$name$strata
    name.cluster <- attr(structure$name$cluster,"original")
    name.nesting <- structure$nesting
    n.strata <- length(U.strata)
    
    ## ** identify parameters
    if(length(name.cov)==0){
    }
    diff.time <- abs(lp2data[lp[,1],name.time]-lp2data[lp[,2],name.time])

    if(structure$twin){ ## shared correlation parameters for all blocks
        index.rho <- paste(diff.time, strata, sep = sep["rho.strata"])
    }else{ ## specific correlation parameters for each block
        index.rho <- paste(diff.time,
                           collapse(lp2data[lp[,1],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                           collapse(lp2data[lp[,2],name.cov,drop=FALSE], sep = sep["lp"], as.factor = FALSE),
                           strata, sep = sep["rho.strata"])
    }
    ## internal identifier for the parameter
    code.rho <- tapply(paste(block,lp[,1],lp[,2], sep=sep["rho.name"]), index.rho, FUN = identity, simplify = FALSE)
    ## strata
    strata.rho <- tapply(strata, index.rho, unique)
    ## index of the corresponding observations
    indexObs.rho <- tapply(1:NROW(index), index.rho, function(iVec){index[iVec,,drop=FALSE]})

    ## provide a name for the parameter
    if(structure$twin){
        if(is.na(structure$class["correlation.cross"])){
            level.rho <- paste0("(",tapply(diff.time, index.rho, unique),")")
        }else{
            level.rho <- paste0(block,"(",tapply(diff.time, index.rho, unique),")")
        }
    }else{
        level.rho <- tapply(paste0("(",diff.time,"|",collapse(lp2data[lp[,1],name.cov,drop=FALSE], sep = sep["lp"]),",",collapse(lp2data[lp[,2],name.cov,drop=FALSE], sep = sep["lp"]),")"),
                            index.rho, unique)
    }
    if(n.strata>1){
        level.rho <- paste(level.rho, U.strata[strata.rho], sep = sep["rho.strata"])
    }

    ## ** constraint
    constraint.rho <- rep(NA, length(code.rho))

    ## ** export
    return(list(code = code.rho,
                level = level.rho,
                strata = unname(strata.rho),
                nesting = NULL,
                constraint = constraint.rho,
                indexObs = indexObs.rho))

}

## * .skeletonRho.EXP
## .skeletonRho.EXP <- function(structure, data, 
##                            U.cluster, index.cluster,
##                            U.time, index.clusterTime, 
##                            U.strata, index.clusterStrata){

##     ## *** param rho
##     regressor <- colnames(X.cor)[which(attr(X.cor, "assign") == max(attr(X.cor, "assign")))]

##     if(n.strata==1){
##         param.rho <- "lambda"
##         strata.rho <- 1
##         code.rho <- regressor
##         level.rho <- ""
##         if(structure$heterogeneous){
##             param.rho <- c(param.rho,"nugget")
##             strata.rho <- c(strata.rho,1)
##             code.rho <- c(code.rho,NA)
##             level.rho <- c(code.rho,"")
##         }
##     }else{
##         param.rho <- paste0("lambda",U.strata)
##         strata.rho <- 1:n.strata
##         code.rho <- regressor
##         level.rho <- U.strata
##         if(structure$heterogeneous){
##             param.rho <- c(param.rho,paste0("nugget",U.strata))
##             strata.rho <- c(strata.rho,1:n.strata)
##             code.rho <- c(code.rho,rep(NA, n.strata))
##             level.rho <- c(level.rho,U.strata)
##         }
##     }

## }

## * helper
## ** .pairPatternX
##' @description Generate combinations of pairwise differences from linear predictors
##' 
##' @param object [list] output of model.matrix containing at least the design matrix in element \code{"X"}.
##' @param name.cov [character vector] covariate influencing the correlation, other than strata and time.
##' @param name.time [character] name of the time variable
##' @param U.cluster [character vector] cluster levels.
##' @param index.cluster [list of integer vector] position of each cluster in the dataset.
##' @param U.strata [character vector] strata levels.
##' @param indexCluster.strata [character vector] strata associated to each cluster.
##' @param sep [character] character used to collapse the covariate values into a linear predictor.
##' 
##' @noRd
##' @return When \code{blockwise} is TRUE, a list with the following elements \itemize{
##' \item lpW [numeric matrix, n.lp x 2] index of the linear predictor (diagonal blocks).
##' \item strataW [numeric vector, n.lp] index of the strata.
##' \item indexW [numeric matrix, n.lp x 2] index of the observations forming the pair in the dataset.
##' \item clusterW [numeric vector, n.lp] index of the cluster.
##' \item lpB [numeric matrix, n.lp x 2] index of the linear predictor (off-diagonal blocks).
##' \item diffXB [numeric matrix, n.lp n.X] difference in design matrix value in the  pair.
##' \item strataB [numeric vector, n.lp] index of the strata.
##' \item indexB [numeric matrix, n.lp x 2] index of the observations forming the pair in the dataset.
##' \item clusterB [numeric vector, n.lp] index of the cluster.
##' }
##' where n.lp is the number of distinct value that the linear predictor can take
##' and n.X the number of columns in the design matrix.
##' Otherwise, a list with the following elements \itemize{
##' \item lp [numeric matrix, n.lp x 2] index of the linear predictor (all blocks).
##' \item diffX [numeric matrix, n.lp] difference in design matrix value in the  pair.
##' \item strata [numeric vector, n.lp] index of the strata.
##' \item index [numeric matrix, n.lp x 2] index of the observations forming the pair in the dataset.
##' \item cluster [numeric vector, n.lp] index of the cluster.
##' }
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##'
##' X.cor <- model.matrix(~1, data = gastricbypassL)
##' .pairPatternX(list(X=X.cor),
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))),
##'               blockwise = TRUE)
##' 
##' X.cor <- model.matrix(~visit, data = gastricbypassL)
##' .pairPatternX(list(X=X.cor),
##'               U.cluster = unique(gastricbypassL$id),
##'               index.cluster = tapply(1:NROW(gastricbypassL), gastricbypassL$id, identity),
##'               U.strata = "1", index.clusterStrata = rep(1, length(unique(gastricbypassL$id))),
##'               blockwise = TRUE)
.pairPatternX <- function(object, name.cov, name.time,
                          U.cluster, index.cluster,
                          U.strata, index.clusterStrata,
                          sep = c("")){

    ## ** normalize user input
    n.strata <- length(U.strata)
    
    ## ** extract unique pairs of linear predictor 
    lp.pair <- list(within = NULL, between = NULL)
    index <- list(within = NULL, between = NULL)
    cluster <- list(within = NULL, between = NULL)
    strata <- list(within = NULL, between = NULL)

    if("X.cross" %in% names(object)){
        blocks <- c("within","between")
    }else{
        blocks <- "within"
    }
    
    for(iBlock in blocks){

        iSuffix <- switch(iBlock,
                          "within" = "",
                          "between" = ".cross")

        for(iS in 1:n.strata){ ## iS <- 1

            ## index of the clusters in the strata
            iIndex.clusterStrata <- which(index.clusterStrata==iS)
            ## index of the clusters with unique linear predictor pattern
            iIndex.clusterU <- iIndex.clusterStrata[.patternINtable12(pattern = object[[paste0("pattern",iSuffix)]][iIndex.clusterStrata],
                                                                      pattern2lp = object[[paste0("pattern2lp",iSuffix)]])
                                                    ]
            ## linear predictor for each cluster
            iLs.lp.clusterU <- object[[paste0("pattern2lp",iSuffix)]][object[[paste0("pattern",iSuffix)]][iIndex.clusterU]]

            ## remove triplicates
            iLs.lpU.clusterU <- lapply(1:length(iIndex.clusterU),function(iC){ ## iC <- 1
                iIndex <- which(!triplicated(iLs.lp.clusterU[[iC]]))
                iIndex2 <- iIndex[order(iLs.lp.clusterU[[iC]][iIndex])]
                iOut <- iLs.lp.clusterU[[iC]][iIndex2]
                attr(iOut,"index") <- index.cluster[[iIndex.clusterU[[iC]]]][iIndex2]
                return(iOut)
            })

            ## form all pairs
            iLS.pairs <- lapply(iLs.lpU.clusterU, unorderedPairs, distinct = TRUE)
            
            iNpair.lpU.clusterU <- sapply(iLS.pairs,NCOL)
            iM.pairs <- base::t(do.call(cbind,iLS.pairs))
            ## combine while remove pairs having the same linear predictors
            iTest.duplicated <- !duplicated(iM.pairs)
            ## retrieve cluster corresponding to the linear predictor
            iCluster <- do.call(c,mapply(x = iIndex.clusterU, y = iNpair.lpU.clusterU, function(x,y){rep(x, times = y)}, SIMPLIFY = FALSE))
            ## retrieve index of the observation corresponding to the linear predictor
            iIndex <- do.call(cbind,lapply(1:length(iIndex.clusterU), function(iC){
                iiIndex <- unorderedPairs(attr(iLs.lpU.clusterU[[iC]],"index"), distinct = TRUE)
                iiTest <- iTest.duplicated[cumsum(c(1,iNpair.lpU.clusterU))[iC]:cumsum(iNpair.lpU.clusterU)[iC]]
                return(iiIndex[,iiTest,drop=FALSE])
            }))

            ## update
            lp.pair[[iBlock]] <- rbind(lp.pair[[iBlock]], iM.pairs[iTest.duplicated,,drop=FALSE])
            index[[iBlock]] <- rbind(index[[iBlock]], t(iIndex))
            cluster[[iBlock]] <- c(cluster[[iBlock]], iCluster[iTest.duplicated])
            strata[[iBlock]] <- c(strata[[iBlock]],rep(iS, sum(iTest.duplicated)))
        }
    }

    ## ** index of the within and between block pairs + contrast in linear predictor
    if(length(name.cov)==0){ ## single block
        diffX.pair <- list(within = object$X[index$within[,2],,drop=FALSE] - object$X[index$within[,1],,drop=FALSE],
                           between = NULL)
        index.within <- 1:NROW(diffX.pair$within)
        index.between <- NULL
        
    }else{ ## block diagonal matrix
        if(length(blocks)==1){ ## duplicate in the case of same design matrix in diagonal and off diagonal blocks
            lp.pair$between <- lp.pair$within
            index$between <- index$within
            cluster$between <- cluster$within
            strata$between <- strata$within
            object$X.cross <- object$X
        }
        
        if(all(is.na(name.time)) || name.time %in%  attr(object$X,"variable") == FALSE){ ## time not in the design matrix
            indexX.cov <- 1:NCOL(object$X)
        }else{
            indexX.cov <- which(attr(object$X,"M.level")[,name.time] == attr(object$X,"reference.level")[,name.time]) ## exclude time when distinguishing within vs between block
        }
        index.within <- which(rowSums(object$X[index$within[,2],indexX.cov,drop=FALSE] - object$X[index$within[,1],indexX.cov,drop=FALSE]!=0)==0)
        index.between <- which(rowSums(object$X[index$between[,2],indexX.cov,drop=FALSE] - object$X[index$between[,1],indexX.cov,drop=FALSE]!=0)>0)
        
        diffX.pair <- list(within = object$X[index$within[index.within,2],,drop=FALSE] - object$X[index$within[index.within,1],,drop=FALSE],
                           between = object$X.cross[index$between[index.between,2],,drop=FALSE] - object$X.cross[index$between[index.between,1],,drop=FALSE])
        rownames(diffX.pair$within) <- NULL
        rownames(diffX.pair$between) <- NULL
    }

    ## ** export    
    return(list(lpW = lp.pair$within[index.within,,drop=FALSE],
                diffXW = diffX.pair$within,
                strataW = strata$within[index.within],
                indexW = index$within[index.within,,drop=FALSE],
                clusterW = cluster$within[index.within],
                lpB = lp.pair$between[index.between,,drop=FALSE],
                diffXB = diffX.pair$between,
                strataB = strata$between[index.between],
                indexB = index$between[index.between,,drop=FALSE],
                clusterB = cluster$between[index.between]
                ))
}

## ** .patternINtable12
##' @description Find the subset with unique pairs from a table
##'
##' @return column/cluster index for which unique pairs are observed
##' @noRd
##' @examples
##'
##' M1 <- cbind(c(1,0,0,0),c(1,1,1,1),c(1,0,1,0))
##' .patternINtable12(M1)
##' 
##' M2 <- cbind(c(1,0,0,0,0,0),c(1,1,1,0,0,0),c(0,0,0,1,1,1),c(0,1,0,0,1,1))
##' .patternINtable12(M2)
##' 
##' M3 <- cbind(c(1,0,0,0,0,0),c(1,1,1,0,0,0),c(0,0,0,1,1,1),c(0,2,0,0,1,1))
##' .patternINtable12(M3)
##' 
##' M3.bis <- cbind(c(0,0,0,2,1,1),c(1,0,0,0,0,0),c(1,1,1,0,0,0),c(0,0,0,1,1,1))
##' .patternINtable12(M3.bis)
##'
##' M4 <- cbind(c(0,0,0,1,1,1),c(2,0,0,0,0,0),c(1,1,1,0,0,0),c(0,0,0,1,1,1))
##' .patternINtable12(M4)
.patternINtable12 <- function(object, pattern = NULL, pattern2lp = NULL){
    
    ## *** initialization
    if(missing(object) && !is.null(pattern) && !is.null(pattern2lp)){

        index.unique <- which(!duplicated(pattern)) ## only keep one cluster per pattern
        Upattern <- pattern[index.unique] ## find unique patterns
        Upattern2lp <- pattern2lp[Upattern] ## find linear predictor per pattern
        Ulp <- unique(unlist(Upattern2lp)) ## find unique linear predictor

        tablePatternLp <- do.call(cbind,lapply(Upattern2lp, function(iVec){table(factor(iVec, levels = Ulp))})) ## table of linear predictors per patterns
        object <- pmin(tablePatternLp,2) ## observing more than twice a linear predictor does not lead to different parameter
                                         ## but observing twice may lead to \rho(lp1,lp1) instead of \rho(lp1,lp2)
    }else{
        index.unique <- 1:NCOL(object)
    }

    ## *** special case 1: single pattern
    n.pattern <- NCOL(object)
    if(n.pattern==1){
        return(1)
    }

    ## *** special case 2: single linear predictor
    n.lp <- NROW(object)
    if(n.lp==1){
        return(which.max(object[1,]))
    }


    ## *** find unique patterns
    index.current <- 1:n.pattern
    index.keep <- NULL
    ## max.row <- apply(object,1,max)
    ## max.current <- rep(0,n.lp)

    for(iP in 1:n.pattern){ ## iP <- 1
        iNewKeep <- index.current[which.max(colSums(object[,index.current,drop=FALSE]))]
        index.keep <- c(index.keep,iNewKeep)

        ## update remaining patterns
        index.current <- setdiff(index.current,iNewKeep)
        if(length(index.current)>0){
            ## also remove patterns nested in the current one
            iDiff <- sweep(object[,index.current,drop=FALSE], FUN = "-", MARGIN = 1, STATS = object[,iNewKeep])
            index.current <- index.current[which(colSums(iDiff>0)>0)]
        }
        if(length(index.current)==0){ ## no other pattern remaining
            break
        }
    }

    ## *** export
    return(index.unique[index.keep])
}


##----------------------------------------------------------------------
### structure-skeletonRho.R ends here
