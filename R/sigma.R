### sigma.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: maj 20 2022 (18:07) 
##           By: Brice Ozenne
##     Update #: 313
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * sigma.lmm (documentation)
##' @title Extract The Residuals Variance-Covariance Matrix From a Linear Mixed Model
##' @description Extract the unique set of residuals variance-covariance matrices or the one relative to specific clusters.
##' @name sigma
##' 
##' @param object a \code{lmm} object.
##' @param cluster [character, data.frame, NULL] identifier of the cluster(s) for which to extract the residual variance-covariance matrix.
##' For new clusters, a dataset containing the information (cluster, time, strata, ...) to be used to generate the residual variance-covariance matrices.
##' When \code{NULL}, will output complete data covariance patterns.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the residual variance-covariance matrix. Only relevant if differs from the fitted values.
##' @param inverse [logical] Output the matrix inverse of the variance-covariance matrix.
##' @param simplifies [logical] When there is only one variance-covariance matrix, output a matrix instead of a list of matrices.
##' @param ... Not used. For compatibility with the generic method.
##'
##'
##' @return A list where each element contains a residual variance-covariance matrix.
##' Can also be directly a matrix when argument is \code{simplifies=TRUE} and there is a single residual variance-covariance matrix. 
##'
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' ## extract residuals variance covariance matrix
##' sigma(eUN.lmm) ## unique patterns
##' sigma(eUN.lmm, cluster = c("1","5")) ## existing clusters
##' sigma(eUN.lmm, cluster = dL[1:7,,drop=FALSE]) ## new clusters

## * sigma.lmm
##' @rdname sigma
##' @export
sigma.lmm <- function(object, cluster = NULL, p = NULL, inverse = FALSE, simplifies = TRUE, ...){

    ## ** extract from object
    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.level <- stats::setNames(object$design$param$level,param.name)
    param.sigma <- stats::setNames(object$design$param$sigma,param.name)
    param.k.x <- stats::setNames(object$design$param$k.x,param.name)
    param.k.y <- stats::setNames(object$design$param$k.y,param.name)

    strata <- object$strata$levels
    n.strata <- length(strata)

    U.cluster <- object$design$cluster$levels ## use object$design$cluster instead object$cluster to remove clusters with missing values
    n.cluster <- object$design$cluster$n
    name.cluster <- object$design$cluster$var
    index.cluster <- object$design$index.cluster
    index.clusterTime <- object$design$index.clusterTime
    
    n.strata <- object$strata$n
    U.strata <- object$strata$levels

    structure <- object$design$vcov
    pattern.cluster <- structure$X$pattern.cluster
    Upattern <- structure$X$Upattern
    object.Omega <- object$Omega

    ## ** normalize user imput
    ## dot
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## p
    if(!is.null(p) && any(names(which(param.type %in% c("sigma","k","rho"))) %in% names(p) == FALSE)){
        stop("Incorrect argument \'p\' - it should be a vector with names containing all variance and correlation parameters. \n")
    }

    ## cluster
    if(!is.null(cluster)){
        if(inherits(cluster, "data.frame")){
            structure <- stats::model.matrix(object, data = cluster, effect = "variance")
        }else{
            if(any(duplicated(cluster))){
                stop("Argument \'cluster\' should contain duplicates. \n")
            }
            if(is.numeric(cluster)){
                if(any(cluster %in% 1:n.cluster == FALSE)){ 
                    stop("When numeric, elements in argument \'cluster\' should index the clusters, i.e. be between 1 and ",n.cluster,". \n", sep = "")
                }
                cluster <- U.cluster[cluster]
            }else if(is.character(U.cluster)){
                if(any(cluster %in% U.cluster == FALSE)){
                    stop("When character, elements in argument \'cluster\' should refer to clusters used to fit the model \n", sep = "")
                }
                cluster <- match.arg(as.character(cluster), U.cluster, several.ok = TRUE)
            }else{
                stop("Incorrect value for argument \'cluster\'. Should be a numeric vector or a character vector. \n")
            }
        }
    }

    ## ** rebuild residual variance-covariance matrix
    if(inherits(cluster, "data.frame")){ ## for new clusters/times
        if(is.null(p)){
            p <- coef(object, effects = "all")
        }
        Omega <- .calc_Omega(object = structure,
                             param = p,
                             keep.interim = TRUE)
    }else{ ## for existing clusters and time
        if(!is.null(p)){
            Omega <- .calc_Omega(object = structure,
                                 param = p,
                                 keep.interim = TRUE)
        }else{
            Omega <- object.Omega
        }
    }

    ## ** inverse
    if(inverse){
        Omega <- lapply(Omega, function(iO){
            iOut <- solve(iO)
            attr(iOut,"time") <- attr(iO,"time")
            attr(iOut,"sd") <- attr(iO,"sd")
            attr(iOut,"cor") <- attr(iO,"cor")
            return(iOut)
        })
    }

    ## ** subset
    if(is.null(cluster)){ ## unique covariance patterns
        if(n.strata==1){            
            out <- stats::setNames(list(.getUVarCov(object, Omega = Omega)),U.strata)
        }else{
            out <- stats::setNames(vector(mode = "list", length = n.strata),strata)
            for(iStrata in 1:n.strata){ ## iStrata <- 1
                out[[iStrata]] <- .getUVarCov(object, Omega = Omega[strata[Upattern$strata]==strata[iStrata]])
            }
        }
    }else if(inherits(cluster,"data.frame")){ ## cluster specific covariance patterns (new)
        
        Ucluster <- as.character(unique(cluster[[name.cluster]]))
        out <- stats::setNames(Omega[pattern.cluster[Ucluster]],Ucluster)

        for(iO in 1:length(out)){ ## iO <- 6
            dimnames(out[[iO]]) <- list(U.time[attr(out[[iO]],"time")],U.time[attr(out[[iO]],"time")])
            attr(out[[iO]],"time") <- NULL
            attr(out[[iO]],"sd") <- NULL
            attr(out[[iO]],"cor") <- NULL
        }

    }else{ ## cluster specific covariance patterns (existing)

        out <- stats::setNames(Omega[pattern.cluster[cluster]], cluster)

        for(iO in 1:length(out)){ ## iO <- 6
            ## DO NOT USE
            ## dimnames(out[[iO]]) <- list(object$time$levels[attr(out[[iO]],"time")],object$time$levels[attr(out[[iO]],"time")])
            ## as this is incorrect with CS structure and missing data (indeed CS can be for time 1,2 but also work for 2,3. However the previous line would incorrectly label the times)

            iC <- which(U.cluster == cluster[iO])
            iO.sort <- index.cluster[[iC]]
            iO.time <- U.time[index.clusteTime[[iC]]]
            dimnames(out[[iO]]) <- list(iO.time,iO.time)
            
            attr(out[[iO]],"time") <- NULL
            attr(out[[iO]],"sd") <- NULL
            attr(out[[iO]],"cor") <- NULL
        }

    }

    ## ** export
    if(is.list(out) && length(out)==1 && simplifies){
        out <- out[[1]]
        attr(out,"pattern") <- NULL
        return(out)
    }else{
        if(simplifies && is.null(cluster)){
            for(iStrata in 1:length(out)){
                attr(out[[iStrata]],"pattern") <- NULL
            }
        }
        return(out)
    }
}

## * getVarCov.lmm
##' @title Depreciated Extractor of the Residual Variance-Covariance Matrix
##' @description Depreciated extractor of the residual variance-covariance matrix.
##' @param obj  a \code{lmm} object.
##' @param ... other arguments.
##' @return Nothing
##' @seealso \code{\link{sigma.lmm}}
#' @export
getVarCov.lmm <- function(obj, ...) {
  .Deprecated("sigma.lmm")
  NULL
}

## * .getUVarPattern
##' @title Find unique variance-covariance patterns 
##' @noRd
.getUVarCov <- function(object, Omega, sep = c(":","::X::")){

    ## ** extract from object
    n.time <- object$time$n
    U.time <- object$time$level
    Upattern <- object$design$vcov$X$Upattern
    Xpattern.var <- object$design$vcov$X$Xpattern.var
    Xpattern.cor <- object$design$vcov$X$Xpattern.cor
    
    ## ** restrict to current strata
    iUpattern <- Upattern[Upattern$name %in% names(Omega),]
    iXpattern.var <- Xpattern.var[iUpattern$var]
    if(!is.null(Xpattern.cor)){
        iXpattern.cor <- Xpattern.cor[iUpattern$cor]
    }

    ## ** summary statistic of each pattern
    iXCpattern.var <- lapply(iXpattern.var, function(iVar){as.character(interaction(as.data.frame(iVar),sep=sep[1]))})
    if(!is.null(object$design$vcov$X$Xpattern.cor)){
        iXCpattern.cor <- lapply(iXpattern.cor, function(iCor){as.character(interaction(as.data.frame(iCor),sep=sep[1]))})
        iXCpattern <- mapply(x = iXCpattern.var, y = iXCpattern.cor, function(x,y){paste(x,y,sep=sep[2])})
        names(iXCpattern) <- iUpattern$name
    }else{
        iXCpattern <- iXCpattern.var
    }
    labels <- unique(unlist(iXCpattern))
    table.Xpattern <- do.call(rbind,lapply(iXCpattern, function(iPattern){table(factor(iPattern, levels = labels))}))
    table.Xpattern <- table.Xpattern[order(rowSums(table.Xpattern), decreasing = TRUE),,drop=FALSE]

    ## ** identify patterns with unique set of parameters
    iUpattern.sort <- iUpattern[match(rownames(table.Xpattern), iUpattern$name),]
   
    keep.pattern <- iUpattern.sort$name[!duplicated(iUpattern.sort$param)]
    table.Xpattern.keep <- table.Xpattern[rownames(table.Xpattern) %in% keep.pattern,,drop=FALSE]

    ## ** check whether other patterns are not nested in the set of unique patterns
    ## if so add them
    possibleNested.pattern <- setdiff(Upattern$name, keep.pattern)
    for(iPattern in possibleNested.pattern){ ## iPattern <- possibleNested.pattern[1]
       iTable <- table.Xpattern[rownames(table.Xpattern) == iPattern,,drop=FALSE]

       iTest <- rowSums(sweep(table.Xpattern.keep, MARGIN = 1, FUN = "-", STATS = iTable)<0)
       if(any(iTable>0)){
           keep.pattern <- c(keep.pattern, iPattern)
           table.Xpattern.keep <- table.Xpattern[rownames(table.Xpattern) %in% keep.pattern,,drop=FALSE]
       }
    }

    ## ** export Omega
    out <- Omega[[keep.pattern]]

    return(out)
}


##----------------------------------------------------------------------
### sigma.R ends here
