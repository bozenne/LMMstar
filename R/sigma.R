### sigma.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: May 28 2022 (16:34) 
##           By: Brice Ozenne
##     Update #: 386
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
    U.original.cluster <- object$design$cluster$levels.original
    n.cluster <- object$design$cluster$n
    name.cluster <- object$design$cluster$var
    index.cluster <- object$design$index.cluster
    index.clusterTime <- object$design$index.clusterTime
    
    n.strata <- object$strata$n
    U.strata <- object$strata$levels
    name.time <- object$time$var
    U.time <- object$time$levels

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
            cluster <- unique(stats::model.matrix(object, data = cluster, effect = "variance")[["XXindex.clusterXX"]])
        }else{
            if(any(duplicated(cluster))){
                stop("Argument \'cluster\' should contain duplicates. \n")
            }
            if(is.numeric(cluster)){
                if(any(cluster %in% 1:n.cluster == FALSE)){ 
                    stop("When numeric, elements in argument \'cluster\' should index the clusters, i.e. be between 1 and ",n.cluster,". \n", sep = "")
                }
                cluster <- cluster
                U.output <- U.cluster                
            }else if(is.character(cluster)){
                if(any(cluster %in% U.original.cluster == FALSE)){
                    stop("When character, elements in argument \'cluster\' should refer to clusters used to fit the model \n", sep = "")
                }
                cluster <- match(cluster, U.original.cluster)
                U.output <- U.original.cluster                
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
        
        out <- .getUVarCov(object, Omega = Omega)
        for(iO in 1:length(out)){ ## iO <- 6
            attr(out[[iO]],"time") <- NULL
            attr(out[[iO]],"sd") <- NULL
            attr(out[[iO]],"cor") <- NULL
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
        out.pattern <- pattern.cluster$pattern[match(cluster,pattern.cluster$index.cluster)]
        out <- stats::setNames(Omega[out.pattern], U.output[cluster])
        attr(out, "pattern") <- as.character(out.pattern)

        for(iO in 1:length(out)){ ## iO <- 6
            ## DO NOT USE
            ## dimnames(out[[iO]]) <- list(object$time$levels[attr(out[[iO]],"time")],object$time$levels[attr(out[[iO]],"time")])
            ## as this is incorrect with CS structure and missing data (indeed CS can be for time 1,2 but also work for 2,3. However the previous line would incorrectly label the times)
            if(!is.numeric(object$data[[name.time]])){
                iC <- cluster[[iO]]
                iO.sort <- index.cluster[[iC]]
                iO.time <- U.time[index.clusterTime[[iC]]]
                dimnames(out[[iO]]) <- list(iO.time,iO.time)
            }
            attr(out[[iO]],"sd") <- NULL
            attr(out[[iO]],"cor") <- NULL
        }
    }

    ## ** export
    if(is.list(out) && length(out)==1 && simplifies){
        out <- out[[1]]
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
    heterogeneous <- object$design$vcov$heterogeneous
    n.time <- object$time$n
    U.time <- object$time$level
    Upattern <- object$design$vcov$X$Upattern
    Xpattern.var <- object$design$vcov$X$Xpattern.var
    Xpattern.cor <- object$design$vcov$X$Xpattern.cor
    data <- object$data
    formula <- object$design$vcov$formula
    Vindex.cluster <- attr(object$design$index.cluster, "vectorwise")
    
    ## ** summary statistic of each pattern
    XCpattern.var <- lapply(Xpattern.var, function(iVar){as.character(interaction(as.data.frame(iVar),sep=sep[1]))})
    if(!is.null(object$design$vcov$X$Xpattern.cor)){
        XCpattern.cor <- lapply(Xpattern.cor, function(iCor){as.character(interaction(as.data.frame(iCor),sep=sep[1]))})
        XCpattern <- mapply(x = XCpattern.var[Upattern$var], y = XCpattern.cor[Upattern$cor], function(x,y){paste(x,y,sep=sep[2])})
    }else{
        XCpattern <- XCpattern.var[Upattern$var]
    }
    names(XCpattern) <- Upattern$name

    labels <- unique(unlist(XCpattern))
    table.Xpattern <- do.call(rbind,lapply(XCpattern, function(iPattern){table(factor(iPattern, levels = labels))}))
    pcObs.Xpattern <- rowSums(table.Xpattern)/n.time
    
    ## ** identify patterns with unique set of parameters
    all.param <- object$design$vcov$param$name
    table.param <- do.call(rbind,lapply(Upattern$param, function(iParam){table(factor(iParam, levels = all.param))}))
    
    keep.pattern <- rownames(table.param)[which.max(rowSums(table.param) + pcObs.Xpattern/10)[1]]
    iTable.param <- table.param[setdiff(rownames(table.param),keep.pattern),colSums(table.param[keep.pattern,,drop=FALSE])==0,drop=FALSE]
    iter.max <- NROW(Upattern)-1
    while(NCOL(iTable.param)>0 && any(1 %in% iTable.param) && iter.max>0){
        iMax <- which.max(rowSums(iTable.param) + pcObs.Xpattern[rownames(iTable.param)]/10)[1] ## favor patterns with more timepoints
        keep.pattern <- c(keep.pattern, rownames(iTable.param)[iMax])
        iTable.param <- table.param[setdiff(rownames(table.param),keep.pattern),colSums(table.param[keep.pattern,,drop=FALSE])==0,drop=FALSE]
        iter.max <- iter.max-1
    }
    table.Xpattern.keep <- table.Xpattern[rownames(table.Xpattern) %in% keep.pattern,,drop=FALSE]

    ## ** check whether other patterns are not nested in the set of unique patterns
    ## if so add them
    if(heterogeneous){
        possibleNested.pattern <- setdiff(Upattern$name, keep.pattern)
        for(iPattern in possibleNested.pattern){ ## iPattern <- possibleNested.pattern[1]
            iTable <- table.Xpattern[rownames(table.Xpattern) == iPattern,,drop=FALSE]

            iTest <- rowSums(sweep(table.Xpattern.keep, MARGIN = 2, FUN = "-", STATS = iTable)<0)
            if(all(iTest>0)){
                keep.pattern <- c(keep.pattern, iPattern)
                table.Xpattern.keep <- table.Xpattern[rownames(table.Xpattern) %in% keep.pattern,,drop=FALSE]
            }
        }
    }

    ## ** recover time
    if(!is.numeric(object$data[[object$time$var]])){
        Upattern.var <- Upattern$var[match(keep.pattern,Upattern$name)]
        iIndex.time <- lapply(Xpattern.var[Upattern.var], function(iX){object$design$index.clusterTime[[attr(iX,"index.cluster")[1]]]})
        out <- mapply(x = Omega[keep.pattern], y = iIndex.time, function(x,y){
            if(!is.null(attr(x,"time"))){
                dimnames(x) <- list(U.time[attr(x,"time")],U.time[attr(x,"time")])
            }else{
                dimnames(x) <- list(U.time[y],U.time[y])
            }
            dimnames(attr(x,"sd")) <- list(U.time[y],NULL)
            dimnames(attr(x,"cor")) <- list(U.time[y],U.time[y])
            return(x)
        }, SIMPLIFY = FALSE)
    }else{
        out <- Omega[keep.pattern]
    }
    ## ** rename patterns
    Upattern.strata <- unlist(Upattern[match(keep.pattern,Upattern$name),"index.strata"])

    all.cov <- union(all.vars(formula$var), all.vars(formula$cor))
    Ucov.cluster <- do.call(rbind,by(data[,all.cov,drop=FALSE],
                                     INDICES = Vindex.cluster,
                                     FUN = function(iData){apply(iData,MARGIN=2,FUN = function(iCol){iU <- unique(iCol); if(length(iU)==1){iU}else{NA}})}, simplify = FALSE))
    Ucov.pattern <- which(colSums(is.na(Ucov.cluster))==0)
    if(length(Ucov.pattern)>0){
        iIndex.cluster <- sapply(Upattern[match(names(out),Upattern$name),"index.cluster"],"[",1)
        covname.pattern <- as.character(interaction(Ucov.cluster[iIndex.cluster,Ucov.pattern,drop=FALSE], sep = sep[1], drop=TRUE))
        names(out) <- covname.pattern
    }

    ## ** export Omega
    return(out)
}


##----------------------------------------------------------------------
### sigma.R ends here
