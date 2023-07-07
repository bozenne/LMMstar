### sigma.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: jul  7 2023 (15:22) 
##           By: Brice Ozenne
##     Update #: 570
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
##' 
##' @param object a \code{lmm} object.
##' @param cluster [character, data.frame, NULL] identifier of the cluster(s) for which to extract the residual variance-covariance matrix.
##' For new clusters, a dataset containing the information (cluster, time, strata, ...) to be used to generate the residual variance-covariance matrices.
##' When \code{NULL}, will output complete data covariance patterns.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the residual variance-covariance matrix. Only relevant if differs from the fitted values.
##' @param chol [logical] Output the cholesky factorization of the variance-covariance matrix.
##' @param inverse [logical] Output the matrix inverse of the variance-covariance matrix.
##' @param simplifies [logical] When there is only one variance-covariance matrix, output a matrix instead of a list of matrices.
##' @param ... Not used. For compatibility with the generic method.
##'
##'
##' @return A list where each element contains a residual variance-covariance matrix.
##' Can also be directly a matrix when argument is \code{simplifies=TRUE} and there is a single residual variance-covariance matrix. 
##'
##' @keywords methods
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' dL$id.fac <- paste0("id",dL$id)
##' 
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id.fac,
##'                structure = "UN", data = dL, df = FALSE)
##'
##' ## extract residuals variance covariance matrix
##' sigma(eUN.lmm) ## unique patterns
##' sigma(eUN.lmm, cluster = c("id1","id5")) ## existing clusters
##' sigma(eUN.lmm, cluster = dL[1:7,,drop=FALSE]) ## new clusters

## * sigma.lmm
##' @export
sigma.lmm <- function(object, cluster = NULL, p = NULL, chol = FALSE, inverse = FALSE, simplifies = TRUE, ...){

    ## ** extract from object
    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.level <- stats::setNames(object$design$param$level,param.name)
    param.sigma <- stats::setNames(object$design$param$sigma,param.name)
    param.k.x <- stats::setNames(object$design$param$k.x,param.name)
    param.k.y <- stats::setNames(object$design$param$k.y,param.name)

    strata <- object$strata$levels
    n.strata <- length(strata)

    n.strata <- object$strata$n
    strata.var <- object$strata$var
    U.strata <- object$strata$levels
    if(!is.null(attr(object$time$var,"original"))){
        time.var <- attr(object$time$var,"original")
    }else{
        time.var <- object$time$var
    }
    U.time <- object$time$levels
    if(is.null(attr(U.time,"original"))){
        U.time.original <- U.time
    }else{
        U.time.original <- attr(U.time,"original")
    }
    n.time <- object$time$n
    if(!is.null(attr(object$cluster$var,"original"))){
        cluster.var <- attr(object$cluster$var,"original")
    }else{
        cluster.var <- object$cluster$var
    }
        
    object.structure <- object$design$vcov
    object.Omega <- object$Omega
    object.n.cluster <- object$design$cluster$n
    object.cluster <- object$design$cluster$levels
    object.index.cluster <- object$design$index.cluster
    object.index.clusterTime <- object$design$index.clusterTime
    
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
        test.clusterDF <- inherits(cluster, "data.frame")
        if(test.clusterDF){
            newdesign <- stats::model.matrix(object, data = cluster, effect = "variance", simplifies = FALSE)
            cluster <- 1:newdesign$cluster$n
            U.output <- newdesign$cluster$levels
        }else{
            if(any(duplicated(cluster))){
                stop("Argument \'cluster\' should contain duplicates. \n")
            }
            newdesign <- NULL
            if(is.numeric(cluster)){
                if(any(cluster %in% 1:object.n.cluster == FALSE)){ 
                    stop("When numeric, elements in argument \'cluster\' should index the clusters, i.e. be between 1 and ",object.n.cluster,". \n", sep = "")
                }
                cluster <- cluster
                U.output <- object.cluster[cluster] ## use object$design$cluster instead object$cluster to remove clusters with missing valuesb
            }else if(is.character(cluster)){
                if(any(cluster %in% object.cluster == FALSE)){
                    stop("When character, elements in argument \'cluster\' should refer to clusters used to fit the model \n", sep = "")
                }
                cluster <- match(cluster, object.cluster)
                U.output <- object.cluster[cluster]               
            }else{
                stop("Incorrect value for argument \'cluster\'. Should be a numeric vector or a character vector. \n")
            }
        }
    }else{
        newdesign <- NULL
        test.clusterDF <- FALSE
    }
    
    ## ** rebuild residual variance-covariance matrix
    if(test.clusterDF){ ## for new clusters/times
        if(is.null(p)){
            p <- coef(object, effects = "all")
        }
        Omega <- .calc_Omega(object = newdesign$vcov,
                             param = p,
                             keep.interim = FALSE)

        pattern.cluster <- as.character(newdesign$vcov$X$pattern.cluster$pattern[match(cluster,newdesign$vcov$X$pattern.cluster$index.cluster)])
        index.cluster <- newdesign$index.cluster[cluster]
        index.clusterTime <- newdesign$index.clusterTime[cluster]

    }else{ ## for existing clusters and time
        vec.cov <- c(setdiff(object.structure$name$var[[1]], c(time.var,object.structure$name$time)),
                     setdiff(object.structure$name$cor[[1]], c(time.var,object.structure$name$time)))
        vec.ntime <- tapply(object.structure$X$Upattern$n.time,unlist(object.structure$X$Upattern$index.strata), max)

        if(is.null(cluster) && any(vec.ntime < n.time) && length(stats::na.omit(vec.cov))==0){
            ## Agregate patterns in presence of missing values
            ## Only work when no covariate
            keep.index.strata <- which(vec.ntime < n.time)
            df.fulltime <- do.call(rbind,lapply(keep.index.strata, function(iStrata){
                if(!is.na(attr(strata.var,"original"))){
                    ## NOTE: use U.time.original instead of U.time in case multiple time variables
                    iDF <- data.frame(U.strata[iStrata],
                                      U.time.original, 
                                      object$cluster$levels[iStrata])
                    names(iDF) <- c(time.var,cluster.var,attr(strata.var,"original"))
                }else{
                    iDF <- data.frame(U.time.original,object$cluster$levels[iStrata])
                    names(iDF) <- c(time.var,cluster.var)
                }
                return(iDF)
            }))
            object$design$vcov <- model.matrix(object, data = df.fulltime, effects = "variance", simplifies = TRUE)
            Omega <- .calc_Omega(object = object$design$vcov,
                                 param = if(is.null(p)){object$param}else{p},
                                 keep.interim = FALSE)

        }else if(!is.null(p)){
            Omega <- .calc_Omega(object = object.structure,
                                 param = p,
                                 keep.interim = FALSE)
        }else{
            Omega <- object.Omega
            for(iO in 1:length(Omega)){
                attr(Omega[[iO]], "sd") <- NULL
                attr(Omega[[iO]], "cor") <- NULL
                attr(Omega[[iO]], "time") <- NULL
            }
        }
        
        pattern.cluster <- as.character(object.structure$X$pattern.cluster$pattern[match(cluster,object.structure$X$pattern.cluster$index.cluster)])
        index.cluster <- object.index.cluster[cluster]
        index.clusterTime <- object.index.clusterTime[cluster]
    }

    ## ** inverse
    if(chol){
        Omega <- lapply(Omega, chol)
    }
    if(inverse){
        Omega <- lapply(Omega, solve)
    }

    ## ** subset
    if(is.null(cluster)){ ## unique covariance patterns
        out <- .getUVarCov(object, Omega = Omega)
    }else{ ## cluster specific covariance patterns
        out <- stats::setNames(Omega[pattern.cluster],U.output)
        attr(out, "pattern") <- stats::setNames(pattern.cluster,U.output)
        attr(out,"index.cluster") <- stats::setNames(index.cluster,U.output)
        attr(out,"index.clusterTime") <- stats::setNames(index.clusterTime,U.output)

        for(iC in 1:length(out)){ ## iC <- 1
            ## DO NOT USE
            ## dimnames(out[[iO]]) <- list(object$time$levels[attr(out[[iO]],"time")],object$time$levels[attr(out[[iO]],"time")])
            ## as this is incorrect with CS structure and missing data (indeed CS can be for time 1,2 but also work for 2,3. However the previous line would incorrectly label the times)
            iC.time <- U.time[index.clusterTime[[iC]]]
            dimnames(out[[iC]]) <- list(iC.time,iC.time)
        }
        
    }

    ## ** export
    attr(out,"structure") <- newdesign
    if(is.list(out) && length(out)==1 && simplifies){
        out <- out[[1]]
        return(out)
    }else{
        if(simplifies){
            if(is.null(cluster)){
                for(iStrata in 1:length(out)){
                    attr(out[[iStrata]],"pattern") <- NULL
                }
            }else{
                attr(out,"index.clusterTime") <- NULL
                attr(out,"index.cluster") <- NULL
                attr(out,"pattern") <- NULL
            }
        }
        return(out)
    }
}

## * sigma.clmm
##' @export
sigma.clmm <- function(object, ...){

    object$Omega <- .calc_Omega(object$design$vcov, param = object$param, keep.interim = FALSE)
    out <- sigma.lmm(object, ...)
    return(out)

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
    if(length(object$index.na)){
        data <- data[-object$index.na,,drop=FALSE]
        data$XXcluster.indexXX <- attr(object$design$index.cluster,"vectorwise")
        data$XXtime.indexXX <- attr(object$design$index.clusterTime,"vectorwise")
        data$XXstrata.indexXX <- attr(object$design$index.clusterStrata,"vectorwise")
    }
    formula <- object$design$vcov$formula
    Vindex.cluster <- attr(object$design$index.cluster, "vectorwise")

    ## ** special case (no covariance structure)
    if(all(sapply(Upattern$param,length)==0)){
        index.patternMax <- which.max(Upattern$n.time)
        time.patternMax <- U.time[attr(Xpattern.var[[index.patternMax]],"index.time")]
        OmegaMax <- Omega[index.patternMax]
        dimnames(OmegaMax[[1]]) <- list(time.patternMax,time.patternMax)
        return(OmegaMax)
    }

    ## ** percentage of repetitions at each pattern
    pcObs.Xpattern <- stats::setNames(Upattern$n.time/n.time, Upattern$name)

    ## ** identify patterns with unique set of parameters
    all.param <- object$design$vcov$param$name
    table.param <- do.call(rbind,lapply(Upattern$param, function(iParam){table(factor(iParam, levels = all.param))}))
    rownames(table.param) <- Upattern$name

    keep.pattern <- rownames(table.param)[which.max(rowSums(table.param) + pcObs.Xpattern/10)[1]]
    iTable.param <- table.param[setdiff(rownames(table.param),keep.pattern),colSums(table.param[keep.pattern,,drop=FALSE])==0,drop=FALSE]
    iter.max <- NROW(Upattern)-1
    while(NCOL(iTable.param)>0 && any(1 %in% iTable.param) && iter.max>0){
        iMax <- which.max(rowSums(iTable.param) + pcObs.Xpattern[rownames(iTable.param)]/10)[1] ## favor patterns with more timepoints
        keep.pattern <- c(keep.pattern, rownames(iTable.param)[iMax])
        iTable.param <- table.param[setdiff(rownames(table.param),keep.pattern),colSums(table.param[keep.pattern,,drop=FALSE])==0,drop=FALSE]
        iter.max <- iter.max-1
    }

    ## ** recover time
    out <- mapply(x = Omega[keep.pattern], y = Upattern$index.time[match(keep.pattern,Upattern$name)], function(x,y){
        dimnames(x) <- list(U.time[y],U.time[y])
        return(x)
    }, SIMPLIFY = FALSE)

    ## ** rename patterns
    all.cov <- union(all.vars(formula$var), all.vars(formula$cor))
    if(length(all.cov)>0){
        Ucov.cluster <- do.call(rbind,by(data[,all.cov,drop=FALSE],
                                         INDICES = Vindex.cluster,
                                         FUN = function(iData){apply(iData,MARGIN=2,FUN = function(iCol){iU <- unique(iCol); if(length(iU)==1){iU}else{NA}})}, simplify = FALSE))
        Ucov.pattern <- which(colSums(is.na(Ucov.cluster))==0)
        if(length(Ucov.pattern)>0){
            iIndex.cluster <- sapply(Upattern[match(names(out),Upattern$name),"index.cluster"],"[",1)
            covname.pattern <- as.character(interaction(Ucov.cluster[iIndex.cluster,Ucov.pattern,drop=FALSE], sep = sep[1], drop=TRUE))
            names(out) <- covname.pattern
        }
    }

    ## ** export Omega
    return(out)
}


##----------------------------------------------------------------------
### sigma.R ends here
