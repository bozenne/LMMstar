### sigma.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: May  9 2024 (15:20) 
##           By: Brice Ozenne
##     Update #: 764
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
##' @param simplify [logical] When there is only one variance-covariance matrix, output a matrix instead of a list of matrices.
##' @param ... Not used. For compatibility with the generic method.
##'
##'
##' @return A list where each element contains a residual variance-covariance matrix.
##' Can also be directly a matrix when argument is \code{simplify=TRUE} and there is a single residual variance-covariance matrix. 
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
sigma.lmm <- function(object, cluster = NULL, p = NULL, chol = FALSE, inverse = FALSE, simplify = TRUE, ...){

    ## ** extract from object
    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.level <- stats::setNames(object$design$param$level,param.name)
    param.sigma <- stats::setNames(object$design$param$sigma,param.name)
    param.k.x <- stats::setNames(object$design$param$k.x,param.name)
    param.k.y <- stats::setNames(object$design$param$k.y,param.name)

    outcome.var <- object$outcome$var

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
    Upattern <- object.structure$Upattern
    object.Omega <- object$Omega
    object.cluster <- object$cluster$levels
    object.cluster.num <- object$cluster$index
    object.index.cluster <- object$design$index.cluster
    object.index.clusterTime <- object$design$index.clusterTime
    
    ## find cluster index associated to each pattern
    pattern.index.cluster1 <- sapply(attr(object.structure$pattern,"list")[Upattern$name],"[",1)
    ## find time associated to each pattern
    pattern.Utime <- stats::setNames(lapply(object.index.clusterTime[pattern.index.cluster1], function(iIndex){U.time[iIndex]}),
                                    Upattern$name)

    ## ** normalize user imput

    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    ## p
    if(!is.null(p)){
        init <- .init_transform(p = p, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, 
                                x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                                table.param = object$design$param)
        theta <- init$p
    }else{
        theta <- object$param
    }

    ## cluster    
    if(!is.null(cluster)){
        test.clusterDF <- inherits(cluster, "data.frame")
        if(test.clusterDF){
            if(outcome.var %in% names(cluster) == FALSE){
                cluster[[outcome.var]] <- NA
            }
            newdesign <- stats::model.matrix(object, newdata = cluster, effect = "variance", simplify = FALSE, na.rm = FALSE)            
            cluster.num <- 1:length(newdesign$index.cluster)
            if(!is.null(attr(object$cluster$var,"original"))){
                cluster.level <- cluster[sapply(newdesign$index.cluster, "[", 1), attr(object$cluster$var,"original")]                
            }else{
                cluster.level <- cluster.num
            }
            cluster <- cluster.level
        }else{
            if(any(duplicated(cluster))){
                stop("Argument \'cluster\' should contain duplicates. \n")
            }
            newdesign <- NULL
            if(is.numeric(cluster)){ ## numeric matching the XXcluster.indexXX variable
                if(any(cluster %in% stats::na.omit(object.cluster.num) == FALSE)){ 
                    stop("When numeric, elements in argument \'cluster\' should index the clusters, i.e. be between 1 and ",max(object.cluster.num, na.rm = TRUE),". \n", sep = "")
                }
                cluster.num <- cluster
                cluster.level <- object.cluster[match(cluster,object.cluster.num)] 
            }else if(is.character(cluster)){
                if(any(cluster %in% object.cluster == FALSE)){
                    stop("When character, elements in argument \'cluster\' should refer to clusters used to fit the model \n", sep = "")
                }
                cluster.num <- match(cluster, object.cluster)
                cluster.level <- object.cluster[cluster.num]               
            }else{
                stop("Incorrect value for argument \'cluster\'. Should be a numeric vector or a character vector. \n")
            }
        }
    }else{
        newdesign <- NULL
        test.clusterDF <- FALSE
    }


    ## ** rebuild residual variance-covariance matrix
    if(is.null(cluster)){ ## representative covariance patterns

        vec.ntime <- tapply(Upattern$n.time,Upattern$index.strata, max)
        if(all(vec.ntime == n.time)){
            ## Retrieve covariance from existing pattern

            if(!is.null(p)){
                Omega <- .calc_Omega(object = object.structure,
                                     param = theta,
                                     keep.interim = FALSE)
            }else{
                Omega <- object.Omega
                for(iO in 1:length(Omega)){
                    attr(Omega[[iO]], "sd") <- NULL
                    attr(Omega[[iO]], "cor") <- NULL
                }
            }
            Omega.time <- pattern.Utime[names(Omega)]

        }else{
            ## find group pattern (ignoring the time variable)
            Upattern$group <- .nameUpatterns(object$design$vcov, xfactor = object$xfactor, ignore.time = TRUE, sep = LMMstar.options()$sep[c("Gpattern.var","Gpattern.level")])
            Ugroup <- unique(Upattern$group)
            cluster.var <- attr(object$design$vcov$name$cluster,"original")
            vcov.var <- unique(stats::na.omit(c(attr(object$design$vcov$name$time,"original"),
                                                object$design$vcov$name$strata,
                                                object$design$vcov$name$var[[1]],
                                                object$design$vcov$name$cor[[1]])))
            
            
            ## Create artifical clusters agregating patterns
            ## (typically usefull in presence of missing values, e.g. only observe time A-B or time A-C but not A-B-C together)
            df.fulltime <- do.call(rbind,lapply(Ugroup, function(iGroup){ ## iGroup <- Ugroup[1]
                iPattern <- Upattern[Upattern$group==iGroup,c("name","n.time")]
                iCluster <- unlist(attr(object$design$vcov$pattern,"list")[iPattern$name])
                iDF <- object$data[object$data$XXcluster.indexXX %in% iCluster,c("XXtime.indexXX",vcov.var),drop=FALSE]
                if(!is.na(cluster.var)){
                    iDF[[cluster.var]] <- iGroup
                }
                return(unique(iDF)[c(cluster.var,vcov.var)])
            }))
            ## update structure
            object$design <- stats::model.matrix(object, newdata = df.fulltime, effects = "variance", simplify = FALSE)
            ## evaluate residual varince covariance matrix
            Omega <- .calc_Omega(object = object$design$vcov,
                                 param = theta,
                                 keep.interim = FALSE)
            ## identify timepoints
            Upattern <- object$design$vcov$Upattern
            object.index.clusterTime <- object$design$index.clusterTime
            pattern.index.cluster1 <- sapply(attr(object$design$vcov$pattern,"list")[names(Omega)],"[",1)
            pattern.Utime <- stats::setNames(lapply(object.index.clusterTime[pattern.index.cluster1], function(iIndex){U.time[iIndex]}),
                                             Upattern$name)
            Omega.time <- pattern.Utime[names(Omega)]
        }

        ## add groups
        Upattern$group <- .nameUpatterns(object$design$vcov, xfactor = object$xfactor, ignore.time = FALSE,
                                         sep = LMMstar.options()$sep[c("Gpattern.var","Gpattern.level")])

    }else if(test.clusterDF){ ## for new clusters/times

        Omega <- .calc_Omega(object = newdesign$vcov,
                             param = theta,
                             keep.interim = FALSE)
        ## identify timepoints
        Unewpattern <- newdesign$vcov$Upattern
        newdesign.index.clusterTime <- newdesign$index.clusterTime
        newpattern.index.cluster1 <- sapply(attr(newdesign$vcov$pattern,"list")[names(Omega)],"[",1)
        newpattern.Utime <- stats::setNames(lapply(newdesign.index.clusterTime[newpattern.index.cluster1], function(iIndex){U.time[iIndex]}),
                                         Unewpattern$name)
        Omega.time <- newpattern.Utime[names(Omega)]
        ## identify index
        pattern.cluster <- newdesign$vcov$pattern
        index.cluster <- newdesign$index.cluster[cluster.num]
        index.clusterTime <- newdesign$index.clusterTime[cluster.num]

    }else{ ## for existing clusters and time

        if(!is.null(p)){
            Omega <- .calc_Omega(object = object.structure,
                                 param = theta,
                                 keep.interim = FALSE)
        }else{
            Omega <- object.Omega
            for(iO in 1:length(Omega)){
                attr(Omega[[iO]], "sd") <- NULL
                attr(Omega[[iO]], "cor") <- NULL
            }
        }
        Omega.time <- pattern.Utime[names(Omega)]
        pattern.cluster <- object.structure$pattern[cluster.num]
        
    }

    ## ** add time names    
    if(!is.null(Omega.time)){
        for(iP in 1:length(Omega)){ ## iP <- 1
            dimnames(Omega[[iP]]) <- list(Omega.time[[iP]], Omega.time[[iP]])
        }
    }
    
    ## ** inverse
    if(chol){
        Omega <- lapply(Omega, chol)
    }
    if(inverse){
        Omega <- lapply(Omega, solve)
    }

    ## ** subset residual variance-covariance matrix
    if(is.null(cluster)){ ## find unique covariance patterns

        vec.Upattern <- unlist(by(Upattern,Upattern$group,function(iDf){
            iDf$name[which.max(iDf$n.time)]
        }, simplify = FALSE))
        
        ## subset
        out <- Omega[vec.Upattern]
        if(!is.null(names(vec.Upattern))){
            names(out) <- names(vec.Upattern)
        }

        ## add possibly missing times
        for(iPattern  in 1:length(vec.Upattern)){
            if(!identical(colnames(out[[iPattern]]),U.time) || !identical(rownames(out[[iPattern]]),U.time)){
                iOmega.save <- out[[iPattern]]
                out[[iPattern]] <- matrix(NA, nrow = n.time, ncol = n.time, dimnames = list(U.time,U.time))
                out[[iPattern]][rownames(iOmega.save),colnames(iOmega.save)] <- iOmega.save
            }
        }

        ## prepare for export
        if(!simplify){
            attr(out,"pattern") <- vec.Upattern
        }

    }else{ ## cluster specific covariance patterns

        out <- stats::setNames(Omega[pattern.cluster],cluster)
        if(!simplify){
            attr(out, "pattern") <- stats::setNames(pattern.cluster, cluster.level)
        }

    }

    ## ** export
    if(!simplify){
        attr(out,"design") <- newdesign
    }else if(length(out)==1){
        out <- out[[1]]
    }
    return(out)
}

## * sigma.clmm
##' @export
sigma.clmm <- function(object, ...){

    object$Omega <- .calc_Omega(object$design$vcov, param = object$param, keep.interim = FALSE)
    out <- sigma.lmm(object, ...)
    return(out)

}



##----------------------------------------------------------------------
### sigma.R ends here
