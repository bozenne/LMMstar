### getVarCov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: Dec 15 2021 (17:37) 
##           By: Brice Ozenne
##     Update #: 275
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getVarCov.lmm (documentation)
##' @title Extract The Residuals Variance-Covariance Matrix From a Linear Mixed Model
##' @description Extract the unique set of residuals variance-covariance matrices or the one relative to specific clusters.
##' @name getVarCov
##' 
##' @param obj a \code{lmm} object.
##' @param individual [character, data.frame, NULL] identifier of the cluster(s) for which to extract the residual variance-covariance matrix.
##' For new clusters, a dataset containing the information (cluster, time, strata, ...) to be used to generate the residual variance-covariance matrices.
##' When \code{NULL}, will output complete data covariance patterns.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the residual variance-covariance matrix. Only relevant if differs from the fitted values.
##' @param strata [character vector] When not \code{NULL} and argument \code{individual} is not specified, only output the residual variance-covariance matrix relative to specific levels of the variable used to stratify the mean and covariance structure.
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
##' getVarCov(eUN.lmm) ## unique patterns
##' getVarCov(eUN.lmm, individual = c("1","5")) ## existing individuals
##' getVarCov(eUN.lmm, individual = dL[1:7,,drop=FALSE]) ## new individuals

## * getVarCov.lmm
##' @rdname getVarCov
##' @export
getVarCov.lmm <- function(obj, individual = NULL, p = NULL, simplifies = TRUE, strata = NULL, ...){
    object <- obj

    ## ** normalize user imput
    ## dot
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## p
    if(!is.null(p) && any(names(which(object$param$type %in% c("sigma","k","rho"))) %in% names(p) == FALSE)){
        stop("Incorrect argument \'p\' - it should be a vector with names containing all variance and correlation parameters. \n")
    }

    ## strata
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }else{
        strata <- object$strata$levels
    }
    n.strata <- length(strata)

    ## individual
    if(!is.null(individual)){
        if(inherits(individual, "data.frame")){
            object$design$vcov <- model.matrix(object, data = individual, effect = "variance")
        }else{
            if(any(duplicated(individual))){
                stop("Argument \'individual\' should contain duplicates. \n")
            }
            if(is.numeric(individual)){
                if(any(individual %in% 1:length(object$design$cluster$levels) == FALSE)){ ## use object$design$cluster instead object$cluster to remove clusters with missing values
                    stop("When numeric, elements in argument \'individual\' should index the clusters, i.e. be between 1 and ",object$design$cluster$n,". \n", sep = "")
                }
                individual <- object$design$cluster$levels[individual]
            }else if(is.character(object$design$cluster$levels)){
                if(any(individual %in% object$design$cluster$levels == FALSE)){
                    stop("When character, elements in argument \'individual\' should refer to clusters used to fit the model \n", sep = "")
                }
                individual <- match.arg(as.character(individual), object$design$cluster$levels, several.ok = TRUE)
            }else{
                stop("Incorrect value for argument \'individual\'. Should be a numeric vector or a character vector. \n")
            }
        }
    }

    ## ** rebuild residual variance-covariance matrix
    if(inherits(individual, "data.frame")){ ## for new individuals/times
        if(is.null(p)){
            p <- coef(object, effects = "all")
        }
        Omega <- .calc_Omega(object = object$design$vcov,
                             param = p,
                             keep.interim = TRUE)
    }else{ ## for existing individuals and time
        if(!is.null(p)){
            Omega <- .calc_Omega(object = object$design$vcov,
                                 param = p,
                                 keep.interim = TRUE)
        }else{
            Omega <- object$Omega
        }
    }

    ## ** subset
    if(is.null(individual)){ ## unique covariance patterns
        if(object$strata$n==1){            
            out <- stats::setNames(list(.getUVarCov(object, Omega = Omega)),object$strata$levels)
        }else{
            out <- stats::setNames(vector(mode = "list", length = n.strata),strata)
            for(iStrata in 1:n.strata){ ## iStrata <- 1
                out[[iStrata]] <- .getUVarCov(object, Omega = Omega[strata[object$design$vcov$X$Upattern$strata]==strata[iStrata]])
            }
        }
    }else if(inherits(individual,"data.frame")){ ## individual specific covariance patterns (new)
        
        Ucluster <- as.character(unique(individual[[object$cluster$var]]))
        out <- stats::setNames(Omega[object$design$vcov$X$pattern.cluster[Ucluster]],Ucluster)

        for(iO in 1:length(out)){ ## iO <- 6
            dimnames(out[[iO]]) <- list(object$time$levels[attr(out[[iO]],"time")],object$time$levels[attr(out[[iO]],"time")])
            attr(out[[iO]],"time") <- NULL
            attr(out[[iO]],"sd") <- NULL
            attr(out[[iO]],"cor") <- NULL
        }

    }else{ ## individual specific covariance patterns (existing)

        out <- stats::setNames(Omega[object$design$vcov$X$pattern.cluster[individual]], individual)

        for(iO in 1:length(out)){ ## iO <- 6
            ## DO NOT USE
            ## dimnames(out[[iO]]) <- list(object$time$levels[attr(out[[iO]],"time")],object$time$levels[attr(out[[iO]],"time")])
            ## as this is incorrect with CS structure and missing data (indeed CS can be for time 1,2 but also work for 2,3. However the previous line would incorrectly label the times)
            
            iO.sort <- attr(object$design$index.cluster,"sorted")[[which(object$design$cluster$levels == individual[iO])]]
            iO.time <- object$time$levels[object$design$index.time[iO.sort]]
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
        if(simplifies && is.null(individual)){
            for(iStrata in 1:length(out)){
                attr(out[[iStrata]],"pattern") <- NULL
            }
        }
        return(out)
    }
}

## * .getUVarPattern
## get residual variance covariance matrix at all timepoints
.getUVarCov <- function(object, Omega){
    ntime <- object$time$n
    time.level <- object$time$level
    time.n <- object$time$n
    Upattern <- object$design$vcov$X$Upattern
    Upattern.strata <- Upattern[names(Upattern$time) %in% names(Omega),]
    index.time <- Upattern.strata$time ## subset when strata
    varPattern.ntime <- sapply(index.time,length)

    index.maxPattern <- which(varPattern.ntime==ntime)
    if(length(index.maxPattern)==1){ ## exactly one covariance structure cover all times
        out <- Omega[[index.maxPattern]]
        dimnames(out) <- list(time.level[attr(out,"time")],time.level[attr(out,"time")])
        attr(out,"time") <- NULL
        attr(out,"sd") <- NULL
        attr(out,"cor") <- NULL
        attr(out,"pattern") <- Upattern.strata$name[index.maxPattern]
    }else if(sum(varPattern.ntime==ntime)>1){ ## more than one covariance structure cover all times
        out <- lapply(Omega[index.maxPattern], function(iO){
            dimnames(iO) <- list(time.level[attr(iO,"time")],time.level[attr(iO,"time")])
            attr(iO,"time") <- NULL
            attr(iO,"sd") <- NULL
            attr(iO,"cor") <- NULL
            return(iO)
        })
        ## rename each pattern according to the covariates
        out <- .patternName(Omega = out, data = object$data, name = object$design$vcov$name, Upattern = Upattern.strata, time.n = time.n)        
    }else{ ## no covariance structure covering all times
        out <- matrix(NA, nrow = time.n, ncol = time.n, dimnames = list(time.level, time.level))

        ## 1- get all covariance matrix
        ls.Omega <- getVarCov(object, individual = unlist(object$design$vcov$X$cluster.pattern[Upattern.strata$name]))
        
        ## 2- get all unique covariance matrix
        ls.UOmega <- ls.Omega[!duplicated(ls.Omega)]

        ## 3- assemble
        out <- matrix(NA, nrow = time.n, ncol = time.n, dimnames = list(time.level, time.level))
        warn <- FALSE
        for(iC in 1:length(ls.UOmega)){ ## iC <- 1
            diff <- stats::na.omit(out[rownames(ls.UOmega[[iC]]),colnames(ls.UOmega[[iC]])] - ls.UOmega[[iC]])
            if(length(diff)>0 && any(abs(diff)>1e-10)){warn <- TRUE}
            out[rownames(ls.UOmega[[iC]]),colnames(ls.UOmega[[iC]])] <- ls.UOmega[[iC]]
        }
        if(warn){
            warning("Issue when trying to extract unique covariance patterns: heterogeneity of the variance-covariance values within the same pattern \n",
                    "Contact the package manager. \n")
        }
        
    }
    return(out)
}


## * .patternName
## reduce data to the clusters and variables involved in the unique patterns
.patternName <- function(Omega, data, name, Upattern, time.n){
    Upattern.red <- Upattern[Upattern$name %in% names(Omega),,drop=FALSE]
    keep.col <- stats::na.omit(c(name$var[[1]],name$cor[[1]]))
    data$XXcluster.indexXX <- as.numeric(data$XXclusterXX)
    data <- data[data$XXcluster.indexXX %in% Upattern.red$example, c("XXcluster.indexXX","XXtime.indexXX",keep.col),drop=FALSE]
    ## remove column corresponding to the time variable
    if(time.n>1){
        test.time <- sapply(data[-(1:2)], function(iCol){sum(diag(table(iCol,data$XXtime.indexXX)))==NROW(data)})
        if(any(test.time)){
            data <- data[test.time==FALSE]
        }
    }

    ## for each subject, combine all combinations of variables
    if(NCOL(data)==0){
        OmegaSave <- Omega
        Omega <- as.matrix(Matrix::bdiag(OmegaSave))
        dimnames(Omega) <- list(unname(sapply(OmegaSave, rownames)),
                                unname(sapply(OmegaSave, colnames)))
        attr(Omega,"pattern") <- names(OmegaSave)
    }else{
        Omega.name <- tapply(1:NROW(data), data$XXcluster.index, function(iRow){paste(interaction(unique(data[-(1:2)][iRow,,drop=FALSE])), collapse = ", ")})
        names(Omega) <- unname(Omega.name)[order(Upattern.red$example)]
        attr(Omega,"pattern") <- Upattern.red$name
    }
    return(Omega)
}
##----------------------------------------------------------------------
### getVarCov.R ends here
