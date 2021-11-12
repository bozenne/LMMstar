### getVarCov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: nov 12 2021 (15:46) 
##           By: Brice Ozenne
##     Update #: 194
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
##' @param individual [character] identifier of the cluster for which to extract the residual variance-covariance matrix.
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
##' getVarCov(eUN.lmm)
##' getVarCov(eUN.lmm, individual = c("1","5"))

## * getVarCov.lmm
##' @rdname getVarCov
##' @export
getVarCov.lmm <- function(obj, individual = NULL, p = NULL, simplifies = TRUE, strata = NULL, ...){
    object <- obj

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.null(p) && any(names(which(object$param$type %in% c("sigma","k","rho"))) %in% names(p) == FALSE)){
        stop("Incorrect argument \'p\' - it should be a vector with names containing all variance and correlation parameters. \n")
    }
    
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }else{
        strata <- object$strata$levels
    }
    n.strata <- length(strata)

    if(!is.null(individual)){
        if(is.character(object$design$cluster$levels)){
            individual <- match.arg(as.character(individual), object$design$cluster$levels, several.ok = TRUE)
        }else if(any(individual %in% object$design$cluster$levels == FALSE) ){
            stop("Unknown values for argument \'individual\'. Should correspond to cluster id from the training dataset. \n")
        }
    }

    ## ** rebuild residual variance-covariance matrix
    if(!is.null(p)){
        Omega <- .calc_Omega(object = object$design$vcov,
                             param = p,
                             keep.interim = TRUE)
    }else{
        Omega <- object$Omega
    }

    ## ** subset
    if(is.null(individual)){
        if(object$strata$n==1){
            out <- stats::setNames(list(.getUVarCov(object, Omega = Omega)),object$strata$levels)
        }else{
            out <- stats::setNames(vector(mode = "list", length = n.strata),strata)
            for(iStrata in 1:n.strata){ ## iStrata <- 1
                out[[iStrata]] <- .getUVarCov(object, Omega = Omega[strata[object$design$vcov$X$Upattern$strata]==strata[iStrata]])
            }
        }
    }else{
        out <- Omega[stats::setNames(object$design$vcov$X$pattern.cluster,object$design$cluster$levels)[individual]]
        for(iO in 1:length(out)){
            dimnames(out[[iO]]) <- list(object$time$levels[attr(out[[iO]],"time")],object$time$levels[attr(out[[iO]],"time")])
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
        if(simplifies){
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
    }else{
        index.maxtime <- which(varPattern.ntime==max(varPattern.ntime))
        M.patterns <- matrix(0, nrow = time.n, ncol = length(index.maxtime),
                             dimnames = list(time.level, names(index.maxtime)))
        iMindex <- 1
        for(iPattern in index.maxtime){ ## iPattern <- 2
            M.patterns[index.time[[iPattern]],iMindex] <- 1
            iMindex <- iMindex+1
        }
        n.UvarPattern <- length(index.time)
        if(any(rowSums(M.patterns)==0)){ ## all covariance structure with max timepoints do not cover all times - add remaining structures one by one if they do not overlap previous timepoints
            for(iPattern in setdiff(1:n.UvarPattern,index.maxtime)){ ## iPattern <- 1
                if(all(rowSums(M.patterns)[index.time[[iPattern]]]==0)){
                    M.patterns <- cbind(M.patterns,NA)
                    M.patterns[index.time[[iPattern]],iMindex] <- 1
                    colnames(M.patterns)[iMindex] <- names(index.time)[iPattern]
                    iMindex <- iMindex+1
                }
            }
        }
        out <- lapply(Omega[colnames(M.patterns)], function(iO){
            dimnames(iO) <- list(time.level[attr(iO,"time")],time.level[attr(iO,"time")])
            attr(iO,"time") <- NULL
            attr(iO,"sd") <- NULL
            attr(iO,"cor") <- NULL
            return(iO)
        })
        ## rename each pattern according to the covariates
        out <- .patternName(Omega = out, data = object$data, name = object$design$vcov$name, Upattern = Upattern.strata, time.n = time.n)        

        ## reduce data to the clusters and variables involved in the unique patterns
        Upattern.strataRed <- Upattern.strata[Upattern.strata$name %in% names(out),,drop=FALSE]
        keep.col <- stats::na.omit(c(object$design$vcov$name$var[[1]],object$design$vcov$name$cor[[1]]))
        object$data$XXcluster.indexXX <- as.numeric(object$data$XXclusterXX)
        iData <- object$data[object$data$XXcluster.indexXX %in% Upattern.strataRed$example, c("XXcluster.indexXX","XXtime.indexXX",keep.col),drop=FALSE]
        

        dimnames(out) <- list(time.level[sapply(Omega[colnames(M.patterns)], attr, "time")],
                              time.level[sapply(Omega[colnames(M.patterns)], attr, "time")])
        attr(out,"pattern") <- colnames(M.patterns)

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
