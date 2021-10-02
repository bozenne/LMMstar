### getVarCov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:57) 
## Version: 
## Last-Updated: okt  1 2021 (16:53) 
##           By: Brice Ozenne
##     Update #: 153
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
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
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
getVarCov.lmm <- function(obj, individual = NULL, p = NULL, type.object = c("lmm","gls"), simplifies = TRUE, strata = NULL, ...){
    object <- obj

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.null(p) && any(names(which(object$param$type %in% c("sigma","k","rho"))) %in% names(p) == FALSE)){
        stop("Incorrect argument \'p\' - it should be a vector with names containing all variance and correlation parameters. \n")
    }
    
    type.object <- match.arg(type.object, c("lmm","gls"))
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

    if(type.object == "lmm"){

        if(!is.null(p)){
            Omega <- .calc_Omega(object = object$design$vcov,
                                 param = p,
                                 keep.interim = TRUE)
        }else{
            Omega <- object$Omega
        }
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
        if(is.list(out) && length(out)==1 && simplifies){
            return(out[[1]])
        }else{
            return(out)
        }
        
    }else if(type.object == "gls"){
        if(object$strata$n==1){
            if(is.null(individual)){
                return(nlme::getVarCov(object$gls[[1]]))
            }else{
                return(nlme::getVarCov(object$gls[[1]], individual = individual))
            }
        }else{
            if(is.null(individual)){
                return(lapply(object$gls, nlme::getVarCov))
            }else{
                out <- stats::setNames(vector(mode = "list", length = length(individual)),individual)
                for(iStrata in 1:object$strata$n){ ## iStrata <- 1
                    iIndiv <- intersect(individual,names(object$design$index.strata[object$design$index.strata==iStrata]))
                    if(length(iIndiv)>0){
                        out[match(iIndiv,names(out))] <- nlme::getVarCov(object$gls[[iStrata]], individual = iIndiv)
                    }
                }
                return(out)
            }
        }

    }
        
}

## * .getUVarPattern
## get residual variance covariance matrix at all timepoints
.getUVarCov <- function(object, Omega){
    ntime <- object$time$n
    time.level <- object$time$level
    time.n <- object$time$n
    index.time <- object$design$vcov$X$Upattern$time[names(Omega)] ## subset when strata
    varPattern.ntime <- sapply(index.time,length)
    if(any(varPattern.ntime==ntime)){ ## one covariance structure cover all times
        out <- Omega[[which(varPattern.ntime==ntime)]]
        dimnames(out) <- list(time.level[attr(out,"time")],time.level[attr(out,"time")])
        attr(out,"time") <- NULL
        attr(out,"sd") <- NULL
        attr(out,"cor") <- NULL
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
            if(any(rowSums(M.patterns)==0)){
                stop("Something went wrong when trying to identify the residual variance-covariance matrix at all timepoints. \n")
            }

        }

        out <- as.matrix(Matrix::bdiag(Omega[colnames(M.patterns)]))
        dimnames(out) <- list(time.level[sapply(Omega[colnames(M.patterns)], attr, "time")],
                              time.level[sapply(Omega[colnames(M.patterns)], attr, "time")])

    }
    return(out)
}
##----------------------------------------------------------------------
### getVarCov.R ends here
