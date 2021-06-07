### residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:40) 
## Version: 
## Last-Updated: Jun  7 2021 (17:01) 
##           By: Brice Ozenne
##     Update #: 128
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * residuals.lmm (documentation)
##' @title Extract The Residuals From a Multivariate Gaussian Model
##' @description Extract or compute the residuals of a multivariate gaussian model.
##' @name residuals
##' 
##' @param object a \code{lmm} object.
##' @param type.residual [character] Should the raw residuals be output (\code{"response"}), or the Pearson residuals (\code{"pearson"}),  or normalized residuals (\code{"normalized"} or \code{"tnormalized"}).
##' @param format [character] Should the residuals be output relative as a vector (\code{"long"}), or as a matrix with in row the clusters and in columns the outcomes (\code{"wide"}).
##' @param data [data.frame] dataset relative to which the residuals should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the residuals. Only relevant if differs from the fitted values.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \itemize{
##' \item \code{"raw"}: observed outcome minus fitted value.
##' \item \code{"pearson"}: each raw residual is divided by its modeled standard deviation.
##' \item \code{"normalized"}: raw residuals are multiplied, within clusters, by the inverse of the (lower) Cholesky factor.
##' \item \code{"tnormalized"}: raw residuals are multiplied, within clusters, by the inverse of the (upper) Cholesky factor.
##' }
##'
##' @return
##' When argument format is \code{"long"} and type.oobject is \code{"lmm"}, a vector containing the value of the residual realtive to each observation.
##' When argument format is \code{"wide"} and type.oobject is \code{"lmm"}, a data.frame with the value of the residual relative to each cluster (in rows) at each timepoint (in columns).
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Multivariate Gaussian Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' ## prediction
##' residuals(eUN.lmm, format = "wide")
##' residuals(eUN.lmm, format = "wide", type.residual = "normalized")


## * residuals.lmm (code)
##' @rdname residuals
##' @export
residuals.lmm <- function(object, type.residual = "response", format = "long",
                          data = NULL, p = NULL, type.object = "lmm", ...){
    options <- LMMstar.options()

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    type.object <- match.arg(type.object, c("lmm","gls"))
    format <- match.arg(format, c("wide","long"))
    type.residuals <- match.arg(type.residual, c("response","pearson","normalized","tnormalized"))

    ## ** extract
    if(type.object == "lmm"){

        if(!is.null(data)){
            ff.allvars <- c(all.vars(object$formula$mean), all.vars(object$formula$var))
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            design <- .model.matrix.lmm(formula.mean = object$formula$mean.design,
                                        formula.var = object$formula$var.design,
                                        data = data,
                                        var.outcome = object$outcome$var,
                                        var.strata = object$strata$var, U.strata = object$strata$levels,
                                        var.time = object$time$var, U.time = object$time$levels,
                                        var.cluster = object$cluster$var,
                                        structure = object$structure
                                        )
        }else{
            design <- object$design
        }
        Y <- design$Y
        X <- design$X.mean
        X.var <- design$X.var
        n.cluster <- design$cluster$n
        index.cluster <- design$index.cluster
        index.variance <- design$X.var$cluster
        index.time <- design$index.time

        if(!is.null(p)){
            if(any(duplicated(names(p)))){
                stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
            }
            if(any(names(object$param$type) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(object$param$type)[names(object$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            beta <- p[names(object$param$type=="mu")]
            if(type.residuals %in% c("pearson","normalized","tnormalized")){
                Omega <- .calc_Omega(object = X.var, param = p)
                precision <- lapply(Omega, solve)
            }
        }else{
            beta <- object$param$value[object$param$type=="mu"]
            if(type.residuals %in% c("pearson","normalized","tnormalized")){
                Omega <- object$Omega
                precision <- object$OmegaM1
            }
        }
        if(type.residuals == "pearson"){
            sqrtPrecision <- lapply(precision,function(iM){sqrt(diag(iM))})
        }else if(type.residuals == "normalized"){
            sqrtPrecision <- lapply(precision,function(iP){t(chol(iP))})
        }else if(type.residuals == "tnormalized"){
            sqrtPrecision <- lapply(precision,chol)
        }

        if(!is.null(data) || !is.null(p)){
            res <-  Y - X %*% beta
        }else{
            res <- object$residuals
        }

        ## normalization
        if(type.residuals %in% c("pearson","normalized","tnormalized")){
            for(iId in 1:n.cluster){ ## iId <- 7
                iIndex <- which(index.cluster==iId)
                iOrder <- order(index.time[iIndex])
                iResidual <- res[iIndex[iOrder],,drop=FALSE]
                if(type.residuals == "pearson"){
                    resnorm <- res[index.cluster==iId] * sqrtPrecision[[index.variance[iId]]]
                }else if(type.residuals %in% c("normalized","tnormalized")){
                    resnorm <- as.double(res[index.cluster==iId] %*% sqrtPrecision[[index.variance[iId]]])
                }
                res[iIndex] <- resnorm[order(iOrder)]
            }
        }

        ## add NA
        if(is.null(data) && length(object$index.na)>0){
            inflateNA <-  .addNA(index.na = object$index.na, design = design, time = object$time)
            res.save <- res
            res <- matrix(NA, nrow = inflateNA$n.allobs, ncol = 1)
            res[-object$index.na,] <- res.save

            level.cluster <- inflateNA$level.cluster
            level.time <- inflateNA$level.time
        }else{
            level.cluster <- object$design$cluster$levels[index.cluster]
            level.time <- object$time$levels[index.time]
        }

        ## 
        if(format=="wide"){
            res <- reshape2::dcast(data = data.frame(residuals = res, cluster = level.cluster, time = level.time),
                                   formula = cluster~time, value.var = "residuals")
        }else{
            res <- as.vector(res)
        }
        return(res)

    }else if(type.object == "gls"){
        if(!is.null(data)){
            stop("Cannot handle argument \'data\' when argument \'type.object\' is \"gls\". \n")
        }
        if(!is.null(p)){
            stop("Cannot handle argument \'p\' when argument \'type.object\' is \"gls\". \n")
        }
        if(object$strata$n == 1){
            return(stats::residuals(object$gls, type = type.residual))
        }else {
            return(lapply(object$gls, residuals, type = type.residual))
        }
    }

}

## * .addNA
.addNA <- function(index.na, design, time){

    attr.cluster <- attr(index.na,"cluster")
    attr.time <- attr(index.na,"time")

    ## ** if no missing value or missing information about cluster or time returns nothing
    if(length(index.na)==0 || any(is.na(attr.cluster)) || any(is.na(attr.time))){
        return(NULL)
    }

    ## ** identify all clusters
    allcluster <- sort(union(design$cluster$levels, attr.cluster))
    n.allcluster <- length(allcluster)

    n.allobs <- length(design$index.cluster)+length(index.na)
    
    ## ** identify missing clusters
    missing.cluster <- setdiff(attr.cluster, design$cluster$levels)

    ## ** create extended vector of observations
    level.cluster <- rep(NA, n.allobs)
    level.cluster[-index.na] <- design$cluster$levels[design$index.cluster]
    level.cluster[index.na] <- attr.cluster

    level.time <- rep(NA, n.allobs)
    level.time[-index.na] <- time$levels[design$index.time]
    level.time[index.na] <- time$levels[attr.time]

    ## ** export
    return(list(allcluster = allcluster,
                missing.cluster = missing.cluster,
                n.allcluster = length(allcluster),
                n.allobs = length(level.cluster),
                level.cluster = level.cluster,
                level.time = level.time))
    
}
##----------------------------------------------------------------------
### residuals.R ends here
