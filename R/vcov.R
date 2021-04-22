### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: Apr 22 2021 (18:09) 
##           By: Brice Ozenne
##     Update #: 118
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * vcov.lmm (documentation)
##' @title Extract The Variance-Covariance Matrix From a Linear Mixed Model
##' @description Extract the variance-covariance matrix of the model coefficients of a linear mixed model
##' @name vcov
##' 
##' @param object a \code{lmm} object.
##' @param effects [character] Should only variance/covariance relative to mean parameters (\code{"mean"})
##' or only relative to parameters for the variance-covariance structure (\code{"variance"}) be output, or both (\code{all}).
##' @param df [logical] Should degree of freedom, computed using Satterthwaite approximation, for the model parameters be output.
##' @param data [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the information. Only relevant if differs from the fitted values.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param transform [0,1,2] Transformation used on the variance coefficient. See details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \bold{transform}: \cr
##' \itemize{
##' \item 0 means no transformation i.e. ouput stanrdard error, ratio of standard errors, and correlations.
##' \item 1 means log/atanh transformation i.e. ouput log(stanrdard error), log(ratio of standard errors), and atanh(correlations).
##' \item 2 ouput variance coefficients and correlations.
##' }
##'
##' @return A matrix with an attribute \code{"df"} when argument df is set to \code{TRUE}.
##' 

## * vcov.lmm (code)
##' @rdname vcov
##' @export
vcov.lmm <- function(object, effects = "all", df = FALSE, type.object = "lmm", strata = NULL, data = NULL, p = NULL,
                     transform = NULL, transform.names = TRUE, type.information = NULL, ...){
    options <- LMMstar.options()
    object.transform <- object$transform

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    if(is.null(transform)){transform <- options$transform}
    if(is.null(type.information)){type.information <- options$type.information}

    ## ** extract or recompute variance covariance matrix
    if(type.object=="lmm"){
        keep.name <- names(coef(object, effects = effects, type.object = type.object, strata = strata))

        if(!is.null(data) || !is.null(p) || (transform != object.transform)){
            if(df){
                attr(type.information,"detail") <- TRUE
            }
            infoFull <- information(object, data = data, p = p, transform = transform, type.information = type.information, transform.name = FALSE)
            vcovFull <- solve(infoFull)
            vcov <- vcovFull[keep.name,keep.name]
            if(df){
                param <- object$param
                param$mu <- attr(infoFull,"detail")$param[param$type=="mu"]
                param$sigma <- attr(infoFull,"detail")$param[param$type=="sigma"]
                param$k <- attr(infoFull,"detail")$param[param$type=="k"]
                param$rho <- attr(infoFull,"detail")$param[param$type=="rho"]
                attr(vcov,"df") <- .df(param = param, Y = attr(infoFull,"detail")$Y, X.mean = attr(infoFull,"detail")$X.mean, X.var = attr(infoFull,"detail")$X.var,
                                       index.variance = attr(infoFull,"detail")$index.variance, time.variance = attr(infoFull,"detail")$time.variance, index.cluster = attr(infoFull,"detail")$index.cluster, 
                                       pair.meanvarcoef = attr(infoFull,"detail")$pair.meanvarcoef, pair.varcoef = attr(infoFull,"detail")$pair.varcoef,
                                       REML = attr(infoFull,"detail")$REML, type.information = type.information, transform = attr(infoFull,"detail")$transform, object.transform = object.transform,
                                       vcov = vcovFull, diag = TRUE)[keep.name]

            }

            if(transform>0 && transform.names){
                newnames <- names(coef(e.lmm, transform = transform, effects = effects, transform.names = transform.names))
                dimnames(vcov) <- list(newnames,newnames)
                names(attr(vcov,"df")) <- newnames
            }

        }else{
            vcov <- object$vcov[keep.name,keep.name]
            if(df){
                attr(vcov,"df") <- object$df[keep.name]
            }
        }

        return(vcov)

    }else if(type.object=="gls"){
        if(!is.null(data)){
            stop("Cannot handle argument \'data\' when argument \'type.object\' is \"gls\". \n")
        }
        if(!is.null(p)){
            stop("Cannot handle argument \'p\' when argument \'type.object\' is \"gls\". \n")
        }
        if(transform>0){
            stop("Cannot handle argument \'transform\' when argument \'type.object\' is \"gls\". \n")
        }
        if(df){
            stop("Cannot handle argument \'df\' when argument \'type.object\' is \"gls\". \n")
        }
        
        .getVcov <- function(oo, effects){
            if(identical(effects,"mean")){
                return(vcov(oo))
            }else if(identical(effects,"variance")){
                return(oo$apVar)
            }else{
                out.names <- c(colnames(vcov(oo)),colnames(oo$apVar))
                out <- as.matrix(Matrix::bdiag(vcov(oo), oo$apVar))
                dimnames(out) <- list(out.names,out.names)
                return(out)
            }
        }

        if(is.null(strata) && (object$strata$n == 1)){
            .getVcov(object$gls[[1]], effects = effects)
        }else{
            return(lapply(object$gls[which(object$strata$levels %in% strata)], .getVcov, effects = effects))
        }

    }
}

## * .dinformation
.df <- function(param, Y, X.mean, X.var,
                index.variance, time.variance, index.cluster, 
                pair.meanvarcoef, pair.varcoef, REML, type.information, transform, object.transform, vcov, diag){

    ## ** prepare vector of parameters
    param.type <- param$type
    param.namemu <- names(which(param.type=="mu"))
    param.namesigma <- names(which(param.type=="sigma"))
    param.namek <- names(which(param.type=="k"))
    param.namerho <- names(which(param.type=="rho"))
    param.value <- unlist(unname(param[c("mu","sigma","k","rho")]))
    n.param <- length(param.type)
    name.param <- names(param.type)
    
    ## ** warper for computing information
    FUN_information <- function(p){

        sigma <- p[param.namesigma]
        k <- p[param.namek]
        rho <- p[param.namerho]
        
        Omega <- .calc_Omega(object = X.var, sigma = sigma, k = k, rho = rho, keep.interim = TRUE)
        OmegaM1 <- lapply(Omega,solve)
        dOmega <- .calc_dOmega(object = X.var, sigma = sigma, k = k, rho = rho, Omega = Omega, transform = transform)

        if(type.information == "observed"){
            mu <- p[param.namemu]
            residuals <- Y - X.mean %*% mu
        }else{
            residuals <- NULL
        }
        if(REML || type.information == "observed"){
            d2Omega <- .calc_d2Omega(object = X.var, sigma = sigma, k = k, rho = rho, Omega = Omega, dOmega = dOmega, pair = pair.varcoef, transform = transform)
        }else{
            d2Omega <- NULL
        }

        info <- .information(X = X.mean, residuals = residuals, precision = OmegaM1, dOmega = dOmega, d2Omega = d2Omega,
                             index.variance = index.variance, time.variance = time.variance, index.cluster = index.cluster,
                             pair.meanvarcoef = pair.meanvarcoef, pair.varcoef = pair.varcoef, indiv = FALSE, REML = REML, type.information = type.information)
        return(as.double(info))
    }

    ## ** derivative of the information using numerical derivative
    ## FUN_information(param.value)
    if(transform == object.transform){
        if(type.information == "observed"){
            M.dInfo <- numDeriv::jacobian(func = FUN_information, x = param.value)
            colnames(M.dInfo) <- name.param
        }else{
            M.dInfo <- numDeriv::jacobian(func = function(p){FUN_information(c(param.value[param.namemu],p))}, x = param.value[c(param.namesigma,param.namek,param.namerho)])
            colnames(M.dInfo) <- c(param.namesigma,param.namek,param.namerho)
        }
    }else{
        ## transform
        param.valueT  <- c(param.value[param.namemu],
                           .coefVar(sigma = param.value[param.namesigma], k = param.value[param.namek], rho = param.value[param.namerho], transform = transform, transform.names = FALSE,
                                    param.type = param.type, param.strata = param$strata))
        FUN_informationT <- function(p){
            ## back-transform
            pT <- c(p[param.namemu],
                    .coefVar(sigma = p[param.namesigma], k = p[param.namek], rho = p[param.namerho], transform = -transform, transform.names = FALSE,
                             param.type = param.type, param.strata = param$strata))
            return(FUN_information(pT))
        }
        if(type.information == "observed"){
            M.dInfo <- numDeriv::jacobian(func = FUN_informationT, x = param.valueT)
            colnames(M.dInfo) <- name.param
        }else{
            M.dInfo <- numDeriv::jacobian(func = function(p){FUN_informationT(c(param.valueT[param.namemu],p))}, x = param.valueT[c(param.namesigma,param.namek,param.namerho)])
            colnames(M.dInfo) <- c(param.namesigma,param.namek,param.namerho)
        }
    }

    A.dVcov <- array(0, dim = rep(n.param,3), dimnames = list(name.param,name.param,name.param))
    for(iParam in 1:NCOL(M.dInfo)){
        iName <- colnames(M.dInfo)[iParam]
        A.dVcov[,,iName] <- M.dInfo[,iName]
        A.dVcov[,,iName] <- - vcov %*% A.dVcov[,,iName] %*% vcov
    }

    ## ** degrees of freedom
    if(diag){
        df <- setNames(sapply(1:n.param, function(iP){
            2 * vcov[iP,iP]^2 / (A.dVcov[iP,iP,] %*% vcov %*% A.dVcov[iP,iP,])
        }), name.param)
    }else{
        df <- matrix(NA, nrow = n.param, ncol = n.param, dimnames = list(name.param, name.param))
        for(iParam in 1:n.param){
            for(iiParam in 1:iParam){
                df[iParam,iiParam] <- 2 * vcov[iParam,iiParam]^2 / (A.dVcov[iParam,iiParam,] %*% vcov %*% A.dVcov[iiParam,iParam,])
                if(iParam != iiParam){
                    df[iiParam,iParam] <- df[iParam,iiParam]
                }
            }
        }
    }
    
    ## ** export
    attr(df,"dVcov") <- A.dVcov
    return(df)
}
##----------------------------------------------------------------------
### vcov.R ends here
