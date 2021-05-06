### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: May  4 2021 (10:59) 
##           By: Brice Ozenne
##     Update #: 198
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
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef} function.
##'
##' @return A matrix with an attribute \code{"df"} when argument df is set to \code{TRUE}.
##' 

## * vcov.lmm (code)
##' @rdname vcov
##' @export
vcov.lmm <- function(object, effects = "all", df = FALSE, type.object = "lmm", strata = NULL, data = NULL, p = NULL,
                     type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    options <- LMMstar.options()
    x.transform.sigma <- object$reparametrize$transform.sigma
    x.transform.k <- object$reparametrize$transform.k
    x.transform.rho <- object$reparametrize$transform.rho

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

    if(is.null(type.information)){
        type.information <- options$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    if(is.null(transform.sigma)){
        transform.sigma <- options$transform.sigma
    }else{
        transform.sigma <- match.arg(transform.sigma, c("none","log","square","logsquare"))
    }

    if(is.null(transform.k)){
        transform.k <- options$transform.k
    }else{
        transform.k <- match.arg(transform.k, c("none","log","square","logsquare","sd","logsd","var","logvar"))
    }

    if(is.null(transform.rho)){
        transform.rho <- options$transform.rho
    }else{
        transform.rho <- match.arg(transform.rho, c("none","atanh","cov"))
    }
    test.notransform <- (transform.sigma==x.transform.sigma) && (transform.k==x.transform.k) && (transform.rho==x.transform.rho)

    ## ** extract or recompute variance covariance matrix
    if(type.object=="lmm"){
        keep.name <- setNames(names(coef(object, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.name = TRUE)),
                              names(coef(object, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.name = TRUE)))
        
        if(is.null(data) && is.null(p) && test.notransform){
            vcov <- object$vcov[keep.name,keep.name]
            if(df){
                attr(vcov,"df") <- object$df[keep.name]
            }
        }else{
            if(df){
                attr(type.information,"detail") <- TRUE
            }
            infoFull <- information(object, data = data, p = p, type.information = type.information,
                                    transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.name = TRUE)
            detail <- attr(infoFull,"detail")
            vcovFull <- solve(infoFull)
            vcov <- vcovFull[keep.name,keep.name]
            if(df){
                param <- object$param
                param$mu <- detail$param[param$type=="mu"]
                param$sigma <- detail$param[param$type=="sigma"]
                param$k <- detail$param[param$type=="k"]
                param$rho <- detail$param[param$type=="rho"]
                df <- .df(param = param, reparametrize = detail$reparametrize, Y = detail$Y, X.mean = detail$X.mean, X.var = detail$X.var,
                          index.variance = detail$index.variance, time.variance = detail$time.variance, index.cluster = detail$index.cluster, 
                          pair.meanvarcoef = detail$pair.meanvarcoef, pair.varcoef = detail$pair.varcoef,
                          REML = detail$REML, type.information = type.information,
                          type = object$param$type, strata = object$param$strata, transform.sigma = detail$transform.sigma, transform.k = detail$transform.k, transform.rho = detail$transform.rho, 
                          vcov = vcovFull, diag = TRUE)
                attr(vcov,"df") <- setNames(df[names(keep.name)],keep.name)
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
.df <- function(param, reparametrize, Y, X.mean, X.var,
                index.variance, time.variance, index.cluster, 
                pair.meanvarcoef, pair.varcoef, REML, type.information,
                type, strata, transform.sigma, transform.k, transform.rho, 
                vcov, diag){

    ## ** prepare vector of parameters
    param.type <- param$type
    param.namemu <- names(which(param.type=="mu"))
    param.namesigma <- names(which(param.type=="sigma"))
    param.namek <- names(which(param.type=="k"))
    param.namerho <- names(which(param.type=="rho"))
    param.nameVar <- c(param.namesigma,param.namek,param.namerho)
    n.param <- length(param.type)
    name.param <- names(param.type)
    test.transform <- (transform.sigma != "none") || (transform.k != "none") || (transform.rho != "none")

    ## param.value <- unlist(unname(param[c("mu","sigma","k","rho")]))
    param.trans.value <- c(param$mu,reparametrize$p)
    
    ## ** warper for computing information
    FUN_information <- function(p){

        sigma <- p[param.namesigma]
        k <- p[param.namek]
        rho <- p[param.namerho]

        if(test.transform){ ## back-transform
            backp <- .reparametrize(p = c(sigma,k,rho), type = type[param.nameVar], strata = strata[param.nameVar], 
                                    Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                    transform.sigma = transform.sigma,
                                    transform.k = transform.k,
                                    transform.rho = transform.rho,
                                    transform.names = FALSE)
            sigma <- backp$p[param.namesigma]
            k <- backp$p[param.namek]
            rho <- backp$p[param.namerho]
        }

        ## get jacobian
        newp <- .reparametrize(p = c(sigma,k,rho), type = type[param.nameVar], strata = strata[param.nameVar], 
                               Jacobian = TRUE, dJacobian = (REML || type.information == "observed"), inverse = FALSE,
                               transform.sigma = transform.sigma,
                               transform.k = transform.k,
                               transform.rho = transform.rho,
                               transform.names = FALSE)

        if(newp$transform){
            Jacobian <- newp$Jacobian
            dJacobian <- newp$Jacobian
        }else{
            Jacobian <- NULL
            dJacobian <- NULL
        }

        Omega <- .calc_Omega(object = X.var, sigma = sigma, k = k, rho = rho, keep.interim = TRUE)
        OmegaM1 <- lapply(Omega,solve)
        dOmega <- .calc_dOmega(object = X.var, sigma = sigma, k = k, rho = rho, Omega = Omega, Jacobian = Jacobian)

        if(type.information == "observed"){
            mu <- p[param.namemu]
            residuals <- Y - X.mean %*% mu
        }else{
            residuals <- NULL
        }
        if(REML || type.information == "observed"){
            d2Omega <- .calc_d2Omega(object = X.var, sigma = sigma, k = k, rho = rho, Omega = Omega, dOmega = dOmega, pair = pair.varcoef, Jacobian = Jacobian, dJacobian = dJacobian)
        }else{
            d2Omega <- NULL
        }

        info <- .information(X = X.mean, residuals = residuals, precision = OmegaM1, dOmega = dOmega, d2Omega = d2Omega,
                             index.variance = index.variance, time.variance = time.variance, index.cluster = index.cluster,
                             pair.meanvarcoef = pair.meanvarcoef, pair.varcoef = pair.varcoef, indiv = FALSE, REML = REML, type.information = type.information)
        return(as.double(info))
    }

    ## ** derivative of the information using numerical derivative
    ## matrix(FUN_information(param.trans.value), nrow = sqrt(n.param), ncol = sqrt(n.param))
    if(type.information == "observed"){
        M.dInfo <- numDeriv::jacobian(func = FUN_information, x = param.trans.value)
        colnames(M.dInfo) <- name.param
    }else{
        M.dInfo <- numDeriv::jacobian(func = function(p){FUN_information(c(param.trans.value[param.namemu],p))}, x = param.trans.value[c(param.namesigma,param.namek,param.namerho)])
        colnames(M.dInfo) <- c(param.namesigma,param.namek,param.namerho)
    }

    A.dVcov <- array(0, dim = rep(n.param,3), dimnames = list(name.param,name.param,name.param))
    for(iParam in 1:NCOL(M.dInfo)){
        iName <- colnames(M.dInfo)[iParam]
        A.dVcov[,,iName] <- M.dInfo[,iName]
        A.dVcov[,,iName] <- - vcov %*% A.dVcov[,,iName] %*% vcov
    }
    ## solve(crossprod(model.matrix(e.lmm, effects = "mean")))
    ## 4*coef(e.lmm)["sigma"]^2/nobs(e.lmm)[1]
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
