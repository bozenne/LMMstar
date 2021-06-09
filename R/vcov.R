### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: Jun  9 2021 (11:43) 
##           By: Brice Ozenne
##     Update #: 401
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * vcov.lmm (documentation)
##' @title Extract The Variance-Covariance Matrix From a Multivariate Gaussian Model
##' @description Extract the variance-covariance matrix of the model coefficients of a multivariate gaussian model.
##' @name vcov
##' 
##' @param object a \code{lmm} object.
##' @param effects [character] Should the variance-covariance matrix for all coefficients be output (\code{"all"}),
##' or only for coefficients relative to the mean (\code{"mean"}),
##' or only for coefficients relative to the variance structure (\code{"variance"}),
##' or only for coefficients relative to the correlation structure (\code{"correlation"}).
##' @param robust [logical] Should robust standard error (aka sandwich estimator) be output instead of the model-based standard errors. Not feasible for variance or correlation coefficients estimated by REML.
##' @param df [logical] Should degree of freedom, computed using Satterthwaite approximation, for the model parameters be output.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param data [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the information. Only relevant if differs from the fitted values.
##' @param strata [character vector] When not \code{NULL}, only output the variance-covariance matrix for the estimated parameters relative to specific levels of the variable used to stratify the mean and covariance structure.
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
vcov.lmm <- function(object, effects = "all", robust = FALSE, df = FALSE, type.object = "lmm", strata = NULL, data = NULL, p = NULL,
                     type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    
    options <- LMMstar.options()

    ## ** normalize user imput
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)

    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }

    if(is.null(type.information)){
        type.information <- attr(object$information,"type.information")
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    ## ** extract or recompute variance covariance matrix
    if(type.object=="lmm"){

        keep.name <- stats::setNames(names(coef(object, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                              names(coef(object, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))

        if(is.null(data) && is.null(p) && test.notransform && (df == FALSE || !is.null(object$df)) && (robust == FALSE)){
            vcov <- object$vcov[keep.name,keep.name,drop=FALSE]
            if(transform.names){
                dimnames(vcov) <- list(names(keep.name),names(keep.name))
            }
            if(df>0){
                attr(vcov,"df") <- object$df[keep.name]
                if(transform.names){
                    names(attr(vcov,"df")) <- names(keep.name)
                }
            }
            if(df>1){
                attr(vcov,"dVcov") <- object$dVcov[keep.name,keep.name,keep.name,drop=FALSE]
                if(transform.names){
                    dimnames(attr(vcov,"dVcov")) <- list(names(keep.name),names(keep.name),names(keep.name))
                }
            }
        }else{
            attr(type.information,"robust") <- robust
            if(df>0){
                attr(type.information,"detail") <- TRUE
            }
            infoFull <- lava::information(object, data = data, p = p, type.information = type.information,
                                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
            detail <- attr(infoFull,"detail")
            vcovFull <- solve(infoFull)
            vcov <- vcovFull[keep.name,keep.name,drop=FALSE]
            if(transform.names){
                dimnames(vcov) <- list(names(keep.name), names(keep.name))
            }
            if(df>0){
                outdf <- .df(param = detail$param, reparametrize = detail$reparametrize, Y = detail$Y, X.mean = detail$X.mean, X.var = detail$X.var,
                             index.variance = detail$index.variance, time.variance = detail$time.variance, index.cluster = detail$index.cluster, name.varcoef = detail$name.varcoef,
                             time.k = object$design$param$time.k, time.rho = object$design$param$time.rho,
                             pair.meanvarcoef = detail$pair.meanvarcoef, pair.varcoef = detail$pair.varcoef,
                             REML = detail$REML, type.information = as.vector(type.information), effects = effects, 
                             transform.sigma = detail$transform.sigma, transform.k = detail$transform.k, transform.rho = detail$transform.rho, 
                             vcov = vcovFull, diag = TRUE, method.numDeriv = options$method.numDeriv, robust = robust)
                attr(vcov,"df") <- stats::setNames(outdf[keep.name],names(keep.name))
                if(df>1){
                    attr(vcov,"dVcov") <- attr(outdf,"dVcov")[keep.name,keep.name,keep.name]
                    dimnames(attr(vcov,"dVcov")) <- list(names(keep.name),names(keep.name),names(keep.name))
                }
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
            }else if("mean" %in% effects == FALSE){
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
                index.variance, time.variance, index.cluster, name.varcoef, 
                time.k, time.rho,
                pair.meanvarcoef, pair.varcoef, REML, type.information, effects, 
                transform.sigma, transform.k, transform.rho, 
                vcov, diag, method.numDeriv, robust){

    ## ** prepare vector of parameters
    param.type <- param$type
    param.strata <- param$strata
    name.allcoef <- names(param$value)
    n.allcoef <- length(param$type)

    param.nameVar <- name.allcoef[param$type %in% c("sigma","k","rho")]
    param.nameMean <- name.allcoef[param$type %in% c("mu")]

    test.transform <- (transform.sigma != "none") || (transform.k != "none") || (transform.rho != "none")

    param.trans.value <- c(param$value[param.nameMean],reparametrize$p)[name.allcoef]

    ## ** warper for computing information
    FUN_information <- function(p, as.double){

        paramVar <- p[param.nameVar]
        paramMean <- p[param.nameMean]

        if(test.transform){ ## back-transform
            backp <- .reparametrize(p = paramVar, type = param.type[param.nameVar], strata = param.strata[param.nameVar], 
                                    Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                    transform.sigma = transform.sigma,
                                    transform.k = transform.k,
                                    transform.rho = transform.rho,
                                    transform.names = FALSE)
            paramVar <- backp$p
        }

        ## get jacobian
        newp <- .reparametrize(p = paramVar, type = param.type[param.nameVar], strata = param.strata[param.nameVar], 
                               time.k = time.k, time.rho = time.rho,
                               Jacobian = TRUE, dJacobian = 2*(REML || type.information == "observed"), inverse = FALSE,
                               transform.sigma = transform.sigma,
                               transform.k = transform.k,
                               transform.rho = transform.rho,
                               transform.names = FALSE)

        if(newp$transform){
            Jacobian <- newp$Jacobian
            dJacobian <- newp$dJacobian
        }else{
            Jacobian <- NULL
            dJacobian <- NULL
        }

        Omega <- .calc_Omega(object = X.var, param = paramVar, keep.interim = TRUE)
        OmegaM1 <- lapply(Omega,solve)
        dOmega <- .calc_dOmega(object = X.var, param = paramVar, type = param.type[names(paramVar)], Omega = Omega, Jacobian = Jacobian)

        if(type.information == "observed"){
            residuals <- Y - X.mean %*% paramMean
        }else{
            residuals <- NULL
        }
        if(REML || type.information == "observed"){
            d2Omega <- .calc_d2Omega(object = X.var, param = paramVar, type = param.type[names(paramVar)], Omega = Omega, dOmega = dOmega, pair = pair.varcoef, Jacobian = Jacobian, dJacobian = dJacobian)
        }else{
            d2Omega <- NULL
        }
        info <- .information(X = X.mean, residuals = residuals, precision = OmegaM1, dOmega = dOmega, d2Omega = d2Omega, robust = robust,
                             index.variance = index.variance, time.variance = time.variance, index.cluster = index.cluster, name.varcoef = name.varcoef, name.allcoef = name.allcoef,
                             pair.meanvarcoef = pair.meanvarcoef, pair.varcoef = pair.varcoef, indiv = FALSE, REML = REML, type.information = type.information, effects = effects)

        if(as.double){
            return(as.double(info))
        }else{
            return(info)
        }
    }

    ## ** derivative of the information using numerical derivative
    ## matrix(FUN_information(param.trans.value, as.double = TRUE), nrow = n.allcoef, ncol = n.allcoef)
    name.effects <- colnames(FUN_information(param.trans.value, as.double = FALSE))
    if(type.information == "observed"){
        M.dInfo <- numDeriv::jacobian(func = function(p){FUN_information(p,as.double = TRUE)}, x = param.trans.value, method = method.numDeriv)
        colnames(M.dInfo) <- name.allcoef
    }else{
        M.dInfo <- numDeriv::jacobian(func = function(p){FUN_information(c(param$value[param.nameMean],p)[name.allcoef], as.double = TRUE)}, x = param.trans.value[param.nameVar], method = method.numDeriv)
        colnames(M.dInfo) <- param.nameVar
    }
    A.dVcov <- array(0, dim = rep(n.allcoef,3), dimnames = list(name.allcoef,name.allcoef,name.allcoef))
    for(iParam in 1:NCOL(M.dInfo)){ ## iParam <- 1
        iName <- colnames(M.dInfo)[iParam]
        A.dVcov[name.effects,name.effects,iName] <- M.dInfo[,iName]
        A.dVcov[,,iName] <- - vcov %*% A.dVcov[,,iName] %*% vcov
    }
    ## solve(crossprod(model.matrix(e.lmm, effects = "mean")))
    ## 4*coef(e.lmm)["sigma"]^2/stats::nobs(e.lmm)[1]
    ## ** degrees of freedom
    if(diag){
        df <- stats::setNames(sapply(1:n.allcoef, function(iP){
            2 * vcov[iP,iP]^2 / (A.dVcov[iP,iP,] %*% vcov %*% A.dVcov[iP,iP,])
        }), name.allcoef)
    }else{
        df <- matrix(NA, nrow = n.allcoef, ncol = n.allcoef, dimnames = list(name.allcoef, name.allcoef))
        for(iParam in 1:n.allcoef){
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
