### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: nov 13 2021 (16:35) 
##           By: Brice Ozenne
##     Update #: 471
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
##' @description Extract the variance-covariance matrix of the model coefficients of a linear mixed model.
##' @name vcov
##' 
##' @param object a \code{lmm} object.
##' @param effects [character] Should the variance-covariance matrix for all coefficients be output (\code{"all"}),
##' or only for coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only for coefficients relative to the variance structure (\code{"variance"}),
##' or only for coefficients relative to the correlation structure (\code{"correlation"}).
##' @param robust [logical] Should robust standard error (aka sandwich estimator) be output instead of the model-based standard errors. Not feasible for variance or correlation coefficients estimated by REML.
##' @param df [logical] Should degree of freedom, computed using Satterthwaite approximation, for the model parameters be output.
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
vcov.lmm <- function(object, effects = "mean", robust = FALSE, df = FALSE, strata = NULL, data = NULL, p = NULL,
                     type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    
    options <- LMMstar.options()

    ## ** normalize user imput
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(is.null(effects)){
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"

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

        if(is.null(data) && is.null(p) && test.notransform && (df == FALSE || !is.null(object$df)) && (robust == FALSE) && attr(object$information,"type.information")==type.information){
            keep.name <- stats::setNames(names(coef(object, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                                                     names(coef(object, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))    

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
            test.precompute <- !is.null(object$design$precompute.XX)
         
            if(!is.null(data)){
                ff.allvars <- c(all.vars(object$formula$mean), all.vars(object$formula$var))
                if(any(ff.allvars %in% names(data) == FALSE)){
                    stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
                }

                design <- .model.matrix.lmm(formula.mean = object$formula$mean.design,
                                            structure = object$design$vcov,
                                            data = data,
                                            U.strata = object$strata$levels,
                                            U.time = object$time$levels,
                                            stratify.mean = object$opt$name=="gls",
                                            precompute.moments = test.precompute)
            }else{
                design <- object$design
            }

            if(!is.null(p)){
                if(any(duplicated(names(p)))){
                    stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
                }
                if(any(names(object$param$type) %in% names(p) == FALSE)){
                    stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(object$param$type)[names(object$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
                }
                p <- p[names(object$param$value)]
            }else{
                p <- object$param$value
            }

            outMoments <- .moments.lmm(value = p, design = design, time = object$time, method.fit = object$method.fit, type.information = type.information,
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                       logLik = FALSE, score = FALSE, information = FALSE, vcov = TRUE, df = df, indiv = FALSE, effects = effects, robust = robust,
                                       trace = FALSE, precompute.moments = test.precompute, method.numDeriv = options$method.numDeriv, transform.names = transform.names)

            if("variance" %in% effects && transform.k %in% c("sd","var","logsd","logvar") && object$strata$n>1 && transform.names){
                ## re-order values when converting to sd with strata (avoid sd0:0 sd0:1 sd1:0 sd1:1 sd2:0 sd2:1 ...)
                out.name <- names(stats::coef(object, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = TRUE))
                vcov <- outMoments$vcov[out.name,out.name,drop=FALSE]
                if(df>0){
                    attr(vcov,"df") <- outMoments$df[out.name]
                }
                if(df>1){
                    attr(vcov,"dVcov") <- outMoments$dVcov[out.name,out.name,out.name,drop=FALSE]
                }
            }else{
                vcov <- outMoments$vcov
                if(df>0){
                    attr(vcov,"df") <- outMoments$df
                }
                if(df>1){
                    attr(vcov,"dVcov") <- outMoments$dVcov
                }
            }
        }
    return(vcov)    
}


##----------------------------------------------------------------------
### vcov.R ends here
