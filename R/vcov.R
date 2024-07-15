### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: jul 15 2024 (16:49) 
##           By: Brice Ozenne
##     Update #: 658
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
##' 
##' @param object a \code{lmm} object.
##' @param effects [character] Should the variance-covariance matrix for all coefficients be output (\code{"all"}),
##' or only for coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only for coefficients relative to the variance structure (\code{"variance"}),
##' or only for coefficients relative to the correlation structure (\code{"correlation"}).
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors.
##' Not feasible for variance or correlation coefficients estimated by REML.
##' @param df [logical] Should degree of freedom, computed using Satterthwaite approximation, for the model parameters be output.
##' @param newdata [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the variance-covariance matrix. Only relevant if differs from the fitted values.
##' @param strata [character vector] When not \code{NULL}, only output the variance-covariance matrix for the estimated parameters relative to specific levels of the variable used to stratify the mean and covariance structure.
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef.lmm} function.
##'
##' @return A matrix with an attribute \code{"df"} when argument df is set to \code{TRUE}.
##'
##' @keywords methods 

## * vcov.lmm (code)
##' @export
vcov.lmm <- function(object, effects = "mean", robust = FALSE, df = FALSE, strata = NULL,
                     newdata = NULL, p = NULL,
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
        type.information <- object$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }
    if(df && robust && object$args$method.fit == "REML"){
        stop("Cannot compute degrees of freedom under REML for robust standard errors. \n",
             "Consider setting the argument df to FALSE",
             " \n or using ML estimation by setting the argument method.fit=\"ML\" when calling lmm.")
    }

    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho,
                            table.param = object$design$param)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    if(is.null(p)){
        theta <- object$param
    }else{
        theta <- init$p
    }    

    ## ** extract or recompute variance covariance matrix
    if(is.null(newdata) && is.null(p) && test.notransform && (df == FALSE || !is.null(object$df)) && (robust == FALSE) && object$args$type.information==type.information){
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
         
        if(!is.null(newdata)){
            design <- stats::model.matrix(object, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- object$design
        }

        outMoments <- .moments.lmm(value = theta, design = design, time = object$time, method.fit = object$args$method.fit, type.information = type.information,
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

## * vcov.mlmm
##' @export
vcov.Wald_lmm <- function(object, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    return(object$vcov)
    
}

## * vcov.mlmm (documentation)
##' @title Extract The Variance-Covariance Matrix From Multiple Linear Mixed Models
##' @description Extract the variance-covariance matrix of the model coefficients from multiple linear mixed models.
##'
##' @param object a \code{mlmm} object.
##' @param effects [character] By default will output the estimates relative to the hypotheses being tested (\code{"contrast"}).
##' But can also output all model coefficients (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' @param p [list of numeric vector] list of model coefficients to be used. Only relevant if differs from the fitted values.
##' @param newdata [NULL] Not used. For compatibility with the generic method.
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors.
##' Not feasible for variance or correlation coefficients estimated by REML.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param simplify [logical] Should the column names contain the level of the by variable?
##' Not relevant when \code{effects=\"contrast\"}.
##' @param ... Not used. For compatibility with the generic method.

## * vcov.mlmm (code)
##' @export
vcov.mlmm <- function(object, effects = "contrast", robust = object$args$robust, type.information = object$object$type.information, 
                      newdata = NULL, p = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, simplify = TRUE, ...){

    options <- LMMstar.options()
    pool.method <- options$pool.method
    adj.method <- options$adj.method

    ## ** normalize use input

    ## *** dots
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character vector")
    }
    valid.effects <- c("contrast","mean","fixed","variance","correlation","all")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }    
    if("contrast" %in% effects && length(effects)>1){
        stop("Argument \'effects\' must have length 1 when containing the element \'effects\'. \n")
    }
    if("all" %in% effects && length(effects)>1){
        stop("Argument \'effects\' must have length 1 when containing the element \'all\'. \n")
    }

    ## *** transformation
    test.sigma <- (is.null(transform.sigma) || transform.sigma == object$args$transform.sigma)
    test.k <- (is.null(transform.k) || transform.k == object$args$transform.k)
    test.rho <- (is.null(transform.rho) || transform.rho == object$args$transform.rho)
    test.notransform <- test.sigma & test.k & test.rho

    ## *** newdata
    if(!is.null(newdata)){
        message("Argument \'newdata\' is being ignored. \n")
    }

    ## ** extract
    if(all(effects=="contrast")){
        if((length(unlist(p))==0) && (robust == object$args$robust) && (type.information == object$object$type.information) && test.notransform){

            out <- object$vcov

        }else{

            e.iid <- iid.mlmm(object, effects = "contrast", p = p, robust = robust, type.information = type.information,
                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
            out <- crossprod(e.iid)
            attr(out,"original.name") <- attr(e.iid,"original.name")
            attr(out,"by") <- attr(e.iid,"by")
            attr(out,"message") <- attr(e.iid,"message")
            
        }
    }else{

        e.iid <- iid.mlmm(object, effects = effects, p = p, robust = robust, type.information = type.information,
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names, simplify = simplify)
        if(is.matrix(e.iid)){
            out <- crossprod(e.iid)
            attr(out,"original.name") <- attr(e.iid,"original.name")
            attr(out,"by") <- attr(e.iid,"by")
            attr(out,"message") <- attr(e.iid,"message")
        }else if(is.list(e.iid)){

            out <- lapply(e.iid, function(iIID){
                iOut <- crossprod(iIID)
                attr(iOut,"message") <- attr(iIID,"message")
                return(iOut)
            })
        }
        
    }

    ## ** export
    return(out)
}


##----------------------------------------------------------------------
### vcov.R ends here
