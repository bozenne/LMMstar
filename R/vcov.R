### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: jul 31 2024 (15:37) 
##           By: Brice Ozenne
##     Update #: 842
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
##' @param effects [character vector] Should the variance-covariance matrix for all coefficients be output (\code{"all"}),
##' or only for coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only for coefficients relative to the variance structure (\code{"variance"}),
##' or only for coefficients relative to the correlation structure (\code{"correlation"}).
##' Can also contain \code{"gradient"} to also output the gradient of the Variance Covariance Matrix .
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
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef.lmm} function. \cr
##' 
##' @return A matrix with one column and column per parameter.
##' An attribute \code{"df"} is added when argument df is set to \code{TRUE}, containing a numeric vector with one element per parameter.
##' An attribute \code{"dVcov"} is added when argument df is greater than 1, containing a 3 dimensional array with with dimension being the number of parameters.
##'
##' @keywords methods 

## * vcov.lmm (code)
##' @export
vcov.lmm <- function(object, effects = NULL, robust = FALSE, df = FALSE, strata = NULL,
                     newdata = NULL, p = NULL,
                     type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    
    options <- LMMstar.options()

    ## ** normalize user imput
    ## *** dots
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(is.null(effects)){
        keep.grad <- FALSE
        if((is.null(transform.sigma) || identical(transform.sigma,"none")) && (is.null(transform.k) || identical(transform.k,"none")) && (is.null(transform.rho) || identical(transform.rho,"none"))){
            effects <- options$effects
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else{
        if(!is.character(effects) || !is.vector(effects)){
            stop("Argument \'effects\' must be a character vector. \n")
        }
        valid.effects <- c("mean","fixed","variance","correlation","all","gradient")
        if(any(effects %in% valid.effects == FALSE)){
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
        }
        if("gradient" %in% effects){
            keep.grad <- TRUE
            effects <- setdiff(effects, "gradient")
            if(object$args$df==0){
                stop("Argument \'effects\' cannot contain \"gradient\" when no degrees of freedom have been stored. \n",
                     "Consider setting the argument \'df\' to TRUE when calling lmm. \n")
            }
            if(length(effects)==0){
                stop("Argument \'effects\' cannot only contain \"gradient\". \n",
                     "Consider setting adding \"mean\" or \"all\". \n")
            }
        }else{
            keep.grad <- FALSE
        }
        if(all("all" %in% effects)){
            if(length(unique(effects))>1){
                stop("When containing the element \"all\" the argument \'effects\' should not contain \"mean\", \"fixed\", \"variance\", or \"correlation\" . \n")
            }else{
                effects <- c("mean","variance","correlation")
            }
        }else{
            effects[effects == "fixed"] <- "mean"
        }
    }

    ## *** strata
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }

    ## *** type information
    if(is.null(type.information)){
        type.information <- object$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## *** transformation & p
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

    ## *** df
    if(length(df)>1){
        stop("Argument \'df\' must have length 1. \n")
    }
    if(((is.numeric(df) && df %in% 0:1 == FALSE) && !is.logical(df)) || !is.vector(df)){
        stop("Argument \'df\' must be a logical value. \n")
    }
    if(df && object$args$df==0){
        stop("Argument \'df\' cannot be TRUE when no degrees of freedom have been stored. \n",
             "Consider setting the argument \'df\' to TRUE when calling lmm. \n")
    }
    
    ## ** extract or recompute variance covariance matrix
    if(is.null(newdata) && is.null(p) && test.notransform && (df == FALSE || object$args$df) && (robust == FALSE) && object$args$type.information==type.information){
        keep.name <- stats::setNames(names(coef(object, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                                     names(coef(object, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))    

        vcov <- object$vcov[keep.name,keep.name,drop=FALSE]
        if(transform.names){
            dimnames(vcov) <- list(names(keep.name),names(keep.name))
        }
        if(df){
            attr(vcov,"df") <- object$df[keep.name]
            if(transform.names){
                names(attr(vcov,"df")) <- names(keep.name)
            }
        }
        if(keep.grad){
            attr(vcov,"dVcov") <- object$dVcov[keep.name,keep.name,keep.name,drop=FALSE]
            if(transform.names){
                dimnames(attr(vcov,"dVcov")) <- list(names(keep.name),names(keep.name),names(keep.name))
            }
        }

    }else{
         
        if(!is.null(newdata)){
            design <- stats::model.matrix(object, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- object$design
        }

        outMoments <- .moments.lmm(value = theta, design = design, time = object$time, method.fit = object$args$method.fit, type.information = type.information,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   logLik = FALSE, score = FALSE, information = FALSE, vcov = TRUE, df = df, indiv = FALSE, effects = effects, robust = robust,
                                   trace = FALSE, precompute.moments = !is.null(object$design$precompute.XX), method.numDeriv = options$method.numDeriv, transform.names = transform.names)

        if("variance" %in% effects && transform.k %in% c("sd","var","logsd","logvar") && object$strata$n>1 && transform.names){
            ## re-order values when converting to sd with strata (avoid sd0:0 sd0:1 sd1:0 sd1:1 sd2:0 sd2:1 ...)
            out.name <- names(stats::coef(object, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = TRUE))
            vcov <- outMoments$vcov[out.name,out.name,drop=FALSE]
            if(df){
                attr(vcov,"df") <- outMoments$df[out.name]
            }
            if(keep.grad){
                attr(vcov,"dVcov") <- outMoments$dVcov[out.name,out.name,out.name,drop=FALSE]
            }
        }else{
            vcov <- outMoments$vcov
            if(df){
                attr(vcov,"df") <- outMoments$df
            }
            if(keep.grad){
                attr(vcov,"dVcov") <- outMoments$dVcov
            }
        }

    }

    ## ** export
    return(vcov)    
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
##' @param ordering [character] should the output be ordered by type of parameter (\code{"parameter"}) or by model (\code{"by"}).
##' Not relevant when \code{effects="contrast"}.
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
vcov.mlmm <- function(object, effects = "contrast", p = NULL, newdata = NULL, ordering = "by",
                      robust = object$args$robust, type.information = object$object$type.information, 
                      transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, simplify = TRUE, ...){

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

    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1 or FALSE/0. \n")
    }

    ## ** extract
    if(all(effects=="contrast")){
        if((length(unlist(p))==0) && (robust == object$args$robust) && (type.information == object$object$type.information) && test.notransform){

            out <- object$vcov

        }else{

            e.iid <- iid.mlmm(object, effects = "contrast", p = p, REML2ML = REML2ML, robust = robust, type.information = type.information,
                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
            out <- crossprod(e.iid)
            attr(out,"original.name") <- attr(e.iid,"original.name")
            attr(out,"by") <- attr(e.iid,"by")
            attr(out,"message") <- attr(e.iid,"message")
            
        }
    }else if(ordering == "by" && !simplify){
        out <- lapply(object$model, FUN = vcov, effects = effects, p = p, robust = robust, type.information = type.information,
                      transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    }else{
        e.iid <- iid.mlmm(object, effects = effects, p = p, ordering = ordering, robust = robust, type.information = type.information,
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


## * vcov.rbindWald_lmm
##' @title Extract  the Variance-Covariance From Wald Tests for Linear Mixed Models
##' @description Extract  the variance-covariance matrix from Wald tests applied to a linear mixed models.
##'
##' @param object a \code{rbindWald_lmm} object.
##' @param effects [character] should the variance-covariance matrix involved in the Wald tests be output (\code{"contrast"}),
##' or the variance-covariance matrix of the linear mixed model parameters (\code{"all"})?
##' @param method [character vector] type of adjustment for multiple comparisons across the linear contrasts (one of \code{"none"}, \code{"bonferroni"}, ..., \code{"single-step2"}).
##' Only relevant when \code{effects = "contrast"}.
##' @param df [logical] Should degree of freedom, computed using Satterthwaite approximation, for the model parameters be output.
##' Also output the first derivative of the variance-covariance matrix whenever the argument is stricly greater than 1.
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' Only relevant when \code{effects="contrast"}, \code{effects="all"}, or \code{effects="all.original"}.
##' @param simplify [logical] should the output be a vector or a list with one element specific to each possible ordering (i.e. contrast or model).
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A matrix with one column and column per parameter.
##' An attribute \code{"df"} is added when argument df is set to \code{TRUE}, containing a numeric vector with one element per parameter.
##' An attribute \code{"dVcov"} is added when argument df is greater than 1, containing a 3 dimensional array with with dimension being the number of parameters.
##'
##' @export
vcov.rbindWald_lmm <- function(object, effects = "Wald", method = "none", df = FALSE, ordering = NULL, simplify = TRUE, ...){

    options <- LMMstar.options()
    adj.method <- options$adj.method

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character.")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1.")
    }
    valid.effects <- c("Wald","all")
    if(effects %in% valid.effects == FALSE){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }

    ## *** method
    if(!is.character(method) || !is.vector(method)){
        stop("Argument \'method\' must be a character.")
    }
    if(length(method)!=1){
        stop("Argument \'method\' must have length 1.")
    }
    valid.method <- c("none",adj.method)
    if(any(method %in% valid.method == FALSE)){
        stop("Incorrect value for argument \'method\': \"",paste(setdiff(method,valid.method), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.method, collapse ="\", \""),"\". \n")
    }

    ## *** ordering
    if(!is.null(ordering)){
        ordering <- match.arg(ordering, c("contrast","model"))
    }

    ## ** extract
    out <- vcov.Wald_lmm(object, effects = effects, df = df)

    if(!is.null(ordering)){
        if(effects=="Wald"){
            ordering.out <- object$univariate[[c(model = "model", contrast = "term")[ordering]]]
        }else if(effects=="all"){
            ordering.out <- object$glht[[1]][[c(model = "model", contrast = "term")[ordering]]]
        }
    }else if(simplify == FALSE){
        if(effects=="Wald"){
            attr(out,"parameter") <- object$univariate$term
            attr(out,"model") <- object$univariate$model
        }else if(effects=="all"){
            attr(out,"parameter") <- object$glht[[1]]$term
            attr(out,"model") <- object$glht[[1]]$model
        }
    }

    ## ** export
    if(simplify == FALSE && !is.null(ordering)){
        out <- tapply(1:length(ordering.out), factor(ordering.out, levels = unique(ordering.out)), FUN = function(iIndex){
            return(out[iIndex,iIndex,drop=FALSE])
        }, simplify = FALSE)
    }else if(!is.null(ordering)){
        neworder <- unlist(tapply(1:length(ordering.out), factor(ordering.out, levels = unique(ordering.out)), base::identity, simplify = FALSE))
        out <- out[neworder,neworder,drop=FALSE]
    }
    return(out)
}

## * vcov.Wald_lmm
##' @title Extract the Variance-Covariance Matrix from Wald Tests
##' @description Extract the variance-covariance matrix of the linear contrasts involved in the Wald test.
##'
##' @param object a \code{Wald_lmm} object.
##' @param effects [character vector] should the variance-covariance of the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or of the linear mixed model parameters (\code{"all"})?
##' Can also contain \code{"gradient"} to also output the gradient of the Variance Covariance Matrix.
##' @param df [logical] Should degree of freedom, computed using Satterthwaite approximation, for the model parameters be output.
##' Also output the first derivative of the variance-covariance matrix whenever the argument is stricly greater than 1.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A matrix with one column and column per parameter.
##' An attribute \code{"df"} is added when argument \code{df} is set to \code{TRUE}, containing a numeric vector with one element per parameter.
##' An attribute \code{"dVcov"} is added when argument \code{effects} contain \code{"gradient"}, containing a 3 dimensional array with with dimension being the number of parameters.
##' 
##' @export
vcov.Wald_lmm <- function(object, effects = "Wald", df = FALSE, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character. \n")
    }
    valid.effects <- c("Wald","all","gradient")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }
    if("gradient" %in% effects){
        keep.grad <- TRUE
        effects <- setdiff(effects, "gradient")
        if(object$args$df==0){
            stop("Argument \'effects\' cannot contain \"gradient\" when no degrees of freedom have been stored. \n",
                 "Consider setting the argument \'df\' to TRUE when calling lmm and anova. \n")
        }            
        if(length(effects)==0){
            stop("Argument \'effects\' cannot only contain \"gradient\". \n",
                 "Consider setting adding \"mean\" or \"all\". \n")
        }
    }else{
        keep.grad <- FALSE
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' should not contain \"Wald\" and \"all\". \n")
    }

    ## *** df
    if(length(df)>1){
        stop("Argument \'df\' must have length 1. \n")
    }
    if(((is.numeric(df) && df %in% 0:1 == FALSE) && !is.logical(df)) || !is.vector(df)){
        stop("Argument \'df\' must be a logical value. \n")
    }
    if(df && object$args$df==0){
        stop("Argument \'df\' cannot be TRUE when no degrees of freedom have been stored. \n",
             "Consider setting the argument \'df\' to TRUE when calling lmm and anova. \n")
    }
    
    ## ** extract
    out <- lapply(object$glht,"[[","vcov")
    if(effects=="Wald"){
        ls.contrast <- model.tables(object, effects = "contrast", simplify = FALSE)
        out <- mapply(iVcov = out, iC = ls.contrast, FUN = function(iVcov,iC){iC %*% iVcov %*% t(iC)},
                      SIMPLIFY = FALSE)
    }else if(effects == "all"){
        trans.names <- names(coef(object, effects = "all"))
        dimnames(out[[1]]) <- list(trans.names,trans.names)
    }

    if(object$args$type=="user"){
        out <- out[[1]]
    }

    ## ** add keep.grad
    if(keep.grad || (df && effects == "all")){
        
        if(effects == "all"){
            original.names <- names(coef(object, effects = "all", backtransform = TRUE))
            dVcov <- object$glht[[1]]$dVcov
            dimnames(dVcov) <- list(trans.names[match(dimnames(dVcov)[[1]], original.names)], trans.names[match(dimnames(dVcov)[[1]], original.names)], trans.names)
        }else if(effects == "Wald"){
            trans.names <- names(coef(object, effects = "all"))
            n.contrast <- NROW(ls.contrast[[1]])
            name.contrast <- rownames(ls.contrast[[1]])
            
            dVcov <- array(NA, dim = c(n.contrast,n.contrast,length(trans.names)), dimnames = list(name.contrast,name.contrast,trans.names))
            for(iParam in 1:length(trans.names)){
                dVcov[,,iParam] <- ls.contrast[[1]] %*% object$glht[[1]]$dVcov[,,iParam] %*% t(contrast[[1]])
            }
        }

        if(keep.grad){
            attr(out,"dVcov") <- dVcov
        }
    }
    
    ## ** add degrees of freedom
    if(df){
        if(object$args$type=="user"){
            if(effects == "all"){
                iDF <- .dfX(X.beta = NULL, vcov = out, dVcov.param = dVcov)
                attr(out,"df") <- stats::setNames(rep(NA, length(trans.names)), trans.names) ## handle the case where simplify is used, i.e., not all dVcov is stored
                attr(out,"df")[names(iDF)] <- iDF
            }else if(effects == "Wald"){
                attr(out, "df") <- setNames(object$univariate$df,rownames(object$univariate))
            }
        }else{
            for(iName in names(out)){ ## iName <- names(out)[1]
                attr(out[[iName]], "df") <- setNames(object$univariate[object$univariate$name==iName,"df"],rownames(object$univariate)[object$univariate$name==iName])
            }        
        }
    }

    ## ** export
    return(out)
    
}

##----------------------------------------------------------------------
### vcov.R ends here
