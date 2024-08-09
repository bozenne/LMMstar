### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: aug  9 2024 (17:39) 
##           By: Brice Ozenne
##     Update #: 1037
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
##' Can also contain \code{"gradient"} to also output the gradient of the Variance-Covariance matrix .
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors.
##' Can also be \code{2} compute the degrees-of-freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
##' @param df [logical] Should degrees-of-freedom, computed using Satterthwaite approximation, for the model parameters be output.
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
##' @return A matrix with one column and column per parameter. \itemize{
##' \item \code{df=TRUE}: with an attribute \code{"df"} containing a numeric vector with one element per parameter. \cr
##' \item \code{effects} includes \code{"gradient"}: with an attribute \code{"gradient"} containing a 3 dimensional array with dimension the number of parameters.
##' }
##' 
##' @keywords methods 

## * vcov.lmm (code)
##' @export
vcov.lmm <- function(object, effects = NULL, robust = FALSE, df = FALSE, strata = NULL,
                     newdata = NULL, p = NULL,
                     type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    
    ## ** normalize user imput
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
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
                stop("Argument \'effects\' cannot contain \"gradient\" when no degrees-of-freedom have been stored. \n",
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
    table.param <- stats::model.tables(object, effects = c("param",effects),
                                        transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)

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

    ## *** robust
    if(length(df)>1){
        stop("Argument \'robust\' must have length 1. \n")
    }
    if(((is.numeric(robust) && df %in% 0:2 == FALSE) && !is.logical(df)) || !is.vector(df)){
        stop("Argument \'df\' must be a logical value or an integer 0, 1, or 2. \n")
    }

    ## *** df
    if(length(df)>1){
        stop("Argument \'df\' must have length 1. \n")
    }
    if(((is.numeric(df) && df %in% 0:1 == FALSE) && !is.logical(df)) || !is.vector(df)){
        stop("Argument \'df\' must be a logical value. \n")
    }
    if(df && object$args$df==0){
        stop("Argument \'df\' cannot be TRUE when no degrees-of-freedom have been stored. \n",
             "Consider setting the argument \'df\' to TRUE when calling lmm. \n")
    }

    ## ** extract or recompute variance covariance matrix
    test.new <- !is.null(newdata) || !is.null(p) || !test.notransform || object$args$type.information!=type.information || (df == TRUE && !object$args$df) || (robust == 2)
    ## need for new computation if: new dataset, param, transform, information, or df not computed before, or robust df instead of model-based df 
    if((test.new == FALSE) && (robust == FALSE)){
        keep.name <- stats::setNames(table.param$name, table.param$trans.name)

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
            attr(vcov,"gradient") <- object$dVcov[keep.name,keep.name,,drop=FALSE]
            if(transform.names){
                newname.all <- dimnames(object$dVcov)[[3]]
                newname.all[match(names(object$reparametrize$p),newname.all)] <- object$reparametrize$newname
                dimnames(attr(vcov,"gradient")) <- list(names(keep.name),names(keep.name),newname.all)
            }
        }

    }else{

        ## *** re-compute variance (and possibly df and gradient)
        if(!is.null(newdata)){
            design <- stats::model.matrix(object, newdata = newdata, effects = "all", simplify = FALSE, options = options)
        }else{
            design <- object$design
        }
        outMoments <- .moments.lmm(value = theta, design = design, time = object$time, method.fit = object$args$method.fit, type.information = type.information,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   logLik = FALSE, score = FALSE, information = FALSE,  vcov = TRUE, df = (df || keep.grad) && test.new, indiv = FALSE, effects = effects, robust = robust,
                                   trace = FALSE, precompute.moments = !is.null(object$design$precompute.XX), method.numDeriv = options$method.numDeriv, transform.names = transform.names)

        trans2original.names <- colnames(outMoments$vcov) ## output names without transformation
        if(any(trans2original.names %in% outMoments$reparametrize$newname)){
            index.trans2original <- match(outMoments$reparametrize$newname,trans2original.names)
            trans2original.names[stats::na.omit(index.trans2original)] <- names(outMoments$reparametrize$p)[!is.na(index.trans2original)]
        }

        if(test.new == FALSE && df){
            outMoments$df <- stats::setNames(object$df[trans2original.names], colnames(outMoments$vcov))
        }
        if(test.new == FALSE && keep.grad){
            outMoments$dVcov <- object$dVcov[trans2original.names,trans2original.names,,drop=FALSE]
            dimnames(outMoments$dVcov) <- list(colnames(outMoments$vcov),colnames(outMoments$vcov),dimnames(outMoments$dVcov)[[3]])
        }

        ## *** re-order and prepare export
        if("variance" %in% effects && transform.k %in% c("sd","var","logsd","logvar") && object$strata$n>1 && transform.names){
            ## re-order values when converting to sd with strata (avoid sd0:0 sd0:1 sd1:0 sd1:1 sd2:0 sd2:1 ...)a
            vcov <- outMoments$vcov[table.param$trans.name,table.param$trans.name,drop=FALSE]
            if(df){
                attr(vcov,"df") <- outMoments$df[table.param$trans.name]
            }
            if(keep.grad){
                attr(vcov,"gradient") <- outMoments$dVcov[table.param$trans.name,table.param$trans.name,table.param$trans.name,drop=FALSE]
            }
        }else{
            vcov <- outMoments$vcov
            if(df){
                attr(vcov,"df") <- outMoments$df
            }
            if(keep.grad){
                attr(vcov,"gradient") <- outMoments$dVcov
            }
        }

        ## *** retrieve model-based vcov for gradient when robust is 1
        if(robust == 1 && keep.grad){
            if(test.new){
                iVcov <-  .moments.lmm(value = theta, design = design, time = object$time, method.fit = object$args$method.fit, type.information = type.information,
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                       logLik = FALSE, score = FALSE, information = FALSE,  vcov = TRUE, df = FALSE, indiv = FALSE, effects = effects, robust = FALSE,
                                       trace = FALSE, precompute.moments = !is.null(object$design$precompute.XX), method.numDeriv = options$method.numDeriv, transform.names = transform.names)$vcov
                attr(attr(vcov,"gradient"),"vcov") <- iVcov[rownames(vcov),colnames(vcov)] ## possible re-ordering to account for transformation + drop attributes by subsetting
            }else{
                attr(attr(vcov,"gradient"),"vcov") <- object$vcov[trans2original.names,trans2original.names] ## drop attributes by subsetting
                dimnames(attr(attr(vcov,"gradient"),"vcov")) <- list(rownames(vcov),colnames(vcov))
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
##' Can also be \code{2} compute the degrees-of-freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
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


    ## ** normalize use input

    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    pool.method <- options$pool.method
    adj.method <- options$adj.method
    
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
##' @title Extract the Variance-Covariance From Wald Tests For Linear Mixed Models
##' @description Extract the variance-covariance matrix from Wald tests applied to a linear mixed models.
##'
##' @param object a \code{rbindWald_lmm} object.
##' @param effects [character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or the value of the linear mixed model parameters (\code{"all"})?
##' Can also contain \code{"gradient"} to also output the gradient of the Variance-Covariance matrix.
##' @param method [character vector] type of adjustment for multiple comparisons across the linear contrasts (one of \code{"none"}, \code{"bonferroni"}, ..., \code{"single-step2"}).
##' Only relevant when \code{effects = "Wald"}.
##' @param df [logical] Should degrees-of-freedom, computed using Satterthwaite approximation, for the model parameters be output.
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' Only relevant when \code{effects="all"}.
##' @param simplify [logical] should the output be a vector or a list with one element specific to each possible ordering (i.e. contrast or model).
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A matrix with one column and column per parameter. 
##' 
##' @export
vcov.rbindWald_lmm <- function(object, effects = "Wald", method = "none", df = FALSE, ordering = NULL, transform.names = TRUE, simplify = TRUE, ...){

    table.param <- stats::model.tables(object, effects = "param")
    user.call <- match.call()

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    adj.method <- options$adj.method

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** effects
    ## initialized in vcov.Wald_lmm

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

    ## *** transform.names
    if(!is.numeric(transform.names) && !is.logical(transform.names)){
        stop("Argument \'transform.names\' must be numeric or logical. \n")
    }
    if(length(transform.names)!=1){
        stop("Argument \'transform.names\' must have length 1. \n")
    }
    if(transform.names %in% c(0,1) == FALSE){
        stop("Argument \'transform.names\' must be TRUE/1, 0.5, or FALSE/0. \n")
    }

    ## ** extract

    ## *** select value 
    value.out <- vcov.Wald_lmm(object, effects = effects, df = df, transform.names = FALSE)
    
    ## *** rename
    if(transform.names == FALSE){
        if(simplify==FALSE && !is.null(ordering)){
            newname <- table.param[[switch(ordering, "model" = "name", "contrast" = "model")]]                
        }else{
            newname <- table.param$Uname
        }
    }else{
        if(simplify==FALSE && !is.null(ordering)){
            newname <- table.param[[switch(ordering, "model" = "trans.name", "contrast" = "model")]]
        }else{
            newname <- table.param$trans.Uname
        }
    }
       
    ## vcov
    if("all" %in% effects){
        dimnames(value.out) <- list(newname,newname)
    }

    ## gradient
    if("gradient" %in% effects){
        index.name3 <- match(dimnames(attr(iOut,"gradient"))[[3]], table.param$Uname)

        if("all" %in% effects){
            index.name12 <- match(dimnames(attr(iOut,"gradient"))[[1]], table.param$Uname)
            dimnames(attr(value.out,"gradient")) <- list(newname[[index.name12]],newname[[index.name12]],newname[[index.name3]])
        }else if("Wald" %in% effects){
            dimnames(attr(value.out,"gradient")) <- c(dimnames(attr(value.out,"gradient"))[1:2],list(newname[[index.name3]]))
        }       
    }
            
    ## *** ordering
    if(!is.null(ordering)){        
        if("Wald" %in% effects){ ## effects can have an extra element "gradient", i.e. be a vector
            table.univariate <- object$univariate
            ordering.var <- switch(ordering, "model" = "model", "contrast" = "term")
            ordering.out <- factor(table.univariate[[ordering.var]], unique(table.univariate[[ordering.var]])) 
        }else if("all" %in% effects){
            if(transform.names && simplify){
                ordering.var <- switch(ordering, "model" = "model", "contrast" = "trans.name")
            }else{
                ordering.var <- switch(ordering, "model" = "model", "contrast" = "name")
            }
            ordering.out <- factor(table.param[[ordering.var]], unique(table.param[[ordering.var]]))
        }
    }else if(simplify == FALSE){
        if("Wald" %in% effects){
            attr(value.out,"parameter") <- object$univariate$term
            attr(value.out,"model") <- object$univariate$model
        }else if("all" %in% effects){
            attr(value.out,"parameter") <- table.param$name
            attr(value.out,"model") <- table.param$model
        }
    }

    ## ** split
    if(simplify == FALSE && !is.null(ordering)){
        out <- tapply(1:length(ordering.out), ordering.out, FUN = function(iIndex){ ## iIndex <- 1:2
            iOut.grad <- attr(value.out,"gradient")
            iOut <- value.out[iIndex,iIndex,drop=FALSE]
            if("gradient" %in% effects){
                attr(iOut,"gradient") <- iOut.grad[,,iIndex,drop=FALSE]
            }
            return(iOut)
        }, simplify = FALSE)
    }else{
        if(!is.null(ordering)){
            neworder <- order(ordering.out)
            out <- value.out[neworder,neworder,drop=FALSE]
            if("gradient" %in% effects){
                attr(out,"gradient") <- attr(value.out,"gradient")[,,neworder,drop=FALSE]
            }
        }else{
            out <- value.out

        }
    }

    ## ** export
    return(out)
}

## * vcov.Wald_lmm
##' @title Extract the Variance-Covariance Matrix From Wald Tests For Linear Mixed Model
##' @description Extract the variance-covariance matrix of the linear contrasts involved in the Wald test.
##'
##' @param object a \code{Wald_lmm} object.
##' @param effects [character vector] should the variance-covariance of the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or of the linear mixed model parameters (\code{"all"})?
##' Can also contain \code{"gradient"} to also output the gradient of the Variance-Covariance matrix.
##' @param df [logical] Should degrees-of-freedom, computed using Satterthwaite approximation, for the model parameters be output.
##' Also output the first derivative of the variance-covariance matrix whenever the argument is stricly greater than 1.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A matrix with one column and column per parameter. \itemize{
##' \item \code{df=TRUE}: with an attribute \code{"df"} containing a numeric vector with one element per parameter. \cr
##' \item \code{effects} includes \code{"gradient"}: with an attribute \code{"gradient"} containing a 3 dimensional array with dimension the number of parameters.
##' }
##' 
##' @export
vcov.Wald_lmm <- function(object, effects = "Wald", df = FALSE, transform.names = TRUE, ...){

    table.param <- stats::model.tables(object, effects = "param")

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- NULL
    }
    dots$options <- NULL
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
            stop("Argument \'effects\' cannot contain \"gradient\" when no degrees-of-freedom have been stored. \n",
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
        stop("Argument \'df\' cannot be TRUE when no degrees-of-freedom have been stored. \n",
             "Consider setting the argument \'df\' to TRUE when calling lmm and anova. \n")
    }

    ## ** extract
    if(effects=="Wald"){
        out <- lapply(object$glht,function(iGlht){
            iGlht$linfct %*% iGlht$vcov %*% t(iGlht$linfct)
        })
        if(object$args$type=="user"){
            out <- out[[1]]
        }
    }else if(effects == "all"){
        out <- object$glht[[1]]$vcov
        if(transform.names){
            trans.names <- table.param$trans.name[match(rownames(out),table.param$name)]
            dimnames(out) <- list(trans.names,trans.names)
        }
    }

    ## ** add keep.grad
    if(keep.grad || (df && effects == "all")){
        
        if(transform.names){
            dparam.names <- table.param$trans.name[match(dimnames(object$glht[[1]]$dVcov)[[3]],table.param$name)]
        }else if(effects == "Wald"){
            dparam.names <- table.param$name[match(dimnames(object$glht[[1]]$dVcov)[[3]],table.param$name)]
        }

        if(effects == "Wald"){
            n.contrast <- NROW(ls.contrast[[1]])
            name.contrast <- rownames(ls.contrast[[1]])
            dVcov <- array(NA, dim = c(n.contrast,n.contrast,length(dparam.names)),
                           dimnames = list(name.contrast,name.contrast,dparam.names))
            contrast <- ls.contrast[[1]][,rownames(object$glht[[1]]$dVcov),drop=FALSE] ## need to subset in the case dVcov is only partially stored (only first two dimension relevant to the contrast)
            for(iParam in 1:length(dparam.names)){ ## iParam <- 1
                dVcov[,,iParam] <- contrast %*% object$glht[[1]]$dVcov[,,iParam] %*% t(contrast)
            }
        }else if(effects == "all"){
            dVcov <- object$glht[[1]]$dVcov
            if(transform.names){
                ## trans.names[match(dimnames(dVcov)[[1]], table.param$name)] may differ from trans.names when only a subset of dVcov is exported (simplify = 0.5) 
                dimnames(dVcov) <- list(table.param$trans.name[match(dimnames(dVcov)[[1]], table.param$name)],
                                        table.param$trans.name[match(dimnames(dVcov)[[2]], table.param$name)],
                                        dparam.names)
            }
        }

        if(keep.grad){
            attr(out,"gradient") <- dVcov
        }
    }
    
    ## ** add degrees-of-freedom
    if(df){
        if(effects == "Wald"){
            if(object$args$type=="user"){
                attr(out, "df") <- setNames(object$univariate$df,rownames(object$univariate))
            }else{
                for(iName in names(out)){ ## iName <- names(out)[1]
                    attr(out[[iName]], "df") <- setNames(object$univariate[object$univariate$name==iName,"df"],rownames(object$univariate)[object$univariate$name==iName])
                }        
            }
        }else if(effects == "all"){
            iDF <- .df_contrast(contrast = NULL, vcov = out, dVcov.param = dVcov)
            attr(out,"df") <- stats::setNames(rep(NA, NCOL(out)), colnames(out)) ## handle the case where simplify is used, i.e., not all dVcov is storeda
            attr(out,"df")[names(iDF)] <- iDF
        }
        
    }

    ## ** export
    return(out)
    
}

##----------------------------------------------------------------------
### vcov.R ends here
