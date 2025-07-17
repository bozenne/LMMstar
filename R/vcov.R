### vcov.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:28) 
## Version: 
## Last-Updated: jul 17 2025 (10:20) 
##           By: Brice Ozenne
##     Update #: 1230
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
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
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
vcov.lmm <- function(object, effects = NULL, robust = FALSE, df = FALSE, 
                     newdata = NULL, p = NULL,
                     type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    
    ## ** normalize user imput
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){ ## hidden argument for internal use
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
        attr(vcov,"type.information") <- attr(object$vcov,"type.information")
        attr(vcov,"robust") <- attr(object$vcov,"robust")
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


## * vcov.Wald_lmm
##' @title Extract the Variance-Covariance Matrix From Wald Tests
##' @description Extract the variance-covariance matrix of the linear contrasts involved in Wald tests.
##'
##' @param object a \code{Wald_lmm} object.
##' @param effects [character vector] should the variance-covariance of the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or of the linear mixed model parameters (\code{"all"})?
##' Can also contain \code{"gradient"} to also output the gradient of the Variance-Covariance matrix.
##' @param df [logical] Should degrees-of-freedom be output?
##' @param transform.sigma,transform.k,transform.rho [character] for internal use (delta-method via estimate).
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A matrix with one column and column per parameter. \itemize{
##' \item \code{df=TRUE}: with an attribute \code{"df"} containing a numeric vector with one element per parameter. \cr
##' \item \code{effects} includes \code{"gradient"}: with an attribute \code{"gradient"} containing a 3 dimensional array with dimension the number of parameters.
##' }
##' 
##' @export
vcov.Wald_lmm <- function(object, effects = "Wald", df = FALSE,
                          transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

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
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling anova. \n")
        return(invisible(NULL))
    }
    table.param <- stats::model.tables(object, effects = "param")
    if(object$args$p.null || ("Wald" %in% effects && "gradient" %in% effects == FALSE)){
        p <- NULL
    }else{
        p <- stats::setNames(table.param$value, table.param$name)
    }

    ## *** transform
    if("Wald" %in% effects){
        if(!is.null(transform.sigma) || !is.null(transform.k) || !is.null(transform.rho)){
            txt.print <- paste0("transform.",c("sigma","k","rho"))[c(test.sigma,test.k,test.rho)]
            message("Argument(s) \'",paste(txt.print, collapse = "\', \'"),"\' are ignored when argument \'effects\' equals to \"Wald\". \n")    
        }
        transform.sigma <- object$args$transform.sigma
        transform.k <- object$args$transform.k
        transform.rho <- object$args$transform.rho

        if(!is.numeric(transform.names) && !is.logical(transform.names)){
            stop("Argument \'transform.names\' must be numeric or logical. \n")
        }
        if(length(transform.names)!=1){
            stop("Argument \'transform.names\' must have length 1. \n")
        }
        if(transform.names %in% c(0,1) == FALSE){
            stop("Argument \'transform.names\' must be TRUE/1, 0.5, or FALSE/0. \n")
        }
    }else{
        if(!is.null(transform.sigma)){ transform.sigma <- object$args$transform.sigma}
        if(!is.null(transform.k)){ transform.k <- object$args$transform.k}
        if(!is.null(transform.rho)){ transform.rho <- object$args$transform.rho}
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
    if("all" %in% effects){
        if("Wald" %in% effects){
            stop("Incorrect value for argument \'effect\': should not simultaneously contain \"all\" and \"Wald\". \n")
        }
        out <- vcov(lmm(object), effects = effects, df = df, transform.names = transform.names,
                    robust = object$args$robust, type.information = object$args$type.information, p = p,
                    transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        return(out)
        
    }else if("gradient" %in% effects){
        if(object$args$df==0){
            stop("Argument \'effects\' cannot contain \"gradient\" when no degrees-of-freedom have been stored. \n",
                 "Consider setting the argument \'df\' to TRUE when calling lmm and anova. \n")
        }            
        if(length(unique(effects))==1){
            stop("Argument \'effects\' cannot only contain \"gradient\". \n",
                 "Consider setting adding \"Wald\" or \"all\". \n")
        }
        if(object$args$type == "auto"){
            stop("Argument \'effect\' should not contain \"gradient\" w.r.t. automatically generated Wald tests. \n",
                 "(the derivative of the variance-covariance matrix is not stored) \n",
                 "Consider specifying the linear hypothesis (argument \'effect\') when calling anova. \n")
        }        
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
    out <- lapply(object$glht,function(iGlht){
        iGlht$linfct %*% iGlht$vcov %*% t(iGlht$linfct)
    })
    if(object$args$type=="user"){
        out <- out[[1]]
    }

    ## ** add gradient
    if("gradient" %in% effects){
        object.dVcov <- attr(stats::vcov(lmm(object), df = FALSE, effects = c("all","gradient"), type.information = object$args$type.information, robust = object$args$robust, p = p,
                                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE,
                                         options = options),"gradient")

        if(transform.names){
            dparam.names <- table.param$trans.name[match(dimnames(object.dVcov)[[3]],table.param$name)]
        }else{
            dparam.names <- table.param$name[match(dimnames(object.dVcov)[[3]],table.param$name)]
        }

        ## need to subset in the case dVcov is only partially stored (only first two dimensions relevant to the contrast)
        contrast <- model.tables(object, effects = "contrast", simplify = FALSE, transform.names = FALSE)[[1]][,rownames(object.dVcov),drop=FALSE]
        n.contrast <- NROW(contrast)
        name.contrast <- rownames(contrast)
        dVcov <- array(NA, dim = c(n.contrast,n.contrast,length(dparam.names)),
                       dimnames = list(name.contrast,name.contrast,dparam.names))
        for(iParam in 1:length(dparam.names)){ ## iParam <- 1
            dVcov[,,iParam] <- contrast %*% object.dVcov[,,iParam] %*% t(contrast)
        }
        
        attr(out,"gradient") <- dVcov
    }
    
    ## ** add degrees-of-freedom
    if(df){
        if(object$args$type=="user"){
            ## NOTE: to retrieve the degrees-of-freedom from the gradient (assuming identity contrast)
            ## object.vcov <- stats::vcov(lmm(object), effects = "all")
            ## 2*diag(out)^2 / sum((object.vcov %*% apply(dVcov,3,identity)) * apply(dVcov,3,identity))
            attr(out, "df") <- stats::setNames(object$univariate$df,rownames(object$univariate))
        }else{
            for(iName in names(out)){ ## iName <- names(out)[1]
                attr(out[[iName]], "df") <- stats::setNames(object$univariate[object$univariate$name==iName,"df"],rownames(object$univariate)[object$univariate$name==iName])
            }        
        }
    }

    ## ** export
    return(out)
    
}

## * vcov.rbindWald_lmm
##' @title Extract the Variance-Covariance Matrix From Combined Wald Tests
##' @description Extract the variance-covariance matrix from Wald tests applied to a linear mixed models.
##'
##' @param object a \code{rbindWald_lmm} object.
##' @param effects [character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or the value of the linear mixed model parameters (\code{"all"})?
##' Can also contain \code{"gradient"} to also output the gradient of the Variance-Covariance matrix.
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' @param transform.sigma,transform.k,transform.rho [character] for internal use (delta-method via estimate).
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' Only relevant when \code{effects="all"}.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @keywords methods
##' @return A matrix with one column and column per parameter.  \itemize{
##' \item \code{effects} includes \code{"gradient"}: with an attribute \code{"gradient"} containing a 3 dimensional array with dimension the number of parameters.
##' }
##' 
##' @export
vcov.rbindWald_lmm <- function(object, effects = "Wald", ordering = NULL,
                               transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

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

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }
    if(object$args$p.null || ("Wald" %in% effects && "gradient" %in% effects == FALSE)){
        p <- NULL
    }else{
        table.param <- stats::model.tables(object, effects = "param")
        p <- stats::setNames(table.param$value, table.param$name)
    }

    ## *** ordering
    if(!is.null(ordering)){
        ordering <- match.arg(ordering, c("contrast","model"))
        table.param.order <- object$param[order(object$param[[switch(ordering, "model" = "model", "contrast" = "name")]]),,drop=FALSE]
    }

    ## *** transform
    if("Wald" %in% effects){
        if(!is.null(transform.sigma) || !is.null(transform.k) || !is.null(transform.rho)){
            txt.print <- paste0("transform.",c("sigma","k","rho"))[c(test.sigma,test.k,test.rho)]
            message("Argument(s) \'",paste(txt.print, collapse = "\', \'"),"\' are ignored when argument \'effects\' equals to \"Wald\". \n")    
        }
        transform.sigma <- object$args$transform.sigma
        transform.k <- object$args$transform.k
        transform.rho <- object$args$transform.rho

        if(!is.numeric(transform.names) && !is.logical(transform.names)){
            stop("Argument \'transform.names\' must be numeric or logical. \n")
        }
        if(length(transform.names)!=1){
            stop("Argument \'transform.names\' must have length 1. \n")
        }
        if(transform.names %in% c(0,1) == FALSE){
            stop("Argument \'transform.names\' must be TRUE/1, 0.5, or FALSE/0. \n")
        }
    }else{
        if(is.null(transform.sigma)){ transform.sigma <- object$args$transform.sigma}
        if(is.null(transform.k)){ transform.k <- object$args$transform.k}
        if(is.null(transform.rho)){ transform.rho <- object$args$transform.rho}
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

    if("all" %in% effects && "Wald" %in% effects){
        stop("Incorrect value for argument \'effect\': should not simultaneously contain \"all\" and \"Wald\". \n")
    }

    if("gradient" %in% effects){
        if(length(unique(effects))==1){
            stop("Argument \'effects\' cannot only contain \"gradient\". \n",
                 "Consider setting adding \"Wald\" or \"all\". \n")
        }
        if(object$args$type == "auto"){
            stop("Argument \'effect\' should not contain \"gradient\" w.r.t. automatically generated Wald tests. \n",
                 "(the derivative of the variance-covariance matrix is not stored) \n",
                 "Consider specifying the linear hypothesis (argument \'effect\') when calling anova. \n")
        }
    }

    ## ** extract full variance-covariance matrix
    if("all" %in% effects || "gradient" %in% effects){
        ls.model <- lmm(object)
        seq.cluster <- .rbind.cluster(ls.model)
        all.vcov <- .rbind.vcov(ls.model, robust = object$args$robust, type.information = object$args$type.information, keep.grad = "gradient" %in% effects, p = p,
                                transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                                seq.cluster = seq.cluster, n.cluster = length(seq.cluster), independence = object$args$independence,
                                all.table.param = object$param, all.coefUnames = object$param$trans.Uname, all.coefUnamesO = object$param$Uname, options = options)
        
    }

    ## ** process variance-covariance matrix
    if("all" %in% effects){
        if(transform.names){
            out <- all.vcov[object$param$Uname,object$param$Uname,drop=FALSE]
            dimnames(out) <- list(object$param$trans.Uname,object$param$trans.Uname)
            if("gradient" %in% effects){
                attr(out,"gradient") <- attr(all.vcov,"gradient")[object$param$Uname,object$param$Uname,object$param$Uname,drop=FALSE]
                dimnames(attr(out,"gradient")) <- list(object$param$trans.Uname,object$param$trans.Uname,object$param$trans.Uname)
            }
        }else{
            out <- all.vcov
            attr(out,"iid") <- NULL
        }
        if(!is.null(ordering)){
            if(transform.names){
                out <- out[table.param.order$trans.Uname,table.param.order$trans.Uname,drop=FALSE]
                if("gradient" %in% effects){
                    attr(out,"gradient") <- attr(out,"gradient")[table.param.order$trans.Uname,table.param.order$trans.Uname,table.param.order$trans.Uname,drop=FALSE]
                }
            }else{
                out <- out[table.param.order$Uname,table.param.order$Uname,drop=FALSE]
                if("gradient" %in% effects){
                    attr(out,"gradient") <- attr(out,"gradient")[table.param.order$Uname,table.param.order$Uname,table.param.order$Uname,drop=FALSE]
                }
            }
        }

    }else{

        object$param$name <- object$param$Uname
        object$param$trans.name <- object$param$Uname.trans
        if("gradient" %in% effects){
            object$glht[[1]]$dVcov <- attr(all.vcov, "gradient")[colnames(object$glht[[1]]$vcov),colnames(object$glht[[1]]$vcov),,drop=FALSE]
        }
        out <- vcov.Wald_lmm(object, effects = setdiff(effects,"gradient"), df = FALSE, transform.names = FALSE)
        if("gradient" %in% effects){

            if(transform.names){
                table.param <- stats::model.tables(object, effects = "param")
                dparam.names <- table.param$trans.Uname[match(dimnames(attr(all.vcov,"gradient"))[[3]],table.param$Uname)]
            }else{
                dparam.names <- dimnames(attr(all.vcov,"gradient"))[[3]]
            }

            ## need to subset in the case dVcov is only partially stored (only first two dimensions relevant to the contrast)
            contrast <- model.tables(object, effects = "contrast", simplify = FALSE, transform.names = FALSE)[[1]][,rownames(attr(all.vcov,"gradient")),drop=FALSE]
            n.contrast <- NROW(contrast)
            name.contrast <- rownames(contrast)
            attr(out,"gradient") <- array(NA, dim = c(n.contrast,n.contrast,length(dparam.names)),
                                          dimnames = list(name.contrast,name.contrast,dparam.names))
            for(iParam in 1:length(dparam.names)){ ## iParam <- 1
                attr(out,"gradient")[,,iParam] <- contrast %*% attr(all.vcov,"gradient")[,,iParam] %*% t(contrast)
            }
            
        }

        if(!is.null(ordering)){
            
            if(ordering=="model"){
                reorder <- order(object$univariate$model)
            }else if(is.list(object$univariate$term)){
                reorder <- order(object$univariate$type,sapply(object$univariate$term, paste, collapse = ";"))
            }else{
                reorder <- order(object$univariate$type,object$univariate$term)
            }
            if("gradient" %in% effects){
                out.dVcov <- attr(out,"gradient")[reorder,reorder,table.param.order$Uname,drop=FALSE]
            }
            out <- out[reorder,reorder,drop=FALSE]
            if("gradient" %in% effects){
                attr(out,"gradient") <- out.dVcov
            }
            
        }

    }

    ## ** export
    return(out)
}

## * vcov.mlmm (documentation)
##' @title Extract The Variance-Covariance Matrix From Multiple Linear Mixed Models
##' @description Extract the variance-covariance matrix of the model coefficients from multiple linear mixed models.
##'
##' @param object a \code{mlmm} object.
##' @inheritParams vcov.rbindWald_lmm
##' 
##' @keywords methods
##' @return A matrix with one column and column per parameter.  \itemize{
##' \item \code{effects} includes \code{"gradient"}: with an attribute \code{"gradient"} containing a 3 dimensional array with dimension the number of parameters.
##' }
##' 
## * vcov.mlmm (code)
##' @export
vcov.mlmm <- vcov.rbindWald_lmm


##----------------------------------------------------------------------
### vcov.R ends here
