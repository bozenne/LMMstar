### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  4 2021 (10:04) 
## Version: 
## Last-Updated: jul 12 2024 (17:16) 
##           By: Brice Ozenne
##     Update #: 91
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * iid.lmm (documentation)
##' @title Extract the Influence Function From a Linear Mixed Model
##' @description Extract the influence function from a linear mixed model.
##' @rdname iid.lmm
##' @rdname influence.lmm
##' 
##' @param x a \code{lmm} object.
##' @param effects [character] Should the variance-covariance matrix for all coefficients be output (\code{"all"}),
##' or only for coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only for coefficients relative to the variance structure (\code{"variance"}),
##' or only for coefficients relative to the correlation structure (\code{"correlation"}).
##' @param p [numeric vector] value of the model coefficients at which to evaluate the influence function. Only relevant if differs from the fitted values.
##' @param robust [logical] If \code{FALSE} the influence function is rescaled to match the model-based standard errors.
##' The correlation however will not (necessarily) match the model-based correlation.
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @keywords methods

## * iid.lmm (code)
##' @export
iid.lmm <- function(x,
                    effects = "mean",
                    p = NULL,
                    robust = TRUE,
                    type.information = NULL,
                    transform.sigma = NULL,
                    transform.k = NULL,
                    transform.rho = NULL,
                    transform.names = TRUE,
                    ...){

    ## ** normalize user imput

    ## dots
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## effects
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)

    ## type.information
    if(is.null(type.information)){
        type.information <- x$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## p and transform
    init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    ## ** get information and score
    x.vcov <- vcov.lmm(x, effects = effects, p = p, robust = FALSE, type.information = type.information, df = FALSE,
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    x.score <- score.lmm(x, effects = effects, p = p, indiv = TRUE, 
                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)

    ## ** compute iid
    out <- x.score %*% x.vcov
    if(robust==FALSE){
        out <- sweep(out, MARGIN = 2, FUN = "*", STATS = sqrt(diag(x.vcov))/sqrt(colSums(out^2, na.rm = TRUE)))
    }

    ## ** export
    return(out)
}

## * iid.Wald_lmm (code)
##' @export
iid.Wald_lmm <- function(x, ...){

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** export
    return(x$iid)
}

## * iid.mlmm (code)
##' @export
iid.mlmm <- function(x,
                     effects = "contrast",
                     p = NULL,
                     robust = TRUE,
                     type.information = NULL,
                     transform.sigma = NULL,
                     transform.k = NULL,
                     transform.rho = NULL,
                     transform.names = TRUE,
                     simplify = TRUE,
                     ...){

    ## ** normalize use input
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(!is.null(effects)){
        effects <- match.arg(effects, c("contrast","mean","fixed","variance","correlation","all"), several.ok = TRUE)
    }
    if(!is.null(effects) && length(effects)==1 && effects=="contrast"){
        transform.names  <- FALSE
        if(!is.null(x$univariate) & all(x$univariate$type=="mu")){
            effects2 <- "mean"
        }else{
            effects2 <- "all"
        }
    }else{
        effects2 <- effects
    }

    ## *** type.information
    if(is.null(type.information)){
        type.information <- x$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## *** p and transform
    if(!is.null(p)){
        if(!is.list(p)){
            stop("Argument \'p\' should either be NULL or a list. \n")
        }
        if(is.null(names(p))){
            stop("Argument \'p\' should either be NULL or a named list. \n")
        }
        if(any(names(p) %in% names(x$model) == FALSE)){
            stop("Incorrect names for argument \'p\': \"",paste(setdiff(names(p),names(x$model)), collapse = "\", \""),"\". \n", 
                 "Should be among \"",paste(names(x$model), collapse = "\", \""),"\". \n")
        }
    }else{
        p <- stats::setNames(vector(mode = "list", length = length(x$model)), names(x$model))
    }

    ls.init <- lapply(names(x$model),function(iM){ ## iM <- names(x$model)[1]
        .init_transform(p = p[[iM]], transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                        x.transform.sigma = x$args$transform.sigma, x.transform.k = x$args$transform.k, x.transform.rho = x$args$transform.rho,
                        table.param = x$model[[iM]]$design$param)
    })    
    theta <- setNames(lapply(ls.init, "[[","p"),names(x$model))
    test.notransform <- unique(sapply(ls.init,"[[","test.notransform"))
    transform.sigma <- unique(sapply(ls.init,"[[","transform.sigma"))
    transform.k <- unique(sapply(ls.init,"[[","transform.k"))
    transform.rho <- unique(sapply(ls.init,"[[","transform.rho"))
    if(length(test.notransform)>1 || length(transform.sigma)>1 || length(transform.k)>1 || length(transform.rho)>1){
        stop("Something went wrong when initializing the transformation parameters. \n",
             "Not the same transformation for all models. \n")
    }

    ## ** prepare
    if(!is.null(effects) && length(effects)==1 && effects=="contrast"){

        M.contrast <- stats::coef(x, type = "contrast")
        lsM.contrast <- stats::coef(x, type = "ls.contrast")

    }

    ## ** get iid from each model
    cluster <- x$object$cluster
    n.cluster <- length(cluster)
    
    ls.iid <- lapply(names(x$model), FUN = function(iBy){ ## iBy <- names(x$model)[[1]]
        iO <- x$model[[iBy]]
        if(iO$args$method.fit == "REML" && any(c("variance", "correlation", "all") %in% effects2)){
            iO$args$method.fit <- "ML"
            message <- "iid.REML2ML" ## the influence function is computed under the ML loss function plugging in the REML estimates
        }else{
            message <- NULL
        }
        iIID <- iid(iO, p = p[[iBy]], effects = effects2, robust = robust, type.information = type.information,
                    transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        if(!is.null(effects) && length(effects)==1 && effects=="contrast"){
            iIID <- iIID %*% t(lsM.contrast[[iBy]][,colnames(iIID),drop=FALSE])

        }
        if(simplify){        
            iOut <- matrix(0, nrow = n.cluster, ncol = NCOL(iIID), dimnames = list(cluster, paste0(iBy,": ", colnames(iIID))))
        }else{
            iOut <- matrix(0, nrow = n.cluster, ncol = NCOL(iIID), dimnames = list(cluster, colnames(iIID)))
        }
        iOut[rownames(iIID),] <- iIID
        if(simplify){       
            attr(iOut,"name") <- colnames(iIID)
            attr(iOut,"by") <- rep(iBy, NCOL(iIID))
        }
        attr(iOut,"message") <- message
        return(iOut)
    })

    if(!is.null(effects) && length(effects)==1 && effects=="contrast"){
        out <- do.call(cbind,ls.iid) %*% M.contrast
        attr(out,"message") <- unique(unlist(lapply(ls.iid,attr,"message")))
    }else if(simplify){        
        out <- do.call(cbind,ls.iid)
        attr(out,"original.name") <- unlist(lapply(ls.iid,attr,"name"))
        attr(out,"by") <- unlist(lapply(ls.iid,attr,"by"))
        attr(out,"message") <- unique(unlist(lapply(ls.iid,attr,"message")))
    }else{
        out <- ls.iid
    }

    ## ** export
    return(out)
    
}

##----------------------------------------------------------------------
### iid.R ends here
