### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  4 2021 (10:04) 
## Version: 
## Last-Updated: jul  4 2024 (12:31) 
##           By: Brice Ozenne
##     Update #: 64
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
    dots <- list(...)
    dots$complete <- NULL ## for multcomp which passes an argument complete when calling vcov
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)

    if(is.null(type.information)){
        type.information <- x$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    ## ** get information and score
    x.vcov <- vcov.lmm(x, effects = effects, p = p, robust = FALSE, type.information = type.information, df = FALSE,
                          transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
    x.score <- score.lmm(x, effects = effects, p = p, indiv = TRUE, 
                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)

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
iid.mlmm <- function(x, effects = "contrast", robust = TRUE, ...){

    ## ** normalize use input

    ## effects
    if(!is.null(effects)){
        effects <- match.arg(effects, c("contrast","mean","fixed","variance","correlation","all"), several.ok = TRUE)
    }
    if(!is.null(effects) && length(effects)==1 && effects=="contrast"){
        if(!is.null(x$univariate) & all(x$univariate$type=="mu")){
            effects2 <- "mean"
        }else{
            effects2 <- "all"
        }
    }else{
        effects2 <- effects
    }

    ## ** get information and score
    x.score <- score.mlmm(x, effects = effects2, indiv = TRUE, ordering = "by", ...)
    x.vcov <- vcov.mlmm(x, effects = effects2, robust = FALSE, ...)

    ## ** compute iid
    ls.out <- mapply(iScore = x.score, iVcov = x.vcov, function(iScore, iVcov){
        iIID <- iScore %*% iVcov
        if(robust==FALSE){
            return(sweep(iIID, MARGIN = 2, FUN = "*", STATS = sqrt(diag(iVcov))/sqrt(colSums(iIID^2, na.rm = TRUE))))
        }else{
            return(iIID)
        }
    }, SIMPLIFY = FALSE)

    ## ** export
    if(!is.null(effects) && length(effects)==1 && effects=="contrast"){

        M.contrast <- stats::coef(x, type = "contrast")
        lsM.contrast <- stats::coef(x, type = "ls.contrast")

        ls.Cout <- mapply(iIID = ls.out, iC = lsM.contrast, function(iIID,iC){ ## iIID <- ls.out[[1]] ; iC <- lsM.contrast[[1]]
            iIID %*% t(iC[,colnames(iIID),drop=FALSE])
        }, SIMPLIFY = FALSE)

        out <- do.call(cbind,ls.Cout) %*% M.contrast
        colnames(out) <- names(stats::coef(x))
        
    }else{
        out <- ls.out
    }

    return(out)
    
}

##----------------------------------------------------------------------
### iid.R ends here
