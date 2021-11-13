### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: nov 13 2021 (17:01) 
##           By: Brice Ozenne
##     Update #: 261
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * coef.lmm (documentation)
##' @title Extract Coefficients From a Linear Mixed Model
##' @description Extract coefficients from a linear mixed model.
##' @name coef
##'
##' @param object a \code{lmm} object.
##' @param effects [character] Should all coefficients be output (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' @param strata [character vector] When not \code{NULL}, only output coefficient relative to specific levels of the variable used to stratify the mean and covariance structure.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##' 
##'
##' @details \bold{transform.sigma}: \cr
##' \itemize{
##' \item \code{"none"} ouput residual standard error.
##' \item \code{"log"} ouput log-transformed residual standard error.
##' \item \code{"square"} ouput residual variance.
##' \item \code{"logsquare"} ouput log-transformed residual variance.
##' }
##'
##'  \bold{transform.k}: \cr
##' \itemize{
##' \item \code{"none"} ouput ratio between the residual standard error of the current level and the reference level.
##' \item \code{"log"} ouput log-transformed ratio between the residual standard errors.
##' \item \code{"square"} ouput ratio between the residual variances.
##' \item \code{"logsquare"} ouput log-transformed ratio between the residual variances.
##' \item \code{"sd"} ouput residual standard error of the current level.
##' \item \code{"logsd"} ouput residual log-transformed standard error of the current level.
##' \item \code{"var"} ouput residual variance of the current level.
##' \item \code{"logvar"} ouput residual log-transformed variance of the current level.
##' }
##' 
##'  \bold{transform.rho}: \cr
##' \itemize{
##' \item \code{"none"} ouput correlation coefficient.
##' \item \code{"atanh"} ouput correlation coefficient after tangent hyperbolic transformation.
##' \item \code{"cov"} ouput covariance coefficient.
##' }
##'
##' @return A vector with the value of the model coefficients.
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit linear mixed model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' ## output coefficients
##' coef(eUN.lmm)
##' coef(eUN.lmm, effects = "mean")
##' coef(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none")

## * coef.lmm (code)
##' @rdname coef
##' @export
coef.lmm <- function(object, effects = NULL, strata = NULL,
                     transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE, ...){

    ## ** normalize user imput
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    if(is.null(effects)){
        if(transform.sigma == "none" && transform.k == "none" && transform.rho == "none"){
            effects <- options$effects
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    
    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    if(transform.rho == "cov" && ("variance" %in% effects == FALSE || "correlation" %in% effects == FALSE)){
        stop("Cannot use the argument \'transform.rho\' set to \"cov\" when \"variance\" or \"correlation\" is not in argument \'effect\'. \n")
    }
    
    ## ** extract
    out <- NULL
    if("mean" %in% effects){
        out <- c(out, object$param$value[object$param$type=="mu"])
    }

    if(any(c("variance","correlation") %in% effects)){
        pVar <- NULL
        if("variance" %in% effects){
            if(test.notransform){
                index.sigmak <- names(object$param$type)[object$param$type %in% c("sigma","k")]
                if(transform.names && !is.null(object$reparametrize$newname)){
                    pVar <- c(pVar, stats::setNames(object$reparametrize$p[index.sigmak],object$reparametrize$newname[match(index.sigmak,names(object$reparametrize$p))]))
                }else{
                    pVar <- c(pVar, object$reparametrize$p[index.sigmak])
                }                    
            }else{
                pVar <- c(pVar, object$param$value[object$param$type %in% c("sigma","k")])
            }
        }
        if("correlation" %in% effects){
            if(test.notransform){
                index.rho <- names(object$param$type)[object$param$type %in% c("rho")]
                if(transform.names && !is.null(object$reparametrize$newname)){
                    pVar <- c(pVar, stats::setNames(object$reparametrize$p[index.rho],object$reparametrize$newname[match(index.rho,names(object$reparametrize$p))]))
                }else{
                    pVar <- c(pVar, object$reparametrize$p[index.rho])
                }                    
            }else{
                pVar <- c(pVar, object$param$value[object$param$type %in% c("rho")])
            }
        }
        if(!test.notransform){
            ls.reparam <- .reparametrize(p = pVar,  
                                         type = object$param$type[names(pVar)], strata = object$param$strata[names(pVar)],
                                         time.k = object$design$param$time.k, time.rho = object$design$param$time.rho,
                                         name2sd = stats::setNames(object$design$vcov$param$name2,object$design$vcov$param$name),
                                         Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
            outVar <- ls.reparam$p
            if(ls.reparam$transform){
                newname <- stats::setNames(ls.reparam$newname,names(pVar))
            }else{
                newname <- NULL
            }
        }else{
            outVar <- pVar
            newname <- NULL
        }
        out <- c(out,outVar)

        }else{
            newname <- NULL
        }
    
    ## ** post process
    ## re-order values when converting to sd with strata (avoid sd0:0 sd0:1 sd1:0 sd1:1 sd2:0 sd2:1 ...)
    if("variance" %in% effects && transform.k %in% c("sd","var","logsd","logvar") && object$strata$n>1 && transform.names){
        out.strata <- object$param$strata[names(pVar)]
        out.type <- object$param$type[names(pVar)]
        index.sd <- which(out.type %in% c("sigma","k"))
        savenames <- names(out)
        savenames[savenames %in% names(out.type)[index.sd]] <- names(out.type)[index.sd][order(out.strata[index.sd])]
        out <- out[savenames]
        
    }
    if(!is.null(strata)){
        out <- out[object$param$strata[names(out)] %in% strata]
    }
    if(length(newname)>0){
        names(out)[match(names(newname),names(out))] <- as.character(newname)
    }
    return(out)
}
##----------------------------------------------------------------------
### coef.R ends here
