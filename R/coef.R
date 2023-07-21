### coef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:30) 
## Version: 
## Last-Updated: jul 21 2023 (17:03) 
##           By: Brice Ozenne
##     Update #: 602
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
##'
##' @param object a \code{lmm} object.
##' @param effects [character] Should all coefficients be output (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance structure (\code{"variance"}),
##' or only coefficients relative to the correlation structure (\code{"correlation"}).
##' Can also be \code{"ranef"} to output random effect (only for \code{CS} structure).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ordering [character] According to which variable the coefficient should be ordered in the output: \code{"by"} or \code{"parameter"}.
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
##' When using a (pure) compound symmetry covariance structure (\code{structure = "CS"}),
##' estimated random effects can be extracted by setting argument \code{effects} to \code{"ranef"}.
##'
##' @return A vector with the value of the model coefficients.
##' 
##' @seealso
##' \code{\link{confint.lmm}} or \code{\link{model.tables.lmm}} for a data.frame containing estimates with their uncertainty. \cr
##' 
##' @keywords methods
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
##' @export
coef.lmm <- function(object, effects = NULL, p = NULL,
                     transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE, ...){

    ## ** extract from object
    param.name <- object$design$param$name
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.level <- stats::setNames(object$design$param$level,param.name)
    param.sigma <- stats::setNames(object$design$param$sigma,param.name)
    param.k.x <- stats::setNames(object$design$param$k.x,param.name)
    param.k.y <- stats::setNames(object$design$param$k.y,param.name)

    object.reparametrize.name <- names(object$reparametrize$p)
    object.reparametrize.value <- object$reparametrize$p
    object.reparametrize.newname <- object$reparametrize$newname

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
    if("ranef" %in% effects){
        if(length(effects)>1){
            stop("Argument \'effects\' should be of length 1 when it contains \"ranef\". \n")
        }
        return(ranef(object, p = p))
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation","ranef"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    
    effects2 <- effects
    if(transform.rho == "cov"){
        if(all("correlation" %in% effects == FALSE)){
            stop("Cannot use the argument \'transform.rho\' set to \"cov\" when \"correlation\" is not in argument \'effect\'. \n")
        }
        if(all("variance" %in% effects == FALSE)){
            effects2 <- c("variance",effects2)
        }
    }
    if(!is.null(p)){
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(param.name %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(param.name[param.name %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        p <- p[param.name]
        if(object$reparametrize$transform){
            reparametrize.p <- .reparametrize(p = p[object.reparametrize.name],  
                                              type = param.type[object.reparametrize.name],
                                              sigma = param.sigma[object.reparametrize.name],
                                              k.x = param.k.x[object.reparametrize.name],
                                              k.y = param.k.y[object.reparametrize.name],
                                              level = param.level[object.reparametrize.name],                                              
                                              Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                              transform.sigma = transform.sigma,
                                              transform.k = transform.k,
                                              transform.rho = transform.rho,
                                              transform.names = FALSE)$p
        }else{
            reparametrize.p <- p[object.reparametrize.name]
        }
    }else{
        p <- object$param
        reparametrize.p <- object.reparametrize.value
    }

    ## ** extract
    out <- NULL
    if("mean" %in% effects2){
        out <- c(out, p[param.type=="mu"])
    }
    if(any(c("variance","correlation") %in% effects2)){
        pVar <- NULL
        if("variance" %in% effects2){
            if(test.notransform){                
                index.sigmak <- names(param.type)[param.type %in% c("sigma","k")]
                if(transform.names && !is.null(object.reparametrize.newname)){
                    pVar <- c(pVar, stats::setNames(reparametrize.p[index.sigmak],object.reparametrize.newname[match(index.sigmak,names(reparametrize.p))]))
                }else{
                    pVar <- c(pVar, reparametrize.p[index.sigmak])
                }                    
            }else{
                pVar <- c(pVar, p[param.name[param.type %in% c("sigma","k")]])
            }
        }
        if("correlation" %in% effects2){
            if(test.notransform){
                index.rho <- names(param.type)[param.type %in% c("rho")]
                if(transform.names && !is.null(object.reparametrize.newname)){
                    pVar <- c(pVar, stats::setNames(reparametrize.p[index.rho],object.reparametrize.newname[match(index.rho,names(reparametrize.p))]))
                }else{
                    pVar <- c(pVar, reparametrize.p[index.rho])
                }                    
            }else{
                pVar <- c(pVar, p[param.name[param.type %in% c("rho")]])
            }
        }
        if(!test.notransform){
            ls.reparam <- .reparametrize(p = pVar,  
                                         type = param.type[names(pVar)], 
                                         sigma = param.sigma[names(pVar)], 
                                         k.x = param.k.x[names(pVar)], 
                                         k.y = param.k.y[names(pVar)], 
                                         level = param.level[names(pVar)], 
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
    if(("variance" %in% effects2) && ("variance" %in% effects == FALSE)){
        index.rm <- which(names(newname) %in% param.name[param.type %in% c("sigma","k")])
        newname <- newname[-index.rm]
        out <- out[-index.rm]
    }
    if(length(newname)>0){
        ## rename
        names(out)[match(names(newname),names(out))] <- as.character(newname)
    }
    return(out)
}

## * coef.lmmCC (code)
##' @export
coef.lmmCC <- function(object, effects = NULL, ...){

    if(object$time$n==4 && (is.null(effects) || effects == "change")){

        dots <- list(...)
        if(length(dots)>0){
            stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
        }
    
        Mcon <- cbind(c(-1,1,0,0),c(0,0,-1,1))
        Sigma.change <- t(Mcon) %*% stats::sigma(object) %*% Mcon
        out <- c(cor = stats::cov2cor(Sigma.change)[1,2],
                 beta = Sigma.change[1,2]/Sigma.change[1,1])
        
    }else{

        class(object) <- setdiff(class(object),"lmmCC")
        out <- coef(object, effects = effects, ...)

    }

    ## ** export
    return(out)

}

## * coef.Wald_lmm
##' @export
coef.Wald_lmm <- function(object, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    if(is.null(object$univariate)){
        return(NULL)
    }else{
        return(stats::setNames(object$univariate$estimate, rownames(object$univariate)))
    }
    
}
## * coef.LRT_lmm
##' @export
coef.LRT_lmm <- function(object, ...){
    message("No effect size available for likelihood ratio tests.")
    return(NULL)
}

## * coef.mlmm
##' @export
coef.mlmm <- function(object, effects = "contrast", ordering = "parameter", ...){
    ordering <- match.arg(ordering, c("by","parameter"))
    
    if(!is.null(effects) && effects=="contrast"){
        table.order <- object$univariate[order(object$univariate[[ordering]]),,drop=FALSE]
        return(stats::setNames(table.order$estimate,rownames(table.order)))
    }else{
        ls.out <- lapply(object$model, coef, effects = effects, ...)
        if(ordering == "by"){
            return(ls.out)
        }else if(ordering == "parameter"){
            Uname <- unique(unlist(lapply(ls.out,names)))
            ls.out2 <- stats::setNames(lapply(Uname, function(iName){ ## iName <- "X1"
                unlist(lapply(ls.out,function(iVec){unname(iVec[iName])}))
            }), Uname)
            return(ls.out2)
        }
    }
    
}



##----------------------------------------------------------------------
### coef.R ends here
