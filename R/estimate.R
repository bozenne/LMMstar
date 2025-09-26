### estimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 20 2021 (23:25) 
## Version: 
## Last-Updated: jul 25 2025 (13:38) 
##           By: Brice Ozenne
##     Update #: 1442
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * estimate.lmm (documentation)
##' @title Delta Method for a Linear Mixed Model
##' @description Estimate standard errors, confidence intervals, and p-values for a smooth transformation of parameters from a linear mixed model.
##'
##' @param x  a \code{lmm} object.
##' @param f [function] function taking as input \code{p}, the linear mixed model parameters, which outputs the parameter(s) of interest.
##' Can accept extra-arguments, such as \code{object} representing the linear mixed model.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. 
##' Can also be \code{2} compute the degrees-of-freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
##' @param df [logical] Should degree-of-freedom, computed using Satterthwaite approximation, for the parameter of interest be output.
##' Can also be a numeric vector providing providing the degrees-of-freedom relative to each estimate.
##' @param type.information [character] Should the expected information be used  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param average [logical] is the estimand the average output of argument \code{f}?
##' Otherwise consider each individual output of argument \code{f}.
##' @param method.numDeriv [character] method used to approximate the gradient: either \code{"simple"} or \code{"Richardson"}.
##' Passed to \code{numDeriv::jacobian}.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param ... extra arguments passed to \code{f}.
##'
##' @keywords mhtest
##' 
##' @details Based a first order delta method to evaluate the variance of the estimate.
##' The derivative of the transformation is evaluated using numerical differentiation (\code{numDeriv::jacobian}). \cr \cr
##' 
##' Argument \bold{robust}: the Satterhwaite approximation for the degrees-of-freedom of robust standard errors is often unreliable. This is why the default is to use the degrees-of-freedom of the modeled based standard error instead.
##' 
##' @examples
##' 
##' if(require(lava) && require(nlme)){
##' 
##' #### Random effect ####
##' set.seed(10)
##' dL <- sampleRem(1e2, n.times = 3, format = "long")
##' e.lmm1 <- lmm(Y ~ X1+X2+X3 + (1|id), repetition = ~visit|id, data = dL)
##' nlme::ranef(e.lmm1, se = TRUE)
##' e.ranef <- estimate(e.lmm1, f  = function(object, p){nlme::ranef(object, p = p)})
##' e.ranef
##'
##' if(require(ggplot2)){
##' df.gg <- cbind(index = 1:NROW(e.ranef), e.ranef)
##' gg.ranef <- ggplot(df.gg, aes(x = index, y=estimate, ymin=lower, ymax = upper))
##' gg.ranef + geom_point() + geom_errorbar() + ylab("estimated random effect") + xlab("id")
##' }
##' 
##' #### ANCOVA via mixed model ####
##' set.seed(10)
##' d <- sampleRem(1e2, n.time = 2)
##' e.ANCOVA1 <- lm(Y2~Y1+X1, data = d)
##'
##' dL2 <- reshape(d, direction = "long", idvar = c("id","X1"), 
##'                timevar = "time", times = c("1","2"), varying = c("Y1","Y2"), 
##'                v.names = "Y")
##'
##' ## estimated treatment effect (no baseline constraint)
##' e.lmm <- lmm(Y ~ time + time:X1, data = dL2, repetition = ~time|id)
##' 
##' e.delta <- estimate(e.lmm, function(p){
##'      c(Y1 = p["rho(1,2)"]*p["k.2"],
##'       X1 = p["time2:X1"]-p["k.2"]*p["rho(1,2)"]*p["time1:X1"])
##' }) ## same estimate and similar standard errors. 
##' e.delta ## Degrees-of-freedom are a bit off though
##' cbind(summary(e.ANCOVA1)$coef, df = df.residual(e.ANCOVA1))
##'
##' ## estimated treatment effect (baseline constraint)
##' dL2$time2 <- as.numeric(dL2$time=="2")
##' e.lmmC <- lmm(Y ~ time2 + time2:X1, data = dL2, repetition = ~time|id)
##' e.deltaC <- estimate(e.lmmC, function(p){
##'        c(Y1 = p["rho(1,2)"]*p["k.2"],
##'          X1 = p["time2:X1"])
##' })
##' e.deltaC ## Degrees-of-freedom are a bit more accurate
##' }

## * estimate.lmm (code)
##' @export
estimate.lmm <- function(x, f, df = !is.null(x$df) & !robust, robust = FALSE, type.information = NULL, level = 0.95,
                         method.numDeriv = NULL, average = FALSE,
                         transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, ...){


    ## ** normalize arguments
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL

    ## *** function
    f.formals <- formals(f)
    names.formals <- names(f.formals)
    is.empty <- sapply(f.formals,inherits,"name")
    if("p" %in% names.formals == FALSE){
        stop("Incorrect argument \'f\': the function should take \'p\' as argument. \n",
             "It refers to the model parameters, coef(x, effects = \"all\"), that will be varied to assess uncertainty. \n")
    }
    if(any(names(dots) %in% names.formals == FALSE)){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    ls.args <- stats::setNames(vector(mode = "list", length = length(names.formals)-1), setdiff(names.formals,"p"))
    ls.args[names(which(!is.empty))] <- f.formals[names(which(!is.empty))]
    ls.args[names(dots)] <- dots
    if("object" %in% names(ls.args) && is.null(ls.args$object)){
        ls.args <- c(list(x),ls.args[setdiff(names(ls.args),"object")])
    }

    ## *** df
    ## e.df: [numeric] value for the degrees-of-freedom (estimated or given by the user)
    ## df: [logical] whether degrees-of-freedom should be estimated
    if(length(df)>1){
        stop("Argument \'df\' must have length 1. \n")
    }
    if((!is.numeric(df) && !is.logical(df)) || !is.vector(df)){
        stop("Argument \'df\' must be a numeric or logical value. \n")
    }
    if(is.numeric(df) && df %in% 0:1){
        df <- as.logical(df)
    }
    if(is.numeric(df)){
        e.df <- df
        df <- FALSE ## no need to estimate the degrees-of-freedom
    }else if(is.logical(df)){
        if(df==TRUE){
            if(x$args$df==0){
                stop("Argument \'df\' cannot be TRUE when no degrees-of-freedom have been stored. \n",
                     "Consider setting the argument \'df\' to TRUE when calling lmm. \n")
            }
            e.df <- NULL
        }else{
            e.df <- Inf
        }
    }
    
    ## *** average
    if(!is.logical(average)){
        stop("Argument \'average\' must be TRUE or FALSE. \n")
    }
    
    ## *** derivative approximation
    if(is.null(method.numDeriv)){
        method.numDeriv <- options$method.numDeriv
    }

    ## ** warper
    if(average){
        ff <- function(p, keep.indiv = FALSE){
            iRes <- do.call(f, args = c(list(p = p),ls.args))
            iOut <- mean(iRes)
            if(keep.indiv){
                attr(iOut,"indiv") <- iRes
            }
            return(iOut)
        }
    }else{
        ff <- function(p, keep.indiv = FALSE){
            do.call(f, args = c(list(p = p),ls.args))
        }
    }

    ## ** delta-method
    out <- .estimate(x, ff = ff, average = average, df = df, level = level, 
                     robust = robust, type.information = type.information, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                     method.numDeriv = method.numDeriv, options = options)

    ## ** export
    return(out)
}

## * estimate.Wald (documentation)
##' @title Delta Method for Wald Tests
##' @description Estimate standard errors, confidence intervals, and p-values for a smooth transformation of parameters involved in Wald tests.
##'
##' @param x  a \code{Wald_lmm} object.
##' @param f [function] function taking as input \code{object}, the \code{Wald_lmm} object, which outputs the parameter(s) of interest.
##' Can accept extra-arguments.
##' @param df [logical] Should degree-of-freedom, computed using Satterthwaite approximation, for the parameter of interest be output.
##' Can also be a numeric vector providing providing the degrees-of-freedom relative to each estimate.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param method.numDeriv [character] method used to approximate the gradient: either \code{"simple"} or \code{"Richardson"}.
##' Passed to \code{numDeriv::jacobian}.
##' @param ... extra arguments passed to \code{f}.
##'
##' @keywords mhtest
##' 
##' @details Based a first order delta method to evaluate the variance of the estimate.
##' The derivative of the transformation is evaluated using numerical differentiation (\code{numDeriv::jacobian}). \cr \cr
##' 
##' Compared to \code{estimate.lmm}, the \code{f} argument cannot take \code{p} as argument.
##' One should instead explicitely request the estimated contrasts in the function \code{f} using \code{coef(object)}.

## * estimate.Wald_lmm (code)
##' @export
estimate.Wald_lmm <- function(x, f, df = x$args$df, level = 0.95, method.numDeriv = NULL, ...){

    ## ** normalize arguments
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL

    ## *** function
    f.formals <- formals(f)
    names.formals <- names(f.formals)
    is.empty <- sapply(f.formals,inherits,"name")
    if("object" %in% names.formals == FALSE){
        stop("Incorrect argument \'f\': the function should take \'object\' as argument. \n",
             "It refers to the output of anova based on model parameters that will be varied to assess uncertainty. \n")
    }
    if("p" %in% names.formals){
        stop("Incorrect argument \'f\': the function should not take \'p\' as argument. \n",
             "Consider using instead argument \'object\', refering to the output of anova based on model parameters that will be varied to assess uncertainty. \n")
    }
    if(any(names(dots) %in% names.formals == FALSE)){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    ls.args <- stats::setNames(vector(mode = "list", length = length(names.formals)-1), setdiff(names.formals,"object"))
    ls.args[names(which(!is.empty))] <- f.formals[names(which(!is.empty))]
    ls.args[names(dots)] <- dots
    
    ## *** df
    ## e.df: [numeric] value for the degrees-of-freedom (estimated or given by the user)
    ## df: [logical] whether degrees-of-freedom should be estimated
    if(length(df)>1){
        stop("Argument \'df\' must have length 1. \n")
    }
    if((!is.numeric(df) && !is.logical(df)) || !is.vector(df)){
        stop("Argument \'df\' must be a numeric or logical value. \n")
    }
    if(is.numeric(df) && df %in% 0:1){
        df <- as.logical(df)
    }
    if(is.numeric(df)){
        e.df <- df
        df <- FALSE ## no need to estimate the degrees-of-freedom
    }else if(is.logical(df)){
        if(df==TRUE){
            if(x$args$df==0){
                stop("Argument \'df\' cannot be TRUE when no degrees-of-freedom have been stored. \n",
                     "Consider setting the argument \'df\' to TRUE when calling lmm. \n")
            }
            e.df <- NULL
        }else{
            e.df <- Inf
        }
    }
    
    ## *** derivative approximation
    if(is.null(method.numDeriv)){
        method.numDeriv <- options$method.numDeriv
    }

    ## ** wraper
    ff <- function(p, keep.indiv){ ## x <- NULL
        do.call(f, args = c(list(object = stats::update(x, p = p)), ls.args))
    }

    ## ** delta-method
    out <- .estimate(x, ff = ff, average = FALSE, df = df, level = level, 
                     transform.sigma = x$args$transform.sigma, transform.k = x$args$transform.k, transform.rho = x$args$transform.rho,
                     method.numDeriv = method.numDeriv, options = options)

    ## ** export
    return(out)
}

## * estimate.rbindWald (documentation)
##' @title Delta Method for Combined Wald Tests
##' @description Estimate standard errors, confidence intervals, and p-values for a smooth transformation of parameters involved in combined Wald tests.
##'
##' @param x  a \code{rbindWald_lmm} object.
##' @param f [function] function taking as input \code{object}, the \code{rbindWald_lmm} object, which outputs the parameter(s) of interest.
##' Can accept extra-arguments.
##' @param df [logical] Should degree-of-freedom, computed using Satterthwaite approximation, for the parameter of interest be output.
##' Can also be a numeric vector providing providing the degrees-of-freedom relative to each estimate.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param method.numDeriv [character] method used to approximate the gradient: either \code{"simple"} or \code{"Richardson"}.
##' Passed to \code{numDeriv::jacobian}.
##' @param ... extra arguments passed to \code{f}.
##'
##' @keywords mhtest
##' 
##' @details Based a first order delta method to evaluate the variance of the estimate.
##' The derivative of the transformation is evaluated using numerical differentiation (\code{numDeriv::jacobian}). \cr \cr
##' 
##' Compared to \code{estimate.lmm}, the \code{f} argument cannot take \code{p} as argument.
##' One should instead explicitely request the estimated contrasts in the function \code{f} using \code{coef(object)}.
##' 
##' @examples
##' if(require(lava)){
##' 
##' data(gastricbypassL, package = "LMMstar")
##' e.reg <- by(gastricbypassL, gastricbypassL$time, function(iData){
##'    iLMM <- lmm(glucagonAUC ~ weight, data = iData, repetition = ~1|id)
##'    anova(iLMM, "weight=0", simplify = FALSE)
##' })
##'
##' e.Wald14 <- rbind(e.reg[[1]], e.reg[[2]], e.reg[[3]], e.reg[[4]])
##'
##' estimate(e.Wald14, function(object){ 
##'    p <- coef(object)
##'    c(p, average = mean(p))
##' })
##'
##' estimate(e.Wald14, function(object){ 
##'    coef(object, method = c("none","average"))
##' })
##'
##' }
##' @export
estimate.rbindWald_lmm <- estimate.Wald_lmm

## * estimate.mlmm (code)
##' @title Delta Method for Multiple Linear Mixed Models
##' @description Estimate standard errors, confidence intervals, and p-values for a smooth transformation of parameters from multiple linear mixed models.
##'
##' @param x  a \code{rbindWald_lmm} object.
##' @param f [function] function taking as input \code{object}, the \code{rbindWald_lmm} object, which outputs the parameter(s) of interest.
##' Can accept extra-arguments.
##' @param df [logical] Should degree-of-freedom, computed using Satterthwaite approximation, for the parameter of interest be output.
##' Can also be a numeric vector providing providing the degrees-of-freedom relative to each estimate.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param method.numDeriv [character] method used to approximate the gradient: either \code{"simple"} or \code{"Richardson"}.
##' Passed to \code{numDeriv::jacobian}.
##' @param ... extra arguments passed to \code{f}.
##'
##' @keywords mhtest
##' 
##' @details Based a first order delta method to evaluate the variance of the estimate.
##' The derivative of the transformation is evaluated using numerical differentiation (\code{numDeriv::jacobian}). \cr \cr
##' 
##' Compared to \code{estimate.lmm}, the \code{f} argument cannot take \code{p} as argument.
##' One should instead explicitely request the estimated contrasts in the function \code{f} using \code{coef(object)}.
##' 
##' @keywords mhtest
##' 
##' @export
estimate.mlmm <- estimate.Wald_lmm

## * .estimate
##' @description Perform delta method based on a refit function
##' @noRd 
.estimate <- function(x, ff, average, df, level, 
                      robust = NULL, type.information = NULL, transform.sigma, transform.k, transform.rho,
                      method.numDeriv, options){

    ## ** initialize transform (for lmm only, otherwise transform.sigma/k/rho are provided, i.e. not NULL)
    if(is.null(transform.sigma)){
        transform2.sigma <- "none"
    }else{
        transform2.sigma <- transform.sigma
    }
    if(is.null(transform.k)){
        transform2.k <- "none"
    }else{
        transform2.k <- transform.k
    }
    if(is.null(transform.rho)){
        transform2.rho <- "none"
    }else{
        transform2.rho <- transform.rho
    }

    ## ** estimate
    beta <- stats::coef(x, effects = "all", transform.sigma = transform2.sigma, transform.k = transform2.k, transform.rho = transform2.rho, transform.names = FALSE, simplify = FALSE, options = options)
    type.beta <- attr(beta,"type")
    attr(beta, "type") <- NULL
    attr(beta, "sigma") <- NULL
    attr(beta, "k.x") <- NULL
    attr(beta, "k.y") <- NULL
    attr(beta, "transform.sigma") <- transform2.sigma
    attr(beta, "transform.k") <- transform2.k
    attr(beta, "transform.rho") <- transform2.rho

    ## ** partial derivative
    fbeta <- ff(beta, keep.indiv = TRUE)
    if(average){
        fbeta.indiv <- attr(fbeta,"indiv")
        attr(fbeta,"indiv") <- NULL
    }
    ## remove attributes otherwise is.vector(fbeta) is FALSE
    if(!is.null(attr(fbeta,"transform.sigma"))){attr(fbeta,"transform.sigma") <- NULL}
    if(!is.null(attr(fbeta,"transform.k"))){attr(fbeta,"transform.k") <- NULL}
    if(!is.null(attr(fbeta,"transform.rho"))){attr(fbeta,"transform.rho") <- NULL}

    if(!is.vector(fbeta) || (!is.numeric(fbeta) && !is.logical(fbeta))){
        stop("The output of the function defined in the argument \'FUN\' must be a numeric vector. \n")
    }
    if(!is.null(names(fbeta)) && any(duplicated(names(fbeta)))){
        stop("The output of the function defined in the argument \'FUN\' should not contain duplicated names. \n")
    }

    grad <- numDeriv::jacobian(func = ff, x = beta, method = method.numDeriv)
    colnames(grad) <- names(beta)

    ## revert back to usual transformation if vcov parameter not used (for lmm only, otherwise transform.sigma/k/rho are provided, i.e. not NULL)
    if(all(!is.na(grad))){
        if(all(colSums(grad[,type.beta=="sigma",drop=FALSE]!=0)==0) && is.null(transform.sigma)){
            transform2.sigma <- x$reparametrize$transform.sigma
        }
        if(all(colSums(grad[,type.beta=="k",drop=FALSE]!=0)==0) && is.null(transform.k)){
            transform2.k <- x$reparametrize$transform.k
        }
        if(all(colSums(grad[,type.beta=="rho",drop=FALSE]!=0)==0) && is.null(transform.rho)){
            transform2.rho <- x$reparametrize$transform.rho
        }
    }

    ## ** extract variance-covariance
    if(is.null(robust) && is.null(type.information)){ ## Wald_lmm, rbindWald_lmm, mlmm
        Sigma <- stats::vcov(x, effects = list("all",c("all","gradient"))[[df+1]],
                             transform.sigma = transform2.sigma, transform.k = transform2.k, transform.rho = transform2.rho, options = options)
    }else{ ## lmm
        Sigma <- stats::vcov(x, effects = list("all",c("all","gradient"))[[df+1]],
                             robust = robust, type.information = type.information, 
                             transform.sigma = transform2.sigma, transform.k = transform2.k, transform.rho = transform2.rho, options = options)
    }

    ## ** delta-method
    C.Sigma.C <- grad %*% Sigma %*% t(grad)

    ## second order?
    ## g(\thetahat) = g(\theta) + (\thetahat-\theta)grad + 0.5(\thetahat-\theta)lap(\thetahat-\theta) + ...
    ## Var[g(\thetahat)] = grad\Var[\thetahat-\theta]grad + 0.25\Var[(\thetahat-\theta)lap(\thetahat-\theta)] + \Cov((\thetahat-\theta)grad,(\thetahat-\theta)lap(\thetahat-\theta)) + ...
    ## https://stats.stackexchange.com/questions/427332/variance-of-quadratic-form-for-multivariate-normal-distribution
    ## Var[g(\thetahat)] = grad\Var[\thetahat-\theta]grad + 0.25*2*tr((lap\Var[\thetahat-\theta])^2) + 0 + ...
    
    ## lap <- numDeriv::hessian(func = f, x = beta) ## laplacian
    ## 2 * sum(diag(Sigma %*% lap %*% Sigma %*% lap))
    if(average){
        C.sigma.C <- sqrt(diag(C.Sigma.C) + sum((fbeta.indiv - fbeta)^2)/(length(fbeta.indiv)-1))
    }else{
        C.sigma.C <- sqrt(diag(C.Sigma.C))
    }
    
    ## ** df
    if(df){
        colnames(grad) <- colnames(Sigma)
        e.df <- .df_contrast(contrast = grad, vcov.param = Sigma, dVcov.param = attr(Sigma, "gradient"))
    }else if(length(df)==1){
        if(identical(FALSE,df)){
            e.df <- rep(Inf, length(fbeta))
        }else{
            e.df <- rep(df,length(fbeta))
        }
    }else if(length(e.df) != length(fbeta)){
        stop("Incorrect length of argument \'df\': when a numeric vector it should have lenght either 1 or the output of argument \'f\'. \n",
             "Valid length: ",length(fbeta),". \n")
    }else{
        e.df <- df
    }

    ## ** export
    alpha <- 1-level
    out <- data.frame(estimate = as.double(fbeta),
                      se = as.double(C.sigma.C),
                      df = as.double(e.df),
                      lower = as.double(fbeta + stats::qt(alpha/2, df = e.df) * C.sigma.C),
                      upper = as.double(fbeta + stats::qt(1-alpha/2, df = e.df) * C.sigma.C),
                      p.value = as.double(2*(1-stats::pt(abs(fbeta/C.sigma.C), df = e.df))))
    attr(out,"gradient") <- grad
    if(!is.null(names(fbeta))){
        rownames(out) <- names(fbeta)
        rownames(attr(out,"gradient")) <- names(fbeta)
    }
    colnames(attr(out,"gradient")) <- names(beta)
    return(out)
}
                
##----------------------------------------------------------------------
### estimate.R ends here
