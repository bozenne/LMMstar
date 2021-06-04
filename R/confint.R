### confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: Jun  4 2021 (09:32) 
##           By: Brice Ozenne
##     Update #: 231
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.lmm (documentation)
##' @title Statistical Inference for Multivariate Gaussian Models.
##' @description Compute confidence intervals (CIs) and p-values for the coefficients of a multivariate gaussian model.
##' @name confint
##' 
##' @param object a \code{lmm} object.
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param effects [character] Should the CIs/p-values for all coefficients be output (\code{"all"}),
##' or only for mean coefficients (\code{"mean"}),
##' or only for variance coefficients (\code{"variance"}),
##' or only for correlation coefficients (\code{"correlation"}).
##' @param null [numeric vector] the value of the null hypothesis relative to each coefficient.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficient. Otherwise a normal distribution is used.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param strata [character vector] When not \code{NULL}, only output coefficient relative to specific levels of the variable used to stratify the mean and covariance structure.
##' @param type.information,transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @seealso the function \code{multcomp::glht} to perform inference about linear combinations of coefficients and adjust for multiple comparisons.
##'
##' 
##' @return A data.frame containing for each coefficient (in rows): \itemize{
##' \item column estimate: the estimate.
##' \item column se: the standard error.
##' \item column statistic: the test statistic.
##' \item column df: the degree of freedom.
##' \item column lower: the lower bound of the confidence interval.
##' \item column upper: the upper bound of the confidence interval.
##' \item column null: the null hypothesis.
##' \item column p.value: the p-value relative to the null hypothesis.
##' }
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Multivariate Gaussian Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##' 
##' ## based on normal distribution with transformation
##' confint(eUN.lmm)
##' ## based on normal distribution without transformation
##' confint(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none")
##' ## based on Student's t-distribution with transformation
##' \dontrun{
##' confint(eUN.lmm, df = TRUE)
##' }
##' 

## * confint.lmm (code)
##' @export
confint.lmm <- function (object, parm = NULL, level = 0.95, effects = "all", null = NULL, type.object = "lmm", strata = NULL, 
                         df = !is.null(object$df), type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE,
                          ...){

    options <- LMMstar.options()

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.null(parm)){
        stop("Argument \'parm\' should not be used. It is here for compatibility with the generic method. \n",
             "Use \'effects\' instead. \n")
    }
    type.object <- match.arg(type.object, c("lmm","gls"))
    type.param <- object$param$type
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)
    if(!is.null(strata)){
        strata <- match.arg(strata, object$strata$levels, several.ok = TRUE)
    }
    if(type.object=="gls"){
        if(!is.null(strata)){
            return(lapply(object$gls[strata], intervals))
        }else{
            return(lapply(object$gls, intervals))
        }
    }
    ## used to decide on the null hypothesis of k parameters
    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, options = options,
                            x.transform.sigma = NULL, x.transform.k = NULL, x.transform.rho = NULL)
    
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    transform <- init$transform

    if(is.null(type.information)){
        type.information <- options$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## ** get estimate
    beta <- coef(object, effects = effects, type.object = "lmm", strata = strata,
                 transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    p <- length(beta)
    nameNoTransform.beta <- names(coef(object, effects = effects, type.object = "lmm", strata = strata,
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE))
   
    ## ** get uncertainty 
    vcov.beta <- vcov(object, effects = effects, df = df, type.object = "lmm", strata = strata,
                      type.information = type.information, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)

    if(df){
        df <- attr(vcov.beta,"df")
        attr(vcov.beta,"df") <- NULL
        if((type.information != "observed") && ("mean" %in% effects)){
            warning("when using REML with expected information, the degree of freedom of the mean parameters may depend on the parametrisation of the variance parameters. \n")
        }
    }else{
        df <- stats::setNames(rep(Inf,p),names(beta))
    }
    
    ## ** get null
    if(is.null(null)){
        if(transform.k %in% c("sd","logsd","var","logvar")){
            null <- stats::setNames(sapply(type.param[nameNoTransform.beta],switch,
                           "mu" = 0,
                           "sigma" = NA,
                           "k" = NA,
                           "rho" = 0),nameNoTransform.beta)
        }else{
            null <- stats::setNames(sapply(type.param[nameNoTransform.beta],switch,
                           "mu" = 0,
                           "sigma" = NA,
                           "k" = 1,
                           "rho" = 0), nameNoTransform.beta)
        }
               
    }else{
        if(length(null)!=p){
            stop("Incorrect argument \'null\': there should be exactly one null hypothesis per coefficients. \n",
                 "Number of coefficients: ",p,"\n",
                 "Number of null hypothesis: ",length(null),"\n")
        }
        if(!is.null(names(null)) && any(sort(names(null)) != sort(names(beta)))){
            stop("Incorrect argument \'null\': when specified, the names of the null hypotheses should match the name of the coefficients. \n",
                 "Incorrect name(s): \"",paste(names(null)[names(null) %in% names(beta) == FALSE], collapse ="\" \""),"\" \n",
                 "Missing name(s): \"",paste(names(beta)[names(beta) %in% names(null) == FALSE], collapse ="\" \""),"\" \n")
        }
        null <- null[nameNoTransform.beta]
    }

    ## ** combine
    out <- data.frame(estimate = beta, se = sqrt(diag(vcov.beta)), statistic = as.numeric(NA), df = df, lower = as.numeric(NA), upper = as.numeric(NA), null = null, p.value = as.numeric(NA))
    out$statistic <- (out$estimate-null)/out$se
    out$p.value <- 2*(1-stats::pt(abs(out$statistic), df = out$df))

    alpha <- 1-level
    out$lower <- out$estimate + stats::qt(alpha/2, df = out$df) * out$se
    out$upper <- out$estimate + stats::qt(1-alpha/2, df = out$df) * out$se

    
    ## ** export
    attr(out, "transform") <- list(sigma = transform.sigma,
                                   k = transform.k,
                                   rho = transform.rho)
    attr(out, "type") <- type.param
    if(transform.k %in% c("sd","var","logsd","logvar")){
        attr(out, "type")[attr(out, "type")=="sigma"] <- "k"
    }
    attr(out, "old2new") <-  stats::setNames(nameNoTransform.beta, rownames(out))
    attr(out, "backtransform.names") <- names(coef(object, effects = effects, type.object = "lmm", strata = strata,
                                                   transform.sigma = gsub("log","",transform.sigma), transform.k = gsub("log","",transform.k), transform.rho = gsub("atanh","",transform.rho), transform.names = transform.names))

    attr(out, "backtransform") <-  FALSE
    class(out) <- append("confint_lmm", class(out))
    return(out)
}

## * print.confint_lmm
##' @export
print.confint_lmm <- function(x, ...){

    print(as.data.frame(x))

    backtransform <- attr(x,"backtransform")
    if(identical(backtransform, TRUE)){
        transform <- attr(x,"transform")
        type <- attr(x,"type")

        missing.type <- stats::setNames(c("sigma","k","rho") %in% type == FALSE,c("transform.sigma","transform.k","transform.rho"))

        transform2 <- unlist(transform[c("transform.sigma","transform.k","transform.rho")])
        transform2[names(missing.type)[missing.type]] <- "none"

        cat("Note: estimates and confidence intervals for ",paste(names(transform)[transform!="none"], collapse = ", ")," have been back-transformed. \n",
            "      standard errors are not back-transformed.\n", sep="")
    }
    return(NULL)
}

## * backtransform.confint_lmm
##' @title BackTransformation for Outputs from Multivariate Gaussian Models.
##' @description Back-transform estimates and confidence intervals (CIs).
##' @name confint
##' 
##' @param object a \code{confint_lmm} object, i.e. the output of the confint function applied to a \code{lmm} object.
##'
##' @details if the option \code{transform.sigma} and/or  \code{transform.k} is one of \code{"log"}, \code{"logsd"}, \code{"logvar"}, \code{"logsqaure"},
##' the estimate and CIs are transformed back to the original scale by applying the exponential function.
##' 
##' If the option \code{transform.rho} is \code{"atanh"}, the estimate and CIs are transformed back to the original scale by applying the tangent hyperbolic function.
##' 
##' @export
backtransform <-   function(object,...) UseMethod("backtransform")

##' @rdname confint
##' @export
backtransform.confint_lmm <- function(object, ...){

    type.param <- attr(object,"type")
    transform <- attr(object,"transform")
    old2new <- attr(object,"old2new")
    if(attr(object,"backtransform")){
        message("Estimates and confidence intervals have already been backtransformed.")
        return(object)
    }

    for(iType in names(transform)){ ## iType <- names(transform)[1]

        if(transform[[iType]] %in%  c("log","logsd","logvar","logsquare")){
            iBacktransform <- exp
        }else if(identical("exp",transform[[iType]])){
            iBacktransform <- log
        }else if(identical("tanh",transform[[iType]])){
            iBacktransform <- atanh
        }else if(identical("atanh",transform[[iType]])){
            iBacktransform <- tanh
        }else if(transform[[iType]] %in% c("none","sd","var","square","cov")){
            next
        }else{
            stop("Unknown transformation. \n")
        }
        iIndex.type <- which(type.param[old2new]==iType)
        ## ** back-transform
        if(length(iIndex.type)>0){
            object$estimate[iIndex.type] <- sapply(object$estimate[iIndex.type],iBacktransform)
            object$lower[iIndex.type] <- sapply(object$lower[iIndex.type],iBacktransform)
            object$upper[iIndex.type] <- sapply(object$upper[iIndex.type],iBacktransform)
        }
    }
        

    ## ** rename
    backtransform.names <- attr(object,"backtransform.names")
    rownames(object) <- backtransform.names
    
    ## ** export
    attr(object, "backtransform") <-  TRUE    
    return(object)
}
##----------------------------------------------------------------------
### confint.R ends here
