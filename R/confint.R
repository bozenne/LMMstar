### confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: okt  1 2021 (17:24) 
##           By: Brice Ozenne
##     Update #: 320
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.lmm (documentation)
##' @title Statistical Inference for Linear Mixed Model
##' @description Compute confidence intervals (CIs) and p-values for the coefficients of a multivariate gaussian model. 
##' @name confint
##' 
##' @param object a \code{lmm} object.
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param effects [character] Should the CIs/p-values for all coefficients be output (\code{"all"}),
##' or only for mean coefficients (\code{"mean"} or \code{"fixed"}),
##' or only for variance coefficients (\code{"variance"}),
##' or only for correlation coefficients (\code{"correlation"}).
##' @param robust [logical] Should robust standard error (aka sandwich estimator) be output instead of the model-based standard errors. Not feasible for variance or correlation coefficients estimated by REML.
##' @param null [numeric vector] the value of the null hypothesis relative to each coefficient.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficient. Otherwise a normal distribution is used.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param strata [character vector] When not \code{NULL}, only output coefficient relative to specific levels of the variable used to stratify the mean and covariance structure.
##' @param columns [character vector] Columns to be output. Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param type.information,transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param backtransform [logical] should the variance/covariance/correlation coefficient be backtransformed?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @seealso the function \code{anova} to perform inference about linear combinations of coefficients and adjust for multiple comparisons.
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
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL)
##' 
##' ## based on a Student's t-distribution with transformation
##' confint(eUN.lmm)
##' ## based on a Student's t-distribution without transformation
##' confint(eUN.lmm, transform.sigma = "none", transform.k = "none", transform.rho = "none")
##' ## based on a Normal distribution with transformation
##' confint(eUN.lmm, df = FALSE)
##' 

## * confint.lmm (code)
##' @export
confint.lmm <- function (object, parm = NULL, level = 0.95, effects = NULL, robust = FALSE, null = NULL,
                         type.object = "lmm", strata = NULL, columns = NULL,
                         df = NULL, type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE,
                         backtransform = NULL, ...){

    ## ** normalize user imput
    dots <- list(...)
    options <- LMMstar.options()
    
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.null(parm)){
        stop("Argument \'parm\' should not be used. It is here for compatibility with the generic method. \n",
             "Use \'effects\' instead. \n")
    }
    type.object <- match.arg(type.object, c("lmm","gls"))
    type.param <- object$param$type
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
    if(type.object=="gls"){
        if(!is.null(strata)){
            return(lapply(object$gls[strata], intervals))
        }else{
            return(lapply(object$gls, intervals))
        }
    }
    if(is.null(df)){
        df <- (!is.null(object$df)) && (robust==FALSE)
    }
    if(is.null(backtransform)){
        if(is.null(transform.sigma) && is.null(transform.k) && is.null(transform.rho)){
            backtransform <- options$backtransform.confint
        }else{
            backtransform <- FALSE
        }
    }else if(is.character(backtransform)){
        backtransform <-  eval(parse(text=backtransform))
    }
    ## used to decide on the null hypothesis of k parameters
    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    transform <- init$transform

    if(is.null(type.information)){
        type.information <- attr(object$information,"type.information")
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }
    if(!is.null(columns)){
        columns  <- match.arg(columns, c("estimate","se","statistic","df","lower","upper","null","p.value"), several.ok = TRUE)
    }else{
        columns <- options$columns.confint
    }

    ## ** get estimate
    beta <- coef(object, effects = effects, type.object = "lmm", strata = strata,
                 transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    p <- length(beta)
    nameNoTransform.beta <- names(coef(object, effects = effects, type.object = "lmm", strata = strata,
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE))
   
    ## ** get uncertainty
    vcov.beta <- vcov(object, effects = effects, df = df, type.object = "lmm", strata = strata, robust = robust,
                      type.information = type.information, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)

    if(df){
        df <- pmax(attr(vcov.beta,"df"), options$min.df)
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
    name.beta <- names(beta)
    out <- data.frame(estimate = beta, se = sqrt(diag(vcov.beta[name.beta,name.beta,drop=FALSE])),
                      statistic = as.numeric(NA), df = df[name.beta], lower = as.numeric(NA), upper = as.numeric(NA), null = null, p.value = as.numeric(NA))
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
    if(is.function(backtransform)){
        out.save <- out
        out$estimate <- do.call(backtransform, list(out.save$estimate))
        out$se  <- numDeriv::grad(func = exp, x = out.save$estimate) * out.save$se
        out$lower <- do.call(backtransform, list(out.save$lower))
        out$upper <- do.call(backtransform, list(out.save$upper))
        attr(out, "backtransform") <-  2    
    }else if(backtransform){
        out <- .backtransform(out)
    }
    out[names(out)[names(out) %in% columns == FALSE]] <- NULL
    return(out)
}

## * print.confint_lmm
##' @export
print.confint_lmm <- function(x, digit = 3, ...){
    print(as.data.frame(x), digits = digit)

    backtransform <- attr(x,"backtransform")
    if(identical(backtransform, TRUE)){
        transform <- attr(x,"transform")
        type <- attr(x,"type")

        transform2 <- unlist(transform[c("sigma","k","rho")])
        ## missing.type <- stats::setNames(c("sigma","k","rho") %in% type == FALSE,c("transform.sigma","transform.k","transform.rho"))
        ## transform2[names(missing.type)[missing.type]] <- "none"
        iType <- attr(x,"type")[rownames(x)]
        
        if(length(intersect(iType,names(transform2)))>0){
            if(all(c("estimate","se","lower","upper") %in% names(x))){
                cat("Note: estimates and confidence intervals for ",paste(intersect(iType,names(transform2)), collapse = ", ")," have been back-transformed. \n",
                    "      standard errors are not back-transformed.\n", sep="")
            }else if(all(c("estimate","lower","upper") %in% names(x))){
                cat("Note: estimates and confidence intervals for ",paste(intersect(iType,names(transform2)), collapse = ", ")," have been back-transformed. \n", sep="")
            }else if(all(c("estimate","se") %in% names(x))){
                cat("Note: estimates for ",paste(intersect(iType,names(transform2)), collapse = ", ")," have been back-transformed. \n",
                    "      standard errors are not back-transformed.\n", sep="")
            }else if("estimate" %in% names(x)){
                cat("Note: estimates for ",paste(intersect(iType,names(transform2)), collapse = ", ")," have been back-transformed. \n")
            }
        }
    }else if(identical(backtransform, 2)){
        txt <- unique(c("estimates","standard errors","confidence intervals","confidence intervals")[c("estimate","se","lower","upper") %in% names(x)])
        cat("Note: ",paste(txt,collapse = ", ")," have been back-transformed. \n", sep ="")
    }
    return(NULL)
}

## * backtransform
##' @title BackTransformation for Outputs from Linear Mixed Models
##' @description Back-transform estimates and confidence intervals (CIs).
##' @noRd
##' 
##' @param object a \code{confint_lmm} object, i.e. the output of the confint function applied to a \code{lmm} object.
##'
##' @details if the option \code{transform.sigma} and/or  \code{transform.k} is one of \code{"log"}, \code{"logsd"}, \code{"logvar"}, \code{"logsqaure"},
##' the estimate and CIs are transformed back to the original scale by applying the exponential function.
##' 
##' If the option \code{transform.rho} is \code{"atanh"}, the estimate and CIs are transformed back to the original scale by applying the tangent hyperbolic function.
##'
##' @keywords internal
.backtransform <-   function(object,...) UseMethod(".backtransform")

.backtransform.confint_lmm <- function(object, ...){

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
