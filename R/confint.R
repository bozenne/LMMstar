### confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: May 15 2021 (19:04) 
##           By: Brice Ozenne
##     Update #: 71
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.lmm (documentation)
##' @title Statistical Inference for Linear Mixed Model.
##' @description Compute confidence intervals (CIs) and p-values for the coefficients of a linear mixed model.
##' @name confint
##' 
##' @param object a \code{lmm} object.
##' @param param Not used. For compatibility with the generic method.
##' @param effects [character] Should the CIs/p-values for all coefficients be output (\code{"all"}),
##' or only for mean coefficients (\code{"mean"}),
##' or only for variance coefficients (\code{"variance"}),
##' or only for correlation coefficients (\code{"correlation"}).
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param null [numeric vector] the value of the null hypothesis relative to each coefficient.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficient. Otherwise a normal distribution is used.
##' @param type.object [character] Set this argument to \code{"gls"} to obtain the output from the gls object and related methods.
##' @param strata [character vector] When not \code{NULL}, only output coefficient relative to specific levels of the variable used to stratify the mean and covariance structure.
##' @param type.information,transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method.
##' @param backtransform.mu,backtransform.sigma,backtransform.k,backtransform.rho [function] possible backtransformation for the estimate, lower and upper bounds of the confidence interval.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @seealso the function \code{multcomp::glht} to perform inference about linear combinations of coefficients and adjust for multiple comparisons.
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

## * confint.lmm (code)
##' @export
confint.lmm <- function (object, parm = NULL, effects = "all", level = 0.95, null = NULL, type.object = "lmm", strata = NULL, 
                         df = !is.null(object$df), type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE,
                         backtransform.mu = NULL, backtransform.sigma = NULL, backtransform.k = NULL, backtransform.rho = NULL, ...){

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
    if(is.null(transform.k)){
        transform.k <- options$transform.k
    }else{ 
        transform.k <- match.arg(transform.k, c("none","log","square","logsquare","sd","logsd","var","logvar"))
    }
    if(is.null(type.information)){
        type.information <- options$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## ** get estimate
    beta <- coef(object, effects = effects, type.object = "lmm", strata = strata,
                 transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    p <- length(beta)
    nameNoTransform.beta <- names(coef(object, effects = effects, type.object = "lmm", strata = strata))
    
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
        df <- setNames(rep(Inf,p),names(beta))
    }
    
    ## ** get null
    if(is.null(null)){
        if(transform.k %in% c("sd","logsd","var","logvar")){
            null <- sapply(object$param$type[nameNoTransform.beta],switch,
                           "mu" = 0,
                           "sigma" = NA,
                           "k" = NA,
                           "rho" = 0)
        }else{
            null <- sapply(object$param$type[nameNoTransform.beta],switch,
                           "mu" = 0,
                           "sigma" = NA,
                           "k" = 1,
                           "rho" = 0)
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
        null <- null[names(beta)]
    }

    ## ** combine
    out <- data.frame(estimate = beta, se = sqrt(diag(vcov.beta)), statistic = as.numeric(NA), df = df, lower = as.numeric(NA), upper = as.numeric(NA), null = null, p.value = as.numeric(NA))
    out$statistic <- (out$estimate-null)/out$se
    out$p.value <- 2*(1-stats::pt(abs(out$statistic), df = out$df))

    alpha <- 1-level
    out$lower <- out$estimate + stats::qt(alpha/2, df = out$df) * out$se
    out$upper <- out$estimate + stats::qt(1-alpha/2, df = out$df) * out$se

    ## ** back-transform
    backtransform <- list(mu = backtransform.mu,
                          sigma = backtransform.sigma,
                          k = backtransform.k,
                          rho = backtransform.rho)

    for(iType in names(backtransform)){

        if(!is.null(backtransform[[iType]])){
            if(identical("exp",backtransform[[iType]])){
                backtransform[[iType]] <- "exp"
            }else if(identical("log",backtransform[[iType]])){
                backtransform[[iType]] <- "log"
            }else if(identical("tanh",backtransform[[iType]])){
                backtransform[[iType]] <- "tanh"
            }else if(identical("atanh",backtransform[[iType]])){
                backtransform[[iType]] <- "atanh"
            }else if(is.function(backtransform.sigma)==FALSE){
                stop("Argument \'backtransform.sigma\' must be a function. \n")
            }
            iIndex.type <- which(object$param$type[nameNoTransform.beta]==iType)
            out$estimate[iIndex.type] <- sapply(out$estimate[iIndex.type],backtransform.sigma)
n            ## out$se[iIndex.type] <- sapply(iIndex.type,function(iBeta){numDeriv::jacobian(func = backtransform.sigma, x = out$estimate[iBeta])*out$se[iBeta]})
            out$lower[iIndex.type] <- sapply(out$lower[iIndex.type],backtransform.sigma)
            out$upper[iIndex.type] <- sapply(out$upper[iIndex.type],backtransform.sigma)
        }
        
    }
    
    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### confint.R ends here
