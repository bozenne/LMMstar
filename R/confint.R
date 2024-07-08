### confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: jul  5 2024 (11:03) 
##           By: Brice Ozenne
##     Update #: 926
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.effect_lmm
##' @export
confint.effect_lmm <- function(object, parm, level = 0.95, method = "none", ...){
    return(confint.Wald_lmm(object, parm, level = 0.95, method = method,  ...))
}


## * confint.lmm (documentation)
##' @title Confidence Intervals for Linear Mixed Model
##' @description Compute confidence intervals (CIs) and p-values for the coefficients of a linear mixed model. 
##' 
##' @param object a \code{lmm} object.
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param effects [character] Should the CIs/p-values for all coefficients be output (\code{"all"}),
##' or only for mean coefficients (\code{"mean"} or \code{"fixed"}),
##' or only for variance coefficients (\code{"variance"}),
##' or only for correlation coefficients (\code{"correlation"}).
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. Not feasible for variance or correlation coefficients estimated by REML.
##' @param null [numeric vector] the value of the null hypothesis relative to each coefficient.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficient. Otherwise a normal distribution is used.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param type.information,transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param backtransform [logical] should the variance/covariance/correlation coefficient be backtransformed?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @seealso the function \code{anova} to perform inference about linear combinations of coefficients and adjust for multiple comparisons.
##' 
##' @return A data.frame containing some of the following coefficient (in rows): \itemize{
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
##' @seealso
##' \code{\link{coef.lmm}} for a simpler output (e.g. only estimates). \cr
##' \code{\link{model.tables.lmm}} for a more detailed output (e.g. with p-value). \cr
##'
##' @keywords methods
##' 
##' @examples
##' #### simulate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' #### fit Linear Mixed Model ####
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL)
##'
##' #### Confidence intervals ####
##' ## based on a Student's t-distribution with transformation
##' confint(eUN.lmm, effects = "all")
##' ## based on a Student's t-distribution without transformation
##' confint(eUN.lmm, effects = "all",
##'         transform.sigma = "none", transform.k = "none", transform.rho = "none")
##' ## based on a Student's t-distribution transformation but not backtransformed
##' confint(eUN.lmm, effects = "all", backtransform = FALSE)
##' ## based on a Normal distribution with transformation
##' confint(eUN.lmm, df = FALSE)
##' 

## * confint.lmm (code)
##' @export
confint.lmm <- function (object, parm = NULL, level = 0.95, effects = NULL, robust = FALSE, null = NULL,
                         columns = NULL,
                         df = NULL, type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE,
                         backtransform = NULL, ...){


    ## ** extract from object
    name.param <- object$design$param$name
    type.param <- stats::setNames(object$design$param$type, name.param)

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
    if(is.null(effects)){
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"
    if(is.null(df)){
        df <- (!is.null(object$df)) && (robust==FALSE)
    }
    if(is.null(backtransform)){
        backtransform <- rep(as.logical(options$backtransform.confint),4)
        if(!is.null(transform.sigma)){
            backtransform[1] <- FALSE
        }
        if(!is.null(transform.k)){
            backtransform[2] <- FALSE
        }
        if(!is.null(transform.rho)){
            backtransform[3] <- FALSE
        }
    }else if(any(is.character(backtransform))){
        backtransform <-  eval(parse(text=backtransform))
    }else if(all(is.numeric(backtransform))){
        backtransform <- as.logical(backtransform)
    }

    ## used to decide on the null hypothesis of k parameters
    init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    transform <- init$transform

    if(is.null(type.information)){
        type.information <- object$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }
    valid.columns <- c("estimate","se","statistic","df","lower","upper","null","p.value")
    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(!is.null(columns)){
        columns <- tolower(columns)
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse = "\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse = "\" \""),"\"\n")
        }
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            columns <- union(options$columns.confint, unname(columns))
        }
        if(!is.null(names(columns)) && all(names(columns)=="remove")){
            columns <- setdiff(options$columns.confint, unname(columns))
        }
    }else{
        columns <- options$columns.confint
    }

    ## ** get estimate
    beta <- coef(object, effects = effects, 
                 transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    p <- length(beta)
    nameNoTransform.beta <- names(coef(object, effects = effects, 
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE))
    name.beta <- names(beta)
    type.beta <- type.param[name.beta]

    ## ** get uncertainty
    vcov.beta <- vcov(object, effects = effects, df = df, robust = robust,
                      type.information = type.information, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
    if(df){
        df <- pmax(attr(vcov.beta,"df"), options$min.df)
        attr(vcov.beta,"df") <- NULL
        if((object$args$method.fit=="REML") && (type.information != "observed") && ("mean" %in% effects)){
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
        }else if(transform.k %in% c("log")){
            null <- stats::setNames(sapply(type.param[nameNoTransform.beta],switch,
                           "mu" = 0,
                           "sigma" = NA,
                           "k" = 0,
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
                      statistic = as.numeric(NA), df = df[name.beta], lower = as.numeric(NA), upper = as.numeric(NA), null = null, p.value = as.numeric(NA),
                      stringsAsFactors = FALSE)

    out$statistic <- (out$estimate-null)/out$se
    out$p.value <- 2*(1-stats::pt(abs(out$statistic), df = out$df))
    index.cor <- setdiff(which(type.beta=="mu"), which(name.beta=="(Intercept)"))
    
    alpha <- 1-level
    out$lower <- out$estimate + stats::qt(alpha/2, df = out$df) * out$se
    out$upper <- out$estimate + stats::qt(1-alpha/2, df = out$df) * out$se

    ## ** back-transform
    if(!identical(backtransform,FALSE) && !identical(backtransform,c(FALSE,FALSE,FALSE))){

        if(is.function(backtransform) || all(is.character(backtransform))){

            out <- .backtransform(out, type.param = type.param[match(nameNoTransform.beta, names(type.param))],
                                  backtransform = TRUE, backtransform.names = names(beta),
                                  transform.mu = backtransform,
                                  transform.sigma = backtransform,
                                  transform.k = backtransform,
                                  transform.rho = backtransform)

        }else{
            backtransform.names <- names(coef(object, effects = effects, 
                                              transform.sigma = gsub("log","",transform.sigma), transform.k = gsub("log","",transform.k), transform.rho = gsub("atanh","",transform.rho), transform.names = transform.names))

            out <- .backtransform(out,
                                  type.param = type.param[match(nameNoTransform.beta, names(type.param))],
                                  backtransform = backtransform, backtransform.names = backtransform.names,
                                  transform.mu = "none",
                                  transform.sigma = transform.sigma,
                                  transform.k = transform.k,
                                  transform.rho = transform.rho)            
        }
    }

    ## ** export
    out[names(out)[names(out) %in% columns == FALSE]] <- NULL
    class(out) <- append("confint_lmm", class(out))
    return(out)
}

## * confint.lmmCC (code)
##' @export
confint.lmmCC <- function(object, parm = NULL, level = 0.95, effects = NULL, columns = NULL, ...){

    if(object$time$n==4 && (is.null(effects) || effects == "change")){

        Mcon <- cbind(c(-1,1,0,0),c(0,0,-1,1))
        out.estimate <- estimate(object, function(p){
            Sigma.change <- t(Mcon) %*% sigma(object, p = p) %*% Mcon
            c(cor = stats::cov2cor(Sigma.change)[1,2],
              beta = Sigma.change[1,2]/Sigma.change[1,1])
        }, level = level, ...)
        if(!is.null(columns)){
            out <- out.estimate[,intersect(names(out.estimate),columns)]
        }else{
            out <- out.estimate[,c("estimate","lower","upper")]
        }
        
    }else{
        class(object) <- setdiff(class(object),"lmmCC")
        out <- confint(object, parm = parm, effects = effects, level = level, columns = columns, ...)
    }

    ## ** export
    return(out)

}

## * confint.mlmm (documentation)
##' @title Confidence Intervals for Multiple Linear Mixed Model.
##' @description Compute confidence intervals for several linear mixed models.
##' 
##' @param object an \code{mlmm} object, output of \code{mlmm}.
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}, or \code{"pool"}.
##' @param ordering [character] should the output be ordered by type of parameter (\code{parameter}) or by model (\code{by}).
##' Only relevant for \code{mlmm} objects.n
##' @param ... other arguments are passed to \code{\link{confint.Wald_lmm}}.
##'
##' @details Statistical inference following pooling is performed according to Rubin's rule whose validity requires the congeniality condition of Meng (1994).
##'
##' @references
##' Meng X. L.(1994). Multiple-imputation inferences with uncongenial sources of input. Statist. Sci.9, 538–58.
##' 
##' @keywords methods

## * confint.mlmm (code)
##' @export
confint.mlmm <- function(object, parm = NULL, level = 0.95, method = NULL, ordering = "parameter", ...){

    ## ** normalize user input
    if(is.null(method)){
        method <- "none"
    }
    ordering <- match.arg(ordering, c("by","parameter"))
    table.transform <- attr(object$confint.nocontrast,"backtransform")
    options <- LMMstar.options()
    pool.method <- options$pool.method

    ## ** extract confidence intervals
    out.confint <- confint.Wald_lmm(object, parm = parm, level = level, method = method, backtransform = object$args$backtransform, ...)
    if(all(method %in% pool.method == FALSE)){
        if(ordering=="by"){
            reorder <- order(object$univariate$by)
        }else if(is.list(object$univariate$parameter)){
            reorder <- order(object$univariate$type,sapply(object$univariate$parameter, paste, collapse = ";"))
        }else{
            reorder <- order(object$univariate$type,object$univariate$parameter)
        }
        out <- out.confint[reorder,,drop=FALSE]
    }else{
        out <- out.confint
    }

    ## ** export
    return(out)
}

## * confint.LRT_lmm
##' @export
confint.LRT_lmm <- function(object, parm, level = 0.95, ...){
    message("No confidence interval available for likelihood ratio tests.")
    return(NULL)
}

## * confint.resample (documentation)
##' @title Resampling Confidence Intervals for Linear Mixed Model
##' @description Compute confidence intervals for linear mixed model using resampling (permutation or bootstrap).
##'
##' @param object a \code{reample} object.
##' @param parm Not used. For compatibility with the generic method.
##' @param null [numeric vector] the value of the null hypothesis relative to each coefficient.
##' Only relevant for when using bootstrap.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param method [character] method used to compute the confidence intervals and p-values: \code{"percentile"}, \code{"gaussian"}, or \code{"studentized"}.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"sample.estimate"}, \code{"se"}, \code{"sample.se"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param correction [logical] correction to ensure non-0 p-values when using the percentile method,
##' e.g. with permutations the p.value is evaluated as (#more extreme + 1)/(n.sample + 1) instead of (#more extreme)/(n.sample).
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Argument \bold{correction}: if we denote by \code{n.sample} the number of permutations that have been performed,
##' having the correction avoids p-values of 0 by ensuring they are at least \code{as.numeric(correction)/(n.sample+as.numeric(correction))},
##' e.g. at least 0.000999001 for a thousand permutations. \cr \cr
##'
##' Argument \bold{correction} (bootstrap): percentile confidence intervals are computed based on the quantile of the resampling distribution. \cr
##' Gaussian confidence intervals and p-values are computed by assuming a normal distribution and estimating its variance based on the bootstrap distribution. \cr
##' Studentized confidence intervals are computed using quantiles based on the boostrap distribution after it has been centered and rescaled by the point estimate and its standard error
##' Percentile and Studentized p-values are computed by finding the coverage level such that one of the bound of the confidence interval equals the null. \cr \cr
##' 
##' Argument \bold{correction} (bootstrap): 
##' Percentile p-values are computed based on the proportion of times the estimate was more extreme than the permutation estimates. \cr
##' Gaussian p-values are computed by assuming a normal distribution and estimating its variance based on the permutation distribution. \cr
##' Studentized p-values are computed based on the proportion of times the test statistic was more extreme than the permutation test statistic. \cr
##' No confidence intervals are provided.
##' 

## * confint.resample (code)
## '@export
confint.resample <-  function(object, parm = NULL, null = 0, level = 0.95, method = NULL, columns = NULL, correction = TRUE, ...){

    ## ** normalize user input
    if(!missing(parm) && !is.null(parm)){
        stop("Argument \'parm\' is not used - only there for compatibility with the generic method. \n")
    }
    options <- LMMstar.options()
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    alpha <- 1-level
    type <- object$args$type
    studentized <- object$args$studentized

    valid.method <- c("percentile","gaussian")
    if(studentized){
        valid.method <- c("studentized",valid.method)
    }
    if(is.null(method)){
        method <- valid.method[1]
    }else{
        method <- match.arg(method, valid.method)
    }
    param <- names(object$estimate)
    n.param <- length(param)
    if(type %in% c("perm-res","perm-var")){
        if(!is.null(null) && "null" %in% names(match.call())){
            message("Argument \'null\' is disregarded when performing a permutation test. \n")
        }
    }else if(length(null)==1){
        null <- stats::setNames(rep(null,n.param), param)
    }else if(length(null)!=n.param){
        stop("Incorrect argument \'null\': should have length 1 or length the number of parameters to be tested (here ",n.param,"). \n")
    }else if(is.null(names(null))){
        stop("Incorrect argument \'null\': should be named with the parameters names when it has length strictly greater than 1. \n")
    }else if(any(param %in% names(null) == FALSE)){
        stop("Incorrect argument \'null\': incorrect names. \n")
    }

    valid.columns <- c("estimate","sample.estimate","se","sample.se","lower","upper","null","p.value")
    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(!is.null(columns)){
        columns <- tolower(columns)
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse = "\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse = "\" \""),"\"\n")
        }
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            columns <- union(intersect(options$columns.confint,valid.columns), unname(columns))
        }
        if(!is.null(names(columns)) && all(names(columns)=="remove")){
            columns <- setdiff(intersect(options$columns.confint,valid.columns), unname(columns))
        }
    }else{
        columns <- intersect(options$columns.confint,valid.columns)
    }

    ## ** extract information from object 
    sample.estimate <- object$sample.estimate
    sample.se <- object$sample.se
    
    ## ** point estimate
    out <- data.frame(estimate = unname(object$estimate),
                      sample.estimate = NA,
                      se = NA,
                      sample.se = NA,
                      lower = NA,
                      upper = NA,
                      null = null,
                      p.value = NA)
    rownames(out) <- names(object$estimate)
    if(method == "studentized"){
        out$se <- unname(object$se)
    }else if(method == "gaussian" || (studentized & type == "boot")){
        out$se <- apply(sample.estimate, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)
    }

    ## ** evaluate CI and p-value
    out$sample.estimate <- apply(sample.estimate, MARGIN = 2, FUN = mean, na.rm = TRUE)
    out$sample.se <- apply(sample.estimate, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)

    if(type %in% c("perm-var","perm-res") ){

        if(method == "percentile"){
            
            out$p.value <- sapply(param, function(iParam){
                iTest <- abs(sample.estimate[,iParam]) > abs(out[iParam,"estimate"])
                iP <- (correction + sum(iTest, na.rm = TRUE)) / (correction + sum(!is.na(iTest), na.rm = TRUE))
                return(iP)
            })

        }else if(method == "gaussian"){
         
            out$p.value <- 2*(1 - pnorm(abs(out$estimate-out$null)/out$sample.se))
            
        }else if(method == "studentized"){
            
            statistic  <- (out$estimate-null)/out$se
            sample.statistic <- sweep(sample.estimate, MARGIN = 1, FUN = "-", STATS = null)/sample.se
            outTable[index.var,"p.value"] <- sapply(param, function(iParam){
                iTest <- abs(sample.statistic[,iParam]) > abs(statistic[iParam])
                iP <- (correction + sum(iTest, na.rm = TRUE)) / (correction + sum(!is.na(iTest), na.rm = TRUE))
                return(iP)
            })

        }

    }else if(type == "boot"){

        if(method == "percentile"){

            out$lower <- apply(sample.estimate, MARGIN = 2, FUN = stats::quantile, probs = alpha/2, na.rm = TRUE)
            out$upper <- apply(sample.estimate, MARGIN = 2, FUN = stats::quantile, probs = 1-alpha/2, na.rm = TRUE)
            out$p.value <- sapply(param, function(iName){
                boot2pvalue(stats::na.omit(sample.estimate[,iName]),
                            null = null[iName],
                            estimate = out[iName,"estimate"],
                            alternative = "two.sided",
                            add.1 = correction)
            })

        }else if(method == "gaussian"){

            out$lower <- out$estimate + qnorm(alpha/2) * out$se
            out$upper <- out$estimate + qnorm(1-alpha/2) * out$se
            out$p.value <- 2*(1 - pnorm(abs(out$estimate-out$null)/out$se))
            
        }else if(method == "studentized"){

            sample.Wald0 <- sweep(sample.estimate, MARGIN = 2, FUN = "-", STATS = out$estimate)/sample.se  ## center around the null
            out$lower <- out$estimate + out$se * apply(sample.Wald0, MARGIN = 2, FUN = stats::quantile, probs = alpha/2, na.rm = TRUE)
            out$upper <- out$estimate + out$se * apply(sample.Wald0, MARGIN = 2, FUN = stats::quantile, probs = 1-alpha/2, na.rm = TRUE)
            out$p.value <- sapply(param, function(iName){
                boot2pvalue(stats::na.omit(out[iName,"estimate"] + out[iName,"se"] * sample.Wald0[,iName]),
                            null = null[iName],
                            estimate = out[iName,"estimate"],
                            alternative = "two.sided",
                            add.1 = correction)
            })
            
        }

    }

    ## ** export
    out <- out[,columns,drop=FALSE]
    attr(out,"method") <- method
    return(out)
}

## * confint.Wald_lmm (documentation)
##' @title Confidence Intervals for Multivariate Wald Tests
##' @description Compute confidence intervals for linear hypothesis tests, possibly with adjustment for multiple comparisons.
##' 
##' @param object a \code{Wald_lmm} object
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param method [character] type of adjustment for multiple comparisons, one of \code{"none"}, \code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}.
##' Alternatively, a method for combining the estimates, one of \code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the Wald statistic. Otherwise a normal distribution is used.
##' @param backtransform [logical] should the estimates, standard errors, and confidence intervals be backtransformed?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \bold{Adjustment for multiple comparisons}: available methods are:
##' \itemize{
##'  \item \code{"none"}, \code{"bonferroni"}, \code{"single-step2"}
##'  \item \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}: adjustment performed by [stats::p.adjust()], no confidence interval is computed.
##'  \item \code{"single-step"}, \code{"free"}, \code{"Westfall"}, \code{"Shaffer"}: adjustment performed by [multcomp::glht()],  for all but the first method no confidence interval is computed.
##' }
##' Note: method \code{"single-step"} adjusts for multiple comparisons using equicoordinate quantiles of the multivariate Student's t-distribution over all tests, instead of the univariate quantiles. It assumes equal degrees of freedom in the marginal and is described in section 7.1 of Dmitrienko et al. (2013) under the name single-step Dunnett procedure. The name \code{"single-step"} is borrowed from the multcomp package. In the book Bretz et al. (2010) written by the authors of the package, the procedure is refered to as max-t tests which is the terminology adopted in the LMMstar package.  \cr
##' When degrees of freedom differs between individual hypotheses, method \code{"single-step2"} is recommended. It simulates data using copula whose marginal distributions are Student's t-distribution (with possibly different degrees of freedom) and elliptical copula with parameters the estimated correlation between the test statistics (via the copula package). It then computes the frequency at which the simulated maximum exceed the observed maximum and appropriate quantile of simulated maximum for the confidence interval.
##'
##'  \bold{Pooling estimates}: available methods are:
##' \itemize{
##'  \item \code{"average"}: average estimates
##'  \item \code{"pool.fixse"}: weighted average of the estimates, with weights being the inverse of the squared standard error. The uncertainty about the weights is neglected.
##'  \item \code{"pool.se"}: weighted average of the estimates, with weights being the inverse of the squared standard error. The uncertainty about the weights is computed under independence of the variance parameters between models. 
##'  \item \code{"pool.gls"}: weighted average of the estimates, with weights being based on the variance-covariance matrix of the estimates. When this matrix is singular, its spectral decomposition is truncated when the correspodning eigenvalues are below \eqn{10^{-12}}. The uncertainty about the weights is neglected. 
##'  \item \code{"pool.gls1"}: similar to \code{"pool.gls"} but ensure that the weights are at most 1 in absolute value by shrinking toward the average.
##'  \item \code{"pool.rubin"}: average of the estimates and compute the uncertainty according to Rubin's rule (Barnard et al. 1999).
##' }
##'
##' @references Barnard and Rubin, \bold{Small-sample degrees of freedom with multiple imputation}. \emph{Biometrika} (1999), 86(4):948-955. \cr
##' Dmitrienko, A. and D'Agostino, R., Sr (2013), \bold{Traditional multiplicity adjustment methods in clinical trials}. \emph{Statist. Med.}, 32: 5172-5218.
##' Frank Bretz, Torsten Hothorn and Peter Westfall (2010), \bold{Multiple Comparisons Using R}, \emph{CRC Press}, Boca Raton.
##' 
##' @keywords methods

## * confint.Wald_lmm (code)
##' @export
confint.Wald_lmm <- function(object, parm, level = 0.95, df = NULL, method = NULL, columns = NULL, backtransform = NULL, ...){

    ## ** normalize user input
    if(!missing(parm) && !is.null(parm)){
        stop("Argument \'parm\' is not used - only there for compatibility with the generic method. \n")
    }
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    alpha <- 1-level
    pool.method <- options$pool.method
    adj.method <- options$adj.method

    ## handle multiple types of adjustment for multiplicty
    if(!is.null(method)){
        if(length(method)>1){
            method <- match.arg(method, c(pool.method,adj.method), several.ok = TRUE)
            if(sum(method %in% adj.method)>1){
                stop("Incorrect argument \'method\' \n",
                     "confint.Wald_lmm cannot handle several methods to adjust for multiple comparisons.")
            }
            ls.confint <- lapply(method, function(iMethod){
                iOut <- confint.Wald_lmm(object = object, method = iMethod, level = level, columns = columns, ...)
                if(iMethod %in% pool.method){
                    rownames(iOut) <- iMethod
                }
                return(iOut)
            })
            out <- do.call(rbind, ls.confint)
            attr(out,"level") <- level
            attr(out,"method") <- method
            attr(out,"contrast") <- stats::setNames(lapply(ls.confint,attr,"contrast"), method)
            attr(out,"error") <- stats::setNames(lapply(ls.confint,attr,"error"), method)
            return(out)
        }else{
            name.method <- names(method)
            method <- match.arg(method, c(pool.method,adj.method))
        }
    }else{
        name.method <- NULL
    }
    if(identical(method,"pool.se") && !inherits(object,"rbindWald_lmm")){
        method <- "pool.fixse"
        message("Argument \'method\' has been changed from \"pool.se\" to \"pool.fixse\". \n",
                "Consider using the estimate() function to account for the uncertainty of the weights. \n")
    }
    valid.columns <- names(object$univariate)
    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(!is.null(columns)){
        columns <- tolower(columns)
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse = "\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse = "\" \""),"\"\n")
        }
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            columns <- union(options$columns.confint, unname(columns))
        }
        if(!is.null(names(columns)) && all(names(columns)=="remove")){
            columns <- setdiff(options$columns.confint, unname(columns))
        }
    }else{
        columns <- options$columns.confint
    } 

    type <- unique(object$multivariate$type)
    if(is.null(df)){
        df2 <- object$args$df
    }else{
        df2 <- df
    }
    if(object$args$ci==FALSE){
        level <- NA
        if(method %in% adj.method){
            method <- "none"
        }
    }

    transform.sigma <- object$args$transform.sigma
    transform.k <- object$args$transform.k
    transform.rho <- object$args$transform.rho

    if(is.character(backtransform)){
        backtransform <-  eval(parse(text=backtransform))
    }else if(is.numeric(backtransform)){
        backtransform <- as.logical(backtransform)
    }

    test.backtransform <- stats::na.omit(c(sigma = transform.sigma, k = transform.k, rho = transform.rho))    
    if(is.null(backtransform)){            
        if(options$backtransform.confint==FALSE || length(test.backtransform[test.backtransform != "none"])==0){
            backtransform <- FALSE
        }else{
            backtransform <- object$args$backtransform
        }
    }else if(is.logical(backtransform) && length(test.backtransform[test.backtransform != "none"])==0){
        backtransform <- FALSE
    }
    n.model <- length(object$model)

    ## ** normalize df
    out <- object$univariate
    out$method <- "NA"
    if(df2){
        out$df <- pmax(out$df, options$min.df)
    }else{
        out$df <- Inf
    }
    out$statistic <- (out$estimate-out$null) / out$se

    ## ** extract info and compute CI
    n.sample <- options$n.sampleCopula
    grid <- unique(object$univariate[,c("type","test"),drop=FALSE])
    grid$type.original <- object$args$type[[1]]
    n.grid <- NROW(grid)
    attr(out,"error") <- rep(NA,n.grid)

    for(iGrid in 1:n.grid){ ## iGrid <- 1

        if(n.grid>1){
            iIndex.table <- intersect(which(out$type==grid$type[iGrid]),
                                      which(out$test==grid$test[iGrid]))
        }else{
            iIndex.table <- 1:NROW(out)
        }
        iTable <- out[iIndex.table,,drop=FALSE]
        iN.test <- NROW(iTable)

        ## *** method for multiple comparisons adjustment
        if(is.null(method)){
            if(NROW(iTable$df)==1 || all(is.na(iTable$statistic))){
                iMethod  <- "none"
            }else if(df2 == FALSE || all(is.infinite(iTable$df)) || all(abs(iTable$df - round(mean(iTable$df)))<0.1)){
                iMethod <- "single-step"
            }else{
                iMethod <- "single-step2"
            }
        }else{
            if(NROW(iTable$df)==1){
                iMethod  <- "none"
            }else{
                iMethod <- method
            }
        }

        if(iMethod %in% c("Westfall","Shaffer","free","single-step","single-step2")){
            iGlht <- object$glht[[grid[iGrid,"type.original"]]][[grid[iGrid,"test"]]]
            iGlht$df <- max(iGlht$df, min(out$df)) ## update df: cannot be smaller than the smaller df in the table
            ## may happen when df=FALSE and all df in the table have been set to Inf
        }
        out[iIndex.table,"method"] <-  iMethod
   
        ## *** evaluation
        if(iMethod == "none"){
            
            out[iIndex.table,"lower"] <- iTable$estimate + iTable$se * stats::qt(alpha/2, df = iTable$df)
            out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * stats::qt(1-alpha/2, df = iTable$df)
            out[iIndex.table,"p.value"] <- 2*(1-stats::pt(abs(iTable$statistic), df = iTable$df))
            
        }else if(iMethod %in% pool.method){

            outPool <- contrastWald_pool(object = object, index = iIndex.table, 
                                         method = method, name.method = name.method, ci = !is.na(level), df = df, alpha = alpha)
            out <- rbind(outPool,out[-iIndex.table,,drop=FALSE])

        }else if(iMethod == "single-step"){

            iCi <- confint(iGlht)
            iP <- summary(iGlht, test = multcomp::adjusted("single-step"))

            out[iIndex.table,"lower"] <- iCi$confint[,"lwr"]
            out[iIndex.table,"upper"] <- iCi$confint[,"upr"]
            out[iIndex.table,"p.value"] <- as.double(iP$test$pvalues)
            
            if(df2){
                out[iIndex.table,"df"] <- iGlht$df
            }
            attr(out, "error")[iGrid] <-  attr(iP$test$pvalues,"error")
            attr(out, "quantile")[iGrid] <-  attr(iCi$confint,"calpha")

        }else if(iMethod %in% c("Westfall","Shaffer","free")){

            iP <- summary(iGlht, test = multcomp::adjusted(iMethod))
            out[iIndex.table,"p.value"] <- as.double(iP$test$pvalues)
            
            out[iIndex.table,"lower"] <- NA
            out[iIndex.table,"upper"] <- NA
            if(df2){
                out[iIndex.table,"df"] <- iGlht$df
            }

            if(!is.null(attr(iP$test$pvalues,"error"))){
                attr(out, "error")[iGrid] <-  attr(iP$test$pvalues,"error")
            }

        }else if(iMethod == "single-step2"){

            sigma.linfct <- iGlht$linfct %*% iGlht$vcov %*% t(iGlht$linfct)
            index.n0sigma <- which(diag(sigma.linfct)>0) ## handles no variance (e.g. no treatment effect at baseline)
            rho.linfct <- stats::cov2cor(sigma.linfct[index.n0sigma,index.n0sigma,drop=FALSE])

            if(all(rho.linfct>=(1-1e-6))){ ## handles perfectly colinear case (e.g. same treatment effect at all timepoints)
                maxH0 <- abs(stats::rt(n.sample, df = mean(iTable$df[index.n0sigma])))
            }else{
                myMvd <- copula::mvdc(copula = copula::normalCopula(param=rho.linfct[lower.tri(rho.linfct)], dim = NROW(rho.linfct), dispstr = "un"),
                                      margins = rep("t", NROW(rho.linfct)),
                                      paramMargins = as.list(stats::setNames(iTable$df[index.n0sigma],rep("df",NROW(rho.linfct)))))
                maxH0 <- apply(abs(copula::rMvdc(n.sample, myMvd)), 1, max)
            }
            cH0 <- stats::quantile(maxH0, 1-alpha) 

            out[iIndex.table,"p.value"] <- sapply(abs(iTable$statistic), function(iT){(sum(iT <= maxH0)+1)/(n.sample+1)})
            out[iIndex.table,"lower"] <- iTable$estimate - iTable$se * cH0
            out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * cH0                
            attr(out, "n.sample") <-  n.sample
            attr(out, "quantile")[iGrid] <-  cH0
                
        }else if(iMethod == "bonferroni"){
                
            out[iIndex.table,"lower"] <- iTable$estimate + iTable$se * stats::qt(alpha/(2*iN.test), df = iTable$df)
            out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * stats::qt(1-alpha/(2*iN.test), df = iTable$df)
            out[iIndex.table,"p.value"] <- pmin(1,iN.test*(2*(1-stats::pt(abs(iTable$statistic), df = iTable$df))))
                
        }else{
                
            out[iIndex.table,"lower"] <- NA
            out[iIndex.table,"upper"] <- NA
            out[iIndex.table,"p.value"] <- stats::p.adjust(2*(1-stats::pt(abs(iTable$statistic), df = iTable$df)), method = iMethod)
                
        }
    }

    ## ** back-transformation
    if(is.function(backtransform) || identical(backtransform,TRUE)){

        if(is.function(backtransform)){

            out <- .backtransform(out, type.param = out$type,
                                  backtransform = TRUE, backtransform.names = NULL,
                                  transform.mu = backtransform,
                                  transform.sigma = backtransform,
                                  transform.k = backtransform,
                                  transform.rho = backtransform)

        }else{

            out <- .backtransform(out, type.param = out$type,  
                                  backtransform = TRUE, backtransform.names = object$args$backtransform.names[[1]],
                                  transform.mu = "none",
                                  transform.sigma = object$args$transform.sigma,
                                  transform.k = object$args$transform.k,
                                  transform.rho = object$args$transform.rho)

            vec.backtransform <- attr(object$univariate,"backtransform")
            if(!is.null(vec.backtransform)){
                ## case where a contrast is performed on transformed coefficients (e.g. sigma:male vs sigma:female)
                ## the back transformed version exp(log(sigma:male) - log(sigma:female)) differs from the original version sigma:male - sigma:female
                ## thus without further indication the original version is output
                out[names(vec.backtransform),"estimate"] <- unname(vec.backtransform)
                out[names(vec.backtransform),"se"] <- NA
                out[names(vec.backtransform),"df"] <- NA
                out[names(vec.backtransform),"lower"] <- NA
                out[names(vec.backtransform),"upper"] <- NA
            }

        }
    }

    ## ** export
    Umethod <- unique(out$method)
    if(length(Umethod)>1){
        Umethod <- setdiff(Umethod,"none")
    }
    out[names(out)[names(out) %in% columns == FALSE]] <- NULL
    attr(out, "level") <- 0.95
    attr(out, "method") <- Umethod
    class(out) <- append("confint_lmm", class(out))
    return(out)
}

## * helper
## ** contrastWald_pool
contrastWald_pool <- function(object, index, method, name.method, ci, df, alpha){

    ## *** extract information
    table <- object$univariate[index,,drop=FALSE]

    if(ci){
        transform.sigma <-  object$args$transform.sigma
        transform.k <-  object$args$transform.k
        transform.rho <-  object$args$transform.rho

        Sigma <- vcov(object, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)[index,index,drop=FALSE]
        independence <- attr(object$object,"independence")

        if(is.null(df)){
            if(method %in% "pool.rubin"){
                df <- TRUE
            }else if(independence & any(is.infinite(table$df)==FALSE) & method %in% c("average","pool.fixse")){
                df <- TRUE
            }else{
                df <- FALSE
            }
        }
    }

    ## *** name for the pooled estimator
    Wald.name <- rownames(table)
    sep <- object$args$sep
    pool.splitname <- strsplit(Wald.name, split = sep, fixed = TRUE)
    if(!is.null(name.method)){
        pool.name <- name.method
    }else if(all(lengths(pool.splitname)==2)){
        ## ??
        Mpool.splitname <- do.call(rbind,pool.splitname)
        pool.splitname2 <- strsplit(Mpool.splitname[,1],split="=", fixed = TRUE)
        if(all(lengths(pool.splitname2)==2)){
            Mpool.splitname2 <- do.call(rbind,pool.splitname2)
            if(length(unique(Mpool.splitname2[,1]))==1){
                pool.name <- paste0(Mpool.splitname2[1,1],"=<",paste(Mpool.splitname2[1,2],Mpool.splitname2[NROW(Mpool.splitname2),2],sep=":"),">",sep,Mpool.splitname[1,2])                   
            }else{
                pool.name <- paste0("<",paste(Mpool.splitname[1,1],Mpool.splitname[NROW(Mpool.splitname),1],sep=":"),">",sep,Mpool.splitname[1,2])                   
            }
        }else if(length(unique(Mpool.splitname[,2]))==1){
            pool.name <- paste0("<",paste(Mpool.splitname[1,1],Mpool.splitname[NROW(Mpool.splitname),1],sep=":"),">",sep,Mpool.splitname[1,2])                   
        }else{
            pool.name <- paste0("<",paste(Wald.name[1],Wald.name[length(Wald.name)],collapse=":"),">")                   
        }
    }else{
        ## <first test, last test>
        pool.name <- paste0("<",paste(Wald.name[1],Wald.name[length(Wald.name)],sep=", "),">")
    }

    ## *** Point estimate
    n.test <- NROW(table)
    if(method %in% c("average","pool.rubin")){

        pool.contrast <- matrix(1/n.test, nrow = 1, ncol = n.test)

    }else if(method %in% c("pool.fixse","pool.se")){

        pool.contrast <- matrix(1/diag(Sigma)/sum(1/diag(Sigma)), nrow = 1, ncol = n.test)

    }else if(method %in% c("pool.gls","pool.gls1")){

        poolGLS <- function(Sigma, method, tol = 1e-10){
            ## Note: without tuncation the weights can be computed as
            ## rowSums(solve(Sigma))/sum(solve(Sigma))

            ## **** truncate eigen value decomposition
            Sigma.eigen <- eigen(Sigma)
            index.subset <- which(abs(Sigma.eigen$values) > tol)

            if(length(index.subset)==0){
                stop("All eigenvalues  of the variance-covariance matrix are close to 0 (<",tol,"). \n",sep="")
            }else if(any(abs(Sigma.eigen$values) <= tol)){
                error <-  length(Sigma.eigen$values)-length(index.subset)
            }else{
                error <- NULL
            }

            lambda.k <- Sigma.eigen$values[index.subset]
            q.kj <- Sigma.eigen$vector[,index.subset,drop=FALSE]

            ## **** evaluate weigths
            qbar.k <- colSums(q.kj)
            w.k <- qbar.k^2/lambda.k
            out <- rbind(rowSums(sweep(q.kj, FUN = "*", MARGIN = 2, STATS = w.k/qbar.k))/sum(w.k))
            ## shortcut for
            ## out <- rbind(colSums(sweep(t(q.kj), FUN = "*", MARGIN = 1, STATS = w.k/qbar.k))/sum(w.k))

            if(method == "pool.gls1" && max(abs(out))>1){
                ## ensure no weight greater than 1
                iIndex.max <- which.max(abs(out))
                iMax <- abs(out)[iIndex.max]
                iMaxC <- max((1-iN.test*out)/(1-sign(out)*iN.test))
                out <- rbind(rep((1-1/iMaxC)/iN.test,iN.test) + out[1,]/iMaxC)
            }

            ## **** export
            attr(out,"error") <- error
            return(out)
        }

        pool.contrast <- poolGLS(Sigma = Sigma, method = method)        
    } 

    ## *** Variance
    if(!ci){

        pool.se  <- NA

    }else if(method %in% c("average","pool.fixse")){
        
        pool.se <- sqrt(as.double(pool.contrast %*% Sigma %*% t(pool.contrast)))
        
    }else if(independence && method == "pool.se"){
        ## \Var[\hat{pool}] = d\hat{pool}/d\Theta \Var[\Theta] t(d\hat{pool}/d\Theta)
        ##                  = [C1,C2] [A 0; 0 B] [C1;C2] = C1 A t(C1) + C2 B t(C2)
        browser()
        theta <- lapply(coef(object, effects = "all", ordering = "by", transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE),
                        function(iCoef){
                            attr(iCoef,"transform.sigma") <- transform.sigma
                            attr(iCoef,"transform.k") <- transform.k
                            attr(iCoef,"transform.rho") <- transform.rho
                            return(iCoef)
                        })
        
        ls.estimate <- lapply(names(object$model), FUN = function(iM){ ## iM <- "A"
            lava::estimate(object$model[[iM]], f = function(p){ 
                iTheta <- theta
                iTheta[[iM]] <- p
                iSigma <- vcov(object, p = iTheta, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)[index,index,drop=FALSE]
                sum(table$estimate/diag(iSigma))/sum(1/diag(iSigma))
            }, df = FALSE, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        })
        pool.se <- sqrt(sum(sapply(ls.estimate,"[[","se")^2))

    }else if((!independence && method == "pool.se") || (method %in% c("pool.gls","pool.gls1"))){

        browser()
        
    }else if(method == "pool.rubin"){
        pool.U <- mean(diag(Sigma))
        pool.B <- sum((table$estimate - mean(table$estimate))^2)/(n.test-1)
        pool.se <- sqrt(pool.U + (1 + 1/n.test) * pool.B)
    }

    ## *** Degrees of freedom
    if(!ci || df == FALSE){

        pool.df <- Inf

    }else if(method == "pool.rubin"){

        ## MICE's approach: https://stefvanbuuren.name/fimd/sec-whyandwhen.html
        pool.lambda <- (1+1/n.test)*pool.B / out$se
        pool.nu_old <- (n.test-1)/pool.lambda^2 ## formula 2.30 (Rubin 1987b eq 3.1.6)
        pool.nu_obs <- (estimate.df+1)/(estimate.df+3)*estimate.df*(1-pool.lambda) ## formula 2.31 (Barnard and Rubin (1999) )
        pool.df <- mean(pool.nu_old*pool.nu_obs/(pool.nu_old+pool.nu_obs)) ## formula 2.32

    }else if(abs(diff(range(table$df)))<0.1){

        pool.df <- mean(table$df)

    }else if(independence && method %in% c("average","pool.fixse")){

            ## \hat{pool} = \sum_k w_k \hat{beta}_k
            ## \sigma^2_{\hat{pool}} = \Var[\hat{pool}] = \sum_k w^2_k \Var[\hat{beta}_k] = \sum_k w^2_k \sigma^2_{\hat{beta}_k} by independence
            ## \Var[\sigma^2_{\hat{pool}}] = \sum_k w^4_k \Var[\sigma^2_{\hat{beta}_k}] by independence

            ## retrieve each denominator: \Var[\sigma^2_{\hat{beta}_k}]
            lmm.contrast <- coef(object, type = "ls.contrast",
                                 transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
            lmm.vcov <- vcov(object, effects = "all", df = 2,
                             transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
            lmm.df <- mapply(iC = lmm.contrast, iVcov = lmm.vcov, FUN = function(iC,iVcov){
                iDf <- .dfX(X.beta = iC, vcov.param = iVcov, dVcov.param = attr(iVcov,"dVcov"), return.vcov = TRUE)
                return(iDf)
            }, SIMPLIFY = FALSE)

            Vdenum.df <- sapply(lmm.df, attr, "denum")

            if(pool.contrast^4 %*% Vdenum.df <= 0){
                pool.df <- Inf
            }else{
                pool.df <- as.double(2*pool.se^4 / (pool.contrast^4 %*% Vdenum.df))
            }

    }else{ ## average (non-independence), pool.se, pool.gls, pool.gls1

        ## ADD-HOC APPROXIMATION (ignores correlation and uncertainty about the weights)
        lmm.contrast <- coef(object, type = "ls.contrast",
                             transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
        lmm.vcov <- vcov(object, effects = "all", df = 2,
                         transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = FALSE)
        lmm.df <- mapply(iC = lmm.contrast, iVcov = lmm.vcov, FUN = function(iC,iVcov){
            iDf <- .dfX(X.beta = iC, vcov.param = iVcov, dVcov.param = attr(iVcov,"dVcov"), return.vcov = TRUE)
            return(iDf)
        }, SIMPLIFY = FALSE)
        Vdenum.df <- sapply(lmm.df, attr, "denum")

        if(pool.contrast^4 %*% Vdenum.df <= 0){
            pool.df <- Inf
        }else{
            pool.df <- as.double(2*pool.se^4 / (pool.contrast^4 %*% Vdenum.df))
        }

    }

    ## *** export
    out <- table[1,,drop=FALSE]
    out[setdiff(names(out),c("parameter","type","test","estimate","se","df","statistic","lower","upper","null","p.value"))] <- NA
    if(length(unique(out$parameter))!=1){
        out$parameter <- NA
    }
    out$estimate <- as.double(pool.contrast %*% table$estimate)
    out$se <- pool.se
    out$df <- pool.df
    out$statistic <- as.double(pool.contrast %*% (table$estimate-table$null)/out$se)
    out$lower <- out$estimate + out$se * stats::qt(alpha/2, df = out$df)
    out$upper <- out$estimate + out$se * stats::qt(1-alpha/2, df = out$df)
    out$p.value = 2*(1-stats::pt(abs(out$statistic), df = out$df))
    rownames(out) <- pool.name
    attr(out,"contrast") <- pool.contrast
    return(out)
}


##----------------------------------------------------------------------
### confint.R ends here
