### confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: sep 29 2025 (13:15) 
##           By: Brice Ozenne
##     Update #: 1361
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * confint.lmm (documentation)
##' @title Confidence Intervals for a Linear Mixed Model
##' @description Compute confidence intervals (CIs) relative to parameters from a linear mixed model. 
##' 
##' @param object a \code{lmm} object.
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param effects [character] Should the CIs/p-values for all coefficients be output (\code{"all"}),
##' or only for mean coefficients (\code{"mean"} or \code{"fixed"}),
##' or only for variance coefficients (\code{"variance"}),
##' or only for correlation coefficients (\code{"correlation"}).
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors.
##' Can also be \code{2} compute the degrees-of-freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
##' @param null [numeric vector] the value of the null hypothesis relative to each coefficient.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficients. Otherwise a normal distribution is used.
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
##' \item column df: the degrees-of-freedom.
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
    table.param <- stats::model.tables(object, effects = "param")
    name.param <- table.param$name
    type.param <- stats::setNames(table.param$type, name.param)

    ## ** normalize user imput
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

    ## *** parm
    if(!is.null(parm)){
        stop("Argument \'parm\' should not be used. It is here for compatibility with the generic method. \n",
             "Use \'effects\' instead. \n")
    }

    ## *** effects
    if(is.null(effects)){
        if((is.null(transform.sigma) || identical(transform.sigma,"none")) && (is.null(transform.k) || identical(transform.k,"none")) && (is.null(transform.rho) || identical(transform.rho,"none"))){
            effects <- options$effects
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else{
        if(!is.character(effects) || !is.vector(effects)){
            stop("Argument \'effects\' must be a character vector. \n")
        }
        valid.effects <- c("mean","fixed","variance","correlation","all")
        if(any(effects %in% valid.effects == FALSE)){
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
        }
        if(all("all" %in% effects)){
            if(length(effects)>1){
                stop("Argument \'effects\' must have length 1 when containing the element \"all\". \n")
            }else{
                effects <- c("mean","variance","correlation")
            }
        }else{
            effects[effects== "fixed"] <- "mean"
        }
    }

    ## *** df
    if(is.null(df)){
        df <- !is.null(object$df)
    }

    ## *** backtransform
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

    ## *** transform
    ## used to decide on the null hypothesis of k parameters
    init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = object$reparametrize$transform.sigma, x.transform.k = object$reparametrize$transform.k, x.transform.rho = object$reparametrize$transform.rho)
    
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    transform <- init$transform

    ## *** type.information
    if(is.null(type.information)){
        type.information <- object$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## *** columns
    valid.columns <- c("estimate","se","statistic","df","lower","upper","null","p.value")
    save.columns <- columns
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
    beta <- stats::coef(object, effects = effects, 
                        transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names, options = options)
    p <- length(beta)
    nameNoTransform.beta <- stats::model.tables(object, effects = c("param",effects))$name
    name.beta <- names(beta)
    type.beta <- type.param[name.beta]

    ## ** get uncertainty
    vcov.beta <- stats::vcov(object, effects = effects, df = df, robust = robust,
                             type.information = type.information, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names,
                             options = options)

    if(df){
        df <- pmax(attr(vcov.beta,"df"), options$min.df)
        attr(vcov.beta,"df") <- NULL
        if((object$args$method.fit=="REML") && (type.information != "observed") && ("mean" %in% effects)){
            warning("when using REML with expected information, the degrees-of-freedom of the mean parameters may depend on the parametrisation of the variance parameters. \n")
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
                      df = df[name.beta], lower = as.numeric(NA), upper = as.numeric(NA), statistic = as.numeric(NA), null = null, p.value = as.numeric(NA),
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
            backtransform.names <- stats::model.tables(object, effects = c("param",effects), 
                                                       transform.sigma = gsub("log","",transform.sigma), transform.k = gsub("log","",transform.k), transform.rho = gsub("atanh","",transform.rho),
                                                       transform.names = transform.names)$trans.name

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
    if(!identical(save.columns,"all")){
        if(is.null(columns) || !is.null(names(columns))){
            out[names(out)[names(out) %in% columns == FALSE]] <- NULL
        }else{
            attr.out <- attributes(out)[setdiff(names(attributes(out)),c("names","row.names","class"))]
            out <- out[,columns,drop=FALSE]
            attributes(out) <- c(attributes(out),attr.out)
        }
    }
    class(out) <- append("confint_lmm", class(out))
    return(out)
}

## * confint.lmmCC (code)
##' @export
confint.lmmCC <- function(object, parm = NULL, level = 0.95, effects = NULL, columns = NULL, ...){

    if(object$time$n==4 && (is.null(effects) || effects == "change")){

        Mcon <- cbind(c(-1,1,0,0),c(0,0,-1,1))
        out.estimate <- lava::estimate(object, function(p){
            Sigma.change <- t(Mcon) %*% stats::sigma(object, p = p) %*% Mcon
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
        out <- stats::confint(object, parm = parm, effects = effects, level = level, columns = columns, ...)
    }

    ## ** export
    return(out)

}



## * confint.Wald_lmm (documentation)
##' @title Confidence Intervals From Wald Tests
##' @description Compute pointwise or simultaneous confidence
##'     intervals of linear contrasts involved in Wald tests, along
##'     with corresponding p-values. Pointwise confidence intervals
##'     have nominal coverage w.r.t. a single contrast whereas
##'     simultaneous confidence intervals have nominal coverage
##'     w.r.t. to all contrasts.  Can also be used to pool estimates.
##' 
##' @param object a \code{Wald_lmm} object
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the Wald statistic. Otherwise a normal distribution is used.
##' @param method [character] Should pointwise confidence intervals be output (\code{"none"}) or simultaneous confidence intervals (\code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}) and/or confidence intervals for pooled linear contrast estimates (\code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}, \code{"p.rejection"})?
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param backtransform [logical] should the estimates, standard errors, and confidence intervals be backtransformed? Ignored when pooling estimates.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \bold{Multiple testing methods}:
##' \itemize{
##'  \item \code{"none"}: no adjustment
##'  \item \code{"bonferroni"}: bonferroni adjustment
##'  \item \code{"single-step"}: max-test adjustment performed by \code{multcomp::glht()} finding equicoordinate quantiles of the multivariate Student's t-distribution over all tests, instead of the univariate quantiles. It assumes equal degrees-of-freedom in the marginal and is described in section 7.1 of Dmitrienko et al. (2013) under the name single-step Dunnett procedure. The name \code{"single-step"} is borrowed from the multcomp package. In the book Bretz et al. (2010) written by the authors of the package, the procedure is refered to as max-t tests which is the terminology adopted in the LMMstar package.
##'  \item \code{"single-step2"}: max-test adjustment performed using Monte Carlo integration. It simulates data using copula whose marginal distributions are Student's t-distribution (with possibly different degrees-of-freedom) and elliptical copula with parameters the estimated correlation between the test statistics (via the copula package). It then computes the frequency at which the simulated maximum exceed the observed maximum and appropriate quantile of simulated maximum for the confidence interval.
##'  \item \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}: other adjustments performed by \code{stats::p.adjust()}. No confidence interval is computed.
##'  \item \code{"free"}, \code{"Westfall"}, \code{"Shaffer"}: other adjustments performed by \code{multcomp::glht()}. No confidence interval is computed.
##' }
##' When degrees-of-freedom differs between individual hypotheses, \code{"single-step2"} is recommended over \code{"single-step"}. \cr \cr
##'
##' \bold{Pooling methods}:
##' \itemize{
##'  \item \code{"average"}: average estimates
##'  \item \code{"pool.se"}: weighted average of the estimates, with weights being the inverse of the squared standard error. 
##'  \item \code{"pool.gls"}: weighted average of the estimates, with weights being based on the variance-covariance matrix of the estimates. When this matrix is singular, the Moore–Penrose inverse is used which correspond to truncate the spectral decomposition for eigenvalues below \eqn{10^{-12}}. 
##'  \item \code{"pool.gls1"}: similar to \code{"pool.gls"} with weights shrinked toward the average whenever they exceed 1 in absolute value.
##'  \item \code{"pool.rubin"}: average of the estimates and compute the uncertainty according to Rubin's rule (Barnard et al. 1999). Validity requires the congeniality condition of Meng (1994).
##'  \item \code{"p.rejection"}: proportion of null hypotheses where there is evidence for an effect. By default the critical quantile (defining the level of evidence required) is evaluated using a \code{"single-step"} method but this can be changed by adding adjustment method in the argument \code{method}, e.g. \code{effects=c("bonferronin","p.rejection")}.
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
    pool.method <- options$pool.method
    adj.method <- options$adj.method
    n.sample <- options$n.sampleCopula

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling anova. \n")
        return(invisible(NULL))
    }

    ## *** parm
    if(!missing(parm) && !is.null(parm)){
        stop("Argument \'parm\' is not used - only there for compatibility with the generic method. \n")
    }

    ## *** level
    if(object$args$univariate==FALSE){
        level <- NA
        method <- "none"
    }
    alpha <- 1-level

    ## *** df
    if(is.null(df)){
        df <- object$args$df
    }

    ## *** method
    if(!is.null(method)){
        if(!is.character(method) || !is.vector(method)){
            stop("Argument \'method\' must be a character.")
        }
        if(NROW(object$univariate)==1 && any(method!="none")){
            message("Argument \'method\' is ignored when there is a single test. \n")
            method <- "none"
        }
        if(any(method %in% c(adj.method,pool.method) == FALSE)){
            stop("Unknown value for argument \'method\': \"",paste(setdiff(method,c(adj.method,pool.method)),collapse = "\", \""),"\". \n",
                 "Possible values: \"",paste(c(adj.method,pool.method), collapse = "\", \""),"\". \n")
        }
        if(sum(method %in% adj.method)>1){
            stop("Incorrect argument \'method\' \n",
                 "Cannot handle several approaches to adjust for multiple comparisons.")
        }
        ## add method used to adjust the rejection threshold
        if("p.rejection" %in% method){
            if(!is.null(attr(method,"qt"))){
                qt <- attr(method,"qt")
                if(length(qt)>1){
                    stop("Attribute \"qt\" to argument \'method\' should have length 1. \n")
                }
                if(is.na(qt)){
                    stop("Attribute \"qt\" to argument \'method\' should no be NA. \n")
                }
                if(!is.numeric(qt) && !is.character(qt)){
                    stop("Attribute \"qt\" to argument \'method\' should be a numeric value or a character. \n")
                }
                if(is.numeric(qt) && qt<0){
                    stop("When a numeric, attribute \"qt\" to argument \'method\' should be non-negatie. \n")
                }
                if(is.character(qt) && (qt %in% c("none","bonferroni","single-step","single-step2")==FALSE)){
                    stop("When a character, attribute \"qt\" to argument \'method\' should be one of \"none\", \"bonferroni\", \"single-step\", \"single-step2\". \n")
                }
            }
            if(is.null(attr(method,"qt")) & length(intersect(method, adj.method))==1){
                qt <- intersect(method, adj.method)
            }else{
                message("Use \"single-step\" method to define statistical significance with method=\"p.rejection\". \n",
                        "Consider adding an attribute \"qt\" to the \'method\' argument to avoid this message. \n")
                qt <- "single-step"
            }
        }
    }else{
        qt <- NULL
    }

    ## *** columns
    valid.columns <- names(object$univariate)
    save.columns <- columns
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

    ## simplify method
    if(!is.null(method) && "p.value" %in% columns == FALSE && any(method %in% c("Westfall","Shaffer","free"))){
        method[method %in% adj.method] <- "none"
    }
    if(!is.null(method) && (all(c("lower","upper","quantile","p.value") %in% columns == FALSE) && any(method %in% adj.method))){
        method[method %in% adj.method] <- "none"
    }

    ## *** backtransform
    transform.sigma <- object$args$transform.sigma
    transform.k <- object$args$transform.k
    transform.rho <- object$args$transform.rho

    if(is.null(backtransform)){
        backtransform <- any(object$univariate$tobacktransform)
    }else if(is.character(backtransform)){
        backtransform <-  eval(parse(text=backtransform))
    }else if(is.numeric(backtransform)){
        backtransform <- as.logical(backtransform)
    }

    ## ** normalize df
    out <- object$univariate
    out$method <- "NA"
    if(df){
        out$df <- pmax(out$df, options$min.df)
    }else{
        out$df <- Inf
    }

    ## ** prepare pooling
    if("null" %in% columns || "p.value" %in% columns){
        compute.null <- NULL ## compute the null
    }else{
        compute.null <- NA ## do not compute the null
    }
    compute.ci <- any(c("se","lower","upper","p.value") %in%  columns)
    
    ## ** extract info and compute CI
    grid <- unique(object$univariate$name)
    n.grid <- length(grid)
    attr(out,"error") <- rep(NA,n.grid)

    for(iGrid in 1:n.grid){ ## iGrid <- 1

        iIndex.table <- which(out$name == grid[iGrid])
        iTable <- out[iIndex.table,,drop=FALSE]
        iN.test <- NROW(iTable)

        ## *** method 
        if(is.null(method)){
            if(NROW(iTable$df)==1 || all(is.na(iTable$statistic))){
                iMethod  <- "none"
            }else if(df == FALSE || all(is.infinite(iTable$df)) || all(abs(iTable$df - round(mean(iTable$df)))<0.1)){
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
        iMethod.adj <- intersect(iMethod, adj.method)

        ## *** adjustment for multiple comparisons
        if(length(iMethod.adj)>0 && iMethod.adj %in% c("Westfall","Shaffer","free","single-step","single-step2")){
            iGlht <- object$glht[[grid[iGrid]]]
            iGlht$df <- max(iGlht$df, min(out$df)) ## update df: cannot be smaller than the smaller df in the table
            ## may happen when df=FALSE and all df in the table have been set to Inf
        }
        if(length(iMethod.adj)>0){
            out[iIndex.table,"method"] <-  iMethod.adj
        }

        if(length(iMethod.adj)>0 && iMethod.adj == "single-step"){

            if(df){
                iGlht$df <- round(iGlht$df)
                out[iIndex.table,"df"] <- iGlht$df
            }

            if("p.value" %in% columns){
                iP <- summary(iGlht, test = multcomp::adjusted("single-step"))
                out[iIndex.table,"p.value"] <- as.double(iP$test$pvalues)
                attr(out, "error")[iGrid] <-  attr(iP$test$pvalues,"error")
            }
            
            if(any(c("lower","upper","quantile") %in% columns) || is.character(qt)){
                iCi <- stats::confint(iGlht)
                out[iIndex.table,"lower"] <- iCi$confint[,"lwr"]
                out[iIndex.table,"upper"] <- iCi$confint[,"upr"]
                out[iIndex.table,"quantile"] <-  attr(iCi$confint,"calpha")
            }

        }else if(length(iMethod.adj)>0 && iMethod.adj %in% c("Westfall","Shaffer","free")){

            iP <- summary(iGlht, test = multcomp::adjusted(iMethod.adj))
            out[iIndex.table,"p.value"] <- as.double(iP$test$pvalues)
            
            if(df){
                out[iIndex.table,"df"] <- iGlht$df
            }

            if(!is.null(attr(iP$test$pvalues,"error"))){
                attr(out, "error")[iGrid] <-  attr(iP$test$pvalues,"error")
            }

        }else if(length(iMethod.adj)>0 && iMethod.adj == "single-step2"){

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

            if("p.value" %in% columns){
                out[iIndex.table,"p.value"] <- sapply(abs(iTable$statistic), function(iT){(sum(iT <= maxH0)+1)/(n.sample+1)})
            }

            if(any(c("lower","upper","quantile") %in% columns) || is.character(qt)){
                cH0 <- stats::quantile(maxH0, 1-alpha) 
                out[iIndex.table,"lower"] <- iTable$estimate - iTable$se * cH0
                out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * cH0                
                out[iIndex.table,"quantile"] <- cH0
            }

            attr(out, "n.sample") <-  n.sample
            
        }else if(length(iMethod.adj)>0 && iMethod.adj %in% c("none",stats::p.adjust.methods)){
            
            if("p.value" %in% columns){
                out[iIndex.table,"p.value"] <- 2*(1-stats::pt(abs(iTable$statistic), df = iTable$df))
                if(iMethod.adj %in% stats::p.adjust.methods){
                    out[iIndex.table,"p.value"] <- stats::p.adjust(out[iIndex.table,"p.value"], method = iMethod.adj)
                }
            }

            if(any(c("lower","upper","quantile") %in% columns) || is.character(qt)){
                if(iMethod.adj == "none"){
                    cH0 <- stats::qt(p = 1-alpha/2, df = iTable$df)
                }else if(iMethod.adj == "bonferroni"){
                    cH0 <- stats::qt(p = 1-alpha/(2*iN.test), df = iTable$df)
                }else{
                    cH0 <- as.numeric(NA)
                }
                out[iIndex.table,"lower"] <- iTable$estimate - iTable$se * cH0
                out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * cH0
                out[iIndex.table,"quantile"] <- cH0
            }
        }

        ## *** pooling
        if(any(iMethod %in% pool.method)){
            if(is.character(qt) && (qt %in% method)){
                qt <- out$quantile
            }
            ## do not always request null and p.value as it can be numerically expensive for method = "p.rejection"
            iPool <- aggregate(object, method = intersect(method, pool.method), columns = union(columns,c("type","transformed","tobacktransform","method")), qt = qt,
                               level = level, df = df, options = options)            
            out <- rbind(out,matrix(NA, nrow = NROW(iPool), ncol = NCOL(out),
                                    dimnames = list(rownames(iPool),names(out))))
            out[seq(from = NROW(out)-NROW(iPool)+1,to = NROW(out), by = 1),union(columns,c("type","transformed","tobacktransform","method"))] <- iPool
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

            ## force no back-transform, e.g. when comparing two correlation coefficients
            out <- .backtransform(out,
                                  type.param = ifelse(out$tobacktransform,out$type,"mu"),  
                                  backtransform = TRUE, backtransform.names = NULL,
                                  transform.mu = "none",
                                  transform.sigma = object$args$transform.sigma,
                                  transform.k = object$args$transform.k,
                                  transform.rho = object$args$transform.rho)

            if(object$args$transform.rho=="atanh" && any(out$type=="rho" & out$n.param > 1 & out$tobacktransform==FALSE)){
                ## case where a contrast is performed on transformed coefficients (e.g. rho:male vs rho:female)
                ## the back transformed version tanh(atanh(rho:male) - atanh(rho:female)) differs from the original version rho:male - rho:female
                ## thus without further indication the original version is output
                out[out$type=="rho" & out$n.param > 1 & out$tobacktransform==FALSE,"se"] <- NA
                out[out$type=="rho" & out$n.param > 1 & out$tobacktransform==FALSE,"lower"] <- NA
                out[out$type=="rho" & out$n.param > 1 & out$tobacktransform==FALSE,"upper"] <- NA
            }

        }
    }

    ## ** export
    Umethod <- unique(out$method)
    if(length(Umethod)>1){
        Umethod <- setdiff(Umethod,"none")
    }
    if(!is.null(method) && all(method %in% pool.method)){
        out <- out[which(out$method %in% method),,drop=FALSE]
    }
    if(!identical(save.columns,"all")){
        if(is.null(columns) || !is.null(names(columns))){
            out[names(out)[names(out) %in% columns == FALSE]] <- NULL
        }else{
            attr.out <- attributes(out)[setdiff(names(attributes(out)),c("names","row.names","class"))]
            out <- out[,columns,drop=FALSE]
            attributes(out) <- c(attributes(out),attr.out)
        }
    }
    attr(out, "level") <- 0.95
    attr(out, "method") <- Umethod
    class(out) <- append("confint_lmm", class(out))
    return(out)
}

## * confint.LRT_lmm
##' @export
confint.LRT_lmm <- function(object, parm, level = 0.95, ...){
    message("No confidence interval available for likelihood ratio tests.")
    return(NULL)
}
## * confint.rbindWald_lmm
##' @title Confidence Intervals From Combined Wald Tests
##' @description Combine pointwise or simultaneous confidence intervals relative to linear contrasts of parameters from different linear mixed models. 
##'  Can also output p-values (corresponding to pointwise confidence intervals) or adjusted p-values (corresponding to simultaneous confidence intervals).
##'
##' @param object a \code{rbindWald_lmm} object.
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the Wald statistic. Otherwise a normal distribution is used.
##' @param method [character] Should pointwise confidence intervals be output (\code{"none"}) or simultaneous confidence intervals (\code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}) and/or confidence intervals for pooled linear contrast estimates (\code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.gls1"}, \code{"pool.rubin"}, \code{"p.rejection"})?
##' Only relevant when \code{effects = "Wald"}.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' @param backtransform [logical] should the estimates be back-transformed?
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @export
confint.rbindWald_lmm <- function(object, parm, level = 0.95, df = NULL, method = NULL, columns = NULL, ordering = NULL, backtransform = NULL, ...){

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
    pool.method <- options$pool.method
    adj.method <- options$adj.method

    ## *** object
    if(object$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm \n")
        return(invisible(NULL))
    }

    ## *** ordering
    if(!is.null(ordering)){
        ordering <- match.arg(ordering, c("contrast","model"))
        if(any(method %in% pool.method)){
            message("Argument \'ordering\' is ignored when argument \'method\' is \"",intersect(pool.method,method)[1],"\". \n")
            ordering <- NULL
        }
    }

    ## ** confint
    out <- confint.Wald_lmm(object, level = level, df = df, method = method, columns = columns, backtransform = backtransform, options = options)

    ## ** re-order
    if(!is.null(ordering)){
        if(ordering=="model"){
            reorder <- order(object$univariate$model)
        }else if(is.list(object$univariate$term)){
            reorder <- order(object$univariate$type,sapply(object$univariate$term, paste, collapse = ";"))
        }else{
            reorder <- order(object$univariate$type,object$univariate$term)
        }
        out <- out[reorder,,drop=FALSE]
    }

    ## ** export
    return(out)
}


## * confint.mlmm (documentation)
##' @title Confidence Intervals for Multiple Linear Mixed Models
##' @description Compute pointwise or simultaneous confidence intervals relative to parameter or linear contrasts of parameters from group-specific linear mixed models.
##' Can also output p-values (corresponding to pointwise confidence intervals) or adjusted p-values (corresponding to simultaneous confidence intervals).
##' 
##' @param object an \code{mlmm} object, output of \code{mlmm}.
##' @inheritParams confint.rbindWald_lmm
##'
##' @export
confint.mlmm <- confint.rbindWald_lmm

## * confint.effect_lmm
##' @export
confint.effect_lmm <- function(object, parm, level = 0.95, method = "none", ...){
    return(confint.Wald_lmm(object, parm, level = 0.95, method = method,  ...))
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
confint.resample <-  function(object, parm = NULL, level = 0.95, null = NULL, method = NULL, columns = NULL, correction = TRUE, ...){

    mycall <- match.call()
    type <- object$args$type
    studentized <- object$args$studentized
    param <- names(object$estimate)
    n.param <- length(param)

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

    ## *** parm
    if(!missing(parm) && !is.null(parm)){
        stop("Argument \'parm\' is not used - only there for compatibility with the generic method. \n")
    }

    ## *** level
    alpha <- 1-level

    ## *** method
    valid.method <- c("percentile","gaussian")
    if(studentized){
        valid.method <- c("studentized",valid.method)
    }
    if(is.null(method)){
        method <- valid.method[1]
    }else{
        method <- match.arg(method, valid.method)
    }

    ## *** null
    if(type %in% c("perm-res","perm-var")){
        if(!is.null(null) && "null" %in% names(mycall)){
            message("Argument \'null\' is disregarded when performing a permutation test. \n")
        }
    }else if(is.null(null)){
        ## default based on propagating the null from each parameter
        null <- object$null 
    }else if(length(null)==1){
        null <- stats::setNames(rep(null,n.param), param)
    }else if(length(null)!=n.param){
        stop("Incorrect argument \'null\': should have length 1 or length the number of parameters to be tested (here ",n.param,"). \n")
    }else if(is.null(names(null))){
        stop("Incorrect argument \'null\': should be named with the parameters names when it has length strictly greater than 1. \n")
    }else if(any(param %in% names(null) == FALSE)){
        stop("Incorrect argument \'null\': incorrect names. \n")
    }

    ## *** columns
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
         
            out$p.value <- 2*(1 - stats::pnorm(abs(out$estimate-out$null)/out$sample.se))
            
        }else if(method == "studentized"){
            
            statistic  <- (out$estimate-null)/out$se
            sample.statistic <- sweep(sample.estimate, MARGIN = 1, FUN = "-", STATS = null)/sample.se
            outTable[,"p.value"] <- sapply(param, function(iParam){
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

            out$lower <- out$estimate + stats::qnorm(alpha/2) * out$se
            out$upper <- out$estimate + stats::qnorm(1-alpha/2) * out$se
            out$p.value <- 2*(1 - stats::pnorm(abs(out$estimate-out$null)/out$se))
            
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



##----------------------------------------------------------------------
### confint.R ends here
