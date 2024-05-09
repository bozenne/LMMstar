### confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: May  9 2024 (12:55) 
##           By: Brice Ozenne
##     Update #: 721
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
confint.Wald_lmm <- function(object, parm, level = 0.95, method = NULL, columns = NULL, backtransform = NULL, ...){

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
    pool.method <- c("average","pool.fixse","pool.se","pool.gls","pool.gls1","pool.rubin")
    adj.method <- c(stats::p.adjust.methods,"single-step", "Westfall", "Shaffer", "free", "single-step2")

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
    df <- object$args$df
    
    if(object$args$ci==FALSE){
        return(NULL)
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
    out$df <- pmax(out$df, options$min.df)
    out$statistic <- (out$estimate-out$null) / out$se

    ## ** extract info and compute CI
    n.sample <- options$n.sampleCopula
    grid <- object$multivariate[,c("type","test")]
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
        out[iIndex.table,"method"] <-  iMethod
   
        ## *** evaluation
        if(iMethod %in% c("single-step","single-step2","Westfall","free","Shaffer")){
            iGlht <- object$glht[[grid[iGrid,"type.original"]]][[grid[iGrid,"test"]]]
        }
        iTable$p.value <- 2*(1-stats::pt( abs(iTable$statistic), df = iTable$df))

        if(iMethod == "none"){
            
            out[iIndex.table,"lower"] <- iTable$estimate + iTable$se * stats::qt(alpha/2, df = iTable$df)
            out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * stats::qt(1-alpha/2, df = iTable$df)
            out[iIndex.table,"p.value"] <- iTable$p.value
            
        }else if(iMethod %in% c("average","pool.fixse","pool.se","pool.gls","pool.gls1","pool.rubin")){

            ## replace many lines (individual effect) by a single line (for the average)
            out.save <- out
            if(utils::tail(iIndex.table,1)==NROW(out)){
                out <- out.save[1:iIndex.table[1],,drop=FALSE]
            }else{
                out <- out.save[c(1:iIndex.table[1],(iIndex.table[length(iIndex.table)]+1):NROW(out.save)),]
            }

            ## find new name
            pool.name <- rownames(out.save)[iIndex.table]
            sep <- object$args$sep
            pool.splitname <- strsplit(pool.name, split = sep, fixed = TRUE)
            if(!is.null(name.method)){
                rownames(out)[iIndex.table[1]] <- name.method
            }else if(all(lengths(pool.splitname)==2)){
                Mpool.splitname <- do.call(rbind,pool.splitname)
                pool.splitname2 <- strsplit(Mpool.splitname[,1],split="=", fixed = TRUE)
                if(all(lengths(pool.splitname2)==2)){
                    Mpool.splitname2 <- do.call(rbind,pool.splitname2)
                    if(length(unique(Mpool.splitname2[,1]))==1){
                        rownames(out)[iIndex.table[1]] <- paste0(Mpool.splitname2[1,1],"=<",paste(Mpool.splitname2[1,2],Mpool.splitname2[NROW(Mpool.splitname2),2],sep=":"),">",sep,Mpool.splitname[1,2])                   
                    }else{
                        rownames(out)[iIndex.table[1]] <- paste0("<",paste(Mpool.splitname[1,1],Mpool.splitname[NROW(Mpool.splitname),1],sep=":"),">",sep,Mpool.splitname[1,2])                   
                    }
                }else if(length(unique(Mpool.splitname[,2]))==1){
                    rownames(out)[iIndex.table[1]] <- paste0("<",paste(Mpool.splitname[1,1],Mpool.splitname[NROW(Mpool.splitname),1],sep=":"),">",sep,Mpool.splitname[1,2])                   
                }else{
                    rownames(out)[iIndex.table[1]] <- paste0("<",paste(pool.name[1],pool.name[length(pool.name)],collapse=":"),">")                   
                }
            }else{
                ## nchar.hypo <- min(nchar(pool.name))
                ## Mchar.hypo <- do.call(rbind,lapply(strsplit(pool.name, split = "", fixed = TRUE), function(iSplit){iSplit[1:nchar.hypo]}))
                ## indexchar.hypo <- which(colSums(apply(Mchar.hypo,2,duplicated)==FALSE)==1)
                ## if(length(indexchar.hypo)>0){
                ##     ## only keep common letters                    
                ##     rownames(out)[iIndex.table[1]] <- paste(Mchar.hypo[1,indexchar.hypo],collapse = "")
                ## }else{
                rownames(out)[iIndex.table[1]] <- paste0("<",paste(pool.name[1],pool.name[length(pool.name)],sep=", "),">")
                ## }    
            }

            ## find contrast
            iVcov <- object$vcov[iIndex.table,iIndex.table,drop=FALSE]
            if(iMethod %in% c("average","pool.fixse","pool.gls","pool.gls1")){

                if(iMethod=="average"){
                    iC.pool <- matrix(1/iN.test, nrow = 1, ncol = iN.test)
                }else if(iMethod %in% "pool.fixse"){
                    iIvar <- 1/out.save[iIndex.table,"se"]^2                
                    iC.pool <- matrix(iIvar/sum(iIvar), nrow = 1, ncol = iN.test)
                }else if(iMethod %in% c("pool.gls","pool.gls1")){
                    iEigen <- eigen(iVcov)
                    iEigen.subset <- which(abs(iEigen$values) > 1e-10)

                    if(length(iEigen.subset)==0){
                        stop("All eigenvalues  of the variance-covariance matrix are close to 0 (<1e-12). \n")
                    }else if(any(abs(iEigen$values) <= 1e-10)){
                        attr(out, "error")[iGrid] <-  length(iEigen$values)-length(iEigen.subset)
                    }
                    iPsum <- colSums(iEigen$vectors[,iEigen.subset,drop=FALSE])
                    iWeight <- iPsum^2/iEigen$values[iEigen.subset]
                    iWPstar <- rowSums(sweep(iEigen$vectors[,iEigen.subset,drop=FALSE], FUN = "*", MARGIN = 2, STATS = iWeight/iPsum))
                    iC.pool <- rbind(iWPstar/sum(iWeight))
                    if(iMethod %in% "pool.gls1" && max(abs(iC.pool))>1){
                        ## +/-max /(max + a) + 1/p (1-1/(max+a)) = +/-1
                        ## +/-p max + max + a - 1 = +/- p max +/- p a
                        ## max + a - 1 = +/- p a
                        ## a = (-max + 1)/(1+/-p)
                        ## +/-max + a = (1 +/- p max)/(1+/-p)
                        iIndex.max <- which.max(abs(iC.pool))
                        iMax <- abs(iC.pool)[iIndex.max]
                        iMaxC <- max((1-iN.test*iC.pool)/(1-sign(iC.pool)*iN.test))
                        iC.pool <- rbind(rep((1-1/iMaxC)/iN.test,iN.test) + iC.pool[1,]/iMaxC)
                    }                    
                }
                iVcov.pool <- as.double(iC.pool %*% iVcov %*% t(iC.pool))

                ## degree of freedom
                if(df == FALSE){
                    pool.df <- Inf
                }else if(abs(diff(range(out.save[iIndex.table,"df"])))<0.1){
                    pool.df <- mean(out.save[iIndex.table,"df"])
                }else if(is.null(attr(out.save$df,"vcov"))){
                    pool.df <- Inf
                }else{
                    iCvar.pool <- colSums(sweep(attr(out.save$df,"contrast"), FUN = "*", MARGIN = 1, STATS = iC.pool)) ## iCvar.pool[iCvar.pool!=0]
                    iVcov.vcov <- attr(out.save$df,"vcov")
                    iPair <- expand.grid(names(iCvar.pool),names(iCvar.pool))

                    idXTX.dT <- matrix(NA, nrow = 1, ncol = NROW(iPair), dimnames = list(NULL,nlme::collapse(iPair,sep=".")))
                    for(iP in 1:NROW(iPair)){ ## iP <- 35
                        idXTX.dT[,iP] <- iCvar.pool[iPair[iP,1]] * iCvar.pool[iPair[iP,2]]
                    }
                    denum <- sum((idXTX.dT %*% iVcov.vcov) * idXTX.dT)
                    if(denum==0){
                        out[denum==0] <- Inf
                    }else{
                        pool.df <- 2*iVcov.pool^2 / denum
                    }
                }
            
            }else if(iMethod %in% "pool.se"){

                iIvar <- 1/out.save[iIndex.table,"se"]^2                
                iC.pool <- matrix(iIvar/sum(iIvar), nrow = 1, ncol = iN.test)
                
                ## uncertainty relative to \beta in C\beta
                iVcov.pool <- iC.pool %*% iVcov %*% t(iC.pool)
                
                ## uncertainty relative to C in C\beta
                iLS.model <- object$model
                iLS.meancoef <- lapply(object$model, stats::coef, effects = "mean")
                iLS.vcovcoef <- lapply(object$model, stats::coef, effects = c("variance","correlation"))
                iLS.vcovvcov <- lapply(object$model, vcov, effects = c("variance","correlation"), transform.sigma = "none", transform.k = "none", transform.rho = "none")
                iLS.Coriginal <- stats::setNames(object$glht$all[[1]]$linfct.original, names(iLS.model))
                iLS.tCoriginal <- lapply(iLS.Coriginal,t)
                iVEC.beta <- cbind(out.save[iIndex.table,"estimate"]-out.save[iIndex.table,"null"])

                if(all(lengths(iLS.meancoef) == sapply(iLS.Coriginal,NCOL))){
                    effects.vcov <- "mean"
                }else{
                    effects.vcov <- "all"
                }
                                        
                icalcCB <- function(x, index){
                    iIvar.new <- iIvar
                    iIvar.new[index] <- 1/(iLS.Coriginal[[index]] %*% vcov(iLS.model[[index]], p = c(iLS.meancoef[[index]], x), effects = effects.vcov) %*% iLS.tCoriginal[[index]])
                    return((iIvar.new/sum(iIvar.new)) %*% iVEC.beta)
                }

                ## assumes independence between models
                iVcov.pool2 <- lapply(1:n.model, function(iModel){ ## iModel <- 1
                    d_calciCB <- numDeriv::jacobian(icalcCB, x = iLS.vcovcoef[[iModel]], index = iModel)
                    return(as.double(d_calciCB %*% unname(iLS.vcovvcov)[[iModel]] %*% t(d_calciCB)))
                })
                iVcov.pool <- iVcov.pool + sum(unlist(iVcov.pool2))
            
                if(length(unique(out.save[iIndex.table,"df"]))!=1){
                    pool.df <- Inf                                
                }else{
                    pool.df <- out.save[iIndex.table[1],"df"]
                }
            
            }else if(iMethod == "pool.rubin"){
                iC.pool <- matrix(1/iN.test, nrow = 1, ncol = iN.test)

                pool.U <- mean(out.save[iIndex.table,"se"]^2)
                pool.B <- sum((out.save[iIndex.table,"estimate"]-mean(out.save[iIndex.table,"estimate"]))^2)/(iN.test-1)
                iVcov.pool <- pool.U+(1+1/iN.test)*pool.B
                
                ## trying to get dfct
                if(df){
                    ## ## VERSION 1: MICE's approach ## ##
                    ## https://stefvanbuuren.name/fimd/sec-whyandwhen.html
                        
                    ## proportion of variance attributable to the missing data
                    pool.lambda <- (1+1/iN.test)*pool.B / iVcov.pool
                    pool.nu_old <- (iN.test-1)/pool.lambda^2 ## formula 2.30 (Rubin 1987b eq 3.1.6)
                    pool.nu_obs <- (out.save[iIndex.table,"df"]+1)/(out.save[iIndex.table,"df"]+3)*out.save[iIndex.table,"df"]*(1-pool.lambda) ## formula 2.31 (Barnard and Rubin (1999) )
                    pool.df <- mean(pool.nu_old*pool.nu_obs/(pool.nu_old+pool.nu_obs)) ## formula 2.32

                    ## ## VERSION 2: SATTERWAITE APPROXIMATION ## ##
                    ## see section EXTRA at the end of the file

                }else{
                    pool.df <- Inf
                }

            }
            
            out[iIndex.table[1],"estimate"]  <- iC.pool %*% out.save[iIndex.table,"estimate"]
            out[iIndex.table[1],"se"] <- sqrt(iVcov.pool)
            out[iIndex.table[1],"df"] <- pool.df
            out[iIndex.table[1],"statistic"] <- (iC.pool %*% (out.save[iIndex.table,"estimate"]-out.save[iIndex.table,"null"]))/out[iIndex.table[1],"se"]
            out[iIndex.table[1],"lower"] <- out[iIndex.table[1],"estimate"] + out[iIndex.table[1],"se"] * stats::qt(alpha/2, df = pool.df)
            out[iIndex.table[1],"upper"] <- out[iIndex.table[1],"estimate"] + out[iIndex.table[1],"se"] * stats::qt(1-alpha/2, df = pool.df)
            out[iIndex.table[1],"p.value"] <- 2*(1-stats::pt( abs(out[iIndex.table[1],"statistic"]), df = pool.df ))
            attr(out,"contrast") <- iC.pool

        }else if(iMethod == "single-step"){
            iCi <- confint(iGlht)
            iP <- summary(iGlht, test = multcomp::adjusted("single-step"))

            out[iIndex.table,"lower"] <- iCi$confint[,"lwr"]
            out[iIndex.table,"upper"] <- iCi$confint[,"upr"]
            out[iIndex.table,"p.value"] <- as.double(iP$test$pvalues)
            
            if(df){
                out[iIndex.table,"df"] <- iGlht$df
            }
            attr(out, "error")[iGrid] <-  attr(iP$test$pvalues,"error")
            attr(out, "quantile")[iGrid] <-  attr(iCi$confint,"calpha")

        }else if(iMethod %in% c("Westfall","Shaffer","free")){

            iP <- summary(iGlht, test = multcomp::adjusted(iMethod))
            out[iIndex.table,"p.value"] <- as.double(iP$test$pvalues)
            
            out[iIndex.table,"lower"] <- NA
            out[iIndex.table,"upper"] <- NA
            if(df){
                out[iIndex.table,"df"] <- iGlht$df
            }
            attr(out, "error")[iGrid] <-  attr(iP$test$pvalues,"error")

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
            cH0 <- stats::quantile(maxH0, 1-alpha)  ## attr(confint(iGlht)$confint,"calpha")

            out[iIndex.table,"p.value"] <- sapply(abs(iTable$statistic), function(iT){(sum(iT <= maxH0)+1)/(n.sample+1)})
            out[iIndex.table,"lower"] <- iTable$estimate - iTable$se * cH0
            out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * cH0                
            attr(out, "n.sample") <-  n.sample
            attr(out, "quantile")[iGrid] <-  cH0
                
            }else if(iMethod == "bonferroni"){
                
                out[iIndex.table,"lower"] <- iTable$estimate + iTable$se * stats::qt(alpha/(2*iN.test), df = iTable$df)
                out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * stats::qt(1-alpha/(2*iN.test), df = iTable$df)
                out[iIndex.table,"p.value"] <- pmin(1,iN.test*iTable$p.value)
                
            }else{
                
                out[iIndex.table,"lower"] <- NA
                out[iIndex.table,"upper"] <- NA
                out[iIndex.table,"p.value"] <- stats::p.adjust(iTable$p.value, method = iMethod)
                
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

## * confint.LRT_lmm
##' @export
confint.LRT_lmm <- function(object, parm, level = 0.95, ...){
    message("No confidence interval available for likelihood ratio tests.")
    return(NULL)
}

## * confint.effect_lmm
##' @export
confint.effect_lmm <- function(object, parm, level = 0.95, method = "none", ...){
    return(confint.Wald_lmm(object, parm, level = 0.95, method = method,  ...))
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
##' Only relevant for \code{mlmm} objects.
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

    ## ** extract confidence intervals
    out.confint <- confint.Wald_lmm(object, parm = parm, level = level, method = method, backtransform = object$args$backtransform, ...)
    if(all(method %in% c("average","pool.se","pool.gls","pool.gls1","pool.rubin") == FALSE)){
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

## * Extra
## ## VERSION 2: SATTERWAITE APPROXIMATION ## ## 
## -> weird results ## 
## note df = num/denum i.e. denum = num/df = 2*se^4/df
## and Var[mean] = mean(Var)/n = mean(2*se^4/(df*n))
## so df(mean) = 2*mean(Var)^2 / mean(2*se^4/(df*n))
##             = mean(Var)^2 / mean(se^4/(df*n))
## here the second term in the pooled variance is ignored
## pool.df <- mean(iOut.save$se^2)^2 / mean((iOut.save$se^4/(iOut.save$df)))
                        
## ## get all models
## ls.model <- estimate(object)
## ## get all coef
## ls.coef <- lapply(ls.model,coef, effects = "all")
## index.coef <- unlist(lapply(1:n.pool, function(iPool){rep(iPool, length(ls.coef[[iPool]]))}))
## vec.coef <- stats::setNames(unlist(ls.coef, use.names = FALSE),unlist(lapply(ls.coef,names)))
## ## get contrast matrix
## ls.C <- attr(iC,"model-specific")
## M.C <- as.matrix(do.call(Matrix::bdiag, ls.C))
## ## get variance-covariance of all coef
## ls.vcov <- lapply(attr(iO,"glht")[[iTest]]$model, vcov, effects = "all")
## M.vcov <- as.matrix(do.call(Matrix::bdiag, ls.vcov))
## ## derivative of the variance with respect to the model parameters
## warperPoolVar <- function(x){ ## x <- vec.coef
##     term1 <- mean(sapply(1:n.pool, function(iPool){ ## iPool <- 1
##         ls.C[[iPool]] %*% vcov(ls.model[[iPool]], p = x[index.coef==iPool], effects = "all") %*% t(ls.C[[iPool]])
##     }))
##     term2 <- var(M.C %*% x)*(1+1/n.pool)
##     return(as.double(term1+term2)) ## WARNING: only consider the first term
## }
## ## sanity check
## ## warperPoolVar(vec.coef) - pool.variance
## nabla <- numDeriv::jacobian(warperPoolVar, x = vec.coef, method = "Richardson")
## ## Satterthwaite approximation
## pool.df <- as.double(2 * warperPoolVar(vec.coef)^2/ (nabla %*% M.vcov %*% t(nabla)))

##----------------------------------------------------------------------
### confint.R ends here
