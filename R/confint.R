### confint.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: jul 20 2022 (16:16) 
##           By: Brice Ozenne
##     Update #: 262
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
##' @name confint
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
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}, \code{"partial.R"}.
##' @param type.information,transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param backtransform [logical] should the variance/covariance/correlation coefficient be backtransformed?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @seealso the function \code{anova} to perform inference about linear combinations of coefficients and adjust for multiple comparisons.
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
    valid.columns <- c("estimate","se","statistic","df","lower","upper","null","p.value","partial.r")
    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(!is.null(columns)){
        columns <- tolower(columns)
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse = "\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse = "\" \""),"\"\n")
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
        if((object$method.fit=="REML") && (type.information != "observed") && ("mean" %in% effects)){
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
                      partial.R = as.numeric(NA),
                      stringsAsFactors = FALSE)

    out$statistic <- (out$estimate-null)/out$se
    out$p.value <- 2*(1-stats::pt(abs(out$statistic), df = out$df))
    index.cor <- setdiff(which(type.beta=="mu"), which(name.beta=="(Intercept)"))
    if(length(index.cor)>0){
        out[index.cor,"partial.R"] <- sign(out$statistic[index.cor])*sqrt(out$statistic[index.cor]^2/(out$df[index.cor]+out$statistic[index.cor]^2))
        ## from "An R2 statistic for fixed effects in the linear mixed model" by Lloyd J. Edwards et al. 2008 (Statistic in medicine)
        ## Equation 19
        ## DOI: 10.1002/sim.3429
    }
    
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

## * confint.Wald_lmm (documentation)
##' @title Confidence Intervals for Multivariate Wald Tests
##' @description Compute confidence intervals for linear hypothesis tests, possibly with adjustment for multiple comparisons.
##' 
##' @param object a \code{anova_lmm} object
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}, \code{"partial.r"}.
##' @param backtransform [logical] should the estimates, standard errors, and confidence intervals be backtransformed?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Method \code{"single-step"} adjust for multiple comparisons using quantiles of the multivariate Student's t-distribution, assuming equal degrees of freedom in the marginal.
##' This is performed by the multcomp package.
##'
##' When degrees of freedom differs between individual hypotheses, method \code{"single-step2"} is recommended. It simulates data using copula whose marginal distributions are Student's t-distribution (with possibly different degrees of freedom) and elliptical copula with parameters the estimated correlation between the test statistics. This is performed by the copula package.
##' 

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
    if(!is.null(method)){
        method <- match.arg(method, c(stats::p.adjust.methods,"single-step", "single-step2","pool"))
    }

    valid.columns <- c("type","test","method","estimate","se","statistic","df","lower","upper","null","p.value","partial.r")
    if(identical(columns,"all")){
        columns <- valid.columns
    }else if(!is.null(columns)){
        columns <- tolower(columns)
        if(any(columns %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns[columns %in% valid.columns == FALSE], collapse = "\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns), collapse = "\" \""),"\"\n")
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
    if((is.null(backtransform) || is.logical(backtransform)) && (is.na(transform.sigma) && is.na(transform.k) && is.na(transform.rho))){
        backtransform <- FALSE
    }else if(is.null(backtransform)){
        if(all(object$univariate$type %in% c("mu","all"))){
            backtransform <- FALSE ## no need for back-transformation
        }else{
            backtransform <- options$backtransform.confint
        }
    }

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
            if(NROW(iTable$df)==1){
                iMethod  <- "none"
            }else if(length(unique(round(iTable$df)))>1){
                iMethod <- "single-step2"
            }else{
                iMethod <- "single-step"
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
        if(iMethod %in% c("single-step","single-step2")){
            iGlht <- object$glht[[grid[iGrid,"type.original"]]][[grid[iGrid,"test"]]]
        }
        iTable$p.value <- 2*(1-stats::pt( abs(iTable$statistic), df = iTable$df))

        if(iMethod == "none"){
            
            out[iIndex.table,"lower"] <- iTable$estimate + iTable$se * stats::qt(alpha/2, df = iTable$df)
            out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * stats::qt(1-alpha/2, df = iTable$df)
            out[iIndex.table,"p.value"] <- iTable$p.value
            
        }else if(iMethod == "pool"){
            browser()
            iOut.save <- iOut[[iTest]]
            pool.name <- rownames(iOut.save)
            n.pool <- NROW(iOut.save)
            df <- attr(object,"df")
            iC <- attr(iO,"glht")[[iTest]]$linfct
                
                pool.estimate  <- mean(iOut.save$estimate-iOut.save$null)
                
                if(inherits(object,"mlmm")){
                    sep <- attr(object,"sep")
                    pool.name <- sapply(strsplit(pool.name,split=sep,fixed =TRUE), function(iVec){paste(iVec[-1],collapse="")})
                    pool.U <- mean(iOut.save$se^2)
                    pool.B <- sum((iOut.save$estimate-pool.estimate)^2)/(n.pool-1)
                    pool.variance <- pool.U+(1+1/n.pool)*pool.B
                
                    ## trying to get dfct
                    if(df){
                        ## ## VERSION 1: MICE's approach ## ##
                        ## https://stefvanbuuren.name/fimd/sec-whyandwhen.html
                        
                        ## proportion of variance attributable to the missing data
                        pool.lambda <- (1+1/n.pool)*pool.B / pool.variance
                        pool.nu_old <- (n.pool-1)/pool.lambda^2 ## formula 2.30 (Rubin 1987b eq 3.1.6)
                        pool.nu_obs <- (iOut.save$df+1)/(iOut.save$df+3)*iOut.save$df*(1-pool.lambda) ## formula 2.31 (Barnard and Rubin (1999) )
                        pool.df <- mean(pool.nu_old*pool.nu_obs/(pool.nu_old+pool.nu_obs)) ## formula 2.32

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
                    }else{
                        pool.df <- Inf
                    }
                }else{
                    iC.pool <- rep(1/n.pool, n.pool) %*% iC
                    pool.variance <- iC.pool %*% attr(iO,"glht")[[iTest]]$vcov %*% t(iC.pool)
                    if(df){
                        pool.df <- .dfX(X.beta = iC.pool, vcov.param = attr(iO,"glht")[[iTest]]$model$vcov, dVcov.param = attr(iO,"glht")[[iTest]]$model$dVcov)
                    }else{
                        pool.df <- Inf
                    }
                }
                iOut[[iTest]] <- data.frame(estimate = pool.estimate,
                                            se = sqrt(pool.variance),
                                            df = pool.df,
                                            statistic = NA,
                                            lower = NA,
                                            upper = NA,
                                            null = 0,
                                            p.value = NA)

            ## try to find common pattern
            nchar.hypo <- min(nchar(pool.name))
            Mchar.hypo <- do.call(rbind,lapply(strsplit(pool.name, split = "", fixed = TRUE), function(iSplit){iSplit[1:nchar.hypo]}))
            indexchar.hypo <- which(colSums(apply(Mchar.hypo,2,duplicated)==FALSE)==1)
            if(length(indexchar.hypo)>0){
                rownames(iOut[[iTest]]) <- paste(Mchar.hypo[1,indexchar.hypo],collapse = "")
            }
            
            out[iIndex.table,"lower"] <- iTable$estimate + iTable$se * stats::qt(alpha/2, df = iTable$df)
            out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * stats::qt(1-alpha/2, df = iTable$df)
            out[iIndex.table,"p.value"] <- 2*(1-stats::pt( abs(iTable$statistic), df = iTable$df))
            
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

            }else if(iMethod == "single-step2"){

                rho <- stats::cov2cor(iGlht$linfct %*% iGlht$vcov %*% t(iGlht$linfct))

                myMvd <- copula::mvdc(copula = copula::normalCopula(param=rho[lower.tri(rho)], dim = NROW(rho), dispstr = "un"),
                                      margins = rep("t", iN.test),
                                      paramMargins = as.list(stats::setNames(iTable$df,rep("df",iN.test))))
                maxH0 <- apply(abs(copula::rMvdc(n.sample, myMvd)), 1, max)
                cH0 <- stats::quantile(maxH0, 1-alpha)  ## attr(confint(iGlht)$confint,"calpha")
                
                out[iIndex.table,"p.value"] <- sapply(abs(iTable$statistic), function(iT){(sum(iT <= maxH0)+1)/(n.sample+1)})
                out[iIndex.table,"lower"] <- iTable$estimate - iTable$se * cH0
                out[iIndex.table,"upper"] <- iTable$estimate + iTable$se * cH0                
                attr(out, "n.sample") <-  n.sample
                
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

## * confint.mlmm (documentation)
##' @title Confidence Intervals for Multiple Linear Mixed Model.
##' @description Compute confidence intervals for several linear mixed models.
##' 
##' @param object an \code{mlmm} object, output of \code{mlmm}.
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}, or \code{"pool"}.
##' @param ... other arguments are passed to \code{\link{confint.anova_lmm}}.
##'
##' @details Statistical inference following pooling is performed according to Rubin's rule whose validity requires the congeniality condition of Meng (1994).
##'
##' @references
##' Meng X. L.(1994). Multiple-imputation inferences with uncongenial sources of input. Statist. Sci.9, 538–58.

## * confint.mlmm (code)
##' @export
confint.mlmm <- function(object, parm = NULL, level = 0.95, method = NULL, ...){

    if(is.null(method)){
        method <- "none"
    }

    return(confint.anova_lmm(object, parm = parm, level = level, method = method, ...))
    
}

##----------------------------------------------------------------------
### confint.R ends here
