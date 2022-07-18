### confint.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: jul 18 2022 (16:39) 
##           By: Brice Ozenne
##     Update #: 240
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.anova_lmm
##' @title Confidence Intervals for Multivariate Wald Tests
##' @description Compute confidence intervals for linear hypothesis tests, possibly with adjustment for multiple comparisons.
##' 
##' @param object a \code{anova_lmm} object
##' @param parm Not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}, \code{"partial.r"}.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Method \code{"single-step"} adjust for multiple comparisons using quantiles of the multivariate Student's t-distribution, assuming equal degrees of freedom in the marginal.
##' This is performed by the multcomp package.
##'
##' When degrees of freedom differs between individual hypotheses, method \code{"single-step2"} is recommended. It simulates data using copula whose marginal distributions are Student's t-distribution (with possibly different degrees of freedom) and elliptical copula with parameters the estimated correlation between the test statistics. This is performed by the copula package.
##' 
##' @export
confint.anova_lmm <- function(object, parm, level = 0.95, method = NULL, columns = NULL, ...){

    ## ** normalize user input
    if(attr(object,"test") == "LRT"){
        message("No confidence interval available for likelihood ratio tests.")
        return(NULL)
    }
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
    backtransform <- object$args$backtransform

    if(is.function(backtransform) || all(is.character(backtransform)) || identical(backtransform,TRUE) || identical(backtransform,1)){

        if(is.function(backtransform) || all(is.character(backtransform))){

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

##----------------------------------------------------------------------
### confint.anova_lmm.R ends here
