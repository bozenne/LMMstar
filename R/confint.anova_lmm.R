### confint.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: jun 27 2022 (11:53) 
##           By: Brice Ozenne
##     Update #: 127
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
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}, \code{"partial.R"}.
##' @param simplify [logical] Return a data.frame instead of a list containing a data.frame when possible.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Method \code{"single-step"} adjust for multiple comparisons using quantiles of the multivariate Student's t-distribution, assuming equal degrees of freedom in the marginal.
##' This is performed by the multcomp package.
##'
##' When degrees of freedom differs between individual hypotheses, method \code{"single-step2"} is recommended. It simulates data using copula whose marginal distributions are Student's t-distribution (with possibly different degrees of freedom) and elliptical copula with parameters the estimated correlation between the test statistics. This is performed by the copula package.
##' 
##' @export
confint.anova_lmm <- function(object, parm, level = 0.95, method = NULL, columns = NULL, simplify = TRUE, ...){

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

    if(!is.null(columns)){
        columns  <- match.arg(columns, c("estimate","se","statistic","df","lower","upper","null","p.value","partial.R"), several.ok = TRUE)
    }else{
        columns <- options$columns.confint
    }
    ## ** extract info and compute CI
    out <- lapply(object[setdiff(names(object),"call")], function(iO){ ## iO <- object[[1]]

        iTable <- attr(iO,"CI")
        if(is.null(iTable) || all(sapply(iTable,is.null))){return(NULL)}
        iOut <- stats::setNames(vector(mode = "list", length = length(iTable)),names(iTable))

        for(iTest in 1:length(iTable)){ ## iTest <- 1
            iOut[[iTest]] <- iTable[[iTest]]
            iOut[[iTest]]$df <- pmax(iOut[[iTest]]$df, options$min.df)

            if(is.null(method)){
                if(length(unique(round(iOut[[iTest]]$df)))>1){
                    iMethod <- "single-step2"
                }else{
                    iMethod <- "single-step"
                }
            }else{
                iMethod <- method
            }
            attr(iOut[[iTest]], "method") <-  iMethod

            if(iMethod == "none" || NROW(iOut[[iTest]])==1){
                iOut[[iTest]]$lower <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(alpha/2, df = iOut[[iTest]]$df)
                iOut[[iTest]]$upper <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(1-alpha/2, df = iOut[[iTest]]$df)
                iOut[[iTest]]$p.value <- 2*(1-stats::pt( abs((iOut[[iTest]]$estimate-iOut[[iTest]]$null) / iOut[[iTest]]$se), df = iOut[[iTest]]$df))
            }else if(iMethod == "pool"){
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
                iOut[[iTest]]$statistic <- iOut[[iTest]]$estimate/iOut[[iTest]]$se
                iOut[[iTest]]$lower <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(alpha/2, df = iOut[[iTest]]$df)
                iOut[[iTest]]$upper <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(1-alpha/2, df = iOut[[iTest]]$df)
                iOut[[iTest]]$p.value <- 2*(1-stats::pt( abs((iOut[[iTest]]$estimate-iOut[[iTest]]$null) / iOut[[iTest]]$se), df = iOut[[iTest]]$df))
            }else if(iMethod == "single-step"){

                iGlht <- attr(iO,"glht")[[iTest]]
                iCi <- confint(iGlht)
                iOut[[iTest]]$lower <- iCi$confint[,"lwr"]
                iOut[[iTest]]$upper <- iCi$confint[,"upr"]
                iOut[[iTest]]$p.value <- summary(iGlht, test = multcomp::adjusted("single-step"))$test$pvalues
                iOut[[iTest]]$df <- iGlht$df

            }else if(iMethod == "single-step2"){
                iGlht <- attr(iO,"glht")[[iTest]]
                n.sample <- options$n.sampleCopula

                rho <- stats::cov2cor(iGlht$linfct %*% iGlht$vcov %*% t(iGlht$linfct))
                n.marginal <- length(iOut[[iTest]]$df)

                myMvd <- copula::mvdc(copula = copula::normalCopula(param=rho[lower.tri(rho)], dim = NROW(rho), dispstr = "un"),
                                      margins = rep("t", n.marginal),
                                      paramMargins = as.list(stats::setNames((iOut[[iTest]]$df),rep("df",n.marginal))))
                maxH0 <- sort(apply(abs(copula::rMvdc(n.sample, myMvd)), 1, max))
                
                iOut[[iTest]]$p.value <- sapply(abs(iOut[[iTest]]$statistic), function(iT){(sum(iT <= maxH0)+1)/(n.sample+1)})

                cH0 <- stats::quantile(maxH0, 0.95) ## attr(confint(iGlht)$confint,"calpha")
                iOut[[iTest]]$lower <- iOut[[iTest]]$estimate - iOut[[iTest]]$se * cH0
                iOut[[iTest]]$upper <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * cH0
                attr(iOut[[iTest]], "n.sample") <-  n.sample
            }else if(iMethod == "bonferroni"){
                p <- NROW(iOut[[iTest]])
                iOut[[iTest]]$lower <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(alpha/(2*p), df = iOut[[iTest]]$df)
                iOut[[iTest]]$upper <- iOut[[iTest]]$estimate + iOut[[iTest]]$se * stats::qt(1-alpha/(2*p), df = iOut[[iTest]]$df)
                iOut[[iTest]]$p.value <- pmin(1,2*p*(1-stats::pt( abs((iOut[[iTest]]$estimate-iOut[[iTest]]$null) / iOut[[iTest]]$se), df = iOut[[iTest]]$df)))
            }else{
                iOut[[iTest]]$lower <- NA
                iOut[[iTest]]$upper <- NA
                iOut[[iTest]]$p.value <- stats::p.adjust(2*(1-stats::pt( abs((iOut[[iTest]]$estimate-iOut[[iTest]]$null) / iOut[[iTest]]$se), df = iOut[[iTest]]$df)), method = iMethod)
            }
            iOut[[iTest]][names(iOut[[iTest]])[names(iOut[[iTest]]) %in% columns == FALSE]] <- NULL
        }        
        return(iOut)
    })
    if(simplify && length(out)==1){
        out <- out[[1]]
        if(simplify && length(out)==1){
            out <- out[[1]]
        }
    }
    return(out)
}


##----------------------------------------------------------------------
### confint.anova_lmm.R ends here
