### confint.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:51) 
## Version: 
## Last-Updated: mar 14 2022 (13:27) 
##           By: Brice Ozenne
##     Update #: 57
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
##' @param simplify [logical] Return a data.frame instead of a list containing a data.frame when possible.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Method \code{"single-step"} adjust for multiple comparisons using quantiles of the multivariate Student's t-distribution, assuming equal degrees of freedom in the marginal.
##' This is performed by the multcomp package.
##'
##' When degrees of freedom differs between individual hypotheses, method \code{"single-step2"} is recommended. It simulates data using copula whose marginal distributions are Student's t-distribution (with possibly different degrees of freedom) and elliptical copula with parameters the estimated correlation between the test statistics. This is performed by the copula package.
##' 
##' @export
confint.anova_lmm <- function(object, parm, level = 0.95, method = NULL, simplify = TRUE, ...){

    ## ** normalize user input
    if(attr(object,"test") == "LRT"){
        message("No confidence interval available for likelihood ratio tests.")
        return(NULL)
    }
    if(!missing(parm)){
        stop("Argument \'parm\' is not used - only there for compatibility with the generic method. \n")
    }
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    alpha <- 1-level
    if(!is.null(method)){
        method <- match.arg(method, c(stats::p.adjust.methods,"single-step", "single-step2"))
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
            }else if(iMethod == "single-step"){
                iGlht <- attr(iO,"glht")[[iTest]]
                iCi <- confint(iGlht)
                iOut[[iTest]]$lower <- iCi$confint[,"lwr"]
                iOut[[iTest]]$upper <- iCi$confint[,"upr"]
                iOut[[iTest]]$p.value <- summary(iGlht, test = multcomp::adjusted("single-step"))$test$pvalues
                iOut[[iTest]]$df <- iGlht$df
            }else if(iMethod == "single-step2"){
                iGlht <- attr(iO,"glht")[[iTest]]
                requireNamespace("copula")
                n.sample <- options$n.sampleCopula

                rho <- stats::cov2cor(iGlht$linfct %*% iGlht$vcov %*% t(iGlht$linfct))
                n.marginal <- length(iOut[[iTest]]$df)

                myMvd <- copula::mvdc(copula = copula::normalCopula(param=rho[lower.tri(rho)], dim = NROW(rho), dispstr = "un"),
                                      margins = rep("t", n.marginal),
                                      paramMargins = as.list(stats::setNames(iOut[[iTest]]$df,rep("df",n.marginal))))
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
