### summary.Manova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 24 2022 (11:09) 
## Version: 
## Last-Updated: jan 24 2022 (11:11) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.Manova_lmm
#' @title Summary of Testing Across Linear Mixed Models
#' @description Estimates, p-values, and confidence intevals for linear hypothesis testing, possibly adjusted for multiple comparisons.
#' 
#' @param object a \code{Manova_lmm} object, output of \code{glht}.
#' @param confint [logical] should confidence intervals be output
#' @param conf.level [numeric 0-1] level of the confidence intervals.
#' @param transform [function] function to backtransform the estimates, standard errors, null hypothesis, and the associated confidence intervals
#' (e.g. \code{exp} if the outcomes have been log-transformed).
#' @param seed [integer] value that will be set before adjustment for multiple comparisons to ensure reproducible results.
#' Can also be \code{NULL}: in such a case no seed is set.
#' @param rowname.rhs [logical] when naming the hypotheses, add the right-hand side (i.e. "X1-X2=0" instead of "X1-X2").
#' @param ... argument passed to \code{multcomp:::summary.glht}, e.g. argument \code{test} to choose the type of adjustment for multiple comparisons.
#' @export
summary.Manova_lmm <- function(object, confint = TRUE, conf.level = 0.95, transform = NULL, seed = NULL, rowname.rhs = TRUE, ...){

    if(!is.null(seed)){
        old.seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit( assign(".Random.seed", old.seed, envir = .GlobalEnv, inherits = FALSE) )
        set.seed(seed)
    }

    keep.class <- class(object)
    object$test <- NULL
    object$confint <- NULL
    class(object) <- setdiff(keep.class, "Manova_lmm")
    keep.df <- object$df
    test.df <- any( (keep.df>0) * (!is.infinite(keep.df)) == 1 )
    object$df <- round(stats::median(object$df))
    output <- summary(object, ...)
    ## restaure df when possible
    method.adjust <- output$test$type
    if(NROW(object$linfct)==1){method.adjust <- "none"}
    if(test.df && method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none","univariate")){
        output$df <- keep.df
        output$test$pvalues <- stats::p.adjust(2*(1-stats::pt(abs(output$test$tstat), df = keep.df)), method = method.adjust)
    }
    
    name.hypo <- rownames(output$linfct)
    n.hypo <- length(name.hypo)
    if(confint && method.adjust %in% c("univariate","none","bonferroni","single-step")){
        if(method.adjust %in% c("none","univariate","bonferroni")){
            alpha <- switch(method.adjust,
                            "none" = 1-conf.level,
                            "univariate" = 1-conf.level,
                            "bonferroni" = (1-conf.level)/n.hypo)
            if(test.df){
                q <- stats::qt(1-alpha/2, df = output$df)
            }else{
                q <- stats::qnorm(1-alpha/2)
            }
            output$confint <- data.frame(matrix(NA, ncol = 3, nrow = n.hypo,
                                                dimnames = list(name.hypo, c("Estimate","lwr","upr"))))
            output$confint$Estimate <- as.double(output$test$coef)
            output$confint$lwr <- as.double(output$test$coef - q * output$test$sigma)
            output$confint$upr <- as.double(output$test$coef + q * output$test$sigma)
            ## range(confint(output, level = 1-alpha, calpha = univariate_calpha())$confint-output$confint)
        }else if(method.adjust == "single-step"){
            output <- confint(output, level = conf.level, calpha = multcomp::adjusted_calpha())
        }else{
            output$confint <- matrix(NA, nrow = n.hypo, ncol = 3,
                                     dimnames = list(name.hypo, c("Estimate","lwr","upr")))
        }
    }
    if(rowname.rhs){
        table2.rownames <- paste0(name.hypo, " == ", output$rhs)
    }else{
        table2.rownames <- name.hypo
    }
    output$table2 <- data.frame(matrix(NA, nrow = n.hypo, ncol = 7,
                                       dimnames = list(table2.rownames,
                                                       c("estimate","se","df","lower","upper","statistic","p.value"))
                                       ), stringsAsFactors = FALSE)
    output$table2$estimate <- output$test$coefficients
    output$table2$se <- output$test$sigma
    output$table2$df <- output$df
    output$table2$df[output$table2$df==0] <- Inf
    output$table2$lower <- output$confint[,"lwr"]
    output$table2$upper <- output$confint[,"upr"]
    output$table2$statistic <- output$test$tstat
    output$table2$p.value <- output$test$pvalues
    output$seed <- seed
    
    ## ** transformation
    output$transform <- transform
    output$table2 <- transformSummaryTable(output$table2,
                                           transform = transform)

    ## ** export    
    class(output) <- append(c("summary.Manova_lmm","summary.glht"),keep.class)
    return(output)
}

## * print.summary.Manova_lmm
#' @export
print.summary.Manova_lmm <- function(x,
                                     digits = max(3L, getOption("digits") - 2L),
                                     digits.p.value = max(3L, getOption("digits") - 2L),
                                     columns = c("estimate","se","df","lower","upper","statistic","p.value"),
                                     ...){
    
    columns <- match.arg(columns, choices = c("estimate","se","df","lower","upper","statistic","p.value"), several.ok = TRUE)
    type <- x$type
    call <- if(isS4(x$model)){x$model@call}else{x$model$call}
    alternative <- x$alternativ
    type <- x$test$type
    txt.type <- switch(type,
                       "univariate" = "(CIs/p-values not adjusted for multiple comparisons)", 
                       "none" = "(CIs/p-values not adjusted for multiple comparisons)", 
                       "single-step" = paste0("(CIs/p-values adjusted for multiple comparisons -- single step max-test)"), 
                       "free" = paste0("(CIs/p-values adjusted for multiple comparisons -- step down max-test)"), 
                       "Westfall" = paste0("(CIs/p-values adjusted for multiple comparisons -- step down max-test with logical restrictions)"), 
                       paste0("(CIs/p-values adjusted for multiple comparisons -- ", type, " method)")
                       )
    txt.robust <- switch(as.character(x$robust),
                         "TRUE" = "Robust",
                         "FALSE" = "Model-based"
                         )

    ## txt.correction <- switch(as.character(x$ssc),
    ##                          "Cox" = " corrected for small sample bias (Cox correction)",
    ##                          "residuals" = " corrected for small sample bias (residual correction)",
    ##                          "NA" = ""
    ##                          )
    
    txt.alternative <- switch(alternative,
                              "less" = "one sided tests - inferiority",
                              "greater" = "one sided tests - superiority",
                              "two.sided" = "two sided tests")

    ## display
    cat("\n\t", "Simultaneous Tests for General Linear Hypotheses\n\n")
    if (!is.null(type)) {
        cat("Multiple Comparisons of Means (",txt.alternative,") \n\n", sep = "")
    }
    cat("Standard errors: ",txt.robust,"\n",sep="")
    cat("Linear Hypotheses:\n")
    stats::printCoefmat(x$table2[,columns[columns %in% names(x$table2)],drop=FALSE], digits = digits,
                        has.Pvalue = "p.value" %in% columns,
                        P.values = "p.value" %in% columns,
                        eps.Pvalue = 10^{-digits.p.value})

    if(NROW(x$table2)>1){
        cat(txt.type,"\n")
    }
    error <- attr(x$test$pvalues,"error")
    if(!is.null(error) && error > 1e-12 && "p.value" %in% columns){
        txt.error <- paste0("Error when computing the adjusted p-value by numerical integration: ", signif(error, digits = digits))
        if(!is.null(x$seed)){
            txt.error <- paste0(txt.error," (seed ",x$seed,")")
        }
        cat(txt.error,"\n")
    }
    
    if(!is.null(x$global)){
        cat("\nGlobal test: p.value=",format.pval(x$global["p.value"], digits = digits, eps = 10^(-digits.p.value)),
            " (statistic=",round(x$global["statistic"], digits = digits),
            ", df=",round(x$global["df"], digits = digits),")\n",sep="")
    }
    ## if(nchar(txt.correction)>0){cat("(",txt.correction,")\n",sep="")}
    cat("\n")
    return(invisible(x))
}



##----------------------------------------------------------------------
### summary.Manova_lmm.R ends here
