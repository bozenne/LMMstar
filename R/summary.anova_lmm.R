### summary.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:50) 
## Version: 
## Last-Updated: mar  4 2022 (15:42) 
##           By: Brice Ozenne
##     Update #: 144
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.anova_lmm (documentation)
##' @title Summary of Testing for a Linear Mixed Models
##' @description Estimates, p-values, and confidence intevals for linear hypothesis testing, possibly adjusted for multiple comparisons.
##' 
##' @param object an \code{anova_lmm} object, output of \code{anova}.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}.
##' @param transform [function] function to backtransform the estimates, standard errors, null hypothesis, and the associated confidence intervals
##' (e.g. \code{exp} if the outcomes have been log-transformed).
##' @param level [numeric 0-1] level of the confidence intervals.
##' @param print.nulls [logical] should the estimates for the individual null hypotheses be displayed?
##' @param seed [integer] value that will be set before adjustment for multiple comparisons to ensure reproducible results.
##' Can also be \code{NULL}: in such a case no seed is set.
##' @param columns [character vector] Columns to be displayed for each null hypothesis.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"df.denom"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param digits [interger] number of digits used to display estimates.
##' @param digits.p.value [interger] number of digits used to display p-values.
##' @param ... Not used. For compatibility with the generic method.
 
## * summary.anova_lmm (code)
##' @export
summary.anova_lmm <- function(object, method = "single-step", transform = NULL, level = 0.95, print.nulls = TRUE, seed = NULL, columns = NULL,
                              digits = max(3L, getOption("digits") - 2L),
                              digits.p.value = max(3L, getOption("digits") - 2L),
                              ...){

    if(!is.null(seed)){
        old.seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit( assign(".Random.seed", old.seed, envir = .GlobalEnv, inherits = FALSE) )
        set.seed(seed)
    }

    dots <- list(...)
    if(length(dots)>0){
        if(identical(names(dots),"test")){
            stop("Unknown argument \'test\'. Use argument \'method\' instead, e.g. method=\"single-step\". \n")
        }else{
            stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
        }
    }
    if(!is.null(match.call()$method) && is.null(match.call()$print.nulls)){
        print.nulls <- 2
    }

    options <- LMMstar.options()
    if(is.null(columns)){
        columns <- c(options$columns.anova, "df.denom", "df.num")
    }else{
        columns <- match.arg(columns, choices = c("null","estimate","se","statistic","df","lower","upper","p.value"), several.ok = TRUE)
        if("df" %in% columns){
            columns <- c(columns,"df.num","df.denom")
        }
    }

    out <- list()
    
    if(attr(object,"test")=="Wald"){
        type <- setdiff(names(object),"call")
        ci <- stats::confint(object, level = level, method = method, simplify = FALSE)
        for(iType in type){

            if(is.null(object[[iType]])){next}
            object.print <- object[[iType]]
            object.print <- cbind(object.print,
                                  " " = stats::symnum(object.print$p.value, corr = FALSE, na = FALSE, 
                                                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                                       symbols = c("***", "**", "*", ".", " "))
                                  )
            object.print$p.value <- as.character(signif(object.print$p.value, digits = digits.p.value))
            iNoDf <- is.infinite(object.print$df.denom)
            txt.test <- "Multivariate Wald test (global null hypothesis)"
            if(iType == "all"){
                cat("\n\t", "|| User-specified linear hypotheses || \n", sep="")
                ## print(object$call)
            }else{
                cat("\n\t","|| ",iType," coefficients || \n", sep="")
            }
            if(print.nulls>-1){
                if(iType == "all"){
                    if(is.na(object.print$statistic) && !is.null(attr(object.print$statistic,"error"))){
                        cat("\n - ",txt.test,": ",attr(object.print$statistic,"error"),"\n", sep="")
                    }else{
                        cat("\n - ",txt.test,"\n", sep="")
                        print(object.print[,names(object.print) %in% c(columns,"statistic"," "),drop=FALSE], digits = digits, row.names = FALSE)
                        out <- c(out, list(object[[iType]]))
                    }
                }else{
                    cat("\n - ",txt.test,"\n",sep="")
                    print(object.print[,names(object.print) %in% c(columns,"statistic"," "),drop=FALSE], digits = digits)
                    out <- c(out, list(object[[iType]]))
                }
                if(print.nulls==FALSE && "p.value" %in% columns && iType == utils::tail(type,1)){
                    cat("---\n",
                        "Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n")
                }
            }
            if(!is.null(ci[[iType]]) && abs(print.nulls)>=1){
                if(!is.null(transform)){
                    ci[[iType]] <- lapply(ci[[iType]], function(iCI){transformSummaryTable(iCI, transform = transform)})
                }

                cat("\n - Univariate Wald test (individual null hypotheses) \n", sep="")
                object.print <- do.call(rbind,unname(ci[[iType]]))
                stats::printCoefmat(object.print[,names(object.print) %in% columns,drop=FALSE], digits = digits,
                                    has.Pvalue = "p.value" %in% columns,
                                    P.values = "p.value" %in% columns,
                                    eps.Pvalue = 10^{-digits.p.value},
                                    signif.legend = TRUE)
                out <- c(out, list(do.call(rbind,unname(ci[[iType]]))))
                if(!is.null(transform)){
                    cat("\n Columns ",paste(intersect(columns,c("estimate","se","lower","upper")), collapse =", ")," have been back-transform. \n",sep="")
                }

                if(attr(object,"robust")){
                    cat("Standard errors: robust\n")
                }else{
                    cat("Standard errors: model-based\n")
                }
                if(all(sapply(ci[[iType]],NROW)==1) || method == "none"){ ## always only one hypothesis in each global test
                    cat("(CIs/p-values not adjusted for multiple comparisons) \n", sep="")
                }else if(length(ci[[iType]])==1){ ## only one global test
                    if(method=="bonferroni"){
                        cat("(CIs/p-values adjusted for multiple comparisons -- Bonferroni)\n", sep="")
                    }else if(method == "single-step"){
                        cat("(CIs/p-values adjusted for multiple comparisons -- single step max-test)\n", sep="")
                    }else{
                        cat(paste0("(CIs/p-values adjusted for multiple comparisons -- ",method,")\n", sep=""),sep="")
                    }
                }else{
                    if(method=="bonferroni"){
                        cat("(CIs/p-values adjusted for multiple comparisons within each global test -- bonferroni) \n", sep="")
                    }else if(method == "single-step"){
                        cat("(CIs/p-values adjusted for multiple comparisons within each global test -- single step max-test) \n", sep="")
                    }else{
                        cat(paste0("(CIs/p-values adjusted for multiple comparisons within each global test -- ",method,") \n", sep=""),sep="")
                    }
                }

                if(method == "single-step" && "p.value" %in% columns){
                    error <- max(c(0,unlist(lapply(ci[[iType]],function(iO){attr(iO$p.value,"error")}))))
                    if(error > 1e-12){
                        txt.error <- paste0("Error when computing the adjusted p-value by numerical integration: ", signif(error, digits = 5))
                        if(!is.null(seed)){
                            txt.error <- paste0(txt.error," (seed ",seed,")")
                        }
                        cat(txt.error,"\n")
                    }
                }


            }
        }
        cat("\n")
    }else if(attr(object,"test")=="LRT"){
        cat(" - Likelihood ratio test \n")
        out <- as.data.frame(object)
        out.print <- out
        if("null" %in% columns == FALSE){
            out.print[["null"]] <- NULL
        }
        print(out.print)
    }

    return(invisible(out))
}


##----------------------------------------------------------------------
### summary.anova_lmm.R ends here
