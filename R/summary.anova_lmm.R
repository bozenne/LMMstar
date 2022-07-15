### summary.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  9 2022 (14:50) 
## Version: 
## Last-Updated: jul 15 2022 (18:28) 
##           By: Brice Ozenne
##     Update #: 284
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
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}.
##' @param transform [function] function to backtransform the estimates, standard errors, null hypothesis, and the associated confidence intervals
##' (e.g. \code{exp} if the outcomes have been log-transformed).
##' @param level [numeric 0-1] level of the confidence intervals.
##' @param print [logical] should the output be printed in the console.
##' Can be a vector of length 2 where the first element refer to the global tests and the second to the individual tests.
##' @param seed [integer] value that will be set before adjustment for multiple comparisons to ensure reproducible results.
##' Can also be \code{NULL}: in such a case no seed is set.
##' @param columns [character vector] Columns to be displayed for each null hypothesis.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param digits [interger] number of digits used to display estimates.
##' @param digits.p.value [interger] number of digits used to display p-values.
##' @param ... Not used. For compatibility with the generic method.
##'
##'
##' @details By default adjustment for multiple comparisons via a single step max-test adjustment,
##' either using the multcomp package (equal degrees of freedom) or the copula package (unequal degrees of freedom).
 
## * summary.anova_lmm (code)
##' @export
summary.anova_lmm <- function(object, method = NULL, transform = NULL, level = 0.95, print = TRUE, seed = NULL, columns = NULL,
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
    if(!is.null(match.call()$method) && is.null(match.call()$print)){
        print <- TRUE
    }
    if(length(print)==1){
        print.indiv <- print
        print.global <- print
    }else if(length(print)>2){
        stop("Argument \'print\' should have length at most 2. \n",
             "The first element refering to global test and the second to individual hypotheses. \n")
    }else{
        print.global <- print[1]
        print.indiv <- print[2]
    }

    options <- LMMstar.options()
    valid.columns <- c("null","estimate","se","statistic","df","lower","upper","p.value","partial.r","")
    if(identical(columns,"all")){
        columns.global <- setdiff(valid.columns, c("estimate", "se", "lower", "upper"))
        columns.indiv <- valid.columns
    }else  if(is.null(columns)){
        columns.indiv <- options$columns.anova
        columns.global <- union("statistic", setdiff(options$columns.anova, c("estimate", "se", "lower", "upper")))
    }else{
        columns.indiv <- tolower(columns)
        if(any(columns.indiv %in% valid.columns == FALSE)){
            stop("Incorrect value(s) \"",paste(columns.indiv[columns.indiv %in% valid.columns == FALSE], collapse ="\" \""),"\" for argument \'columns\'. \n",
                 "Valid values: \"",paste(setdiff(valid.columns, columns.indiv), collapse ="\" \""),"\".\n")
        }
        columns.global <- setdiff(columns.indiv, c("estimate", "se", "lower", "upper"))
    }
    if("df" %in% columns.global){
        index.df <- which(columns.global == "df")
        if(index.df == 1){
            columns.global <- c("df.num", "df.denom", columns.global[(index.df+1):length(columns.global)])
        }else if(index.df == length(columns.global)){
            columns.global <- c(columns.global[1:(index.df-1)], "df.num", "df.denom")
        }else{
            columns.global <- c(columns.global[1:(index.df-1)], "df.num", "df.denom", columns.global[(index.df+1):length(columns.global)])
        }
    }
    columns.global <- gsub("^partial.r$","partial.r2", columns.global)
    object.df <- object$args$df
    object.robust <- object$args$robust
    object.ci <- object$args$ci
    
    if(attr(object,"test")=="Wald"){

        table.multivariate <- object$multivariate
        table.multivariate$type.original <- object$args$type[[1]]
        type <- unique(table.multivariate$type.original)

        if(object.ci){
            table.univariate <- confint(object, level = level, method = method, columns = union(c("type","test","method"),columns.indiv))

            typetest2type.original <- stats::setNames(table.multivariate$type.original,paste(table.multivariate$type,table.multivariate$test,sep="|"))
            table.univariate$type.original <- typetest2type.original[paste(table.univariate$type,table.univariate$test,sep="|")]
            univariate.method <- attr(table.univariate,"method")
        }

        for(iType in type){ ## iType <- type[1]

            ## ** type of test
            if(iType == "all" && (print.global>0.5 || print.indiv>0.5)){
                cat("\n\t", "|| User-specified linear hypotheses || \n", sep="")
                ## print(object$call)
            }else if(print.global>0.5 || print.indiv>0.5){
                cat("\n\t","|| ",iType," coefficients || \n", sep="")
            }

            ## ** Global tests
            object.print <- table.multivariate[table.multivariate$type.original==iType,,drop=FALSE]
            object.print <- cbind(object.print,
                                  " " = stats::symnum(object.print$p.value, corr = FALSE, na = FALSE, 
                                                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                                      symbols = c("***", "**", "*", ".", " "))
                                  )
            names(object.print)[NCOL(object.print)] <- ""
            object.print$p.value <- as.character(signif(object.print$p.value, digits = digits.p.value))
            
            if(print.global){

                txt.test <- "Multivariate Wald test (global null hypothesis)"
                if(all(is.na(object.print$statistic)) && !is.null(attr(object.print$statistic,"error"))){
                    cat("\n - ",txt.test,": ",attr(object.print$statistic,"error"),"\n", sep="")
                }else{
                    cat("\n - ",txt.test,"\n", sep="")
                    object.print <-  object.print[,names(object.print) %in% columns.global,drop=FALSE]
                    if(object.df){
                        names(object.print) <- gsub("^statistic","F-statistic",names(object.print))
                    }else{
                        names(object.print) <- gsub("^statistic","chi2-statistic",names(object.print))
                    }
                    print(object.print, digits = digits, row.names = (iType != "all"))
                }
                if(print.indiv==FALSE && "" %in% columns.global && iType == utils::tail(type,1)){
                    cat("---\n",
                        "Signif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n")
                }
            }

            ## ** individual specific tests
            if(object.ci){
                                
                if(print.indiv){

                    cat("\n - Univariate Wald test (individual null hypotheses) \n", sep="")
                    object.print <- table.univariate[table.univariate$type.original==iType,,drop=FALSE]
                    n.hypoPerTest <- table(paste(object.print$type,object.print$test,sep="|"))
                    if(length(univariate.method)>1){
                        warning("Different methods have been used to adjust for multiple comparisons - text describing the adjustment will not be accurate.")
                    }
                    object.print <- object.print[,names(object.print) %in% columns.indiv,drop=FALSE]
                    if(object.df){
                        names(object.print) <- gsub("^statistic","t-statistic",names(object.print))
                    }else{
                        names(object.print) <- gsub("^statistic","z-statistic",names(object.print))
                    }
                    stats::printCoefmat(object.print, digits = digits,
                                        has.Pvalue = "p.value" %in% columns.indiv,
                                        P.values = "p.value" %in% columns.indiv,
                                        eps.Pvalue = 10^{-digits.p.value},
                                        signif.legend = TRUE)
                    if(!is.null(transform)){
                        cat("\n Columns ",paste(intersect(columns.indiv,c("estimate","se","lower","upper")), collapse =", ")," have been back-transform. \n",sep="")
                    }

                    if(object.robust && "se" %in% columns.indiv){
                        cat("Standard errors: robust\n")
                    }else if("se" %in% columns.indiv){
                        cat("Standard errors: model-based\n")
                    }
                    
                    if(any(c("p.value", "lower", "upper") %in% columns.indiv)){
                        if("p.value" %in% columns.indiv == FALSE){
                            txt.cip <- "P-values"
                        }else if("lower" %in% columns.indiv == FALSE && "upper" %in% columns.indiv == FALSE){
                            txt.cip <- "CIs"
                        }else{
                            txt.cip <- "CIs/p-values"
                        }
                        if(univariate.method[1] == "none"){ ## always only one hypothesis in each global test
                            cat("(",txt.cip," not adjusted for multiple comparisons) \n", sep="")
                        }else if(length(n.hypoPerTest)==1){ ## only one global test
                            if(univariate.method[1] == "bonferroni"){
                                cat("(",txt.cip," adjusted for multiple comparisons -- Bonferroni)\n", sep="")
                            }else if(univariate.method[1] %in% c("single-step", "single-step2")){
                                cat("(",txt.cip," adjusted for multiple comparisons -- max-test adjustment)\n", sep="")
                            }else{
                                cat(paste0("(",txt.cip," adjusted for multiple comparisons -- ",univariate.method[1],")\n", sep=""),sep="")
                            }
                        }else{
                            if(univariate.method[1] == "bonferroni"){
                                cat("(",txt.cip," adjusted for multiple comparisons within each global test -- bonferroni) \n", sep="")
                            }else if(univariate.method[1] %in% c("single-step","single-step2")){
                                cat("(",txt.cip," adjusted for multiple comparisons within each global test -- max-test adjustment) \n", sep="")
                            }else{
                                cat(paste0("(",txt.cip," adjusted for multiple comparisons within each global test -- ",univariate.method[1],") \n", sep=""),sep="")
                            }
                        }

                        if(univariate.method[1] == "single-step"){
                            error <- max(c(0,abs(attr(table.univariate,"error")[table.multivariate$type.original==iType])), na.rm = TRUE)
                            if(error > 1e-12){
                                txt.error <- paste0("Error when computing the adjusted ",txt.cip," by numerical integration: ", signif(error, digits = 5))
                                if(!is.null(seed)){
                                    txt.error <- paste0(txt.error," (seed ",seed,")")
                                }
                                cat(txt.error,".\n",sep="")
                            }
                        }else if(univariate.method[1] == "single-step2"){
                            txt.sample <- paste("Adjusted ",txt.cip," computed using ",attr(table.univariate,"n.sample")," samples", sep = "")
                            if(!is.null(seed)){
                                txt.sample <- paste0(txt.error," (seed ",seed,")")
                            }
                            cat(txt.sample,".\n")                            
                        }
                    }
                }
            }
        }

        if(print.global || print.indiv){ cat("\n") }

    }else if(attr(object,"test")=="LRT"){
        cat(" - Likelihood ratio test \n")
        out <- as.data.frame(object)
        out.print <- out
        if("null" %in% columns.indiv == FALSE){
            out.print[["null"]] <- NULL
        }
        if(print.global>0){
            print(out.print)
        }
    }

    return(invisible(NULL))
}


##----------------------------------------------------------------------
### summary.anova_lmm.R ends here
