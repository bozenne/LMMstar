### print.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 24 2022 (10:05) 
## Version: 
## Last-Updated: jan 24 2022 (16:45) 
##           By: Brice Ozenne
##     Update #: 19
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.anova_lmm
##' @rdname anova
##' @export
print.anova_lmm <- function(x, ...){
    dots <- list(...)
    dots$print.null <- NULL
    return(do.call(summary, c(list(object = x, print.null = -1), dots)))
}

## * summary.anova_lmm
##' @rdname anova
##' @export
summary.anova_lmm <- function(object, level = 0.95, method = "single-step", transform = NULL, print.null = FALSE, columns = NULL, ...){

    if(attr(object,"test")=="Wald"){
        type <- setdiff(names(object),"call")
        ci <- stats::confint(object, level = level, method = method)
        for(iType in type){

            if(is.null(object[[iType]])){next}

            if(!print.null){
                object[[iType]][["null"]] <- NULL
            }
            iNoDf <- is.infinite(object[[iType]]$df.denom)
            txt.test <- ifelse(all(iNoDf),"Chi-square test","F-test")
            if(iType == "all"){
                if(print.null>=0){cat("                     ** User-specified hypotheses ** \n", sep="")}
                if(is.na(object[[iType]]$statistic) && !is.null(attr(object[[iType]]$statistic,"error"))){
                    cat(" - ",txt.test,": ",attr(object[[iType]]$statistic,"error"),"\n", sep="")
                }else{
                    if(print.null>=0){cat(" - ",txt.test,"\n", sep="")}
                    print(object[[iType]], row.names = FALSE)
                }
            }else{
                if(print.null>=0){
                    cat("                     ** ",iType," coefficients ** \n", sep="")
                    cat(" - ",txt.test,"\n",sep="")
                }
                print(object[[iType]])
            }
            if(!is.null(ci[[iType]]) && print.null>=0){
                options <- LMMstar.options()
                if(is.null(columns)){
                    columns <- options$columns.anova
                }
                if(all(sapply(ci[[iType]],NROW)==1) || method == "none"){ ## always only one hypothesis in each global test
                    cat("\n - P-values and confidence interval \n", sep="")
                }else if(length(ci[[iType]])==1){ ## only one global test
                    cat("\n - P-values and confidence interval (adjusted for multiplicity) \n", sep="")
                }else{
                    cat("\n - P-values and confidence interval (adjusted for multiplicity within each global test) \n", sep="")
                }
                if(!is.null(transform)){
                    ci[[iType]] <- lapply(ci[[iType]], function(iCI){transformSummaryTable(iCI, transform = transform)})
                }                    
                print(do.call(rbind,unname(ci[[iType]]))[,columns,drop=FALSE])
                if(!is.null(transform)){
                    cat("\n Columns ",paste(setdiff(columns,c("df","p.value")), collapse =", ")," have been back-transform. \n",sep="")
                }                    
            }
            cat("\n")
        }
    }else if(attr(object,"test")=="LRT"){
        cat(" - Likelihood ratio test \n")
        object.print <- as.data.frame(object)
        if(print.null==FALSE){object.print[["null"]] <- NULL}
        print(object.print)
    }

    return(invisible(NULL))
}
##----------------------------------------------------------------------
### print.anova_lmm.R ends here
