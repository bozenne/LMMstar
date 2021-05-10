### summary.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: May 10 2021 (19:27) 
##           By: Brice Ozenne
##     Update #: 76
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.lmm (documentation)
##' @title Summary Output for a Linear Mixed Model 
##' @description Summary output for a linear mixed model fitted with \code{lmm}.
##' This is a modified version of the \code{nlme::summary.gls} function.
##' @name summary
##'
##' @param object [lmm] output of the \code{lmm} function.
##' @param digit [integer,>0] number of digit used to display numeric values.
##' @param level [numeric,0-1] confidence level for the confidence intervals.
##' @param print [logical] should the output be printed in the console.
##' @param hide.fit [logical] should information about the model fit not be printed.
##' @param hide.cor [logical] should information about the correlation structure not be printed.
##' @param hide.var [logical] should information about the variance structure not be printed.
##' @param hide.mean [logical] should information about the mean structure not be printed.
##' @param ... not used. For compatibility with the generic function.

## * summary.lmm (code)
##' @rdname summary
##' @export
summary.lmm <- function(object, digit = 3, level = 0.95, print = TRUE,
                        hide.fit = FALSE, hide.cor = FALSE, hide.var = FALSE, hide.sd = FALSE, hide.mean = FALSE, ...){

    ## ** welcome message
    if(print){
        if(length(object$param$rho) > 0){
            if(length(c(object$param$sigma,object$param$k))==1){
                cat("  Linear model \n")
            }else{
                cat("  Linear model with heterogeneous residual variance \n")
            }
        }else{
            if(object$structure=="UN"){
                cat("  Linear mixed effect model with an unstructured covariance matrix \n")
            }else if(object$structure=="CS"){
                cat("  Linear mixed effect model with a compound symmetry covariance matrix \n")
            }
        }
    }

    
    ## ** fit message
    if(print && !hide.fit){
        if(object$method == "REML"){
            cat("  - fitted using Restricted Maximum Likelihood (REML) \n")
        }else{
            cat("  - fitted using Maximum Likelihood (ML) \n")
        }
        cat("  - log-likelihood :", as.double(object$logLik), " (parameters: mean = ",length(object$param$mu),", variance = ",length(object$param[c("sigma","k")]),", correlation = ",length(object$param[c("rho")]),")\n",sep="")    
        cat(" \n")

        cat("Dataset:", deparse(object$call$data), "\n")
        cat(" - ", object$design$cluster$n, " clusters \n" , sep = "")
        cat(" - ", sum(object$design$cluster$nobs), " observations \n",  sep = "")
        cat(" - ", max(object$design$cluster$nobs), " maximum number of observations per cluster \n", sep = "")
        
        data.X <- object$data[all.vars(delete.response(terms(object$formula$mean)))]
        C <- lapply(data.X, function(iCol){
            if(inherits(iCol,"factor")){contrasts(iCol)}else if(inherits(iCol,"character")){contrasts(as.factor(iCol))}
        })
        C <- C[!unlist(lapply(C, is.null))]
        if(length(C)>0){
            cat(" - levels of the categorical variables \n", sep = "")
            if(attr(stats::terms(object$formula$mean),"intercept") == 1){
                ref.level <- paste(unlist(lapply(names(C), function(iC){
                    paste0(iC,"=",rownames(C[[iC]])[1])
                })), collapse = " ; ")
                cat(" - reference level: ",ref.level," \n", sep = "")
                cat(" \n")
            }
            print(C)
        }
    }

    ## ** correlation structure
    ## handle no correlation structure.
    if(any("rho" %in% object$param$type)){
        if(print && !hide.cor){
            cat("Correlation structure:",deparse(object$call$correlation),"\n")
        }
        vec.corcoef <- object$param$rho
        browser()
        ## TOFIX
        n.rep <- attr(object$modelStruct$corStruct,"Dim")$maxLen
        if(!is.null(object$modelStruct$varStruct)){
            name.rep <- attr(object$modelStruct$varStruct,"groupNames")
        }else{
            name.rep <- NULL
        }
        M.corcoef <- matrix(NA, nrow = n.rep, ncol = n.rep,
                            dimnames = list(name.rep,name.rep))
        M.corcoef[lower.tri(M.corcoef)] <- vec.corcoef
        M.corcoef[upper.tri(M.corcoef)] <- t(M.corcoef)[upper.tri(M.corcoef)]
        diag(M.corcoef) <- 1
        if(print && !hide.cor){
            print(M.corcoef, digit = digit)
            cat("\n")
        }
    }else{
        M.corcoef <- NULL
    }
    

    ## ** variance structure
    if(print && !hide.var && !hide.sd){
        cat("Variance structure:",deparse(object$formula$var.design),"\n")
    }

    M.varcoef <- rbind("variance" = stats::coef(object, effects = "variance", transform.k = "sd"),
                       "relative" = stats::coef(object, effects = "variance", transform.sigma = "one", transform.k = "none"))

    M.sdcoef <- rbind("standard deviation" = stats::coef(object, effects = "variance", transform.k = "var"), ## implicit transformation of sigma parameters when using var for the k parameters
                       "relative" = stats::coef(object, effects = "variance", transform.sigma = "one", transform.k = "square"))

    if(print){
        if(!hide.var && !hide.sd){
            print(rbind(M.varcoef,M.sdcoef), digit = digit)
            cat("\n")
        }else if(!hide.var){
            print(M.varcoef, digit = digit)
            cat("\n")
        }else if(!hide.sd){
            print(M.sdcoef, digit = digit)
            cat("\n")
        }
    }

    ## ** mean structure
    if(print && !hide.mean){
        cat("Mean structure:",deparse(object$call$formula),"\n")
    }
    tTable <- confint(object, level = level, effects = "mean")[,c("estimate","se","df","lower","upper","p.value")]
    starSymbol <- stats::symnum(tTable[,"p.value"], corr = FALSE, na = FALSE,
                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                symbols = c("***", "**", "*", ".", " "))
    tTable.print <- cbind(tTable, starSymbol)
    names(tTable.print)[7] <- ""
    tTable.print[["p.value"]] <- format.pval(tTable[["p.value"]], digits = digit, eps = 10^(-digit))
    tTable.print[["estimate"]] <- as.character(round(tTable[["estimate"]], digits = digit))
    tTable.print[["df"]] <- as.character(round(tTable[["df"]], digits = digit))
    tTable.print[["lower"]] <- as.character(round(tTable[["lower"]], digits = digit))
    tTable.print[["upper"]] <- as.character(round(tTable[["lower"]], digits = digit))
    tTable.print[["se"]] <- as.character(round(tTable[["se"]], digits = digit))
    if(print && !hide.mean){
        print(tTable.print)
        cat("\n")
        cat("The columns lower and upper correspond to the ",100*level,"% confidence interval of the estimated coefficient\n", sep = "")
        if(NROW(tTable.print)>1){
            cat("Note: p-values and confidence intervals are not adjusted for multiple comparisons\n")
        }
    }
    ## ** export
    return(invisible(list(correlation = M.corcoef,
                          variance = M.varcoef,
                          mean = tTable)))
}


######################################################################
### summary.lmm.R ends here
