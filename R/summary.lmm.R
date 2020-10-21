### summary.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: okt 21 2020 (15:32) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.lmm
##' @title Summary Output for a Linear Mixed Model 
##' @description Summary output for a linear mixed model fitted with \code{lmm}.
##' This is a modified version of the \code{nlme::summary.gls} function.
##'
##' @param object [lmm] output of the \code{lmm} function.
##' @param digit [integer,>0] number of digit used to display numeric values.
##' @param conf.level [numeric,0-1] confidence level for the confidence intervals.
##' @param print [logical] should the output be printed in the console.
##' @param ... not used. For compatibility with the generic function.
##' @export
summary.lmm <- function(object, digit = 3, conf.level = 0.95, print = TRUE, ...){

    ## ** welcome message
    if(print){
        if(!is.null(object$modelStruct$varStruct)){
            cat("  Linear mixed effect model with an unstructured covariance matrix \n")
        }else{
            cat("  Linear mixed effect model with a compound symmetry covariance matrix \n")
        }
        if(object$dim$REML){
            cat("  - fitted using Restricted Maximum Likelihood (REML) \n")
        }else{
            cat("  - fitted using Maximum Likelihood (REML) \n")
        }
        cat("  - likelihood :", as.double(object$logLik), " (df = ",object$dim$p,")\n",sep="")    
        cat(" \n")

        cat("Dataset:", deparse(object$call$data), "\n")
        cat(" - ", attr(object$modelStruct$corStruct,"Dim")$M, " clusters \n" , sep = "")
        cat(" - ", attr(object$modelStruct$corStruct,"Dim")$N, " observations \n",  sep = "")
        cat(" - ", attr(object$modelStruct$corStruct,"Dim")$maxLen, " maximum number of observations per cluster \n", sep = "")
        if(length(object$contrasts)>0){
            cat(" - levels of the categorical variables \n", sep = "")
            print(object$contrasts)
        
            if(attr(stats::terms(eval(object$call$model)),"intercept") == 1){
                ref.level <- paste(unlist(lapply(names(object$contrasts), function(iC){
                    paste0(iC,"=",rownames(object$contrasts[[iC]])[1])
                })), collapse = " ; ")
                cat(" - reference level: ",ref.level," \n", sep = "")
                cat(" \n")
            }
        }else{
            cat(" \n")

        }
    }

    ## ** correlation structure
    if(print){
        cat("Correlation structure:",deparse(object$call$correlation),"\n")
    }
    vec.corcoef <- stats::coef(object$modelStruct$corStruct, unconstrained = FALSE)

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
    if(print){
        print(M.corcoef, digit = digit)
        cat("\n")
    }

    ## ** variance structure
    if(print){
        cat("Variance structure:",deparse(object$call$weights),"\n")
    }
    M.varcoef <- matrix(NA, nrow = 2, ncol = n.rep,
                        dimnames = list(c("variance","relative variance"),name.rep))
    if(!is.null(object$modelStruct$varStruct)){
        vec.varcoef <- stats::coef(object$modelStruct$varStruct, unconstrained = FALSE)^2
        M.varcoef[1,] <- stats::sigma(object)^2*c(1,vec.varcoef)
        M.varcoef[2,] <- c(1,vec.varcoef)
    }else{
        M.varcoef[1,] <- stats::sigma(object)^2
        M.varcoef[2,] <- 1
    }
    if(print){
    print(M.varcoef, digit = digit)
    cat("\n")
    }

    ## ** mean structure
    if(print){
        cat("Mean structure:",deparse(object$call$model),"\n")
    }
    object2 <- object
    class(object2) <- setdiff(class(object2),"lmm")
    tTable <- summary(object2, verbose=FALSE)$tTable
    starSymbol <- stats::symnum(tTable[,4], corr = FALSE, na = FALSE,
                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                symbols = c("***", "**", "*", ".", " "))
    tTable <- data.frame(tTable[,1], nlme::intervals(object, which = "coef", level = conf.level)$coef[,c("lower","upper")],tTable[,2:4],starSymbol)
    names(tTable) <- c("estimate","lower","upper","se","t-value","p-value","")
    tTable.print <- tTable[,-5]
    tTable.print[["p-value"]] <- format.pval(tTable[["p-value"]], digits = digit, eps = 10^(-digit))
    tTable.print[["estimate"]] <- as.character(round(tTable[["estimate"]], digits = digit))
    tTable.print[["lower"]] <- as.character(round(tTable[["lower"]], digits = digit))
    tTable.print[["upper"]] <- as.character(round(tTable[["lower"]], digits = digit))
    tTable.print[["se"]] <- as.character(round(tTable[["se"]], digits = digit))
    if(print){
        print(tTable.print)
        cat("\n")
        cat("The columns lower and upper correspond to the ",100*conf.level,"% confidence interval of the estimated coefficient\n", sep = "")
        cat("Note: p-value(s) and confidence interval(s) are not adjusted for multiple comparisons\n")
    }
    ## ** export
    return(invisible(list(correlation = M.corcoef,
                          variance = M.varcoef,
                          mean = tTable)))
}


######################################################################
### summary.lmm.R ends here
