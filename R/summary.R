### summary.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: mar 22 2021 (23:48) 
##           By: Brice Ozenne
##     Update #: 51
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
##' @param hide.fit [logical] should information about the model fit not be printed.
##' @param hide.cor [logical] should information about the correlation structure not be printed.
##' @param hide.var [logical] should information about the variance structure not be printed.
##' @param hide.mean [logical] should information about the mean structure not be printed.
##' @param ... not used. For compatibility with the generic function.
##' @export
summary.lmm <- function(object, digit = 3, conf.level = 0.95, print = TRUE,
                        hide.fit = FALSE, hide.cor = FALSE, hide.var = FALSE, hide.sd = FALSE, hide.mean = FALSE, ...){

    ## ** welcome message
    if(print){
        if(max(object$design$cluster$nobs)==1){
            if(length(c(object$param$sigma,object$param$k))==1){
                cat("  Linear model \n")
            }else if(object$structure=="CS"){
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
        cat("  - likelihood :", as.double(object$logLik), " (df: mean = ",length(object$param$mu),", variance = ",length(object$param[c("sigma","k","cor")]),")\n",sep="")    
        cat(" \n")

        cat("Dataset:", deparse(object$call$data), "\n")
        cat(" - ", object$design$cluster$n, " clusters \n" , sep = "")
        cat(" - ", sum(object$design$cluster$nobs), " observations \n",  sep = "")
        cat(" - ", max(object$design$cluster$nobs), " maximum number of observations per cluster \n", sep = "")
        if(length(object$contrasts)>0){
            cat(" - levels of the categorical variables \n", sep = "")

            data.X <- object$data[all.vars(delete.response(terms(object$formula$mean)))]
            C <- lapply(data.X, function(iCol){
                if(inherits(iCol,"factor")){contrasts(iCol)}else if(inherits(iCol,"character")){contrasts(as.factor(iCol))}
            })
            C <- C[!unlist(lapply(C, is.null))]
            print(C)

            if(attr(stats::terms(object$formula$mean),"intercept") == 1){
                ref.level <- paste(unlist(lapply(names(C), function(iC){
                    paste0(iC,"=",rownames(C[[iC]])[1])
                })), collapse = " ; ")
                cat(" - reference level: ",ref.level," \n", sep = "")
                cat(" \n")
            }
        }else{
            cat(" \n")

        }
    }

    ## ** correlation structure
    ## handle no correlation structure.
    if(print && !hide.cor){
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
    if(print && !hide.cor){
        print(M.corcoef, digit = digit)
        cat("\n")
    }

    ## ** variance structure
    if(print && !hide.var && !hide.sd){
        cat("Variance structure:",deparse(object$call$weights),"\n")
    }
    M.varcoef <- matrix(NA, nrow = 2, ncol = n.rep,
                        dimnames = list(c("variance","relative"),name.rep))
    if(!is.null(object$modelStruct$varStruct)){
        vec.varcoef <- stats::coef(object$modelStruct$varStruct, unconstrained = FALSE)^2
        M.varcoef[1,] <- stats::sigma(object)^2*c(1,vec.varcoef)
        M.varcoef[2,] <- c(1,vec.varcoef)
    }else{
        M.varcoef[1,] <- stats::sigma(object)^2
        M.varcoef[2,] <- 1
    }

    M.sdcoef <- matrix(NA, nrow = 2, ncol = n.rep,
                       dimnames = list(c("standard deviation","relative"),name.rep))
    if(!is.null(object$modelStruct$varStruct)){
        vec.sdcoef <- stats::coef(object$modelStruct$varStruct, unconstrained = FALSE)
        M.sdcoef[1,] <- stats::sigma(object)*c(1,vec.varcoef)
        M.sdcoef[2,] <- c(1,vec.sdcoef)
    }else{
        M.sdcoef[1,] <- stats::sigma(object)
        M.sdcoef[2,] <- 1
    }
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
        cat("Mean structure:",deparse(object$call$model),"\n")
    }
    object2 <- object
    class(object2) <- setdiff(class(object2),"lmm")
    tTable <- summary(object2, verbose=FALSE)$tTable
    tTable <- data.frame(estimate = tTable[,"Value"],
                         se = tTable[,"Std.Error"],
                         nlme::intervals(object, which = "coef", level = conf.level)$coef[,c("lower","upper")],
                         t.value = tTable[,"t-value"],
                         p.value = tTable[,"p-value"])
    starSymbol <- stats::symnum(tTable[,"p.value"], corr = FALSE, na = FALSE,
                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                symbols = c("***", "**", "*", ".", " "))
    tTable.print <- cbind(tTable, starSymbol)
    names(tTable.print)[7] <- ""
    tTable.print[["p.value"]] <- format.pval(tTable[["p.value"]], digits = digit, eps = 10^(-digit))
    tTable.print[["estimate"]] <- as.character(round(tTable[["estimate"]], digits = digit))
    tTable.print[["lower"]] <- as.character(round(tTable[["lower"]], digits = digit))
    tTable.print[["upper"]] <- as.character(round(tTable[["lower"]], digits = digit))
    tTable.print[["se"]] <- as.character(round(tTable[["se"]], digits = digit))
    if(print && !hide.mean){
        print(tTable.print)
        cat("\n")
        cat("The columns lower and upper correspond to the ",100*conf.level,"% confidence interval of the estimated coefficient\n", sep = "")
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
