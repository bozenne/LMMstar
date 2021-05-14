### summary.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: May 14 2021 (17:10) 
##           By: Brice Ozenne
##     Update #: 172
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
##' @param ci [logical] should the confidence intervals be displayed.
##' @param hide.fit [logical] should information about the model fit not be printed.
##' @param hide.cor [logical] should information about the correlation structure not be printed.
##' @param hide.var [logical] should information about the variance structure not be printed.
##' @param hide.mean [logical] should information about the mean structure not be printed.
##' @param ... not used. For compatibility with the generic function.

## * summary.lmm (code)
##' @rdname summary
##' @export
summary.lmm <- function(object, digit = 3, level = 0.95, print = TRUE, ci = FALSE,
                        hide.fit = FALSE, hide.cor = FALSE, hide.var = TRUE, hide.sd = FALSE, hide.mean = FALSE, ...){

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
    if(any("rho" %in% object$param$type)){
        if(print && !hide.cor){
            cat("Correlation structure:",deparse(object$formula$cor),"\n")
        }
        vec.corcoef <- object$param$rho
        table.cor <- lapply(getVarCov(object, simplifies = FALSE),cov2cor)
        if(length(table.cor)==1){
            table.cor <- table.cor[[1]]
        }
        if(print && !hide.cor){
            print(table.cor, digit = digit)
            cat("\n")
        }
    }else{
        table.cor <- NULL
    }
    

    ## ** variance structure
    if(print && (!hide.var || !hide.sd)){
        cat("Variance structure:",deparse(object$formula$var.design),"\n")
        name.sigma <- names(coef(object, transform.k = "sd", effects = "variance"))
        index.ref <- which(names(coef(object, effects = "variance", transform.names = FALSE)) %in% names(object$param$sigma))
    }

    if(!hide.var){
        M.varcoef <- confint(object, level = level, df = ci, transform.k = "var", effects = "variance")
        M.varcoefRe <- confint(object, level = level, df = ci, transform.sigma = "none", transform.k = "square", effects = "variance", transform.names = FALSE)
        M.varcoefRe[index.ref,c("estimate","lower","upper")] <- 1
        M.varcoefRe[index.ref,c("se","sd")] <- NA

        if(ci){
            table.var <- data.frame(estimate = M.varcoef[,"estimate"], lower = M.varcoef[,"lower"], upper = M.varcoef[,"upper"],
                                    estimate.ratio = M.varcoefRe[,"estimate"], lower.ratio = M.varcoefRe[,"lower"], upper.ratio = M.varcoefRe[,"upper"])
        }else{
            table.var <- data.frame(estimate = M.varcoef[,"estimate"], estimate.ratio = M.varcoefRe[,"estimate"])
        }
        rownames(table.var) <- name.sigma
        test.k <- NROW(table.var) > length(index.ref)
    }else{
        table.var <- NULL
    }
    if(!hide.sd){
        M.sdcoef <- confint(object, level = level, df = ci, transform.k = "sd", effects = "variance")
        M.sdcoefRe <- confint(object, level = level, df = ci, transform.sigma = "none", transform.k = "none", effects = "variance", transform.names = FALSE)
        M.sdcoefRe[index.ref,c("estimate","lower","upper")] <- 1
        M.sdcoefRe[index.ref,c("se","sd")] <- NA

        if(ci){
            table.sd <- data.frame(estimate = M.sdcoef[,"estimate"], lower = M.sdcoef[,"lower"], upper = M.sdcoef[,"upper"],
                                   estimate.ratio = M.sdcoefRe[,"estimate"], lower.ratio = M.sdcoefRe[,"lower"], upper.ratio = M.sdcoefRe[,"upper"])
        }else{
            table.sd <- data.frame(estimate = M.sdcoef[,"estimate"], estimate.ratio = M.sdcoefRe[,"estimate"])
        }
        rownames(table.sd) <- name.sigma
        test.k <- NROW(table.sd) > length(index.ref)
    }else{
        table.sd <- NULL
    }
    
    if(print && (!hide.var || !hide.sd)){
        if(!ci){
            printtable <- matrix(NA, ncol = 0, nrow = length(name.sigma))
            if(!hide.var){
                printtable <- cbind(printtable, data.frame(variance = unname(table.var[,"estimate"])))
                if(test.k){
                    printtable <- cbind(printtable, data.frame(ratio = unname(table.var[,"estimate.ratio"])))
                }
            }
            if(!hide.sd){
                printtable <- cbind(printtable, data.frame("standard deviation" = unname(table.sd[,"estimate"])))
                if(test.k){
                    printtable <- cbind(printtable, data.frame(ratio = unname(table.sd[,"estimate.ratio"])))
                }
            }

            rownames(printtable) <- name.sigma
            print(printtable, digit = digit)
            cat("\n")
        }else{
            if(!hide.var){
                printtable.var <- table.var
                printtable.var <- apply(printtable.var, 2, round, digit = digit)
                if(test.k){
                    cat("Variance estimates (relative to reference): \n")
                    printtable.var[index.ref,"estimate.ratio"] <- "reference"
                    printtable.var[index.ref,c("lower.ratio","upper.ratio")] <- ""
                    printtable.var[-index.ref,"estimate"] <- paste0(printtable.var[-index.ref,"estimate"], " (",printtable.var[-index.ref,"estimate.ratio"],")")
                    printtable.var[-index.ref,"lower"] <- paste0(printtable.var[-index.ref,"lower"], " (",printtable.var[-index.ref,"lower.ratio"],")")
                    printtable.var[-index.ref,"upper"] <- paste0(printtable.var[-index.ref,"upper"], " (",printtable.var[-index.ref,"upper.ratio"],")")
                }else{
                    cat("Variance estimates: \n")
                }
                print(printtable.var[,c("estimate","lower","upper")], quote = FALSE)
                cat("\n")
            }
            if(!hide.sd){
                printtable.sd <- table.sd
                printtable.sd <- apply(printtable.sd,2,round,digit=digit)
                if(test.k){
                    cat("Standard deviation estimates (relative to reference): \n")
                    printtable.sd[index.ref,"estimate.ratio"] <- "reference"
                    printtable.sd[index.ref,c("lower.ratio","upper.ratio")] <- ""
                    printtable.sd[-index.ref,"estimate"] <- paste0(printtable.sd[-index.ref,"estimate"], " (",printtable.sd[-index.ref,"estimate.ratio"],")")
                    printtable.sd[-index.ref,"lower"] <- paste0(printtable.sd[-index.ref,"lower"], " (",printtable.sd[-index.ref,"lower.ratio"],")")
                    printtable.sd[-index.ref,"upper"] <- paste0(printtable.sd[-index.ref,"upper"], " (",printtable.sd[-index.ref,"upper.ratio"],")")
                }else{
                    cat("Standard deviation estimates: \n")
                }
                print(printtable.sd[,c("estimate","lower","upper")], quote = FALSE)
                cat("\n")
            }
        }
    }

    ## ** mean structure
    if(print && !hide.mean){
        cat("Mean structure:",deparse(object$call$formula),"\n")
    }
    table.mean <- confint(object, level = level, effects = "mean")[,c("estimate","se","df","lower","upper","statistic","p.value")]
    starSymbol <- stats::symnum(table.mean[,"p.value"], corr = FALSE, na = FALSE,
                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                symbols = c("***", "**", "*", ".", " "))
    printtable.mean <- cbind(table.mean, starSymbol)
    names(printtable.mean)[8] <- ""
    printtable.mean[["estimate"]] <- as.character(round(table.mean[["estimate"]], digits = digit))
    printtable.mean[["se"]] <- as.character(round(table.mean[["se"]], digits = digit))
    printtable.mean[["df"]] <- as.character(round(table.mean[["df"]], digits = digit))
    printtable.mean[["lower"]] <- as.character(round(table.mean[["lower"]], digits = digit))
    printtable.mean[["upper"]] <- as.character(round(table.mean[["lower"]], digits = digit))
    printtable.mean[["statistic"]] <- as.character(round(table.mean[["statistic"]], digits = digit))
    printtable.mean[["p.value"]] <- format.pval(table.mean[["p.value"]], digits = digit, eps = 10^(-digit))
    if(print && !hide.mean){
        if(ci){
            printtable.mean$statistic <- NULL
            print(printtable.mean)
            cat("\n")
            cat("The columns lower and upper correspond to the ",100*level,"% confidence interval of the estimated coefficient\n", sep = "")
            if(NROW(printtable.mean)>1){
                cat("Note: p-values and confidence intervals are not adjusted for multiple comparisons\n")
            }
        }else{
            printtable.mean$lower <- NULL
            printtable.mean$upper <- NULL
            print(printtable.mean)
            cat("\n")
            if(NROW(printtable.mean)>1){
                cat("Note: p-values are not adjusted for multiple comparisons\n")
            }
        }
    }
    ## ** export
    return(invisible(list(correlation = table.cor,
                          variance = table.var,
                          sd = table.sd,
                          mean = table.mean)))
}


######################################################################
### summary.lmm.R ends here
