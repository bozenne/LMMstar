### summary.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: May 27 2021 (17:06) 
##           By: Brice Ozenne
##     Update #: 195
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
##' @param hide.sd [logical] should information about the standard deviation not be printed.
##' @param hide.var [logical] should information about the variance not be printed.
##' @param hide.mean [logical] should information about the mean structure not be printed.
##' @param ... not used. For compatibility with the generic function.

## * summary.lmm (code)
##' @rdname summary
##' @export
summary.lmm <- function(object, digit = 3, level = 0.95, print = TRUE, ci = FALSE,
                        hide.fit = FALSE, hide.cor = FALSE, hide.var = TRUE, hide.sd = FALSE, hide.mean = FALSE, ...){

    param.mu <- object$param$value[object$param$type=="mu"]
    param.sigma <- object$param$value[object$param$type=="sigma"]
    param.k <- object$param$value[object$param$type=="k"]
    param.rho <- object$param$value[object$param$type=="rho"]
    data <- object$data
    call <- object$call
    structure <- object$structure
    logLik <- stats::logLik(object)
    nobs <- stats::nobs(object)
    method.fit <- object$method
    nobsByCluster <- object$design$cluster$nobs
    formula <- object$formula
    Omega <- getVarCov(object, simplifies = FALSE)
    df <- !is.null(object$df)
    
    ## ** welcome message
    if(print){
        if(length(param.rho) > 0){
            if(length(c(param.sigma,param.k))==1){
                cat("  Linear model \n")
            }else{
                cat("  Linear model with heterogeneous residual variance \n")
            }
        }else{
            if(structure=="UN"){
                cat("  Linear mixed effect model with an unstructured covariance matrix \n")
            }else if(structure=="CS"){
                cat("  Linear mixed effect model with a compound symmetry covariance matrix \n")
            }
        }
    }

    
    ## ** fit message
    if(print && !hide.fit){
        if(method.fit == "REML"){
            cat("  - fitted using Restricted Maximum Likelihood (REML) \n")
        }else{
            cat("  - fitted using Maximum Likelihood (ML) \n")
        }
        cat("  - log-likelihood :", as.double(logLik), " (parameters: mean = ",length(param.mu),", variance = ",length(c(param.sigma,param.k)),", correlation = ",length(param.rho),")\n",sep="")    
        cat(" \n")

        cat("Dataset:", deparse(call$data), "\n")
        cat(" - ", nobs["cluster"], " clusters \n" , sep = "")
        cat(" - ", sum(nobsByCluster), " observations \n",  sep = "")
        cat(" - ", max(nobsByCluster), " maximum number of observations per cluster \n", sep = "")

        data.X <- data[all.vars(delete.response(terms(formula$mean)))]
        C <- lapply(data.X, function(iCol){
            if(inherits(iCol,"factor")){contrasts(iCol)}else if(inherits(iCol,"character")){contrasts(as.factor(iCol))}
        })
        C <- C[!unlist(lapply(C, is.null))]
        if(length(C)>0){
            cat(" - levels of the categorical variables \n", sep = "")
            if(attr(stats::terms(formula$mean),"intercept") == 1){
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
    if(length(param.rho)>0){
        if(print && !hide.cor){
            cat("Correlation structure:",deparse(formula$cor),"\n")
        }
        table.cor <- lapply(Omega,cov2cor)
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
        cat("Variance structure:",deparse(formula$var.design),"\n")
        name.sigma <- names(coef(object, transform.k = "sd", effects = "variance"))
        index.ref <- which(names(coef(object, effects = "variance", transform.names = FALSE)) %in% names(param.sigma))
    }

    if(!hide.var){
        M.varcoef <- confint(object, level = level, df = df, transform.k = "var", effects = "variance")
        M.varcoefRe <- confint(object, level = level, df = df, transform.sigma = "none", transform.k = "square", effects = "variance", transform.names = FALSE)
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
        M.sdcoef <- confint(object, level = level, df = df, transform.k = "logsd", effects = "variance")
        M.sdcoefRe <- confint(object, level = level, df = df, transform.sigma = "log", effects = "variance", transform.names = FALSE)

        M.sdcoefRe[index.ref,c("estimate","lower","upper")] <- 1
        M.sdcoefRe[index.ref,c("se")] <- NA

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
                for(iCol in 1:NCOL(printtable.sd))
                    printtable.sd[,iCol] <- round(printtable.sd[,iCol], digits = digit)
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
                if(object$reparametrize$transform.sigma=="log" && (all(object$param$type!="k") || (object$reparametrize$transform.k=="log"))){
                    cat("\n Note: confidence intervals have been computed on the log-scale and then back-transformed.\n",sep="")
                }
                cat("\n")
            }
        }
    }

    ## ** mean structure
    if(print && !hide.mean){
        cat("Mean structure:",deparse(call$formula),"\n")
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
