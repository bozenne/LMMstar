### summary.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:13) 
## Version: 
## Last-Updated: okt 20 2021 (15:58) 
##           By: Brice Ozenne
##     Update #: 356
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
##' @param columns [character vector] Columns to be output for the fixed effects. Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param robust [logical] Should robust standard error (aka sandwich estimator) be output instead of the model-based standard errors. 
##' @param hide.fit [logical] should information about the model fit not be printed.
##' @param hide.data [logical] should information about the dataset not be printed.
##' @param hide.cor [logical] should information about the correlation structure not be printed.
##' @param hide.sd [logical] should information about the standard deviation not be printed.
##' @param hide.var [logical] should information about the variance not be printed.
##' @param hide.mean [logical] should information about the mean structure not be printed.
##' @param ... not used. For compatibility with the generic function.
##'
##' @return A list containing elements displayed in the summary: \itemize{
##' \item \code{correlation}: the correlation structure.
##' \item \code{variance}: the variance structure.
##' \item \code{sd}: the variance structure expressed in term of standard deviations.
##' \item \code{mean}: the mean structure.
##' }

## * summary.lmm (code)
##' @rdname summary
##' @export
summary.lmm <- function(object, digit = 3, level = 0.95, robust = FALSE, print = TRUE, columns = NULL,
                        hide.fit = FALSE, hide.data = FALSE, hide.cor = FALSE, hide.var = TRUE, hide.sd = FALSE, hide.mean = FALSE, ...){

    ## ** normalize user input
    param.mu <- object$param$value[object$param$type=="mu"]
    param.sigma <- object$param$value[object$param$type=="sigma"]
    param.k <- object$param$value[object$param$type=="k"]
    param.rho <- object$param$value[object$param$type=="rho"]
    data <- object$data
    call <- object$call
    structure <- object$design$vcov
    logLik <- stats::logLik(object)
    nobs <- stats::nobs(object)
    method.fit <- object$method
    nobsByCluster <- object$design$cluster$nobs
    formula <- object$formula
    Omega <- nlme::getVarCov(object, simplifies = FALSE)
    df <- !is.null(object$df)
    options <- LMMstar.options()
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    if(!is.null(columns)){
        columns  <- match.arg(columns, c("estimate","se","statistic","df","lower","upper","null","p.value",""), several.ok = TRUE)
    }else{
        columns <- options$columns.summary
    }


    ## ** welcome message
    if(print){
        if(length(param.rho) == 0){
            cat("           Linear regression \n")
        }else{
            cat("           Linear Mixed Model \n")
        }
        cat(" \n")
    }

    
    ## ** data message    
    if(print && !hide.data){
        cat("Dataset:", deparse(call$data), "\n\n")

        if(nobs["missing"]>0){
            missing.cluster <- .addNA(object$index.na, design = object$design, time = object$time)$missing.cluster
            if(length(missing.cluster)>0){
                cat("  - ", nobs["cluster"], " clusters were analyzed, ",length(missing.cluster)," were excluded because of missing values \n" , sep = "")
            }else{
                cat("  - ", nobs["cluster"], " clusters \n" , sep = "")
            }
            cat("  - ", sum(nobsByCluster), " observations were analyzed, ",nobs["missing"]," were excluded because of missing values \n",  sep = "")
        }else{
            cat("  - ", nobs["cluster"], " clusters \n" , sep = "")
            cat("  - ", sum(nobsByCluster), " observations \n",  sep = "")
        }
        if(length(unique(nobsByCluster))==1){
            cat("  - ", nobsByCluster[1], " observations per cluster \n", sep = "")
        }else{
            cat("  - between ", min(nobsByCluster), " and ",max(nobsByCluster)," observations per cluster \n", sep = "")
        }

        cat("\nSummary of the outcome and covariates: \n\n")
        data.XY <- data[all.vars(stats::terms(formula$mean))]
        str.XY <- utils::capture.output(utils::str(data.XY))[-1]
        str.XY[1] <- paste0(" ",str.XY[1])
        cat(paste0("  ",str.XY,"\n"))

        reference.level <- levels(object)$reference
        if(!is.null(reference.level)){
            cat("    reference level: ",paste(paste(names(reference.level),reference.level,sep="="),collapse=";")," \n", sep = "")
        }
        cat("\n")
    }
    
    ## ** optim message    
    if(print && !hide.fit){
        cat("Estimation procedure \n\n")
        if(method.fit == "REML"){
            cat("  - Restricted Maximum Likelihood (REML) \n")
        }else{
            cat("  - Maximum Likelihood (ML) \n")
        }
        cat("  - log-likelihood :", as.double(logLik), "\n",sep="")
        cat("  - parameters: mean = ",length(param.mu),", variance = ",length(c(param.sigma,param.k)),", correlation = ",length(param.rho),"\n", sep = "")
        if(!is.null(object$opt)){
            index.score <- which.max(abs(object$score))
            cat("  - convergence: ",object$opt$cv," (",object$opt$n.iter," iterations, largest |score|=",object$score[index.score]," is for ",names(index.score),")\n", sep = "")
        }
        cat(" \n")
    }

    ## ** vcov structure
    if(print && (!hide.cor || !hide.var || !hide.sd)){
        cat("Residual variance-covariance: ")
        if(length(param.rho)==0){
            if(length(c(param.sigma,param.k))==1){
                cat("identity \n\n")
            }else{
                cat("diagonal \n\n")
            }
        }else if(structure$type == "UN"){
            cat("unstructured \n\n")
        }else if(structure$type == "CS"){
            cat("compound symmetry \n\n")
        }
    }
    ## *** correlation
    if(object$time$n>1){
        if(print && !hide.cor){
            cat("  - correlation structure:",deparse(formula$cor),"\n")
        }
        table.cor <- lapply(Omega,stats::cov2cor)
        if(length(table.cor)==1){
            table.cor <- table.cor[[1]]
        }
        if(print && !hide.cor){
            if(is.list(table.cor)){
                table.cor  <- lapply(table.cor, function(iCor){
                    rownames(iCor) <- paste0("    ",rownames(iCor))
                    return(iCor)
                })
            }else{
                rownames(table.cor) <- paste0("    ",rownames(table.cor))
            }
            print(table.cor, digit = digit)
            cat("\n")
        }
    }else{
        table.cor <- NULL
    }
    

    ## *** variance
    if(!hide.var || !hide.sd){
        name.sigma <- gsub("sigma:","",names(coef(object, transform.sigma = "none", transform.k = "sd", effects = "variance")))
        index.ref <- which(names(coef(object, effects = "variance", transform.names = FALSE)) %in% names(param.sigma))
        if(print){
            cat("  - variance structure:",deparse(formula$var.design),"\n")
        }
    }
    if(!hide.var){
        table.var <- cbind(estimate = coef(object, transform.sigma = "none", transform.k = "var", effects = "variance"),
                           estimate.ratio = coef(object, transform.sigma = "none", transform.k = "square", effects = "variance", transform.names = FALSE))
        table.var[index.ref,"estimate.ratio"] <- 1
        test.k <- NROW(table.var) > length(index.ref)
        rownames(table.var) <- name.sigma
    }else{
        table.var <- NULL
    }
    if(!hide.sd){
        table.sd <- cbind(estimate = coef(object, transform.sigma = "none", transform.k = "sd", effects = "variance"),
                          estimate.ratio = coef(object, transform.sigma = "none", transform.k = "none", effects = "variance", transform.names = FALSE))
        table.sd[index.ref,"estimate.ratio"] <- 1
        test.k <- NROW(table.sd) > length(index.ref)
        rownames(table.sd) <- name.sigma
    }else{
        table.sd <- NULL
    }
    if(print && (!hide.var || !hide.sd)){
            printtable <- matrix(NA, ncol = 0, nrow = length(name.sigma))
            if(!hide.var){
                printtable <- cbind(printtable, data.frame(variance = unname(table.var[,"estimate"]),stringsAsFactors = FALSE))
                if(test.k){
                    printtable <- cbind(printtable, data.frame(ratio = unname(table.var[,"estimate.ratio"]),stringsAsFactors = FALSE))
                }
            }
            if(!hide.sd){
                printtable <- cbind(printtable, data.frame("standard deviation" = unname(table.sd[,"estimate"]),stringsAsFactors = FALSE))
                if(test.k){
                    printtable <- cbind(printtable, data.frame(ratio = unname(table.sd[,"estimate.ratio"]),stringsAsFactors = FALSE))
                }
            }

            rownames(printtable) <- paste0("    ",name.sigma)
            print(printtable, digit = digit)
    }
    if(print && (!hide.cor || !hide.var || !hide.sd)){
        cat("\n")
    }
    
    ## ** mean structure
    if(print && !hide.mean){
        cat("Fixed effects:",deparse(call$formula),"\n\n")
    }
    table.mean <- confint(object, level = level, robust = robust, effects = "mean", columns = c("estimate","se","df","lower","upper","statistic","p.value"))
    starSymbol <- stats::symnum(table.mean[,"p.value"], corr = FALSE, na = FALSE,
                                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                symbols = c("***", "**", "*", ".", " "))
    printtable.mean <- cbind(table.mean, starSymbol)
    names(printtable.mean)[8] <- ""
    printtable.mean[["estimate"]] <- as.character(round(table.mean[["estimate"]], digits = digit))
    printtable.mean[["se"]] <- as.character(round(table.mean[["se"]], digits = digit))
    printtable.mean[["df"]] <- as.character(round(table.mean[["df"]], digits = digit))
    printtable.mean[["lower"]] <- as.character(round(table.mean[["lower"]], digits = digit))
    printtable.mean[["upper"]] <- as.character(round(table.mean[["upper"]], digits = digit))
    printtable.mean[["statistic"]] <- as.character(round(table.mean[["statistic"]], digits = digit))
    printtable.mean[["p.value"]] <- format.pval(table.mean[["p.value"]], digits = digit, eps = 10^(-digit))
    if(any(names(printtable.mean) %in% columns == FALSE)){
        printtable.mean <- printtable.mean[,-which(names(printtable.mean) %in% columns == FALSE),drop=FALSE]
    }
    
    if(print && !hide.mean){
        print(printtable.mean)
        cat("\n")
        if(robust){
            cat("Uncertainty was quantified using robust standard errors (column se). \n", sep = "")
        }else{
            cat("Uncertainty was quantified using model-based standard errors (column se). \n", sep = "")
        }
        if(!is.null(object$df)){
            cat("Degrees of freedom were computed using a Satterthwaite approximation (column df). \n", sep = "")
        }
        if("lower" %in% columns && "upper" %in% columns){
            cat("The columns lower and upper indicate a ",100*level,"% confidence interval for each coefficient.\n", sep = "")
        }else if("lower" %in% columns){
            cat("The column lower indicates a ",100*level,"% confidence interval for each coefficient.\n", sep = "")
        }else if("upper" %in% columns){
            cat("The column upper indicate a ",100*level,"% confidence interval for each coefficient.\n", sep = "")
        }
        cat("\n")
    }
    ## ** export
    return(invisible(list(correlation = table.cor,
                          variance = table.var,
                          sd = table.sd,
                          mean = table.mean)))
}


######################################################################
### summary.lmm.R ends here
