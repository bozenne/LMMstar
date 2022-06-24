### summary.mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2022 (14:45) 
## Version: 
## Last-Updated: jun 24 2022 (18:01) 
##           By: Brice Ozenne
##     Update #: 46
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.mlmm
##' @title Summary of Multiple Linear Mixed Models
##' @description Estimates, p-values, and confidence intevals for multiple linear mixed models.
##' 
##' @param object an \code{mlmm} object, output of \code{mlmm}.
##' @param method [character] type of adjustment for multiple comparisons: one of \code{"none"}, \code{"bonferroni"}, \code{"single-step"}, \code{"single-step2"}.
##' @param print [logical] should the output be printed in the console.
##' Can be a vector of length 2 where the first element refer to the global tests and the second to the individual tests.
##' @param hide.data [logical] should information about the dataset not be printed.
##' @param hide.fit [logical] should information about the model fit not be printed.
##' @param ... other arguments are passed to \code{\link{summary.anova_lmm}}.
##'
##' @export
summary.mlmm <- function(object, method = NULL, print = NULL, hide.data = FALSE, hide.fit = FALSE, ...){

    if(is.null(method)){
        method <- "none"
    }
    if(is.null(print)){
        print <- c(FALSE,TRUE)
    }


    ## extract models
    ls.model <- estimate(object)
    method.fit <- ls.model[[1]]$method.fit
    optimizer <- ls.model[[1]]$opt$name
    logLik <- sapply(ls.model, logLik)
    cv <- sapply(ls.model, function(iM){iM$opt$cv})
    n.iter <- sapply(ls.model, function(iM){iM$opt$n.iter})
    param.value <- coef(object, effects = "all")

    nparam.mu  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="mu")})
    nparam.sigma  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="sigma")})
    nparam.k  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="k")})
    nparam.rho  <- sapply(ls.model, function(iM){sum(iM$design$param$type=="rho")})

    M.nobs <- do.call(rbind,lapply(ls.model, stats::nobs))
    vec.nobs <- apply(M.nobs,2,paste,collapse = ", ")
    nobsByCluster <- lapply(ls.model,function(iM){iM$design$cluster$nobs})
    n.cluster.original <- sapply(ls.model,function(iM){iM$design$cluster$n})
    n.cluster.design <- sapply(ls.model,function(iM){iM$design$cluster$n})

    ## ** welcome message
    if(any(print)){
        cat("	Linear Mixed Models stratified according to \"",object$call$by,"\" \n\n",sep="")
    }

    ## ** data message    
    if(!hide.data){
        cat("Dataset:", deparse(object$call$data), "\n")
        cat("Strata : \"",paste(names(ls.model),collapse = "\", \""),"\"\n\n",sep="")
        if(any(M.nobs[,"missing"]>0)){
            if(any(n.cluster.original-n.cluster.design>0)){
                cat("  - ", vec.nobs["cluster"], " clusters were analyzed, ",paste(n.cluster.original-n.cluster.design,collapse=", ")," were excluded because of missing values \n" , sep = "")
            }else{
                cat("  - ", vec.nobs["cluster"], " clusters \n" , sep = "")
            }
            cat("  - ", paste(sapply(nobsByCluster,sum), collapse = ", "), " observations were analyzed, ",vec.nobs["missing"]," were excluded because of missing values \n",  sep = "")
        }else{
            cat("  - ", vec.nobs["cluster"], " clusters \n" , sep = "")
            cat("  - ", paste(sapply(nobsByCluster,sum), collapse = ", "), " observations \n",  sep = "")
        }
        cat("\n")
    }
    
    ## ** optim message    
    if(!hide.fit){
        cat("Estimation procedure \n\n")
        if(method.fit == "REML"){
            cat("  - Restricted Maximum Likelihood (REML) \n")
        }else{
            cat("  - Maximum Likelihood (ML) \n")
        }
        cat("  - log-likelihood :", paste(as.double(logLik),collapse = ", "), "\n",sep="")
        cat("  - parameters: mean = ",paste(nparam.mu,collapse = ", "),", variance = ",paste(nparam.sigma+nparam.k,collapse = ", "),", correlation = ",paste(nparam.rho,collapse = ", "),"\n", sep = "")
        if(optimizer!="gls"){
            cat("  - convergence: ",paste(cv>0,collapse = ", ")," (",paste(n.iter,collapse = ", ")," iterations) \n", sep = "")
        }
        cat(" \n")
    }


    ## ** extract test
    if(any(print)){
        cat("Statistical inference \n")
        options <- LMMstar.options()
        out <- summary.anova_lmm(object, method = method, columns = options$columns.summary, print = print/2, ...)
    }

    ## ** export
    return(invisible(out))
}

##----------------------------------------------------------------------
### summary.mlmm.R ends here
