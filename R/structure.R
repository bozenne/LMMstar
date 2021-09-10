### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: sep  9 2021 (17:21) 
##           By: Brice Ozenne
##     Update #: 154
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .formulaStructure
##' @examples
##' \dontrun{
##' .formulaStructure(strata ~ time|id)
##' .formulaStructure( ~ time|id)
##' .formulaStructure(list( ~ gender+time|id,  ~ time|id))
##' .formulaStructure(strata ~ 1|id, missing.time.ok = TRUE)
##' .formulaStructure(strata ~ time, missing.id.ok = TRUE)
##' .formulaStructure( ~ time|id)
##' }
.formulaStructure <- function(formula, missing.time.ok = FALSE, missing.id.ok = FALSE){
    
    if(is.list(formula) && length(formula)==2 && all(sapply(formula,inherits,"formula"))){
        init.formula.var <- .formulaStructure(formula[[1]], missing.time.ok = missing.time.ok, missing.id.ok = missing.id.ok)
        init.formula.cor <- .formulaStructure(formula[[2]], missing.time.ok = missing.time.ok, missing.id.ok = missing.id.ok)
        if(!identical(init.formula.cor$cluster,init.formula.var$cluster)){
            stop("Inconsistent cluster variable between the formula.\n")
        }
        if(!identical(init.formula.cor$strata,init.formula.var$strata)){
            stop("Inconsistent strata variable between the formula.\n")
        }
        out <- list(cluster = init.formula.var$cluster,
                    strata = init.formula.var$strata,
                    X.var = init.formula.var$X.var,
                    X.cor = init.formula.cor$X.var,
                    formula.var = init.formula.var$formula.var,
                    formula.cor = init.formula.cor$formula.var) ## note: this is not a mistake that $formula.var (and not $formula.cor)
        return(out)
    }else if(!inherits(formula,"formula")){
        stop("Incorrect argument \'formula\': should be a formula or a list of 2 formula (var, cor).\n")
    }

    ## ** right hand side
    res.split <- strsplit(deparse(formula),"|", fixed = TRUE)[[1]]

    if(length(res.split)>2){
        stop("Incorrect formula for the residual variance-covariance structure. \n",
             "The symbol | should appear at most once, something like: ~ time|id or group ~ time|id. \n")
    }

    if(!grepl("|",deparse(formula),fixed = TRUE)){
        if(missing.id.ok){
            var.cluster <- NULL
            formula.var <- formula
            var.time <- all.vars(stats::update(formula.var,0~.))
        }else{
            stop("Incorrect formula for the residual variance-covariance structure. \n",
                 "No | symbol found so no grouping variable could be defined. \n",
                 "Shoud be something like: ~ time|id or group ~ time|id. \n")
        }
    }else{
        var.cluster <- trimws(res.split[2], which = "both")
        if(length(var.cluster)!=1){
            stop("Incorrect formula for the residual variance-covariance structure. \n",
                 "Should have exactly one variable after the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
        }

        formula.var <- stats::as.formula(res.split[1])
        var.time <- all.vars(stats::update(formula.var,0~.))
        if(length(var.time)>1){
            var.time <- NULL
        }
    }
    
    if(length(var.time)==0){
        if(missing.time.ok){
            var.time <- NULL
        }else{ 
            stop("Incorrect formula for the residual variance-covariance structure. \n",
                 "There should be at least one variable before the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
        }
    }

    ## if(any(attr(terms(formula.var),"order")>1)){
    ##     stop("Does not handle interactions in the formula. \n")
    ## }

    ## ** left hand side
    var.strata <- lhs.vars(formula)
    if(length(var.strata)==0){
        var.strata <- NULL
    }else if(length(lhs.vars(formula))!=1){
        stop("Incorrect formula for the residual variance-covariance structure. \n",
             "There should be at most one variable on the left hand side, something like: ~ time|id or group ~ time|id. \n")
    }

    ## ** export
    out <- list(pattern = NULL,
                cluster = var.cluster,
                strata = var.strata,
                X.var = var.time,
                X.cor = var.time,
                formula.var = formula.var,
                formula.cor = formula.var
                )
    return(out)
}




## * IND (independence)
##' @title Independence Structure
##' @description Variance-covariance structure where the residuals are independent.
##'
##' @param formula formula indicating factors influencing the residual variance.
##'
##' @details A typical formula would be either \code{~1} indicating constant variance
##' or \code{~time} indicating a time dependent variance.
##'
##' @examples
##' IND(~1)
##' IND(~time)
##' 
##' @export
IND <- function(formula,...){
    out0 <- .formulaStructure(formula, missing.time.ok = TRUE, missing.id.ok = TRUE)
    if(!is.null(out0$strata)){
        stop("Cannot stratify when using structure IND. \n")
    }
    out <- list(pattern = NULL,
                cluster = out0$cluster,
                strata = NULL,
                X.var = out0$X.var,
                X.cor = NULL,
                time.var = out0$X.var,
                formula.var = out0$formula.var,
                formula.cor = NULL,
                type = "IND")
    
    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("IND",class(out))
    return(out)
}


## * CS (compound symmetry)
##' @title Compound Symmetry Structure
##' @description Variance-covariance structure where the residuals have constant variance and correlation.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating the cluster and a possible stratification.
##' @param var.cluster [character] used to check the cluster variable in the formula.
##' @param var.time [character] used to check the time variable in the formula.
##'
##' @details A typical formula would be \code{~1|id}, indicating a variance constant over time and the same correlation between all pairs of times.
##'
##' @examples
##' CS(~1|id, var.time = "time", var.cluster = "id")
##' CS(group~1|id, var.time = "time", var.cluster = "id")
##' CS(group~time|id, var.time = "time", var.cluster = "id")
##' 
##' @export
CS <- function(formula, var.cluster, var.time, ...){
    out0 <- .formulaStructure(formula, missing.time.ok = TRUE)

    ## check cluster
    if(!missing(var.cluster) && var.cluster!=out0$cluster){
        stop("Mismatch between cluster variables. \n")
    }

    ## remove time variable that may have been propagated via the argument repetition
    if(!missing(var.time) && !is.null(out0$X.var) && var.time!=out0$X.var){
        stop("Mismatch between time variables for the variance structure. \n")
    }
    if(!missing(var.time) && !is.null(out0$X.cor) && var.time!=out0$X.cor){
        stop("Mismatch between time variables for the correlation structure. \n")
    }
    out <- list(pattern = NULL,
                cluster = out0$cluster,
                strata = out0$strata,
                X.var = NULL,
                X.cor = NULL,
                formula.var = out0$formula.var,
                formula.cor = out0$formula.cor,
                type = "CS")
    if(identical(out0$X.var,out0$X.cor) && length(out0$X.var)==1){
        out$time.var <- out0$X.var
    }

    ##  add strata to the formula (and possibly remove time effect)
    if(is.null(out$strata)){
        out$formula.var <- ~1
        out$formula.cor <- ~1
    }else{
        out$formula.var <- stats::as.formula(paste0("~0+",out$strata))
        out$formula.cor <- stats::as.formula(paste0("~0+",out$strata))
    }
    
    
    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("CS",class(out))
    return(out)
}

## * UN (unstructured)
##' @title Unstructured Structure 
##' @description Variance-covariance structure where the residuals have time-specific variance and correlation.
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating the cluster, factors influencing the variance and the correlation, and a possible stratification.
##' @param var.cluster [character] used to check the cluster variable in the formula.
##'
##' @details A typical formula would be \code{~time} or \code{~time|id}, indicating a time-specific variance parameter and a correlation parameter specific to each pair of times.
##'
##' @examples
##' UN(~time|id, var.cluster = "id")
##' UN(group~time|id, var.cluster = "id")
##' 
##' @export
UN <- function(formula, var.cluster, ...){
    out <- .formulaStructure(formula, missing.time.ok = TRUE)

    ## check cluster
    if(!missing(var.cluster) && var.cluster!=out$cluster){
        stop("Mismatch between cluster variables. \n")
    }
    if(identical(out$X.var,out$X.cor) && length(out$X.var)==1){
        out$time.var <- out$X.var
    }
    out$type <- "UN"

    ##  add strata to the formula
    if(!is.null(out$strata)){

        if(!identical(out$time,all.vars(update(out$formula.var,"0~.")))){
            stop("When using stratification with unstructured covariance pattern, the regressor in the variance formula must the time variable. \n")
        }

        terms.var <- stats::delete.response(stats::terms(out$formula.var))
        ## out$formula.var <- stats::update(terms.var, paste0("~0+",out$strata,"+",out$strata,":."))
        out$formula.var <- stats::update(terms.var, paste0("~0+",out$strata,"+",out$strata,":."))

        terms.cor <- stats::delete.response(stats::terms(out$formula.cor))
        if(any(attr(terms.cor,"order")>1)){
            stop("Does not handle interactions in the correlation formula. \n")
        }
        out$formula.cor <- stats::update(terms.cor, paste0("~0+",out$strata,"+",out$strata,":."))

        ## using ".:var.strata" does not work (it gives the same formula - does not invert . var.strata around the : symbol)
    }
 
    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("UN",class(out))
    return(out)
}

## * EXP (exponential)



##----------------------------------------------------------------------
### structure.R ends here
