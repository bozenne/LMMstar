### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: okt  1 2021 (17:08) 
##           By: Brice Ozenne
##     Update #: 293
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .formulaStructure
##' @title Extract Variable From Formula For VCOV Structure
##' @description Extract the variables from the variance and correlation formula to be used to initialize the variance-covariance structure.
##' @noRd
##'
##' @param formula A formula or a list of two formulas.
##' @param missing.time.ok [logical] If FALSE an error is triggered when the function could not identify the time variable. 
##' @param missing.id.ok [logical] If FALSE an error is triggered when the function could not identify the cluster variable.
##' 
##' @keywords internal
##' @examples
##' .formulaStructure(strata ~ time|id)
##' .formulaStructure( ~ time|id)
##' .formulaStructure(list( ~ gender+time|id,  ~ time|id))
##' .formulaStructure(strata ~ 1|id, missing.time.ok = TRUE)
##' .formulaStructure(strata ~ time, missing.id.ok = TRUE)
##' .formulaStructure( ~ time|id)
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
                    time.var = if(identical(init.formula.var$time.var,init.formula.cor$time.var)){init.formula.var$time.var}else{NULL},
                    X.var = init.formula.var$X.var,
                    X.cor = init.formula.cor$X.var,
                    formula.var = init.formula.var$formula.var,
                    formula.cor = init.formula.cor$formula.var) ## note: this is not a mistake that $formula.var (and not $formula.cor)
        return(out)
    }else if(!inherits(formula,"formula")){
        stop("Incorrect argument \'formula\': should be a formula or a list of 2 formula (var, cor).\n")
    }

    ## ** left hand side
    var.strata <- lhs.vars(formula)
    if(length(var.strata)==0){
        var.strata <- NULL
    }else if(length(var.strata)!=1){
        stop("Incorrect formula for the residual variance-covariance structure. \n",
             "There should be at most one variable on the left hand side, something like: ~ time|id or group ~ time|id. \n")
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
            var.rhs <- rhs.vars(formula)
            if(length(var.rhs)==0){
                var.time <- NULL
                formula.var <- ~1 ## no need to stratify on anything
            }else if(length(var.rhs)==1){
                var.time <- var.rhs
                if(length(var.strata)==0){
                    formula.var <- stats::as.formula(paste("~",var.rhs))
                }else{
                    formula.var <- stats::as.formula(paste(var.strata,"~",var.rhs))
                }
            }else{
                var.time <- NULL
                if(length(var.strata)==0){
                    formula.var <- stats::as.formula(paste("~",paste(var.rhs,collapse=":")))
                }else{
                    formula.var <- stats::as.formula(paste(var.strata,"~",paste(var.rhs,collapse=":")))
                }
            }
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
    }
    
    formula.var <- stats::as.formula(res.split[1])
    terms.var <- stats::terms(formula.var)
    if(any(attr(terms.var,"order")>1)){
        stop("Does not handle interactions in the formula. \n")
    }
    var.rhs <- rhs.vars(formula.var)
    if(length(var.rhs)==0){
        var.time <- NULL
        formula.var <- ~1 ## no need to stratify on anything
    }else if(length(var.rhs)==1){
        var.time <- var.rhs
        if(length(var.strata)==0){
            formula.var <- stats::as.formula(paste("~",var.rhs))
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata,"+",var.rhs,":",var.strata))
        }
    }else{
        var.time <- NULL
        if(length(var.strata)==0){
            formula.var <- stats::as.formula(paste("~",paste(var.rhs,collapse=":")))
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata,"+",paste(paste(var.rhs,collapse=":"),var.strata,sep=":")))
            ## stats::update(terms.var, paste0("~0+",out$name$strata,"+",out$name$strata,":."))
            ## using ".:var.strata" does not work (it gives the same formula - does not invert . var.strata around the : symbol)
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

    ## ** export
    out <- list(pattern = NULL,
                cluster = var.cluster,
                strata = var.strata,
                var.time = var.time,
                X.var = all.vars(formula.var),
                X.cor = all.vars(formula.var),
                formula.var = formula.var,
                formula.cor = formula.var
                )
    return(out)
}




## * ID (identity)
##' @title identity Structure
##' @description Variance-covariance structure where the residuals are independent and identically distribution.
##'
##' @param formula formula indicating the time and cluster variables.
##' @param var.time [character] name of the time variable.
##' @param ... not used.
##'
##' @details A typical formula would be either \code{~1}.
##'
##' @return An object of class \code{IND} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @examples
##' ID(~1)
##' ID(~time)
##' ID(~time+gender)
##' ID(~time+gender,var.time="time")
##' ID(gender~time,var.time="time")
##' 
##' @export
ID <- function(formula, var.time, ...){
    out0 <- .formulaStructure(formula, missing.time.ok = TRUE, missing.id.ok = TRUE)
    out <- list(call = match.call(),
                name = data.frame(cluster = if(length(out0$cluster)==1){out0$cluster}else{NA},
                                  strata = if(!is.null(out0$strata)){out0$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else if(length(out0$var.time)==1){out0$var.time}else{NA},
                                  var = NA,
                                  cor = NA),
                formula = list(var = ~1,
                               cor = NULL),
                type = "IND")


    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("IND",class(out))
    return(out)
}

## * IND (independence)
##' @title Independence Structure
##' @description Variance-covariance structure where the residuals are independent.
##'
##' @param formula formula indicating factors influencing the residual variance.
##' @param var.time [character] name of the time variable.
##' @param ... not used.
##'
##' @details A typical formula would be either \code{~1} indicating constant variance
##' or \code{~time} indicating a time dependent variance.
##' 
##' @return An object of class \code{IND} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @examples
##' IND(~1)
##' IND(~time)
##' IND(~time+gender)
##' IND(~time+gender,var.time="time")
##' IND(gender~time,var.time="time")
##' 
##' @export
IND <- function(formula, var.time, ...){
    out0 <- .formulaStructure(formula, missing.time.ok = TRUE, missing.id.ok = TRUE)
    out <- list(call = match.call(),
                name = data.frame(cluster = if(length(out0$cluster)==1){out0$cluster}else{NA},
                                  strata = if(!is.null(out0$strata)){out0$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else if(length(out0$var.time)==1){out0$var.time}else{NA},
                                  var = if(length(out0$X.var)>0){I(list(out0$X.var))}else{NA},
                                  cor = NA),
                formula = list(var = out0$formula.var,
                               cor = NULL),
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
##' @param ... not used.
##'
##' @details A typical formula would be \code{~1|id}, indicating a variance constant over time and the same correlation between all pairs of times.
##'
##' @return An object of class \code{CS} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' CS(~1|id)
##' CS(~1|id, var.time = "time", var.cluster = "id")
##' CS(group~1|id)
##' CS(group~time|id, var.time = "time", var.cluster = "id")
##' 
##' @export
CS <- function(formula, var.cluster, var.time, ...){
    out0 <- .formulaStructure(formula, missing.time.ok = TRUE)

    ## check cluster
    if(!missing(var.cluster) && var.cluster!=out0$cluster){
        stop("Mismatch between cluster variables. \n")
    }
    if(!missing(var.time) && !is.null(out0$time.var) && var.time!=out0$time.var){
        stop("Mismatch between time variables for the variance structure. \n")
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = out0$cluster,
                                  strata = if(!is.null(out0$strata)){out0$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else if(length(out0$var.time)==1){out0$var.time}else{NA},
                                  var = NA,
                                  cor = NA),
                formula = list(var = NULL,
                               cor = NULL),
                type = "CS")

    ## check variable in formula
    if(length(setdiff(out$X.var,c(out$var.time,out$strata)))>1){
        stop("Should be at no covariate in the variance formula, except thoses indicating time and strata. \n")
    }
    if(length(setdiff(out$X.cor,c(out$var.time,out$strata)))>1){
        stop("Should be at no covariate in the correlation formula, except thoses indicating time and strata. \n")
    }
    ##  add strata to the formula (and possibly remove time effect)
    if(is.na(out$name$strata)){
        out$formula$var <- ~1
        out$formula$cor <- ~1
    }else{
        out$name$var <- I(list(out0$strata))
        out$name$cor <- I(list(out0$strata))
        out$formula$var <- stats::as.formula(paste0("~0+",out0$strata))
        out$formula$cor <- stats::as.formula(paste0("~0+",out0$strata))
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
##' @param var.time [character] used to check the time variable in the formula.
##' @param ... not used.
##'
##' @details A typical formula would be \code{~time} or \code{~time|id}, indicating a time-specific variance parameter and a correlation parameter specific to each pair of times.
##'
##' @return An object of class \code{UN} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' UN(~time|id)
##' UN(~time+gender|id)
##' UN(group~time|id, var.cluster = "id")
##' UN(group~time|id, var.cluster = "id", var.time = "time")
##' 
##' @export
UN <- function(formula, var.cluster, var.time, ...){
    out0 <- .formulaStructure(formula, missing.time.ok = TRUE)

    ## check cluster
    if(!missing(var.cluster) && var.cluster!=out0$cluster){
        stop("Mismatch between cluster variables. \n")
    }
    if(!missing(var.time) && !is.null(out0$time.var) && var.time!=out0$time.var){
        stop("Mismatch between time variables for the variance structure. \n")
    }
    out <- list(call = match.call(),
                name = data.frame(cluster = out0$cluster,
                                  strata = if(length(out0$strata)==1){out0$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else if(length(out0$var.time)==1){out0$var.time}else{NA},
                                  var = I(list(out0$X.var)),
                                  cor = I(list(out0$X.cor))),
                formula = list(var = out0$formula.var,
                               cor = out0$formula.cor),
                type = "UN")

    ##  check covariates when stratifying
    if(!is.na(out$name$strata)){

        if(!identical(out$name$time,setdiff(out$name$var[[1]],out$name$strata))){
            stop("When using stratification with unstructured covariance pattern, the regressor in the variance formula must the time variable. \n")
        }
        if(!identical(out$name$time,setdiff(out$name$cor[[1]],out$name$strata))){
            stop("When using stratification with unstructured covariance pattern, the regressor in the correlation formula must the time variable. \n")
        }

    }
 
    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("UN",class(out))
    return(out)
}

## * EXP (exponential)



##----------------------------------------------------------------------
### structure.R ends here
