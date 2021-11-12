### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: nov 12 2021 (15:12) 
##           By: Brice Ozenne
##     Update #: 399
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
##' @param add.X additional covariates to be added to the variance and correlation structure.
##' 
##' @keywords internal
##' @examples
##' .formulaStructure(strata ~ time)
##' .formulaStructure( ~ time)
##' .formulaStructure(list( ~ gender+time,  ~ time))
##' .formulaStructure(strata ~ 1)
.formulaStructure <- function(formula, add.X = NULL){

    if(is.list(formula) && length(formula)==2 && all(sapply(formula,inherits,"formula"))){
        init.formula.var <- .formulaStructure(formula[[1]])
        init.formula.cor <- .formulaStructure(formula[[2]])

        if(!identical(init.formula.var$strata,init.formula.cor$strata)){
            stop("Incorrect argument \'formula\': strata variable differ between the correlation and variance structure. \n")
        }
        out <- list(strata = init.formula.var$strata,
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
    }else if(length(var.strata)>1){
        stop("There should be at most one strata variable. \n")
    }

    ## ** right hand side
    var.X <- unique(c(rhs.vars(formula),add.X))
    if(length(var.strata)==0){
        var.strata <- NULL
    }    
    terms.var <- stats::delete.response(stats::terms(formula))
    if(any(attr(terms.var,"order")>1)){
        stop("Does not handle interactions in the formula. \n")
    }

    ## ** combine left and right hand side
    if(length(var.X)==0){
        if(length(var.strata)==0){
            formula.var <- ~1 
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata))
        }
    }else if(length(var.X)==1){
        if(length(var.strata)==0){
            formula.var <- stats::as.formula(paste("~",var.X))
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata,"+",var.X,":",var.strata))
        }
    }else{
        if(length(var.strata)==0){
            formula.var <- stats::as.formula(paste("~",paste(var.X,collapse=":")))
        }else{
            formula.var <- stats::as.formula(paste("~0+",var.strata,"+",paste(paste(var.X,collapse=":"),var.strata,sep=":")))
            ## stats::update(terms.var, paste0("~0+",out$name$strata,"+",out$name$strata,":."))
            ## using ".:var.strata" does not work (it gives the same formula - does not invert . var.strata around the : symbol)
        }
    }
    
    ## ** export
    out <- list(strata = var.strata,
                X.var = unname(var.X),
                X.cor = unname(var.X),
                formula.var = formula.var,
                formula.cor = formula.var
                )
    return(out)
}




## * ID (identity)
##' @title identity Structure
##' @description Variance-covariance structure where the residuals are independent and identically distribution.
##' Can be stratified on a categorical variable.
##' 
##' @param formula formula indicating on which variable to stratify the residual variance (left hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param ... Not used. For compatibility with other structures.
##'
##' @details A typical formula would be \code{~1}.
##'
##' @return An object of class \code{IND} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @examples
##' ID(NULL, var.cluster = "id", var.time = "time")
##' ID(~1, var.cluster = "id", var.time = "time")
##' ID(gender~1, var.cluster = "id", var.time = "time")
##' @export
ID <- function(formula, var.cluster, var.time, ...){

    if(is.null(formula)){
        outCov <- .formulaStructure(~1)
    }else{
        outCov <- .formulaStructure(formula)
    }
    if(length(outCov$X.var)>0 || length(outCov$X.cor)>0){
        stop("Structure \"ID\" cannot handle covariates, consider using structure \"IND\" instead. \n")
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = NA,
                                  cor = NA,
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
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
##' Can be stratified on a categorical variable.
##'
##' @param formula formula indicating variables influencing the residual variance,
##' using either as a multiplicative factor (right hand side) or stratification (left hand side) to model their effect.
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param add.time Should the default formula (i.e. when \code{NULL}) contain a time effect.
##'
##' @details A typical formula would be either \code{~1} indicating constant variance
##' or \code{~time} indicating a time dependent variance.
##' 
##' @return An object of class \code{IND} that can be passed to the argument \code{structure} of the \code{lmm} function.
##'
##' @examples
##' IND(NULL, var.cluster = "id", var.time = "time", add.time = TRUE)
##' IND(~1, var.cluster = "id", var.time = "time")
##' IND(gender~1, var.cluster = "id", var.time = "time")
##' 
##' IND(~time, var.cluster = "id", var.time = "time")
##' IND(gender~time, var.cluster = "id", var.time = "time")
##' IND(~time+gender, var.cluster = "id", var.time = "time")
##' @export
IND <- function(formula, var.cluster, var.time, add.time){
    if(is.null(formula)){
        outCov <- .formulaStructure(~1, add.X = if(add.time && !is.na(var.time)){var.time}else{NULL})
    }else{
        outCov <- .formulaStructure(formula)
    }
    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = NA,
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
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
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance (right hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param ... Not used. For compatibility with other structures.
##'
##' @details A typical formula would be \code{~1}, indicating a variance constant over time and the same correlation between all pairs of times.
##'
##' @return An object of class \code{CS} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' CS(~1, var.cluster = "id", var.time = "time")
##' CS(gender~1, var.cluster = "id", var.time = "time")
##' CS(list(~time,~1), var.cluster = "id", var.time = "time")
##' CS(list(gender~time,gender~1), var.cluster = "id", var.time = "time")
##' 
##' @export
CS <- function(formula, var.cluster, var.time, ...){
    if(is.null(formula)){
        outCov <- .formulaStructure(~1)
    }else{
        outCov <- .formulaStructure(formula)
    }

    if(length(outCov$X.cor)>0){
        stop("Structure \"CS\" cannot handle covariates in the correlation structure. \n")
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = NA,
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                type = "CS")

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
##' @param formula formula indicating on which variable to stratify the residual variance and correlation (left hand side)
##' and variables influencing the residual variance (right hand side).
##' @param var.cluster [character] cluster variable.
##' @param var.time [character] time variable.
##' @param add.time Should the default formula (i.e. when \code{NULL}) contain a time effect.
##'
##' @details A typical formula would be \code{~1}, indicating a time-specific variance parameter and a correlation parameter specific to each pair of times.
##'
##' @return An object of class \code{UN} that can be passed to the argument \code{structure} of the \code{lmm} function.
##' 
##' @examples
##' UN(NULL, var.cluster = "id", var.time = "time", add.time = TRUE)
##' UN(list(~gender,~time), var.cluster = "id", var.time = "time")
##' UN(gender~time, var.cluster = "id", var.time = "time")
##' UN(list(gender~time,gender~time), var.cluster = "id", var.time = "time")
##' 
##' @export
UN <- function(formula, var.cluster, var.time, add.time){
    if(is.null(formula)){
        outCov <- .formulaStructure(~1, add.X = if(add.time){var.time}else{NULL})
    }else{
        outCov <- .formulaStructure(formula)
    }

    out <- list(call = match.call(),
                name = data.frame(cluster = if(!missing(var.cluster)){var.cluster}else{NA},
                                  strata = if(!is.null(outCov$strata)){outCov$strata}else{NA},
                                  time = if(!missing(var.time)){var.time}else{NA},
                                  var = if(length(outCov$X.var)>0){I(list(outCov$X.var))}else{NA},
                                  cor = I(list(outCov$X.cor)),
                                  stringsAsFactors = FALSE),
                formula = list(var = outCov$formula.var,
                               cor = outCov$formula.cor),
                type = "UN")

    ## export
    class(out) <- append("structure",class(out))
    class(out) <- append("UN",class(out))
    return(out)
}

## * EXP (exponential)



##----------------------------------------------------------------------
### structure.R ends here
