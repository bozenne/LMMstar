### table.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (10:48) 
## Version: 
## Last-Updated: Oct 10 2022 (12:23) 
##           By: Brice Ozenne
##     Update #: 30
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.tables.lmm (documentation)
##' @title Statistical Inference for Linear Mixed Model
##' @description Export estimates, standard errors, degrees of freedom, confidence intervals (CIs) and p-values for the mean coefficients of a linear mixed model. 
##'
##' @param x a \code{lmm} object.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}, \code{"partial.r"}.
##' @param method [character] type of adjustment for multiple comparisons, one of \code{"none"}, \code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}.
##' Alternatively, method for combining the estimates, one of \code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.rubin"}.
##' @param ... arguments to be passed to the \code{confint} method. Should not contain the argument \code{column}.
##' 
##' @details This function simply calls \code{\link{confint}} with a specific value for the argument \code{column}.
##' 
##' @export
model.tables.lmm <- function(x, columns, ...){

    newcolumns <- c("estimate","se","df","lower","upper","p.value")

    if(!missing(columns)){
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            newcolumns <- union(newcolumns, unname(columns))
        }else if(!is.null(names(columns)) && all(names(columns)=="remove")){
            newcolumns <- setdiff(newcolumns, unname(columns))
        }else{
            newcolumns <- columns
        }
    }

    out <- confint(x, ..., columns = newcolumns)
    attr(out, "backtransform") <- NULL
    class(out) <- "data.frame"
    return(out)
}

## * model.tables.Wald_lmm (documentation)
##' @export
model.tables.Wald_lmm <- function(x, columns, ...){

    newcolumns <- c("estimate","se","df","lower","upper","p.value")

    if(!missing(columns)){
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            newcolumns <- union(newcolumns, unname(columns))
        }else if(!is.null(names(columns)) && all(names(columns)=="remove")){
            newcolumns <- setdiff(newcolumns, unname(columns))
        }else{
            newcolumns <- columns
        }
    }

    out <- confint(x, ..., columns = newcolumns)
    attr(out, "backtransform") <- NULL
    attr(out, "error") <- NULL
    attr(out, "level") <- NULL
    attr(out, "method") <- NULL
    class(out) <- "data.frame"
    return(out)
}

## * model.tables.mlmm (documentation)
##' @export
model.tables.mlmm <- function(x, columns, method = NULL, ...){

    if(!is.null(method) && method %in% c("average","pool.se","pool.gls","pool.rubin")){
        newcolumns <- c("estimate","se","df","lower","upper","p.value")
    }else{
        newcolumns <- c("by","parameter","estimate","se","df","lower","upper","p.value")
    }

    if(!missing(columns)){
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            newcolumns <- union(newcolumns, unname(columns))
        }else if(!is.null(names(columns)) && all(names(columns)=="remove")){
            newcolumns <- setdiff(newcolumns, unname(columns))
        }else{
            newcolumns <- columns
        }
    }

    out <- confint(x, method = method, ..., columns = newcolumns)
    attr(out, "backtransform") <- NULL
    class(out) <- "data.frame"
    rownames(out) <- NULL
    return(out)
}

##----------------------------------------------------------------------
### model.tables.R ends here
