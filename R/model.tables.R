### table.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (10:48) 
## Version: 
## Last-Updated: maj 10 2024 (16:12) 
##           By: Brice Ozenne
##     Update #: 62
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
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param method [character] type of adjustment for multiple comparisons, one of \code{"none"}, \code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}.
##' Alternatively, method for combining the estimates, one of \code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.rubin"}.
##' @param ... arguments to be passed to the \code{confint} method. Should not contain the argument \code{column}.
##' 
##' @details This function simply calls \code{\link{confint}} with a specific value for the argument \code{column}.
##' 
##' @keywords methods
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

    if(!is.null(method) && all(method %in% c("average","pool.fixse","pool.se","pool.gls","pool.gls1","pool.rubin"))){
        newcolumns <- c("estimate","se","df","lower","upper","p.value")
        rm.rownames <- FALSE
    }else{
        newcolumns <- c("by","parameter","estimate","se","df","lower","upper","p.value")
        rm.rownames <- TRUE
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
    if(rm.rownames){
        rownames(out) <- NULL
    }
    return(out)
}

## * model.tables.effect_lmm
##' @export
model.tables.effect_lmm <- function(x, columns, ...){

    extra.var <- c(x$args$variable,unlist(x$args$time),x$args$strata)
    
    ## ** usual model.tables
    newcolumns <- c("estimate","se","df","lower","upper","p.value")

    if(!missing(columns)){
        
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            newcolumns <- union(c(newcolumns, extra.var), unname(columns))
        }else if(!is.null(names(columns)) && all(names(columns)=="remove")){
            newcolumns <- setdiff(c(newcolumns, extra.var), unname(columns))
        }else{
            newcolumns <- c(newcolumns, extra.var)
        }
        if(any(newcolumns %in% extra.var)){
            add <- x$univariate[intersect(newcolumns,extra.var)]
            newcolumns <- setdiff(newcolumns,extra.var)
        }else{
            add <- NULL
        }
    }else{
        add <- x$univariate[extra.var]
    }

    out <- cbind(add,
                 confint(x, ..., columns = newcolumns))
    rownames(out) <- NULL
    attr(out, "backtransform") <- NULL
    attr(out, "error") <- NULL
    attr(out, "level") <- NULL
    attr(out, "method") <- NULL
    class(out) <- "data.frame"

    return(out)
}
##----------------------------------------------------------------------
### model.tables.R ends here
