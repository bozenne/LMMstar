### table.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (10:48) 
## Version: 
## Last-Updated: jul 26 2024 (17:55) 
##           By: Brice Ozenne
##     Update #: 147
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.tables.effect_lmm
##' @export
model.tables.effect_lmm <- function(x, columns, ...){

    extra.var <- c(x$args$variable,unlist(x$args$time),x$args$strata)
    
    ## ** usual model.tables
    newcolumns <- c("estimate","se","df","lower","upper")
    if(x$args$effect[[1]][1]=="difference"){
        newcolumns <- c(newcolumns,"p.value")
    }

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

    out <- cbind(add, confint(x, ..., columns = newcolumns))
    attr(out, "backtransform") <- NULL
    attr(out, "error") <- NULL
    attr(out, "level") <- NULL
    attr(out, "method") <- NULL
    class(out) <- "data.frame"

    return(out)
}
## * model.tables.lmm
##' @title Statistical Inference and parametrization of Linear Mixed Models
##' @description Export estimated parameters with their uncertainty (standard errors, degrees of freedom, confidence intervals and p-values) from a linear mixed model
##' or a table describing each parameter (type, associated sigma or k parameter, ...).
##'
##' @param x a \code{lmm} object.
##' @param effects [character] Should the CIs/p-values for all coefficients be output (\code{"all"}),
##' or only for mean coefficients (\code{"mean"} or \code{"fixed"}),
##' or only for variance coefficients (\code{"variance"}),
##' or only for correlation coefficients (\code{"correlation"}).
##' Alternatively can be \code{"param"} to output the name and characteristics of each parameter (type, strata, ...).
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param method [character] type of adjustment for multiple comparisons, one of \code{"none"}, \code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}.
##' Alternatively, method for combining the estimates, one of \code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.rubin"}.
##' @param ... arguments to be passed to the \code{confint} method. Should not contain the argument \code{column}.
##' 
##' @details When \code{effects} is not \code{"param"}, this function simply calls \code{\link{confint}} with a specific value for the argument \code{column}.
##' 
##' @keywords methods
##' 
##' @return A \code{data.frame} object.
##' 
##' @export
model.tables.lmm <- function(x, effects = NULL, columns, ...){

    options <- LMMstar.options()

    ## ** normalize user input
    ## *** effects
    if(is.null(effects)){
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }else{
        valid.effects <- c("param","mean","fixed","variance","correlation")
        if(any(effects %in% valid.effects == FALSE)){
            stop("Incorrect value for argument \'effects\'. \n",
                 "Possible values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
        }else if("param" %in% effects & length(effects)>1){
            stop("When argument \'effects\' contains \"param\" it should have length 1. \n")
        }
        effects[effects== "fixed"] <- "mean"
    }

    ## *** columns
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
        
    ## ** extract
    if(identical(effects,"param")){
        out <- x$design$param
    }else{
        out <- confint(x, effects = effects, ..., columns = newcolumns)
        attr(out, "backtransform") <- NULL
        class(out) <- "data.frame"
    }

    ## ** export
    return(out)
}

## * model.tables.mlmm
##' @export
model.tables.mlmm <- function(x, columns, method = NULL, ...){

    options <- LMMstar.options()
    pool.method <- options$pool.method

    if(!is.null(method) && all(method %in% pool.method)){
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
        out$parameter <- rownames(out)
        rownames(out) <- NULL
    }
    return(out)
}

## * model.tables.resample
##' @export
model.tables.resample <- function(x, columns, ...){

    newcolumns <- c("estimate","sample.estimate","se","sample.se","lower","upper","p.value")

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
    class(out) <- "data.frame"
    return(out)
}

## * model.tables.Wald_lmm
##' @title Statistical Inference for Wald test
##' @description Export estimates, standard errors, degrees of freedom, confidence intervals (CIs) and p-values for the mean coefficients of a linear mixed model. 
##'
##' @param x a \code{lmm} object.
##' @param method [character] type of adjustment for multiple comparisons, one of \code{"none"}, \code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"} (see \code{\link{confint.Wald_lmm}} for details).
##' Can also be \code{"contrast"} to output the contrast matrix.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param simplify [logical] 
##' @param ... arguments to be passed to the \code{confint} method. 
##' 
##' @keywords methods
##'
##' @details Argument \bold{simplify}: \itemize{
##' \item \code{method="none"}: attributes (backtransform, error, level, method) are removed from the data.frame.
##' \item \code{method="method"}: the output will be converted into a matrix (instead of a list of matrix) whenever possible.
##' }
##' @return A \code{data.frame} object or (list of) contrast matrices.
##' 
##' @export
model.tables.Wald_lmm <- function(x, method = "none", columns, simplify = TRUE, ...){

    options <- LMMstar.options()
    adj.method <- options$adj.method

    ## ** normalize user input

    ## *** object
    if(x$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** method
    if(!is.character(method) || !is.vector(method)){
        stop("Argument \'method\' must be a character.")
    }
    if(length(method)!=1){
        stop("Argument \'method\' must have length 1.")
    }    
    valid.method <- c("none","contrast",adj.method)
    if(any(method %in% valid.method == FALSE)){
        stop("Unknown value for argument \'type\': \"",paste(setdiff(method, valid.method),collapse = "\", \""),"\". \n",
             "Possible values: \"",paste(valid.method, collapse = "\", \""),"\". \n")
    }

    ## *** columns
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

    ## ** extract from object
    if(method == "contrast"){
        ls.out <- lapply(x$glht, "[[","linfct")
        if(simplify){
            if(x$args$type=="auto"){
                ## combine matrices that are type specific
                out <- as.matrix(Matrix::bdiag(ls.out))
                rownames(out) <- do.call(base::c,lapply(ls.out,rownames))
                colnames(out) <- do.call(base::c,lapply(ls.out,colnames))
            }else{
                ## remove columns with only 0
                out <- ls.out[[1]][,colSums(ls.out[[1]]!=0)>0,drop=FALSE]
            }            
        }else{
            out <- ls.out
        }

    }else{
        out <- confint(x, method = method, columns = newcolumns, ...)
        if(simplify){
            attr(out, "backtransform") <- NULL
            attr(out, "error") <- NULL
            attr(out, "level") <- NULL
            attr(out, "method") <- NULL
        }
        class(out) <- "data.frame"
    }

    ## ** export
    return(out)
}


##----------------------------------------------------------------------
### model.tables.R ends here

