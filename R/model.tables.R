### table.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (10:48) 
## Version: 
## Last-Updated: Jul 28 2024 (20:04) 
##           By: Brice Ozenne
##     Update #: 165
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
##' @param method [character] type of adjustment for multiple comparisons, one of \code{"none"}, \code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}.
##' Alternatively, method for combining the estimates, one of \code{"average"}, \code{"pool.se"}, \code{"pool.gls"}, \code{"pool.rubin"}.
##' @param simplify [logical] omit from the output the backtransform attribute.
##' Not relevant when the argument \code{effects="param"}, 
##' @param ... arguments to be passed to \code{\link{confint.lmm}}.
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
        if(simply){
            attr(out, "backtransform") <- NULL
        }
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
##' @param effects [character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' or the contrast matrix (\code{"contrast"})?
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param simplify [logical] with argument \code{effects="Wald"}, omit from the output attributes containing additional information (e.g. approximation error made when adjusting p-values).
##' with argument \code{effects="contrast"} the output will be converted into a matrix (instead of a list of matrix) whenever possible.
##' @param ... arguments to be passed to \code{\link{confint.Wald_lmm}}. 
##' 
##' @keywords methods
##'
##' @export
model.tables.Wald_lmm <- function(x, effects = "Wald", columns, simplify = TRUE, ...){

    ## ** normalize user input

    ## *** object
    if(x$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character.")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1.")
    }
    effects <- match.arg(effects, c("Wald","contrast"))
    
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
    if(effects == "contrast"){
        ls.out <- lapply(x$glht, "[[","linfct")
        if(simplify){
            if(x$args$type=="auto"){
                ## combine matrices that are type specific
                lsType.out <- tapply(names(ls.out), INDEX = sapply(strsplit(names(ls.out), split = "_"),"[[",1), FUN = function(iName){do.call(rbind,ls.out[iName])}, simplify = FALSE)
                out <- as.matrix(Matrix::bdiag(lsType.out))
                rownames(out) <- do.call(base::c,lapply(lsType.out,rownames))
                colnames(out) <- do.call(base::c,lapply(lsType.out,colnames))
            }else{
                ## remove columns with only 0
                out <- ls.out[[1]][,colSums(ls.out[[1]]!=0)>0,drop=FALSE]
            }            
        }else{
            out <- ls.out
        }

    }else if(effects == "Wald"){
        out <- confint(x, columns = newcolumns, ...)
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

