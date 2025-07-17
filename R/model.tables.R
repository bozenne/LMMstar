### table.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (10:48) 
## Version: 
## Last-Updated: jul 17 2025 (10:16) 
##           By: Brice Ozenne
##     Update #: 398
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.tables.lmm
##' @title Statistical Inference and parametrization of a Linear Mixed Model
##' @description Export estimated parameters with their uncertainty (standard errors, degrees-of-freedom, confidence intervals and p-values) from a linear mixed model
##' or a table describing each parameter (type, associated sigma or k parameter, ...).
##'
##' @param x a \code{lmm} object.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param effects [character] Should the CIs/p-values for all coefficients be output (\code{"all"}),
##' or only for mean coefficients (\code{"mean"} or \code{"fixed"}),
##' or only for variance coefficients (\code{"variance"}),
##' or only for correlation coefficients (\code{"correlation"}).
##' Alternatively can be \code{"param"} to output the name and characteristics of each parameter (type, strata, ...).
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors.
##' Can also be \code{2} compute the degrees-of-freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
##' @param null [numeric vector] the value of the null hypothesis relative to each coefficient.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficient. Otherwise a normal distribution is used.
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"statistic"}, \code{"df"}, \code{"null"}, \code{"lower"}, \code{"upper"}, \code{"p.value"}.
##' @param type.information,transform.sigma,transform.k,transform.rho,transform.names are passed to the \code{vcov} method. See details section in \code{\link{coef.lmm}}.
##' @param backtransform [logical] should the variance/covariance/correlation coefficient be backtransformed?
##' @param simplify [logical] omit from the output the backtransform attribute.
##' Not relevant when the argument \code{effects="param"}, 
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @details When \code{effects} differs from \code{"param"}, this function is a wrapper for \code{\link{confint}} with different default value for the argument \code{column}.
##' 
##' @keywords methods
##' 
##' @return A \code{data.frame} object.
##' 
##' @export
model.tables.lmm <- function(x, level = 0.95, effects = NULL, robust = FALSE, null = NULL,
                             columns = NULL,
                             df = NULL, type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE,
                             backtransform = NULL, simplify = TRUE, ...){

    ## ** normalize user input
    ## *** dots
    ## will be initialized in confint (e.g. to keep the options)
    
    ## *** effects
    if(is.null(effects)){
        ## do nothing, let confint initialize
    }else{ 
        if(!is.character(effects) || !is.vector(effects)){
            stop("Argument \'effects\' must be a character vector. \n")
        }
        valid.effects <- c("param","mean","fixed","variance","correlation","all")
        if(any(effects %in% valid.effects == FALSE)){
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
        }
        ## otherwise let anova initialize, e.g. all or fixed
    }

    ## *** columns
    newcolumns <- c("estimate","se","df","lower","upper","p.value")
    if(!is.null(columns)){
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            newcolumns <- union(newcolumns, unname(columns))
        }else if(!is.null(names(columns)) && all(names(columns)=="remove")){
            newcolumns <- setdiff(newcolumns, unname(columns))
        }else{
            newcolumns <- columns
        }
    }

    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1 or FALSE/0. \n")
    }

    ## ** extract
    if(!is.null(effects) && "param" %in% effects){

        out <- cbind(trans.name = x$design$param$name, x$design$param)        

        ## *** add transformed names
        if(is.null(transform.sigma) && is.null(transform.k) && is.null(transform.rho)){
            if(!is.null(x$reparametrize$newname) && transform.names){
                out[match(names(x$reparametrize$p),out$name),"trans.name"] <- x$reparametrize$newname
            }else{
                ## do nothing as name should not be changed or no reparametrisation was performed
            }
        }else{

            init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                                    x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho,
                                    simplify = FALSE)

            if((init$transform.sigma=="none") & (init$transform.k=="none" | "k" %in% out$type == FALSE) & (init$transform.rho=="none" | "rho" %in% out$type == FALSE) | transform.names == FALSE){
                ## do nothing as there is no transformation
            }else{
                index.reparametrize <- match(names(x$reparametrize$p),out$name)
                out[index.reparametrize,"trans.name"] <- .reparametrize(p = rep(NA, length(index.reparametrize)),  
                                                                        type = out[index.reparametrize,"type"],
                                                                        sigma = out[index.reparametrize,"sigma"],
                                                                        k.x = out[index.reparametrize,"k.x"],
                                                                        k.y = out[index.reparametrize,"k.y"],
                                                                        level = out[index.reparametrize,"level"],                                              
                                                                        Jacobian = FALSE, dJacobian = FALSE, inverse = FALSE,
                                                                        transform.sigma = init$transform.sigma,
                                                                        transform.k = init$transform.k,
                                                                        transform.rho = init$transform.rho,
                                                                        transform.names = TRUE)$newname
            }
        }

        ## *** subset
        effects2 <- setdiff(effects,"param")
        if(length(effects2)>0){
            if(all("all" %in% effects2)){
                effects2 <- c("mean","variance","correlation")
            }else{
                effects2[effects2== "fixed"] <- "mean"
            }
            keep.type <- unlist(lapply(effects2, switch,
                                       "mean" = "mu",
                                       "variance" = c("sigma","k"),
                                       "correlation" = "rho"))
            out <- out[out$type %in% keep.type,,drop=FALSE]
        }
    }else{
        out <- stats::confint(x, level = 0.95, effects = effects, robust = robust, null = null,
                              columns = newcolumns,
                              df = df, type.information = type.information, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names,
                              backtransform = backtransform, ...)
        if(simplify){
            attr(out, "backtransform") <- NULL
        }
        class(out) <- "data.frame"
    }

    ## ** export
    return(out)
}

## * model.tables.Wald_lmm
##' @title Statistical Inference for Wald tests
##' @description Export estimates, standard errors, degrees-of-freedom, confidence intervals (CIs) and p-values
##' relative to linear contrasts involved in Wald tests.
##'
##' @param x a \code{Wald_lmm} object.
##' @param effects [character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' the contrast matrix (\code{"contrast"}),
##' or the name/value/type of the underlying mixed model parameters (\code{"param"})?
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the Wald statistic. Otherwise a normal distribution is used.
##' @param method [character] Should pointwise confidence intervals be output (\code{"none"}) or simultaneous confidence intervals (\code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"})? 
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"df"}, \code{"quantile"}, \code{"lower"}, \code{"upper"}, \code{"statistic"}, \code{"null"}, \code{"p.value"}.
##' @param backtransform [logical] should the estimates, standard errors, and confidence intervals be backtransformed?
##' @param transform.names [logical] should the name of the coefficients be updated to reflect the transformation that has been used?
##' Only relevant when \code{effects="contrast"}.
##' @param simplify [logical] with argument \code{effects="Wald"}, omit from the output attributes containing additional information (e.g. approximation error made when adjusting p-values).
##' with argument \code{effects="contrast"} the output will be converted into a matrix (instead of a list of matrix) whenever possible.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @details When \code{effects} is \code{"Wald"}, this function is a wrapper for \code{\link{confint}} with different default value for the argument \code{column}.
##' 
##' @keywords methods
##'
##' @export
model.tables.Wald_lmm <- function(x, effects = "Wald",
                                  level = 0.95, df = NULL, method = NULL, columns = NULL, backtransform = NULL, transform.names = TRUE,
                                  simplify = TRUE, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- NULL
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** object
    if(x$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character. \n")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1. \n")
    }
    valid.effects <- c("Wald","contrast","param")
    if(effects %in% valid.effects == FALSE){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }

    ## *** columns
    if(effects == "Wald"){
        newcolumns <- c("estimate","se","df","lower","upper","p.value")
        if(!is.null(columns)){
            if(!is.null(names(columns)) && all(names(columns)=="add")){
                newcolumns <- union(newcolumns, unname(columns))
            }else if(!is.null(names(columns)) && all(names(columns)=="remove")){
                newcolumns <- setdiff(newcolumns, unname(columns))
            }else{
                newcolumns <- columns
            }
        }
    }

    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1 or FALSE/0. \n")
    }

    ## ** extract from object
    if(effects == "param"){

        out <- x$param

    }else if(effects == "contrast"){
        table.param <- stats::model.tables(x, effects = "param")
        ls.out <- lapply(x$glht, function(iGlht){ ## iGlht <- x$glht[[1]]
            if(simplify){
                iOut <- iGlht$linfct
            }else{
                iOut <- matrix(0, nrow = NROW(iGlht$linfct), ncol = NROW(table.param),
                               dimnames = list(rownames(iGlht$linfct), table.param$name))
                iOut[,colnames(iGlht$linfct)] <- iGlht$linfct
            }
            if(transform.names){
                colnames(iOut) <- table.param[match(colnames(iOut), table.param$name),"trans.name"]
            }
            return(iOut)
        })

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
        out <- stats::confint(x, level = level, df = df, method = method, columns = newcolumns, backtransform = backtransform, options = options)
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

## * model.tables.rbindWald_lmm
##' @title Statistical Inference From Combined Wald Tests
##' @description Combine estimates, standard errors, degrees-of-freedom, confidence intervals (CIs) and p-values
##' relative to linear contrasts of parameters from different linear mixed models. 
##'
##' @param x a \code{mlmm} object.
##' @param effects [character] should the linear contrasts involved in the Wald test be output (\code{"Wald"}),
##' the contrast matrix (\code{"contrast"}),
##' or the name/value/type of the underlying mixed model parameters (\code{"param"})?
##' @param level [numeric, 0-1] nominal coverage of the confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the Wald statistic. Otherwise a normal distribution is used.
##' @param method [character] Should pointwise confidence intervals be output (\code{"none"}) or simultaneous confidence intervals (\code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"})? 
##' @param columns [character vector] Columns to be output.
##' Can be any of \code{"estimate"}, \code{"se"}, \code{"df"}, \code{"quantile"}, \code{"lower"}, \code{"upper"}, \code{"statistic"}, \code{"null"}, \code{"p.value"}.
##' @param ordering [character] should the output be ordered by name of the linear contrast (\code{"contrast"}) or by model (\code{"model"}).
##' @param backtransform [logical] should the estimates, standard errors, and confidence intervals be backtransformed?
##' @param transform.names [logical] should the name of the coefficients be updated to reflect the transformation that has been used?
##' Only relevant when \code{effects="contrast"}.
##' @param simplify [logical] with argument \code{effects="Wald"}, omit from the output attributes containing additional information (e.g. approximation error made when adjusting p-values).
##' with argument \code{effects="contrast"} the output will be converted into a matrix (instead of a list of matrix) whenever possible.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @details When \code{effects} is \code{"Wald"}, this function simply calls \code{\link{confint}} with a specific value for the argument \code{column}.
##' 
##' @keywords methods
##' @return A \code{data.frame} object.
##'
##' @export
model.tables.rbindWald_lmm <- function(x, effects = "Wald",
                                       level = 0.95, df = NULL, method = NULL, columns = NULL, ordering = NULL, backtransform = NULL,
                                       transform.names = TRUE, simplify = TRUE, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- NULL
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** object
    if(x$args$univariate == FALSE){
        message("Nothing to return: consider setting argument \'univariate\' to TRUE when calling rbind.Wald_lmm. \n")
        return(invisible(NULL))
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character. \n")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1. \n")
    }
    valid.effects <- c("Wald","contrast","param")
    if(effects %in% valid.effects == FALSE){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }

    ## *** columns
    if(effects == "Wald"){
        newcolumns <- c("estimate","se","df","lower","upper","p.value")
        if(!is.null(columns)){
            if(!is.null(names(columns)) && all(names(columns)=="add")){
                newcolumns <- union(newcolumns, unname(columns))
            }else if(!is.null(names(columns)) && all(names(columns)=="remove")){
                newcolumns <- setdiff(newcolumns, unname(columns))
            }else{
                newcolumns <- columns
            }
        }
    }

    ## *** ordering
    table.param <- x$param
    if(!is.null(ordering)){
        ordering <- match.arg(ordering, c("contrast","model"))
        ordering.var <- switch(ordering, "model" = "model", "contrast" = "name")
        ## make sure that (Intercept) weight sigma is not ordered (Intercept) sigma weight 
        table.param.order <- table.param[order(factor(table.param[[ordering.var]], levels = unique(table.param[[ordering.var]]))),,drop=FALSE]
    }
    
    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1 or FALSE/0. \n")
    }

    ## ** extract from object
    if(effects == "param"){

        if(is.null(ordering)){
            out <- table.param
        }else{
            out <- table.param.order
        }

        
    }else if(effects == "contrast"){
        glht.linfct <- x$glht[[1]]$linfct ## by design single glht when rbindWald object

        if(simplify){ ## remove columns with only 0
            if(transform.names){
                colnames(glht.linfct) <- table.param[match(colnames(glht.linfct), table.param$Uname),"trans.Uname"]
            }
            out <- glht.linfct[,colSums(glht.linfct!=0)>0,drop=FALSE]
        }else{
            out <- matrix(0, nrow =  NROW(glht.linfct), ncol = NROW(table.param),
                          dimnames = list(rownames(glht.linfct),table.param$Uname))
            out[,colnames(glht.linfct)] <- glht.linfct
            if(transform.names){
                colnames(out) <- table.param$trans.Uname
            }
        }

        if(!is.null(ordering)){
            ordering.var2 <- switch(ordering, "model" = "model", "contrast" = "term")
            ordering.row <- rownames(x$univariate)[order(x$univariate[[ordering.var2]])]
            if(transform.names){                
                out <- out[ordering.row,intersect(table.param.order$trans.Uname,colnames(out)),drop=FALSE]
            }else{
                out <- out[ordering.row,intersect(table.param.order$Uname,colnames(out))]
            }
        }

        if(!simplify){
            out <- list(user = out)
        }


    }else if(effects == "Wald"){

        out <- stats::confint(x, level = level, df = df, method = method, columns = newcolumns, ordering = ordering, backtransform = backtransform, options = options)
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

## * model.tables.mlmm
##' @title Statistical Inference and parametrization of Multiple Linear Mixed Model
##' @description Combine estimated parameters with their uncertainty (standard errors, degrees-of-freedom, confidence intervals and p-values) from group-specific linear mixed models
##' or a table describing each parameter (type, associated sigma or k parameter, ...).
##'
##' @param x a \code{mlmm} object.
##' @inheritParams model.tables.rbindWald_lmm
##' 
##' @keywords methods
##' @return A \code{data.frame} object.
##' 
##' @export
model.tables.mlmm <- model.tables.rbindWald_lmm

## * model.tables.effect_lmm
##' @export
model.tables.effect_lmm <- function(x, effects = "Wald", columns, ...){

    ## ** normalize user input
    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character. \n")
    }
    if(length(effects)!=1){
        stop("Argument \'effects\' must have length 1. \n")
    }
    valid.effects <- c("Wald","contrast","param")
    if(effects %in% valid.effects == FALSE){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }
    if(effects %in% c("contrast","param")){
        return(model.tables.Wald_lmm(x, ...))
    }
    
    ## ** usual model.tables
    extra.var <- c(x$args$variable,unlist(x$args$time),x$args$strata)

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

    out <- cbind(add, stats::confint(x, ..., columns = newcolumns))
    attr(out, "backtransform") <- NULL
    attr(out, "error") <- NULL
    attr(out, "level") <- NULL
    attr(out, "method") <- NULL
    class(out) <- "data.frame"

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

    out <- stats::confint(x, ..., columns = newcolumns)
    class(out) <- "data.frame"
    return(out)
}

##----------------------------------------------------------------------
### model.tables.R ends here

