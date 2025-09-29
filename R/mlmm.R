### mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 14 2022 (09:45) 
## Version: 
## Last-Updated: sep 29 2025 (13:14) 
##           By: Brice Ozenne
##     Update #: 591
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mlmm (documentation)
##' @title Fit Multiple Linear Mixed Model
##' @description Fit several linear mixed models, extract relevant coefficients, and combine them into a single table. 
##'
##' @param ... arguments passed to \code{\link{lmm.formula}}.
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param by [character] variable used to split the dataset. On each split a seperate linear mixed model is fit.
##' @param effects [character or numeric matrix] Linear combinations of coefficients relative to which Wald test should be computed.
##' Argument passed to \code{\link{anova.lmm}}.
##' Right hand side can be specified via an attribute \code{"rhs"}.
##' @param contrast.rbind [character or numeric matrix] Contrast to be be applied to compare the groups.
##' Argument passed to the argument \code{effects} of \code{\link{rbind.Wald_lmm}}.
##' Right hand side can be specified via an attribute \code{"rhs"}.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors.
##' Can also be \code{2} compute the degrees-of-freedom w.r.t. robust standard errors instead of w.r.t. model-based standard errors.
##' Argument passed to \code{\link{anova.lmm}}.
##' @param df [logical] Should the degrees-of-freedom be computed using a Satterthwaite approximation?
##' Argument passed to \code{\link{anova.lmm}}.
##' @param p [list of numeric vector] parameter values for each model. For internal use (delta-method via \code{\link{estimate.mlmm}}).
##' @param name.short [logical] use short names for the output coefficients (omit the name of the by variable).
##' @param trace [interger, >0] Show the progress of the execution of the function.
##' @param transform.sigma,transform.k,transform.rho,transform.names [character] transformation used on certain type of parameters.
##' 
##' @details
##' \bold{Grouping variable} in argument repetition: when numeric, it will be converted into a factor variable, possibly adding a leading 0 to preserve the ordering.
##' This transformation may cause inconsistency when combining results between different \code{lmm} object. 
##' This is why the grouping variable should preferably be of type character or factor.
##' 
##' @seealso
##' \code{\link{confint.mlmm}} for a data.frame containing estimates with their uncertainty. \cr
##' \code{\link{summary.mlmm}} for a summary of the model and estimates. \cr
##' \code{\link{autoplot.Wald_lmm}} for a graphical display. \cr
##' 
##' @keywords models
##' 
##' @examples
##' #### univariate regression ####
##' if(require(lava) && require(multcomp)){
##' 
##' set.seed(10)
##' d1 <- cbind(sim(lvm(Y~0.5*X1), 25), group = "A")
##' d2 <- cbind(sim(lvm(Y~0.1*X1), 100), group = "B")
##' d3 <- cbind(sim(lvm(Y~0.01*X1), 1000), group = "C")
##' d1$id <- 1:NROW(d1)
##' d2$id <- 1:NROW(d2)
##' d3$id <- 1:NROW(d3)
##' 
##' d <- rbind(d1,d2,d3)
##' 
##' e.mlmm <- mlmm(Y~X1, data = d, repetition = ~1|id,
##'                by = "group", effects = "X1=0")
##' summary(e.mlmm)
##' summary(e.mlmm, method = "single-step")
##' summary(e.mlmm, method = "single-step2")
##'
##' ## re-work contrast
##' summary(anova(e.mlmm, effects = mcp(X1 = "Dunnett")), method = "none")
##' ## summary(mlmm(Y~X1, data = d, by = "group", effects = mcp(X1="Dunnett")))
##' }
##' 
##' #### multivariate regression ####
##' set.seed(10)
##' dL <- sampleRem(250, n.times = 3, format = "long")
##'
##' e.mlmm <- mlmm(Y~X1+X2+X6, repetition = ~visit|id, data = dL,
##'                by = "X4", structure = "CS")
##' summary(e.mlmm)
##' 
##' e.mlmmX1 <- mlmm(Y~X1+X2+X6, repetition = ~visit|id, data = dL,
##'                by = "X4", effects = "X1=0", structure = "CS")
##' summary(e.mlmmX1)
##' summary(e.mlmmX1, method = "single-step")
##' 

## * mlmm (code)
##' @export
mlmm <- function(..., data, by, contrast.rbind = NULL, effects = NULL, robust = FALSE, df = NULL, p = NULL,
                 name.short = TRUE, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, trace = TRUE){

    options <- LMMstar.options()

    ## ** normalizer user input
    ## *** data
    if(!inherits(data,"data.frame")){
        stop("Argument \'data\' must inherit from \"data.frame\". \n")
    }
    if(any(by %in% names(data) == FALSE)){
        stop("Mismatch between argument \'by\' and \'data\'.\n",
             "Could not find column \"",paste(by[by %in% names(data) == FALSE], collapse = "\" \""),"\" in data \n")
    }
    reserved.names <- c("by","type","test","estimate","se","df","statistic","lower","upper","null","p.value")
    if(any(names(data) %in% reserved.names)){
        stop("Argument \'data\' should not contain a column named \"",paste(names(data)[names(data) %in% reserved.names], collapse = "\" \""),
             "\" as this name is used internally by the mlmm function. \n")
    }

    ## *** by
    if(length(by)>1){
        by.keep <- by
        by <- paste(by.keep,collapse=",")            
        if(by %in% names(data)){
            stop("Argument \'data\' should not contain a column named \"",by,"\" as this name is used internally by the mlmm function. \n")
        }
        if(any(duplicated(by.keep))){
            stop("Argument \'by\' should not contain duplicated values. \n")
        }
        data[[by]] <- nlme::collapse(data[by.keep], sep=",", as.factor = TRUE)
    }else{
        by.keep <- by
        if(is.factor(data[[by]])){
            data[[by]] <- droplevels(data[[by]])
        }
    }

    ## *** df
    if(is.null(df)){
        df <- options$df
    }
    
    ## *** p
    if(!is.null(p) && !is.list(p)){
        stop("Argument \'p\' should either be NULL or a list with the parameter values for each model. \n")
    }
    if(!is.null(p) && is.null(names(p))){
        stop("Argument \'p\' should either be NULL or a named list with the parameter values for each model. \n")
    }

    ## ** fit mixed models
    repetition <- list(...)$repetition
    if(is.null(repetition)){
        ## add cluster variable so that each individual is in a different cluster
        ## other the first from each split may be considered from the same cluster
        data$XXclusterXX <- addLeading0(1:NROW(data), as.factor = FALSE, code = NULL)
    }else if(!is.null(repetition)){
        ## avoid that cluster variable is transformed differently across linear mixed models
        ## 1 --> "01" (when 10+ but 100- individuals) vs. 1 --> "001" (when 100+ but 1000- individuals)
        detail.repetition <- formula2var(repetition)
        if(detail.repetition$special=="repetition"){
            var.cluster <- detail.repetition$var$cluster
            if(var.cluster %in% names(data) && is.numeric(data[[var.cluster]])){
                if(all(data[[var.cluster]]>0) && all(data[[var.cluster]] %% 1 == 0)){
                    data[[var.cluster]] <- addLeading0(data[[var.cluster]], as.factor = FALSE, code = attr(var.cluster,"code"))
                }else{
                    data[[var.cluster]] <- as.character(data[[var.cluster]])
                }
            }
        }
    }

    ls.data <- by(data, data[[by]], function(iDF){
        if(is.factor(iDF[[by]])){ ## ensure that by variable is a character and no a factor
            iDF[[by]] <- as.character(iDF[[by]])
        }
        return(iDF)
    })

    if(trace>0){
        dots <- list(...)
        if("repetition" %in% names(dots) && all.vars(stats::delete.response(stats::terms(dots[["repetition"]])))[1] != by){
            cat("Fitting linear mixed models:\n")
        }else{
            cat("Fitting linear regressions:\n")
        }
    }
    ls.lmm <- lapply(ls.data, function(iData){ ## iData <- ls.data[[2]]
        if(trace>0.5){
            cat(" - ",by,"=",unique(iData[[by]]),"\n", sep = "")
        }
        return(lmm(..., data = iData, df = df, trace = max(0,trace-1)))
    })
    name.lmm <- names(ls.lmm)
    n.lmm <- length(name.lmm)
    variable.lmm <- stats::variable.names(ls.lmm[[1]])

    if(!is.null(p) && length(p) != n.lmm){
        stop("If not NULL, the argument \'p\' should be a list whose length equal the number of strata (here ",n.lmm,"). \n")
    }
    if(!is.null(p) && any(sort(names(p) != name.lmm))){
        stop("If not NULL, invalid names in argument \'p\': \"",paste(setdiff(names(p),name.lmm), collapse ="\", \""),"\"\n",
             "Example of valid names: \"",paste(setdiff(name.lmm,names(p)), collapse ="\", \""),"\"\n")
    }

    ## ** test linear combinations
    if(trace>0){
        cat("\nHypothesis test:\n")
    }
    ls.name.allCoef <- lapply(ls.lmm, function(iLMM){stats::model.tables(iLMM, effects = "param")$trans.name})

    ## *** identify contrast
    if(trace>0.5){
        cat(" - generate contrast matrix\n")
    }
    if(inherits(effects,"mcp")){
        if(length(effects)!=1){
            stop("Argument \'effects\' must specify a single hypothesis test when being of class \"mcp\". \n",
                 "Something like mcp(group = \"Dunnett\") or mcp(group = \"Tukey\") \n")
        }
        if(!is.null(contrast.rbind)){
            message("Argument \'contrast.rbind\' ignored when argument \'effects\' inherits from class \"mcp\". \n")
        }
        effects.save <- effects
        contrast.rbind <- effects.save[[1]]
        effects <- names(effects.save)
        if(!grepl("=",effects)){
            effects <- paste0(effects,"=0")
        }
    }

    if(is.null(effects)){ ## Case 1: no user-input so export all coefficients, each in a separate test

        ls.contrast <- lapply(ls.name.allCoef, function(iName){matrix(0, nrow = length(iName), ncol = length(iName), dimnames = list(iName,iName))})

        for(iName in name.lmm){## iName <- name.lmm[[1]]
            iNameCoef <- stats::model.tables(ls.lmm[[iName]], effects = c("param",options$effects))$trans.name
            if(length(iNameCoef)==1){
                ls.contrast[[iName]][iNameCoef,iNameCoef] <- 1
            }else{
                diag(ls.contrast[[iName]][iNameCoef,iNameCoef]) <- 1
            }
            ls.contrast[[iName]] <- ls.contrast[[iName]][rowSums(abs(ls.contrast[[iName]]))!=0,,drop=FALSE] ## remove lines with only 0
        }
        ls.rhs <- NULL

    }else{

        if(length(unique(lengths(ls.name.allCoef)))>1){
            stop("Cannot use matrix interface for argument \'effects\' when the number of model parameters varies over splits. \n")
        }
        if(any(sapply(ls.name.allCoef[-1],identical,ls.name.allCoef[[1]])==FALSE)){
            stop("Cannot use matrix interface for argument \'effects\' when the model parameters varies over splits. \n")
        }
        ls.e2c <- effects2contrast(ls.lmm[[1]], effects = effects, rhs = attr(effects,"rhs"), options = options)
        
        ls.contrast <- stats::setNames(lapply(1:n.lmm, function(iI){ls.e2c$contrast$user$user}), name.lmm)
        ls.rhs <- stats::setNames(lapply(1:n.lmm, function(iI){ls.e2c$null$user$user}), name.lmm)
        
    }

    ## *** run
    if(trace>0.5){
        cat(" - univariate test\n")
    }

    ls.anova <- stats::setNames(lapply(name.lmm, function(iName){ ## iName <- name.lmm[1]
        iWald <- stats::anova(ls.lmm[[iName]], effects = ls.contrast[[iName]], rhs = ls.rhs[[iName]], robust = robust, df = df,
                              univariate = TRUE, multivariate = FALSE, p = p[[iName]],
                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                              simplify = FALSE)

        return(iWald)
    }), name.lmm)

    ## ** joint inference
    if(trace>0.5){
        cat(" - combine tests\n")
    }
    if(length(by.keep)==1){
        if(name.short[1]){
            name.model <- name.lmm
        }else{
            name.model <- paste0(by,"=",name.lmm)
        }
        
    }else{
        name.model <- unlist(lapply(ls.data, function(iData){
            if(name.short[1]){
                return(paste(iData[1,by.keep],collapse=","))
            }else{
                return(paste(paste0(by.keep,"=",iData[1,by.keep]),collapse=","))
            }
        }))
        
    }

    if(!is.null(contrast.rbind)){
        rhs.by <- attr(contrast.rbind,"rhs")
        attr(contrast.rbind,"rhs") <- NULL
        sep <- ":"
    }else{
        rhs.by <- NULL
        sep <- ": "
    }
    out <- rbind.Wald_lmm(model = ls.anova, effects = contrast.rbind, rhs = rhs.by, name = name.model, name.short = name.short, sep = sep)
    names(out$univariate)[names(out$univariate)=="model"] <- "by"

    ## *** add covariate values
    if(is.null(contrast.rbind)){
        keep.by.level <- do.call(rbind,lapply(ls.data, function(iData){iData[1,by.keep,drop=FALSE]}))
        keep.rownames <- rownames(out$univariate)
        out$univariate <- merge(x = out$univariate, y = cbind(by = name.model, keep.by.level), by = "by", all.x = TRUE, sort = FALSE)
        rownames(out$univariate) <- keep.rownames
    }

    ## ** export
    if(trace>0){
        cat("\n")
    }
    out$args$by <- NA ## deal with the case of multiple by
    out$args$by <- list(by.keep)
    out$call <- match.call()
    class(out) <- append("mlmm", class(out))
    return(out)
}
##----------------------------------------------------------------------
### mlmm.R ends here
