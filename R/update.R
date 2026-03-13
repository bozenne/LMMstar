### update.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 17 2025 (11:39) 
## Version: 
## Last-Updated: mar 13 2026 (14:03) 
##           By: Brice Ozenne
##     Update #: 217
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * update.ID
##' @description update structure according to information from the repetition argument
##' @noRd
update.ID <- function(object, var.cluster, var.time, var.strata, ...){

    ## ** update

    ## *** cluster 
    object$name$cluster <- var.cluster

    ## *** ordering
    object$name$ordering <- var.time
    
    ## *** prepare strata
    update.strata <- sum(!is.na(setdiff(var.strata,object$name$strata)))>0

    ## *** formula/name and strata
    if(update.strata){
        ## if both repetition and structure contained a strata, it ought to be the same (see check in .lmmNormalizeArgs)
        ## so no need to add to the existing strata
        object$name$strata <- var.strata
        ## use .formulaStructure to obtain ~0+strata instead of 'just' adding strata to the formula
        outCov <- .formulaStructure(list(variance = formula.variance, correlation = NULL, correlation.cross = NULL),
                                    var.time = list(variance = var.time, correlation = NULL, correlation.cross = NULL),
                                    correlation = FALSE)
        object$formula$variance <- outCov$formula$variance
        object$name$variance <- list(outCov$name$variance)
    }
    
    ## ** export
    return(object)
}

## * update.IND
##' @description update structure according to information from the repetition argument
##' @details May add argument \code{var.time} as covariate.
##' @noRd
update.IND <- function(object, var.cluster, var.time, var.strata, ...){

    ## ** update

    ## *** cluster 
    object$name$cluster <- var.cluster

    ## *** ordering
    object$name$ordering <- var.time
    
    ## *** prepare strata
    update.strata <- sum(!is.na(setdiff(var.strata,object$name$strata)))>0

    ## *** prepare formula
    update.time <- attr(object$formula,"update.time")    
    if(update.time["variance"]){
        if(!is.null(attr(var.time,"original")) && any(!is.na(attr(var.time,"original")))){
            var.time <- stats::na.omit(attr(var.time,"original"))            
        }else if(any(!is.na(var.time))){
            var.time <- stats::na.omit(var.time)
        }        
        
    }else{
        var.time <- NULL
    }

    ## *** formula/name and strata
    if(update.strata || update.time["variance"]){
        ## if both repetition and structure contained a strata, it ought to be the same (see check in .lmmNormalizeArgs)
        ## so no need to add to the existing strata
        if(update.strata){
            object$name$strata <- var.strata
            formula.variance <- stats::update(object$formula$variance, paste0(paste(var.strata,collapse="+"),"~."))
        }else{
            formula.variance <- object$formula$variance
        }
        ## use .formulaStructure to obtain ~0+strata instead of 'just' adding strata to the formula
        ## If both repetition and structure contained a strata, it ought to be the same (see check in .lmmNormalizeArgs)
        outCov <- .formulaStructure(list(variance = formula.variance, correlation = NULL, correlation.cross = NULL),
                                    ls.time = list(variance = var.time, correlation = NULL, correlation.cross = NULL),
                                    correlation = FALSE)
        object$formula$variance <- outCov$formula$variance
        object$name$variance <- list(outCov$name$variance)
    }

    ## ** export
    attr(object$formula,"update.time") <- NULL
    return(object)
}


## * update.CS
##' @description update structure according to information from the repetition argument
##' @noRd
update.CS <- function(object, var.cluster, var.time, var.strata, ...){

    ## ** update

    ## *** cluster 
    object$name$cluster <- var.cluster

    ## *** ordering
    object$name$ordering <- var.time

    ## *** prepare strata
    update.strata <- sum(!is.na(setdiff(var.strata,object$name$strata)))>0

    ## *** time
    update.time <- attr(object$formula,"update.time")    
    if(update.time["correlation"] || update.time["correlation.cross"]){
        if(!is.null(attr(var.time,"original")) && any(!is.na(attr(var.time,"original")))){
            object$name$time <- stats::na.omit(attr(var.time,"original"))
        }else if(any(!is.na(var.time))){
            object$name$time <- stats::na.omit(var.time)
        }        
    }

    ls.time <- stats::setNames(vector(mode = "list",length = 3), c("variance","correlation","correlation.cross"))
    if(!is.null(attr(var.time,"original")) && any(!is.na(attr(var.time,"original")))){
        var.time <- stats::na.omit(attr(var.time,"original"))
    }
    for(iMoment in c("variance","correlation","correlation.cross")){
        if(update.time[iMoment]){
            if(any(var.time %in% object$name[[iMoment]][[1]])){
                stop("Time variable \"",paste(intersect(var.time,object$name[[iMoment]][[1]]),collapse="\", \""),"\" already present in the ",iMoment," structure. \n",
                     ifelse(iMoment=="variance","Consider setting the argument \'heterogeneous\' to FALSE when creating the covariance structure. \n",""))
            }
            ls.time[[iMoment]] <- var.time
        }
    }
 
    ## *** strata and formula
    if(update.strata || any(update.time)){
        ## if both repetition and structure contained a strata, it ought to be the same (see check in .lmmNormalizeArgs)
        ## so no need to add to the existing strata
        if(update.strata){
            object$name$strata <- var.strata
            formula.variance <- stats::update(object$formula$variance, paste0(paste(var.strata,collapse="+"),"~."))
            formula.correlation <- stats::update(object$formula$correlation, paste0(paste(var.strata,collapse="+"),"~."))
            if(!is.null(object$formula$correlation.cross)){
                formula.correlation.cross <- stats::update(object$formula$correlation.cross, paste0(paste(var.strata,collapse="+"),"~."))
            }else{
                formula.correlation.cross <- NULL
            }
        }else{
            formula.variance <- object$formula$variance
            formula.correlation <- object$formula$correlation
            formula.correlation.cross <- object$formula$correlation.cross
        }
        ## use .formulaStructure to obtain ~0+strata instead of 'just' adding strata to the formula
        ## If both repetition and structure contained a strata, it ought to be the same (see check in .lmmNormalizeArgs)
        outCov <- .formulaStructure(list(variance = formula.variance, correlation = formula.correlation, correlation.cross = formula.correlation.cross),
                                    ls.time = ls.time, correlation = TRUE)
   
        if(update.time["variance"] || update.strata){
            object$formula$variance <- outCov$formula$variance
            object$name$variance <- list(outCov$name$variance)
        }
        if(update.time["correlation"] || update.strata){
            object$formula$correlation <- outCov$formula$correlation
        }
        if(update.time["correlation.cross"] || update.strata){
            object$formula$correlation.cross <- outCov$formula$correlation.cross
            
        }
        if(update.time["correlation"] || update.time["correlation.cross"] || update.strata){
            object$name$correlation <- list(outCov$name$correlation)
        }
    }

    ## ** export
    attr(object$formula,"update.time") <- NULL
    return(object)
}

## * update.RE
##' @description update structure according to information from the repetition argument
##' @details May add argument \code{var.time} as covariate.
##' @noRd
##' @param ranef Random effect structure identified via the formula argument of lmm (mean structure).
update.RE <- function(object, var.cluster, var.time, var.strata, ranef, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** what to update
    if(!identical(sort(object$name$strata),sort(var.strata))){
        update.strata <- TRUE
        rm.strata <- unique(stats::na.omit(c(object$name$strata,var.strata)))
    }else{
        update.strata <- FALSE
    }

    if(!missing(ranef) && !is.null(ranef) && is.null(object$ranef)){
        ## handle the case where structure = RE(~strata) or RE(strata~1) whereas formula = Y ~ (1|id/session)
        ## one needs to update the correlation formula with the ranef
        add.RE <- TRUE
    }else{
        add.RE <- FALSE
    }

    ## ** update
    if(update.strata || add.RE){
        call.structure <- object$call
        ranef.structure <- object$ranef

        ls.call.structure <- as.list(call.structure)
        fct.structure <- eval(ls.call.structure[[1]])
        args.structure <- lapply(ls.call.structure[-1], eval)

        if("var.cluster" %in% names(args.structure) == FALSE){
            args.structure$var.cluster <- var.cluster
        }
        if("var.time" %in% names(args.structure) == FALSE){
            args.structure$var.time <- var.time
        }

        if(update.strata){
            if(is.list(args.structure$formula)){
                args.structure$formula <- list(updateFormula(args.structure$formula[[1]], drop.y = TRUE, drop.x = rm.strata, add.y = var.strata),
                                               updateFormula(args.structure$formula[[2]], drop.y = TRUE, drop.x = rm.strata, add.y = var.strata, add.x = ranef$term))
            }else if(inherits(args.structure$formula,"formula")){
                args.structure$formula <- updateFormula(args.structure$formula, drop.y = TRUE, drop.x = rm.strata, add.y = var.strata, add.x = ranef$term)
            }
        }

        object <- do.call(fct.structure, args = args.structure)
        object$call <- call.structure
        if(add.RE){
            object$ranef <- ranef
        }
    }else{
        if(is.na(object$name$cluster) && !is.na(var.cluster)){
            object$name$cluster <- var.cluster
        }
        if(length(object$name$time)==1 && is.na(object$name$time) && all(!is.na(var.time))){
            object$name$time <- var.time
        }
    }

    ## ** export
    return(object)

}

## * update.TOEPLITZ
##' @description update structure according to information from the repetition argument
##' @details May add argument \code{var.time} as covariate.
##' @noRd
update.TOEPLITZ <- update.CS

## * update.UN
##' @noRd
update.UN <- update.CS

## * update.DUN
##' @noRd
update.DUN <- update.CS

## * update.EXP
##' @noRd
update.EXP <- update.CS

## * update.AR1
##' @noRd
update.AR1 <- update.CS

## * update.CUSTOM
##' @noRd
update.CUSTOM <- update.CS

## * update.lmm
##' @export
update.lmm <- function(object, ...){

    out <- try(stats::update.default(object, ...), silent = TRUE)
    if(inherits(out,"try-error")){ ## try to handle not finding the dataset by using what is stored in the object
        if(inherits(object$call$data,"name") && inherits(try(get(object$call$data), silent = TRUE), "try-error")){
            object$call$data <- object$data.original
            out <- stats::update.default(object, ...)
        }else{
            stop(out)
        }
    }
    return(out)

}


## * update.rbindWald_lmm
##' @export
update.rbindWald_lmm <- function(object, formula., ..., evaluate){

    ## ** normalize input
    dots <- list(...)
    if(length(dots)!=1 || "p" %in% names(dots) == FALSE){
        stop("Unknonwn input: can only update the rbindWald_lmm object w.r.t. to \'p\'. \n")
    }
    if(!missing(formula.)){
        message("Argument \'formula.\' is ignored when updating the rbindWald_lmm object. \n")
    }
    if(!missing(evaluate)){
        message("Argument \'evaluate\' is ignored when updating the rbindWald_lmm object. \n")
    }
    p <- dots$p
    table.param <- stats::model.tables(object, effects = "param")
    if(length(p) != NROW(table.param)){
        stop("Incorrect length for argument \'p\': should have length ",NROW(table.param)," instead of ",length(p),".\n",
             "There are ",NROW(table.param),"=",paste(table(table.param$model),collapse="+")," parameters w.r.t. the linear mixed models. \n")
    }
    if(any(names(p) %in% table.param$Uname == FALSE)){
        stop("Incorrect names for argument \'p\': \"",paste(setdiff(names(p),table.param$Uname), collapse = "\", \""),"\".\n",
             "Example of valid names: ",paste(setdiff(table.param$Uname,names(p)), collapse = "\", \""),"\n")
    }
    if(any(table.param$Uname %in% names(p) == FALSE)){
        stop("Missing parameter in argument \'p\': \"",paste(setdiff(table.param$Uname, names(p)), collapse = "\", \""),"\".\n")
    }

    ## ** update anova
    ## retrieve call
    call.anova <- object$call[setdiff(names(object$call),"rbind")]

    ## retrieve lmm
    ls.lmm <- lmm(object)
    n.lmm <- length(ls.lmm)

    ## update
    ls.anova <- lapply(1:n.lmm, function(iLMM){ ## iLMM <-  1
        iTable <- table.param[table.param$model==iLMM,]

        iArgs <- as.list(call.anova[[iLMM]])[-1]
        iArgs$object <- ls.lmm[[iLMM]]
        iArgs$p <- stats::setNames(p[iTable$Uname],iTable$name)
        
        return(do.call(anova, args = iArgs))
    })

    ## ** combine anova
    ## retrieve call
    call.rbind <- object$call[["rbind"]]

    ## reformat arguments
    keep.argsRBIND <- setdiff(names(formals(fun = rbind.Wald_lmm)),c("model","..."))
    orginal.argsRBIND <- as.list(call.rbind)[-1][]
    if(any(names(orginal.argsRBIND) %in% keep.argsRBIND)){
        ls.argsRBIND <- orginal.argsRBIND[intersect(names(orginal.argsRBIND),keep.argsRBIND)]
    }else{
        ls.argsRBIND <- NULL
    }
    out <- do.call(rbind.Wald_lmm, args = c(list(model = ls.anova), ls.argsRBIND))

    ## ** export
    return(out)
}

## * update.rbindWald_lmm
##' @export
update.mlmm <- function(object, formula., ..., evaluate){

    ## ** normalize input
    dots <- list(...)
    if(is.null(names(dots)) || "" %in% names(dots) || any(is.na(names(dots)))){
        stop("Argument used to update the mlmm object should be named. \n")
    }

    ## ** update fit
    ls.args <- c(as.list(object$call)[-1],
                 dots[setdiff(names(dots),names(object$call))])
    ls.args[intersect(names(ls.args),names(dots))] <- dots[names(dots)]

    out <- do.call(mlmm, args = ls.args)

    ## ** export
    return(out)
}
##----------------------------------------------------------------------
### update.R ends here
