### structure-update.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2023 (14:01) 
## Version: 
## Last-Updated: May 12 2024 (19:16) 
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

## * update.structure
##' @description update structure according to information from the repetition argument
##'
##' @details update.ID and update.IND differs in that update.IND may add the time as covariate
##' 
##' @noRd

## * update.ID
##' @noRd
update.ID <- function(object, var.cluster, var.time, var.strata, n.time, ...){

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

    ## ** update
    if(update.strata){
        call.structure <- object$call
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
                                               updateFormula(args.structure$formula[[2]], drop.y = TRUE, drop.x = rm.strata, add.y = var.strata))
            }else if(inherits(args.structure$formula,"formula")){
                args.structure$formula <- updateFormula(args.structure$formula, drop.y = TRUE, drop.x = rm.strata, add.y = var.strata)
            }
        }

        object <- do.call(fct.structure, args = args.structure)
        object$call <- call.structure
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

## * update.IND
##' @noRd
update.IND <- function(object, var.cluster, var.time, var.strata, n.time, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** what to update
    if(n.time>1 && length(object$name$time)==1 && is.na(object$name$time) && !is.null(attr(var.time,"original")) && all(!is.na(attr(var.time,"original")))){
        add.time <- TRUE
    }else{
        add.time <- FALSE
    }

    if(!identical(sort(object$name$strata),sort(var.strata))){
        update.strata <- TRUE
        rm.strata <- unique(stats::na.omit(c(object$name$strata,var.strata)))
    }else{
        update.strata <- FALSE
    }

    ## ** update
    if(add.time || update.strata){
        call.structure <- object$call
        ls.call.structure <- as.list(call.structure)
        fct.structure <- eval(ls.call.structure[[1]])
        args.structure <- lapply(ls.call.structure[-1], eval)

        if("var.cluster" %in% names(args.structure) == FALSE){
            args.structure$var.cluster <- var.cluster
        }
        if("var.time" %in% names(args.structure) == FALSE){
            args.structure$var.time <- var.time
        }
        if("add.time" %in% names(args.structure) == FALSE && !is.list(args.structure$formula)){
            args.structure$add.time <- attr(var.time,"original")
        }

        if(update.strata){
            if(is.list(args.structure$formula)){
                args.structure$formula <- list(updateFormula(args.structure$formula[[1]], drop.y = TRUE, drop.x = rm.strata, add.y = var.strata),
                                               updateFormula(args.structure$formula[[2]], drop.y = TRUE, drop.x = rm.strata, add.y = var.strata))
            }else if(inherits(args.structure$formula,"formula")){
                args.structure$formula <- updateFormula(args.structure$formula, drop.y = TRUE, drop.x = rm.strata, add.y = var.strata)
            }
        }
        object <- do.call(fct.structure, args = args.structure)
        object$call <- call.structure
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

## * update.CS
##' @noRd
update.CS <- update.ID

## * update.RE
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
##' @noRd
update.TOEPLITZ <- update.IND

## * update.UN
##' @noRd
update.UN <- update.IND

## * update.CUSTOM
##' @noRd
update.CUSTOM <- update.ID

##----------------------------------------------------------------------
### structure-update.R ends here
