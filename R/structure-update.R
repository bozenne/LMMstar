### structure-update.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2023 (14:01) 
## Version: 
## Last-Updated: jul  7 2023 (12:01) 
##           By: Brice Ozenne
##     Update #: 99
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
update.ID <- function(object, var.cluster, var.time, var.strata, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** what to update
    if(is.na(object$name$strata) && !is.null(var.strata) && !is.na(var.strata)){
        add.strata <- TRUE
    }else{
        add.strata <- FALSE
    }

    ## ** update
    if(add.strata){
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

        if(add.strata && is.list(args.structure$formula)){
            args.structure$formula <- list(stats::update(args.structure$formula[[1]], stats::as.formula(paste0(var.strata,"~."))),
                                           stats::update(args.structure$formula[[2]], stats::as.formula(paste0(var.strata,"~."))))
        }else if(add.strata && inherits(args.structure$formula,"formula")){
            args.structure$formula <- stats::update(args.structure$formula, stats::as.formula(paste0(var.strata,"~.")))
        }

        object <- do.call(fct.structure, args = args.structure)
        object$call <- call.structure
    }else{
        if(is.na(object$name$cluster) && !is.na(var.cluster)){
            object$name$cluster <- var.cluster
        }
        if(is.na(object$name$time) && !is.na(var.time)){
            object$name$time <- var.time
        }
    }

    ## ** export
    return(object)
}

## * update.IND
update.IND <- function(object, var.cluster, var.time, var.strata, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** what to update
    if(is.na(object$name$time) && !is.null(attr(var.time,"original")) && !is.na(attr(var.time,"original"))){
        add.time <- TRUE
    }else{
        add.time <- FALSE
    }

    if(is.na(object$name$strata) && !is.null(var.strata) && !is.na(var.strata)){
        add.strata <- TRUE
    }else{
        add.strata <- FALSE
    }

    ## ** update
    if(add.time || add.strata){
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
        if("add.time" %in% names(args.structure) == FALSE){
            args.structure$add.time <- attr(var.time,"original")
        }
        if(add.strata && is.list(args.structure$formula)){
            args.structure$formula <- list(stats::update(args.structure$formula[[1]], stats::as.formula(paste0(var.strata,"~."))),
                                           stats::update(args.structure$formula[[2]], stats::as.formula(paste0(var.strata,"~."))))
        }else if(add.strata && inherits(args.structure$formula,"formula")){
            args.structure$formula <- stats::update(args.structure$formula, stats::as.formula(paste0(var.strata,"~.")))
        }

        object <- do.call(fct.structure, args = args.structure)
        object$call <- call.structure
    }else{
        if(is.na(object$name$cluster) && !is.na(var.cluster)){
            object$name$cluster <- var.cluster
        }
        if(is.na(object$name$time) && !is.na(var.time)){
            object$name$time <- var.time
        }
    }

    ## ** export
    return(object)
}

## * update.CS
update.CS <- update.ID

## * update.RE
##' @param ranef Random effect structure identified via the formula argument of lmm (mean structure).
##' @noRd
update.RE <- function(object, var.cluster, var.time, var.strata, ranef, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** what to update
    if(is.na(object$name$strata) && !is.null(var.strata) && !is.na(var.strata)){
        add.strata <- TRUE
    }else{
        add.strata <- FALSE
    }

    if(!missing(ranef) && !is.null(ranef) && is.null(object$ranef)){
        ## handle the case where structure = RE(~strata) or RE(strata~1) whereas formula = Y ~ (1|id/session)
        ## one needs to update the correlation formula with the ranef
        add.RE <- TRUE
    }else{
        add.RE <- FALSE
    }

    ## ** update
    if(add.strata || add.RE){
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
        if(add.RE){
            if(add.strata==FALSE){
                if(is.na(object$name$strata)){
                    var.strata <- NULL
                }else{
                    var.strata <- object$name$strata
                }
            }
            if(is.list(args.structure$formula)){
                args.structure$formula <- list(stats::as.formula(paste0(var.strata,"~",paste(ranef$term,collapse=" + "))),
                                               stats::as.formula(paste0(var.strata,"~",paste(ranef$term,collapse=" + "))))
            }else if(inherits(args.structure$formula,"formula")){
                args.structure$formula <- stats::as.formula(paste0(var.strata,"~",paste(ranef$term,collapse=" + ")))
            }
        }else if(add.strata){
            if(is.list(args.structure$formula)){
                args.structure$formula <- list(stats::update(args.structure$formula[[1]], stats::as.formula(paste0(var.strata,"~."))),
                                               stats::update(args.structure$formula[[2]], stats::as.formula(paste0(var.strata,"~."))))
            }else if(inherits(args.structure$formula,"formula")){
                args.structure$formula <- stats::update(args.structure$formula, stats::as.formula(paste0(var.strata,"~.")))
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
        if(is.na(object$name$time) && !is.na(var.time)){
            object$name$time <- var.time
        }
    }

    ## ** export
    return(object)

}

## * update.TOEPLITZ
update.TOEPLITZ <- update.IND

## * update.UN
update.UN <- update.IND

## * update.CUSTOM
update.CUSTOM <- update.ID

##----------------------------------------------------------------------
### structure-update.R ends here
