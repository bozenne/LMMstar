### structure-update.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2023 (14:01) 
## Version: 
## Last-Updated: jul  6 2023 (18:03) 
##           By: Brice Ozenne
##     Update #: 37
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
    if(is.na(object$name$strata) && !is.null(attr(var.strata,"original")) && !is.na(attr(var.strata,"original"))){
        add.strata <- TRUE
    }else{
        add.strata <- FALSE
    }
        
    ## ** update
    if(add.strata){
        call.structure <- as.list(object$call)
        fct.structure <- eval(call.structure[[1]])
        args.structure <- eval(call.structure[-1])

        if("var.cluster" %in% names(args.structure) == FALSE){
            args.structure$var.cluster <- var.cluster
        }
        if("var.time" %in% names(args.structure) == FALSE){
            args.structure$var.time <- var.time
        }

        if(add.strata && is.list(args.structure$formula)){
            args.structure$formula <- list(stats::update(stats::as.formula(args.structure$formula[[1]], stats::as.formula(paste0(var.strata,"~.")))),
                                           stats::update(stats::as.formula(args.structure$formula[[2]], stats::as.formula(paste0(var.strata,"~.")))))
        }else if(add.strata && inherits(args.structure$formula,"formula")){
            args.structure$formula <- stats::update(stats::as.formula(args.structure$formula), stats::as.formula(paste0(var.strata2,"~.")))
        }

        object <- do.call(fct.structure, args = args.structure)

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

    if(is.na(object$name$strata) && !is.null(attr(var.strata,"original")) && !is.na(attr(var.strata,"original"))){
        add.strata <- TRUE
    }else{
        add.strata <- FALSE
    }

    ## ** update
    if(add.time || add.strata){
        call.structure <- as.list(object$call)
        fct.structure <- eval(call.structure[[1]])
        args.structure <- eval(call.structure[-1])

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
            args.structure$formula <- list(stats::update(stats::as.formula(args.structure$formula[[1]], stats::as.formula(paste0(var.strata,"~.")))),
                                           stats::update(stats::as.formula(args.structure$formula[[2]], stats::as.formula(paste0(var.strata,"~.")))))
        }else if(add.strata && inherits(args.structure$formula,"formula")){
            args.structure$formula <- stats::update(stats::as.formula(args.structure$formula), stats::as.formula(paste0(var.strata2,"~.")))
        }

        object <- do.call(fct.structure, args = args.structure)

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
update.RE <- function(object, var.cluster, var.time, var.strata, ranef, ...){
    browser()
    if(!missing(ranef)){
        object$name$strata <- NA ## will be updated via var.strata
        object$name$cor <- list(setdiff(ranef$vars,attr(var.cluster,"original")))
        object$formula$var <- ~1
        object$formula$cor <- stats::as.formula(paste("~0 +", paste(object$name$cor[[1]], collapse = "+")))
    }

    object <- update.ID(object = object, var.cluster = var.cluster, var.time = var.time, var.strata = var.strata, ...)
    

}

## * update.TOEPLITZ
update.TOEPLITZ <- update.IND

## * update.UN
update.UN <- update.IND

## * update.CUSTOM
update.CUSTOM <- update.ID

##----------------------------------------------------------------------
### structure-update.R ends here
