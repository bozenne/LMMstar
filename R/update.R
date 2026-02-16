### update.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 17 2025 (11:39) 
## Version: 
## Last-Updated: Feb 13 2026 (15:21) 
##           By: Brice Ozenne
##     Update #: 45
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
update.ID <- function(object, var.cluster, var.time, var.strata, n.time, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** what to update
    ## *** strata
    if(length(var.strata)==0){
        update.strata <- FALSE
    }else if(identical(object$name$strata,NA)){
        update.strata <- TRUE
    }else{
        stop("Argument \'repetition\' should not contain a left hand side to indicate the strata as a strata variable has already been defined in the \'structure\' argument. \n")
    }
    
    ## *** time
    ## nothing

    ## ** update
    ## *** strata
    if(update.strata){

        object$formula <- .formulaStructure(list(variance = stats::update(object$formula$variance, paste0(paste(var.strata,collapse="+"),"~.")),
                                                 correlation = NULL),
                                            correlation = FALSE)

    }

    ## *** time
    ## nothing

    ## *** cluster
    if(is.na(object$name$cluster) && !is.na(var.cluster)){
        object$name$cluster <- var.cluster
    }

    ## *** time
    if(length(object$name$time)==1 && is.na(object$name$time) && any(!is.na(var.time))){
        object$name$time <- na.omit(var.time)
    }
    
    ## ** export
    return(object)
}

## * update.IND
##' @description update structure according to information from the repetition argument
##' @details May add argument \code{var.time} as covariate.
##' @noRd
update.IND <- function(object, var.cluster, var.time, var.strata, n.time, ...){

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** what to update
    ## *** strata
    if(length(var.strata)==0){
        update.strata <- FALSE
    }else if(identical(object$name$strata,NA)){
        update.strata <- TRUE
    }else{
        stop("Argument \'repetition\' should not contain a left hand side to indicate the strata as a strata variable has already been defined in the \'structure\' argument. \n")
    }
    browser()
    ## *** time
    if(xx){
        add.time <- n.time>1 && length(object$name$time)==1 && is.na(object$name$time) && !is.null(attr(var.time,"original")) && all(!is.na(attr(var.time,"original")))
    }else{
        add.time <- FALSE
    }

    ## ** update
    ## *** strata
    if(update.strata){

        object$formula <- .formulaStructure(list(variance = stats::update(object$formula$variance, paste0(paste(var.strata,collapse="+"),"~.")),
                                                 correlation = NULL),
                                            correlation = FALSE)

    }

    ## *** cluster
    if(is.na(object$name$cluster) && !is.na(var.cluster)){
        object$name$cluster <- var.cluster
    }

    ## *** time
    if(length(object$name$time)==1 && is.na(object$name$time) && any(!is.na(var.time))){
        object$name$time <- na.omit(var.time)
    }
    
    ## ** export
    return(object)
}

update.IND <- function(object, var.cluster, var.time, var.strata, n.time, ...){



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
update.CS <- function(object, var.cluster, var.time, var.strata, n.time, ...){
browser()
    if(object$cross %in% c("TOEPLITZ","DUN","UN")){
        ## add var.time as covariate
        out <- update.IND(object = object, var.cluster = var.cluster, var.time = var.time, var.strata = var.strata, n.time = n.time, ...)
    }else{
        ## does not add var.time as covariate
        out <- update.ID(object = object, var.cluster = var.cluster, var.time = var.time, var.strata = var.strata, n.time = n.time, ...)
    }

    ## ** export
    return(out)

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
##' @noRd
update.TOEPLITZ <- update.IND

## * update.UN
##' @noRd
update.UN <- update.IND

## * update.CUSTOM
##' @noRd
update.CUSTOM <- update.ID

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
