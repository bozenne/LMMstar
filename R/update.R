### update.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 17 2025 (11:39) 
## Version: 
## Last-Updated: jul 24 2025 (17:30) 
##           By: Brice Ozenne
##     Update #: 27
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

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
