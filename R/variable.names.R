### variable.names.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 31 2022 (15:05) 
## Version: 
## Last-Updated: aug  8 2024 (13:32) 
##           By: Brice Ozenne
##     Update #: 128
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * variable.names.lmm (documentation)
##' @title Variables Involved in a Linear Mixed Model
##' @description Extract the variables used by the linear mixed model.
##'
##' @param object a \code{lmm} object.
##' @param effects [character] Should all variable be output (\code{"all"}),
##' or only those related to the outcome (\code{"outcome"}), mean (\code{"mean"}), variance (\code{"variance"}),
##' correlation (\code{"correlation"}), time (\code{"time"}), cluster (\code{"cluster"}), strata (\code{"strata"})?
##' @param original [logical] Should only the variables present in the original data be output?
##' When \code{FALSE}, variables internally created are output instead of the original variable for time, cluster, and strata.
##' @param simplify [logical] Should the list be converted into a vector if a single \code{effects} is requested?
##' @param ... not used. For compatibility with the generic function
##'
##' @return A list of character vectors or a character vector.
##' 
##' @keywords methods

## * variable.names.lmm 
##' @export
variable.names.lmm <- function(object, effects = "all", original = TRUE, simplify = TRUE, ...){

    ## ** extract variables
    var.outcome <- object$outcome$var
    var.cluster <- object$cluster$var
    var.time <- object$time$var
    var.strata <- object$strata$var
    n.time <- object$time$n

    ## ** check user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character vector. \n")
    }
    valid.effects <- c("outcome","mean","variance","correlation","time","cluster","strata",
                       "all","mean.type")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }
    if(all("all" %in% effects)){
        if(length(effects)>1){
            stop("Argument \'effects\' must have length 1 when containing the element \"all\". \n")
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

    ## ** prepare output
    ls.out <- list(outcome = var.outcome,
                   mean = attr(object$design$mean, "variable"),
                   variance = all.vars(object$design$vcov$formula$var),
                   correlation = all.vars(object$design$vcov$formula$cor),
                   time = var.time,
                   cluster = var.cluster,
                   strata = var.strata)
    if(original){
        ls.out$time <- attr(ls.out$time,"original")
        ls.out$cluster <- attr(ls.out$cluster,"original")
        ls.out$strata <- attr(ls.out$strata,"original")
    }
    ## normalize numeric(0) and NA into NULL
    ls.out[sapply(ls.out, function(iE){sum(!is.na(iE))==0})] <- list(NULL)

    ## ** time-invariant variables
    if("mean.type" %in% effects){
        ## hidden: no returned by effects = "all" to save computation time
        ## only internal function need effects="mean.type" (predict,model.frame) for reshaping data
        n.time <- object$time$n
        vars.time <- stats::na.omit(unique(c(var.time,attr(var.time,"original"))))
        vars.cluster <- stats::na.omit(unique(c(var.cluster,attr(var.cluster,"original"))))
        vars.mean_notime <- ls.out$mean
        data <- object$data

        ls.out$mean.type <- stats::setNames(rep("baseline", length(ls.out$mean)), ls.out$mean)
        if(length(vars.mean_notime)>0){
            ## find covariates that are a simple transformation of time
            cov.levels <- apply(data[vars.mean_notime], MARGIN = 2, function(iCol){length(unique(iCol))})
            if(any(cov.levels <= n.time)){ ## needs to have less levels than time
                possibleTime.vars <- vars.mean_notime[cov.levels <= n.time]

                ls.test.trans <- lapply(possibleTime.vars, function(iVar){ ## iVar <- possibleTime.vars[1]
                    ## table to find the transformation from the time variable to the covariate
                    iTable <- table(data[[iVar]],data$XXtime.indexXX)
                    iOut <- all(colSums(iTable>0)==1) ## no time value leading to two different covariate values
                    if(iOut){ ## explicit the transformation
                        iIndex.level <- apply(iTable>0, MARGIN = 2, FUN = which)
                        if(iVar %in% names(object$xfactor$mean)){
                            attr(iOut,"table") <- object$xfactor$mean[[iVar]][iIndex.level]
                        }else{
                            attr(iOut,"table") <- unique(data[[iVar]])[iIndex.level]
                        }                        
                    }
                    return(iOut)
                })
                index.newtime <- which(unlist(ls.test.trans)==TRUE)
                if(length(index.newtime)>0){
                    ls.out$mean.type[names(ls.out$mean.type) %in% possibleTime.vars[index.newtime]] <- "time"
                    attr(ls.out$mean.type, "table") <- stats::setNames(lapply(ls.test.trans[index.newtime],attr,"table"), possibleTime.vars[index.newtime])
                    vars.mean_notime <- setdiff(vars.mean_notime, possibleTime.vars[index.newtime])
                }
            }

            ## find time varying covariates
            if(length(vars.mean_notime)>0){
                test.timevar <- sapply(vars.mean_notime, function(iVar){
                    any(tapply(data[[iVar]],data$XXcluster.indexXX, function(iVec){sum(!duplicated(iVec))>1}))
                })            
                ls.out$mean.type[names(ls.out$mean.type) %in% names(which(test.timevar))] <- "timevar"
            }
        }
    }

    ## ** export        
    if(simplify && length(effects)==1){
        if(effects=="all"){
            out <- unname(sort(unique(unlist(ls.out))))
            attributes(out) <- c(attributes(out),ls.out[lengths(ls.out)>0])
        }else{
            out <- ls.out[[effects]]
        }
        
    }else{
        if(effects=="all"){
            out <- ls.out
        }else{
            out <- ls.out[effects]
        }
    }
    return(out)
}

## * variable.names.mlmm (documentation)
##' @title Variables Involved in Multiple Linear Mixed Models
##' @description Extract the variables used to obtain multiple linear mixed models.
##'
##' @param object a \code{mlmm} object.
##' @param effects [character] Should all variable be output (\code{"all"}),
##' or only those related to the outcome (\code{"outcome"}), mean (\code{"mean"}), variance (\code{"variance"}),
##' correlation (\code{"correlation"}), time (\code{"time"}), cluster (\code{"cluster"}), strata (\code{"strata"}),
##' or variable defining the split of the dataset on which a separate linear mixed model is fit (\code{"by"})?
##' @param original [logical] Should only the variables present in the original data be output?
##' When \code{FALSE}, variables internally created are output instead of the original variable for time, cluster, and strata.
##' @param simplify [logical] Should the list be converted into a vector if a single \code{effects} is requested?
##' @param ... not used. For compatibility with the generic function
##'
##' @return A list of character vectors or a character vector.
##' 
##' @keywords methods

## * manifest.mlmm 
##' @export
variable.names.mlmm <- function(object, effects = "all", original = TRUE, simplify = TRUE, ...){

    ## ** check user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character vector. \n")
    }
    valid.effects <- c("outcome","mean","variance","correlation","time","cluster","strata",
                       "all","mean.type","by")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }
    if(all("all" %in% effects)){
        if(length(effects)>1){
            stop("Argument \'effects\' must have length 1 when containing the element \"all\". \n")
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

    ## ** extract lmm specific variables
    ls.manifest <- lapply(object$model, variable.names.lmm, effects = setdiff(effects,"by"), original = original, simplify = simplify)
    
    if(simplify){
        if(any(sapply(ls.manifest,identical,ls.manifest[[1]])==FALSE)){
            stop("Difference in manifest variables between the LMM. \n",
                 "Cannot provide a single output for all models.")
        }
        if(any(c("all","by") %in% effects)){
            out <- union(ls.manifest[[1]], object$object$by)
            attributes(out) <- c(attributes(ls.manifest[[1]]), list(by = object$object$by))
        }else{
            out <- ls.manifest[[1]]
        }
    }else{
        out <- ls.manifest
        if(any(c("all","by") %in% effects)){
            attr(out,"by") <- object$object$by
        }
    }

    ## ** export
    return(out)
}


##----------------------------------------------------------------------
### variable.names.R ends here
