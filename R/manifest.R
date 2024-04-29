### manifest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 31 2022 (15:05) 
## Version: 
## Last-Updated: Mar 27 2024 (12:22) 
##           By: Brice Ozenne
##     Update #: 93
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * manifest (documentation)
##' @title Variables Involved in a Linear Mixed Model
##' @description Extract the variables used by the linear mixed model.
##'
##' @param x a \code{lmm} object.
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

## * manifest.lmm 
##' @export
manifest.lmm <- function(x, effects = "all", original = TRUE, simplify = TRUE, ...){

    ## ** extract variables
    var.outcome <- x$outcome$var
    var.cluster <- x$cluster$var
    var.time <- x$time$var
    var.strata <- x$strata$var
    n.time <- x$time$n

    ls.out <- list(outcome = var.outcome,
                   mean = attr(x$design$mean, "variable"),
                   variance = all.vars(x$design$vcov$formula$var),
                   correlation = all.vars(x$design$vcov$formula$cor),
                   time = var.time,
                   cluster = var.cluster,
                   strata = var.strata)
    if(original){
        ls.out$time <- attr(ls.out$time,"original")
        ls.out$cluster <- attr(ls.out$cluster,"original")
        ls.out$strata <- attr(ls.out$strata,"original")
    }
    valid.effects <- c(names(ls.out),"mean.type")
    ## normalize numeric(0) and NA into NULL
    ls.out[sapply(ls.out, function(iE){sum(!is.na(iE))==0})] <- list(NULL)

    ## ** time-invariant variables
    if("mean.type" %in% effects){
        ## hidden: no returned by effects = "all" to save computation time
        ## only internal function need effects="mean.type" (predict,model.frame) for reshaping data
        n.time <- x$time$n
        vars.time <- stats::na.omit(unique(c(var.time,attr(var.time,"original"))))
        vars.cluster <- stats::na.omit(unique(c(var.cluster,attr(var.cluster,"original"))))
        vars.mean_notime <- ls.out$mean
        data <- x$data

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
                        if(iVar %in% names(x$xfactor$mean)){
                            attr(iOut,"table") <- x$xfactor$mean[[iVar]][iIndex.level]
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
    
    ## ** check user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    effects <-  match.arg(effects, c("all",valid.effects), several.ok = TRUE)
        
    if(simplify && length(effects)==1){
        if(effects=="all"){
            out <- unname(sort(unique(unlist(ls.out))))
            attributes(out) <- c(attributes(out),ls.out[lengths(ls.out)>0])
        }else{
            out <- ls.out[[effects]]
        }
        
    }else{
        out <- ls.out[effects]
    }
    return(out)
}

## * manifest.mlmm 
##' @export
manifest.mlmm <- function(x, ...){
    ls.manifest <- lapply(x$model, manifest)

    if(any(sapply(ls.manifest,identical,ls.manifest[[1]])==FALSE)){
        stop("Difference in manifest variables between the LMM. \n",
             "Cannot provide a single output for all models.")
    }
    out <- union(ls.manifest[[1]], x$object$by)
    attributes(out) <- c(attributes(ls.manifest[[1]]), list(by = x$object$by))
    return(out)
}


##----------------------------------------------------------------------
### manifest.R ends here
