### manifest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 31 2022 (15:05) 
## Version: 
## Last-Updated: jul 26 2023 (11:07) 
##           By: Brice Ozenne
##     Update #: 38
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
##' When \code{NULL}, variables internally created to fill absent variables will be added to the output.
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
    if(!is.null(original) && original){
        ls.out <- list(outcome = x$outcome$var,
                       mean = attr(x$design$mean, "variable"),
                       variance = all.vars(x$design$vcov$formula$var),
                       correlation = all.vars(x$design$vcov$formula$cor),
                       time = attr(x$time$var,"original"),
                       cluster = attr(x$cluster$var,"original"),
                       strata = attr(x$strata$var,"original"))
    }else{
        ls.out <- list(outcome = x$outcome$var,
                       mean = attr(x$design$mean, "variable"),
                       variance = all.vars(x$design$vcov$formula$var),
                       correlation = all.vars(x$design$vcov$formula$cor),
                       time = x$time$var,
                       cluster = x$cluster$var,
                       strata = x$strata$var)
        if(is.null(original)){
            if(is.na(ls.out$time)){  ls.out$time <- x$time$var  }
            if(is.na(ls.out$cluster)){  ls.out$cluster <- x$cluster$var  }
            if(is.na(ls.out$strata)){  ls.out$strata <- x$strata$var  }
        }
    }
    valid.effects <- names(ls.out)
    ## normalize numeric(0) and NA into NULL
    ls.out[sapply(ls.out, function(iE){sum(!is.na(iE))==0})] <- list(NULL)

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
            out <- unname(sort(unique(unlist(ls.out[effects]))))
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
