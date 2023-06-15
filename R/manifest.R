### manifest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 31 2022 (15:05) 
## Version: 
## Last-Updated: jun 15 2023 (16:16) 
##           By: Brice Ozenne
##     Update #: 19
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * manifest
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
##' @param ... not used. For compatibility with the generic function
##'
##' @return A character vector
##' 
##' @keywords methods
##' 
##' @export
manifest.lmm <- function(x, effects = "all", original = TRUE, ...){

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

    ## ** check user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(length(effects)==1 && effects == "all"){
        out <- unlist(ls.out)
    }else{
        effects <-  match.arg(effects, valid.effects, several.ok = TRUE)
        out <- unlist(ls.out[effects])
    }
    return(out[!is.na(out)])
}

##----------------------------------------------------------------------
### manifest.R ends here
