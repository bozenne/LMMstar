### manifest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 31 2022 (15:05) 
## Version: 
## Last-Updated: okt 31 2022 (15:29) 
##           By: Brice Ozenne
##     Update #: 14
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
##' @param object a \code{lmm} object.
##' @param effects [character] Should all variable be output (\code{"all"}),
##' or only those related to the outcome (\code{"outcome"}), mean (\code{"mean"}), variance (\code{"variance"}),
##' correlation (\code{"correlation"}), time (\code{"time"}), cluster (\code{"cluster"}), strata (\code{"strata"})?
##' @param original [logical] Should only the variables present in the original data be output?
##' When \code{NULL}, variables internally created to fill absent variables will be added to the output.
##' When \code{FALSE}, variables internally created are output instead of the original variable for time, cluster, and strata.
##'
##' @return A character vector
##' @export
manifest.lmm <- function(object, effects = "all", original = TRUE){


    ## ** extract variables
    if(!is.null(original) && original){
        ls.out <- list(outcome = object$outcome$var,
                       mean = attr(object$design$mean, "variable"),
                       variance = all.vars(object$design$vcov$formula$var),
                       correlation = all.vars(object$design$vcov$formula$cor),
                       time = attr(object$time$var,"original"),
                       cluster = attr(object$cluster$var,"original"),
                       strata = attr(object$strata$var,"original"))
    }else{
        ls.out <- list(outcome = object$outcome$var,
                       mean = attr(object$design$mean, "variable"),
                       variance = all.vars(object$design$vcov$formula$var),
                       correlation = all.vars(object$design$vcov$formula$cor),
                       time = object$time$var,
                       cluster = object$cluster$var,
                       strata = object$strata$var)
        if(is.null(original)){
            if(is.na(ls.out$time)){  ls.out$time <- object$time$var  }
            if(is.na(ls.out$cluster)){  ls.out$cluster <- object$cluster$var  }
            if(is.na(ls.out$strata)){  ls.out$strata <- object$strata$var  }
        }
    }
    valid.effects <- names(ls.out)

    ## ** check user input
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
