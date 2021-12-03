### terms.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  3 2021 (10:05) 
## Version: 
## Last-Updated: dec  3 2021 (10:59) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## used by mu
##' @title Model Terms For Linear Mixed Models
##' @description Model terms for linear mixed models. Used by \code{multcomp::glht}.
##' 
##' @param x a \code{lmm} object
##' @param ... not used, for compatibility with the generic method.
##'
##' @return An object of class \code{terms} giving a symbolic representation of the mean structure.
##' 
##' @export
terms.lmm <- function(x, ...){
    return(attr(x$design$mean,"terms"))
}

##----------------------------------------------------------------------
### terms.R ends here
