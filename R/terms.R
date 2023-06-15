### terms.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  3 2021 (10:05) 
## Version: 
## Last-Updated: jun 15 2023 (16:57) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## used by glht
##' @title Model Terms For Linear Mixed Models
##' @description Model terms for linear mixed models. Used by \code{multcomp::glht}.
##' 
##' @param x a \code{lmm} object
##' @param ... not used, for compatibility with the generic method.
##'
##' @return An object of class \code{terms} giving a symbolic representation of the mean structure.
##'
##' @keywords interface
##' 
##' @export
terms.lmm <- function(x, ...){
    return(attr(x$design$mean,"terms"))
}

##----------------------------------------------------------------------
### terms.R ends here
