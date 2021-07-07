### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: jul  7 2021 (17:19) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * IND (independence)

## * CS (compound symmetry)
##' @title Compound Symmetry Covariance Matrix
##' @description  Compound symmetry covariance matrix
##'
##' @param formula Time and cluster variables. The left hand side of the formula is ignored.
##'
##' @details A typical formula would be \code{~time|id}, indicating a variance constant over time and the same correlation between all pairs of times.
##' 
##' @export
CS <- function(formula){
    if(missing(formula)){formula <- NULL}
    out <- list(formula = formula,
                type = "CS")
    class(out) <- append("structure",class(out))
    return(out)
}

## * UN (unstructured)
##' @title Unstructured Covariance Matrix
##' @description  Unstructured covariance matrix
##'
##' @param formula Time and cluster variables. The left hand side of the formula is ignored.
##'
##' @details A typical formula would be \code{~time} or \code{~time|id}, indicating a time-specific variance parameter and a correlation parameter specific to each pair of times.
##' 
##' @export
UN <- function(formula){
    out <- list(formula = formula,
                type = "UN")
    class(out) <- append("structure",class(out))
    return(out)
}

## * EXP (exponential)

##----------------------------------------------------------------------
### structure.R ends here
