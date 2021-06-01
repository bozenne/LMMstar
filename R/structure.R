### structure.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 31 2021 (15:28) 
## Version: 
## Last-Updated: May 31 2021 (15:52) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Unstructured Covariance Matrix
##' @description  Unstructured covariance matrix
##'
##' @param formula Time and cluster variables. The left hand side of the formula is ignored.
##'
##' @details A typical formula would be \code{~time} or \code{~time|id}, indicating a time-specific variance parameter and a correlation parameter specific to each pair of times.
##' 
## @export
UN <- function(formula){
    out <- list(formula = formula,
                type = "UN")
    class(out) <- append("structure",class(out))
    return(out)
}

##' @title Compound Symmetry Covariance Matrix
##' @description  Compound symmetry covariance matrix
##'
##' @param formula Time and cluster variables. The left hand side of the formula is ignored.
##'
##' @details A typical formula would be \code{~time|id}, indicating a variance constant over time and the same correlation between all pairs of times.
##' 
## @export
CS <- function(formula){
    out <- list(formula = formula,
                type = "CS")
    class(out) <- append("structure",class(out))
    return(out)
}
##----------------------------------------------------------------------
### structure.R ends here
