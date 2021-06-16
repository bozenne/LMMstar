### utils-formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 23 2021 (09:41) 
## Version: 
## Last-Updated: Jun 14 2021 (13:08) 
##           By: Brice Ozenne
##     Update #: 70
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * rhs.vars
##' @title Extract Variables in the Right Hand Side of a Formula
##' 
##' @param formula a formula
##' 
##' @examples
##' rhs.vars(~X)
##' rhs.vars(Y~X)
##' rhs.vars(Y~X1+X2)
##' rhs.vars(Y~1)
##' rhs.vars(Y+Z~1)
##' rhs.vars(~1|id)
##' rhs.vars(Y~X|id)
##' rhs.vars(Y+Z~X|id)
##' @export
rhs.vars <- function(formula){
    group.var <- all.vars(nlme::getGroupsFormula(formula))
    other.var <- all.vars(nlme::getCovariateFormula(formula))
    return(stats::setNames(c(other.var,group.var),
                           c(rep("covariate",length(other.var)),
                             rep("group",length(group.var)))))
}

## * lhs.vars
##' @title Extract Variables in the Left Hand Side of a Formula
##' 
##' @param formula a formula
##' 
##' @examples
##' lhs.vars(~X)
##' lhs.vars(Y~X)
##' lhs.vars(Y~X1+X2)
##' lhs.vars(Y~1)
##' lhs.vars(Y+Z~1)
##' lhs.vars(~1|id)
##' lhs.vars(Y~X|id)
##' lhs.vars(Y+Z~X|id)
##' @export
lhs.vars <- function(formula){
    setdiff(all.vars(formula),rhs.vars(formula))
}

## * tr
## compute the trace
tr <- function(object){
    sum(diag(object))
}

## * tblock
## transpose by block (depend on the dimension of the matrix)
tblock <- function(M){
    out <- matrix(t(M), nrow=nrow(M), byrow = FALSE)[, c(matrix(1:ncol(M), nrow(M), byrow=T)) ]
    return(out)
}

##----------------------------------------------------------------------
### utils-formula.R ends here
