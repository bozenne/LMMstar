### influence.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  4 2024 (11:53) 
## Version: 
## Last-Updated: jul 11 2024 (11:55) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * influence.lmm (code)
##' @export
influence.lmm <- function(model, ...){
    iid.lmm(x = model, ...)
}

## * influence.Wald_lmm (code)
##' @export
influence.Wald_lmm <- function(model, ...){
    iid.Wald_lmm(x = model, ...)
}

## * influence.mlmm (code)
##' @export
influence.mlmm <- function(model, ...){
    iid.mlmm(x = model, ...)
}
##----------------------------------------------------------------------
### influence.R ends here
