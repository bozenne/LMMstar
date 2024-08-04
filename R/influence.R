### influence.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  4 2024 (11:53) 
## Version: 
## Last-Updated: Aug  4 2024 (14:33) 
##           By: Brice Ozenne
##     Update #: 4
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

## * influence.mlmm (code)
##' @export
influence.mlmm <- function(model, ...){
    iid.mlmm(x = model, ...)
}

## * influence.rbindWald_lmm (code)
##' @export
influence.rbindWald_lmm <- function(model, ...){
    iid.rbindWald_lmm(x = model, ...)
}

## * influence.Wald_lmm (code)
##' @export
influence.Wald_lmm <- function(model, ...){
    iid.Wald_lmm(x = model, ...)
}

##----------------------------------------------------------------------
### influence.R ends here
