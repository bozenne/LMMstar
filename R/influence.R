### influence.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  4 2024 (11:53) 
## Version: 
## Last-Updated: jul  4 2024 (11:54) 
##           By: Brice Ozenne
##     Update #: 2
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
##----------------------------------------------------------------------
### influence.R ends here
