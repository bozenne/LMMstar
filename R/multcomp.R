### multcomp.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 10 2021 (15:57) 
## Version: 
## Last-Updated: May 31 2021 (18:01) 
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

## * estfun.lmm
##' @title Extract the Score Function for Multcomp
##' @description Extract the Score Function for Multcomp. For internal use.
##' 
##' @param x  a \code{lmm} object.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @name estfun
##'
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit mixed model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##' 
##' ## test multiple linear hypotheses
##' if(require(multcomp)){
##' multcomp::glht(eUN.lmm)
##' }
#' @method estfun lmm
#' @export
estfun.lmm <- function(x, ...){
    U <- lava::score(x, indiv = TRUE)
    return(U)
}

##----------------------------------------------------------------------
### multcomp.R ends here
