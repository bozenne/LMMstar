### multcomp.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 10 2021 (15:57) 
## Version: 
## Last-Updated: May  5 2024 (20:14) 
##           By: Brice Ozenne
##     Update #: 29
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
##' @return A matrix containing the score function for each model parameter (columns) relative to each cluster (rows).
##' 
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##' 
##' ## test multiple linear hypotheses
##' if(require(multcomp)){
##' LMMstar.options(effects = c("mean"))
##' e.glht <- multcomp::glht(eUN.lmm)
##' e.glht$linfct
##' }
##' 
#' @method estfun lmm
#' 
#' @keywords interface
estfun.lmm <- function(x, ...){
    U <- lava::score(x, indiv = TRUE)
    return(U)
}

##----------------------------------------------------------------------
### multcomp.R ends here
