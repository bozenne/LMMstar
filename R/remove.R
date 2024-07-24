### remove.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 23 2022 (16:59) 
## Version: 
## Last-Updated: jul 24 2024 (11:39) 
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

##' @title Remove Columns from Output
##' @description Auxiliary function that can be used when specifying the argument \code{columns} (e.g. calling \code{confint.lmm}) to remove columns.
##' Not called remove to avoid confusion with \code{base::remove}
##'
##' @param ... [character vector] name of the columns to be removed to the default output.
##'
##' @return A character vector
##' 
##' @keywords utilities
##' 
##' @examples
##' set.seed(10)
##' dW <- sampleRem(25, n.times = 1, format = "long")
##' e.lmm <- lmm(Y~X1, data = dW)
##'
##' confint(e.lmm, columns = rem("estimate"))
##' 
##' @export
rem <- function(...){
    dots <- list(...)
    return(stats::setNames(unlist(dots),rep("remove",length(unlist(dots)))))
}

##----------------------------------------------------------------------
### remove.R ends here
