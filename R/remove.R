### remove.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 23 2022 (16:59) 
## Version: 
## Last-Updated: jun 15 2023 (16:24) 
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

##' @title Remove Columns from Output
##' @description Auxiliary function that can be used when specifying the argument \code{columns} (e.g. calling \code{confint.lmm}) to remove columns.
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
##' confint(e.lmm, columns = remove("estimate"))
##' 
##' @export
remove <- function(...){
    dots <- list(...)
    return(stats::setNames(unlist(dots),rep("remove",length(dots))))
}

##----------------------------------------------------------------------
### remove.R ends here
