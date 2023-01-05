### add.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 23 2022 (16:58) 
## Version: 
## Last-Updated: jan  3 2023 (14:54) 
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

##' @title Add Columns to Output
##' @description Auxiliary function that can be used when specifying the argument \code{columns} (e.g. calling \code{confint.lmm}) to add columns.
##'
##' @param ... [character vector] name of the columns to be added to the default output.
##'
##' @return A character vector
##' 
##' @examples
##' set.seed(10)
##' dW <- sampleRem(25, n.times = 1, format = "long")
##' e.lmm <- lmm(Y~X1, data = dW)
##'
##' confint(e.lmm, columns = add("statistic"))

##' @export
add <- function(...){
    dots <- list(...)
    return(stats::setNames(unlist(dots),rep("add",length(dots))))
}
    

##----------------------------------------------------------------------
### add.R ends here
