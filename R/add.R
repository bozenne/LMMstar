### add.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 23 2022 (16:58) 
## Version: 
## Last-Updated: sep 28 2022 (12:19) 
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

##' @title Add Columns to Output
##' @description Auxiliary function that can be used when specifying the argument \code{columns} (e.g. calling \code{confint.lmm}) to add columns.
##'
##' @param ... [character vector] name of the columns to be added to the default output.
##'
##' @return A character vector
##' 
##' @export
add <- function(...){
    dots <- list(...)
    return(stats::setNames(unlist(dots),rep("add",length(dots))))
}
    

##----------------------------------------------------------------------
### add.R ends here
