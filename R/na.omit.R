### na.omit.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 16 2026 (15:09) 
## Version: 
## Last-Updated: okt 16 2026 (15:10) 
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

##' @description Extract the index of rows removed from the original dataset
##' For internal use.
##' @noRd
##' @keywords internal

na.omit.lmm <- function(object, ...){
    return(object$na.omit)
}

##----------------------------------------------------------------------
### na.omit.R ends here
