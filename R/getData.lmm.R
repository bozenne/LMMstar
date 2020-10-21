### getData.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (15:53) 
## Version: 
## Last-Updated: okt 21 2020 (15:54) 
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

##' @export
getData.lmm <- function (object) {
    return(attr(object,"data"))
}


######################################################################
### getData.lmm.R ends here
