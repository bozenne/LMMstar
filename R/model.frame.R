### model.frame.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (14:57) 
## Version: 
## Last-Updated: Jun  7 2021 (14:58) 
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

## * model.frame.lmm
##' @export
model.frame.lmm <- function(formula, ...){
    return(formula$data)
}

##----------------------------------------------------------------------
### model.frame.R ends here
