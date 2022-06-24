### coef.mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2022 (15:38) 
## Version: 
## Last-Updated: jun 22 2022 (15:44) 
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

## * coef.mlmm
##' @export
coef.mlmm <- function(object, ...){

    ls.model <- lava::estimate(object)
    return(lapply(ls.model,coef, ...))

}

##----------------------------------------------------------------------
### coef.mlmm.R ends here
