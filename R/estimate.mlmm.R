### estimate.mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2022 (15:42) 
## Version: 
## Last-Updated: jun 22 2022 (15:45) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @export
estimate.mlmm <- function(x, ...){

    return(attr(x$all,"glht")[[1]]$model)

}


##----------------------------------------------------------------------
### estimate.mlmm.R ends here
