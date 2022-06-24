### print.mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2022 (15:11) 
## Version: 
## Last-Updated: jun 22 2022 (15:18) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.mlmm
##' @export
print.mlmm <- function(x, ...){

    print(attr(x$all,"glht")[[1]]$model)

    return(invisible(NULL))
}



##----------------------------------------------------------------------
### print.mlmm.R ends here
