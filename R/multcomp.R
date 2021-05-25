### multcomp.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 10 2021 (15:57) 
## Version: 
## Last-Updated: May 21 2021 (09:59) 
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

## * estfun.lmm
#' @title Extract the Score Function for Multcomp
#' @description Extract the Score Function for Multcomp. For internal use.
#' @rdname estfun
#' @export
estfun.lmm <- function(x, ...){
    U <- lava::score(x, indiv = TRUE)
    return(U)
}


##----------------------------------------------------------------------
### multcomp.R ends here
