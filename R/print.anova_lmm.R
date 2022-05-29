### print.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 24 2022 (10:05) 
## Version: 
## Last-Updated: May 29 2022 (11:14) 
##           By: Brice Ozenne
##     Update #: 30
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.anova_lmm
##' @export
print.anova_lmm <- function(x, ...){
    dots <- list(...)
    dots$print <- c(1,0)
    return(do.call(summary, c(list(object = x), dots)))
}

##----------------------------------------------------------------------
### print.anova_lmm.R ends here
