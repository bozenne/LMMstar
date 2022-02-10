### print.anova_lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 24 2022 (10:05) 
## Version: 
## Last-Updated: feb 10 2022 (13:18) 
##           By: Brice Ozenne
##     Update #: 27
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
    dots$print.null <- NULL
    return(do.call(summary, c(list(object = x, print.nulls = FALSE), dots)))
}

##----------------------------------------------------------------------
### print.anova_lmm.R ends here
