### remove.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 23 2022 (16:59) 
## Version: 
## Last-Updated: sep 23 2022 (17:00) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

remove <- function(...){
    dots <- list(...)
    return(stats::setNames(unlist(dots),rep("remove",length(dots))))
}

##----------------------------------------------------------------------
### remove.R ends here
