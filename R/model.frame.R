### model.frame.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (14:57) 
## Version: 
## Last-Updated: apr 14 2023 (17:06) 
##           By: Brice Ozenne
##     Update #: 8
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
    data <- formula$data
    rownames(data) <- NULL
    return(data)
}

## * model.frame.lmmCC
##' @export
model.frame.lmmCC <- function(formula, ...){
    data <- formula$data
    keep.var <- unique(manifest(formula, original = TRUE))
    rownames(data) <- NULL
    return(data[keep.var])
}

##----------------------------------------------------------------------
### model.frame.R ends here
