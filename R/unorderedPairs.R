### unorderedPairs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:07) 
## Version: 
## Last-Updated: apr 13 2022 (10:10) 
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

## ** .unorderedPairs
##' @title Form All Pairs
##' @description Form all pairs of values
##' @noRd
##'
##' @param x vector of values
##' @param distinct [logical] should pairs containing the same value be removed?
##' 
##' @details adapted from RecordLinkage package
##'
##' @examples
##' .unorderedPairs(1:5, distinct = TRUE)
##' .unorderedPairs(1:5, distinct = FALSE)
##' 
.unorderedPairs <- function(x, distinct = FALSE){
    n <- length(x)
    ls <- lapply(1:n, function(k){ rbind(x[k], x[k:n])})
    out <- do.call(cbind,ls)##array(unlist(ls), dim = c(2, n * (n + 1)/2))
    if(distinct){
        out <- out[,apply(out,2,function(iCol){all(!duplicated(iCol))}),drop=FALSE]
    }
    return(out)
}


##----------------------------------------------------------------------
### unorderedPairs.R ends here
