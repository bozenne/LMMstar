### unorderedPairs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 13 2022 (10:07) 
## Version: 
## Last-Updated: jun 27 2022 (12:16) 
##           By: Brice Ozenne
##     Update #: 17
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
##' .unorderedPairs(1:5, distinct = TRUE) - utils::combn(5, m = 2)
##' .unorderedPairs(1:5, distinct = FALSE)
##' 
.unorderedPairs <- function(x, distinct = FALSE){
    n.x <- length(x)
    out <- do.call(cbind,lapply(1:n.x, function(iK) {
        rbind(x[iK], x[iK:n.x])
    }))
    
    if(distinct){## same as combn but faster when x is large
        return(out[,out[1,]!=out[2,],drop=FALSE])
    }else{
        return(out)
    }
}


##----------------------------------------------------------------------
### unorderedPairs.R ends here
