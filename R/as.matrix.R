### as.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2024 (18:56) 
## Version: 
## Last-Updated: jul  5 2024 (18:58) 
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

## * as.matrix.correlate
## ' @export
as.matrix.correlate <- function(x, ...){
    if(length(x)>1 || length(x[[1]])>1){
        message("Only the first correlation matrix was extracted. \n")
    }    
    return(x[[1]][[1]])
}


##----------------------------------------------------------------------
### as.matrix.R ends here
