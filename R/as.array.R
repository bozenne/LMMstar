### as.array.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  9 2024 (10:48) 
## Version: 
## Last-Updated: jul 31 2024 (11:28) 
##           By: Brice Ozenne
##     Update #: 9
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * as.array.correlate
##' @export
as.array.correlate <- function(x, ...){

    ## ** check and normalize user imput
    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown arguments \'",paste(names(dots), collapse = "\' \'"),"\'. \n")
    }

    ## ** convert to array
    if(length(x)==1){
        ls.x <- x[[1]]
    }else{
        ls.x <- unlist(x, recursive = FALSE)
    }
    out <- array(as.numeric(unlist(ls.x)), dim=c(dim(ls.x[[1]]),length(ls.x)),
                 dimnames = c(dimnames(ls.x[[1]]), list(names(ls.x))))

    ## ** export
    return(out)
}


##----------------------------------------------------------------------
### as.array.R ends here
