### as.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2024 (18:56) 
## Version: 
## Last-Updated: jul 10 2024 (15:31) 
##           By: Brice Ozenne
##     Update #: 14
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
as.matrix.correlate <- function(x, index, ...){

    ## ** check and normalize user imput
    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown arguments \'",paste(names(dots), collapse = "\' \'"),"\'. \n")
    }

    ## ** subset
    if(!missing(index)){
        x.array <- as.array(x)
        if(length(index)>1){
            stop("Incorrect argument \'index\': should have length 1. \n")
        }else if(is.numeric(index) && index %in% 1:dim(x.array)[3] == FALSE){
            stop("Incorrect argument \'index\': when numeric it should take integer value between 1 and the number of correlation matrix (here ",dim(x.array)[3],"). \n")
        }else if(is.character(index) && index %in% dimnames(x.array)[[3]] == FALSE){
            stop("Incorrect argument \'index\': when character it should refer to the name of the correlation matrix. \n",
                 "Possible names: \"",paste(dimnames(x.array)[[3]], collapse = "\", \""),"\". \n")
        }
        out <- x[,,index]
    }else{
        if(length(x)>1 || length(x[[1]])>1){
            message("Only the first correlation matrix was extracted. \n")
        }
        out <- x[[1]][[1]]
    }

    ## ** export
    return(out)
}


##----------------------------------------------------------------------
### as.matrix.R ends here
