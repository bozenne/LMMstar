### utils-formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 23 2021 (09:41) 
## Version: 
## Last-Updated: jul  5 2023 (16:12) 
##           By: Brice Ozenne
##     Update #: 127
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * tr
##' @description compute the trace
##' @noRd
tr <- function(object){
    sum(diag(object))
}

## * tblock
##' @description transpose by block (depend on the dimension of the matrix)
##' @noRd
tblock <- function(M){
    if(NROW(M)==1){
        return(M)
    }else{
        out <- matrix(t(M), nrow=nrow(M), byrow = FALSE)[, c(matrix(1:ncol(M), nrow(M), byrow=T)) ]
        return(out)
    }
}

## * ncharTable
##' @description compute the width of a table 
##' @noRd
ncharTable <- function(object, digits){

    xplus <- cbind(rownames(object),formatC(as.matrix(object), digits = digits, format = "f"))
    nchar.colnames <- c(0,nchar(colnames(object)))
    width <- apply(xplus,1, function(iRow){ ## iRow <- xplus[2,]
        sum(pmax(nchar.colnames,nchar(trimws(iRow)))+1)-1
    })
    return(max(width))
   
}

## * is.invertible
##' @description check whether a matrix is invertibel
##' @noRd
is.invertible <- function(object, cov2cor, tol = 10^(-10*sqrt(NCOL(object)))){
    if(any(is.na(object))){
        return(FALSE)
    }
    if(cov2cor){
        if(any(diag(object)<=0)){
            return(FALSE)
        }else{
            object <- stats::cov2cor(object)
            if(any(is.na(object))){
                return(FALSE)
            }
        }
    }
    return(abs(det(object))>tol)
}

## * countChar
##' @description Count the number of time a character appears
##' @noRd
##' @examples
##' countChar("~1|id","|")
##' countChar("~(1|id) + (time|id)","|")
##' countChar("~(1|id) + (time|id)","time")
countChar <- function(value, pattern, fixed = TRUE){
    lengths(regmatches(value, gregexpr(pattern, value, fixed = fixed)))
}


## * .unorderedPairs
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
##' .unorderedPairs(rep(1,5), distinct = TRUE) - (utils::combn(5, m = 2)>0)
##' 
.unorderedPairs <- function(x, distinct = FALSE){
    n.x <- length(x)
    ## work on integers
    y <- 1:n.x
    out <- do.call(cbind,lapply(1:n.x, function(iK) {
        rbind(y[iK], y[iK:n.x])
    }))
    
    ## remove 'diagonal' pairs (e.g. (1,1) or (2,2))
    if(distinct){## same as combn but faster when x is large
        if(all(out[1,]==out[2,])){
            return(NULL)
        }else{
            out <- out[,out[1,]!=out[2,],drop=FALSE]
        }
    }

    ## restaure original values
    out[] <- x[as.vector(out)]
    return(out)
}

## * triplicated
##' @description identify third (or more) occurence of a value in a vector.
##' @noRd
##' @examples
##' triplicated(c(1:3,1:3))
##' triplicated(c(1:3,1:3,1:3))
##' triplicated(c(NA, 1:3, 3, 4:6, 3, NA, 4))
##' triplicated(c(NA, 1:3, 3, 4:6, 3, NA, 4, 3:4))
triplicated <- function (x){
    out <- rep(FALSE,length(x))
    index.duplicated <- which(base::duplicated(x))
    if(length(index.duplicated)==0){return(out)}
    x.duplicated <- x[index.duplicated]
    out[index.duplicated[base::duplicated(x.duplicated)]] <- TRUE
    return(out)
}
##----------------------------------------------------------------------
### utils-formula.R ends here
