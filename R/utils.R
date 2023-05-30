### utils-formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 23 2021 (09:41) 
## Version: 
## Last-Updated: maj 30 2023 (18:22) 
##           By: Brice Ozenne
##     Update #: 119
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * tr
## compute the trace
tr <- function(object){
    sum(diag(object))
}

## * tblock
## transpose by block (depend on the dimension of the matrix)
tblock <- function(M){
    if(NROW(M)==1){
        return(M)
    }else{
        out <- matrix(t(M), nrow=nrow(M), byrow = FALSE)[, c(matrix(1:ncol(M), nrow(M), byrow=T)) ]
        return(out)
    }
}

## * ncharTable
## compute the width of a table 
ncharTable <- function(object, digits){

    xplus <- cbind(rownames(object),formatC(as.matrix(object), digits = digits, format = "f"))
    nchar.colnames <- c(0,nchar(colnames(object)))
    width <- apply(xplus,1, function(iRow){ ## iRow <- xplus[2,]
        sum(pmax(nchar.colnames,nchar(trimws(iRow)))+1)-1
    })
    return(max(width))
   
}

## * is.invertible
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
## countChar("~1|id","|")
## countChar("~(1|id) + (time|id)","|")
## countChar("~(1|id) + (time|id)","time")
countChar <- function(value, pattern, fixed = TRUE){
    lengths(regmatches(value, gregexpr(pattern, value, fixed = fixed)))
}


## * updateFormula
## NOTE: when updating formula not using stats::drop.terms or stats::udpate as it re-write the interaction X1*X2 ---> X1 + X2 + X1:X2
## stats::drop.terms(terms(Y~X1*X2+(1|id)), dropx = 3)
## stats::update(terms(Y~X1*X2+(1|id)), .~.-(1|id))
## updateFormula(Y~X1*X2+(1|id), drop.x = "(1|id)")
## updateFormula(Y~0+X1*X2+(1|id), drop.x = "(1|id)")
##
## updateFormula(Y~X1*X2, add.x = "(1 | id)")
updateFormula <- function(formula, add.x = NULL, drop.x = NULL){

    drop.x <- gsub(" ","",drop.x)
    if(!inherits(formula,"formula") || (length(formula) %in% 2:3 == FALSE)){
        stop("Argument \'formula\' should be a formula with length 2 or 3. \n")
    }
    test.response <- length(formula) == 3
    
    txt.formula <- as.character(utils::tail(formula,1))
    term.formula <- gsub(" ","",strsplit(txt.formula, split = "+", fixed = TRUE)[[1]])
    if(any(drop.x %in% term.formula == FALSE)){
        stop("Mismatch between argument \'formula\' and \'drop.x\', \n",
             "Could not find \"",paste(drop.x[drop.x %in% term.formula == FALSE], collapse = "\" \""),"\". \n")
    }

    if(!is.null(drop.x)){
        term.formula <- term.formula[term.formula %in% drop.x == FALSE]
    }
    if(!is.null(add.x)){
        term.formula <- c(term.formula, add.x)
    }
    if(test.response){
        txt.new <- paste0(deparse(formula[[2]]),deparse(formula[[1]]),paste0(term.formula, collapse="+"))
    }else{
        txt.new <- paste0(deparse(formula[[1]]),paste0(term.formula, collapse="+"))
    }
    return(stats::as.formula(txt.new))
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

##----------------------------------------------------------------------
### utils-formula.R ends here
