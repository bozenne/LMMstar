### utils-formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 23 2021 (09:41) 
## Version: 
## Last-Updated: mar 11 2024 (11:29) 
##           By: Brice Ozenne
##     Update #: 241
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:



## * addLeading0
##' @description Add leading 0 to a numeric vector
##' @noRd
##' @examples
##' addLeading0(c(1,9,19))
##' addLeading0(c(1,9,100))
addLeading0 <- function(object, as.factor = FALSE, code = NULL){
    if(is.null(code)){
        n.0 <- ceiling(log10(max(abs(object)))+0.1)
        code <- paste0("%0",n.0,"d")
    }
    out <- sprintf(code, object)
    if(as.factor){
        out <- as.factor(out)
    }
    attr(out,"code") <- code
    return(out)
}

## * collapse (collapse.data.frame, collapse.list, collapse.matrix)
##' @description Aggregate columns.
##' Alternative to interaction which can be quite slow 
##' @noRd
##' @examples
##' df1 <- data.frame(var1 = 1:10)
##' collapse(df1)
##' 
##' df2 <- data.frame(var1 = 1:10, var2 = -(1:10))
##' collapse(df2)
##' interaction(df2, sep = ".", drop = TRUE)
##' 
##' df3 <- data.frame(var1 = 1:1000, var2 = -(1:1000), var3 = (1:1000)/1000)
##' test <- collapse(df3)
##' GS <- interaction(df3, sep = ".", drop = TRUE)
collapse.data.frame <- function(value, sep = ".", as.factor = TRUE, ...){

    ## case of single column (fast)
    if(length(value)==1){
        if(is.character(as.factor)){
            return(factor(value[[1]], as.factor))
        }else if(as.factor){
            return(as.factor(value[[1]]))
        }else{
            return(as.character(value[[1]]))
        }
    }

    ## case of multiple columns    
    out <- do.call(paste, c(value, sep = unname(sep)))
    if(is.character(as.factor)){
        return(factor(out, as.factor))
    }else if(as.factor){
        return(as.factor(out))
    }else{
        return(out)
    }
}

collapse.list <- collapse.data.frame
collapse.matrix <- function(value, ...){
    return(nlme::collapse(as.data.frame(value), ...))
}

## * countChar
##' @description Count the number of time a character appears
##'
##' @param value [character] 
##' @param pattern [character]
##' 
##' @noRd
##' @examples
##' countChar("~1|id","|")
##' countChar("~(1|id) + (time|id)","|")
##' countChar("~(1|id) + (time|id)","time")
countChar <- function(value, pattern, fixed = TRUE){
    lengths(regmatches(value, gregexpr(pattern, value, fixed = fixed)))
}


## * is.invertible
##' @description check whether a matrix is invertible
##' 
##' @examples
##' M <- matrix(1:9,3,3)
##' is.invertible(M, cov2cor = FALSE)
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


## * ncharTable
##' @description compute the width of a table
##' 
##' @param object [data.frame or matrix] table
##' @param digits [integer, >0] digits used to round the numeric value
##' @param format [character] how to display numeric value (100, 1e2, ...). See \code{\link{base::formatC}}.
##' @noRd
##'
##' @examples
##' df <- data.frame(name = c("ah","ohhh"), value = c(101.501,200.510))
##' ncharTable(df) ## formatC(1.1, format = "g")
##' ncharTable(df, digits = 7)
##' ncharTable(df, digits = 1) 
##'
##' df2 <- data.frame(name = c("ah","ohhh"), value = c("101.501","200.510"))
##' ncharTable(df2)
##' 
ncharTable <- function(object, digits, format = "f"){
    if(is.matrix(object)){
        if(missing(digits)){
            xplus <- cbind(rownames(object),formatC(object, format = format))
        }else{
            xplus <- cbind(rownames(object),formatC(object, digits = digits, format = format))
        }
    }else if(is.list(object)){
        test.num <- sapply(object,is.numeric)
        object.num <- object[test.num]
        if(missing(digits)){
            object.num <- formatC(as.matrix(object.num), format = format)
        }else{
            object.num <- formatC(as.matrix(object.num), digits = digits, format = format)
        }
        object.char <- object[!test.num]

        xplus <- cbind(rownames(object))
        if(length(object.char)>0){
            xplus <- cbind(xplus, object.char)
        }
        if(length(object.num)>0){
            xplus <- cbind(xplus, object.num)
        }
    }
    nchar.colnames <- c(0,nchar(colnames(object)))
    width <- apply(xplus,1, function(iRow){ ## iRow <- xplus[1,]
        sum(pmax(nchar.colnames,nchar(trimws(iRow)))+1, na.rm = TRUE)-1
    })
    return(max(width))
   
}


## * orderLtoR
##' @description order the lines of a matrix according to the column value from left to right
##' @noRd
##'
##' @examples
##' df <- rbind(data.frame(time = c("a","b","c","d"), gender = "M"),
##'             data.frame(time = c("a","b","c","d"), gender = "F"))
##'
##' X0 <- unique(model.matrix(~1, df))
##' X1 <- unique(model.matrix(~time, df))
##' X2 <- .model.matrix_regularize(~0+gender + time:gender, data = df, augmodel = TRUE)
##'
##' orderLtoR(X0)
##' orderLtoR(X1)
##' orderLtoR(X2, strata = attr(X2,"M.level")$gender)
##' orderLtoR(X2[8:1,], strata = attr(X2,"M.level")$gender)
orderLtoR <- function(object, strata = NULL){
    if(is.null(strata)){
        new.order <- rev(1:NCOL(object))
    }else{
        new.order <- unlist(rev(tapply(colnames(object),strata,rev, simplify = FALSE)))
    }
    out <- as.numeric(nlme::collapse(object[,new.order,drop=FALSE], sep=""))
    ## object[order(out),]
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

## * tr
##' @description compute the trace
##' @noRd
##'
##' @examples
##' M <- matrix(1:9,3,3)
##' tr(M)
tr <- function(object){
    sum(diag(object))
}

## * unorderedPairs
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
##' unorderedPairs(1, distinct = FALSE)
##' unorderedPairs(1, distinct = TRUE)
##' 
##' unorderedPairs(1:5, distinct = TRUE) - utils::combn(5, m = 2)
##' unorderedPairs(1:5, distinct = FALSE)
##' unorderedPairs(rep(1,5), distinct = TRUE) - (utils::combn(5, m = 2)>0)
##' 
unorderedPairs <- function(x, distinct = FALSE){
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


## * groupSet
##' @title Group Sets
##' @description Find groups of sets sharing the same values
##' @noRd
##'
##' @param object a list
##'
##'
##' @examples
##' mylistA <- list(1:2,3:4,1:5)
##' groupSet(mylistA)
##'
##' mylistB <- list(1:2,2:3,2,3)
##' groupSet(mylistB)
##' 
##' mylistC <- list(1:2,1,2,2,3,2:3)
##' groupSet(mylistC, strata = c(1,1,1,2,2,2))
##' 
##' mylistD <- list(1:2,1,2,2,3,2:3,4,4:5,6)
##' groupSet(mylistD, strata = c(1,1,1,2,2,2,2,2,3))
groupSet <- function(object, strata = NULL){

    n.set <- length(object)
    out <- rep(NA, n.set)

    ## ** handle strata
    if(!is.null(strata)){
        previous.group <- 0
        U.strata <- unique(strata)
        n.strata <- length(U.strata)
        index.strata <- tapply(1:length(object), strata, identity)
        for(iS in 1:n.strata){ ## iS <- 1
            iOut <- groupSet(object[index.strata[[iS]]])
            out[index.strata[[iS]]] <- iOut + previous.group
            previous.group <- previous.group + max(iOut)
        }
        return(out)
    }

    ## ** initialize with the largest element
    order.set <- order(lengths(object), decreasing = TRUE)
    out[order.set[1]] <- 1
    
    ## ** deal with special case
    if(n.set==1){
        return(out)
    }

    ## ** find representative patterns
    Uvalue <- unique(unlist(object))

    M.group <- rbind(Uvalue %in% object[[order.set[1]]])
    index.group <- order.set[1]
    if(all(M.group)){ ## another special case
        out[] <- 1
        return(out)
    }

    for(iSet in order.set[-1]){ ## iSet <- 2
        iTest <- Uvalue %in% object[[iSet]]
        iContrast <- sweep(M.group, MARGIN = 2, FUN = "-", STATS = iTest)

        if(any(rowSums(iContrast!=0)==0)){ ## since patterns are looped over with decreasing length, equality can only mean pattern with fewer elements
            out[iSet] <- out[index.group[which(rowSums(iContrast!=0)==0)[1]]]
        }else if(any(rowSums(iContrast<0)==0)){ ## subset of another pattern (possibly several but label as first one)
            out[iSet] <- out[index.group[which(rowSums(iContrast<0)==0)[1]]]
        }else if(any(rowSums(iContrast>0)==0)){ ## contain other patterns 
            iMerge <- which(rowSums(iContrast>0)==0)

            out[c(which(out %in% iMerge),iSet)] <- out[index.group[iMerge[1]]]
            index.group[iMerge[1]] <- iSet
            M.group[iMerge[1],] <- iTest
            if(length(iMerge)>1){ ## handle the case where merging two covaraince structures into a bigger one, e.g. A,B and A,C into A,B,C
                index.group <- index.group[-iMerge[-1]]
                M.group <- M.group[-iMerge[-1],,drop=FALSE]
                out <- as.numeric(as.factor(out))
            }
        }else{ ## new pattern
            out[iSet] <- length(index.group)+1
            index.group <- c(index.group, iSet)
            M.group <- rbind(M.group, iTest)
        }
    }

    ## ** export
    return(out)
    
}

##----------------------------------------------------------------------
### utils-formula.R ends here
