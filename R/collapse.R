### collapse.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 18 2024 (12:17) 
## Version: 
## Last-Updated: sep 30 2024 (13:44) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * collapse (documentation)
##' @title Aggregate Dataset by Columns
##' @description Alternative to \code{\link{interaction}} (which can be quite slow) to aggregate columns of a dataset. 
##' @name collapse
##'
##' @param object dataset
##' @param sep [character] character strings used to concatenate values from different columns.
##' @param as.factor [logical] should the result be converted into a factor or be numeric or character.
##' @param ... Not used. For compatibility with the generic method.½ 
##' 
##' @keywords internal
##' 
##' @examples
##' df1 <- data.frame(var1 = 1:10)
##' nlme::collapse(df1)
##' 
##' df2 <- data.frame(var1 = 1:10, var2 = -(1:10))
##' nlme::collapse(df2)
##' interaction(df2, sep = ".", drop = TRUE)
##' 
##' df3 <- data.frame(var1 = 1:1000, var2 = -(1:1000), var3 = (1:1000)/1000)
##' test <- nlme::collapse(df3)
##' GS <- interaction(df3, sep = ".", drop = TRUE)

## * collapse.data.frame (code)
##' @rdname collapse
collapse.data.frame <- function(object, sep = ".", as.factor = TRUE, ...){

    ## case of single column (fast)
    if(length(object)==1){
        if(is.character(as.factor)){
            return(factor(object[[1]], as.factor))
        }else if(as.factor){
            return(as.factor(object[[1]]))
        }else{
            return(as.character(object[[1]]))
        }
    }

    ## case of multiple columns    
    out <- do.call(paste, c(object, sep = unname(sep)))
    if(is.character(as.factor)){
        return(factor(out, as.factor))
    }else if(as.factor){
        return(as.factor(out))
    }else{
        return(out)
    }
}

## * collapse.list (code)
##' @rdname collapse
collapse.list <- collapse.data.frame

## * collapse.matrix (code)
##' @rdname collapse
collapse.matrix <- function(object, ...){
    return(nlme::collapse(as.data.frame(object), ...))
}


##----------------------------------------------------------------------
### collapse.R ends here
