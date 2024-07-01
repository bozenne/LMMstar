### collapse.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 18 2024 (12:17) 
## Version: 
## Last-Updated: May 18 2024 (12:18) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * collapse (documentation)
##' @description Aggregate columns.
##' Alternative to interaction which can be quite slow 
##' @noRd
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

## * collapse.list (code)
collapse.list <- collapse.data.frame

## * collapse.matrix (code)
collapse.matrix <- function(value, ...){
    return(nlme::collapse(as.data.frame(value), ...))
}


##----------------------------------------------------------------------
### collapse.R ends here
