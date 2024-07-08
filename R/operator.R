### operator.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2023 (15:05) 
## Version: 
## Last-Updated: jul  5 2024 (19:03) 
##           By: Brice Ozenne
##     Update #: 30
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * [.summarize
##' @export
`[.summarize` <- function(x, i, j, ...){
    x.df <- as.data.frame(x)
    if(!missing(j)){
        if(is.character(j) && any(j %in% names(x.df) == FALSE)){
            stop("Incorrect column argument: cannot find \"",paste(j[j %in% names(x.df) == FALSE], collapse = "\" \""),"\" among the summary statistics. \n",
                 "Available summary statistics: \"",paste(setdiff(names(x.df),j), collapse = "\" \""),"\". \n")
        }
    }
    y <- x.df[i,j,...]

    return(y)
}

## * [.correlate
##' @export
`[.correlate` <- function(x, i, j, k, ...){

    flat.x <- unlist(x, recursive = FALSE)

    ## ** find correlation matrix
    if(missing(k)){
        k <- 1
    }else if(is.character(k) && any(k %in% names(flat.x) == FALSE)){
        stop("Incorrect slice argument: cannot find \"",paste(k[k %in% names(flat.x) == FALSE], collapse = "\" \""),"\" among the outcome.covariate levels. \n",
             "Slice names: \"",paste(setdiff(names(flat.x),k), collapse = "\" \""),"\". \n")
    }
    
    x.mat <- flat.x[[k]]

    if(!missing(j)){
        if(is.character(j) && any(j %in% names(x.mat) == FALSE)){
            stop("Incorrect column argument: cannot find \"",paste(j[j %in% colnames(x.mat) == FALSE], collapse = "\" \""),"\" among the columns. \n",
                 "Column names: \"",paste(setdiff(colnames(x.mat),j), collapse = "\" \""),"\". \n")
        }
    }
    if(!missing(i)){
        if(is.character(i) && any(i %in% names(x.mat) == FALSE)){
            stop("Incorrect row argument: cannot find \"",paste(i[i %in% names(x.mat) == FALSE], collapse = "\" \""),"\" among the rows. \n",
                 "Row names: \"",paste(setdiff(names(x.mat),i), collapse = "\" \""),"\". \n")
        }
    }
    y <- x.mat[i,j,...]

    return(y)
}

##----------------------------------------------------------------------
### operator.R ends here
