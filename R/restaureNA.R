### restaureNA.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 20 2023 (15:31) 
## Version: 
## Last-Updated: feb 15 2024 (15:54) 
##           By: Brice Ozenne
##     Update #: 46
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * restaureNA (documentation)
##' @description Restaure NA in the user output that have been excluded when fitting the LMM.
##' @noRd
`restaureNA` <-
  function(object, index.na, level, cluster) UseMethod("restaureNA")

## * restaureNA.vector
restaureNA.numeric <- function(object, index.na, level, cluster){

    if(length(object)==0){
        return(object)
    }else if(level == "cluster" && any(is.na(cluster$index))){
        ## expand to store NAs
        out <- stats::setNames(rep(NA, times = cluster$n), cluster$levels)
        out[names(object)] <- unname(object)
        ## put back attributes
        restaure.attributes <- setdiff(names(attributes(object)), names(attributes(out)))
        if(length(restaure.attributes)>0){
            attributes(out) <- c(attributes(out), attributes(object)[restaure.attributes])
        }

        return(out)
    }else if(level == "obs" && length(index.na)>0){
        ## expand to store NAs
        nobs <- NROW(object) + length(index.na)
        out <- rep(NA, times = nobs)
        out[-index.na] <- unname(object)
        ## put back attributes
        restaure.attributes <- setdiff(names(attributes(object)), names(attributes(out)))
        if(length(restaure.attributes)>0){
            attributes(out) <- c(attributes(out), attributes(object)[restaure.attributes])
        }

        return(out)
    }else{
        return(object)
    }

}
restaureNA.character <- restaureNA.numeric
restaureNA.factor <- restaureNA.numeric
restaureNA.integer <- restaureNA.numeric
restaureNA.logical <- restaureNA.numeric

## * restaureNA.matrix
restaureNA.matrix <- function(object, index.na, level, cluster){

    if(level == "cluster" && any(is.na(cluster$index))){
        ## expand to store NAs
        out <- matrix(NA, nrow = cluster$n, ncol = NCOL(object),
                      dimnames = list(cluster$levels, colnames(object)))
        out[rownames(object),] <- unname(object)
        ## put back attributes
        restaure.attributes <- setdiff(names(attributes(object)), names(attributes(out)))
        if(length(restaure.attributes)>0){
            attributes(out) <- c(attributes(out), attributes(object)[restaure.attributes])
        }

        return(out)
    }else if(level == "obs" && length(index.na)>0){
        ## expand to store NAs
        nobs <- NROW(object) + length(index.na)
        out <- matrix(NA, nrow = nobs, ncol = NCOL(object),
                      dimnames = list(NULL, colnames(object)))
        out[-index.na,] <- unname(object)
        ## put back attributes
        restaure.attributes <- setdiff(names(attributes(object)), names(attributes(out)))
        if(length(restaure.attributes)>0){
            attributes(out) <- c(attributes(out), attributes(object)[restaure.attributes])
        }

        return(out)
    }else{
        return(object)
    }

}

## * restaureNA.array
restaureNA.array <- function(object, index.na, level, cluster){

    if(level == "cluster" && any(is.na(cluster$index))){
        ## expand to store NAs
        out <- array(NA, dim = c(cluster$n, dim(object)[-1]),
                     dimnames = c(list(cluster$levels), dimnames(object)[-1]))
        out[dimnames(object)[[1]],,] <- unname(object)
        ## put back attributes
        restaure.attributes <- setdiff(names(attributes(object)), names(attributes(out)))
        if(length(restaure.attributes)>0){
            attributes(out) <- c(attributes(out), attributes(object)[restaure.attributes])
        }
        
        return(out)
    }else if(level == "obs" && length(index.na)>0){
        stop("Not yet implemented - contact the package maintainer")
    }else{
        return(object)
    }

}

##----------------------------------------------------------------------
### restaureNA.R ends here
