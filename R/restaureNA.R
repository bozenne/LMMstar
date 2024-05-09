### restaureNA.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 20 2023 (15:31) 
## Version: 
## Last-Updated: May  9 2024 (10:38) 
##           By: Brice Ozenne
##     Update #: 65
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * restaureNA (documentation)
##' @title Restaure NA
##' @description Restaure NA in the user output that have been excluded when fitting the LMM.
##' @param object results when NA should be added to match the user input.
##' @param index.na [integer vector] index of the missing values.
##' @param level [character] Should missing observations, indicated by the argument \code{index.na}, be restaured (\code{"obs"})
##' or should missing clusters, indicated by the discrepancy between the name of the object and the argument \code{cluster}, be restaured (\code{"cluster"}).
##' @param cluster [list] list containing the number of cluster (\code{n}), name of the clusters (\code{levels}), and integer associated with each clusters (\code{index}).
##'
##' @keywords internal
##' 
##' @export
`restaureNA` <-
  function(object, index.na, level, cluster) UseMethod("restaureNA")

## * restaureNA.numeric
##' @export
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
        if(!is.null(names(object))){
            names(out)[-index.na] <- names(object)
            names(out)[index.na] <- "NA"
        }
        return(out)
    }else{
        return(object)
    }

}

## * restaureNA.character
##' @export
restaureNA.character <- restaureNA.numeric

## * restaureNA.factor
##' @export
restaureNA.factor <- restaureNA.numeric

## * restaureNA.integer
##' @export
restaureNA.integer <- restaureNA.numeric

## * restaureNA.logical
##' @export
restaureNA.logical <- restaureNA.numeric

## * restaureNA.matrix
##' @export
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
##' @export
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

        ## expand to store NAs
        nobs <- dim(object)[1] + length(index.na)
        out <- array(NA, dim = c(nobs, dim(object)[-1]),
                      dimnames = c(list(NULL), dimnames(object)[-1]))
        out[-index.na,,] <- as.double(object)
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

##----------------------------------------------------------------------
### restaureNA.R ends here
