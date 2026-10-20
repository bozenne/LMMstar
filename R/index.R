### index.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 16 2026 (15:26) 
## Version: 
## Last-Updated: okt 19 2026 (12:22) 
##           By: Brice Ozenne
##     Update #: 83
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @description Extract position of the clusters, times, strata
##' For internal use.
##' @param format observation: vector with one value per observation (i.e. per time and cluster)
##' cluster: list with one vector per cluster which has as many values as observation within cluster.
##' both: list with vector as attribute names \code{"vectorwise"}
##' @param type numeric: index of the timepoint/strata/cluster (1,2,3)
##' character: character for the timepoint/strata/cluster
##' factor: factor for the timepoint/strata/cluster
##' original: try to match the original variable so can be character or factor or if the original variable was numeric re-create the original numeric values (not necessarily 1,2,3... can be 1,5,7...)
##' index: for cluster return the index of the observation instead of the index of the cluster.
##' 
##' @noRd
##' @keywords internal

## * index.lmm (code)
index.lmm <- function(x, variable, format = "both", type = "original", ...){

    ## ** normalize user input
    table.variables <- stats::model.tables(x, effects = "variables")
    variable <- match.arg(variable, c("outcome","strata","time","cluster"))
    table.variables <- table.variables[variable,]
    format <- match.arg(format, c("observation","cluster","both"))
    type <- match.arg(type, c("numeric","original","factor","character","index"))
    
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    ## ** extract information
    ## *** variables
    if(type == "numeric"){
        newlevel <- 1:length(table.variables[["levels"]][[1]])
    }else if(type != "numeric"){
        if(type == "original" && is.character(x$data[[table.variables$var]])){
            type <- "character"
        }else if(type == "original" && is.factor(x$data[[table.variables$var]])){
            type <- "factor"
        } ## original means it was numeric

        if(type %in% c("character","factor") || is.null(attr(table.variables[["levels"]][[1]],"original"))){
            newlevel <- table.variables[["levels"]][[1]]
        }else{ ## original means it was numeric
            newlevel <- attr(table.variables[["levels"]][[1]],"original")[,1]
        }
    }

    ## *** index
    if(format %in% c("observation","both")){
        if(variable == "cluster"){
            index.obs <- attr(x$design$index.cluster,"vectorwise")
        }else if(variable == "strata"){
            index.obs <- x$design$index.clusterStrata
        }else  if(variable == "time"){
            index.obs <- attr(x$design$index.clusterTime,"vectorwise")
        }
    }
    if(format %in% c("cluster","both")){
        if(variable == "cluster"){
            index.cluster <- x$design$index.cluster
            attr(index.cluster,"vectorwise") <- NULL
        }else if(variable == "strata"){
            index.cluster <- mapply(iValue = x$design$index.clusterStrata, iRep = lengths(x$design$index.cluster), function(iValue,iRep){rep(iValue, times = iRep)},  SIMPLIFY = FALSE)
        }else  if(variable == "time"){
            index.cluster <- x$design$index.clusterTime
            attr(index.cluster,"vectorwise") <- NULL
        }
    }  

    ## ** long format: vector with 1 value per observation and cluster
    if(format %in% c("observation","both") && type != "index"){
        if(type %in% c("original","character")){
            index.obs <- newlevel[index.obs]
        }else if(type == "factor"){
            index.obs <- factor(index.obs, labels = newlevel)
        }
    }

    ## ** wide format: 1 vector per cluster
    if(format %in% c("cluster","both") && type != "index"){
        if(variable == "time"){
            if(type %in% c("original","character")){
                index.cluster <- lapply(index.cluster, FUN = function(iVec){newlevel[iVec]})
            }else if(type == "factor"){
                index.cluster <- lapply(index.cluster, FUN = function(iVec){factor(iVec, levels = newlevel)})
            }
        }else{
            if(type %in% c("original","character","numeric")){
                index.cluster <- mapply(iT = lengths(index.cluster), iV =  newlevel, FUN = function(iT,iV){rep(iV,times = iT)},SIMPLIFY = FALSE)
            }else if(type == "factor"){
                index.cluster <- mapply(iT = lengths(index.cluster), iV =  factor(newlevel, levels = newlevel), FUN = function(iT,iV){rep(iV,times = iT)},SIMPLIFY = FALSE)
            }
        }
    }

    ## ** export according to format
    if(format == "observation"){
        attr(out, "vectorwise") <- index.obs
    }else if(format %in% c("cluster","both")){
        out <- index.cluster
        if(format == "both"){
            attr(out, "vectorwise") <- index.obs
        }
    }
    return(out)
    
}

##----------------------------------------------------------------------
### index.R ends here
