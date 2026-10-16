### index.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 16 2026 (15:26) 
## Version: 
## Last-Updated: okt 16 2026 (16:11) 
##           By: Brice Ozenne
##     Update #: 38
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
##' @param format observation, cluster, both
##' @noRd
##' @keywords internal

## * index.lmm (code)
index.lmm <- function(x, variable, format = "both", type = "original", ...){


    if(type != "numeric"){
        table.variables <- stats::model.tables(x, effects = "variables")[variable,]
        if(is.character(x$data.original[[table.variables$var]])){
            type <- "character"
        }else if(is.character(x$data.original[[table.variables$var]])){
            type <- "factor"
        } ## original means it was numeric
    }
browser()
    ## ** extract information
    if(variable == "strata"){
        if(format == "observation"){
            if(type == "numeric"){
                return(x$design$index.clusterStrata)
            }else if(type == "factor"){
                return(factor(x$design$index.clusterStrata, labels = table.variables["strata","levels"][[1]]))
            }else if(type == "character"){
                return(as.character(factor(x$design$index.clusterStrata, labels = table.variables["strata","levels"][[1]])))
            }
        }
        out <- lapply(x$design$index.cluster, function(iIndex){x$design$index.clusterStrata[[iIndex]]})
        attr(out,"vectorwise") <- x$design$index.clusterStrata
    }else{
        element.design <- switch(variable,
                                 cluster = "index.cluster",
                                 time = "index.clusterTime",
                                 strata = "index.clusterStrata")
        out <- x$design[[element.design]]
    }

    ## ** re-type
    if(type != "numeric"){

        if(format %in% c("cluster","both")){
            out.max <- max(attr(out,"vectorwise"))
        }

        if(format %in% c("observation","both")){
            out2 <- attr(out,"vectorwise")
        }
        if(format %in% c("observation","both")){
            out2 <- factor(out2, labels = table.variables[variable,"levels"][[1]])
            if(type=="character"){
                out2 <- as.character(out2)
            }
        }
        if(format %in% c("cluster","both")){
            if(type == "factor"){
                out <- lapply(out, function(iIndex){factor(iIndex, levels = 1:out.max, labels = table.variables[variable,"levels"][[1]])})
            }else if(type == "character"){
                out <- lapply(out, function(iIndex){as.character(factor(iIndex, levels = 1:out.max, labels = table.variables[variable,"levels"][[1]]))})
            }
            
        }
    }

    ## ** export according to format
    if(format %in% c("observation","both")){
        if(format == "both"){
            attr(out, "vectorwise") <- out2
        }        
        return(out)
    }else{
        return(out2)
    }
    
}

##----------------------------------------------------------------------
### index.R ends here
