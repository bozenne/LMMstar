### correlate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2024 (16:41) 
## Version: 
## Last-Updated: jul  8 2024 (12:39) 
##           By: Brice Ozenne
##     Update #: 51
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' attr(e.S, "correlation")
##' 
##' ## compute correlations (composite time variable)
##' dL$time2 <- dL$time == 2
##' dL$time3 <- dL$time == 3
##' e.S <- summarize(Y ~ time2 + time3 + X1 | id, data = dL, na.rm = TRUE)
##' e.S
##' attr(e.S, "correlation")


## * correlate (code)
##' @export
correlate <- function(formula, data, repetition, use = "everything", method = "pearson", keep.name = TRUE, na.rm){

    mycall <- match.call()
    options <- LMMstar.options()

    ## ** normalize user input
    
    ## formula & repetition
    ls.init <- formula2repetition(formula = formula, data = data, repetition = repetition, keep.time = FALSE)
    detail.formula <- ls.init$detail.formula
    detail.repetition <- ls.init$detail.repetition
    if(is.null(detail.repetition)){
        stop("Argument \'repetition\' is missing. \n",
             "Should be a formula such as ~time|cluster. \n")
    }
    name.cluster <- detail.repetition$var$cluster
    name.time <- detail.repetition$var$time
    name.Y <- detail.formula$var$response
    name.X <- detail.formula$var$regressor
    n.Y <- length(name.Y)

    if(length(name.time)>1 && paste(name.time, collapse = "_X_XX_X_") %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"",paste(name.time, collapse = "_X_XX_X_"),"\" as this name is used internally. \n")
    }

    ## na.rm
    if(!missing(na.rm)){
        message("Argument \'na.rm\' is ignored, specify instead argument \'use\'. \n",
                "with \"everything\" for na.rm=FALSE and \"pairwise.complete.obs\" for na.rm = TRUE. \n")
    }
    
    ## ** evaluate correlation
    out <- stats::setNames(vector(mode = "list", length = length(name.Y)), name.Y)
    if(is.null(name.X)){
        ls.cluster <- list(unique(data[[name.cluster]]))
    }else{
        ls.cluster <- tapply(data[[name.cluster]],data[,name.X],unique)
        if(any(duplicated(unlist(ls.cluster)))){
            stop("Argument \'formula\' incompatible with argument \'repetition\': covariate values should not change within cluster. \n")
        }
    }

    for(iY in 1:n.Y){ ## iY <- 1
        
        out[[iY]] <- stats::setNames(lapply(ls.cluster, function(iId){ ## iId <- cluster[1]
            iDataL <- data[data[[name.cluster]] %in% iId,,drop = FALSE]
            if(length(name.time)>1){
                iDataL[[paste(name.time, collapse = "_X_XX_X_")]] <- nlme::collapse(iDataL[,name.time, drop=FALSE], as.factor = TRUE)
                Utime <- paste(name.time, collapse = "_X_XX_X_")
            }else{
                Utime <- name.time
            }
            iDataW <- stats::reshape(data = iDataL[,c(name.cluster, Utime, name.Y[iY])],
                                     direction = "wide", timevar = Utime, idvar = name.cluster, v.names = name.Y[iY])

            iCor <- stats::cor(iDataW[,-1,drop=FALSE], use = use, method =  method)
            iLevels <- levels(as.factor(iDataL[[Utime]]))
            if(NROW(iCor)==length(iLevels)){
                dimnames(iCor) <- list(iLevels,iLevels)
            }
            return(iCor)
        }), names(ls.cluster))

    }

    ## ** export
    attr(out,"call") <- mycall
    attr(out,"method") <- method
    attr(out,"use") <- use
    attr(out,"name.Y") <- name.Y
    attr(out,"name.X") <- name.X
    attr(out,"name.time") <- name.time
    attr(out,"name.cluster") <- name.cluster
    class(out) <- append("correlate",class(out))
    return(out)
}

##----------------------------------------------------------------------
### correlate.R ends here
