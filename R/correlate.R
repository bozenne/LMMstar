### correlate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2024 (16:41) 
## Version: 
## Last-Updated: jul  9 2024 (15:17) 
##           By: Brice Ozenne
##     Update #: 88
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * correlate (documentation)
##' @title Compute correlation matrix
##' @description Compute correlation matrix for multiple variables and/or multiple groups.
##'
##' @param formula [formula] on the left hand side the outcome(s) and on the right hand side the grouping variables.
##' E.g. Y1+Y2 ~ Gender will compute for each gender the correlation matrix for Y1 and the correaltion matrix for Y2.
##' @param data [data.frame] dataset containing the observations.
##' @param repetition [formula] Specify the structure of the data: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param use [character] method for computing correlation in the presence of missing values: \code{"everything"}, \code{"all.obs"}, \code{"complete.obs"}, \code{"na.or.complete"}, or \code{"pairwise.complete.obs"}.
##' Passed to \code{\link{stats::cor}}.
##' @param method [character] type of correlation coefficient: \code{"pearson"}, \code{"kendall"}, \code{"spearman"}.
##' Passed to \code{\link{stats::cor}}.
##' @param collapse.value [character] symbol used to combine covariate values when using multiple grouping variables.
##' @param collapse.var [character] symbol used to combine variable names to be pasted left of the covariate values when using multiple grouping variables.
##' Can be disabled setting it to \code{NULL} or \code{FALSE}.
##' @param na.rm [logical] not used. The user may expect this argument though so it is added to help the user with a message pointing toward the argument \code{use}.
##'
##' @seealso
##' \code{\link{as.array}} to convert the output to an array or \code{\link{as.matrix}} to convert the output to a matrix (when a single outcome and no grouping variable). \cr
##' \code{\link{summarize}} for other summary statistics (e.g. mean, standard deviation, ...). 
##' 
##' @examples 
##' #### simulate data (wide format) ####
##' data(gastricbypassL, package = "LMMstar")
##' 
##' ## compute correlations (composite time variable)
##' e.S <- correlate(weight ~ 1, data = gastricbypassL, repetition = ~time|id)
##' e.S
##' as.matrix(e.S)
##' 
##' e.S21 <- correlate(glucagonAUC + weight ~ 1, data = gastricbypassL,
##'                   repetition = ~time|id, use = "pairwise")
##' e.S21
##' as.array(e.S21)
##' as.matrix(e.S21, index = "weight")
##'
##' gastricbypassL$sex <- as.numeric(gastricbypassL$id) %% 2
##' e.S12 <- correlate(weight ~ sex, data = gastricbypassL,
##'                    repetition = ~time|id, use = "pairwise")
##' e.S12
##' as.array(e.S12)
##' 
##' e.S22 <- correlate(glucagonAUC + weight ~ sex, data = gastricbypassL,
##'                    repetition = ~time|id, use = "pairwise")
##' e.S22
##' as.array(e.S22)


## * correlate (code)
##' @export
correlate <- function(formula, data, repetition, use = "everything", method = "pearson",
                      collapse.value = ".", collapse.var = ".", na.rm){

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
        ls.cluster <- tapply(data[[name.cluster]], INDEX = interaction(data[name.X],drop=FALSE,sep=collapse.value), FUN = unique, simplify = FALSE)
        if(any(duplicated(unlist(ls.cluster)))){
            stop("Argument \'formula\' incompatible with argument \'repetition\': covariate values should not change within cluster. \n")
        }
        if(!is.null(collapse.var) & !identical(collapse.var,FALSE)){
            names(ls.cluster) <- paste(paste(name.X, collapse = collapse.var), names(ls.cluster), sep = "=")
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
