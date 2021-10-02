### summarize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: okt  1 2021 (17:10) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summarize (documentation)
##' @title Compute summary statistics
##' @description Compute summary statistics (similar to the SAS macro procmean).
##' This is essentially an interface to the \code{stats::aggregate} function.
##'
##' @param formula [formula] on the left hand side the outcome(s) and on the right hand side the grouping variables.
##' E.g. Y1+Y2 ~ Gender + Gene will compute for each gender and gene the summary statistics for Y1 and for Y2.
##' Passed to the \code{stats::aggregate} function.
##' @param data [data.frame] dataset (in the wide format) containing the observations.
##' @param na.action [function] a function which indicates what should happen when the data contain 'NA' values.
##' Passed to the \code{stats::aggregate} function.
##' @param na.rm [logical] Should the summary statistics be computed by omitting the missing values.
##' @param which [character vector] name of the summary statistics to kept in the output.
##' Can be any of, or a combination of: \code{"observed"} (number of observations with a measurement),
##' \code{"missing"} (number of observations with a missing value), \code{"mean"}, \code{"sd"}, \code{"min"}, \code{"median"}, \code{"max"}.
##'
##' @return a data frame containing summary statistics (in columns) for each outcome and value of the grouping variables (rows).


## * summarize (examples)
##' @examples
##' ## simulate data in the wide format
##' set.seed(10)
##' d <- sampleRem(1e2, n.times = 3)
##'
##' ## add a missing value
##' d2 <- d
##' d2[1,"Y2"] <- NA
##'
##' ## run summarize
##' summarize(Y1+Y2 ~ 1, data = d)
##' summarize(Y1+Y2 ~ X1, data = d)
##' 
##' summarize(Y1+Y2 ~ X1, data = d2)
##' summarize(Y1+Y2 ~ X1, data = d2, na.rm = TRUE)
##' 
##' ## End of examples


## * summarize (code)
##' @export
summarize <- function(formula, data, na.action = stats::na.pass, na.rm = FALSE,
                     which = c("observed","missing","mean","sd","min","median","max")){

    ## ** check and normalize user imput
    valid.which <- c("observed","missing","mean","sd","min","median","max")
    
    which <- match.arg(which, choices = valid.which, several.ok = TRUE)
    name.all <- all.vars(formula)
    if(any(name.all %in% c(valid.which, "outcome"))){
        invalid <- name.all[name.all %in% c(valid.which, "outcome")]
        stop("Name(s) \"",paste(invalid, collapse = "\" \""),"\" are used internally. \n",
             "Consider renaming the variables in the dataset and updating the formula. \n",
             sep = "")
    }
    if(any(name.all %in% names(data) == FALSE)){
        invalid <- name.all[name.all %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }

    name.Y <- lhs.vars(formula)
    n.Y <- length(name.Y)
    if(n.Y==0){
        stop("Wrong specification of argument \'formula\'. \n",
             "there need to be at least one variable in the left hand side of the formula. \n")
    }
    
    name.X <- rhs.vars(formula)
    n.X <- length(name.Y)

    
    ## ** compute summary statistics
    out <- NULL
    for(iY in 1:n.Y){ ## iY <- 1
        iFormula <- stats::update(formula, paste0(name.Y[iY],"~."))
        iAggregate <- stats::aggregate(iFormula, data=data, function(x){
            c("observed" = sum(!is.na(x)),
              "missing" = sum(is.na(x)),
              "mean" = mean(x, na.rm = na.rm),
              "sd" = stats::sd(x, na.rm = na.rm),
              "min" = min(x, na.rm = na.rm),
              "median" = stats::median(x, na.rm = na.rm),
              "max" = max(x, na.rm = na.rm))},
            na.action=na.action)
        iDF <- cbind(outcome = name.Y[iY], iAggregate[name.X], iAggregate[[name.Y[iY]]][,which,drop=FALSE])
        out <- rbind(out,iDF)
    }

    if(na.rm){
        attr(out,"correlation") <- stats::cor(data[,name.Y,drop=FALSE], use = "pairwise")
    }else{
        attr(out,"correlation") <- stats::cor(data[,name.Y,drop=FALSE])
    }
    
    ## ** export
    return(out)
}


######################################################################
### summarize.R ends here
