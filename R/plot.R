### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (11:00) 
## Version: 
## Last-Updated: jan 23 2023 (16:04) 
##           By: Brice Ozenne
##     Update #: 109
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * plot (code)
##' @describeIn autoplot.lmm Graphical Display For Linear Mixed Models
##' @export
plot.lmm <- function(x, ...){

    out <- autoplot.lmm(x, ...)
    if(!is.null(out$plot)){
        ## if ggplot then the graph is stored in out and should be displayed here
        ## otherwise the graph has already been displayed but could not be stored in the object
        print(out$plot)
    }
    return(invisible(out))

}

## * plot (code)
##' @describeIn autoplot.partialCor Graphical Display For Partial Correlation
##' @export
plot.partialCor <- function(x, ...){
    out <- autoplot.partialCor(x, ...)
    print(out$plot)
    return(invisible(out))

}

## * plot (code)
##' @describeIn autoplot.profile_lmm Display Contour of the log-Likelihood
##' @export
plot.profile_lmm <- function(x, ...){
    out <- autoplot.profile_lmm(x, ...)
    print(out$plot)
    return(invisible(out))

}

## * plot (code)
##' @describeIn autoplot.residuals_lmm Graphical Display of the Residuals
##' @export
plot.residuals_lmm <- function(x, ...){
    out <- autoplot.residuals_lmm(x, ...)
    print(out$plot)
    return(invisible(out))

}

## * plot (code)
##' @describeIn autoplot.summarizeNA Graphical Display of Missing Data Pattern
##' @export
plot.summarizeNA <- function(x, ...){
    out <- autoplot.summarizeNA(x, ...)
    print(out$plot)
    return(invisible(out))
}

## * plot (code)
##' @describeIn autoplot.Wald_lmm Graphical Display For Linear Hypothesis Test
##' @export
plot.Wald_lmm <- function(x, ...){
    out <- autoplot.Wald_lmm(x, ...)
    print(out$plot)
    return(invisible(out))
}

##----------------------------------------------------------------------
### plot.R ends here
