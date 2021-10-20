### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (11:00) 
## Version: 
## Last-Updated: okt 20 2021 (14:56) 
##           By: Brice Ozenne
##     Update #: 21
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Graphical Display For Linear Mixed Models
##' @description Display fitted values or residual plot for the mean, variance, and correlation structure.
##' Can also display quantile-quantile plot relative to the normal distribution.
##' @name plot
##' 
##' @param x a \code{lmm} object.
##' @param type [character] the type of plot: \code{"fit"}, \code{"qqplot"}, \code{"correlation"}, \code{"scatterplot"}, \code{"scatterplot2"}.
##' @param type.residual [character] the type of residual to be used. Not relevant for \code{type="fit"}. By default, normalized residuals are used.
##' @param by.time [logical] should a separate plot be made at each repetition or a single plot over all repetitions be used?
##' Only relevant for \code{type="qqplot"}, \code{type="scatterplot"}, and \code{type="scatterplot2"}.
##' @param ... additional argument passed to \code{residuals.lmm} or \code{autoplot.lmm}.
##'
##' @details Call \code{link(autoplot.lmm)} when code{type=="fit"} and \code{link(residuals.lmm)} for the other types.
##'
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to create the graphical display.
##' \item \code{plot}: ggplot object.
##' }
##' 
##' @export
plot.lmm <- function(x, type = "fit", type.residual = "normalized", by.time = TRUE, ...){
    type <- match.arg(type, c("qqplot","correlation","scatterplot","scatterplot2","fit"))
    if(by.time){
        format <- "wide"
    }else{
        format <- "long"
    }
    if(type=="fit"){
        ## requireNamespace("ggplot2")
        out <- autoplot.lmm(x, ...)
    }else{
        outRes <- residuals(x, plot = type, type = type.residual, format = format, ...)
        out <- list(data = outRes,
                    plot = attr(outRes,"plot"))
        attr(out$data,"plot") <- NULL
    }
    return(invisible(out))
}

##----------------------------------------------------------------------
### plot.R ends here
