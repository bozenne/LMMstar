### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  8 2021 (00:01) 
## Version: 
## Last-Updated: Jun  8 2021 (00:27) 
##           By: Brice Ozenne
##     Update #: 27
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * autoplot.lmm (documentation)
##' @title Graphical Display For Multivariate Gaussian Model
##' @name autoplot
##'
##' @param object a \code{lmm} object.
##' @param color [character] name of the variable in the dataset used to color the curve.
##' @param ci [logical] should confidence intervals be displayed?
##' @param alpha [numeric, 0-1] When not NA, transparency parameter used to display the confidence intervals.
##' @param plot [logical] should the plot be displayed?
##' @param size.point [numeric, >0] the size of the point on the plot.
##' @param size.line [numeric, >0] the size of the line on the plot.
##' @param position.error [character] relative position of the errorbars.
##' @param ... Not used. For compatibility with the generic method.

## * autoplot.lmm (code)
##' @rdname autplot
##' @export
autoplot.lmm <- function(object, color = object$cluster$var, ci = TRUE, alpha = NA, plot = TRUE,
                         size.point = 3, size.line = 1, position.errorbar = "identity", ...){

    var.cluster <- object$cluster$var
    var.time <- object$time$var
    
    ## ** find representative individuals
    reorder <- order(object$data[[var.cluster]],object$data[[var.time]])
    data <- object$data[reorder,,drop=FALSE]
    X <- model.matrix(object, effects = "mean")[reorder,,drop=FALSE]
    test.duplicated <- duplicated(X)
    keep.id <- unique(data[test.duplicated==FALSE,var.cluster])
    newdata <- data[data[[var.cluster]] %in% keep.id,,drop=FALSE]

    ## ** compute fitted curve
    preddata <- cbind(newdata, stats::predict(object, newdata = newdata))

    ## ** generate plot
    gg <- ggplot2::ggplot(preddata, ggplot2::aes_string(x = var.time, y = "estimate", group = var.cluster))
    if(ci){
        if(is.na(alpha)){
            gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), position = position.errorbar)
        }else{
            if(!is.null(color)){
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lower", ymax = "upper", fill = color), alpha = alpha)
            }else{
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = alpha)
            }
        }
    }
    if(!is.null(color)){
        gg <- gg + ggplot2::geom_point(ggplot2::aes_string(color = color), size = size.point) + ggplot2::geom_line(ggplot2::aes_string(color = color), size = size.line)
    }else{
        gg <- gg + ggplot2::geom_point(size = size.point) + ggplot2::geom_line(size = size.line)
    }
    gg  <- gg + ggplot2::ylab(object$outcome$var)


    ## ** display
    if(plot){
        print(gg)
    }
    
    ## ** export
    return(invisible(list(data = preddata,
                          plot = gg)))
}

##----------------------------------------------------------------------
### autoplot.R ends here
