### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  8 2021 (00:01) 
## Version: 
## Last-Updated: Jun  8 2021 (12:21) 
##           By: Brice Ozenne
##     Update #: 41
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
##' @param at [data.frame] values for the covariates at which to evaluate the fitted values.
##' @param color [character] name of the variable in the dataset used to color the curve.
##' @param ci [logical] should confidence intervals be displayed?
##' @param alpha [numeric, 0-1] When not NA, transparency parameter used to display the confidence intervals.
##' @param plot [logical] should the plot be displayed?
##' @param size.point [numeric, >0] the size of the point on the plot.
##' @param size.line [numeric, >0] the size of the line on the plot.
##' @param position.errorbar [character] relative position of the errorbars.
##' @param ... Not used. For compatibility with the generic method.

## * autoplot.lmm (code)
##' @rdname autplot
##' @export
autoplot.lmm <- function(object, at = NULL, color = object$cluster$var, ci = TRUE, alpha = NA, plot = TRUE,
                         size.point = 3, size.line = 1, position.errorbar = "identity", ...){

    var.cluster <- object$cluster$var
    var.time <- object$time$var
    
    ## ** find representative individuals
    reorder <- order(object$data[[var.cluster]],object$data[[var.time]])
    data <- object$data[reorder,,drop=FALSE]
    if(!is.null(at)){
        if(is.vector(at)){at <- as.data.frame(as.list(at))}
        if(!is.data.frame(at)){
            stop("Argument \'at\' must be a data.frame or NULL. \n")
        }
        if(NROW(at)!=1){
            stop("Argument \'at\' must have exactly one row. \n")
        }
        
        if(any(names(at) %in% names(data) == FALSE)){
            stop("Argument \'at\' contains variables not used by the model. \n",
                 "Incorrect variable: \"",paste(names(at)[names(at) %in% names(data) == FALSE], collapse = "\" \""),"\" \n")
        }
        for(iCol in names(at)){
            data[[iCol]] <- at[[iCol]]
        }
    }
    ff.mean <- stats::formula(object, effects = "mean")
    X.name <- colnames(model.matrix(object, effects = "mean"))
    X <- model.matrix(ff.mean, data)[,X.name,drop=FALSE]
    test.duplicated <- duplicated(X)
    keep.id <- unique(data[test.duplicated==FALSE,var.cluster])
    newdata <- data[data[[var.cluster]] %in% keep.id,,drop=FALSE]

    ## ** compute fitted curve
    preddata <- cbind(newdata, stats::predict(object, newdata = newdata))

    ## ** generate plot
    gg <- ggplot2::ggplot(preddata, ggplot2::aes_string(x = var.time, y = "estimate", group = var.cluster))
    if(ci){
        if(is.na(alpha)){
            gg <- gg + ggplot2::geom_errorbar(ggplot2::aes_string(ymin = "lower", ymax = "upper"), position = position.errorbar)
        }else{
            if(!is.null(color)){
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lower", ymax = "upper", fill = color), alpha = alpha)
            }else{
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lower", ymax = "upper"), alpha = alpha)
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
