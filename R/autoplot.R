### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  8 2021 (00:01) 
## Version: 
## Last-Updated: nov  4 2021 (16:28) 
##           By: Brice Ozenne
##     Update #: 69
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * autoplot.lmm (documentation)
##' @title Graphical Display For Linear Mixed Models
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
##' @param size.text [numeric, >0] size of the font used to displayed text when using ggplot2.
##' @param position.errorbar [character] relative position of the errorbars.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to create the graphical display.
##' \item \code{plot}: ggplot object.
##' }


## * autoplot.lmm (code)
##' @rdname autplot
##' @export
autoplot.lmm <- function(object, at = NULL, color = TRUE, ci = TRUE, alpha = NA, plot = TRUE,
                         size.point = 3, size.line = 1, size.text = 16, position.errorbar = "identity", ...){

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

    ## design matrix
    X.beta <- model.matrix(object, data = data, effects = "mean")

    ## only keep one representant per type of design matrix
    test.duplicated <- duplicated(X.beta)
    keep.id <- unique(data[test.duplicated==FALSE,var.cluster])
    newdata <- data[data[[var.cluster]] %in% keep.id,,drop=FALSE]

    if(identical(color,TRUE)){
        mean.var <- all.vars(stats::delete.response(stats::terms(stats::formula(object, effects = "mean"))))
        newdataRed <- newdata[order(newdata[[var.cluster]]),mean.var,drop=FALSE]
        order.cluster <- droplevels(newdata[[var.cluster]][order(newdata[[var.cluster]])])

        ## iCol <- newdataRed[,1]
        M.duplicated <- apply(newdataRed, 2, function(iCol){unlist(tapply(iCol, order.cluster, function(iColCluster){duplicated(iColCluster)[-1]}))})
        color <- names(which(colSums(M.duplicated)==NROW(M.duplicated)))
        if(length(color)>1){
            if(paste(color,collapse=".") %in% names(newdataRed)){
                stop("Cannot use argument \'color\'=TRUE when the dataset contain a column ",paste(color,collapse="."),". \n",
                     "This name is used internally. \n")
            }
            newdata[[paste(color,collapse=".")]] <- interaction(newdata[,color,drop=FALSE])
            color <- paste(color,collapse=".")
        }else if(length(color)==0){
            color <-  NULL
        }
    }
    
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
    gg  <- gg + ggplot2::ylab(object$outcome$var) + ggplot2::theme(text = ggplot2::element_text(size=size.text))


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
