### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (11:00) 
## Version: 
## Last-Updated: nov 13 2021 (17:59) 
##           By: Brice Ozenne
##     Update #: 64
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * plot (documentation)
##' @title Graphical Display For Linear Mixed Models
##' @description Display fitted values or residual plot for the mean, variance, and correlation structure.
##' Can also display quantile-quantile plot relative to the normal distribution.
##' @name plot
##' 
##' @param x a \code{lmm} object.
##' @param type [character] the type of plot: \code{"fit"}, \code{"qqplot"}, \code{"correlation"}, \code{"scatterplot"}, \code{"scatterplot2"}, \code{"partial"}.
##' @param type.residual [character] the type of residual to be used. Not relevant for \code{type="fit"}.
##' By default, normalized residuals are used except when requesting a partial residual plot.
##' @param by.time [logical] should a separate plot be made at each repetition or a single plot over all repetitions be used?
##' Only relevant for \code{type="qqplot"}, \code{type="scatterplot"}, and \code{type="scatterplot2"}.
##' @param ci [logical] should confidence intervals be displayed?
##' @param plot [logical] should the plot be displayed?
##' @param ci.alpha [numeric, 0-1] Transparency parameter used to display the confidence intervals.
##' @param mean.size [numeric vector of length 2] size of the point and line for the mean trajectory.
##' @param size.text [numeric, >0] size of the font used to displayed text when using ggplot2.
##' @param ... additional argument passed to \code{residuals.lmm} or \code{autoplot.lmm}.
##'
##' @details Call \code{\link{autoplot.lmm}} when code{type=="fit"} and \code{link(residuals.lmm)} for the other types.
##'
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to create the graphical display.
##' \item \code{plot}: ggplot object.
##' }

 
## * plot (code)
##' @export
plot.lmm <- function(x, type = "fit", type.residual = "normalized", by.time = TRUE, ci = TRUE,
                     plot = TRUE, ci.alpha = 0.2, mean.size = c(3, 1), size.text = 16, ...){
    type <- match.arg(type, c("qqplot","correlation","scatterplot","scatterplot2","fit","partial"))
    if(by.time){
        format <- "wide"
    }else{
        format <- "long"
    }
    
    if(type=="fit"){
        ## requireNamespace("ggplot2")
        out <- autoplot.lmm(x, ci = ci, plot = plot, ci.alpha = ci.alpha, mean.size = mean.size, size.text = size.text, ...)
    }else if(type=="partial"){
        ## prepare
        if(is.null(match.call()$type.residual)){
            if(!is.null(list(...)$var)){
                type.residual <- list(...)$var
            }else{
                type.residual <- attr(x$design$mean,"variable")[1]
            }
        }        
        name.var <- setdiff(type.residual,"(Intercept)")
        if("(Intercept)" %in% type.residual){
            type.predict <- "static"
        }else{
            type.predict <- "static0"
        }
        type.var <- c("numeric","categorical")[name.var %in% names(x$xfactor$mean) + 1]
        if(sum(type.var=="numeric")>1){
            stop("Cannot simulatenously display partial residuals for more 2 numeric variables. \n")
        }

        ## extract partial residuals
        rr <- stats::residuals(x, type = "partial-ref", var = type.residual, keep.data = TRUE)
        ## extract predictions
        gg.data <- stats::predict(x, newdata = rr, keep.newdata = TRUE, type = type.predict)
        ## normalize covariates
        if(sum(type.var=="categorical")==0){
            name.varcat <- NULL
        }else if(sum(type.var=="categorical")==1){
            name.varcat <- paste(name.var[type.var=="categorical"],collapse = ", ")
            if(sum(type.var=="categorical")>1){
                gg.data[[name.varcat]] <- interaction(gg.data[name.var[type.var=="categorical"]], sep = ",")
            }
        }
        if(sum(type.var=="numeric")==0){
            name.varnum <- NULL
        }else{
            name.varnum <- name.var[type.var=="numeric"]
        }
        ## display
        if(all(type.var=="categorical")){
            gg <- ggplot2::ggplot(data = gg.data, mapping = ggplot2::aes_string(x = name.varcat))
            gg <- gg + ggplot2::geom_point(ggplot2::aes_string(y = "r.partial"), color = "gray")
            if(ci){
                gg <- gg + ggplot2::geom_errorbar(ggplot2::aes_string(ymin = "lower", ymax = "upper"))
            }
            gg <- gg + ggplot2::geom_point(ggplot2::aes_string(y = "estimate"), size = mean.size[1], shape = 21, fill = "white")
        }else{
            gg <- ggplot2::ggplot(data = gg.data, mapping = ggplot2::aes_string(x = name.varnum))
            gg <- gg + ggplot2::geom_point(ggplot2::aes_string(y = "r.partial"), color = "gray")
            if(length(type.var)==1){
                gg <- gg + ggplot2::geom_line(ggplot2::aes_string(y = "estimate"), size = mean.size[2])
                if(ci){
                    gg <- gg + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lower", ymax = "upper"), alpha = ci.alpha)
                }
            }else{
                gg <- gg + ggplot2::geom_line(ggplot2::aes_string(y = "estimate", group = name.varcat, color = name.varcat), size = mean.size[2])
                if(ci){
                    gg <- gg + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lower", ymax = "upper", group = name.varcat, color = name.varcat), alpha = ci.alpha)
                }
            }            
        }
        reference <- attr(rr,"reference")[,setdiff(names(attr(rr,"reference")),type.residual)]
        gg <- gg + ggplot2::ggtitle(paste0("Reference: ",paste(paste0(names(reference),"=",reference), collapse = ", ")))
        gg <- gg + ggplot2::ylab(paste0("Partial residuals for ",paste(name.var,collapse=", "))) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
        if(plot){print(gg)}
        out <- list(data = gg.data,
                    plot = gg)
    }else{
        outRes <- residuals(x, plot = type, type = type.residual, format = format, size.text = size.text, ...)
        out <- list(data = outRes,
                    plot = attr(outRes,"plot"))
        attr(out$data,"plot") <- NULL
    }
    return(invisible(out))
}

##----------------------------------------------------------------------
### plot.R ends here
