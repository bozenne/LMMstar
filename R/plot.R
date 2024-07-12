### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 20 2021 (11:00) 
## Version: 
## Last-Updated: jul 11 2024 (11:21) 
##           By: Brice Ozenne
##     Update #: 141
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * plot.correlate (code)
##' @describeIn autoplot.correlate Graphical Display For Correlation Matrix
##' @export
plot.correlate <- function(x, ...){
    out <- autoplot.correlate(x, ...)
    if(!is.null(out$plot)){
        ## if ggplot then the graph is stored in out and should be displayed here
        ## otherwise the graph has already been displayed but could not be stored in the object (e.g. when using qqtest)
        print(out$plot)
    }
    return(invisible(out))

}

## * plot.lmm (code)
##' @describeIn autoplot.lmm Graphical Display For Linear Mixed Models
##' @export
plot.lmm <- function(x, ...){
    out <- autoplot.lmm(x, ...)
    if(!is.null(out$plot)){
        ## if ggplot then the graph is stored in out and should be displayed here
        ## otherwise the graph has already been displayed but could not be stored in the object (e.g. when using qqtest)
        print(out$plot)
    }
    return(invisible(out))

}

## * plot.mlmm (code)
##' @describeIn autoplot.mlmm Graphical Display For Multiple Linear Mixed Model
##' @export
plot.mlmm <- function(x, facet_nrow = NULL, facet_ncol = NULL, labeller = "label_value", ...){

    out <- autoplot.mlmm(x, ...)

    if(inherits(out$plot,"ggplot")){
        print(out$plot)
    }else if(is.list(out$plot) && all(sapply(out$plot, function(iP){inherits(iP,"ggplot")}))){

        if(!is.character(labeller) || length(labeller)>1){
            stop("Argument \'labeller\' should be a character (with length 1). \n")
        }
        labeller <- match.arg(labeller, c("label_value","label_both"))

        ## count plots
        n.plot <- length(out$plot)
        if(labeller == "label_value"){
            names.plot <- names(out$plot)
        }else if(labeller == "label_both"){
            names.plot <- paste0(x$object$by,": ", names(out$plot))
        }

        ## Set up the page
        grid::grid.newpage()
        
        ## Define the layout
        if(!is.null(facet_nrow) & !is.null(facet_ncol)){
            panelVP <- grid::grid.layout(nrow = facet_nrow,
                                         ncol = facet_ncol)
        }else if(!is.null(facet_nrow)){
            panelVP <- grid::grid.layout(nrow = facet_nrow,
                                         ncol = ceiling(n.plot/facet_nrow))
        }else if(!is.null(facet_ncol)){
            panelVP <- grid::grid.layout(nrow = ceiling(n.plot/facet_ncol),
                                         ncol = facet_ncol)
        }else{
            facet_nrow <- round(sqrt(n.plot))
            panelVP <- grid::grid.layout(nrow = facet_nrow,
                                         ncol = ceiling(n.plot/facet_nrow))
        }
        grid.plot <- expand.grid(row = 1:panelVP$nrow, col =  1:panelVP$ncol)

        grid::pushViewport(grid::viewport(layout = panelVP))
        ## grid::showViewport()

        ## Display the plot in the layout
        for(iPlot in 1:n.plot){ ## iPlot <- 1
            print(out$plot[[iPlot]] + ggplot2::ggtitle(names.plot[[iPlot]]), vp = grid::viewport(layout.pos.row = grid.plot[iPlot,"row"], layout.pos.col = grid.plot[iPlot,"col"]))
        }
    
    
    }
    
    return(invisible(out))
}
## * plot.partialCor (code)
##' @describeIn autoplot.partialCor Graphical Display For Partial Correlation
##' @export
plot.partialCor <- function(x, ...){
    out <- autoplot.partialCor(x, ...)
    print(out$plot)
    return(invisible(out))

}

## * plot.profil_lmm (code)
##' @describeIn autoplot.profile_lmm Display Contour of the log-Likelihood
##' @export
plot.profile_lmm <- function(x, ...){
    out <- autoplot.profile_lmm(x, ...)
    print(out$plot)
    return(invisible(out))

}

## * plot.residuals_lmm (code)
##' @describeIn autoplot.residuals_lmm Graphical Display of the Residuals
##' @export
plot.residuals_lmm <- function(x, ...){
    out <- autoplot.residuals_lmm(x, ...)
    print(out$plot)
    return(invisible(out))

}

## * plot.summarize (code)
##' @describeIn autoplot.summarize Graphical Display of Missing Data Pattern
##' @export
plot.summarize <- function(x, ...){
    out <- autoplot.summarize(x, ...)
    print(out$plot)
    return(invisible(out))
}

## * plot.summarizeNA (code)
##' @describeIn autoplot.summarizeNA Graphical Display of Missing Data Pattern
##' @export
plot.summarizeNA <- function(x, ...){
    out <- autoplot.summarizeNA(x, ...)
    print(out$plot)
    return(invisible(out))
}

## * plot.Wald_lmm (code)
##' @describeIn autoplot.Wald_lmm Graphical Display For Linear Hypothesis Test
##' @export
plot.Wald_lmm <- function(x, ...){
    out <- autoplot.Wald_lmm(x, ...)
    print(out$plot)
    return(invisible(out))
}

##----------------------------------------------------------------------
### plot.R ends here
