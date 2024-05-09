### scatterplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2023 (09:39) 
## Version: 
## Last-Updated: May  9 2024 (11:32) 
##           By: Brice Ozenne
##     Update #: 705
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * scatterplot (documentation)
##' @title Scatterplot for Continuous Variables
##' @description Produce a matrix of plot for continuous variables: scatterplots, histograms, correlation and missing values.
##' Inspired from the \code{ggpairs} function of the R package GGally.
##' 
##' @param data [data.frame] dataset containing the variables to be displayed.
##' @param formula [formula] formula indicating the variables to be used (outcome~time|id). Long format only.
##' @param columns [character vector] Columns whose numerical values are to be displayed. Wide format only.
##' @param format [character] Is the dataset in the long (\code{"long"}) or wide (\code{"wide"}) format?
##' @param group [character] optional group variable used to color the points, stratify the histogram/density and correlation.
##' @param transform [character or function] optional transformation to be applied on the outcome.
##' @param facet [character] whether to use \code{ggplot:::facet_grid} (\code{"grid"}) or \code{ggh4x::facet_grid2} (\code{"grid2"}).
##' @param alpha.point [numeric] the transparency level used to display the points in the scatterplot.
##' @param type.diag [character] type of graphical display on the diagonal: \code{"boxplot"},  \code{"histogram"}, or \code{"density"}.
##' @param bins [character or numeric vector] algorithm or values or number of values used to create the histogram cells.
##' When using \code{facet="grid2"} and \code{density=TRUE} a character of length two indicating the bandwith and the kernel to be used.
##' See \code{ggplot2::stat_density}.
##' @param position.bar [character] passed to \code{geom_histogram} (argument \code{position}).
##' Only relevant when having multiple groups and using \code{ggh4x::facet_grid2}.
##' @param linewidth.density [numeric,>0] width of the lines on the density plot.
##' @param alpha.area [numeric, 0-1] the transparency level used to display the area under the density curve or histogram.
##' @param method.cor [character] estimator of the correlation. Argument passed to \code{stats::cor}.
##' When \code{NA}, the correlation is not displayed.
##' @param name.cor [character] character used to represent the correlation. By default \code{"r"} but can be changed to \code{"\u03C1"} to display the greek letter \eqn{\rho}.
##' @param size.cor [numeric,>0] size of the font used to display the correlation or information about missing values.
##' @param digits [numeric of length 2] number of digits used to display the correlation or round the percentage of missing values.
##' @param display.NA [0:2 or "only"] Should the number of missing values be displayed. When taking value 2, will also display the percentage of missing values.
##' @param color [character vector] color used to display the values for each group.
##' @param xlim [numeric,>0 or "common"] range of the x-axis.
##' @param ylim [numeric,>0 or "common"] range of the y-axis.
##' @param size.axis [numeric,>0] size of the font used to display the tick labels.
##' @param size.legend [numeric,>0] size of the font used to display the legend. Can have a second element to control the size of the legend key.
##' @param size.facet [numeric,>0] size of the font used to display the facets (row and column names).
##' 
##' @details In the long format, the outcome variable contains the numerical values to be displayed.
##' The time variable will be used to spit outcome and display each split separately or jointly with one other split.
##' The identifier links the outcome values across time.
##' 
##' @return a list of ggplot objects (\code{facet="grid"}) or a ggplot object (\code{facet="grid2"})
##' 
##' @keywords utilities
##' 
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$group <- as.numeric(gastricbypassL$id) %% 3
##' data(gastricbypassW, package = "LMMstar")
##'
##' ## single group (wide or long format)
##' scatterplot(gastricbypassL, formula = weight~time|id)
##' scatterplot(gastricbypassW, columns = paste0("weight",1:4))
##'
##' \dontrun{
##' ## use histogram instead of boxplot
##' scatterplot(gastricbypassL, formula = weight~time|id, type.diag = "hist")
##' scatterplot(gastricbypassL, formula = weight~time|id, type.diag = "hist", bins = 15)
##'
##' ## same scale
##' scatterplot(gastricbypassL, formula = weight~time|id,
##'             xlim = "common", ylim = "common")
##' 
##' ## transform outcome
##' scatterplot(gastricbypassL, formula = weight~time|id, transform = "log")
##'
##' ## handling missing values
##' scatterplot(gastricbypassL, formula = glucagonAUC~time|id)
##'
##' ## coloring per group
##' scatterplot(gastricbypassL, formula = weight~time|id, group = "group")
##' 
##' ## only display NAs
##' scatterplot(gastricbypassL, formula = glucagonAUC~time|id,
##'             display.NA = "only", group = "group")
##' scatterplot(gastricbypassL, formula = glucagonAUC~time|id,
##'             display.NA = "only", group = "group", size.legend = c(15,2))
##' }


## * scatterplot (code)
##' @export
scatterplot <- function(data, formula, columns, format = NULL, group = NULL, transform = NULL,
                        facet = "grid", 
                        alpha.point = 1,
                        type.diag = "boxplot", bins = NULL, position.bar = "identity", linewidth.density = NULL, alpha.area = NULL,
                        method.cor = "pearson", name.cor = "r", size.cor = NULL, digits = c(3,2), display.NA = NULL,
                        color = NULL, xlim = NULL, ylim = NULL, size.axis = NULL, size.legend = NULL, size.facet = NULL){

    ## ** normalize user input
    ## facet
    type.diag <- match.arg(type.diag, c("hist","histogram","density","boxplot"))
    if(type.diag=="histogram"){
        type.diag <- "hist"
    }

    ## facet
    facet <- match.arg(facet, c("grid","grid2"))

    ## display.NA
    if(identical(display.NA,"only")){
        method.cor <- NA
        display.NA <- TRUE
    }

    ## xlim, ylim
    if(is.character(xlim) && !identical(xlim,"common")){
        stop("When a character, argument \'xlim\' should take value \"common\". \n")
    }
    if(is.character(ylim) && !identical(ylim,"common")){
        stop("When a character, argument \'ylim\' should take value \"common\". \n")
    }

    ## format
    if(is.null(format)){
        if(missing(formula) & !missing(columns)){
            format <- "wide"
        }else if(!missing(formula) & missing(columns)){
            format <- "long"
        }else if(missing(formula) & missing(columns)){
            format <- "wide"
        }else{
            stop("Only one of the arguments \'formula\' and \'columns\' should be specified \n",
                 "(respectively for the long and wide format). \n")
        }
    }
    format <- match.arg(format, c("long","wide"))

    ## data in the right long format
    data <- as.data.frame(data)
    if(!missing(columns)){
        if(any(columns %in% names(data) == FALSE)){
            invalid <- columns[columns %in% names(data) == FALSE]
            stop("Argument \'columns\' is inconsistent with argument \'data\'. \n",
                 "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
                 sep = "")
        }
        data <- data[union(columns,group)]
    }

    if(any(names(data) %in% c("XXtime1XX","XXtime2XX"))){
        invalid <- names(data)[names(data) %in% c("XXtime1XX","XXtime2XX")]
        stop("Name \"",paste(names(data), collapse = "\" \""),"\" is used internally. \n",
             "Consider renaming the variable in the dataset. \n",
             sep = "")
    }


    if(format == "wide"){
        if(!missing(formula)){
            warning("Argument \'formula\' is ignored when using the wide format. \n")
        }
        if(any(names(data) %in% c("XXindexXX"))){
            invalid <- names(data)[names(data) %in% c("XXindexXX")]
            stop("Name \"",paste(names(data), collapse = "\" \""),"\" is used internally. \n",
                 "Consider renaming the variable in the dataset. \n",
                 sep = "")
        }
        name.Y <- "outcome"
        name.time <- "time"
        level.time <- setdiff(names(data),group)
        data$XXindexXX <- 1:NROW(data)
        name.id <- "XXindexXX"

        dataL <- stats::reshape(data, direction = "long",
                                varying = level.time, v.names = name.Y, times = level.time,
                                idvar = c(name.id,group))
        dataL[[name.time]] <- factor(dataL[[name.time]], levels = level.time)
    }else{
        dataL <- data
        if(!missing(columns)){
            warning("Argument \'columns\' is ignored when using the wide format. \n")
        }

        detail.formula <- formula2var(formula)
        name.all <- detail.formula$vars$all
        if(any(name.all %in% names(dataL) == FALSE)){
            invalid <- name.all[name.all %in% names(dataL) == FALSE]
            stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
                 "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
                 sep = "")
        }
        name.Y <- detail.formula$var$response
        n.Y <- length(name.Y)
        if(n.Y!=1){
            stop("Wrong specification of argument \'formula\'. \n",
                 "There need to be exactly one variable in the left hand side of the formula. \n")
        }

        name.id <- detail.formula$var$cluster
        if(length(detail.formula$var$ranef)>0){
            stop("Wrong specification of argument \'formula\'. \n",
                 "Should be something like Y ~ time|id. \n")
        }else if(length(name.id)==0){
            dataL$XXindexXX <- 1:NROW(dataL)
            name.id <- "XXindexXX"
            name.time <- detail.formula$var$regressor
        }else{
            name.time <- detail.formula$var$time
        }
        
        if(length(name.time)==0){
            stop("Wrong specification of argument \'formula\'. \n",
                 "There should be at least one variable on the right hand side of the formula before the symbol |. \n")
        }
        if(length(name.all) == (length(name.Y) + length(name.time))){
            stop("Wrong specification of argument \'formula\'. \n",
                 "There should be exactly one variable on the right hand side of the formula after the symbol |. \n")
        }
        if(name.id %in% name.time){
            stop("Wrong specification of argument \'formula\'. \n",
                 "Should be outcome ~ time | id, where the time variable should not include the id variable. \n")
        }
        if(length(name.time)>1){
            if(any(names(dataL) %in% c("XXtimeXX"))){
                invalid <- names(dataL)[names(dataL) %in% c("XXtimeXX")]
                stop("Name \"",paste(names(dataL), collapse = "\" \""),"\" is used internally. \n",
                     "Consider renaming the variable in the dataset. \n",
                     sep = "")
            }
            dataL$XXtimeXX <- nlme::collapse(dataL[name.time], as.factor = TRUE)
            name.time <- "XXtimeXX"
        }
        if(!is.factor(dataL[[name.time]])){
            dataL[[name.time]] <- factor(dataL[[name.time]])
        }
        level.time <- levels(dataL[[name.time]])
    }
    if(!is.null(transform)){
        dataL[[name.Y]] <- do.call(transform, list(dataL[[name.Y]]))
    }
    if(!is.null(group)){
        if(!is.factor(dataL[[group]])){
            dataL[[group]] <- factor(dataL[[group]])
        }
        level.group <- levels(dataL[[group]])        
    }

    ## bins
    if(is.null(bins)){
        if(type.diag == "hist"){
            bins <- "Sturges"
        }else if(type.diag == "density"){
            bins <- c("nrd0","gaussian")
        }
    }

    ## alpha.area
    if(is.null(alpha.area)){
        if(type.diag == "density" || type.diag == "boxplot"){
            alpha.area <- 0.3
        }else if(position.bar=="dodge" || is.null(group)){
            alpha.area  <- 1
        }else{
            alpha.area  <- 0.7
        }
    }

    ## ** matrix format
    ls.dataL <- split(dataL, dataL[[name.time]])
    grid.time <- expand.grid(XXtime1XX = level.time,
                             XXtime2XX = level.time)
    ls.dataGrid <- lapply(1:NROW(grid.time), function(iGrid){ ## iGrid <- 2
        iData1 <- ls.dataL[[grid.time[iGrid,"XXtime1XX"]]][,c(name.id,name.Y,name.time, group)]
        iData2 <- ls.dataL[[grid.time[iGrid,"XXtime2XX"]]][,c(name.id,name.Y,name.time, group)]
        if(is.null(group)){
            names(iData1) <- c("id","outcome1","time1")
            names(iData2) <- c("id","outcome2","time2")
            iDf <- merge(iData1,iData2,by="id",all=TRUE)
        }else{
            names(iData1) <- c("id","outcome1","time1","group")
            names(iData2) <- c("id","outcome2","time2","group")
            iDf <- merge(iData1,iData2,by=c("id","group"),all=TRUE)
        }
        if(any(is.na(iDf$time1))){
            iDf$time1[is.na(iDf$time1)] <- grid.time[iGrid,"XXtime1XX"]
        }
        if(any(is.na(iDf$time2))){
            iDf$time2[is.na(iDf$time2)] <- grid.time[iGrid,"XXtime1XX"]
        }
        return(iDf)
    })
    dataGrid <- do.call(rbind,ls.dataGrid)
    dataGrid$position <- as.character(NA)
    dataGrid$position[dataGrid$time1==dataGrid$time2] <- "diag"
    dataGrid$position[as.numeric(dataGrid$time1)>as.numeric(dataGrid$time2)] <- "lower"
    dataGrid$position[as.numeric(dataGrid$time1)<as.numeric(dataGrid$time2)] <- "upper"
    dataGrid$time <- as.factor(paste(dataGrid$time1,dataGrid$time2, sep ="."))

    ## ** prepare outcome
    n.time <- length(level.time)
    if(n.time<2){
        stop("The scatterplot function does not handle a single variable or timepoint. \n")
    }else if(n.time==2){
        dataGrid.lower <- dataGrid[dataGrid$position == "lower",]
    }else{
        dataGrid.lower <- dataGrid[dataGrid$position == "lower",]
    }
        
    
    ## ** prepare histograms
    dataGrid.diag <- dataGrid[dataGrid$position=="diag",]

    ## ** prepare correlation
    if(n.time>2){
        dataGrid.upper <- dataGrid[dataGrid$position=="upper",]
        dataCor <- do.call(rbind, by(dataGrid.upper, droplevels(dataGrid.upper$time), function(iData){
            if(is.null(group)){
                iDF <- data.frame(outcome1 = 1,
                                  outcome2 = mean(range(iData$outcome2,na.rm=TRUE)),
                                  time1 = unique(iData$time1),
                                  time2 = unique(iData$time2),
                                  cor = if(is.na(method.cor)){NA}else{stats::cor(iData$outcome1, iData$outcome2, method = method.cor, use = "pairwise")},
                                  n = NROW(is.na(iData)),
                                  n.NNA = sum(rowSums(is.na(iData))==0)
                                  )

            }else{
                iCor <- by(iData, iData$group, function(iiData){
                    data.frame(group = unique(iiData$group),
                               cor = if(is.na(method.cor)){NA}else{stats::cor(iiData$outcome1, iiData$outcome2, method = method.cor, use = "pairwise")},
                               n = NROW(is.na(iiData)),
                               n.NNA = sum(rowSums(is.na(iiData))==0)
                               )
                })
                iDF <- data.frame(outcome1 = 1:length(iCor),
                                  outcome2 = mean(range(iData$outcome2,na.rm=TRUE)),
                                  time1 = unique(iData$time1),
                                  time2 = unique(iData$time2),
                                  do.call(rbind,iCor))
            }
            return(iDF)
        }))
    }

    n.NA <- dataCor$n-dataCor$n.NNA
    pc.NA <- 100*(1-dataCor$n.NNA/dataCor$n)
    if(is.null(display.NA)){
        if(all(n.NA==0)){
            display.NA <- FALSE
        }else{
            display.NA <- TRUE
        }
    }

    if(!is.na(method.cor)){
        dataCor$label <- paste0(name.cor,"=",round(dataCor$cor, digits[1]))
        ## dataCor$label <- paste0("\u03C1=",round(dataCor$cor, digits[1])) ## NOT DONE BY DEFAULT DUE TO ERROR RUNNING CRAN CHECK ON LINUX
        if(display.NA==1){
            dataCor$label <- paste0(dataCor$label, "; ",n.NA," NA")
        }else if(display.NA>1){
            dataCor$label <- paste0(dataCor$label, "; ",n.NA," NA (",round(100*(1-dataCor$n.NNA/dataCor$n), digits[2]),"%)")
        }
    }else if(display.NA==1){
        dataCor$label <- paste0(n.NA," NA")
    }else if(display.NA>1){
        dataCor$label <- paste0(n.NA," NA (",round(100*(1-dataCor$n.NNA/dataCor$n), digits[2]),"%)")
    }

    ## ** graphical display
    ## size.cor
    if(is.null(size.cor)){
        if(display.NA && !is.na(method.cor)){
            size.cor <- 7 
        }else{
            size.cor <- 10
        }
    }
    if(facet=="grid"){
        gg <- .ggscatterplot(dataGrid.diag = dataGrid.diag, dataGrid.lower = dataGrid.lower, dataCor = dataCor,
                             bins = bins, level.time = level.time, n.time = n.time, group = group, 
                             alpha.point = alpha.point, position.bar = position.bar, linewidth.density = linewidth.density,
                             type.diag = type.diag, alpha.area = alpha.area, method.cor = method.cor, size.cor = size.cor, display.NA = display.NA, 
                             color = color, xlim = xlim, ylim = ylim, size.axis = size.axis, size.legend = size.legend, size.facet = size.facet)
        return(invisible(gg))
    }else if(facet == "grid2"){
        gg <- .ggscatterplot2(dataGrid.diag = dataGrid.diag, dataGrid.lower = dataGrid.lower, dataCor = dataCor,
                              bins = bins, n.time = n.time, group = group, 
                              alpha.point = alpha.point, position.bar = position.bar, linewidth.density = linewidth.density,
                              type.diag = type.diag, alpha.area = alpha.area, method.cor = method.cor, size.cor = size.cor, display.NA = display.NA, 
                              color = color, xlim = xlim, ylim = ylim, size.axis = size.axis, size.legend = size.legend, size.facet = size.facet)
        return(gg)
    }

}

## * .ggscatterplot
.ggscatterplot <- function(dataGrid.diag, dataGrid.lower, dataCor,
                           bins, level.time, n.time, group, 
                           alpha.point, position.bar, linewidth.density, type.diag, alpha.area, method.cor, size.cor, display.NA, 
                           color, xlim, ylim, size.axis, size.legend, size.facet){

    if(is.null(size.legend) || length(size.legend) %in% 1:2){
        ratio.legend <- 1/4
    }else{
        ratio.facet <- ratio.legend[3]
    }
    if(is.null(size.facet) || length(size.facet)==1){
        ratio.facet <- 1/8
    }else{
        ratio.facet <- size.facet[2]
    }

    ## ** create each graphical display
    grid <- expand.grid(time1 = 1:n.time, time2 =  1:n.time)
    n.grid <- NROW(grid)
    vec.plot <- vector(mode = "list", length = n.grid)
    if(is.null(group)){
        dataGrid.diag$group <- "1"
        dataGrid.lower$group <- "1"
        dataCor$group <- "1"
        if(is.null(color)){
            color <- "black"
        }
    }


    for(iGrid in 1:n.grid){ ## iGrid <- 5

        iT1 <- grid[iGrid,1]
        iT2 <- grid[iGrid,2]
        iTime1 <- level.time[iT1]
        iTime2 <- level.time[iT2]
        vec.plot[[iGrid]] <- ggplot2::ggplot()

        if(iT1==iT2){
            ## histogram or density plot or boxplot
            iData <- dataGrid.diag[dataGrid.diag$time1 == iTime1,,drop=FALSE]

            if(type.diag == "density"){
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::geom_density(data = iData,
                                                                               mapping = ggplot2::aes(x = .data$outcome1, fill = .data$group, color = .data$group),
                                                                               bw = bins[1], kernel = bins[2], alpha = alpha.area,
                                                                               outline.type = "full", linewidth = linewidth.density, show.legend = FALSE)
                
            }else if(type.diag == "hist"){
                if(is.character(bins)){
                    iBreaks <- graphics::hist(iData$outcome1, breaks = bins, plot = FALSE)$breaks
                    iBins <- NULL
                }else if(length(bins)>1){
                    iBreaks <- bins
                    iBins <- NULL
                }else{
                    iBreaks <- NULL
                    iBins <- bins
                }
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::geom_histogram(data = iData,
                                                                                 mapping = ggplot2::aes(x = .data$outcome1, fill = .data$group, color = .data$group),
                                                                                 bins = iBins, breaks = iBreaks, alpha = alpha.area, position = position.bar, show.legend = FALSE)

            }else if(type.diag == "boxplot"){
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::geom_boxplot(data = iData,
                                                                               mapping = ggplot2::aes(x = .data$outcome1, fill = .data$group, color = .data$group),
                                                                               alpha = alpha.area, show.legend = FALSE)
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::theme(axis.text.y = ggplot2::element_blank(),  #remove y axis labels
                                                                        axis.ticks.y = ggplot2::element_blank()  #remove y axis ticks
                                                                        )
            }


        }else if(iT1>iT2){
            ## scatterplot
            iData <- dataGrid.lower[dataGrid.lower$time1 == iTime1 & dataGrid.lower$time2 == iTime2,,drop=FALSE]

            vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::geom_point(data = iData,
                                                                         mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1, color = .data$group, shape = .data$group),
                                                                         alpha = alpha.point, show.legend = FALSE)

        }else if(iT1<iT2){
            ## correlation
            iData <- dataCor[dataCor$time1 == iTime1 & dataCor$time2 == iTime2,,drop=FALSE]
            
            if(!is.na(method.cor)){
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::geom_text(data = iData,
                                                                            mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                                                   label = .data$label, color = .data$group),
                                                                            size = size.cor, show.legend = FALSE, vjust = "inward")

            }else{
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::geom_text(data = iData,
                                                                            mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                                                   label = .data$label, color = .data$group),
                                                                            size = size.cor, show.legend = FALSE, vjust = "inward")
            }
            if(!is.null(group)){
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::coord_cartesian(ylim = c(min(iData$outcome1)-0.5,max(iData$outcome1)+0.5))
            }

            vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::theme_void() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), #remove x axis labels
                                                                                            axis.ticks.x = ggplot2::element_blank(), #remove x axis ticks
                                                                                            axis.text.y = ggplot2::element_blank(),  #remove y axis labels
                                                                                            axis.ticks.y = ggplot2::element_blank()  #remove y axis ticks
                                                                                            )

        }

        ## ** add to the graphical display
        vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::labs(x = "", y = "") + ggplot2::theme(plot.margin = ggplot2::margin(0.1,0.1,0.1,0.1, "cm"))
        if(!is.null(color)){
            vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::scale_color_manual(values = color) + ggplot2::scale_fill_manual(values = color)           
        }
        if(!is.null(size.axis)){
            vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::theme(text = ggplot2::element_text(size=size.axis))
        }

    }
    
    ## ** generate graphical output
    ## based on https://stackoverflow.com/questions/46893560/manual-facet-position-in-ggplot2
    if(!is.null(group)){ ## add title
        if(is.null(color)){
            color <- unique(ggplot2::ggplot_build(vec.plot[[2]])$data[[1]][["fill"]])
        }
        shape <- unique(ggplot2::ggplot_build(vec.plot[[2]])$data[[1]][["shape"]])
    }

    if(identical(xlim,"common")){
        index.point <- which(grid[,1]>grid[,2])
        ls.x <- lapply(index.point, function(iGG){ggplot2::ggplot_build(vec.plot[[iGG]])$data[[1]][["x"]]})
        xlim <- range(unlist(ls.x))               

        index.hist <- which(grid[,1]==grid[,2])
        if(type.diag == "boxplot"){
            ls.x.hist <- lapply(index.hist, function(iGG){
                iData <- ggplot2::ggplot_build(vec.plot[[iGG]])$data[[1]]
                c(iData$xmin,iData$xmax)
            })
        }else{
            ls.x.hist <- lapply(index.hist, function(iGG){ggplot2::ggplot_build(vec.plot[[iGG]])$data[[1]][["x"]]})
        }
        xlim.hist <- range(unlist(ls.x.hist))
    }else{
        xlim.hist <- xlim
    }
    if(identical(ylim,"common")){
        index.point <- which(grid[,1]>grid[,2])
        ls.y <- lapply(index.point, function(iGG){ggplot2::ggplot_build(vec.plot[[iGG]])$data[[1]][["y"]]})
        ylim <- range(unlist(ls.y))

        if(type.diag == "boxplot"){
            ylim.hist <- NULL
        }else{
            index.hist <- which(grid[,1]==grid[,2])
            ls.y.hist <- lapply(index.hist, function(iGG){ggplot2::ggplot_build(vec.plot[[iGG]])$data[[1]][["y"]]})
            ylim.hist <- range(unlist(ls.y.hist))
        }
    }else{
        ylim.hist <- ylim
    }

    ## Set up the page
    grid::grid.newpage()
    if(!is.null(group)){
        panelVP <- grid::viewport(layout = grid::grid.layout(nrow = n.time+1,
                                                             ncol = n.time+2,
                                                             width = ggplot2::unit(c(ratio.legend,rep(1, n.time),ratio.legend), "null"),
                                                             height = ggplot2::unit(c(rep(1, n.time),ratio.legend), "null")))
        
    }else{
        panelVP <- grid::viewport(layout = grid::grid.layout(nrow = n.time+1,
                                                             ncol = n.time+1,
                                                             width = ggplot2::unit(c(ratio.facet,rep(1, n.time)), "null"),
                                                             height = ggplot2::unit(c(rep(1, n.time),ratio.facet), "null")))
    }
    grid::pushViewport(panelVP)
    ## grid::showViewport()
    for(iTime in 1:n.time){
        grid::grid.rect(gp = grid::gpar(fill="grey"),
                        vp = grid::viewport(layout.pos.row = iTime, layout.pos.col = 1))
        grid::grid.text(level.time[iTime],
                        vp = grid::viewport(layout.pos.row = iTime, layout.pos.col = 1),
                        rot = 90,
                        gp = grid::gpar(fontsize = size.facet[1]))
        grid::grid.rect(gp = grid::gpar(fill="grey"),
                        vp = grid::viewport(layout.pos.row = n.time+1, layout.pos.col = iTime+1))
        grid::grid.text(level.time[iTime],
                        vp = grid::viewport(layout.pos.row = n.time+1, layout.pos.col = iTime+1),
                        rot = 0,
                        gp = grid::gpar(fontsize = size.facet[1]))
    }
    
    ## display plots
    for(iGrid in 1:n.grid){ ## iGrid <- 5
        if(!is.null(xlim) || !is.null(ylim)){
            if(grid[iGrid,1] == grid[iGrid,2]){
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::coord_cartesian(xlim = xlim.hist, ylim = ylim.hist)
            }else if(grid[iGrid,1] > grid[iGrid,2]){
                vec.plot[[iGrid]] <- vec.plot[[iGrid]] + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
            }
        }
        print(vec.plot[[iGrid]],
              vp = grid::viewport(layout.pos.row = grid[iGrid,1], layout.pos.col = grid[iGrid,2]+1))
    }

    ## display legend
    if(!is.null(group)){
        level.group <- levels(dataGrid.lower$group)
        df.legend <- data.frame(x = 1:length(level.group), y = 0, shape = level.group, color = level.group)
        gg.legend <- ggplot2::ggplot(df.legend, ggplot2::aes(x = .data$x, y = .data$y, color = .data$color, shape = .data$color))
        gg.legend <- gg.legend + ggplot2::geom_point() + ggplot2::labs(color = group, shape = group)
        if(!is.null(size.legend)){
            gg.legend <- gg.legend + ggplot2::theme(text = ggplot2::element_text(size=size.legend[1]))
            if(length(size.legend)>1){
                gg.legend <- gg.legend + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=size.legend[2])),
                                                         fill = ggplot2::guide_legend(override.aes = list(size=size.legend[2])))
            }
        }

        ## from https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
        gtable.legend <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(gg.legend)) 
        index.legend <- which(sapply(gtable.legend$grobs, function(x) x$name) == "guide-box") 
        
        print(.grob2ggplot2(gtable.legend$grobs[[index.legend]]),
              vp = grid::viewport(layout.pos.row = 1:n.time, layout.pos.col = n.time+2))

    }

    ## ** export
    return(vec.plot)
}

## * .ggscatterplot2
.ggscatterplot2 <- function(dataGrid.diag, dataGrid.lower, dataCor,
                            bins, n.time, group, 
                            alpha.point, position.bar, linewidth.density, type.diag, alpha.area, method.cor, size.cor, display.NA, 
                            color, xlim, ylim, size.axis, size.legend, size.facet){

    requireNamespace("ggh4x")

    ## ** prepare
    if(is.null(group)){
        dataGrid.diag$group <- "1"
        dataGrid.lower$group <- "1"
        dataCor$group <- "1"
        if(is.null(color)){
            color <- "black"
        }
    }

    ## ** graphical display
    gg <- ggplot2::ggplot()
    gg <- gg + ggh4x::facet_grid2(time1~time2, scales = "free", labeller = "label_value", independent = "y") + ggplot2::labs(x = "", y = "")
    gg <- gg + ggplot2::geom_point(data = dataGrid.lower,
                                   mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1, color = .data$group, shape = .data$group),
                                   alpha = alpha.point, show.legend = !is.null(group))

    if(type.diag == "density"){
        gg <- gg + ggplot2::geom_density(data = dataGrid.diag,
                                         mapping = ggplot2::aes(x = .data$outcome1, fill = .data$group, color = .data$group),
                                         bw = bins[1], kernel = bins[2], alpha = alpha.area,
                                         outline.type = "full", linewidth = linewidth.density, show.legend = !is.null(group))
    }else if(type.diag == "hist"){
        if(is.character(bins)){
            iBreaks <- graphics::hist(dataGrid.diag$outcome1, breaks = bins, plot = FALSE)$breaks
            iBins <- NULL
        }else if(length(bins)>1){
            iBreaks <- bins
            iBins <- NULL
        }else{
            iBreaks <- NULL
            iBins <- bins
        }

        gg <- gg + ggplot2::geom_histogram(data = dataGrid.diag,
                                           mapping = ggplot2::aes(x = .data$outcome1, fill = .data$group),
                                           bins = iBins, breaks = iBreaks, alpha = alpha.area, position = position.bar, show.legend = !is.null(group))
    }else if(type.diag == "boxplot"){
        gg <- gg + ggplot2::geom_boxplot(data = dataGrid.diag,
                                         mapping = ggplot2::aes(x = .data$outcome1, fill = .data$group, color = .data$group),
                                         alpha = alpha.area, show.legend = !is.null(group))
    }

    if(!is.na(method.cor)){
        gg <- gg + ggplot2::geom_text(data = dataCor,
                                      mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                             label = .data$label, color = .data$group),
                                      size = size.cor, show.legend = FALSE, vjust = "inward")
    }else{
        gg <- gg + ggplot2::geom_text(data = dataCor,
                                      mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                             label = .data$label, color = .data$group),
                                      size = size.cor, show.legend = FALSE, vjust = "inward")
    }
            
    ## ** export
    if(!is.null(xlim)){
        if(identical(xlim,"common")){
            xlim <- range(unlist(lapply(ggplot2::ggplot_build(gg)$data, function(iGG){iGG$x})))
        }        
        gg <- gg + ggplot2::coord_cartesian(xlim = xlim)
    }
    if(!is.null(size.axis)){
        gg <- gg + ggplot2::theme(axis.text = ggplot2::element_text(size=size.axis[1]))
    }
    if(!is.null(size.facet)){
        gg <- gg + ggplot2::theme(strip.text.x = ggplot2::element_text(size = size.facet[1]),
                                  strip.text.y = ggplot2::element_text(size = size.facet[1]))
    }
    if(!is.null(size.legend)){
        gg <- gg + ggplot2::theme(legend.text = ggplot2::element_text(size = size.legend[1]),
                                  legend.title = ggplot2::element_text(size = size.legend[1]))
        if(length(size.legend)>1){
            gg <- gg + ggplot2::theme(legend.key.size = ggplot2::unit(size.legend[2], 'cm'))
        }
    }
    if(!is.null(color)){
        gg <- gg + ggplot2::scale_color_manual(values = color) + ggplot2::scale_fill_manual(values = color)           
    }
    return(gg)
}

## * .grob2ggplot2
## from ggplotify:::as.ggplot_internal
.grob2ggplot2 <- function(plot, scale = 1, hjust = 0, vjust = 0){
    ymin <- xmin <- 1 - scale
    xmax <- ymax <- scale
    gg <- ggplot2::ggplot(data.frame(x = 0:1, y = 0:1), ggplot2::aes(x = .data$x, y = .data$y))
    gg <- gg + ggplot2::geom_blank() + ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0))
    gg <- gg + ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
    gg <- gg + ggplot2::annotation_custom(plot, xmin = xmin + hjust, xmax = xmax + hjust, ymin = ymin + vjust, ymax = ymax + vjust)
    gg <- gg + ggplot2::theme_void()
    return(gg)
}
##----------------------------------------------------------------------
### scatterplot.R ends here
