### scatterplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2023 (09:39) 
## Version: 
## Last-Updated: feb 24 2023 (18:34) 
##           By: Brice Ozenne
##     Update #: 371
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
##' @param columns [character vector] Columns whose numerical values are to be displayed. Wide format only.
##' @param formula [formula] formula indicating the variables to be used (outcome~time|id). Long format only.
##' @param format [character] Is the dataset in the long (\code{"long"}) or wide (\code{"wide"}) format?
##' @param group [character] optional group variable used to color the points, stratify the histogram/density and correlation.
##' @param transform [character or function] optional transformation to be applied on the outcome.
##' @param facet [character] whether to use \code{ggplot:::facet_grid} (\code{"grid"}) or \code{ggh4x::facet_grid2} (\code{"grid2"}).
##' @param labeller [character] passed to \code{ggplot2::facet_grid} to modify the strip labels. 
##' @param alpha.point [numeric] the transparency level used to display the points in the scatterplot.
##' @param breaks [character or numeric vector] algorithm or values used to create the histogram cells.
##' When using \code{facet="grid2"} and \code{density=TRUE} a character of length two indicating the bandwith and the kernel to be used.
##' See \code{ggplot2::stat_density}.
##' @param position.bar [character] passed to \code{geom_histogram} (argument \code{position}).
##' Only relevant when having multiple groups and using \code{ggh4x::facet_grid2}.
##' @param size.bar [numeric,>0] width of the bars of the histogram.
##' @param density [logical] should the density be displayed instead of an histogram.
##' @param alpha.area [numeric, 0-1] the transparency level used to display the area under the density curve or histogram.
##' @param method.cor [character] estimator of the correlation. Argument passed to \code{stats::cor}.
##' When \code{NA}, the correlation is not displayed.
##' @param size.text [numeric,>0] size of the font used to display the correlation or information about missing values.
##' @param digits [numeric of length 2] number of digits used to display the correlation or round the percentage of missing values.
##' 
##' @details In the long format, the outcome variable contains the numerical values to be displayed.
##' The time variable will be used to spit outcome and display each split separately or jointly with one other split.
##' The identifier links the outcome values across time.
##' 
##' @return a ggplot object
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
##' ## tune histogram
##' scatterplot(gastricbypassL, formula = weight~time|id, breaks = 15)
##' 
##' ## transform outcome
##' scatterplot(gastricbypassL, formula = weight~time|id, transform = "log")
##'
##' ## handling missing values
##' scatterplot(gastricbypassL, formula = glucagonAUC~time|id)
##'
##' ## coloring per group
##' gg <- scatterplot(gastricbypassL, formula = weight~time|id, group = "group")
##' gg
##' if(require(ggplot2)){
##' gg + scale_color_manual(values = c("orange","blue","purple"))
##' }
##' 
##' ## only display percentage of NAs
##' scatterplot(gastricbypassL, formula = glucagonAUC~time|id, method.cor = NA)
##' scatterplot(gastricbypassL, formula = glucagonAUC~time|id, method.cor = NA,
##' group = "group", size.text = 5)
##' 
##' ## simple scatterplot
##' scatterplot(gastricbypassW, columns = c("weight2","glucagonAUC1"))
##'


## * scatterplot (code)
##' @export
scatterplot <- function(data, formula, columns, format = NULL, group = NULL, transform = NULL,
                        facet = NULL, labeller = "label_value", alpha.point = 1,
                        breaks = NULL, position.bar = "identity", size.bar = NULL, density = FALSE, alpha.area = NULL,
                        method.cor = "pearson", size.text = 10, digits = c(3,2)){

    ## ** normalize user input
    ## package
    if(is.null(facet)){
        test <- try(requireNamespace("ggh4x"), silent = TRUE)
        if(test==FALSE){
            message("Use facet_grid from ggplot2. Consider installing the package ggh4x for a nice graphical display. \n")
            facet <- "grid"
        }else{
            facet <- "grid2"
        }
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

        name.all <- all.vars(formula)
        if(any(name.all %in% names(dataL) == FALSE)){
            invalid <- name.all[name.all %in% names(dataL) == FALSE]
            stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
                 "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
                 sep = "")
        }
        name.Y <- lhs.vars(formula)
        n.Y <- length(name.Y)
        if(n.Y!=1){
            stop("Wrong specification of argument \'formula\'. \n",
                 "There need to be exactly one variable in the left hand side of the formula. \n")
        }
        index.bar <- grep("|",deparse(formula), fixed = TRUE)
        if(length(index.bar)>1){
            stop("Wrong specification of argument \'formula\'. \n",
                 "There should at most one symbol |. \n")
        }else if(length(index.bar)==1){
            formula.split <- strsplit(split = "|",deparse(formula),fixed=TRUE)
            formula2 <- stats::as.formula(formula.split[[1]][1])
            name.time <- rhs.vars(formula2)            
            name.id <- trimws(formula.split[[1]][2], which = "both")
        }else{
            dataL$XXindexXX <- 1:NROW(dataL)
            name.id <- "XXindexXX"
            name.time <- rhs.vars(formula)
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
            dataL$XXtimeXX <- interaction(dataL[name.time], drop = TRUE)
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
        if(is.null(size.bar)){
            if(facet == "grid"){
                size.bar <- 5/length(level.group)
            }else if(facet == "grid2"){
                size.bar <- 1/length(level.group)
            }
        }
    }else if(is.null(size.bar)){
        if(facet == "grid"){
            size.bar <- 5
        }else if(facet == "grid2"){
            size.bar <- 1
        }
    }

    ## facet
    facet <- match.arg(facet, c("grid","grid2"))

    ## breaks
    if(is.null(breaks)){
        breaks.hist <- "Sturges"

        if(facet == "grid2" && density){
            breaks <- c("nrd0","gaussian")           
        }

    }else{
        breaks.hist <- breaks

    }

    ## alpha.area
    if(is.null(alpha.area)){
        if(density){
            alpha.area <- 0.3
        }else if(facet == "grid"){
            alpha.area  <- 1
        }else if(facet == "grid2"){
            if(position.bar=="dodge" || is.null(group)){
                alpha.area  <- 1
            }else{
                alpha.area  <- 0.7
            }
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
    dataGrid$time <- interaction(dataGrid$time1,dataGrid$time2)

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
    if(n.time>2){
        dataHist <- do.call(rbind,by(dataGrid.diag, dataGrid.diag$time1, function(iData){ ## iData <- dataGrid.diag[dataGrid.diag$time1=="12wks",]
            iHist <- graphics::hist(iData$outcome1, breaks = breaks.hist, plot = FALSE)
            if(is.null(group)){
                iDf <- data.frame(start = iHist$breaks[-length(iHist$breaks)],
                                  mids = iHist$mids,
                                  stop = iHist$breaks[-1],
                                  counts = iHist$counts,
                                  density = iHist$density)            
            }else{
                iDf <- do.call(rbind,by(iData,droplevels(iData$group),function(iiData){ ## iiData <- iData[iData$group=="Placebo",]
                    iiHist <- graphics::hist(iiData$outcome1, breaks = iHist$breaks, plot = FALSE)
                    iiDf <- data.frame(start = iiHist$breaks[-length(iiHist$breaks)],
                                       stop = iiHist$breaks[-1],
                                       counts = iiHist$counts,
                                       density = iiHist$density,
                                       group = unique(iiData$group))
                    iiDf$mids <- sapply(1:NROW(iiDf), function(iRow){
                        seq(from = iiDf$start[iRow], to = iiDf$stop[iRow], length.out = length(level.group)+2)[as.numeric(unique(iiData$group))+1]
                    })
                    return(iiDf)
                }))
            }

            iDf$time1 <- unique(iData$time1)
            iDf$time2 <- unique(iData$time2)
            return(iDf)

        }))
            rownames(dataHist) <- NULL
            dataHist$counts.norm <- dataHist$counts
            dataHist$density.norm <- dataHist$density
            dataHist$base <- 0
            for(iTime in level.time[-1]){ ## iTime <- "12wks"
                iRange.outcome <- range(dataGrid[dataGrid$time2 == iTime,"outcome2"], na.rm = TRUE)
                iRange.counts <- range(dataHist[dataHist$time1 == iTime,"counts"])
                iRange.density <- range(dataHist[dataHist$time1 == iTime,"density"])
                dataHist[dataHist$time1 == iTime,"base"] <- iRange.outcome[1]

                iFactor <- diff(iRange.outcome)/diff(iRange.counts)
                dataHist[dataHist$time1 == iTime,"counts.norm"] <- iRange.outcome[1] + dataHist[dataHist$time1 == iTime,"counts"]*iFactor

                iFactor <- diff(iRange.outcome)/diff(iRange.density)
                dataHist[dataHist$time1 == iTime,"density.norm"] <- iRange.outcome[1] + dataHist[dataHist$time1 == iTime,"density"]*iFactor
            }
    }

    ## ** prepare correlation
    if(n.time>2){
        dataGrid.upper <- dataGrid[dataGrid$position=="upper",]
        dataCor <- do.call(rbind, by(dataGrid.upper, droplevels(dataGrid.upper$time), function(iData){
            if(is.null(group)){
                iDF <- data.frame(outcome1 = as.numeric(NA),
                                  outcome2 = as.numeric(NA),
                                  time1 = unique(iData$time1),
                                  time2 = unique(iData$time2),
                                  cor = if(is.na(method.cor)){NA}else{stats::cor(iData$outcome1, iData$outcome2, method = method.cor, use = "pairwise")},
                                  n = NROW(is.na(iData)),
                                  n.NNA = sum(rowSums(is.na(iData))==0)
                                  )

                iDF$outcome2 <- mean(range(iData$outcome2,na.rm=TRUE))
                if(all(iData$time1 == level.time[1])){
                    if(density){
                        iDF$outcome1 <- mean(range(dataHist[dataHist$time1==level.time[1],"density.norm"]))
                    }else{                        
                        iDF$outcome1 <- mean(range(dataHist[dataHist$time1==level.time[1],"counts.norm"]))
                    }
                }else{
                    iDF$outcome1 <- mean(range(iData$outcome1, na.rm = TRUE))
                }
            }else{
                iCor <- by(iData, iData$group, function(iiData){
                    data.frame(group = unique(iiData$group),
                               cor = if(is.na(method.cor)){NA}else{stats::cor(iiData$outcome1, iiData$outcome2, method = method.cor, use = "pairwise")},
                               n = NROW(is.na(iiData)),
                               n.NNA = sum(rowSums(is.na(iiData))==0)
                               )
                })
                iDF <- data.frame(outcome1 = as.numeric(NA),
                                  outcome2 = as.numeric(NA),
                                  time1 = unique(iData$time1),
                                  time2 = unique(iData$time2),
                                  do.call(rbind,iCor))

                iDF$outcome2 <- mean(range(iData$outcome2,na.rm=TRUE))
                if(all(iData$time1 == level.time[1])){
                    if(density){
                        iRange <- range(dataHist[dataHist$time1==level.time[1],"density.norm"])
                    }else{                        
                        iRange <- range(dataHist[dataHist$time1==level.time[1],"counts.norm"])
                    }
                }else{
                    iRange <- range(iData$outcome1, na.rm = TRUE)
                }
                iDF$outcome1 <- seq(from = iRange[1], to = iRange[2], length.out = length(level.group)+2)[2:(1+length(level.group))]
            }
            return(iDF)
        }))
    }

    ## ** graphical display

    if(facet=="grid"){
        gg <- .ggscatterplot(dataHist = dataHist, dataGrid.lower = dataGrid.lower, dataCor = dataCor,
                             n.time = n.time, group = group, 
                             labeller = labeller, alpha.point = alpha.point, size.bar = size.bar,
                             density = density, alpha.area = alpha.area, method.cor = method.cor, size.text = size.text, digits = digits)
    }else if(facet == "grid2"){
        gg <- .ggscatterplot2(dataGrid.diag = dataGrid.diag, dataGrid.lower = dataGrid.lower, dataCor = dataCor,
                              breaks = breaks, n.time = n.time, group = group, 
                              labeller = labeller, alpha.point = alpha.point, position.bar = position.bar, size.bar = size.bar,
                              density = density, alpha.area = alpha.area, method.cor = method.cor, size.text = size.text, digits = digits)
    }

    ## ** export
    return(gg)
}

## * .ggscatterplot
.ggscatterplot <- function(dataHist, dataGrid.lower, dataCor,
                           n.time, group, 
                           labeller, alpha.point, size.bar, density, alpha.area, method.cor, size.text, digits){


    ## ** graphical display
    gg <- ggplot2::ggplot()
    if(n.time>2){
        gg <- gg +  ggplot2::facet_grid(time1~time2, scales = "free", labeller = labeller) + ggplot2::labs(x = "", y = "")
    }else if(n.time==2){
        gg <- gg + ggplot2::labs(x = dataGrid.lower$time2[1], y = dataGrid.lower$time1[1])
    }
    if(is.null(group)){
        gg <- gg + ggplot2::geom_point(data = dataGrid.lower,
                                       mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1),
                                       alpha = alpha.point)

        if(n.time>2){
            if(density){
                gg <- gg + ggplot2::geom_line(data = dataHist,
                                              mapping = ggplot2::aes(x = .data$mids, y = .data$density.norm),
                                              linewidth = 2)
                gg <- gg + ggplot2::geom_area(data = dataHist,
                                              mapping = ggplot2::aes(x = .data$mids, y = .data$density.norm),
                                              alpha = alpha.area)
            }else{
                gg <- gg + ggplot2::geom_segment(data = dataHist,
                                                 mapping = ggplot2::aes(x = .data$mids, xend = .data$mids, y = .data$base, yend = .data$counts.norm),
                                                 stat = "identity", linewidth = size.bar, alpha = alpha.area
                                                 )
            }

            if(is.na(method.cor)){
                gg <- gg + ggplot2::geom_text(data = dataCor,
                                              mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                     label = paste0(.data$n-.data$n.NNA," NA (",round(100*(1-.data$n.NNA/.data$n), digits[2]),"%)")),
                                              size = size.text, vjust = "inward", hjust = "inward")
            }else{
                gg <- gg + ggplot2::geom_text(data = dataCor,
                                              mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                     label = paste0("\u03C1=",round(.data$cor, digits[1]),"\n",round(100*(1-.data$n.NNA/.data$n), digits[2]),"% NA")),
                                              size = size.text, vjust = "inward", hjust = "inward")
            }
        }
    }else{
        gg <- gg + ggplot2::geom_point(data = dataGrid.lower,
                                       mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1, color = .data$group, shape = .data$group),
                                       alpha = alpha.point)
        if(n.time>2){
            if(density){
                gg <- gg + ggplot2::geom_line(data = dataHist,
                                              mapping = ggplot2::aes(x = .data$mids, y = .data$density.norm, group = .data$group, color = .data$group),
                                              linewidth = 2)
                gg <- gg + ggplot2::geom_area(data = dataHist,
                                              mapping = ggplot2::aes(x = .data$mids, y = .data$density.norm, fill = .data$group),
                                              alpha = alpha.area, position = "identity")
            }else{
                gg <- gg + ggplot2::geom_segment(data = dataHist,
                                                 mapping = ggplot2::aes(x = .data$mids, xend = .data$mids, y = .data$base, yend = .data$counts.norm, color = .data$group),
                                                 stat = "identity", linewidth = size.bar)
            }
            if(!is.na(method.cor)){
                gg <- gg + ggplot2::geom_text(data = dataCor,
                                              mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                     label = paste0("\u03C1=",round(.data$cor, digits[1])), color = .data$group),
                                              size = size.text, show.legend = FALSE, vjust = "inward", hjust = "inward")
            }else{
                gg <- gg + ggplot2::geom_text(data = dataCor,
                                              mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                     label = paste0(.data$n-.data$n.NNA," NA (",round(100*(1-.data$n.NNA/.data$n), digits[2]),"%)"), color = .data$group),
                                              size = size.text, show.legend = FALSE, vjust = "inward", hjust = "inward")
            }
            
        }
    }

    ## ** export
    return(gg)
}

## * .ggscatterplot2
.ggscatterplot2 <- function(dataGrid.diag, dataGrid.lower, dataCor,
                            breaks, n.time, group, 
                            labeller, alpha.point, position.bar, size.bar, density, alpha.area, method.cor, size.text, digits){

   ## ** graphical display
    gg <- ggplot2::ggplot()
    if(n.time>2){
        gg <- gg +  ggh4x::facet_grid2(time1~time2, scales = "free", labeller = labeller, independent = "y") + ggplot2::labs(x = "", y = "")
    }else if(n.time==2){
        gg <- gg + ggplot2::labs(x = dataGrid.lower$time2[1], y = dataGrid.lower$time1[1])
    }
    if(is.null(group)){
        gg <- gg + ggplot2::geom_point(data = dataGrid.lower,
                                       mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1),
                                       alpha = alpha.point)

        if(n.time>2){
            if(density){
                gg <- gg + ggplot2::geom_density(data = dataGrid.diag,
                                                 mapping = ggplot2::aes(x = .data$outcome1),
                                                 bw = breaks[1], kernel = breaks[2], alpha = alpha.area,
                                                 outline.type = "full", linewidth = size.bar, fill = "black")
            }else{
               gg <- gg + ggplot2::geom_histogram(data = dataGrid.diag,
                                                 mapping = ggplot2::aes(x = .data$outcome1),
                                                 bins = breaks, alpha = alpha.area, position = position.bar)
            }

            if(is.na(method.cor)){
                gg <- gg + ggplot2::geom_text(data = dataCor,
                                              mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                     label = paste0(.data$n-.data$n.NNA," NA (",round(100*(1-.data$n.NNA/.data$n), digits[2]),"%)")),
                                              size = size.text)
            }else{
                gg <- gg + ggplot2::geom_text(data = dataCor,
                                              mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                     label = paste0("\u03C1=",round(.data$cor, digits[1]),"\n",round(100*(1-.data$n.NNA/.data$n), digits[2]),"% NA")),
                                              size = size.text)
            }
        }
    }else{
        gg <- gg + ggplot2::geom_point(data = dataGrid.lower,
                                       mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1, color = .data$group, shape = .data$group),
                                       alpha = alpha.point)
        if(n.time>2){
            if(density){
                gg <- gg + ggplot2::geom_density(data = dataGrid.diag,
                                                 mapping = ggplot2::aes(x = .data$outcome1, fill = .data$group, color = .data$group),
                                                 bw = breaks[1], kernel = breaks[2], alpha = alpha.area,
                                                 outline.type = "full", linewidth = size.bar)
            }else{
                gg <- gg + ggplot2::geom_histogram(data = dataGrid.diag,
                                                   mapping = ggplot2::aes(x = .data$outcome1, fill = .data$group),
                                                   bins = breaks, alpha = alpha.area, position = position.bar)
            }

            print(dataCor)
            if(!is.na(method.cor)){
                gg <- gg + ggplot2::geom_text(data = dataCor,
                                              mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                     label = paste0("\u03C1=",round(.data$cor, digits[1])), color = .data$group),
                                              size = size.text, show.legend = FALSE, vjust = "inward")
            }else{
                gg <- gg + ggplot2::geom_text(data = dataCor,
                                              mapping = ggplot2::aes(x = .data$outcome2, y = .data$outcome1,
                                                                     label = paste0(.data$n-.data$n.NNA," NA (",round(100*(1-.data$n.NNA/.data$n), digits[2]),"%)"), color = .data$group),
                                              size = size.text, show.legend = FALSE, vjust = "inward")
            }
            
        }
    }

    ## ** export
    return(gg)
}

##----------------------------------------------------------------------
### scatterplot.R ends here
