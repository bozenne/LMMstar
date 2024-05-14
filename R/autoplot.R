### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  8 2021 (00:01) 
## Version: 
## Last-Updated: May 12 2024 (23:04) 
##           By: Brice Ozenne
##     Update #: 1306
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
##' @description Display fitted values or residual plot for the mean, variance, and correlation structure.
##' Can also display quantile-quantile plot relative to the normal distribution.
##'
##' @param object,x a \code{lmm} object.
##' @param type [character] the type of plot \itemize{
##' \item \code{"fit"}: fitted values over repetitions.
##' \item \code{"qqplot"}: quantile quantile plot of the normalized residuals
##' \item \code{"correlation"}: residual correlation over repetitions
##' \item \code{"scatterplot"}: normalized residuals vs. fitted values (diagnostic for missing non-linear effects),
##' \item \code{"scatterplot2"}: square root of the normalized residuals vs. fitted values (diagnostic for heteroschedasticity),
##' \item \code{"partial"}: partial residual plot.
##' }
##' @param type.residual [character] the type of residual to be used. Not relevant for \code{type="fit"}.
##' By default, normalized residuals are used except when requesting a partial residual plot
##' where this argument specify the variable relative to which the partial residuals are computed (argument \code{variable} in \code{\link{residuals.lmm}}).
##' @param at [data.frame] values for the covariates at which to evaluate the fitted values or partial residuals.
##' @param time.var [character] x-axis variable for the plot.
##' @param obs.alpha [numeric, 0-1] When not NA, transparency parameter used to display the original data by cluster.
##' @param obs.size [numeric vector of length 2] size of the point and line for the original data.
##' @param facet [formula] split the plot into a matrix of panels defined by the variables in the formula.
##' Internally it calls \code{ggplot2::facet_wrap} or \code{ggplot2::facet_grid} depending on whether the formula contains a variable on the left hand side.
##' @param scales,labeller [character] Passed to \code{ggplot2::facet_wrap}.
##' @param facet_nrow [integer] number of rows of panels in the graphical display.
##' @param facet_ncol [integer] number of columns of panels  in the graphical display.
##' @param color [character] name of the variable in the dataset used to color the curve. No color is used when set to \code{FALSE}.
##' @param position [character] relative position of the points when colored according to a variable.
##' @param ci [logical] should confidence intervals be displayed?
##' @param ci.alpha [numeric, 0-1] When not NA, transparency parameter used to display the confidence intervals.
##' @param size.text [numeric, >0] size of the font used to display text.
##' @param position.errorbar [character] relative position of the errorbars.
##' @param mean.size [numeric vector of length 2] size of the point and line for the mean trajectory.
##' @param ylim [numeric vector of length 2] the lower and higher value of the vertical axis.
##' @param ... arguments passed to the \code{predict.lmm} or \code{autoplot.residual_lmm} functions.
##'
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to create the graphical display.
##' \item \code{plot}: ggplot object.
##' }
##'
##' @seealso
##' \code{\link{plot.lmm}} for other graphical display (residual plots, partial residual plots).
##'
##' @keywords hplot
##' 
##' @examples
##' if(require(ggplot2)){
##' 
##' #### simulate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' dL$X1 <- as.factor(dL$X1)
##' 
##' #### fit Linear Mixed Model ####
##' eCS.lmm <- lmm(Y ~ visit + X1 + X6,
##'                repetition = ~visit|id, structure = "CS", data = dL, df = FALSE)
##'
##' #### model fit ####
##' plot(eCS.lmm, type = "fit", facet =~X1)
##' ## customize display
##' gg <- autoplot(eCS.lmm, type = "fit", facet =~X1)$plot
##' gg + coord_cartesian(ylim = c(0,6))
##' ## restrict to specific covariate value
##' plot(eCS.lmm, type = "fit", at = data.frame(X6=1), color = "X1")
##'
##' #### qqplot ####
##' plot(eCS.lmm, type = "qqplot")
##' plot(eCS.lmm, type = "qqplot", engine.qqplot = "qqtest")
##'
##' #### residual correlation ####
##' plot(eCS.lmm, type = "correlation")
##'
##' #### residual trend ####
##' plot(eCS.lmm, type = "scatterplot")
##' 
##' #### residual heteroschedasticity ####
##' plot(eCS.lmm, type = "scatterplot2")
##'
##' #### partial residuals ####
##' plot(eCS.lmm, type = "partial", type.residual = "visit") 
##' plot(eCS.lmm, type = "partial", type.residual = c("(Intercept)","X1","visit"))
##' plot(eCS.lmm, type = "partial", type.residual = c("(Intercept)","X1","visit"),
##' facet = ~X1)
##' }

## * autoplot.lmm (code)
##' @export
autoplot.lmm <- function(object, type = "fit", type.residual = NULL, 
                         obs.alpha = 0, obs.size = NULL, facet = NULL, facet_nrow = NULL, facet_ncol = NULL, scales = "fixed", labeller = "label_value", 
                         at = NULL, time.var = NULL, color = NULL, position = NULL, ci = TRUE, ci.alpha = NULL, 
                         ylim = NULL, mean.size = c(3, 1), size.text = 16, position.errorbar = "identity", ...){

    ## use [] to keep attribute reference for partial residuals
    type[] <- match.arg(type, c("fit",
                                "partial","partial-center",
                                "qqplot","correlation","scatterplot","scatterplot2")) 

    if(type=="fit"){ ## model fit
        if(is.null(ci.alpha)){ci.alpha <- 0.25}
        if(is.null(obs.size)){obs.size <- c(2,0.5)}
        if(is.null(color)){color <- TRUE}
        if(is.null(position)){position <- "identity"}
        out <- .autofit(object,
                        facet = facet, facet_nrow = facet_nrow, facet_ncol = facet_ncol, scales = scales, labeller = labeller,
                        obs.alpha = obs.alpha,
                        obs.size = obs.size,
                        at = at,
                        time.var = time.var,
                        color = color,                        
                        position = position,                        
                        ci = ci,
                        ci.alpha = ci.alpha,
                        ylim = ylim,
                        mean.size = mean.size,
                        size.text = size.text,
                        position.errorbar = position.errorbar,
                        ...)
    }else if(type %in% c("partial","partial-center")){ ## partial residual plot

        test <- c(obs.alpha = !is.null(obs.alpha) & abs(obs.alpha)>0)
        if(any(test)){
            message("Arugment(s) \'",paste(names(test[test]), collapse = "\', \'"),"\' diregarded  when displaying partial residuals\n")
        }

        dots <- list(...)
        if(is.null(type.residual)){
            ## handle the case where the user is specifying the argument var (from residuals.lmm) instead of type.residuals
            if(!is.null(dots$var)){
                type.residual <- dots$variable
                dots$variable <- NULL
            }else{
                type.residual <- attr(object$design$mean,"variable")[1]
            }            
        }

        ## extract partial residuals
        outRes <- stats::residuals(object, type = type, format = c("wide","long"), variable = type.residual, at = at, keep.data = TRUE, simplify = FALSE,
                                   fitted.ci = !is.null(ci.alpha) && !is.na(ci.alpha))

        out <- do.call(autoplot.residuals_lmm, args = c(list(outRes,
                                                             type = type,
                                                             time.var = time.var,
                                                             size.text = size.text,
                                                             facet = facet, facet_nrow = facet_nrow, facet_ncol = facet_ncol, 
                                                             scales = scales,
                                                             labeller = labeller,
                                                             color = color,
                                                             position = position,                        
                                                             obs.size = obs.size,
                                                             ci.alpha = ci.alpha),
                                                        dots))

    }else{ ## residual plot
        test <- c(obs.alpha = !is.null(obs.alpha) & abs(obs.alpha)>0,
                  at = !is.null(at))
        if(any(test)){
            message("Arugment(s) \'",paste(names(test[test]), collapse = "\', \'"),"\' diregarded  when displaying residuals\n")
        }
        if(is.null(type.residual)){
            type.residual <- "normalized"
        }
        outRes <- residuals(object, type = type.residual, format = c("wide","long"), keep.data = TRUE, simplify = FALSE)
        out <- autoplot.residuals_lmm(outRes,
                                      type = type,
                                      time.var = time.var,
                                      size.text = size.text,
                                      mean.size = mean.size,
                                      ci.alpha = ci.alpha,
                                      facet = facet, facet_nrow = facet_nrow, facet_ncol = facet_ncol, 
                                      scales = scales,
                                      labeller = labeller,
                                      color = color,
                                      position = position,                        
                                      obs.size = obs.size,
                                      ...)
    }

    ## ** export
    return(out)

}

## ** .autofit (helper to autoplot.lmm)
.autofit <- function(object, facet, facet_nrow, facet_ncol, scales, labeller,
                     obs.alpha, obs.size,
                     at, time.var, color, position, ci, ci.alpha, 
                     ylim, mean.size, size.text, position.errorbar, ...){

    if(is.null(time.var) && object$time$n==1){
        stop("Cannot display the fitted values over time when there only is a single timepoint. \n")
    }

    ## ** extract from object
    object.data <- object$data
        
    Upattern <- object$design$vcov$Upattern
    pattern.cluster <- object$design$vcov$pattern

    outcome.var <- object$outcome$var

    if(is.null(time.var)){
        time.var.plot <- "XXtimeXX" ## nice as it sure to be a categorical variable
        if(!is.null(attr(object$time$var,"original")) && all(!is.na(attr(object$time$var,"original")))){
            xlabel.plot <- attr(object$time$var,"original")
        }else{
            xlabel.plot <- ""
        }
        
    }else{

        if(length(time.var)>1){
            if("sep" %in% names(time.var)){
                sep.time.var <- time.var["sep"]
                time.var <- time.var[names(time.var) != "sep"]
            }else{
                sep.time.var <- ", "
            }
            if(all(time.var %in% names(object$data.original))){
                object.data$XXtimeXX <- nlme::collapse(object$data.original[time.var], sep = sep.time.var)[object.data$XXindexXX]
                time.var <- "XXtimeXX"
            }else{
                stop("Incorrect value for argument \'time.var\'. \n",
                     "No column ",time.var," found in the dataset used to fit the lmm. \n")
            }
        }else if(length(time.var)==1 && time.var %in% names(data) == FALSE){
            if(time.var %in% names(object$data.original)){
                object.data[[time.var]] <- object$data.original[[time.var]][object.data$XXindexXX]
            }else{
                stop("Incorrect value for argument \'time.var\'. \n",
                     "No column ",time.var," found in the dataset used to fit the lmm. \n")
            }
        }
        time.var.plot <- time.var
        xlabel.plot <- time.var
    }

    time.var <- attr(object$time$var,"original") ## need to be after statement on time.var.plot to avoid confusion
    mu.var <- formula2var(object$formula$mean)$var$regressor
    if(length(time.var) == 0 && length(mu.var) == 0){
        message("There is nothing to be displayed: empty time variable and no covariate for the mean structure. \n")
        return(NULL)
    }
    
    if(length(color) == 0){
        color <- NULL
    }else if(length(color)>1){
        stop("Argument \'color\' should either be NULL \n",
             "        or have length 1 and be TRUE or a variable name in the data used to fit the lmm. \n")
    }else if(length(color) == 1 && is.character(color)){
        if(color %in% names(object.data) == FALSE){
            if(color %in% names(object$data.original)){
                object.data[[color]] <- object$data.original[[color]][object.data$XXindexXX]
            }else{
                stop("Incorrect value for argument \'color\'. \n",
                     "No column ",color," found in the dataset used to fit the lmm. \n")
            }
        }
    }else if(!identical(color,FALSE) && !identical(color,TRUE)){
        stop("Argument \'color\' should either be NULL \n",
             "        or have length 1 and be TRUE or a variable name in the data used to fit the lmm. \n")
    }
    
    if(!is.null(facet)){
        if(!inherits(facet,"formula")){
            stop("Argument \'facet\' must either be NULL or a formula. \n",
                 "It will be passed to ggplot2::facet_wrap or ggplot2::facet_grid,\n",
                 " depending on whether there are variables on the left hand side of the formula. \n")
        } 
        if(any(all.vars(facet) %in% names(object.data) == FALSE)){
            hide.name <- c(paste0("XX",c("index","cluster","time","strata"),"XX"),
                           paste0("XX",c("index","cluster","time","strata"),".indexXX"),
                           time.var, color, outcome.var, object$cluster$var)
            stop("When a formula, argument \'facet\' should contain variables available in the dataset used to fit the model. \n",
                 "Unknown variable(s): \"",paste(all.vars(facet)[all.vars(facet) %in% names(object.data) == FALSE], collapse = "\" \""),"\" \n",
                 "Available variable(s): \"",paste(setdiff(names(object.data),c(all.vars(facet),hide.name)), collapse = "\" \""),"\" \n")
        }
    }


    ## ** find representative individuals
    order.nrep <- names(sort(stats::setNames(Upattern$n.time, Upattern$name), decreasing = TRUE))
    col.pattern <- factor(pattern.cluster, order.nrep)[object.data[["XXclusterXX"]]]

    ## put observations with full data first to avoid "holes"  in the plot
    reorder <- order(col.pattern,object.data[["XXclusterXX"]],object.data[["XXtimeXX"]])
    data <- object.data[reorder,,drop=FALSE]
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

    ## design matrix: find unique combinations of covariates
    timemu.var <- stats::na.omit(union(time.var, mu.var))
    X.beta <- stats::model.matrix(object, effects = "mean",
                                  newdata = data[,timemu.var, drop=FALSE])
    IX.beta <- nlme::collapse(X.beta, as.factor = TRUE)
    vec.X.beta <- tapply(IX.beta, data[["XXclusterXX"]],paste, collapse = "_XXX_")
    UX.beta <- unique(vec.X.beta)

    ## remove duplicates due to missing values (unequal number of repetitions)
    UX.ntime <- table(droplevels(data[["XXclusterXX"]][data[["XXclusterXX"]] %in% names(UX.beta)]))
    if(length(UX.beta)>1 && length(unique(UX.ntime))>1){
        test.UX.beta <- rep(TRUE, length(UX.beta))
        for(iUX in 1:length(UX.beta)){ ## iUX <- 2

            iX <- IX.beta[data[["XXclusterXX"]] == names(UX.beta)[iUX]]
            iTest <- sapply(names(UX.beta)[-iUX], function(iId){all(iX %in% IX.beta[data[["XXclusterXX"]] == iId])})
            if(any(iTest)){
                test.UX.beta[iUX] <- FALSE
            }

        }
        UX.beta <- UX.beta[test.UX.beta]
    }

    lsID.beta <- lapply(UX.beta, function(iX.beta){names(iX.beta == vec.X.beta)}) ## cluster(s) within each mean pattern
    newdata <- data[data[["XXclusterXX"]] %in% names(UX.beta),]

    if(identical(color,TRUE)){
        mean.var <- all.vars(stats::delete.response(stats::terms(stats::formula(object, effects = "mean"))))
        if(length(mean.var)>0){
            newdataRed <- newdata[order(newdata[["XXclusterXX"]]),mean.var,drop=FALSE]
            order.cluster <- droplevels(newdata[["XXclusterXX"]][order(newdata[["XXclusterXX"]])])

            M.duplicated <- apply(newdataRed, 2, function(iCol){unlist(tapply(iCol, order.cluster, function(iColCluster){duplicated(iColCluster)[-1]}))})
            if(length(M.duplicated)==0){
                color <- NULL
            }else{
                color <- setdiff(names(which(colSums(M.duplicated)==NROW(M.duplicated))), all.vars(facet))
            }
        }else{
            color <- NULL
        }
        if(length(color)>1){
            if(paste(color,collapse=".") %in% names(newdataRed)){
                stop("Cannot use argument \'color\'=TRUE when the dataset contain a column ",paste(color,collapse="."),". \n",
                     "This name is used internally. \n")
            }
            newdata[[paste(color,collapse=".")]] <- nlme::collapse(newdata[,color,drop=FALSE], as.factor = TRUE)
            color <- paste(color,collapse=".")
        }else if(length(color)==0){
            color <-  NULL
        }
    }else if(identical(color,FALSE)){
        color <- NULL
    }

    if(!is.na(obs.alpha) && obs.alpha>0 && length(color)>1 && color %in% names(data) == FALSE){
        ls.UX <- lapply(as.character(unique(newdata[["XXclusterXX"]])), function(iC){
            iVec <- nlme::collapse(data[data[["XXclusterXX"]] %in% iC,mu.var,drop=FALSE], as.factor = FALSE)
            cbind(repetition = data[data[["XXclusterXX"]] %in% iC,"XXtimeXX",drop=FALSE], lp = iVec)
        })
    
        index.X <- unlist(lapply(as.character(unique(data[["XXclusterXX"]])), function(iC){
            iVec <- nlme::collapse(data[data[["XXclusterXX"]] %in% iC,mu.var,drop=FALSE], as.factor = FALSE)
            iM <- cbind(repetition = data[data[["XXclusterXX"]] %in% iC,"XXtimeXX",drop=FALSE], lp = iVec)
            iScore <- unlist(lapply(ls.UX, function(iUX){sum(iUX[match(iM[,"Days"],iUX[,"Days"]),"lp"]==iM[,"lp"])}))
            which.max(iScore)
        }))
        data[[color]] <- sort(unique(newdata[[color]]))[index.X]
    }

    ## ** compute fitted curve
    if(!is.na(obs.alpha) && obs.alpha>0){
        preddata <- cbind(data, stats::predict(object, newdata = data[,timemu.var, drop=FALSE], simplify = FALSE, ...))
    }else{
        preddata <- cbind(newdata, stats::predict(object, newdata = newdata[,timemu.var, drop=FALSE], simplify = FALSE, ...))
    }
    if("lower" %in% names(preddata) == FALSE){
        preddata$lower <- NA
    }
    if("upper" %in% names(preddata) == FALSE){
        preddata$upper <- NA
    }

    ## ** add missing times (if any)
    if(is.factor(preddata[[time.var.plot]])){
        U.time <- levels(preddata[[time.var.plot]])
    }else{
        U.time <- sort(unique(preddata[[time.var.plot]]))
    }
    keep.col <- c("XXclusterXX",time.var.plot,outcome.var,color,all.vars(facet),"estimate","lower","upper")
    preddata <- do.call(rbind,by(preddata[keep.col], preddata$XXclusterXX, function(iDF){ ## iDF <- preddata[preddata$XXclusterXX==1,]
        if(all(U.time %in% iDF[[time.var.plot]])){
            return(iDF)
        }else{
            iNewTime <- setdiff(U.time,iDF[[time.var.plot]])

            iDFextra <- as.data.frame(lapply(names(iDF), function(iVar){
                if(iVar == time.var.plot){
                    return(iNewTime)
                }else if(iVar %in% c("XXclusterXX",color,all.vars(facet))){
                    return(rep(iDF[[iVar]][1], length(iNewTime)))
                }else{
                    return(rep(NA, length(iNewTime)))
                } 
            }))
            names(iDFextra) <- names(iDF)
            return(rbind(iDF,iDFextra))
        }
    }))

    ## ** generate plot
    gg <- ggplot2::ggplot(data = preddata, mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                  y = .data$estimate,
                                                                  group = .data$XXclusterXX))
    test.line <- all(tapply(preddata[["XXclusterXX"]],preddata[[time.var.plot]], function(iX){any(duplicated(iX))})==FALSE)

    if(!is.na(obs.alpha) && obs.alpha>0){
        if(!is.null(color)){
            gg <- gg + ggplot2::geom_point(mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                  y = .data[[outcome.var]],
                                                                  group = .data$XXclusterXX,
                                                                  color = .data[[color]]),
                                           alpha = obs.alpha,
                                           position = position,
                                           size = obs.size[1])
            
            if(test.line){
                gg <- gg + ggplot2::geom_line(mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                     y = .data[[outcome.var]],
                                                                     group = .data$XXclusterXX,
                                                                     color = .data[[color]]),
                                              alpha = obs.alpha,
                                              position = position,
                                              linewidth = obs.size[2])
            }
            ## gg + facet_wrap(~XXclusterXX)
        }else{
            gg <- gg + ggplot2::geom_point(mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                  y = .data[[outcome.var]],
                                                                  group = .data$XXclusterXX),
                                           alpha = obs.alpha,
                                           size = obs.size[1])
            if(test.line){
                gg <- gg + ggplot2::geom_line(mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                     y = .data[[outcome.var]],
                                                                     group = .data$XXclusterXX),
                                              alpha = obs.alpha,
                                              linewidth = obs.size[2])
            }
        }
    }
    if(ci){
        if(is.na(ci.alpha)){
            gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), position = position.errorbar)
        }else{
            if(!is.null(color)){
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data[[color]]), alpha = ci.alpha, position = position)
            }else{
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = ci.alpha)
            }
        }
    }
    if(!is.null(color)){
        gg <- gg + ggplot2::geom_point(ggplot2::aes(color = .data[[color]]), size = mean.size[1], position = position)
        if(test.line){
            ## exclude NA' to avoid multiple apparent fit (1-2-3-4, 1-4, ...)
            gg <- gg + ggplot2::geom_line(data = preddata[preddata[["XXclusterXX"]] %in% names(UX.beta),], ggplot2::aes(color = .data[[color]]), linewidth = mean.size[2], position = position)
        }
    }else{
        gg <- gg + ggplot2::geom_point(size = mean.size[1])
        if(test.line){
            ## exclude NA' to avoid multiple apparent fit (1-2-3-4, 1-4, ...)
            gg <- gg + ggplot2::geom_line(data = preddata[preddata[["XXclusterXX"]] %in% names(UX.beta),], linewidth = mean.size[2])
        }
    }
    if(!is.null(facet) & length(all.vars(facet))>0){
        if(attr(stats::terms(facet),"response")==0){
            gg  <- gg + ggplot2::facet_wrap(facet, nrow = facet_nrow, ncol = facet_ncol, scales = scales, labeller = labeller)
        }else{
            gg  <- gg + ggplot2::facet_grid(facet, nrow = facet_nrow, ncol = facet_ncol, scales = scales, labeller = labeller)
        }
    }
    gg  <- gg + ggplot2::ylab(outcome.var) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
    if(!is.null(time.var.plot) && any(!is.na(time.var.plot))){
        if(length(xlabel.plot)==1){
            if(xlabel.plot %in% names(object$data.original)){ ## single observed variable 
                gg  <- gg + ggplot2::xlab(xlabel.plot)
            }else if(time.var.plot=="XXtimeXX"){ ## internally made time variable
                gg  <- gg + ggplot2::xlab(NULL)
            }
        }else { ## multiple (observed) variable aggregated into a single name
            gg  <- gg + ggplot2::xlab(paste(stats::na.omit(xlabel.plot), collapse = ", "))
        }
    }
    if(!is.null(ylim)){
        gg <- gg + ggplot2::coord_cartesian(ylim = ylim)
    }

    ## ** export
    return(list(data = preddata, plot = gg))    
}


## * autoplot.partialCor (documentation)
##' @title Graphical Display For Partial Correlation
##' @description Extract and display the correlation modeled via the linear mixed model.
##'
##' @param object,x a \code{partialCor} object.
##' @param size.text [numeric, >0] size of the font used to display text.
##' @param limits [numeric vector of length 2] minimum and maximum value of the colorscale relative to the correlation.
##' @param low,mid,high [character] color for the the colorscale relative to the correlation.
##' @param midpoint [numeric] correlation value associated with the color defined by argument \code{mid}.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to create the graphical display.
##' \item \code{plot}: ggplot object.
##' }
##' 
##' @keywords hplot
##' 
##' @examples
##' if(require(ggplot2)){
##' data(gastricbypassL, package = "LMMstar")
##' 
##' e.pCor <- partialCor(c(weight,glucagonAUC)~time, repetition = ~visit|id,
##'                      data = gastricbypassL)
##' plot(e.pCor)
##' }
##' 

## * autoplot.partialCor (code)
##' @export
autoplot.partialCor <- function(object, size.text = 16,
                                limits = c(-1,1.00001), low = "blue", mid = "white", high = "red", midpoint = 0, ...){

    object.lmm <- attr(object,"lmm")
    Sigma_t <- sigma(object.lmm)
    name.time <- object.lmm$time$levels
    if(!is.matrix(Sigma_t)){
        stop("Could not extract a unique covariance matrix. \n")
    }
    Sigma_t <- stats::cov2cor(Sigma_t)
        
    ## from matrix to long format
    table <- as.data.frame(cbind(which(is.na(NA*Sigma_t), arr.ind = TRUE),value = as.numeric(Sigma_t)))
    rownames(table) <- NULL
    table$col <- factor(colnames(Sigma_t)[table$col], levels = name.time)
    table$row <- factor(rownames(Sigma_t)[table$row], levels = name.time)
    
    gg <- ggplot2::ggplot(table) + ggplot2::geom_tile(ggplot2::aes(x = .data$row, y = .data$col, fill = .data$value))

    if(!is.null(mid)){
        gg <- gg + ggplot2::scale_fill_gradient2(limits = limits, midpoint = midpoint, low = low, mid = mid, high = high)
    }else{
        gg <- gg + ggplot2::scale_fill_gradient(limits = limits, low = low, high = high)
    }
    gg <- gg + ggplot2::labs(x = NULL, y = NULL, fill = "correlation") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
        
    gg <- gg + ggplot2::theme(text = ggplot2::element_text(size=size.text))

    ## ** export
    return(list(data = table,
                plot = gg))
}



## * autoplot.profile_lmm (documentation)
##' @title Graphical Display of Profile Likelihood
##' @description Graphical representation of the profile likelihood from a linear mixed model
##' 
##' @param object,x an object of class \code{profile_lmm}, output of the \code{profile.lmm} function.
##' @param type [character] Should the log-likelihood (\code{"logLik"}) or the ratio to the maximum likelihood (\code{"ratio"}) be displayed?
##' @param quadratic [logical] Should a quadratic approximation of the likelihood be displayed?
##' @param ci [logical] Should a 95\% confidence intervals obtained from the Wald test (vertical lines) and Likelihood ratio test (horizontal line) be displayed?
##' @param size [numeric vector of length 4] Size of the point for the MLE,
##' width of the line representing the likelihood,
##' width of the corresponding quadratic approximation,
##' and width of the line representing the confidence intervals.
##' @param linetype [integer vector of length 2] type of line used to represent the quadratic approximation of the likelihood
##' and the confidence intervals.
##' @param shape [integer, >0] type of point used to represent the MLE.
##' @param scales,nrow,ncol argument passed to \code{ggplot2::facet_wrap}.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A list with three elements \itemize{
##' \item \code{data.fit}: data containing the quadratice approximation of the log-likelihood
##' \item \code{data.ci}: data containing the confidence intervals.
##' \item \code{plot}: ggplot object.
##' }
##' @keywords hplot

## * autoplot.profile_lmm (code)
##' @export
autoplot.profile_lmm <- function(object, type = "logLik", quadratic = TRUE, ci = FALSE,
                                 size = c(3,2,1,1), linetype = c("dashed","dashed","dashed"), shape = 19, scales = "free", nrow = NULL, ncol = NULL,
                                 ...){

    ## ** normalize arguments
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    type <- match.arg(type, c("logLik","ratio"))
    args <- attr(object,"args")
    profile.likelihood <- args$profile.likelihood
    conf.level <- args$conf.level
    object.logLik <- args$logLik
    effects <- args$effects
    transform.names <- args$transform.names
    name.p <- args$name.p
    p.trans2 <- args$ci

    if(ci && profile.likelihood==FALSE){
        stop("Can only display the confidence intervals when performing profile likelihood. \n")
    }

    ## ** build display
    if(type == "ratio"){
        name.y <- "likelihood.ratio"
        reference <- 1
        legend.y <- "Likelihood relative to the maximum likelihood"
        fff <- likelihood.ratio ~ 0+I(value.trans-mean(value.trans)) + I((value.trans-mean(value.trans))^2)
    }else if(type == "logLik"){
        name.y <- "logLik"
        reference <- object.logLik
        legend.y <- "Log-likelihood"
        fff <- logLik ~ 0+I(value.trans-mean(value.trans)) + I((value.trans-mean(value.trans))^2)
    }

    gg <- ggplot2::ggplot(object, ggplot2::aes(x = .data$value.trans, y = .data[[name.y]]))
    gg <- gg + ggplot2::ylab(legend.y)
    if(profile.likelihood>0){
        if(ci){
            gg <- gg + ggplot2::ggtitle(paste("Profile maximum likelihood estimation for parameter (95% CI):"))
        }else{
            gg <- gg + ggplot2::ggtitle(paste("Profile maximum likelihood estimation for parameter:"))
        }
    }else{
        gg <- gg + ggplot2::ggtitle(expression(paste("Varying a ",bold('single')," parameter:")))
    }
    gg <- gg + ggplot2::facet_wrap(~param, scales= scales, nrow = nrow, ncol = ncol)
    if(size[2]>0){
        gg <- gg + ggplot2::geom_line(linewidth = size[2])
    }
    if(size[3]>0 && quadratic){
        df.fit <- do.call(rbind,by(object, object$param, function(iDF){ ## iDF <- object[object$param=="sigma",]
            iDF$myset <- reference
            ## FOR CRAN test
            myset <- NULL
            iLM <- stats::lm(fff, data = iDF, offset = myset)
            iDF[[name.y]] <- stats::predict(iLM, newdata = iDF)
            return(iDF)
        }))
        df.fit$param <- factor(df.fit$param, levels = levels(object$param))
        gg <- gg + ggplot2::geom_line(data = df.fit, linewidth = size[3], linetype = linetype[1], ggplot2::aes(color = "quadratic approximation"))
    }else{
        df.fit <- NULL
    }
    if(size[1]>0){
        gg <- gg  + ggplot2::geom_point(data = object[object$optimum==TRUE,,drop=FALSE], ggplot2::aes(color = "MLE"), size = size[1], shape = shape)
    }
    if(ci){
        if(ci<=1){
            gg <- gg  + ggplot2::geom_abline(slope = 0, intercept = object.logLik - stats::qchisq(0.95, df = 1)/2, size = size[4], linetype = linetype[2])
        }else{
            gg <- gg  + ggplot2::geom_abline(slope = 0, intercept = exp(- stats::qchisq(0.95, df = 1)/2), size = size[4], linetype = linetype[2])
        }

        df.ci <- cbind(param = rownames(p.trans2), p.trans2)[which(name.p %in% effects),,drop=FALSE]
        if(transform.names){
            df.ci$param <- factor(df.ci$param, levels = levels(object$param))
            ## df.ci$param <- factor(name.p.trans[match(df.ci$param,name.p)], levels = levels(df.profile$param))
        }else{
            df.ci$param <- factor(df.ci$param, levels = levels(object$param))
        }
        gg <- gg + ggplot2::geom_vline(data = df.ci, mapping = ggplot2::aes(xintercept = .data$lower), size = size[4], linetype = linetype[2])
        gg <- gg + ggplot2::geom_vline(data = df.ci, mapping = ggplot2::aes(xintercept = .data$upper), size = size[4], linetype = linetype[2])
    }else{
        df.ci <- NULL
    }
    gg <- gg + ggplot2::xlab("") + ggplot2::labs(color = "") + ggplot2::theme(legend.position = "bottom")

    ## ** export
    return(list(data.fit = df.fit,
                data.ci = df.ci,
                plot = gg))
}


## * autoplot.residuals_lmm (documentation)
##' @title Graphical Display of the Residuals
##' @description Graphical representation of the residuals from a linear mixed model.
##' Require a long format (except for the correlation where both format are accepted) and having exported the dataset along with the residual (argument \code{keep.data} when calling \code{residuals.lmm}).
##' 
##' @param object,x an object of class \code{residuals_lmm}, output of the \code{residuals.lmm} function.
##' @param type [character] Should a qqplot (\code{"qqplot"}), or a heatmap of the correlation between residuals  (\code{"correlation"}, require wide format), or a plot of residuals along the fitted values (\code{"scatterplot"}, require long format) be displayed?
##' @param type.residual [character] Type of residual for which the graphical representation should be made.
##' @param time.var [character] x-axis variable for the plot. Only relevant when argument type is one of \code{"scatterplot"}, \code{"scatterplot2"}, \code{"partial"}, \code{"partial-center"},
##' @param engine.qqplot [character] Should ggplot2 or qqtest be used to display quantile-quantile plots?
##' Only used when argument \code{type} is \code{"qqplot"}.
##' @param add.smooth [logical] should a local smoother be used to display the mean of the residual values across the fitted values.
##' Only relevant for when argument \code{type} is \code{"scatterplot"}.
##' @param digits.cor [integer, >0] Number of digit used to display the correlation coefficients?
##' No correlation coefficient is displayed when set to 0. Only used when argument \code{plot} is \code{"correlation"}.
##' @param size.text [numeric, >0] Size of the font used to displayed text when using ggplot2.
##' @param facet [formula] split the plot into a matrix of panels defined by the variables in the formula.
##' Internally it calls \code{ggplot2::facet_wrap} or \code{ggplot2::facet_grid} depending on whether the formula contains a variable on the left hand side.
##' @param facet_nrow [integer] number of rows of panels in the graphical display.
##' @param facet_ncol [integer] number of columns of panels  in the graphical display.
##' @param scales,labeller [character] Passed to \code{ggplot2::facet_wrap}.
##' @param color [character] color of the dots representing the observations.
##' When displaying partial residuals, should contain a second color indicating how to display the model fit.  
##' @param position [character] relative position of the points when colored according to a variable.
##' @param obs.size [numeric vector] size of the dots representing the observations.
##' @param mean.size [numeric vector of length 2] size of the point and line for the mean trajectory.
##' @param ci.alpha [numeric, 0-1] When not NA, transparency parameter used to display the confidence intervals.
##' @param ... Not used. For compatibility with the generic method.
##'  
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to generate the plot.
##' \item \code{plot}: ggplot object.
##' }
##' 
##' @keywords hplot

## * autoplot.residuals_lmm (code)
##' @export
autoplot.residuals_lmm <- function(object, type = NULL, type.residual = NULL, time.var = NULL, facet = NULL, facet_nrow = NULL, facet_ncol = NULL,
                                   engine.qqplot = "ggplot2", add.smooth = TRUE, digits.cor = 2, size.text = 16,
                                   color = NULL, obs.size = NULL, mean.size = c(3, 1), ci.alpha = 0.25, 
                                   position = NULL, scales = "fixed", labeller = "label_value", ...){

    ## ** check arguments
    call <- match.call()

    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## args
    args <- attr(object,"args")
    if(is.null(args)){
        stop("The argument \'simplify\' must be to FALSE when calling residuals() to obtain a graphical display. \n")
    }
    if(args$keep.data == FALSE && type != "correlation"){
        stop("The argument \'keep.data\' must be to TRUE when calling residuals() to obtain a graphical display. \n")
    }
    args.type <- args$type
    n.type <- length(args.type)
    index.time <- attr(object,"index.time") ## save in case no missing time variable in the long format

    ## type of residual
    if(is.null(type.residual) && is.null(type) && "partial" %in% args.type){
        type <- "partial"
        type.residual <- "partial"
    }

    if(is.null(type.residual)){
        if(type == "partial"){
            type.residual <- "partial"
        }else if(n.type==1){
            type.residual <- args.type
        }else{
            stop("Different types of residuals are available: \"",paste(args.type, collapse = "\", \""),"\". \n",
                 "Select one type via the argument \'type.residual\'. \n")
        }
    }else if(n.type==1){
        if(type == "partial" && !identical(type.residual,"partial")){
            message("Argument \'type.residual\' ignored when displaying partial residuals. \n")
            type.residual <- "partial"
        }else if(type.residual %in% args$type == FALSE){
            ## when partial, type.residual encodes covariates names not types of residuals
            stop("Requested type of residual not available. \n",
                 "Available type(s) of residuals: \"",paste(args.type, collapse = "\", \""),"\"\n")
        }
    }else{
        stop("Can only display one type of residual. \n")
    }

    ## type of plot
    if(is.null(type)){
        type <- "qqplot"
    }else{
        type <- match.arg(type, c("qqplot","correlation","scatterplot","scatterplot2","partial","partial-center"))
    }
    if(length(add.smooth)==1){
        add.smooth <- rep(add.smooth,2)
    }

    ## time.var
    if(!is.null(time.var) && type %in% c("qqplot","correlation")){
        message("Argument \'time.var\' is ignored when type is ",type,". \n")
    }
    if(length(time.var)>1){
        stop("Argument \'time.var\' should be NULL or have length 1. \n")
    }
    if(length(time.var) == 1 && time.var %in% names(object) == FALSE){
        stop("Argument \'time.var\' should be refer to an available variable in the output of residual.lmm. \n",
             "Available variables: \"",paste(names(object), collapse="\", \""),"\"\n")
    }

    ## number of timepoints
    n.time <- args$n.time
    name.time <- args$name.time

    ## facet
    if(!is.null(facet)){
        if(!inherits(facet,"formula")){
            stop("Argument \'facet\' must either be NULL or a formula. \n",
                 "It will be passed to ggplot2::facet_wrap or ggplot2::facet_grid,\n",
                 " depending on whether there are variables on the left hand side of the formula. \n")
        }
        if(any(all.vars(facet) %in% names(object) == FALSE) && any(all.vars(facet) %in% names(attr(object,"wide")) == FALSE)){
            hide.name <- c(paste0("XX",c("index","cluster","time","strata"),"XX"),
                           paste0("XX",c("index","cluster","time","strata"),".indexXX"),
                           name.time, args$outcome, args$nameL.cores, args$nameW.cores)
            missing.var1 <- setdiff(all.vars(facet), names(object))
            missing.var2 <- setdiff(all.vars(facet), names(attr(object,"wide")))
            missing.var <- c(missing.var1,missing.var2)[which.min(c(length(missing.var1),length(missing.var2)))]
            stop("When a formula, argument \'facet\' should contain variables available in the dataset used to fit the model. \n",
                 "Unknown variable(s): \"",paste(missing.var, collapse = "\" \""),"\" \n")
        }
        
    }else if(type %in% c("partial","partial-center")){
        facet <- NULL
    }

    ## format    
    format <- args$format
    if(is.null(attr(format,"original"))){
        original.format <- format
    }else{
        original.format <- attr(format,"original")
    }
    

    by.repetition <- FALSE
    if(type == "correlation"){

        if(n.time == 1){
            stop("Cannot display the residual correlation over time when there is only a single timepoint. \n")
        }
        if("wide" %in% original.format == FALSE){
            stop("Residuals must be in the wide format to display the residual correlation. \n",
                 "Consider setting the argument \'format\' to \"wide\" when calling residuals(). \n")
        }
        if(format == "long"){
            object <- attr(object,"wide")
            format <- "wide"
        }
        
        
    }else if(type == "qqplot"){
        by.repetition <- any(all.vars(facet) %in% c(name.time,attr(name.time,"original")))

        if(engine.qqplot == "qqtest"){
            if(any(all.vars(facet) %in% c(name.time,attr(name.time,"original")) == FALSE)){
                stop("Can only stratify the display regarding the time variable when using engine.qqplot=\"qqtest\". \n",
                     "time variable: \"",paste(union(name.time,attr(name.time,"original")), collapse = "\", \""),"\"\n",
                     "Consider using engine.qqplot=\"ggplot2\". \n")
            }
            if(by.repetition){
                if("wide" %in% original.format == FALSE){
                    stop("Residuals must be in the wide format to display qqplots per repetition via the qqtest package. \n",
                         "Consider setting the argument \'format\' to \"wide\" when calling residuals(). \n")
                }
                if(format == "long"){
                    object <- attr(object,"wide")
                    format <- "wide"
                }
            }
        }

    }else{
        if("long" %in% original.format == FALSE){
            stop("Residuals must be in the long format to obtain scatterplots of the residuals \n",
                 "or qqplots of the residuals (except when using the qqtest package for repetition specific qqplots). \n",
                 "Consider setting the argument \'format\' to \"long\" when calling residuals(). \n")
        }
        if(format == "wide"){
            object <- attr(object,"long")
            format <- "long"
        }
    }

    if((format == "long") && (name.time %in% names(object)==FALSE) && !is.null(index.time)){
        object[[name.time]] <- index.time
    }

    if(type %in% c("scatterplot","scatterplot2")){
        by.repetition <- any(all.vars(facet) %in% c(name.time,attr(name.time,"original")))

        if(by.repetition){
            if(name.time %in% names(object) == FALSE){
                if(all(is.na(attr(name.time,"original")))){
                    stop("Cannot display a scatterplot of the residuals per repetition without the repetition variable. \n",
                         "Consider specifying argument \'repetition\' when calling lmm(), something like repetition = ~time|id. \n")
                }
                if(length(attr(name.time,"original"))>1){
                    name.time <- attr(name.time,"original")
                }
            }
        }
    }

    if(is.null(color)){ ## necessary when function called from autplot.lmm
        if(type %in% c("partial","partial-center") && n.time>1 && length(setdiff(args$var,"(Intercept)"))==2 && name.time %in% args$var){
            color <- setdiff(args$var, c(name.time,"(Intercept)"))
        }else{
            color <- switch(type,
                            scatterplot = "black",
                            scatterplot2 = "black",
                            partial = c("gray","white"),
                            "partial-center" = c("gray","white"),
                            qqplot = NA,
                            correlation = NA
                            )
        }
    } 
    if(is.null(obs.size)){obs.size <- 1} ## necessary when function called from autplot.lmm
    if(is.null(position)){position <- ggplot2::position_dodge(width = 0.1)} 

    ## ** process input
    label.residual <- switch(type.residual,
                             "fitted" = "Fitted values",
                             "response" = "Raw residuals",
                             "studentized" = "Studentized residuals",
                             "pearson" = "Pearson residuals",
                             "normalized" = "Normalized residuals",
                             "normalized2" = "Pearson Normalized residuals",
                             "scaled" = "Scaled residuals")

    if(format == "long"){
        name.residual <- args$nameL.colres[which(type.residual == args.type)]
    }else if(format == "wide"){
        name.residual <- args$nameW.colres
    }

    formula.time <- paste("~",paste(name.time,collapse="+"))
    if(format == "long" && by.repetition && any(rowSums(is.na(object[name.time]))>0)){
        object <- object[rowSums(is.na(object[name.time]))==0,,drop=FALSE]
    }

    ## ** build graphical display
    if(type=="qqplot"){ ## overall timepoints

        ## *** qqplot
        df.gg <- object
        if(engine.qqplot=="ggplot2"){
            gg <- ggplot2::ggplot(df.gg, ggplot2::aes(sample = .data[[name.residual]]))
            gg <- gg + ggplot2::stat_qq() + ggplot2::stat_qq_line()
            gg <- gg + ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles")
            gg <- gg + ggplot2::ggtitle(label.residual) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
            if(!is.null(facet) & length(all.vars(facet))>0){
                if(attr(stats::terms(facet),"response")==0){
                    gg <- gg + ggplot2::facet_wrap(facet, scales = scales, labeller = labeller, nrow = facet_nrow, ncol = facet_ncol)
                }else{
                    gg <- gg + ggplot2::facet_grid(facet, scales = scales, labeller = labeller, nrow = facet_nrow, ncol = facet_ncol)
                }
            }
        }else if(engine.qqplot=="qqtest"){

            requireNamespace("qqtest")
            if(by.repetition){
                if(is.null(facet_nrow)){
                    facet_nrow <- ceiling(sqrt(n.time))
                }
                if(is.null(facet_ncol)){
                    facet_ncol <- ceiling(n.time/facet_nrow)
                }

                Utime <- attr(object,"reshapeWide")$times
                if(is.character(labeller) && labeller == "label_both"){
                    Utime <- paste0(name.time,": ",Utime)
                }
                
                oldpar <- graphics::par(no.readonly = TRUE)   
                on.exit(graphics::par(oldpar))            
                graphics::par(mfrow = c(facet_nrow,facet_ncol))                
                tempo <- lapply(1:n.time,function(iCol){
                    qqtest::qqtest(stats::na.omit(object[,iCol+1]), main = Utime[iCol])
                    graphics::mtext(label.residual, side = 3)
                })
                
            }else{
                qqtest::qqtest(stats::na.omit(df.gg[[name.residual]]), main = label.residual)
            }

            df.gg <- NULL
            gg <- NULL
        }
    }else if(type.residual %in% c("partial","partial-center")){ ## must be type="scatterplot"

        ## *** partial residual plot
        name.facet <- all.vars(facet)
        names.time <- union(name.time,attr(name.time,"original"))
        ## only variables varying within panel and color
        name.var <- setdiff(args$var,c("(Intercept)",color,setdiff(name.facet,names.time)))
        type.var <- args$type.var[match(name.var,setdiff(args$var,"(Intercept)"))]
        if(sum(type.var=="numeric")>1){
            stop("Cannot simulatenously display partial residuals for more than 1 numeric variable. \n")
        }
        if(length(type.var)>2){
            stop("Cannot simulatenously display partial residuals for more than 2 variables. \n")
        }

        name.fitted <- name.var
        
        ## dataset
        if("fitted.lower" %in% names(object) && "fitted.lower" %in% names(object)){
            ci <- TRUE
        }else{
            ci <- FALSE
        }
        df.gg <- object

        ## identify continuous and categorical covariates
        if(sum(type.var=="numeric")==0){
            name.varnum <- NULL
        }else{
            name.varnum <- name.var[type.var=="numeric"]
        }
        if(sum(type.var=="categorical")==0){
            name.varcat <- NULL
        }else{
            name.varcat0 <- name.var[type.var=="categorical"]
            name.varcat <- paste(name.varcat0,collapse = ", ")
            if(sum(type.var=="categorical")>1){
                df.gg[[name.varcat]] <- nlme::collapse(df.gg[name.varcat0], sep = ",", as.factor = TRUE)
            }
        }

        ## time.var
        if(is.null(time.var)){
            if(all(type.var == "categorical")){
                time.var <- name.varcat
            }else{
                time.var <- name.varnum
            }
        }
        
        if(any(is.na(df.gg$fitted))){
            if(all(type.var == "categorical")){
                ## do not remove NA as lines are drawn over clusters                
            }else if(length(type.var)==1){
                ## remove NAs as the line is drawn over all observations 
                df.gg <- df.gg[!is.na(df.gg$fitted),,drop=FALSE]                
            }else if(length(type.var)==2){
                ## remove NAs as the line is drawn over all observations within covariate value
                df.gg <- df.gg[!is.na(df.gg$fitted),,drop=FALSE]                
            }            
        }

        ## display
        if(all(type.var == "categorical")){
            gg <- ggplot2::ggplot(data = df.gg, mapping = ggplot2::aes(x = .data[[time.var]]))

            ## observations
            if(identical(color,FALSE)){
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial), size = obs.size)
            }else if(color[1] %in% names(df.gg)){
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial, color = .data[[color[1]]]), size = obs.size, position = position)
            }else{
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial), color = color[1], size = obs.size)
            }

            ## fit
            if(identical(color,FALSE)){
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$fitted), size = mean.size[1], shape = 21)
            }else if(utils::tail(color,1) %in% names(df.gg)){
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$fitted, fill = .data[[utils::tail(color,1)]]), size = mean.size[1], shape = 21, position = position)
            }else{
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$fitted), size = mean.size[1], shape = 21, fill = utils::tail(color,1))
            }

            ## uncertainty about the fit
            if(ci){
                gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper))
            }

            ## connect lines
            ## if cluster variable available, time is an x-variable not in facet
            if((args$name.cluster %in% names(df.gg)) && (any(name.varcat %in% name.time) && all(all.vars(facet) %in% name.time == FALSE))){
                if(utils::tail(color,1) %in% names(df.gg)){
                    gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$fitted, group = .data[[args$name.cluster]], color = .data[[utils::tail(color,1)]]),
                                                  linewidth = mean.size[2], position = position)
                }else{
                    gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$fitted, group = .data[[args$name.cluster]]), linewidth = mean.size[2])
                }            
            }

        }else if(length(type.var)==1){
            gg <- ggplot2::ggplot(data = df.gg, mapping = ggplot2::aes(x = .data[[time.var]]))
            if(identical(color,FALSE)){
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial), size = obs.size)
            }else if(color[1] %in% names(df.gg)){
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial, color = .data[[color[1]]]), size = obs.size, position = position)
            }else{
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial), color = color[1], size = obs.size)
            }
            gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$fitted), linewidth = mean.size[2])
            if(ci){
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$fitted.lower, ymax = .data$fitted.upper), alpha = ci.alpha)
            }
        }else if(length(type.var)==2){
            gg <- ggplot2::ggplot(data = df.gg, mapping = ggplot2::aes(x = .data[[time.var]]))
            gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial, color = .data[[name.varcat]]), color = color[1], size = obs.size)
            if(identical(color,FALSE)){
                gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$fitted, group = .data[[name.varcat]]), linewidth = mean.size[2])
            }else{
                gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$fitted, group = .data[[name.varcat]], color = .data[[name.varcat]]), linewidth = mean.size[2])
            }
            if(ci){
                if(identical(color,FALSE)){
                    gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$fitted.lower, ymax = .data$fitted.upper, group = .data[[name.varcat]]), alpha = ci.alpha)
                }else{
                    gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$fitted.lower, ymax = .data$fitted.upper, group = .data[[name.varcat]], color = .data[[name.varcat]]), alpha = ci.alpha)
                }
            }
        }        
        if(!is.null(facet) & length(all.vars(facet))>0){
            if(attr(stats::terms(facet),"response")==0){
                gg <- gg + ggplot2::facet_wrap(facet, scales = scales, labeller = labeller, nrow = facet_nrow, ncol = facet_ncol)
            }else{
                gg <- gg + ggplot2::facet_grid(facet, scales = scales, labeller = labeller, nrow = facet_nrow, ncol = facet_ncol)
            }
        }
        reference <- attr(object,"reference")[,!is.na(attr(object,"reference")),drop=FALSE]
        if(NCOL(reference)==0){
            if(args$intercept){
                if("(Intercept)" %in% args$var){
                    gg <- gg + ggplot2::ylab(args$outcome)
                }else{
                    gg <- gg + ggplot2::ylab(paste0(args$outcome, " (centered)"))
                }
            }else{
                gg <- gg + ggplot2::ylab(args$outcome)
            }
        }else{
            reference <- lapply(reference, function(iRef){if(is.factor(iRef)){as.character(iRef)}else{iRef}})
            gg <- gg + ggplot2::ggtitle(paste0("Counterfactual: ",paste(paste0(names(reference),"=",reference), collapse = ", ")))
            gg <- gg + ggplot2::labs(x = paste0(time.var, " (observed)"), y = paste0(args$outcome, " (counterfactual)")) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
        }

        gg <- gg + ggplot2::theme(text = ggplot2::element_text(size=size.text))

    }else if(type %in% c("scatterplot","scatterplot2")){ ## overall timepoints

        ## *** residual plot
        name.fitted <- "fitted"
        if(is.null(time.var)){
            time.var <- name.fitted
            xlab.gg <- "Fitted values"
        }else{
            xlab.gg <- time.var
        }
        
        gg <- ggplot2::ggplot(data = object, mapping = ggplot2::aes(x = .data[[time.var]]))
        gg <- gg + ggplot2::xlab(xlab.gg) + ggplot2::theme(text = ggplot2::element_text(size=size.text))

        if(type == "scatterplot"){            
            gg <- gg + ggplot2::geom_abline(slope=0,intercept=0,color ="red")
            if(color %in% names(object)){
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data[[name.residual]], color = .data[[color]]), size = obs.size, position = position)                
            }else{
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data[[name.residual]]), color = color, size = obs.size)
            }
            gg <- gg + ggplot2::ylab(label.residual) 
            if(add.smooth[1]){
                gg <- gg + ggplot2::geom_smooth(ggplot2::aes(y = .data[[name.residual]]), se = add.smooth[2])
            }
        }else if(type == "scatterplot2"){
            label.residual2 <- paste0("|",label.residual,"|")
            if(color %in% names(object)){
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = sqrt(abs(.data[[name.residual]])), color = .data[[color]]), size = obs.size)
            }else{
                gg <- gg + ggplot2::geom_point(ggplot2::aes(y = sqrt(abs(.data[[name.residual]]))), color = color, size = obs.size)
            }
            gg <- gg + ggplot2::ylab(bquote(sqrt(.(label.residual2))))
            if(add.smooth[1]){
                gg <- gg + ggplot2::geom_smooth(ggplot2::aes(y = sqrt(abs(.data[[name.residual]]))), se = add.smooth[2])
            }
        }

        if(!is.null(facet) & length(all.vars(facet))>0){
            if(attr(stats::terms(facet),"response")==0){
                gg <- gg + ggplot2::facet_wrap(facet, scales = scales, labeller = labeller, nrow = facet_nrow, ncol = facet_ncol)
            }else{
                gg <- gg + ggplot2::facet_grid(facet, scales = scales, labeller = labeller, nrow = facet_nrow, ncol = facet_ncol)
            }
        }
        df.gg <- NULL

    }else if(type == "correlation"){

        ## *** residual correlation heatmap
        M.cor  <- stats::cor(object[,-1], use = "pairwise")
        arr.ind.cor <- rbind(which(is.na(M.cor*NA), arr.ind = TRUE))
        arr.ind.cor[] <- name.residual[arr.ind.cor]
        
        df.gg <- data.frame(correlation = as.double(M.cor), arr.ind.cor, stringsAsFactors = FALSE)
        label.time <- names(name.residual)
        df.gg$col <- factor(df.gg$col, levels = name.residual, labels = label.time)
        df.gg$row <- factor(df.gg$row, levels = name.residual, labels = label.time)
        df.gg$row.index <- match(df.gg$row, label.time)
        df.gg$col.index <- match(df.gg$col, label.time)
        dfR.gg <- df.gg[df.gg$col.index>=df.gg$row.index,,drop=FALSE]
        gg <- ggplot2::ggplot(dfR.gg, ggplot2::aes(x = .data$col,
                                                   y = .data$row,
                                                   fill = .data$correlation)) 
        gg <- gg + ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "blue",
                                                                        high = "red",
                                                                        mid = "white",
                                                                        midpoint = 0,
                                                                        limit = c(-1,1),
                                                                        space = "Lab",
                                                                        name="Correlation")
        gg <- gg + ggplot2::labs(x = NULL, y = NULL) + ggplot2::ggtitle(label.residual)
        gg <- gg + ggplot2::theme(text = ggplot2::element_text(size=size.text))
        if(!is.na(digits.cor) && digits.cor>0){
            gg <- gg + ggplot2::geom_text(ggplot2::aes(label = round(.data$correlation,digits.cor)))
        }

    }

    ## ** export
    return(list(data = df.gg,
                plot = gg))

}


## * autoplot.summarize (documentation)
##' @title Graphical Display of the Descriptive Statistics
##' @description Graphical representation of the descriptive statistics.
##' 
##' @param object,x an object of class \code{summarize}, output of the \code{summarize} function.
##' @param type [character] the summary statistic that should be displayed: \code{"mean"}, \code{"sd"}, \ldots
##' @param variable [character] type outcome relative to which the summary statistic should be displayed.
##' Only relevant when multiple variables have been used on the left hand side of the formula when calling \code{summarize}.
##' @param size.text [numeric, >0] size of the text in the legend, x- and y- labels.
##' @param linewidth [numeric, >0] thickness of the line connecting the points.
##' @param size [numeric, >0] width of the points.
##' @param ... additional arguments passed to .ggHeatmap when displaying the correlation: \itemize{
##' \item name.time [character] title for the x- and y- axis.
##' \item digits.cor [integer, >0] number of digits used to display the correlation.
##' \item name.legend [character] title for the color scale.
##' \item title [character] title for the graph.
##' \item scale [function] color scale used for the correlation.
##' \item type.cor [character] should the whole correlation matrix be displayed (\code{"both"}),
##' or only the element in the lower or upper triangle (\code{"lower"}, \code{"upper"}).
##' \item args.scale [list] arguments to be passed to the color scale.
##' }
##'  
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to generate the plot.
##' \item \code{plot}: ggplot object.
##' }
##' 
##' @keywords hplot
##'
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##' dtS <- summarize(weight ~ time, data = gastricbypassL)
##' plot(dtS)
##' dtS <- summarize(glucagonAUC + weight ~ time|id, data = gastricbypassL, na.rm = TRUE)
##' plot(dtS, variable = "glucagonAUC")
##' plot(dtS, variable = "glucagonAUC", type = "correlation", size.text = 1)

## * autoplot.summarize (code)
##' @export
autoplot.summarize <- function(object, type = "mean", variable = NULL,
                               size.text = 16, linewidth = 1.25, size = 3,
                               ...){

    ## ** normalize input
    name.X <- attr(object, "name.X")
    name.Y <- attr(object, "name.Y")
    name.id <- attr(object, "name.id")
    name.time <- attr(object, "name.time")
    if(length(name.time)==0){
        if(length(name.X)>1){
            stop("Unknown time variable: cannot provide graphical display. \n",
                 "Consider indicating the cluster variable in the formula when calling summarize. \n",
                 "Something like Y ~ time | id or Y ~ time + group | id. \n")
        }else{
            name.time <- name.X
        }
    }
    name.stat <- setdiff(names(object), c("outcome",name.X,name.Y))
    correlation <- attr(object,"correlation")

    ## variable
    if(is.null(variable)){
        if(length(name.Y)>1){
            stop("Missing patterns for several variables could be displayed. \n",
                 "Consider specifying which one via the argument \'variable\'. \n")
        }else{
            variable <- name.Y
        }
    }else if(length(variable)!=1){
        stop("Argument \'variable\' should have length 1. \n")
    }else{
        if(variable %in% unique(object$outcome) == FALSE){
            stop("Incorrect value passed to argument \'variable\'. \n",
                 "Possible values: \"",paste(unique(object$outcome), collapse = "\", \""),"\"\n")
        }        
    }

    ## object
    object <- object[object$outcome==variable,]
    if(is.null(correlation)){
        name.stat.all <- name.stat
    }else{
        name.stat.all <- c(name.stat,"correlation")
        correlation <- correlation[[variable]]
    }

    ## type
    type <- match.arg(type, name.stat.all)

    ## ** graph
    if(type == "correlation"){
        gg <- .ggHeatmap(correlation, name.time = name.time, size.text = size.text, ...)
    }else{
        dots <- list(...)
        if(length(dots)>0){
            stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
        }
        name.group <- setdiff(name.X,name.time)
        if(length(name.group)==0){
            gg <- ggplot2::ggplot(object, ggplot2::aes(x = .data[[name.time]], y = .data[[type]], group = ""))
        }else if(length(name.group)==1){
            gg <- ggplot2::ggplot(object, ggplot2::aes(x = .data[[name.time]], y = .data[[type]],
                                                       group = .data[[name.group]], color = .data[[name.group]]))
        }else{
            name.Group <- nlme::collapse(name.group, as.factor = TRUE)
            object[[name.Group]] <- nlme::collapse(object[,name.group], as.factor = TRUE)
            gg <- ggplot2::ggplot(object, ggplot2::aes(x = .data[[name.time]], y = .data[[type]],
                                                       group = .data[[name.Group]], color = .data[[name.Group]]))
        }
        gg <- gg + ggplot2::geom_point(size = size) + ggplot2::geom_line(linewidth = linewidth)
        if(type == "missing"){
            gg <- gg + ggplot2::labs(y = paste0("number of missing values in ",variable))
        }else if(type == "pc.missing"){
            gg <- gg + ggplot2::labs(y = paste0("percentage of missing values in ",variable))
            gg <- gg + ggplot2::scale_y_continuous(labels = scales::percent)
        }else{
            gg <- gg + ggplot2::labs(y = variable)
        }
        if(!is.null(size.text)){
            gg <- gg + ggplot2::theme(text = ggplot2::element_text(size=size.text))
        }
        data <- NULL
    }


    ## ** export
    return(list(data = data,
                plot = gg))
}

## * autoplot.summarizeNA (documentation)
##' @title Graphical Display of Missing Data Pattern
##' @description Graphical representation of the possible missing data patterns in the dataset.
##'
##' @param object,x a \code{summarizeNA} object, output of the \code{\link{summarizeNA}} function.
##' @param variable [character] variable for which the missing patterns should be displayed.
##' Only required when the argument \code{repetition} has been specified when calling \code{summarizeNA}.
##' @param size.text [numeric, >0] size of the font used to display text.
##' @param add.missing [logical] should the number of missing values per variable be added to the x-axis tick labels.
##' @param order.pattern [numeric vector or character] in which order the missing data pattern should be displayed. Can either be a numeric vector indexing the patterns or a character refering to order the patterns per number of missing values (\code{"n.missing"}) or number of observations (\code{"frequency"}).
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to create the graphical display.
##' \item \code{plot}: ggplot object.
##' }
##'
##' @keywords hplot

## * autoplot.summarizeNA (code)
##' @export
autoplot.summarizeNA <- function(object, variable = NULL, size.text = 16,
                                 add.missing = " missing", order.pattern = NULL, ...){

    newnames <- attr(object,"args")$newnames
    keep.data <- attr(object,"args")$keep.data
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.null(object$variable) && length(unique(object$variable))>1){
        if(is.null(variable)){
            stop("Missing patterns for several variables could be displayed. \n",
                 "Consider specifying which one via the argument \'variable\'. \n")
        }else if(length(variable)!=1){
            stop("Argument \'variable\' should have length 1. \n")
        }else{
            if(variable %in% object$variable == FALSE){
                stop("Incorrect value passed to argument \'variable\'. \n",
                     "Possible values: \"",paste(unique(object$variable), collapse = "\", \""),"\"\n")
            }
            object <- object[object$variable == variable,,drop=FALSE]
        }
    }
    if(identical(order.pattern,newnames[4])){
        order.pattern <- order(object[[newnames[4]]])
    }
    if(identical(order.pattern,newnames[2])){
        order.pattern <- order(object[[newnames[2]]])
    }
    
    if(keep.data == FALSE){
        stop("Argument \'keep.data\' should be set to TRUE when calling summarizeNA to obtain a graphical display. \n")
    }

    keep.cols <- setdiff(names(object),c(newnames,all.vars(attr(object,"args")$repetition)))
    data <- as.data.frame(object[,c(newnames[3],keep.cols),drop=FALSE])
    dataL <- stats::reshape(data, direction = "long", idvar = newnames[3], varying = keep.cols,
                            v.names = newnames[[2]],
                            timevar = newnames[[1]])

    nObs.pattern <- stats::setNames(object[[newnames[2]]], object[[newnames[3]]])
    nNA.Var <- colSums(sweep(data[,keep.cols,drop=FALSE], FUN = "*", MARGIN = 1, STATS = nObs.pattern))

    if(!is.null(add.missing) && !is.na(add.missing) && !identical(FALSE,add.missing)){
        dataL[[newnames[1]]] <- factor(dataL[[newnames[1]]], labels = paste0(keep.cols,"\n(",nNA.Var,add.missing,")"))
    }else{
        dataL[[newnames[1]]] <- factor(dataL[[newnames[1]]], labels = keep.cols)
    }
    dataL[[newnames[2]]] <- factor(dataL[[newnames[2]]], levels = 1:0, labels = c("yes","no"))

    if(!is.null(order.pattern)){
        if(length(order.pattern)!=NROW(data)){
            stop("Argument \'order.pattern\' should have length ",NROW(data),".\n",sep="")
        }
        if(any(sort(order.pattern)!=1:NROW(data))){
            stop("Argument \'order.pattern\' should be a vector containing integers from 1 to ",NROW(data),".\n",sep="")
        }
        dataL[[newnames[3]]] <- factor(dataL[[newnames[3]]], levels = data[[newnames[3]]][order.pattern])
    }

    gg.NA <- ggplot2::ggplot(dataL, ggplot2::aes(y = .data[[newnames[3]]], x = .data[[newnames[1]]], fill = .data[[newnames[2]]]))
    gg.NA <- gg.NA + ggplot2::geom_tile(color = "black")
    gg.NA <- gg.NA + ggplot2::scale_y_discrete(breaks = unique(dataL[[newnames[3]]]),
                                               labels = nObs.pattern[unique(dataL[[newnames[3]]])])
    gg.NA <- gg.NA + ggplot2::labs(fill = "missing", x = "", y = "number of observations")
    gg.NA <- gg.NA + ggplot2::theme(text = ggplot2::element_text(size=size.text))

    ## ** export
    return(list(data = dataL,
                plot = gg.NA))
}

## * autoplot.Wald_lmm (documentation)
##' @title Graphical Display For Linear Hypothesis Test
##'
##' @param object,x a \code{Wald_lmm} object.
##' @param type [character] what to display: a forest plot (\code{"forest"}) or a heatmap (\code{"heat"}).
##' @param add.args [list] additional arguments used to customized the graphical display.
##' Must be a named list. See details.
##' @param size.text [numeric, >0] size of the font used to display text.
##' @param ... arguments passed to the confint method.
##'
##' @details Argument \strong{add.args}: parameters specific to the forest plot: \itemize{
##' \item \code{color}: [logical] should the estimates be colored by global null hypothesis, e.g. when testing the effect of a 3 factor covariate, the two corresponding coefficient will have the same color. Alternatively a vector of positive integers giving the color with which each estimator should be displayed.
##' \item \code{color}: [logical] should the estimates be represented by a different shape per global null hypothesis, e.g. when testing the effect of a 3 factor covariate, the two corresponding coefficient will have the same type of point. Alternatively a vector of positive integers describing the shape to be used for each estimator.
##' \item \code{ci}: [logical] should confidence intervals be displayed?
##' \item \code{size.estimate}: [numeric, >0] size of the dot used to display the estimates.
##' \item \code{size.ci}: [numeric, >0] thickness of the line used to display the confidence intervals.
##' \item \code{width.ci}: [numeric, >0] width of the line used to display the confidence intervals.
##' \item \code{size.null}: [numeric, >0] thickness of the line used to display the null hypothesis. 
##' }
##' Parameters specific to the heatmap plot: \itemize{
##' \item \code{limits}: [numeric vector of length 2] minimum and maximum value of the colorscale relative to the correlation.
##' \item \code{low}, \code{mid}, \code{high}: [character] color for the the colorscale relative to the correlation.
##' \item \code{midpoint}: [numeric] correlation value associated with the color defined by argument \code{mid}
##' }
##'
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to create the graphical display.
##' \item \code{plot}: ggplot object.
##' }
##'
##' @examples
##' ## From the multcomp package
##' if(require(datasets) && require(ggplot2)){
##'
##' ## only tests with 1 df
##' ff <- Fertility ~ Agriculture + Examination + Education + Catholic + Infant.Mortality
##' e.lmm <- lmm(ff, data = swiss)
##' e.aovlmm <- anova(e.lmm)
##' 
##' autoplot(e.aovlmm, type = "forest")
##' autoplot(e.aovlmm, type = "heat") ## 3 color gradient
##' autoplot(e.aovlmm, type = "heat", add.args = list(mid = NULL)) ## 2 color gradient
##'
##' ## test with more than 1 df
##' e.lmm2 <- lmm(breaks ~ tension + wool, data = warpbreaks)
##' e.aovlmm2 <- anova(e.lmm2)
##' autoplot(e.aovlmm2)
##' autoplot(e.aovlmm2, add.args = list(color = FALSE, shape = FALSE))
##' }
##'
##' @keywords hplot

## * autoplot.Wald_lmm (code)
##' @export
autoplot.Wald_lmm <- function(object, type = "forest", size.text = 16, add.args = NULL, ...){

    ## ** check user input
    type <- match.arg(type, c("forest","heat"))
    if(!is.null(add.args) && !is.list(add.args)){
        stop("Argument \'add.args\' should be a list. \n")
    }
    if(is.list(add.args) && is.null(names(add.args))){
        stop("Argument \'add.args\' should have names, i.e. names(add.args) should not be NULL. \n")
    }

    if(type=="forest"){
        init.add.args <- list(color = NULL,
                              shape = NULL,
                              ci = TRUE,
                              size.estimate = 3, 
                              size.ci = 1,
                              width.ci = 0.2,
                              size.null = 1)
    }else{
        init.add.args <- list(limits = c(-1,1.00001),
                              low = "blue",
                              mid = "white",
                              high = "red",
                              midpoint = 0,
                              value.text = FALSE,
                              value.round = 2,
                              value.size = 5)
    }
    valid.names <- names(init.add.args)
    if(any(names(add.args) %in% valid.names == FALSE)){
        invalid.names <- names(add.args)[names(add.args) %in% valid.names == FALSE]
        if(length(invalid.names)>1){txt.invalid <- "are not valid arguments"}else{txt.invalid <- "is not a valid argument"}

        possible.names <- setdiff(valid.names,names(add.args))
        if(length(possible.names)>0){
            txt.valid <- paste0("Possible arguments: \"",paste(possible.names, collapse = "\", \""),"\". \n")
        }else{
            txt.valid <- NULL
        }

        stop("Incorrect element in argument \'add.args\'. \n",
             "\"",paste(invalid.names, collapse = "\", \""),"\" ",txt.invalid,". \n",
             txt.valid, sep = "")
    }

    if(is.null(add.args)){
        add.args <- init.add.args
    }else{
        missing.names <- setdiff(valid.names, names(add.args))
        if(length(missing.names)>0){
            add.args[missing.names] <- init.add.args[missing.names]
        }
    }

    ## ** graphical display
    if(type=="forest"){

        color <- add.args$color
        shape <- add.args$shape
        ci <- add.args$ci
        size.estimate <- add.args$size.estimate
        size.ci <- add.args$size.ci
        width.ci <- add.args$width.ci
        size.null <- add.args$size.null
        rhs <- unique(object$univariate$null)

        if(ci){
            table <- confint(object, columns = c("estimate","test","lower","upper"), ...)
        }else{
            table <- object$univariate
        }
        table <- cbind(names = rownames(table), table)
        if(is.null(color)){
            if(length(unique(table$test))==1 || all(duplicated(table$test)==FALSE)){
                color <- FALSE
                color.legend <- FALSE
            }else{
                table$color <- table$test
                color <- TRUE
                color.legend <- FALSE
            }
        }else{
            if(identical(color,FALSE)){
                color <- FALSE
                color.legend <- FALSE
            }else if(length(color)==NROW(table) && all(is.character(color))){
                table$color <- color
                color <- TRUE
                color.legend <- TRUE
            }else{
                stop("Argument \'color\' should have length ",NROW(table)," and be of type character. \n")
            }
        }
        if(is.null(shape)){
            if(length(unique(table$test))==1 || all(duplicated(table$test)==FALSE)){
                shape <- FALSE
                shape.legend <- FALSE
            }else{
                table$shape <- table$test
                shape <- TRUE
                shape.legend <- FALSE
            }
        }else{
            if(identical(shape,FALSE)){
                shape <- FALSE
                shape.legend <- FALSE
            }else if(length(shape)==NROW(table) && all(is.numeric(shape))){
                table$shape <- as.character(shape)
                shape <- TRUE
                shape.legend <- TRUE
            }else{
                stop("Argument \'shape\' should have length ",NROW(table)," and be of type numeric. \n")
            }
        }

        table$test <- as.factor(table$test)
        table$names <- factor(table$names, levels = unique(table$names)) ## ensure same ordering as in the object (instead of alphabetical ordering)
        if(color & shape){
            gg <- ggplot2::ggplot(table, ggplot2::aes(x = .data$names, y = .data$estimate, color = .data$color, shape = .data$shape)) + ggplot2::labs(color = "", shape = "")
        }else if(color){
            gg <- ggplot2::ggplot(table, ggplot2::aes(x = .data$names, y = .data$estimate, color = .data$color)) + ggplot2::labs(color = "")
        }else if(shape){
            gg <- ggplot2::ggplot(table, ggplot2::aes(x = .data$names, y = .data$estimate, shape = .data$shape)) + ggplot2::labs(color = "")
        }else{
            gg <- ggplot2::ggplot(table, ggplot2::aes(x = .data$names, y = .data$estimate))
        }
    
    if(shape.legend){
        gg <- gg + ggplot2::scale_shape_manual(values = as.numeric(unique(table$shape)), breaks = unique(table$shape)) + ggplot2::guides(shape = "none")
    }
    if(color.legend){
        gg <- gg + ggplot2::scale_color_manual(values = unique(table$color), breaks = unique(table$color)) + ggplot2::guides(color = "none")
    }
    gg <- gg + ggplot2::geom_point(size = size.estimate) + ggplot2::labs(x = "", y = "")
    if(ci){
        gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), linewidth = size.ci, width = width.ci)
    }
    if(size.null>0 && length(rhs)==1){
        gg <- gg + ggplot2::geom_hline(yintercept=rhs, lty=2, linewidth = size.null)
    }
    gg <- gg + ggplot2::coord_flip()

    }else if(type=="heat"){

        limits <- add.args$limits
        low <- add.args$low
        mid <- add.args$mid
        high <- add.args$high
        midpoint <- add.args$midpoint
        value.text <- add.args$value.text
        value.round <- add.args$value.round
        value.size <- add.args$value.size

        Sigma_t <- stats::cov2cor(object$vcov)
        ## from matrix to long format
        table <- as.data.frame(cbind(which(is.na(NA*Sigma_t), arr.ind = TRUE), value = as.numeric(Sigma_t)))
        rownames(table) <- NULL

        ## rename
        if(!is.null(object$args$sep) && all(colnames(Sigma_t) == rownames(Sigma_t))){
            splitname <- strsplit(colnames(Sigma_t),split = object$args$sep, fixed = TRUE)
            if(all(lengths(splitname)==2) && length(unique(sapply(splitname,"[",2)))==1){
                name.x <- splitname[[1]][2]
                name.y <- splitname[[1]][2]
                table$col <- sapply(splitname,"[",1)[table$col]
                table$row <- sapply(splitname,"[",1)[table$row]
            }else{
                table$col <- colnames(Sigma_t)[table$col]
                table$row <- rownames(Sigma_t)[table$row]
                name.x <- ""
                name.y <- ""
            }
        }else{
            table$col <- colnames(Sigma_t)[table$col]
            table$row <- rownames(Sigma_t)[table$row]
            name.x <- ""
            name.y <- ""
        }
        table$row <- factor(table$row, levels = unique(table$row))
        table$col <- factor(table$col, levels = rev(levels(table$row)))
        table$rvalue <- round(table$value, digits = value.round)
        gg <- ggplot2::ggplot(table) + ggplot2::geom_tile(ggplot2::aes(x = .data$row, y = .data$col, fill = .data$value))

        if(value.text){
            gg <- gg + ggplot2::geom_text(ggplot2::aes(x = .data$row, y = .data$col, label = .data$rvalue), size = value.size)
        }
        if(!is.null(mid)){
            gg <- gg + ggplot2::scale_fill_gradient2(limits = limits, midpoint = midpoint, low = low, mid = mid, high = high)
        }else{
            gg <- gg + ggplot2::scale_fill_gradient(limits = limits, low = low, high = high)
        }
        gg <- gg + ggplot2::labs(x=name.x,y=name.y, fill = "correlation")
    }
        
    gg <- gg + ggplot2::theme(text = ggplot2::element_text(size=size.text))

    ## ** export
    return(list(data = table,
                plot = gg))
}
## * .ggHeatmap
.ggHeatmap <- function(object, name.time, size.text = 16, digits.cor = 2, name.legend = "Correlation", title = NULL,
                       scale = ggplot2::scale_fill_gradient2, type.cor = "both",
                       args.scale = list(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1))
                       ){

    ## ** from matrix/list format to data.frame
    type.cor <- match.arg(type.cor, c("lower","upper","both"))
    mat2df <- function(mat){
        df <- as.data.frame(mat)
        name.col <- names(df)
        ind.cor <- !is.na(df)
        arr.ind.cor <- which(ind.cor, arr.ind = TRUE)
        arr.ind.cor[] <- name.col[arr.ind.cor]

        df.gg <- data.frame(correlation = df[ind.cor], arr.ind.cor, stringsAsFactors = FALSE)
        df.gg$col <- factor(df.gg$col, levels = name.col, labels = name.col)
        df.gg$row <- factor(df.gg$row, levels = rev(name.col), labels = rev(name.col))
        df.gg$row.index <- match(df.gg$row, name.col)
        df.gg$col.index <- match(df.gg$col, name.col)
        if(type.cor=="lower"){
            df.gg <- df.gg[df.gg$col.index<=df.gg$row.index,,drop=FALSE]
            rownames(df.gg) <- NULL
        }else if(type.cor=="upper"){
            df.gg <- df.gg[df.gg$col.index>=df.gg$row.index,,drop=FALSE]
            rownames(df.gg) <- NULL
        }
        return(df.gg)
    }
    if(is.matrix(object) || is.data.frame(object)){
        name.object <- NULL
        data <- mat2df(object)
    }else if(is.list(object) && length(object) == 1 && is.null(names(object))){
        name.object <- NULL
        data <- mat2df(object[[1]])
    }else if(is.list(object)){
        name.object <- names(object)
        data <- do.call(rbind,lapply(name.object, function(iName){
            cbind(iName, mat2df(object[[iName]]))
        }))        
    }else {
        stop("Unknown data type: should be matrix, data.frame, or list. \n")
    }

    ## ** graphical display
    gg <- ggplot2::ggplot(data, ggplot2::aes(x = .data$col,
                                             y = .data$row,
                                             fill = .data$correlation)) 
    gg <- gg + ggplot2::geom_tile()
    if(!is.null(name.object)){
        gg <- gg + ggplot2::facet_wrap(~iName)
    }
    if(!is.null(scale)){
        gg <- gg + do.call(scale, args.scale)
    }
    if(!is.null(name.time)){
        gg <- gg + ggplot2::labs(x = name.time, y = name.time)
    }
    if(!is.null(name.legend)){
        gg <- gg + ggplot2::labs(fill = name.legend)
    }
    if(!is.null(title)){
        gg <- gg + ggplot2::ggtitle(title)
    }
    if(!is.null(size.text)){
        gg <- gg + ggplot2::theme(text = ggplot2::element_text(size=size.text))
    }
    if(!is.null(digits.cor) && !is.na(digits.cor) && digits.cor>0){
        gg <- gg + ggplot2::geom_text(ggplot2::aes(label = round(.data$correlation,digits.cor)))
    }

    ## ** export
    return(gg)

}


##----------------------------------------------------------------------
### autoplot.R ends here
