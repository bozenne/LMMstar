### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  8 2021 (00:01) 
## Version: 
## Last-Updated: Jun  2 2022 (17:10) 
##           By: Brice Ozenne
##     Update #: 181
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
##' @param time.var [character] x-axis variable for the plot.
##' @param obs.alpha [numeric, 0-1] When not NA, transparency parameter used to display the original data by cluster.
##' @param obs.size [numeric vector of length 2] size of the point and line for the original data.
##' @param color [character] name of the variable in the dataset used to color the curve.
##' @param ci [logical] should confidence intervals be displayed?
##' @param ci.alpha [numeric, 0-1] When not NA, transparency parameter used to display the confidence intervals.
##' @param plot [logical] should the plot be displayed?
##' @param mean.size [numeric vector of length 2] size of the point and line for the mean trajectory.
##' @param size.text [numeric, >0] size of the font used to displayed text when using ggplot2.
##' @param position.errorbar [character] relative position of the errorbars.
##' @param ... arguments passed to the predict method.
##'
##' @return A list with two elements \itemize{
##' \item \code{data}: data used to create the graphical display.
##' \item \code{plot}: ggplot object.
##' }


## * autoplot.lmm (code)
##' @rdname autplot
##' @export
autoplot.lmm <- function(object, obs.alpha = 0, obs.size = c(2,0.5), at = NULL, time.var = NULL, color = TRUE, ci = TRUE, ci.alpha = NA, plot = TRUE,
                         mean.size = c(3, 1), size.text = 16, position.errorbar = "identity", ...){

    if(object$time$n==1){
        stop("Cannot display the fitted values over time when there only is a single timepoint. \n")
    }

    ## ** extract from object
    object.data <- object$data
    Upattern <- object$design$vcov$X$Upattern
    pattern.cluster <- object$design$vcov$X$pattern.cluster

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
            stop("Argument \'time.var\' should have length 1 or be NULL. \n")
        }
        if(time.var %in% names(object$data) == FALSE){
            stop("Could not find the variable \"",time.var,"\" defined by the \'time.var\' argument in the data used to fit the lmm. \n")
        }
        time.var.plot <- time.var
        xlabel.plot <- time.var
    }
    time.var <- attr(object$time$var,"original") ## need to be after statement on time.var.plot to avoid confusion
    mu.var <- attr(object$design$mean,"variable") 
    if(length(time.var) == 0 && length(mu.var) == 0){
        message("There is nothing to be displayed: empty time variable and no covariate for the mean structure. \n")
        return(NULL)
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
    X.beta <- stats::model.matrix(object, data = data[,union(time.var, mu.var), drop=FALSE], effects = "mean")
    IX.beta <- interaction(as.data.frame(X.beta), drop = TRUE)
    vec.X.beta <- tapply(IX.beta, data[["XXclusterXX"]],paste, collapse = "_XXX_")
    UX.beta <- unique(vec.X.beta)
    if(ci && !is.na(ci.alpha)){ ## remove duplicates
        
        test.UX.beta <- rep(TRUE, length(UX.beta))
        for(iUX in 1:length(UX.beta)){

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
    
    ## pattern.alltimes <- which(object$time$n == sapply(object$design$vcov$X$Upattern$time, length))
    ## if(length(pattern.alltimes)>0){
    ##     index.id <- sapply(pattern.alltimes, function(iP){names(which(object$design$vcov$X$pattern.cluster == object$design$vcov$X$Upattern$name[iP]))[1]})
    ##     newdata <- data[data[["XXclusterXX"]] %in% index.id,]
    ## }else{
    ##     warning("No cluster with all timepoints. Generate an artificial cluster. \n")
    ##     X.vars <- all.vars(stats::delete.response(stats::terms(object$formula$mean)))

    ##     data.order <- data[rowSums(is.na(data[X.vars]))==0,] ## rm NA in X
    ##     data.order <- data.order[order(data.order[["XXclusterXX"]],data.order[["XXtimeXX"]]),] ## order data
    ##     newdata <- data.order[!duplicated(data.order[["XXtimeXX"]]),,drop=FALSE]
    ## }

    if(identical(color,TRUE)){
        mean.var <- all.vars(stats::delete.response(stats::terms(stats::formula(object, effects = "mean"))))
        if(length(mean.var)>0){
            newdataRed <- newdata[order(newdata[["XXclusterXX"]]),mean.var,drop=FALSE]
            order.cluster <- droplevels(newdata[["XXclusterXX"]][order(newdata[["XXclusterXX"]])])

            M.duplicated <- apply(newdataRed, 2, function(iCol){unlist(tapply(iCol, order.cluster, function(iColCluster){duplicated(iColCluster)[-1]}))})
            color <- names(which(colSums(M.duplicated)==NROW(M.duplicated)))
        }else{
            color <- NULL
        }
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
    if(!is.na(obs.alpha) && obs.alpha>0 && length(color)>1 && color %in% names(data) == FALSE){
        ls.UX <- lapply(as.character(unique(newdata[["XXclusterXX"]])), function(iC){
            iVec <- as.character(interaction(data[data[["XXclusterXX"]] %in% iC,mu.var,drop=FALSE]))
            cbind(repetition = data[data[["XXclusterXX"]] %in% iC,"XXtimeXX",drop=FALSE], lp = iVec)
        })
    
        index.X <- unlist(lapply(as.character(unique(data[["XXclusterXX"]])), function(iC){
            iVec <- as.character(interaction(data[data[["XXclusterXX"]] %in% iC,mu.var,drop=FALSE]))
            iM <- cbind(repetition = data[data[["XXclusterXX"]] %in% iC,"XXtimeXX",drop=FALSE], lp = iVec)
            iScore <- unlist(lapply(ls.UX, function(iUX){sum(iUX[match(iM[,"Days"],iUX[,"Days"]),"lp"]==iM[,"lp"])}))
            which.max(iScore)
        }))
        data[[color]] <- sort(unique(newdata[[color]]))[index.X]
    }

    ## ** compute fitted curve
    if(!is.na(obs.alpha) && obs.alpha>0){
        preddata <- cbind(data, stats::predict(object, newdata = data[,union(time.var, mu.var), drop=FALSE], ...))
    }else{
        preddata <- cbind(newdata, stats::predict(object, newdata = newdata[,union(time.var, mu.var), drop=FALSE], ...))
    }

    ## ** generate plot
    gg <- ggplot2::ggplot(preddata, ggplot2::aes_string(x = time.var.plot, y = "estimate", group = "XXclusterXX"))
    if(!is.na(obs.alpha) && obs.alpha>0){
        if(!is.null(color)){
            gg <- gg + ggplot2::geom_point(data = data, mapping = ggplot2::aes_string(x = time.var.plot, y = outcome.var, group = "XXclusterXX", color = color),
                                           alpha = obs.alpha, size = obs.size[1])
            gg <- gg + ggplot2::geom_line(data = data, mapping = ggplot2::aes_string(x = time.var.plot, y = outcome.var, group = "XXclusterXX", color = color),
                                          alpha = obs.alpha, size = obs.size[2])
            ## gg + facet_wrap(~XXclusterXX)
        }else{
            gg <- gg + ggplot2::geom_point(data = data, mapping = ggplot2::aes_string(x = time.var.plot, y = outcome.var, group = "XXclusterXX"), alpha = obs.alpha, size = obs.size[1])
            gg <- gg + ggplot2::geom_line(data = data, mapping = ggplot2::aes_string(x = time.var.plot, y = outcome.var, group = "XXclusterXX"), alpha = obs.alpha, size = obs.size[2])
        }
    }
    if(ci){
        if(is.na(ci.alpha)){
            gg <- gg + ggplot2::geom_errorbar(ggplot2::aes_string(ymin = "lower", ymax = "upper"), position = position.errorbar)
        }else{
            if(!is.null(color)){
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lower", ymax = "upper", fill = color), alpha = ci.alpha)
            }else{
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "lower", ymax = "upper"), alpha = ci.alpha)
            }
        }
    }
    if(!is.null(color)){
        gg <- gg + ggplot2::geom_point(ggplot2::aes_string(color = color), size = mean.size[1]) + ggplot2::geom_line(ggplot2::aes_string(color = color), size = mean.size[2])
    }else{
        gg <- gg + ggplot2::geom_point(size = mean.size[1]) + ggplot2::geom_line(size = mean.size[2])
    }
    gg  <- gg + ggplot2::ylab(outcome.var) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
    if(!is.null(time.var.plot) && any(!is.na(time.var.plot))){
        gg  <- gg + ggplot2::xlab(paste(stats::na.omit(xlabel.plot), collapse = ", "))
    }


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
