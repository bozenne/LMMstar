### autoplot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  8 2021 (00:01) 
## Version: 
## Last-Updated: jun 15 2023 (17:27) 
##           By: Brice Ozenne
##     Update #: 720
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
##' By default, normalized residuals are used except when requesting a partial residual plot.
##' @param at [data.frame] values for the covariates at which to evaluate the fitted values.
##' @param time.var [character] x-axis variable for the plot.
##' @param obs.alpha [numeric, 0-1] When not NA, transparency parameter used to display the original data by cluster.
##' @param obs.size [numeric vector of length 2] size of the point and line for the original data.
##' @param color [character] name of the variable in the dataset used to color the curve.
##' @param ci [logical] should confidence intervals be displayed?
##' @param ci.alpha [numeric, 0-1] When not NA, transparency parameter used to display the confidence intervals.
##' @param mean.size [numeric vector of length 2] size of the point and line for the mean trajectory.
##' @param size.text [numeric, >0] size of the font used to display text.
##' @param position.errorbar [character] relative position of the errorbars.
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
##' eCS.lmm <- lmm(Y ~ visit + X1,
##'                repetition = ~visit|id, structure = "CS", data = dL, df = FALSE)
##' 
##' plot(eCS.lmm, type = "fit")
##' autoplot(eCS.lmm, type = "fit")$plot + facet_wrap(~X1)
##' plot(eCS.lmm, type = "qqplot") ## engine.qqplot = "qqtest"
##' plot(eCS.lmm, type = "qqplot", engine.qqplot = "qqtest")
##' plot(eCS.lmm, type = "correlation") 
##' plot(eCS.lmm, type = "scatterplot") 
##' plot(eCS.lmm, type = "scatterplot2") 
##' plot(eCS.lmm, type = "partial", type.residual = "visit") 
##' plot(eCS.lmm, type = "partial", type.residual = "X1") 
##' }

## * autoplot.lmm (code)
##' @export
autoplot.lmm <- function(object, type = "fit", type.residual = "normalized", 
                         obs.alpha = 0, obs.size = c(2,0.5),
                         at = NULL, time.var = NULL, color = TRUE, ci = TRUE, ci.alpha = 0.25, 
                         ylim = NULL, mean.size = c(3, 1), size.text = 16, position.errorbar = "identity", ...){

    attr.ref <- attr(type,"reference")
    type <- match.arg(type, c("qqplot","correlation","scatterplot","scatterplot2","fit","partial"))

    if(type=="fit"){
        out <- .autofit(object,
                        obs.alpha = obs.alpha,
                        obs.size = obs.size,
                        at = at,
                        time.var = time.var,
                        color = color,                        
                        ci = ci,
                        ci.alpha = ci.alpha,
                        ylim = ylim,
                        mean.size = mean.size,
                        size.text = size.text,
                        position.errorbar = position.errorbar,
                        ...)
    }else if(type=="partial"){
        ## prepare
        if(is.null(match.call()$type.residual)){
            if(!is.null(list(...)$var)){
                type.residual <- list(...)$var
            }else{
                type.residual <- attr(object$design$mean,"variable")[1]
            }
        }        
        name.var <- setdiff(type.residual,"(Intercept)")
        if("(Intercept)" %in% type.residual){
            type.predict <- "static"
        }else{
            type.predict <- "static0"
        }
        type.var <- c("numeric","categorical")[name.var %in% names(object$xfactor$mean) + 1]
        if(sum(type.var=="numeric")>1){
            stop("Cannot simulatenously display partial residuals for more 2 numeric variables. \n")
        }

        ## extract partial residuals
        ttt <- "partial"
        attr(ttt,"reference") <- attr.ref
        rr <- stats::residuals(object, type = ttt, var = type.residual, keep.data = TRUE)
        ## extract predictions
        gg.data <- stats::predict(object, newdata = rr, keep.newdata = TRUE, type = type.predict)
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
            gg <- ggplot2::ggplot(data = gg.data, mapping = ggplot2::aes(x = .data[[name.varcat]]))
            gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial), color = "gray")
            if(ci){
                gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper))
            }
            gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$estimate), size = mean.size[1], shape = 21, fill = "white")
        }else{
            gg <- ggplot2::ggplot(data = gg.data, mapping = ggplot2::aes(x = .data[[name.varnum]]))
            gg <- gg + ggplot2::geom_point(ggplot2::aes(y = .data$r.partial), color = "gray")
            if(length(type.var)==1){
                gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$estimate), linewidth = mean.size[2])
                if(ci){
                    gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = ci.alpha)
                }
            }else{
                gg <- gg + ggplot2::geom_line(ggplot2::aes(y = .data$estimate,
                                                           group = .data[[name.varcat]],
                                                           color = .data[[name.varcat]]),
                                              linewidth = mean.size[2])
                if(ci){
                    gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                                                 ymax = .data$upper,
                                                                 group = .data[[name.varcat]],
                                                                 color = .data[[name.varcat]]),
                                                    alpha = ci.alpha)
                }
            }            
        }

        reference <- attr(rr,"reference")[,setdiff(names(attr(rr,"reference")),type.residual),drop=FALSE]
        reference <- lapply(reference, function(iRef){if(is.factor(iRef)){as.character(iRef)}else{iRef}})
        gg <- gg + ggplot2::ggtitle(paste0("Reference: ",paste(paste0(names(reference),"=",reference), collapse = ", ")))
        gg <- gg + ggplot2::ylab(paste0("Partial residuals for ",paste(name.var,collapse=", "))) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
        out <- list(data = gg.data,
                    plot = gg)
    }else{
        outRes <- residuals(object, type = type.residual, format = "long", keep.data = TRUE)
        out <- plot(outRes, type = type, size.text = size.text, ...)
    }

    ## ** export
    return(out)

}

## ** .autofit (helper to autofit.lmm)
.autofit <- function(object,
                     obs.alpha, obs.size,
                     at, time.var, color, ci, ci.alpha, 
                     ylim, mean.size, size.text, position.errorbar, ...){
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
    mu.var <- formula2var(object$formula$mean)$var$regressor
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
    timemu.var <- stats::na.omit(union(time.var, mu.var))
    X.beta <- stats::model.matrix(object, effects = "mean",
                                  data = data[,timemu.var, drop=FALSE])
    IX.beta <- interaction(as.data.frame(X.beta), drop = TRUE)
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
        preddata <- cbind(data, stats::predict(object, newdata = data[,timemu.var, drop=FALSE], ...))
    }else{
        preddata <- cbind(newdata, stats::predict(object, newdata = newdata[,timemu.var, drop=FALSE], ...))
    }

    ## ** generate plot
    gg <- ggplot2::ggplot(preddata, ggplot2::aes(x = .data[[time.var.plot]],
                                                 y = .data$estimate,
                                                 group = .data$XXclusterXX))
    if(!is.na(obs.alpha) && obs.alpha>0){
        if(!is.null(color)){
            gg <- gg + ggplot2::geom_point(data = data,
                                           mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                  y = .data[[outcome.var]],
                                                                  group = .data$XXclusterXX,
                                                                  color = .data[[color]]),
                                           alpha = obs.alpha,
                                           size = obs.size[1])
            gg <- gg + ggplot2::geom_line(data = data,
                                          mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                 y = .data[[outcome.var]],
                                                                 group = .data$XXclusterXX,
                                                                 color = .data[[color]]),
                                          alpha = obs.alpha,
                                          linewidth = obs.size[2])
            ## gg + facet_wrap(~XXclusterXX)
        }else{
            gg <- gg + ggplot2::geom_point(data = data,
                                           mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                  y = .data[[outcome.var]],
                                                                  group = .data$XXclusterXX),
                                           alpha = obs.alpha,
                                           size = obs.size[1])
            gg <- gg + ggplot2::geom_line(data = data,
                                          mapping = ggplot2::aes(x = .data[[time.var.plot]],
                                                                 y = .data[[outcome.var]],
                                                                 group = .data$XXclusterXX),
                                          alpha = obs.alpha,
                                          linewidth = obs.size[2])
        }
    }
    if(ci){
        if(is.na(ci.alpha)){
            gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), position = position.errorbar)
        }else{
            if(!is.null(color)){
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data[[color]]), alpha = ci.alpha)
            }else{
                gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = ci.alpha)
            }
        }
    }
    if(!is.null(color)){
        gg <- gg + ggplot2::geom_point(ggplot2::aes(color = .data[[color]]), size = mean.size[1]) + ggplot2::geom_line(ggplot2::aes(color = .data[[color]]), linewidth = mean.size[2])
    }else{
        gg <- gg + ggplot2::geom_point(size = mean.size[1]) + ggplot2::geom_line(linewidth = mean.size[2])
    }
    gg  <- gg + ggplot2::ylab(outcome.var) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
    if(!is.null(time.var.plot) && any(!is.na(time.var.plot))){
        gg  <- gg + ggplot2::xlab(paste(stats::na.omit(xlabel.plot), collapse = ", "))
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
##' @param by.repetition [logical] Should a seperate graphical display be made for each repetition.
##' @param engine.qqplot [character] Should ggplot2 or qqtest be used to display quantile-quantile plots?
##' Only used when argument \code{type} is \code{"qqplot"}.
##' @param add.smooth [logical] should a local smoother be used to display the mean of the residual values across the fitted values.
##' Only relevant for when argument \code{type} is \code{"scatterplot"}.
##' @param digits.cor [integer, >0] Number of digit used to display the correlation coefficients?
##' No correlation coefficient is displayed when set to 0. Only used when argument \code{plot} is \code{"correlation"}.
##' @param size.text [numeric, >0] Size of the font used to displayed text when using ggplot2.
##' @param scales,labeller [character] Passed to \code{ggplot2::facet_wrap}.
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
autoplot.residuals_lmm <- function(object, type = NULL, type.residual = NULL, by.repetition = TRUE, 
                                   engine.qqplot = "ggplot2", add.smooth = TRUE, digits.cor = 2, size.text = 16,
                                   scales = "free", labeller = "label_value",...){


    ## ** check arguments
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    args <- attr(object,"args")
    args.type <- args$type
    n.type <- args$n.type
    n.time <- args$n.time
    name.time <- args$name.time
    format <- args$format
    if(n.time == 1){
        by.repetition <- FALSE
    }
    
    if(is.null(type.residual)){
        if(length(args.type)==1){
            type.residual <- args.type
        }else{
            stop("Different types of residuals are available: \"",paste(args.type, collapse = "\", \""),"\"\n",
                 "Select one type via the argument \'type.residual\'. \n")
        }
    }else if(length(args.type)==1){
        if(type.residual %in% args$type == FALSE){
            stop("Requested type of residual not available. \n",
                 "Available type(s) of residuals: \"",paste(args.type, collapse = "\", \""),"\"\n")
        }
    }else{
        stop("Can only display one type of residual. \n")
    }
    
    if(is.null(type)){
        if(type.residual %in% c("partial","partial-center")){
            type <- "scatterplot"
        }else{
            type <- "qqplot"
        }
    }else{
        type <- match.arg(type, c("qqplot","correlation","scatterplot","scatterplot2"))
    }
    if(length(add.smooth)==1){
        add.smooth <- rep(add.smooth,2)
    }

    label.residual <- switch(type.residual,
                             "fitted" = "Fitted values",
                             "response" = "Raw residuals",
                             "studentized" = "Studentized residuals",
                             "pearson" = "Pearson residuals",
                             "normalized" = "Normalized residuals",
                             "normalized2" = "Pearson Normalized residuals",
                             "scaled" = "Scaled residuals",
                             "partial" = paste("Partial residuals for ",paste(args$var,collapse=", ")),
                             "partial" = paste("Partial residuals for ",paste(args$var,collapse=", ")))
    
    if(format == "long"){
        name.residual <- args$name.colres[which(type.residual == args.type)]
    }else{
        name.residual <- args$name.colres
    }

    if(type == "correlation"){
        if(n.time == 1){
            stop("Cannot display the residual correlation over time when there is only a single timepoint. \n")
        }
    }else if(args$format == "wide"){
        stop("Residuals must be in the long format to obtain a",type,". \n",
             "Consider setting the argument \'format\' to \"long\" when calling residuals(). \n")
    }
    if(type == "correlation" || (type == "qqplot" && engine.qqplot == "qqtest" && by.repetition)){
        if(format == "long"){
            objectW <- stats::reshape(data = object[,c("XXclusterXX","XXtimeXX",name.residual),drop=FALSE], 
                                      direction = "wide", timevar = "XXtimeXX", idvar = "XXclusterXX", v.names = name.residual,
                                      sep = ".")
            name.residual <- colnames(objectW)[-1]
        }else{
            objectW <- object
        }
    }
    
    if(type %in% c("scatterplot","scatterplot2") && args$keep.data == FALSE){
        stop("Cannot display a scatterplot of the residuals without the original data/fitted values. \n",
             "Consider setting argument \'keep.data\' to TRUE when calling residuals(). \n")
    }
    if(by.repetition>0 && args$keep.data == FALSE){
        stop("Cannot display a scatterplot of the residuals per repetition without the original data/fitted values. \n",
             "Consider setting argument \'keep.data\' to TRUE when calling residuals(). \n")
    }

    ## ** build graphical display
    if(by.repetition>0 && args$keep.data){
        if(length(name.time) == 1 && (name.time %in% names(object) == FALSE)){
            object[[name.time]] <- object[["XXtimeXX"]]
        }
        formula.time <- stats::as.formula(paste("~",paste(name.time, collapse = "+")))
    }

    if(type=="qqplot"){ ## overall timepoints

        if(args$keep.data == FALSE && format == "long"){
            df.gg <- as.data.frame(object)
            names(df.gg) <- name.residual
        }else{
            df.gg <- object
        }

        if(engine.qqplot=="ggplot2"){
            gg <- ggplot2::ggplot(df.gg, ggplot2::aes(sample = .data[[name.residual]]))
            gg <- gg + ggplot2::stat_qq() + ggplot2::stat_qq_line()
            gg <- gg + ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles")
            gg <- gg + ggplot2::ggtitle(label.residual) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
            if(by.repetition>0){
                gg <- gg + ggplot2::facet_wrap(formula.time, scales = scales, labeller = labeller)
            }
        }else if(engine.qqplot=="qqtest"){
            requireNamespace("qqtest")

            if(by.repetition){
                sqrt.round <- ceiling(sqrt(n.time))
                sqrt.round2 <- ceiling(n.time/sqrt.round)
                Utime <- unique(object$XXtimeXX)

                oldpar <- graphics::par(no.readonly = TRUE)   
                on.exit(graphics::par(oldpar))            
                graphics::par(mfrow = c(sqrt.round,sqrt.round2))                
                lapply(1:n.time,function(iCol){
                    qqtest::qqtest(stats::na.omit(objectW[,iCol+1]), main = Utime[iCol])
                    graphics::mtext(label.residual, side = 3)
                })
                
            }else{
                qqtest::qqtest(stats::na.omit(df.gg[[name.residual]]), main = label.residual)
            }

            df.gg <- NULL
            gg <- NULL
        }
    }else if(type %in% c("scatterplot","scatterplot2")){ ## overall timepoints
        if(type.residual %in% c("partial","partial-center")){
            name.fitted <- args$var
            xlab.gg <- args$var
        }else{
            name.fitted <- "fitted"
            xlab.gg <- "Fitted values"
        }
        gg <- ggplot2::ggplot(object) + ggplot2::xlab(xlab.gg) + ggplot2::theme(text = ggplot2::element_text(size=size.text))
        if(type == "scatterplot"){
            gg <- gg + ggplot2::geom_abline(slope=0,intercept=0,color ="red")
            gg <- gg + ggplot2::geom_point(ggplot2::aes(x = .data[[name.fitted]], y = .data[[name.residual]]))
            gg <- gg + ggplot2::ylab(label.residual) 
            if(add.smooth[1]){
                gg <- gg + ggplot2::geom_smooth(ggplot2::aes(x = .data[[name.fitted]], y = .data[[name.residual]]), se = add.smooth[2])
            }
        }else if(type == "scatterplot2"){
            label.residual2 <- paste0("|",label.residual,"|")
            gg <- gg + ggplot2::geom_point(ggplot2::aes(x = .data[[name.fitted]], y = sqrt(abs(.data[[name.residual]]))))
            gg <- gg + ggplot2::ylab(bquote(sqrt(.(label.residual2))))
            if(add.smooth[1]){
                gg <- gg + ggplot2::geom_smooth(ggplot2::aes(x = .data[[name.fitted]], y = sqrt(abs(.data[[name.residual]]))), se = add.smooth[2])
            }
        }
        if(by.repetition>0){
            gg <- gg + ggplot2::facet_wrap(formula.time, scales = scales, labeller = labeller)
        }
        df.gg <- NULL
    }else if(type == "correlation"){ 
        M.cor  <- stats::cor(objectW[,-1], use = "pairwise")
        ind.cor <- !is.na(M.cor)
        arr.ind.cor <- which(ind.cor, arr.ind = TRUE)
        arr.ind.cor[] <- name.residual[arr.ind.cor]

        df.gg <- data.frame(correlation = M.cor[ind.cor], arr.ind.cor, stringsAsFactors = FALSE)
        name.time <- gsub(paste0("^",args$name.colres,"."),"",name.residual, fixed = FALSE)
        df.gg$col <- factor(df.gg$col, levels = name.residual, labels = name.time)
        df.gg$row <- factor(df.gg$row, levels = name.residual, labels = name.time)
        df.gg$row.index <- match(df.gg$row, name.time)
        df.gg$col.index <- match(df.gg$col, name.time)
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
        gg <- gg + ggplot2::labs(x = name.time, y = name.time) + ggplot2::ggtitle(label.residual)
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
##' @param object an object of class \code{summarize}, output of the \code{summarize} function.
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
            name.Group <- as.character(interaction(name.group))
            object[[name.Group]] <- interaction(object[,name.group])
            gg <- ggplot2::ggplot(object, ggplot2::aes(x = .data[[name.time]], y = .data[[type]],
                                                       group = .data[[name.Group]], color = .data[[name.Group]]))
        }
        gg <- gg + geom_point(size = size) + geom_line(linewidth = linewidth)
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

    keep.cols <- setdiff(names(object),newnames)
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
        gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), size = size.ci, width = width.ci)
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
            if(all(sapply(splitname,length)==2) && length(unique(sapply(splitname,"[",2)))==1){
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
