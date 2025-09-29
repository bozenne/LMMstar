### summarize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: sep 29 2025 (14:51) 
##           By: Brice Ozenne
##     Update #: 792
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summarize (documentation)
##' @title Summary Statistics
##' @description Compute summary statistics for multiple variables and/or multiple groups and save them in a data frame.
##'
##' @param formula [formula] on the left hand side the outcome(s) and on the right hand side the grouping variables.
##' E.g. Y1+Y2 ~ Gender + Gene will compute for each gender and gene the summary statistics for Y1 and for Y2.
##' @param data [data.frame] dataset containing the observations.
##' @param repetition [formula] Specify the structure of the data: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' Used in the long format to count the number of missing values (i.e. add the number of missing rows) and evaluate the correlation.
##' @param na.action [function] a function which indicates what should happen when the data contain 'NA' values.
##' Passed to the \code{stats::aggregate} function. 
##' @param na.rm [logical] Should the summary statistics be computed by omitting the missing values.
##' @param columns [character vector] name of the summary statistics to kept in the output.
##' Can be any of, or a combination of:\itemize{
##' \item \code{"observed"}: number of observations with a measurement.
##' \item \code{"missing"}: number of missing observations.
##' When specifying a grouping variable, it will also attempt to count missing rows in the dataset.
##' \item \code{"pc.missing"}: percentage missing observations.
##' \item \code{"mean"}, \code{"mean.lower"} \code{"mean.upper"}: mean with its confidence interval.
##' \item \code{"median"}, \code{"median.lower"} \code{"median.upper"}: median with its confidence interval.
##' \item \code{"sd"}, \code{"sd.lower"}, \code{"sd.upper"}: standard deviation around the mean with its confidence interval.
##' \item \code{"sd0"}, \code{"sd0.lower"}, \code{"sd0.upper"}: standard deviation around 0 with its confidence interval.
##' \item \code{"skewness"}: skewness, as the third standardized moment.
##' \item \code{"kurtosis"}: kurtosis, as the fourth standardized moment.
##' \item \code{"q1"}, \code{"q3"}, \code{"IQR"}: 1st and 3rd quartile, interquartile range.
##' \item \code{"min"}, \code{"max"}: minimum and maximum observation.
##' \item \code{"predict.lower"}, \code{"predict.upper"}: prediction interval for normally distributed outcome.
##' \item \code{"correlation"}: correlation matrix between the outcomes (when feasible, see detail section).
##' }
##' @param FUN [function] user-defined function for computing summary statistics.
##' It should take a vector as an argument and output a named single value or a named vector.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param skip.reference [logical] should the summary statistics for the reference level of categorical variables be omitted?
##' @param digits [integer, >=0] the minimum number of significant digits to be used to display the results. Passed to \code{print.data.frame}
##' @param filter [character] a regular expression passed to \code{grep} to filter the columns of the dataset.
##' Relevant when using \code{.} to indicate all other variables.
##' @param ... additional arguments passed to argument \code{FUN}.
##'
##' @details This function is essentially an interface to the \code{stats::aggregate} function. \cr
##' \bold{WARNING:} it has the same name as a function from the dplyr package. If you have loaded dplyr already, you should use \code{:::} to call summarize i.e. use \code{LMMstar:::summarize}.
##' 
##' Confidence intervals (CI) and prediction intervals (PI) for the mean are computed via \code{stats::t.test}.
##' Confidence intervals (CI) for the standard deviation are computed using a chi-squared approximation.
##' Confidence intervals (CI) for the median are computed via \code{asht::medianTest}.
##' 
##' @return A data frame containing summary statistics (in columns) for each outcome and value of the grouping variables (rows). It has an attribute \code{"correlation"} when it was possible to compute the correlation matrix for each outcome with respect to the grouping variable.
##' 
##' @seealso
##' \code{\link{correlate}} for correlation matrix.
##' 
##' @keywords utilities
 
## * summarize (examples)
##' @examples
##' #### simulate data (wide format) ####
##' set.seed(10)
##' d <- sampleRem(1e2, n.times = 3)
##' d$treat <-  sample(LETTERS[1:3], NROW(d), replace=TRUE, prob=c(0.3, 0.3, 0.4) )
##'
##' ## add a missing value
##' d2 <- d
##' d2[1,"Y2"] <- NA
##'
##' #### summarize (wide format) ####
##' 
##' ## summary statistic (single variable)
##' summarize(Y1 ~ 1, data = d)
##' ## stratified summary statistic (single variable)
##' summarize(Y1 ~ X1, data = d2)
##' ## stratified summary statistic (multiple variable)
##' summarize(Y1+Y2 ~ X1, data = d)
##' ## categorical variable
##' summarize(treat ~ 1, data = d)
##' summarize(treat ~ 1, skip.reference = TRUE, data = d)
##' ## aggregate data
##' summarize( ~ X1 + treat, data = d)
##' ## user defined summary statistic
##' summarize(Y1 ~ 1, data = d, FUN = quantile)
##' summarize(Y1 ~ 1, data = d, FUN = quantile, p = c(0.25,0.75))
##' ## complete case summary statistic
##' summarize(Y1+Y2 ~ X1, data = d2, na.rm = TRUE)
##' ## shortcut to consider all outcomes with common naming 
##' summarize(. ~ treat, data = d2, na.rm = TRUE, filter = "Y")
##' 
##' #### summarize (long format) ####
##' dL <- reshape(d2, idvar = "id", direction = "long",
##'              v.names = "Y", varying = c("Y1","Y2","Y3"))
##' summarize(Y ~ time + X1, data = dL, na.rm  = TRUE)
##' 
##' ## user defined summary statistic (outlier)
##' summarize(Y ~ time + X1, data = dL, FUN = function(x){
##'    c(outlier.down = sum(x<mean(x,na.rm=TRUE)-2*sd(x,na.rm=TRUE), na.rm=TRUE),
##'      outlier.up = sum(x>mean(x,na.rm=TRUE)+2*sd(x,na.rm=TRUE), na.rm=TRUE))
##' }, na.rm = TRUE)
##'
##' ## user defined summary statistic (auc)
##' myAUC <- function(Y,time){approxAUC(x = time, y = Y, from = 1, to = 3)}
##' myAUC(Y = dL[dL$id==1,"Y"], time = dL[dL$id==1,"time"])
##' summarize(Y ~ id, data = dL, FUN = myAUC, na.rm = TRUE)
##'
##' ## add correlation (see correlate function)
##' e.S <- summarize(Y ~ time + X1, data = dL, repetition = ~time|id,
##'                  na.rm = TRUE, columns = add("correlation"))
##' e.S
##' 
##' #### summarize (long format, missing lines) ####
##' ## use repetition argument to count missing lines in the number of missing values
##' dL.NNA <- dL[rowSums(is.na(dL))==0,]
##' summarize(Y ~ time + X1, data = dL.NNA, repetition =~time|id, na.rm  = TRUE)

## * summarize (code)
##' @export
summarize <- function(formula, data, repetition = NULL, columns = NULL, FUN = NULL,
                      na.action = stats::na.pass, na.rm = FALSE,
                      level = 0.95,
                      skip.reference = FALSE,
                      digits = NULL,
                      filter = NULL,
                      ...){

    data <- as.data.frame(data)
    mycall <- match.call()
    options <- LMMstar.options()

    ## ** check and normalize user imput
    ## *** dots
    dots <- list(...)
    if(length(dots)>0 && is.null(FUN)){
        stop("Unknown arguments \'",paste(names(dots), collapse = "\' \'"),"\'. \n")
    }

    ## *** columns
    valid.columns <- c("observed","missing","pc.missing",
                       "mean","mean.lower","mean.upper","predict.lower","predict.upper",
                       "sd","sd.lower","sd.upper", "sd0","sd0.lower","sd0.upper",
                       "skewness",
                       "kurtosis",
                       "min","q1","median","q3","median.upper","median.lower","IQR","max",
                       "correlation")
    default.columns <- options$columns.summarize
    if(identical(columns,"all")){
        newcolumns <- valid.columns
    }else if(is.null(columns)){
        newcolumns <- default.columns
    }else{
        if(!is.null(names(columns)) && all(names(columns)=="add")){
            newcolumns <- union(columns, unname(default.columns))
        }else if(!is.null(names(columns)) && all(names(columns)=="remove")){
            newcolumns <- setdiff(columns, unname(default.columns))
        }else{
            newcolumns <- columns
        }
    }
    if(!is.character(newcolumns)){
        stop("Argument \'columns\' should be a character vector. \n")
    }
    if(any(newcolumns %in% valid.columns == FALSE) ){
        stop("Argument \'columns\' cannot take value(s) \"",paste(setdiff(newcolumns, valid.columns), collapse = "\", \""),"\". \n",
             "Possible values: \"",paste(setdiff(valid.columns, newcolumns), collapse = "\", \""),"\". \n")
    }else{
        columns <- newcolumns
    }

    ## *** data (column names)
    name.all <- names(data)
    if(any(name.all %in% c("XXn.clusterXX","XXindexXX"))){
        invalid <- name.all[name.all %in% c("XXn.clusterXX","XXindexXX")]
        stop("Name(s) \"",paste(invalid, collapse = "\" \""),"\" are used internally. \n",
             "Consider renaming the variables in the dataset and updating the formula. \n",
             sep = "")
    }
    if(any(name.all %in% names(data) == FALSE)){
        invalid <- name.all[name.all %in% names(data) == FALSE]
        stop("Argument \'formula\' is inconsistent with argument \'data\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }

    ## *** formula & repetition
    ls.init <- formula2repetition(formula = formula, data = data, repetition = repetition, keep.time = TRUE, filter = filter)
    detail.formula <- ls.init$detail.formula
    detail.repetition <- ls.init$detail.repetition
    name.Y <- detail.formula$var$response
    name.X <- detail.formula$var$regressor
    name.cluster <- detail.repetition$var$cluster
    name.time <- detail.repetition$var$time
    
    n.Y <- length(name.Y)
    if(n.Y == 0){
        n.Y <- 1
        columns <- "observed"
    }
    n.X <- length(name.X)

    ## *** FUN
    if(!is.null(FUN)){

        formals.fun <- formals(FUN)
        names.FUN <- setdiff(names(formals.fun[sapply(formals.fun,is.symbol)]), c("...",names(dots)))
        
        if(length(names.FUN)==0){

            try.FUN <- try(FUN(), silent = FALSE)
            Mnames.FUN <- NULL

        }else{
            ## match user defined arguments to available variables

            if(any(names.FUN %in% names(data) == FALSE)){
                possible.arg <- c(name.Y[1], name.time, name.X, name.cluster)
                if(length(possible.arg) == 0){
                    stop("Unknown argument(s) \"",paste(setdiff(names.FUN, names(data)), collapse = "\", \""),"\" in argument \'FUN\'. \n")
                }else if(sum(names.FUN %in% names(data) == FALSE) > length(possible.arg)){
                    stop("Unknown argument(s) \"",paste(setdiff(names.FUN, names(data)), collapse = "\", \""),"\" in argument \'FUN\'. \n",
                         "Can use variables \"",paste(possible.arg, collapse = "\", \""),"\" but this is not enough variables. \n")
                }
                Mnames.FUN <- do.call(rbind,lapply(1:n.Y, FUN = function(iY){
                    iPossible.arg <- c(name.Y[iY], name.time, name.X, name.cluster)
                    iNames <- names.FUN
                    iNames[iNames %in% names(data) == FALSE] <- iPossible.arg[1:sum(names.FUN %in% names(data) == FALSE)]
                    return(iNames)
                }))
                
            }else{
                Mnames.FUN <- matrix(names.FUN, byrow = TRUE, nrow = n.Y, ncol = length(names.FUN))
            }
            colnames(Mnames.FUN) <- names.FUN
            if(!is.null(name.Y)){
                rownames(Mnames.FUN) <- name.Y
            }
            if(is.null(name.X)){
                try.FUN <- try(do.call(FUN, args = c(stats::setNames(as.list(data[Mnames.FUN[1,]]),names.FUN),dots)), silent = FALSE)
            }else{
                strata.X <- interaction(data[name.X],drop=TRUE)
                try.FUN <- try(do.call(FUN, args = c(stats::setNames(as.list(data[strata.X==levels(strata.X)[1],Mnames.FUN[1,],drop=FALSE]),names.FUN),dots)), silent = FALSE)
            }
            
        }

        if(inherits(try.FUN,"try-error")){
            return(invisible(try.FUN))
        }else if(!inherits(try.FUN,"try-error")){
            if(is.null(names(try.FUN))){
                stop("Argument \'FUN\' should return values with names. \n")
            }else if(any(names(try.FUN) %in% valid.columns)){
                stop("Argument \'FUN\' should not return values named \"",paste(names(try.FUN)[names(try.FUN) %in% valid.columns], collapse = "\" \""),"\" are used internally. \n")
            }
        }

        name.outFUN <- names(try.FUN)

    }else{
        Mnames.FUN <- NULL
        name.outFUN <- NULL
    }

    ## ** handle categorical variables
    to.rm <- NULL
    to.add <- NULL

    for(iY in 1:n.Y){ ## iY <- 3
        ## Note: name.Y can be NULL when no left hand side of formula argument (n.Y-->1)
        if(!is.null(name.Y) && (is.character(data[[name.Y[iY]]])||is.factor(data[[name.Y[iY]]]))){
            data[[name.Y[iY]]] <- as.factor(data[[name.Y[iY]]])
            iLevel <- levels(data[[name.Y[iY]]])
            if(any(paste(name.Y[iY],iLevel,sep=":") %in% names(data))){
                stop("Name(s) \"",paste(paste(name.Y[iY],iLevel,sep=":")[paste(name.Y[iY],iLevel,sep=":") %in% names(data)], collapse ="\" \""),"\" are being used internally. \n",
                     "Consider rename or dropping columns from argument \'data\'.\n")
            }
            
            for(iL in iLevel){
                if(!skip.reference || iL != iLevel[1]){
                    data[[paste(name.Y[iY],iL,sep=":")]] <- as.numeric(data[[name.Y[iY]]]==iL)
                    to.add <- c(to.add, paste(name.Y[iY],iL,sep=":"))
                }
            }
            to.rm <- c(to.rm, name.Y[iY])
        }        
    }
    if(length(to.add)!=0 || length(to.rm)!=0){
        name.Y <- unique(c(setdiff(name.Y,to.rm),to.add))
        n.Y <- length(name.Y)
    }

    ## ** prepare 
    ## WARNING: the time variable is transformed as factor in iData but kept un-touched in data
    out <- NULL
    level.int <- c((1-level)/2,1-(1-level)/2) ## confidence level (typically 0.025 and 0.975)

    ## duplicate dataset with index and factor time
    formula.index <- stats::update(ls.init$formula, paste0("XXindexXX~."))
    data.index <- cbind(XXindexXX = 1:NROW(data), data)
    if(!is.null(name.time) && !is.null(name.cluster) && ("missing" %in% columns || "pc.missing" %in% columns)){
        if(length(setdiff(name.X,name.time))==0){            
            level.UX <- rep(1, NROW(data))
        }else{
            level.UX <- interaction(data[setdiff(name.X,name.time)])
        }
        data.index$XXn.clusterXX <- tapply(data[[name.cluster]],level.UX, FUN = function(iVec){length(unique(iVec))})[level.UX]
        for(iTime in name.time){ ## necessary to find out missing rows 
            data.index[[iTime]] <- as.factor(data.index[[iTime]])
        }
    }

    ## all output names
    name.out <- setdiff(c(columns, name.outFUN), "correlation")

    ## type.Y
    if(!is.null(name.Y)){
        type.Y <- stats::setNames(c("continuous","binary")[sapply(name.Y, function(iName){
            all(data[[iName]] %in% 0:1)
        })+1], name.Y)
    }

    ## split according to formula
    if(is.null(name.X)){
        ls.indexSplit <- list(XXindexXX = list(1:NROW(data)))
        n.grid <- 1
    }else{
        ls.indexSplit <- stats::aggregate(formula.index, data=data.index, FUN = identity,  na.action=na.action, simplify = FALSE)
        grid <- as.data.frame(ls.indexSplit[name.X])
        n.grid <- NROW(grid)
    }

    ## ** compute summary statistics
    ls.summary <- lapply(1:n.grid, function(iG){

        iIndex <- ls.indexSplit$XXindexXX[[iG]]
        iOut <- as.data.frame(matrix(NA, nrow = n.Y, ncol = length(name.out), dimnames = list(NULL, name.out)))
        if(!is.null(name.X)){
            iOut <- cbind(grid[rep(iG,n.Y),,drop=FALSE], iOut)
        }
        
        ## *** observed
        if(is.null(name.Y)){ ## no outcome variable

            if("observed" %in% columns){
                iOut$observed <- length(iIndex)
            }

        }else{

            iOut <- cbind(outcome = name.Y, iOut)
            iDataY <- data[iIndex,name.Y,drop=FALSE]
            iDataYcont <- iDataY[type.Y=="continuous"]

            iN.obs <- colSums(!is.na(iDataY))
            if("observed" %in% columns){
                iOut$observed <- iN.obs
            }

        }

        ## *** missing
        if("missing" %in% columns || "pc.missing" %in% columns){
            iN.missing <- colSums(is.na(iDataY))

            ## add missing rows in the dataset
            if(!is.null(name.time) && !is.null(name.cluster)){
                iData.index <- data.index[iIndex,c("XXn.clusterXX",name.cluster,name.time),drop=FALSE]
                iN.cluster <- unique(iData.index$XXn.clusterXX)
                if(length(iN.cluster)>1){
                    warning("Something went wrong when identifying the number of missing values. \n")
                }
                if(any(name.time %in% name.X)){
                    for(iTime in intersect(name.time,name.X)){
                        iData.index[[iTime]] <- droplevels(iData.index[[iTime]])
                    }
                }                
                iTable.time <- table(interaction(iData.index[name.time]))
                iN.missing <- iN.missing + sum(iN.cluster[1] - iTable.time)
            }

            if("missing" %in% columns){
                iOut$missing <- iN.missing
            }
            if("pc.missing" %in% columns){
                iOut$pc.missing <- iN.missing/(iN.obs+iN.missing)
            }
        }

        ## *** mean
        if("mean" %in% columns && any(!is.na(iDataY))){ ## otherwise colMeans outputs NaN instead of NA
            iOut$mean <- colMeans(iDataY, na.rm = na.rm)
        }
        
        if("mean.lower" %in% columns || "mean.upper" %in% columns || "predict.lower" %in% columns || "predict.upper" %in% columns){
            iMeanTest <- do.call(rbind,lapply(name.Y, function(iY){
                iTest <- switch(type.Y[iY],
                                "continuous" = stats::t.test(iDataY[[iY]], na.rm = na.rm, conf.level = level, alternative = "two.sided"),
                                "binary" = stats::binom.test(x = sum(iDataY[[iY]]==1, na.rm = TRUE), n = sum(!is.na(iDataY[[iY]])), conf.level = level, alternative = "two.sided"))
                if(type.Y == "continuous" && ("predict.lower" %in% columns || "predict.upper" %in% columns)){
                    iTest$pred.int <- c(iTest$estimate + sqrt(iN.obs+1) * stats::qt(level.int, iTest$parameter) * iTest$stderr[1])
                }
                c(estimate = unname(iTest$estimate), mean.lower = iTest$conf.int[1], mean.upper = iTest$conf.int[2], predict.lower = iTest$pred.int[1], predict.upper = iTest$pred.int[2])
            }))

            iOut[,intersect(columns,colnames(iMeanTest))] <- iMeanTest[,intersect(columns,colnames(iMeanTest))]
        }

        ## *** standard deviation
        if(any(c("sd","sd.lower","sd.upper") %in% columns) && any(type.Y=="continuous")){
            iSD <- apply(iDataYcont, MARGIN = 2, FUN = stats::sd, na.rm = na.rm)
            if("sd" %in% columns){
                iOut$sd[type.Y=="continuous"] <- iSD
            }
            if("sd.lower" %in% columns){
                iOut$sd.lower[type.Y=="continuous"] <- iSD * sqrt((iN.obs-1)/stats::qchisq(level.int[2], iN.obs-1))
            }
            if("sd.upper" %in% columns){
                iOut$sd.upper[type.Y=="continuous"] <- iSD * sqrt((iN.obs-1)/stats::qchisq(level.int[1], iN.obs-1))
            }
        }

        ## *** standard deviation without centering
        if(any(c("sd0","sd0.lower","sd0.upper") %in% columns) && any(type.Y=="continuous")){
            iSD0 <- apply(iDataYcont, MARGIN = 2, FUN = function(x){sqrt(mean(x^2, na.rm = na.rm))})
            if("sd0" %in% columns){
                iOut$sd[type.Y=="continuous"] <- iSD0
            }
            if("sd0.lower" %in% columns){
                iOut$sd0.lower[type.Y=="continuous"] <- iSD0 * sqrt(iN.obs/stats::qchisq(level.int[2], iN.obs))
            }
            if("sd0.upper" %in% columns){
                iOut$sd0.upper[type.Y=="continuous"] <- iSD0 * sqrt(iN.obs/stats::qchisq(level.int[1], iN.obs))
            }
        }

        ## *** higher-order moments
        if("skewness" %in% columns && any(type.Y=="continuous")){
            iOut$skewness[type.Y=="continuous"] <- colMeans(scale(iDataYcont)^3, na.rm = na.rm)
        }
        if("kurtosis" %in% columns && any(type.Y=="continuous")){
            iOut$kurtosis[type.Y=="continuous"] <- colMeans(scale(iDataYcont)^4, na.rm = na.rm)
        }

        ## *** median
        if("median" %in% columns && any(type.Y=="continuous")){
            iOut$median[type.Y=="continuous"] <- apply(iDataYcont, MARGIN = 2, FUN = stats::median, na.rm = na.rm)
        }
        if(("median.lower" %in% columns || "median.upper" %in% columns) && requireNamespace("asht") && any(type.Y=="continuous")){
            if(na.rm){
                iLs.MedianTest <- apply(iDataYcont, MARGIN = 2, FUN = function(x){asht::medianTest(stats::na.omit(x), conf.level = level, alternative = "two.sided")$conf.int},
                                        simplify = FALSE)
            }else{
                iLs.MedianTest <- apply(iDataYcont, MARGIN = 2, FUN = function(x){asht::medianTest(x, conf.level = level, alternative = "two.sided")$conf.int},
                                        simplify = FALSE)
            }
            iMedianTest <- do.call(rbind,iLs.MedianTest)
            colnames(iMedianTest) <- c("median.lower","median.upper")
            iOut[,intersect(columns,colnames(iMedianTest))] <- iMedianTest[,intersect(columns,colnames(iMedianTest))]
            
        }

        ## *** other quantiles
        if("min" %in% columns){
            iOut$min <- apply(iDataY, MARGIN = 2, FUN = function(iVec){
                if(all(is.na(iVec))){NA}else{min(iVec, na.rm = na.rm)} ## otherwise min output a warning
            })
        }
        if("q1" %in% columns && any(type.Y=="continuous")){
            iOut$q1[type.Y=="continuous"] <- apply(iDataYcont, MARGIN = 2, FUN = stats::quantile, prob = 0.25, na.rm = na.rm)
        }
        if("q3" %in% columns && any(type.Y=="continuous")){
            iOut$q3[type.Y=="continuous"] <- apply(iDataYcont, MARGIN = 2, FUN = stats::quantile, prob = 0.75, na.rm = na.rm)
        }
        if("max" %in% columns && any(!is.na(iDataY))){
            iOut$max <- apply(iDataY, MARGIN = 2, FUN = function(iVec){
                if(all(is.na(iVec))){NA}else{max(iVec, na.rm = na.rm)} ## otherwise min output a warning
            })
        }

        ## *** Interquartile range
        if("IQR" %in% columns && any(type.Y=="continuous")){
            iOut$IQR[type.Y=="continuous"] <- apply(iDataYcont, MARGIN = 2, FUN = stats::IQR, na.rm = na.rm)
        }

        ## *** User defined function
        if(!is.null(FUN)){
            if(length(names.FUN)==0){
                iOut[,name.outFUN] <- matrix(FUN(), byrow = TRUE, nrow = n.Y, ncol = length(name.outFUN))
            }else{
                iLS.outfun <- lapply(1:n.Y, FUN = function(iY){
                    do.call(FUN, args = c(stats::setNames(as.list(data[iIndex,Mnames.FUN[iY,],drop=FALSE]),names.FUN),dots))
                })
                iOut[,name.outFUN] <- do.call(rbind,iLS.outfun)
            }                
            
        }

        ## *** export
        return(iOut)
    })

    out <- do.call(rbind,ls.summary)
    if(!is.null(name.Y)){
        out <- out[order(out$outcome),]
    }
    rownames(out) <- NULL

    ## ** correlation
    if("correlation" %in% columns){
        if(all(ls.init$detail.formula$var$regressor %in% ls.init$detail.repetition$var$time)){
            ls.init$formula <- stats::reformulate("1", response = paste(name.Y,collapse="+"))
        }else if(any(ls.init$detail.formula$var$regressor %in% ls.init$detail.repetition$var$time)){
            terms.rm <- intersect(ls.init$detail.formula$var$regressor, ls.init$detail.repetition$var$time)
            formula.terms <- stats::terms(ls.init$formula)            
            ls.init$formula <- stats::drop.terms(formula.terms, dropx = which(attr(formula.terms,"term.labels") %in% terms.rm), keep.response = TRUE)
        }
        use <- ifelse(na.rm, "pairwise.complete.obs", "everything")
        attr(out,"correlation") <- correlate(formula = ls.init$formula, data = data, repetition = ls.init$repetition, use = use)
        
    }

    ## ** export
    if(!is.null(digits)){
        attr(out,"digits") <- digits
    }
    attr(out,"call") <- mycall
    attr(out,"Mnames.FUN") <- Mnames.FUN
    attr(out,"name.Y") <- name.Y
    attr(out,"name.X") <- name.X
    attr(out,"name.time") <- name.time
    attr(out,"name.cluster") <- name.cluster
    class(out) <- append("summarize",class(out))
    return(out)
}


######################################################################
### summarize.R ends here
