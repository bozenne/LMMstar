### summarize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: jul  1 2022 (16:12) 
##           By: Brice Ozenne
##     Update #: 220
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summarize (documentation)
##' @title Compute summary statistics
##' @description Compute summary statistics for multiple variables and/or multiple groups and save them in a data frame.
##'
##' @param formula [formula] on the left hand side the outcome(s) and on the right hand side the grouping variables.
##' E.g. Y1+Y2 ~ Gender + Gene will compute for each gender and gene the summary statistics for Y1 and for Y2.
##' Passed to the \code{stats::aggregate} function.
##' @param data [data.frame] dataset containing the observations.
##' @param na.action [function] a function which indicates what should happen when the data contain 'NA' values.
##' Passed to the \code{stats::aggregate} function. 
##' @param na.rm [logical] Should the summary statistics be computed by omitting the missing values.
##' @param columns [character vector] name of the summary statistics to kept in the output.
##' Can be any of, or a combination of:\itemize{
##' \item \code{"observed"}: number of observations with a measurement.
##' \item \code{"missing"}: number of observations with a missing value.
##' \item \code{"mean"}, \code{"mean.lower"} \code{"mean.upper"}: mean with its confidence interval.
##' \item \code{"median"}, \code{"median.lower"} \code{"median.upper"}: median with its confidence interval.
##' \item \code{"sd"}: standard deviation.
##' \item \code{"q1"}, \code{"q3"}, \code{"IQR"}: 1st and 3rd quartile, interquartile range.
##' \item \code{"min"}, \code{"max"}: minimum and maximum observation.
##' \item \code{"predict.lower"}, \code{"predict.upper"}: prediction interval for normally distributed outcome.
##' \item \code{"correlation"}: correlation matrix between the outcomes (when feasible, see detail section).
##' }
##' @param FUN [function] user-defined function for computing summary statistics.
##' It should take a vector as an argument and output a named single value or a named vector.
##' @param which deprecated, use the argument columns instead.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param skip.reference [logical] should the summary statistics for the reference level of categorical variables be omitted?
##' @param digits [integer, >=0] the minimum number of significant digits to be used to display the results. Passed to \code{print.data.frame}
##' @param ... additional arguments passed to argument \code{FUN}.
##'
##' @details This function is essentially an interface to the \code{stats::aggregate} function.
##' 
##' Confidence intervals (CI) and prediction intervals (PI) for the mean are computed via \code{stats::t.test}.
##' Confidence intervals (CI) for the median are computed via \code{asht::medianTest}.
##'
##' Correlation can be assessed when a grouping and ordering variable are given in the formula interface , e.g. Y ~ time|id.
##' 
##' @return A data frame containing summary statistics (in columns) for each outcome and value of the grouping variables (rows). It has an attribute \code{"correlation"} when it was possible to compute the correlation matrix for each outcome with respect to the grouping variable.

## * summarize (examples)
##' @examples
##' ## simulate data in the wide format
##' set.seed(10)
##' d <- sampleRem(1e2, n.times = 3)
##' d$treat <-  sample(LETTERS[1:3], NROW(d), replace=TRUE, prob=c(0.3, 0.3, 0.4) )
##'
##' ## add a missing value
##' d2 <- d
##' d2[1,"Y2"] <- NA
##'
##' ## run summarize
##' summarize(Y1 ~ 1, data = d)
##' summarize(Y1 ~ 1, data = d, FUN = quantile, p = c(0.25,0.75))
##' summarize(Y1+Y2 ~ X1, data = d)
##' summarize(treat ~ 1, skip.reference = FALSE, data = d)
##' 
##' summarize(Y1 ~ X1, data = d2)
##' summarize(Y1+Y2 ~ X1, data = d2, na.rm = TRUE)
##' 
##' ## long format
##' dL <- reshape(d, idvar = "id", direction = "long",
##'              v.names = "Y", varying = c("Y1","Y2","Y3"))
##' summarize(Y ~ time + X1, data = dL)
##'
##' ## compute correlations (single time variable)
##' e.S <- summarize(Y ~ time + X1 | id, data = dL, na.rm = TRUE)
##' e.S
##' attr(e.S, "correlation")
##' 
##' ## compute correlations (composite time variable)
##' dL$time2 <- dL$time == 2
##' dL$time3 <- dL$time == 3
##' e.S <- summarize(Y ~ time2 + time3 + X1 | id, data = dL, na.rm = TRUE)
##' e.S
##' attr(e.S, "correlation")

## * summarize (code)
##' @export
summarize <- function(formula, data, na.action = stats::na.pass, na.rm = FALSE, level = 0.95,
                      columns = c("observed","missing","mean","sd","min","q1","median","q3","max","correlation"),
                      FUN = NULL,
                      which = NULL,
                      skip.reference = TRUE,
                      digits = NULL,
                      ...){

    data <- as.data.frame(data)

    ## ** check and normalize user imput
    name.all <- all.vars(formula)
    if(any(name.all %in% c("XXindexXX"))){
        invalid <- name.all[name.all %in% c("XXindexXX")]
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

    name.Y <- lhs.vars(formula)
    n.Y <- length(name.Y)
    if(n.Y==0){
        stop("Wrong specification of argument \'formula\'. \n",
             "There need to be at least one variable in the left hand side of the formula. \n")
    }
    if(length(grep("|",deparse(formula), fixed = TRUE))>1){
        stop("Wrong specification of argument \'formula\'. \n",
             "There should at most one symbol |. \n")
    }else if(length(grep("|",deparse(formula), fixed = TRUE))==1){
        formula.split <- strsplit(split = "|",deparse(formula),fixed=TRUE)
        formula2 <- stats::as.formula(formula.split[[1]][1])
        name.X <- rhs.vars(formula2)
        if(length(setdiff(name.all,c(name.Y,name.X)))!=1){
            stop("Wrong specification of argument \'formula\'. \n",
                 "There should be exactly one variable on the right hand side of the formula after the symbol |. \n")
        }
        name.id <- trimws(formula.split[[1]][2], which = "both")

        test.between <- stats::setNames(sapply(name.X, function(iXvar){
            max(tapply(data[[iXvar]],data[[name.id]], FUN = function(x){sum(!duplicated(x))}), na.rm = TRUE)
        })==1,name.X)
        
        if(any(test.between)){
            vec.split <- interaction(data[names(which(test.between))], drop = TRUE)
            ls.id <- tapply(as.character(data[[name.id]]), vec.split, unique)
        }else{
            ls.id <- list(unique(as.character(data[[name.id]])))
        }
    }else{
        name.X <- rhs.vars(formula)
        name.id <- NULL
        formula2 <- formula
    }
    n.X <- length(name.Y)
    
    if("which" %in% names(match.call())){
        warning("Argument \'which\' is deprecated. Consider using argument \'columns\' instead. \n")
        columns <- which
    }
    valid.columns <- c("observed","missing","mean","mean.lower","mean.upper","predict.lower","predict.upper","sd","min","q1","median","q3","median.upper","median.lower","IQR","max","correlation")
    columns <- match.arg(columns, choices = valid.columns, several.ok = TRUE)

    dots <- list(...)
    if(length(dots)>0 && is.null(FUN)){
        stop("Unknown arguments \'",paste(names(dots), collapse = "\' \'"),"\'. \n")
    }
    if(!is.null(FUN)){
        try.FUN <- try(FUN(data[[name.Y[1]]]), silent = TRUE)
        if(!inherits(try.FUN,"try-error")){
            if(is.null(names(try.FUN))){
                stop("Argument \'FUN\' should return values with names. \n")
            }else if(any(names(try.FUN) %in% valid.columns)){
                stop("Argument \'FUN\' should not return values named \"",paste(names(try.FUN)[names(try.FUN) %in% valid.columns], collapse = "\" \""),"\" are used internally. \n")
            }
        }
    }
        
    ## ** handle categorical variables
    to.rm <- NULL
    to.add <- NULL
    for(iY in 1:n.Y){ ## iY <- 3
        if(is.character(data[[name.Y[iY]]])||is.factor(data[[name.Y[iY]]])){
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
    name.Y <- unique(c(setdiff(name.Y,to.rm),to.add))
    n.Y <- length(name.Y)

    ## ** compute summary statistics
    out <- NULL
    iFormula <- stats::update(formula2, paste0("XXindexXX~."))
    iData <- cbind(XXindexXX = 1:NROW(data), data)
    for(iY in 1:n.Y){ ## iY <- 1

        iAggregate <- stats::aggregate(iFormula, data=iData, function(x){
            y <- data[x,name.Y[iY]]
            ## *** missing data
            ## as NA
            n.obs <- sum(!is.na(y))
            n.missing <- sum(is.na(y))
            ## as missing rows in the dataset
            if(!is.null(name.id) && "missing" %in% columns){
                if(any(test.between)){
                    iIndex <- which(names(ls.id)==levels(interaction(data[x,names(which(test.between))], drop = TRUE)))
                    n.missing <- n.missing + sum(ls.id[[iIndex]] %in% unique(as.character(data[x,name.id])) == FALSE)

                }else{
                    n.missing <- n.missing + sum(ls.id[[1]] %in% unique(as.character(data[x,name.id])) == FALSE)
                }
            }

            ## *** gather
            if("mean.lower" %in% columns || "mean.upper" %in% columns || "predict.lower" %in% columns || "predict.upper" %in% columns){
                if(all(y %in% 0:1)){
                    tty <- stats::binom.test(x = sum(y==1, na.rm = TRUE), n = sum(!is.na(y)), conf.level = level, alternative = "two.sided")
                }else{
                    tty <- stats::t.test(y, na.rm = na.rm, conf.level = level, alternative = "two.sided")
                }
            }else{
                tty <- list(estimate = NA, parameter = NA, stderr = NA, conf.int = c(NA, NA))
            }
            if(("median.lower" %in% columns || "median.upper" %in% columns) && requireNamespace("asht") && !all(y %in% 0:1)){
                wty <- asht::medianTest(y, conf.level = level, alternative = "two.sided")
            }else{
                wty <- list(estimate = NA, parameter = NA, stderr = NA, conf.int = c(NA, NA))
            }
            
            if(all(is.na(y))){ ## avoid warning when taking min(), e.g. min(NA, na.rm = TRUE)
                iVec <- c("observed" = n.obs,
                          "missing" = n.missing,
                          "mean" = NA,
                          "mean.lower" = NA,
                          "mean.upper" = NA,
                          "predict.lower" = NA,
                          "predict.upper" = NA,
                          "sd" = NA,
                          "min" = NA,
                          "q1" = NA,
                          "median" = NA,
                          "median.lower" = NA,
                          "median.upper" = NA,
                          "q3" = NA,
                          "IQR" = NA,
                          "max" = NA)
            }else{
                iVec <- c("observed" = sum(!is.na(y)),
                          "missing" = n.missing,
                          "mean" = mean(y, na.rm = na.rm),
                          "mean.lower" = tty$conf.int[1], ## as.double(tty$estimate + stats::qt((1-level)/2,tty$parameter) * tty$stderr[1])  - tty$conf.int[1]
                          "mean.upper" = tty$conf.int[2], ## as.double(tty$estimate + stats::qt(1-(1-level)/2,tty$parameter) * tty$stderr[1])  - tty$conf.int[2]
                          "predict.lower" = as.double(tty$estimate + sqrt(n.obs+1) * stats::qt((1-level)/2,tty$parameter) * tty$stderr[1]), 
                          "predict.upper" = as.double(tty$estimate + sqrt(n.obs+1) * stats::qt(1-(1-level)/2,tty$parameter) * tty$stderr[1]), 
                          "sd" = stats::sd(y, na.rm = na.rm),
                          "min" = min(y, na.rm = na.rm),
                          "q1" = unname(stats::quantile(y, prob = 0.25, na.rm = na.rm)),
                          "median" = stats::median(y, na.rm = na.rm),
                          "median.lower" = wty$conf.int[1],
                          "median.upper" = wty$conf.int[2],
                          "q3" = unname(stats::quantile(y, prob = 0.75, na.rm = na.rm)),
                          "IQR" = stats::IQR(y, na.rm = na.rm),
                          "max" = max(y, na.rm = na.rm))
            }
            if(!is.null(FUN)){
                iVec <- c(iVec,do.call(FUN, c(list(y), dots)))
                if(all(y %in% 0:1)){
                    iVec[c("sd","predict.lower","predict.upper","q1","median","median.lower","median.upper","q3","IQR")] <- NA
                }
            }
            return(iVec)
            
        },
        na.action=na.action)
        
        iDF <- cbind(outcome = name.Y[iY],
                     iAggregate[name.X],
                     iAggregate[["XXindexXX"]][,union(setdiff(columns,"correlation"),setdiff(colnames(iAggregate[["XXindexXX"]]),valid.columns)),drop=FALSE])
        
        out <- rbind(out,iDF)
    }
    ## ** correlation
    if(!is.null(name.id) && any(!test.between) && "correlation" %in% columns){ ## id and time variables

        time <- names(which(!test.between)) ## can be several variables
        table.id.time <- do.call(table,stats::setNames(c(list(data[[name.id]]),data[,name.X,drop=FALSE]),
                                                       c(name.id,name.X)))
        if(length(time)>1 && paste(time, collapse = "_X_XX_X_") %in% names(data)){
            stop("Argument \'data\' should not contain a column named \"",paste(time, collapse = "_X_XX_X_"),"\" as this name is used internally. \n")
        }
        
        if(all(table.id.time %in% 0:1)){

            attr(out,"correlation") <- stats::setNames(vector(mode = "list", length = length(name.Y)),
                                                       name.Y)

            for(iY in 1:n.Y){
                attr(out,"correlation")[[iY]] <- stats::setNames(lapply(ls.id, function(iId){ ## iId <- ls.id[[1]]
                    iDataL <- data[data[[name.id]] %in% iId,,drop = FALSE]
                    if(length(time)>1){
                        iDataL[[paste(time, collapse = "_X_XX_X_")]] <- interaction(iDataL[,time, drop=FALSE])
                        Utime <- paste(time, collapse = "_X_XX_X_")
                    }else{
                        Utime <- time
                    }
                    iDataW <- stats::reshape(data = iDataL[,c(name.id, Utime, name.Y[iY])],
                                             direction = "wide", timevar = Utime, idvar = name.id, v.names = name.Y[iY])
                    
                    if(na.rm){
                        return(stats::cor(iDataW[,-1,drop=FALSE], use = "pairwise"))
                    }else{
                        return(stats::cor(iDataW[,-1,drop=FALSE]))
                    }
                }), names(ls.id))
            }
        }
    }

    ## ** export
    if(!is.null(digits)){
        attr(out,"digits") <- digits
    }
    class(out) <- append("summarize",class(out))
    return(out)
}

## * print.summarize
#' @export
print.summarize <- function(x,...){
    if(!is.null(attr(x,"digits")) && ("digits" %in% names(list(...)) == FALSE)){
        print(as.data.frame(x), digits = attr(x,"digits"), ...)
    }else{
        print(as.data.frame(x), ...)
    }
    if(!is.null(attr(x,"correlation"))){
        cat("\n Pearson's correlation: \n")
        ls.cor <- attr(x,"correlation")
        if(length(ls.cor)==1){ ## outcome
            ls.cor <- ls.cor[[1]]
            if(length(ls.cor)==1){ ## group
                ls.cor <- ls.cor[[1]]
            }
        }else{
            for(iY in 1:length(ls.cor)){ ## group
                if(length(ls.cor[[iY]])==1){
                    ls.cor[[iY]] <- ls.cor[[iY]][[1]]
                }
            }
        }
        print(ls.cor, ...)
    }
    return(invisible(NULL))
}

######################################################################
### summarize.R ends here
