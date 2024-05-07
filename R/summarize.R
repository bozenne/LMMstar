### summarize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: maj  7 2024 (10:19) 
##           By: Brice Ozenne
##     Update #: 392
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
##' \item \code{"missing"}: number of missing observations.
##' When specifying a grouping variable, it will also attempt to count missing rows in the dataset.
##' \item \code{"pc.missing"}: percentage missing observations.
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
##' Confidence intervals (CI) for the median are computed via \code{asht::medianTest}.
##'
##' Correlation can be assessed when a grouping and ordering variable are given in the formula interface , e.g. Y ~ time|id.
##' 
##' @return A data frame containing summary statistics (in columns) for each outcome and value of the grouping variables (rows). It has an attribute \code{"correlation"} when it was possible to compute the correlation matrix for each outcome with respect to the grouping variable.
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
##' summarize(Y1 ~ 1, data = d)
##' summarize(Y1 ~ 1, data = d, FUN = quantile, p = c(0.25,0.75))
##' summarize(Y1+Y2 ~ X1, data = d)
##' summarize(treat ~ 1, data = d)
##' summarize(treat ~ 1, skip.reference = FALSE, data = d)
##' 
##' summarize(Y1 ~ X1, data = d2)
##' summarize(Y1+Y2 ~ X1, data = d2, na.rm = TRUE)
##' summarize(. ~ treat, data = d2, na.rm = TRUE, filter = "Y")
##' 
##' #### summarize (long format) ####
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
##'

## * summarize (code)
##' @export
summarize <- function(formula, data, na.action = stats::na.pass, na.rm = FALSE, level = 0.95,
                      columns = c("observed","missing","pc.missing","mean","sd","min","q1","median","q3","max","correlation"),
                      FUN = NULL,
                      skip.reference = TRUE,
                      digits = NULL,
                      filter = NULL,
                      ...){

    data <- as.data.frame(data)
    mycall <- match.call()

    ## ** check and normalize user imput
    name.all <- setdiff(all.vars(formula), ".")
    if(!is.null(filter)){
        data <- data[union(name.all,grep(filter,names(data),value=TRUE))]
    }

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

    detail.formula <- formula2var(formula, data = data)
    ## handle .~Group or Y~.
    if(identical(detail.formula$vars$response,".") & identical(detail.formula$vars$regressor,".")){
        stop("Argument \'formula\' cannot be .~. as the left or right hand side need to be explicit \n",
             "Consider for instance using .~1. \n",sep = "")
    }else if(identical(detail.formula$vars$response,".")){
        name.Y <- setdiff(names(data),c(detail.formula$vars$regressor,detail.formula$var$repetition))
        if(length(name.Y)==0){
            stop("Incorrect argument \'formula\': no variable on the left hand side of the formula. \n")
        }
        if(!is.null(detail.formula$vars$cluster)){
            termlabels <- ifelse(is.null(detail.formula$vars$time),
                                 paste0("1|",detail.formula$vars$cluster),
                                 paste0(paste(detail.formula$vars$time,collapse="+"),"|",detail.formula$vars$cluster))
        }else{
            termlabels <- ifelse(is.null(detail.formula$vars$regressor),"1",paste(detail.formula$vars$regressor,collapse="+"))
        }
        formula <- stats::reformulate(response = paste0(name.Y,collapse="+"), termlabels = termlabels)
        detail.formula <- formula2var(formula)
    }else{
        name.Y <- detail.formula$var$response
        if(length(name.Y)==0){
            stop("Incorrect argument \'formula\': no variable on the left hand side of the formula. \n")
        }
        if(detail.formula$special=="none"){
            if(utils::tail(detail.formula$vars$all,1)=="."){
                formula <- stats::formula(stats::terms(formula, data = data))
                detail.formula <- formula2var(formula)
            }
        }else if(detail.formula$special=="repetition"){
            if(!is.null(detail.formula$vars$time[1]) && detail.formula$vars$time[1]=="."){
                name.X <- setdiff(names(data),c(detail.formula$vars$response,detail.formula$var$cluster))
                if(length(name.X)==0){
                    formula <- stats::reformulate(response = paste0(name.Y,collapse="+"), termlabels = paste0("1|",detail.formula$var$cluster))
                }else{
                    formula <- stats::reformulate(response = paste0(name.Y,collapse="+"), termlabels = paste0(paste(name.X,collapse="+"),"|",detail.formula$var$cluster))
                }
                detail.formula <- formula2var(formula)
            }
        }
    }

    n.Y <- length(name.Y)
    if(n.Y==0){
        stop("Wrong specification of argument \'formula\'. \n",
             "There need to be at least one variable in the left hand side of the formula. \n")
    }
    name.id <- detail.formula$var$cluster

    if(length(name.id)==0){
        name.X <- detail.formula$var$regressor
        formula <- detail.formula$formula$all
    }else if("ranef" %in% names(detail.formula$vars)){
        stop("Wrong specification of argument \'formula\'. \n",
             "Should be something like Y ~ time or Y ~ time + G | cluster. \n")
    }else if(length(name.id)==1){
        name.X <- detail.formula$var$time
        if(is.null(name.X)){
            formula <- stats::as.formula(paste(paste(name.Y, collapse = "+"), "~1"))
        }else{
            formula <- stats::as.formula(paste(paste(name.Y, collapse = "+"), "~", paste(name.X, collapse = "+")))
        }

        test.between <- stats::setNames(sapply(name.X, function(iXvar){
            max(tapply(data[[iXvar]],data[[name.id]], FUN = function(x){sum(!duplicated(x))}), na.rm = TRUE)
        })==1,name.X)
        if(any(test.between)){
            vec.split <- nlme::collapse(data[names(which(test.between))], as.factor = TRUE)
            ls.id <- tapply(as.character(data[[name.id]]), vec.split, unique)
        }else{
            ls.id <- list(unique(as.character(data[[name.id]])))
        }
    }else{
        stop("Wrong specification of argument \'formula\'. \n",
             "There should be exactly one variable on the right hand side of the formula after the symbol |. \n",
             "Something like Y ~ time + G | cluster. \n")
    }

    valid.columns <- c("observed","missing","pc.missing",
                       "mean","mean.lower","mean.upper","predict.lower","predict.upper",
                       "sd","min","q1","median","q3","median.upper","median.lower","IQR","max","correlation")
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

    ## ** between and within variables
    if(!is.null(name.id)){
        time <- names(which(!test.between)) ## can be several variables
        if(length(time)==0){
            table.id.time <- do.call(table,stats::setNames(list(data[[name.id]],rep(1,NROW(data))),c(name.id,"XXXXXX")))
        }else{
            table.id.time <- do.call(table,stats::setNames(c(list(data[[name.id]]),data[,name.X,drop=FALSE]),
                                                           c(name.id,name.X)))
        }
    }else{
        time <- NULL
        table.id.time <- NULL
    }
        
    
    ## ** compute summary statistics
    out <- NULL
    iFormula <- stats::update(formula, paste0("XXindexXX~."))
    iData <- cbind(XXindexXX = 1:NROW(data), data)
    for(iY in 1:n.Y){ ## iY <- 1

        iAggregate <- stats::aggregate(iFormula, data=iData, function(x){
            y <- data[x,name.Y[iY]]
            ## *** missing data
            ## as NA
            n.obs <- sum(!is.na(y))
            n.missing <- sum(is.na(y))
            ## as missing rows in the dataset
            if(!is.null(table.id.time) && all(table.id.time %in% 0:1) && ("missing" %in% columns || "pc.missing" %in% columns)){

                if(any(test.between)){
                    iIndex <- which(names(ls.id)==unique(nlme::collapse(data[x,names(which(test.between)),drop=FALSE], as.factor = FALSE)))
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
                if(na.rm){
                    wty <- asht::medianTest(stats::na.omit(y), conf.level = level, alternative = "two.sided")
                }else{
                    wty <- asht::medianTest(y, conf.level = level, alternative = "two.sided")
                }
            }else{
                wty <- list(estimate = NA, parameter = NA, stderr = NA, conf.int = c(NA, NA))
            }

            if(all(is.na(y))){ ## avoid warning when taking min(), e.g. min(NA, na.rm = TRUE)
                iVec <- c("observed" = n.obs,
                          "missing" = n.missing,
                          "pc.missing" = n.missing/(n.obs+n.missing),
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
                          "pc.missing" = n.missing/length(y),
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

        if(length(time)>1 && paste(time, collapse = "_X_XX_X_") %in% names(data)){
            stop("Argument \'data\' should not contain a column named \"",paste(time, collapse = "_X_XX_X_"),"\" as this name is used internally. \n")
        }
        
        if(all(table.id.time %in% 0:1)){

            attr(out,"correlation") <- stats::setNames(vector(mode = "list", length = length(name.Y)),
                                                       name.Y)

            for(iY in 1:n.Y){ ## iY <- 1 
                attr(out,"correlation")[[iY]] <- stats::setNames(lapply(ls.id, function(iId){ ## iId <- ls.id[[1]]
                    iDataL <- data[data[[name.id]] %in% iId,,drop = FALSE]
                    if(length(time)>1){
                        iDataL[[paste(time, collapse = "_X_XX_X_")]] <- nlme::collapse(iDataL[,time, drop=FALSE], as.factor = TRUE)
                        Utime <- paste(time, collapse = "_X_XX_X_")
                    }else{
                        Utime <- time
                    }
                    iDataW <- stats::reshape(data = iDataL[,c(name.id, Utime, name.Y[iY])],
                                             direction = "wide", timevar = Utime, idvar = name.id, v.names = name.Y[iY])
                    
                    if(na.rm){
                        iCor <- stats::cor(iDataW[,-1,drop=FALSE], use = "pairwise")
                        
                    }else{
                        iCor <- stats::cor(iDataW[,-1,drop=FALSE])
                    }
                    iLevels <- levels(as.factor(iDataL[[Utime]]))
                    if(NROW(iCor)==length(iLevels)){
                        dimnames(iCor) <- list(iLevels,iLevels)
                    }
                    return(iCor)
                }), names(ls.id))
            }
        }
    }

    ## ** export

    if(!is.null(digits)){
        attr(out,"digits") <- digits
    }
    attr(out,"call") <- mycall
    attr(out,"name.Y") <- name.Y
    attr(out,"name.X") <- name.X
    attr(out,"name.time") <- time
    attr(out,"name.id") <- name.id
    class(out) <- append("summarize",class(out))
    return(out)
}



######################################################################
### summarize.R ends here
