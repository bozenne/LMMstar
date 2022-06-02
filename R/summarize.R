### summarize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: Jun  2 2022 (11:30) 
##           By: Brice Ozenne
##     Update #: 144
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
##' @description Compute summary statistics (similar to the SAS macro procmean).
##' This is essentially an interface to the \code{stats::aggregate} function.
##'
##' @param formula [formula] on the left hand side the outcome(s) and on the right hand side the grouping variables.
##' E.g. Y1+Y2 ~ Gender + Gene will compute for each gender and gene the summary statistics for Y1 and for Y2.
##' Passed to the \code{stats::aggregate} function.
##' @param data [data.frame] dataset (in the wide format) containing the observations.
##' @param na.action [function] a function which indicates what should happen when the data contain 'NA' values.
##' Passed to the \code{stats::aggregate} function.
##' @param na.rm [logical] Should the summary statistics be computed by omitting the missing values.
##' @param which [character vector] name of the summary statistics to kept in the output.
##' Can be any of, or a combination of: \code{"observed"} (number of observations with a measurement),
##' \code{"missing"} (number of observations with a missing value), \code{"mean"}, \code{"mean.lower"}, \code{"mean.upper"},
##' \code{"sd"}, \code{"min"},
##' \code{"median"}, \code{"median.lower"}, \code{"median.upper"},
##' \code{"max"}.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param skip.reference [logical] should the summary statistics for the reference level of categorical variables be omitted?
##'
##' @details Confidence intervals for the mean are computed via \code{stats::t.test}
##' and confidence intervals for the median are computed via \code{asht::medianTest}.
##' 
##' @return a data frame containing summary statistics (in columns) for each outcome and value of the grouping variables (rows). It has an attribute \code{"correlation"} when it was possible to compute the correlation matrix for each outcome with respect to the grouping variable.

## * summarize (examples)
##' @examples
##' ## simulate data in the wide format
##' set.seed(10)
##' d <- sampleRem(1e2, n.times = 3)
##'
##' ## add a missing value
##' d2 <- d
##' d2[1,"Y2"] <- NA
##'
##' ## run summarize
##' summarize(Y1 ~ 1, data = d)
##' summarize(Y1+Y2 ~ X1, data = d)
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
                      which = c("observed","missing","mean","sd","min","median","max","correlation"),
                      skip.reference = TRUE){

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

    valid.which <- c("observed","missing","mean","mean.lower","mean.upper","sd","min","median","median.upper","median.lower","max","correlation")
    which <- match.arg(which, choices = valid.which, several.ok = TRUE)

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
            n.missing <- sum(is.na(y))
            ## as missing rows in the dataset
            if(!is.null(name.id) && "missing" %in% which){
                if(any(test.between)){
                    iIndex <- which(names(ls.id)==levels(interaction(data[x,names(which(test.between))], drop = TRUE)))
                    n.missing <- n.missing + sum(ls.id[[iIndex]] %in% unique(as.character(data[x,name.id])) == FALSE)

                }else{
                    n.missing <- n.missing + sum(ls.id[[1]] %in% unique(as.character(data[x,name.id])) == FALSE)
                }
            }

            ## *** gather
            if("mean.lower" %in% which || "mean.upper" %in% which){
                if(all(y %in% 0:1)){
                    tty <- stats::binom.test(x = sum(y==1, na.rm = TRUE), n = sum(!is.na(y)), conf.level = level, alternative = "two.sided")
                }else{
                    tty <- stats::t.test(y, na.rm = na.rm, conf.level = level, alternative = "two.sided")
                }
            }else{
                tty <- list(conf.int = c(NA, NA))
            }
            if(("median.lower" %in% which || "median.upper" %in% which) && requireNamespace("asht") && !all(y %in% 0:1)){
                wty <- asht::medianTest(y, conf.level = level, alternative = "two.sided")
            }else{
                wty <- list(conf.int = c(NA, NA))
            }

            if(all(is.na(y))){ ## avoid warning when taking min(), e.g. min(NA, na.rm = TRUE)
                iVec <- c("observed" = sum(!is.na(y)),
                          "missing" = n.missing,
                          "mean" = NA,
                          "mean.lower" = NA,
                          "mean.upper" = NA,
                          "sd" = NA,
                          "min" = NA,
                          "median" = NA,
                          "median.lower" = NA,
                          "median.upper" = NA,
                          "max" = NA)
            }else{
                iVec <- c("observed" = sum(!is.na(y)),
                          "missing" = n.missing,
                          "mean" = mean(y, na.rm = na.rm),
                          "mean.lower" = tty$conf.int[1],
                          "mean.upper" = tty$conf.int[2],
                          "sd" = stats::sd(y, na.rm = na.rm),
                          "min" = min(y, na.rm = na.rm),
                          "median" = stats::median(y, na.rm = na.rm),
                          "median.lower" = wty$conf.int[1],
                          "median.upper" = wty$conf.int[2],
                          "max" = max(y, na.rm = na.rm))
            }
            if(all(y %in% 0:1)){
                iVec[c("sd","median","median.lower","median.upper")] <- NA
            }
            return(iVec)
            
        },
        na.action=na.action)

        iDF <- cbind(outcome = name.Y[iY],
                     iAggregate[name.X],
                     iAggregate[["XXindexXX"]][,setdiff(which,"correlation"),drop=FALSE])
        
        out <- rbind(out,iDF)
    }

    ## ** correlation
    if(!is.null(name.id) && any(!test.between) && "correlation" %in% which){ ## id and time variables

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
    class(out) <- append("summarize",class(out))
    return(out)
}

## * print.summarize
#' @export
print.summarize <- function(x,...){
    print(as.data.frame(x), ...)
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
