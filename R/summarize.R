### summarize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2020 (11:12) 
## Version: 
## Last-Updated: Dec 18 2021 (19:53) 
##           By: Brice Ozenne
##     Update #: 83
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
##' \code{"missing"} (number of observations with a missing value), \code{"mean"}, \code{"sd"}, \code{"min"}, \code{"median"}, \code{"max"}.
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
##' summarize(Y1+Y2 ~ 1, data = d)
##' summarize(Y1+Y2 ~ X1, data = d)
##' 
##' summarize(Y1+Y2 ~ X1, data = d2)
##' summarize(Y1+Y2 ~ X1, data = d2, na.rm = TRUE)
##' 
##' ## long format
##' dL <- reshape2::melt(d, id.vars = c("id","X1"), measure.var = c("Y1","Y2","Y3"))
##' dL2 <- dL[-(4:5),]
##' 
##' summarize(value ~ variable + X1, data = dL)
##'
##' ## compute correlations
##' e.S <- summarize(value ~ variable + X1 | id, data = dL2, na.rm = TRUE)
##' e.S
##' attr(e.S, "correlation")


## * summarize (code)
##' @export
summarize <- function(formula, data, na.action = stats::na.pass, na.rm = FALSE,
                     which = c("observed","missing","mean","sd","min","median","max","correlation")){

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

    valid.which <- c("observed","missing","mean","sd","min","median","max","correlation")
    which <- match.arg(which, choices = valid.which, several.ok = TRUE)
    
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
            c("observed" = sum(!is.na(y)),
              "missing" = n.missing,
              "mean" = mean(y, na.rm = na.rm),
              "sd" = stats::sd(y, na.rm = na.rm),
              "min" = min(y, na.rm = na.rm),
              "median" = stats::median(y, na.rm = na.rm),
              "max" = max(y, na.rm = na.rm))},
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
        if(all(table.id.time %in% 0:1)){

            attr(out,"correlation") <- stats::setNames(vector(mode = "list", length = length(name.Y)),
                                                       name.Y)
            for(iY in 1:n.Y){
                attr(out,"correlation")[[iY]] <- stats::setNames(lapply(ls.id, function(iId){ ## iId <- ls.id[[1]]
                    iDataL <- data[data[[name.id]] %in% iId,,drop = FALSE]
                    iDataW <- reshape2::dcast(iDataL,
                                              formula = stats::as.formula(paste0(name.id,"~",paste0(time,collapse="+"))),
                                              value.var = name.Y[iY])
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
