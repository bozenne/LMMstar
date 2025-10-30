### summarizeNA.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  7 2022 (17:13) 
## Version: 
## Last-Updated: okt 24 2025 (12:23) 
##           By: Brice Ozenne
##     Update #: 128
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summarizeNA (documentation)
##' @title Summarize missing data patterns
##' @description Summarize missing data patterns.
##'
##' @param data [data.frame] dataset containing the observations.
##' @param formula [formula] On the left hand side the variable(s) for which the missing data patterns should be evaluated and on the right hand side the grouping variables.
##' E.g. Y1 ~ Gender will compute missing data pattern w.r.t Y1 for each gender.
##' @param repetition [formula] Specify the structure of the data when in the long format: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' When specified the missing data pattern is specific to each variable not present in the formula. 
##' @param sep [character] character used to separate the missing data indicator (0/1) when naming the missing data patterns.
##' @param newnames [character vector of length 4] additional column containing the variable name (only when argument \code{repetition} is used),
##' variables w.r.t. which missing data patterns are identified,
##' frequency of the missing data pattern in the dataset,
##' name of the missing data pattern in the dataset,
##' and number of missing data per pattern.
##' @param keep.data [logical] should the indicator of missing data per variable in the original dataset per pattern be output.
##' @param filter [character] a regular expression passed to \code{grep} to filter the columns of the dataset.
##' Relevant when using \code{.} to indicate all other variables.
##'
##' @return a data frame
##' 
##' @keywords utilities
##'
##' @seealso
##' \code{\link{autoplot.summarizeNA}} for a graphical display.
##' 
##' @examples
##' #### display missing data pattern (wide format) ####
##' data(gastricbypassW, package = "LMMstar")
##' e.SNA <- summarizeNA(gastricbypassW)
##' e.SNA
##' plot(e.SNA)
##'
##' ## only focus on some variables
##' eG.SNA <- summarizeNA(gastricbypassW, filter = "glucagon")
##' eG.SNA
##' plot(eG.SNA)
##' summarizeNA(weight3+glucagonAUC3 ~ 1, data = gastricbypassW)
##'
##' #### display missing data pattern (long format) ####
##' ## example 1 (single group)
##' data(gastricbypassL, package = "LMMstar")
##' e.SNAL <- summarizeNA(gastricbypassL, repetition = ~time|id)
##' e.SNAL
##' plot(e.SNAL, variable = "glucagonAUC")
##' 
##' ## example 2 (two groups)
##' data(calciumL, package = "LMMstar")
##'
##' ## over both groups
##' mp <- summarizeNA(calciumL, repetition = ~visit|girl)
##' plot(mp, variable = "bmd")
##' plot(mp, variable = "bmd", order.pattern = "frequency")
##' plot(mp, variable = "bmd", order.pattern = 5:1)
##'
##' ## per group
##' mp2 <- summarizeNA(bmd ~ grp, data = calciumL, repetition = ~visit|girl)
##' mp2
##' plot(mp2)
##'
##' ## artificially create different patterns in each group
##' calciumL2 <- calciumL[order(calciumL$girl),]
##' calciumL2[calciumL2$girl == 101,"bmd"] <- c(NA,NA,1,1,1)
##' calciumL2[calciumL2$girl == 104,"bmd"] <- c(NA,1,NA,1,NA)
##' mp3 <- summarizeNA(bmd ~ grp, data = calciumL2, repetition = ~visit|girl)
##' mp3
##' plot(mp3)
##' plot(mp3, order.pattern = "n.missing")
##' plot(mp3, order.pattern = "frequency")
##' 

## * summarizeNA (code)
##' @export
summarizeNA <- function(data, formula, repetition = NULL, sep = "",
                        newnames = c("variable","frequency","missing.pattern","n.missing"),
                        filter = NULL, keep.data = TRUE){

    ## ** check and normalize user input

    ## *** newnames
    if(!is.character(newnames) || length(newnames) !=4){
        stop("Argument \'newnames\' should be a character vector of length 4. \n")
    }
    if(any(duplicated(newnames))){
        stop("Argument \'newnames\' should not contain duplicated values. \n")
    }

    ## *** check data
    data <- as.data.frame(data)
    name.all <- names(data)
    if(any(name.all %in% newnames) && keep.data){
        invalid <- name.all[name.all %in% newnames]
        stop("Incorrect argument \'data\': name(s) \"",paste(invalid, collapse = "\" \""),"\" is used internally. \n",
             "Consider renaming the variables in the dataset. \n",
             sep = "")
    }
    internal.name <- c("XXtimeXX") 
    if(any(internal.name %in% name.all)){
        stop("Incorrect argument \'data\': name \"",paste(intersect(internal.name,name.all),collapse = "\", \""),"\" is used internally. \n",
             "Consider renaming the variables in the dataset. \n",
             sep = "")
    }

    ## *** check formula & repetition
    if(missing(formula)){
        
        name.Y0 <- setdiff(names(data), all.vars(repetition))
        if(!is.null(filter)){
            name.Y0 <- grep(filter, name.Y0, value = TRUE)
        }
        formula <- stats::reformulate(termlabels = "1", response = paste(name.Y0, collapse = "+"))
    }
    
    ls.init <- formula2repetition(formula = formula, data = data, repetition = repetition, keep.time = TRUE, filter = filter)
    detail.formula <- ls.init$detail.formula
    detail.repetition <- ls.init$detail.repetition
        
    name.Y <- detail.formula$var$response
    name.X <- detail.formula$var$regressor
    name.cluster <- detail.repetition$var$cluster
    name.time <- detail.repetition$var$time

    if(any(data[name.X]=="")){        
        stop("Variable on the right hand side of the formula should not take value \"\". \n")
    }
    if(!is.null(name.time) && is.null(name.cluster)){
        stop("Missing cluster variable in argument \'repetition\'. \n",
             "Should be something like ~time|cluster. \n")
    }else if(!is.null(name.cluster) && is.null(name.time)){
        stop("Missing repetition variable in argument \'repetition\'. \n",
             "Should be something like ~time|cluster. \n")
    }
    
    ## ** prepare

    ## convert time as factor
    if(!is.null(name.time)){
        for(iVar in name.time){
            if(!is.factor(data[[iVar]])){
                data[[iVar]] <- as.factor(data[[iVar]])
            }
        }
        data.time <- interaction(data[name.time],drop=TRUE)
        Utime <- levels(data.time)
    }        


    ## split data according to grouping levels
    if(!is.null(name.X)){

        ## collapse all grouping variables
        data.UX <- interaction(data[name.X],drop=FALSE)
        ## match unique group variable to original group variables
        test.UX <- !duplicated(data.UX)
        Mlevel.UX <- data[test.UX,name.X,drop=FALSE]
        rownames(Mlevel.UX) <- data.UX[test.UX]
    
        if(is.null(name.time)){
            ls.data <- split(data[name.Y], f = data.UX)
            grid <- Mlevel.UX[names(ls.data),,drop=FALSE]
            
        }else if(!is.null(name.time)){
            ls0.data <- split(cbind(data[c(name.cluster,name.Y)], XXtimeXX = data.time), f = data.UX)
            grid2 <- expand.grid(name.Y, names(ls0.data), stringsAsFactors = FALSE)
            names(grid2) <- c(newnames[1],"group")
            grid <- cbind(grid2[newnames[1]], Mlevel.UX[grid2$group,,drop=FALSE])

            ls.data <- stats::setNames(lapply(1:NROW(grid2), function(iG){ ## iG <- 1
                iData <- stats::reshape(ls0.data[[grid2[iG,"group"]]][c(name.cluster,"XXtimeXX",grid2[iG,newnames[1]])],
                                        direction = "wide", timevar = "XXtimeXX", idvar = name.cluster, varying = Utime)
                iData[[name.cluster]] <- NULL
                return(iData)
            }), paste(grid2[[newnames[1]]],grid2$group,sep=":"))
        }
    }else if(!is.null(name.time)){
        grid <- expand.grid(name.Y, stringsAsFactors = FALSE)
        names(grid) <- newnames[1]
        ls.data <- stats::setNames(lapply(name.Y, function(iY){ ## iY <- name.Y[1]
            iData <- stats::reshape(cbind(data[c(name.cluster,iY)], XXtimeXX = data.time), direction = "wide", timevar = "XXtimeXX", idvar = name.cluster, varying = Utime)
            iData[[name.cluster]] <- NULL
            return(iData)
        }), name.Y)
    }else{
        grid <- NULL
    }

    ## ** warper
    warper.pattern <- function(iData, sep){ ## iData <- ls.data[[1]]
        iMtest <- is.na(iData)*1.0
        iVtest <- nlme::collapse(iMtest, sep = sep, as.factor = TRUE)
        iUpattern <- levels(iVtest)
        iUpattern.nobs <- unname(table(iVtest))
    
        iIndex0.Upattern <- which(duplicated(iVtest)==FALSE)
        iIndex0.Upattern <- stats::setNames(iIndex0.Upattern,iVtest[iIndex0.Upattern])
        iMtest.Upattern <- iMtest[iIndex0.Upattern[iUpattern],,drop=FALSE]

        if(keep.data){
            iOut <- data.frame(as.numeric(iUpattern.nobs), iUpattern, rowSums(iMtest.Upattern), iMtest.Upattern)
            names(iOut) <- c(newnames[2:4], names(iData))
        }else{
            iOut <- data.frame(as.numeric(iUpattern.nobs), iUpattern, rowSums(iMtest.Upattern))
            names(iOut) <- c(newnames[2:4])
        }
        rownames(iOut) <- NULL
        return(iOut)
    }
    
    ## ** identify patterns
    if(is.null(grid)){
        df.pattern <- warper.pattern(data[name.Y], sep = sep)
    }else{
        ls.df.pattern <- lapply(1:NROW(grid), function(iG){ ## iG <- 1
            iOut <- warper.pattern(ls.data[[iG]], sep = sep)
            iOut2 <- cbind(grid[rep(iG,NROW(iOut)),,drop=FALSE], iOut)
            return(iOut2)
        })
        df.pattern <- do.call(rbind,ls.df.pattern)
        rownames(df.pattern) <- NULL
    }

    ## ** export
    attr(df.pattern,"args") <- list(newnames = newnames, keep.data = keep.data, formula = formula, repetition = repetition, sep = sep,
                                    name.Y = name.Y, name.X = name.X, name.cluster = name.cluster, name.time = name.time)

    class(df.pattern) <- append("summarizeNA", class(df.pattern))
    return(df.pattern)
}

##----------------------------------------------------------------------
### summarizeNA.R ends here
