### summarizeNA.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  7 2022 (17:13) 
## Version: 
## Last-Updated: jan  3 2023 (17:27) 
##           By: Brice Ozenne
##     Update #: 46
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
##' @param repetition [formula] Specify the structure of the data when in the long format: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' When specified the missing data pattern is specific to each variable not present in the formula. 
##' @param sep [character] character used to separate the missing data indicator (0/1) when naming the missing data patterns.
##' @param newnames [character vector of length 4] additional column containing the variable name (only when argument \code{repetition} is used),
##' the frequency of the missing data pattern in the dataset, the name of the missing data pattern in the dataset, and the number of missing data per pattern.
##' @param keep.data [logical] should the indicator of missing data per variable in the original dataset per pattern be output.
##'
##' @return a data frame
##'
##' @seealso
##' \code{\link{autoplot.summarizeNA}} for a graphical display.
##' 
##' @examples
##' data(gastricbypassW, package = "LMMstar")
##' summarizeNA(gastricbypassW)
##' summarizeNA(gastricbypassW, keep.data = FALSE)
##' 
##' data(gastricbypassL, package = "LMMstar")
##' summarizeNA(gastricbypassL, repetition = ~time|id)
##' 
##' data(calciumL, package = "LMMstar")
##' summarizeNA(calciumL, repetition = ~visit|girl)
##' summarizeNA(calciumL[,c("visit","girl","bmd")], repetition = ~visit|girl)
##' 
##' data(vasscoresW, package = "LMMstar")
##' summarizeNA(vasscoresW)

## * summarizeNA (code)
##' @export
summarizeNA <- function(data, repetition = NULL, sep = "",
                        newnames = c("variable","frequency","missing.pattern","n.missing"),
                        keep.data = TRUE){

    ## ** check and normalize user input

    ## *** check data
    data <- as.data.frame(data)
    name.all <- names(data)
    if(any(name.all %in% newnames) && keep.data){
        invalid <- name.all[name.all %in% newnames]
        stop("Name(s) \"",paste(invalid, collapse = "\" \""),"\" is used internally. \n",
             "Consider renaming the variables in the dataset. \n",
             sep = "")
    }

    ## *** handle repetition
    if(!is.null(repetition)){

        if(!inherits(repetition,"formula")){
            stop("Argument \'repetition\' should be a formula. \n",
                 "Typcally ~time|cluster. \n")
        }

        res.split <- strsplit(deparse(repetition),"|", fixed = TRUE)[[1]]
        if(length(res.split)!=2){
            stop("Incorrect specification of argument \'repetition\'. \n",
                 "The symbol | should only exacly once, something like: ~ time|cluster. \n")
        }
        var.time <- all.vars(stats::as.formula(paste("~",trimws(res.split[1], which = "both"))))
        if(length(var.time)==0){
            stop("Missing time variable in argument \'repetition\'. \n",
                 "Should be something like: ~time|cluster. \n")
        }
        if(any(var.time %in% name.all == FALSE)){
            stop("Mismatch between argument \'repetition\' and argument \'data\'. \n",
                 "Could not find the time variable in the dataset.\n")
        }

        var.cluster <- all.vars(stats::as.formula(paste("~",trimws(res.split[2], which = "both"))))
        if(length(var.cluster)==0){
            stop("Missing cluster variable in argument \'repetition\'. \n",
                 "Should be something like: ~time|cluster. \n")
        }else if(length(var.cluster)>1){
            stop("Too many cluster variables in argument \'repetition\'. \n",
                 "Should be something like: ~time|cluster. \n")
        }
        if(any(var.cluster %in% name.all == FALSE)){
            stop("Mismatch between argument \'repetition\' and argument \'data\'. \n",
                 "Could not find the cluster variable in the dataset.\n")
        }

        name.Y <- setdiff(name.all, c(var.time,var.cluster))
        if(length(name.Y)==0){
            stop("Missing column in argument \'data\'. \n",
                 "There should be at least one column other than the time and cluster variable. \n")
        }

        if(!is.factor(data[[var.time]])){
            data[[var.time]] <- as.factor(data[[var.time]])
        }
        Utime <- levels(data[[var.time]])
            
        ls.data <- stats::setNames(lapply(name.Y, function(iY){ ## iY <- name.Y[1]
            stats::reshape(data[,c(var.cluster,var.time,iY)], direction = "wide", timevar = var.time, idvar = var.cluster, varying = Utime)
        }), name.Y)
    }

    ## ** warper

    warper.pattern <- function(iData, sep){ ## iData <- ls.data[[1]]
        iMtest <- is.na(iData)*1.0
        iVtest <- interaction(as.data.frame(iMtest), sep = sep, drop = TRUE)
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

        return(iOut)
    }
    
    ## ** identify patterns
    if(is.null(repetition)){
        df.pattern <- warper.pattern(data, sep = sep)
    }else{
        ls.df.pattern <- lapply(name.Y, function(iY){
            iDf <- cbind(iY, warper.pattern(ls.data[[iY]], sep = sep))
            names(iDf)[1] <- newnames[1]
            return(iDf)
        })
        df.pattern <- do.call(rbind,ls.df.pattern)
        rownames(df.pattern) <- NULL
    }


    ## ** export
    attr(df.pattern,"args") <- list(newnames = newnames, keep.data = keep.data, repetition = repetition, sep = sep)

    class(df.pattern) <- append("summarizeNA", class(df.pattern))
    return(df.pattern)
}

##----------------------------------------------------------------------
### summarizeNA.R ends here