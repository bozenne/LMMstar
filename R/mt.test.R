### mt.test.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 17 2023 (17:11) 
## Version: 
## Last-Updated: May  9 2024 (12:54) 
##           By: Brice Ozenne
##     Update #: 44
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mt.test (documentation)
##' @title Multiple Student's t-Test
##' @description Perform multiple Student's t-Test via heteroschedastic linear regression
##' and combine the results, possibly adjusted for multiplicity.
##' 
##' @param formula A formula like \code{Y1+Y2+Y3~X|id} with: \itemize{
##' \item the outcome on the left hand side (separated with +)
##' \item the group variable on the right hand side
##' \item a variable identifying each line in the dataset (optional)
##' } 
##' @param data dataset in the wide format. Should inherit from data.frame.
##' @param method [character] type of adjustment for multiple comparisons, one of \code{"none"}, \code{"bonferroni"}, ..., \code{"fdr"}, \code{"single-step"}, \code{"single-step2"}.
##' See \code{\link{confint.Wald_lmm}} for more details.
##' By default \code{"single-step"} when the test statistics have equal degrees of freedom and otherwise \code{"single-step2"}.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param trace [logical] should a message be displayed in the console when there are missing data.
##'
##' @return A data.frame with the estimates, confidence intervals, and p-values relative to each outcome.
##' Depending on the argument \code{method} confidence intervals and p-values may be adjusted for multiple comparisons.
##' The data.frame has an attribute \code{mlmm} containing the underlying regression models.
##'
##' @details In presence of missing values, performs a outcome specific complete case analysis.
##' 
##' @keywords models

## * mt.test (example)
##' @examples
##' data(calciumW, package = "LMMstar")
##'
##' t.test(bmd1 ~ grp, data = calciumW)
##' 
##' mt.test(bmd1+bmd2+bmd3+bmd4+bmd5 ~ grp, data = calciumW)
##' mt.test(bmd1+bmd2+bmd3+bmd4+bmd5 ~ grp|girl, data = calciumW)
##' mt.test(bmd1+bmd2+bmd3+bmd4+bmd5 ~ grp|girl, data = calciumW, method = "none")

## * mt.test (code)
##' @export
mt.test <- function(formula, data, method = NULL, level = 0.95, trace = TRUE){

    ## ** normalize arguments
    data <- as.data.frame(data)

    if(!inherits(formula,"formula")){
        stop("Incorrect argument \'formula\': should inherit from the class formula. \n")
    }
    name.var <- all.vars(formula)
    name.X <- all.vars(stats::delete.response(stats::terms(formula)))
    name.Y <- setdiff(name.var,name.X)

    formula.txt <- deparse(formula)
    if(grepl(pattern = "|", formula.txt, fixed = TRUE)){
        name.id <- trimws(strsplit(formula.txt, split = "|", fixed = TRUE)[[1]][2])        
        if(name.id %in% names(data) == FALSE){
            stop("Mismatch between argument \'formula\' and argument \'data\'. \n",
                 "Could not find the variable \"",name.id,"\" in the dataset. \n")
        }
        name.X <- setdiff(name.X, name.id)
    }else{
        name.id <- "TTclusterTT"
        if(name.id %in% names(data)){
            stop("Incorrect argument \'data\': should not contain a column name \"",name.id,"\" as this name is used internally. \n")
        }
        if(name.id %in% name.var){
            stop("Incorrect argument \'formula\': should not contain \"",name.id,"\" as this name is used internally. \n")
        }
        data$TTclusterTT <- 1:NROW(data)
        
    }
    if(length(name.X)!=1){
        stop("Incorrect argument \'formula\': should contain exactly one regressor. \n")
    }
    if(name.X %in% names(data) == FALSE){
        stop("Mismatch between argument \'formula\' and argument \'data\'. \n",
             "Could not find the variable \"",name.X,"\" in the dataset. \n")
    }
    if(any(name.Y %in% names(data) == FALSE)){
        stop("Mismatch between argument \'formula\' and argument \'data\'. \n",
             "Could not find the variable(s) \"",paste(name.Y[name.Y %in% names(data) == FALSE], collapse = "\" \""),"\" in the dataset. \n")
    }
    if("TTtimeTT" %in% names(data)){
        stop("Incorrect argument \'data\': should not contain a column name \"TTtimeTT\" as this name is used internally. \n")
    }     
    data <- data[c(name.Y,name.X,name.id)]
    if(any(is.na(data)) && trace){
        cat("Argument \'data\' contains ",sum(is.na(data))," missing values. \n", sep = "")
    }
    if(is.numeric(data[[name.X]])){
        if(any(data[[name.X]] %in% 0:1 == FALSE)){
            stop("The regressor on the right hand side of the formula must be a binary variable or take at most two levels. \n")
        }
        effects <- paste0(name.X,"=0")
    }else{
        if(is.character(data[[name.X]])){
            data[[name.X]] <- as.factor(data[[name.X]])
        }
        effects <- paste0(name.X, utils::tail(levels(data[[name.X]]),1),"=0")
    }
    if(length(unique(data[[name.X]]))!=2){
        stop("The regressor on the right hand side of the formula must have exactly two distinct values. \n")
    }
    

    ## ** reshape
    name.time <- Reduce(x=name.Y, f=.commonString)
    if(is.na(name.time)){        
        name.time <- "rep"
        levels.time <- name.Y
    }else{
        levels.time <- gsub(name.Y, pattern = name.time, replacement = "", fixed = TRUE)
    }
    if(name.time %in% names(data)){
        stop("Incorrect argument \'data\': should not contain a column name \"",name.time,"\" as this name is used internally. \n")
    }
    dataL <- stats::reshape(data, direction = "long", idvar = c(name.X,name.id),
                            varying  = 1:length(name.Y), v.names = "value", timevar = name.time)
    dataL[[name.time]] <- factor(dataL[[name.time]], labels = name.Y)

    ## ** t-test
    ffx.txt <- stats::as.formula(paste("~", name.X))
    fft.txt <- stats::as.formula(paste("~", name.time,"|",name.id))
    out.mlmm <- mlmm(stats::update(value~.,ffx.txt), structure = do.call(IND, list(ffx.txt)),
                     data = dataL, trace = FALSE,
                     effects = effects, by = name.time, repetition = fft.txt)

    ## ** extract estimates + adjustment for multiple comparisons
    test.sameDF <- all(abs(out.mlmm$univariate$df-round(mean(out.mlmm$univariate$df)))<0.1)
    if(is.null(method)){
        if(test.sameDF){
            method <- "single-step"
        }else{
            method <- "single-step2"
        }            
    }
    out <- model.tables(out.mlmm, method = method, level = level)
    
    ## ** export
    attr(out, "mlmm") <- out.mlmm
    return(out)
}


##----------------------------------------------------------------------
### mt.test.R ends here
