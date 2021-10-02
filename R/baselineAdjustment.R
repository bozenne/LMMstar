### baselineAdjustment.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (14:55) 
## Version: 
## Last-Updated: okt  1 2021 (16:35) 
##           By: Brice Ozenne
##     Update #: 30
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * baselineAdjustment (documentation)
##' @title Perform Baseline Adjustment
##' @description Create a new variable based on a time variable and a group variable where groups are constrained to be equal at specific timepoints.  
##' @name baselineAdjustment
##' 
##' @param object [data.frame] dataset
##' @param variable [character] Column in the dataset to be constrained at specific timepoints. 
##' @param repetition [formula] Time and cluster structure, typically \code{~time|id}. See examples below.
##' @param constrain [vector] Levels of the time variable at which the variable is constained.
##' @param new.level [character or numeric] Level used at the constraint. If \code{NULL}, then the first level of the variable argument is used.
##' 
##' @return A vector of length the number of rows of the dataset.
##' 
##' @examples
##' data(ncgsL, package = "LMMstar")
##' 
##' ## baseline adjustment 1
##' ncgsL$treat <- baselineAdjustment(ncgsL, variable = "group",
##'                  repetition= ~ visit|id, constrain = 1)
##' table(treat = ncgsL$treat, visit = ncgsL$visit, group = ncgsL$group)
##'
##' e1.lmm <- suppressWarnings(lmm(cholest~visit*treat,
##'              data=ncgsL, repetition= ~ visit|id,
##'              structure = "CS"))
##' 
##' 
##' ## baseline adjustment 2
##' ncgsL$treat2 <- baselineAdjustment(ncgsL, variable = "group", new.level = "none",
##'                  repetition= ~ visit|id, constrain = 1)
##' table(treat = ncgsL$treat2, visit = ncgsL$visit, group = ncgsL$group)
##'
##' e2.lmm <- suppressWarnings(lmm(cholest~visit*treat2,
##'              data=ncgsL, repetition= ~ visit|id,
##'              structure = "CS"))
##' 
##' 
##' 

## * baselineAdjustment (code)
##' @rdname baselineAdjustment
##' @export
baselineAdjustment <- function(object, variable, repetition, constrain, new.level = NULL){

    ## ** normalize user input
    if(!inherits(object,"data.frame")){
        stop("Argument \'object\' should be a data.frame. \n")
    }
    if(length(variable)!=1){
        stop("Argument \'variable\' must have length 1. \n")
    }
    if(variable %in% names(object) == FALSE){
        stop("Argument \'variable\' does not correspond to any column in argument \'object\'. \n")
    }
    name.repetition <- all.vars(repetition)
    if(any(name.repetition %in% names(object) == FALSE)){
        invalid <- name.repetition[name.repetition %in% names(object) == FALSE]
        stop("Argument \'repetition\' is inconsistent with argument \'object\'. \n",
             "Variable(s) \"",paste(invalid, collapse = "\" \""),"\" could not be found in the dataset. \n",
             sep = "")
    }
    if(!grepl("|",deparse(repetition),fixed = TRUE)){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "No | symbol found so no grouping variable could be defined. \n",
             "Shoud be something like: ~ time|id or group ~ time|id. \n")
    }

    if(length(grepl("|",deparse(repetition),fixed = TRUE))>1){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "The symbol | should only appear once, something like: ~ time|id or group ~ time|id. \n")
    }
    res.split <- strsplit(deparse(repetition),"|", fixed = TRUE)[[1]]
    var.cluster <- trimws(res.split[2], which = "both")
    formula.var <- stats::as.formula(res.split[1])
    var.time <- all.vars(stats::update(formula.var,0~.))

    if(length(var.time)!=1){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "There should be exactly one variable before the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
    }
    if(length(var.cluster)!=1){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "Should have exactly one variable after the grouping symbol (|), something like: ~ time|id or group ~ time|id. \n")
    }
    if(is.factor(object[[var.cluster]])){
        object[[var.cluster]] <- droplevels(object[[var.cluster]])
    }
    test.duplicated <- tapply(object[[var.time]], object[[var.cluster]], function(iT){any(duplicated(iT))})
    if(any(test.duplicated)){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "The time variable (first variable before |) should contain unique values within clusters \n")
    }
    if(is.factor(object[[var.time]])){
        level.time <- levels(object[[var.time]])
    }else{
        level.time <- sort(unique(object[[var.time]]))
    }
    if(is.factor(object[[variable]])){
        level.variable <- levels(object[[variable]])
    }else{
        level.variable <- sort(unique(object[[variable]]))
    }
    if(any(constrain %in% level.time == FALSE)){
        invalid <- constrain[constrain %in% level.time == FALSE]
        stop("Argument \'constrain\' is inconsistent with argument \'object\'. \n",
             "Level(s) \"",paste(invalid, collapse = "\" \""),"\" does not match any value of the variable ",var.time," in the dataset. \n",
             sep = "")
    }
    if(!is.null(new.level) && new.level %in% level.variable){
        stop("Argument \'new.level\' should not correspond to a value of variable ",variable,"in argument \'object\'.")
    }
    
    ## ** perform baseline adjustment
    if(is.null(new.level)){
        out <- factor(object[[variable]], levels = level.variable)
        out[object[[var.time]] %in% constrain] <- level.variable[1]
    }else{
        out <- factor(object[[variable]], levels = c(new.level,level.variable))
        out[object[[var.time]] %in% constrain] <- new.level
    }

    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### baselineAdjustment.R ends here
