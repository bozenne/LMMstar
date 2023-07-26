### lmmCC.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 27 2023 (17:33) 
## Version: 
## Last-Updated: jul 26 2023 (11:07) 
##           By: Brice Ozenne
##     Update #: 184
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * lmmCC (documentation)
##' @title Fit Linear Mixed Model on Complete Case data
##' @description Fit a linear mixed model on the complete case data.
##' Mostly useful as a sanity check, to match the results of a univariate analysis on the change.
##' @name lmmCC
##' 
##' @param object [formula] Specify the model for the mean.
##' On the left hand side the outcome and on the right hand side the covariates affecting the mean value.
##' E.g. Y ~ Gender + Gene.
##' @param repetition [formula] Specify the structure of the data: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param lm.change [logical] Should a linear model on the change in outcome be estimated. Only possible with two repetitions.
##' Will match the mixed model if the later includes repetition-dependent effects for all covariates. 
##' @param name.time [character] name of the time variable.
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##' @param trace [interger, >0] Show the progress of the execution of the function.
##' @param control [list] Control values for the optimization method.
##' The element \code{optimizer} indicates which optimizer to use and additional argument will be pass to the optimizer.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A \code{lmmCC} object, which inherits from \code{lmm}.
##'
##' @keywords models

## * lmmCC (examples)
##' @examples
##' #### 1- simulate data in the wide format ####
##' set.seed(10)
##' dW <- sampleRem(100, n.times = 3, format = "wide")
##' dW$Y3[1:10] <- NA
##' dW$change2 <- dW$Y2 - dW$Y1
##' dW$change3 <- dW$Y3 - dW$Y1
##' 
##' e.lm2 <- lm(change2 ~ X1 + X2, data = dW)
##' summary(e.lm2)$coef
##' e.lm3 <- lm(change3 ~ X1 + X2, data = dW)
##' summary(e.lm3)$coef
##' 
##' #### 2- complete case LMM from LM ####
##' e.lmmCC2 <- lmmCC(e.lm2, repetition = change2~Y2-Y1)
##' model.tables(e.lmmCC2)
##' e.lmmCC3 <- lmmCC(e.lm3, repetition = change3~Y3-Y1)
##' model.tables(e.lmmCC3)
##' 
##' #### 3- complete case LMM ####
##' dL <- reshape(dW[,c("id","X1","X2","Y1","Y2","Y3")],
##'               direction = "long", 
##'               varying = c("Y1","Y2","Y3"), sep = "", idvar = "id")
##' dL$time <- as.character(dL$time)
##' 
##' e.lmm2 <- lmmCC(Y ~ time + time*X1 + time*X2, repetition = ~time|id,
##'                 data = dL[dL$time %in% 1:2,])
##' model.tables(e.lmm2)
##' e.lmm3.bis <- lmmCC(Y ~ time + time*X1 + time*X2, repetition = ~time|id,
##'                 data = dL[dL$time %in% c(1,3),], lm.change = TRUE)
##' model.tables(e.lmm3.bis)
##' e.lmm3.bis$lm.change
##'
##' @export
`lmmCC` <-
    function(object, ...) UseMethod("lmmCC")

## * lmmCC.formula (code)
##' @export
##' @rdname lmmCC
lmmCC.formula <- function(object, repetition, data,
                          lm.change = FALSE,
                          df = NULL, trace = TRUE, control = NULL, ...){

    ## ** normalize arguments
    ## repetition
    if(!inherits(repetition,"formula")){
        stop("Argument \'repetition\' must be of class formula, something like: ~ time|cluster or strata ~ time|cluster. \n")
    }
    res.split <- strsplit(deparse(repetition),"|", fixed = TRUE)[[1]]
    if(length(res.split)!=2){
        stop("Incorrect specification of argument \'repetition\'. \n",
             "The symbol | should only exacly once, something like: ~ time|cluster or strata ~ time|cluster. \n")
    }
    var.cluster <- trimws(res.split[2], which = "both")
    if(var.cluster %in% names(data) == FALSE){
        stop("Argument \'repetition\' incompatible with argument \'data\'. \n",
             "Could not find the cluster variable \"",var.cluster,"\" in argument \'data\'. \n")
    }
    var.time <- all.vars(stats::delete.response(stats::terms(stats::as.formula(res.split[1]))))
    if(var.time %in% names(data) == FALSE){
        stop("Argument \'repetition\' incompatible with argument \'data\'. \n",
             "Could not find the time variable \"",var.time,"\" in argument \'data\'. \n")
    }
    ## data
    data <- as.data.frame(data)
    if(is.factor(data[[var.time]])){
        data[[var.time]] <- droplevels(data[[var.time]])
    }
    if(is.factor(data[[var.cluster]])){
        data[[var.cluster]] <- droplevels(data[[var.cluster]])
    }

    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    ## ** complete case dataset
    Uvars <- unique(c(all.vars(object),all.vars(repetition)))

    ## identify clusters with missing values
    test.NA1 <- rowSums(is.na(data[,Uvars,drop=FALSE]))>0
    id.NA1 <- unique(data[[var.cluster]][test.NA1])

    ## identify clusters with missing timepoints
    test.NA2 <- rowSums(table(data[[var.cluster]],data[[var.time]])==0)>0
    id.NA2 <- unique(data[[var.cluster]][test.NA2])

    ## remove missing values
    id.NA <- union(id.NA1,id.NA2)
    if(length(id.NA)>0){
        if(trace>0){
            cat("Remove ",length(id.NA)," clusters (",sum(data[[var.cluster]] %in% id.NA)," observations) \n",
                " - ",sum(test.NA1)," observations with missing data (",length(id.NA1)," clusters) \n",
                " - ",sum(test.NA2)," missing repetitions (",length(id.NA2)," clusters) \n",
                sep="")
        }
        data.NNA <- data[data[[var.cluster]] %in% id.NA == FALSE,]
    }else{
        data.NNA <- data
    }

    ## update factor levels
    if(is.factor(data.NNA[[var.cluster]])){
        data.NNA[[var.cluster]] <- droplevels(data.NNA[[var.cluster]])
    }
        
    ## ** fit LMM
    out <- lmm(formula = object, repetition = repetition, data = data.NNA,
               structure = "UN", method.fit = "REML", df = df, type.information = "observed",
               trace = trace-1, control = control)

    ## ** fit LM
    if(lm.change){
        var.outcome <- Uvars[1]
        
        if(is.factor(data[[var.time]])){
            Utime <- levels(data[[var.time]])
        }else{
            Utime <- sort(unique(data[[var.time]]))
        }
        if(length(Utime)!=2){
            stop("Can only fit a univariate linear model on the change with two timepoint. \n",
                 "Consider setting argument \'lm.change\' to FALSE. \n")
        }

        var.X <- attr(out$design$mean,"variable")
        n.level <- stats::setNames(rep(2, length(var.X)), var.X)
        if(length(out$xfactor$mean)>0){
            n.level[names(out$xfactor$mean)] <- (lengths(out$xfactor$mean)-1)*2
        }        
        if(NCOL(out$design$mean) != sum(n.level)){
            warning("Number of coefficients in the mixed model is not twice the number of variables. \n",
                    "Results will likely differ between lmmCC and lm due to missing interactions. \n")
        }
        Xvars <- all.vars(stats::delete.response(stats::terms(stats::as.formula(object))))
        Xvars0 <- setdiff(Xvars,var.time)
        if(length(Xvars0)>0){
            test.cstX <- by(data.NNA[Xvars0],data.NNA[[var.cluster]],function(iData){NROW(unique(iData))})
            if(any(test.cstX!=1)){
                stop("Covariates should be constant within cluster to be able to fit a univariate linear model. \n",
                     "Consider setting argument \'lm.change\' to FALSE. \n")
            }
        }
        dataW.NNA <- stats::reshape(data.NNA, direction = "wide", timevar = var.time,
                                    idvar = var.cluster, v.names = var.outcome, sep = ".")
        dataW.NNA$change <- dataW.NNA[[paste(var.outcome,Utime[2],sep=".")]] - dataW.NNA[[paste(var.outcome,Utime[1],sep=".")]]

        if(length(Xvars0)>0){
            formula.lm <- paste0("change~",paste(Xvars0, collapse = "+"))
        }else{
            formula.lm <- "change~1"
        }
        out$lm.change <- stats::lm(formula.lm, data = dataW.NNA)
        out$lm.change$data <- dataW.NNA
    }

    ## ** export
    class(out) <- append("lmmCC",class(out))
    return(out)
}

## * lmmCC.lm (code)
##' @export
##' @rdname lmmCC
lmmCC.lm <- function(object, repetition, data,
                     name.time = "time",
                     df = NULL, trace = TRUE, control = NULL, ...){

    ## ** normalize arguments
    ## data
    if(missing(data)){
        data <- try(eval(object$call$data), silent = TRUE)
        if(inherits(data,"try-error")){
            stop("Could not recover the original dataset. Consider specifying the argument \'data\'. \n")
        }
    }else{
        test <- stats::update(object, data = data)
        if(abs(logLik(test)-logLik(object))>1e-12){
            warning("Argument \'data\' does not match argument \'object\'. \n",
                    "Difference in log-likelihood of the corresponding linear model of ",as.numeric(logLik(test)-logLik(object)),".\n")
        }
    }
    data <- as.data.frame(data)

    ## repetition
    if(!is.list(repetition)){
        repetition <- list(repetition)
    }
    if(length(repetition)>2){
        stop("Cannot handle more than 2 time-varying variables. \n")
    }
    if(any(sapply(repetition, inherits,"formula")==FALSE)){
        stop("Argument \'repetition\' should be formula or list of formula, one for each time varying variable. \n",
             "Typically Ychange ~ Yfollowup - Ybaseline. \n")
    }
    repetition.char <- lapply(repetition, deparse)
    test.id <- sapply(repetition.char, grepl, pattern = "|", fixed = TRUE)
    if(any(test.id)){
        ls.var.cluster <- lapply(repetition.char[test.id], function(iRep){trimws(strsplit(iRep, split = "|", fixed = TRUE)[[1]][2])})
        if(any(lengths(ls.var.cluster)!=1)){
            stop("Incorrect argument \'formula\': there must be a single cluster variable per formula. \n",
                 "Typically Ychange ~ Yfollowup - Ybaseline|id. \n")
        }
        var.cluster <- unique(unlist(ls.var.cluster))
        if(any(lengths(ls.var.cluster)!=1)){
            stop("Incorrect argument \'formula\': the cluster variable must be the same for all formula. \n",
                 "Typically Ychange ~ Yfollowup - Ybaseline|id and Zchange ~ Zfollowup - Zbaseline|id. \n")
        }        
        if(var.cluster %in% names(data) == FALSE){
            stop("Inconsistency between argument \'data\' and argument \'formula\'. \n",
                 "Could not find a column named \"",var.cluster,"\" in \'data\'. \n")
        }
        repetition <- sapply(repetition.char, function(iRep){stats::as.formula(strsplit(iRep, split = "|", fixed = TRUE)[[1]][1])})
    }else{
        var.cluster <- "CCclusterCC"
        if(var.cluster %in% names(data)){
            stop("Argument \'data\' should not contain a column named \"",var.cluster,"\" as this name is used internally by the lmm function. \n")
        }
        data$CCclusterCC <- 1:NROW(data)
    }
    ls.XY <- lapply(repetition, all.vars)
    if(any(lengths(ls.XY)!=3)){
        stop("Argument \'repetition\' should be formula or list of formula, each containing three variables. \n",
             "Typically Ychange ~ Yfollowup - Ybaseline. \n")
    }
    ls.X <- lapply(repetition, function(iFormula){all.vars(stats::delete.response(stats::terms(stats::as.formula(iFormula))))})
    if(any(lengths(ls.X)!=2)){
        stop("Argument \'repetition\' should be formula or list of formula, each containing three variables. \n",
             "Typically Ychange ~ Yfollowup - Ybaseline. \n")
    }
    M.XY <- do.call(rbind,ls.XY)
    if(any(as.character(M.XY) %in% names(data) == FALSE)){
        stop.miss <- as.character(M.XY)[as.character(M.XY) %in% names(data) == FALSE]
        stop("Inconsistency between argument \'repetition\' and argument \'data\'. \n",
             "Could not find variable(s) \"",paste(stop.miss, collapse = "\", \""),"\" in \'data\'. \n")
    }
    object.vars <- all.vars(object$terms)
    if(any(M.XY[,1] %in% object.vars == FALSE)){
        stop.miss <- M.XY[,1][M.XY[,1] %in% object.vars == FALSE]
        stop("Inconsistency between argument \'repetition\' and argument \'object\'. \n",
             "Variable(s) \"",paste(stop.miss, collapse = "\", \""),"\" are not used in \'object\'. \n")

    }
    if(name.time %in% names(data)){
        stop("Argument \'data\' should not contain a column named \"",name.time,"\" as this name is used internally by the lmm function. \n")
    }

    ## dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** move from the wide to the long format
    dataW <- data[,c(var.cluster,setdiff(object.vars,M.XY[,1]),M.XY[,2],M.XY[,3]),drop=FALSE]    
    ## find unique names
    vec.CS <- apply(M.XY,1,function(iVec){.commonString(iVec[2],iVec[3])})
    if(any(is.na(vec.CS))){
        stop.NA <- unlist(vec.CS)[unlist(vec.CS) %in% names(data)]
        stop("Incorrect argument \'repetition\'. \n",
             "Could not find a common name between the time varying variables. \n")
    }
    if(any(c(vec.CS,paste(vec.CS, collapse=".")) %in% names(data))){
        stop.com <- c(vec.CS,paste(vec.CS, collapse="."))[c(vec.CS,paste(vec.CS, collapse=".")) %in% names(data)]
        stop("Column name(s) \"",paste(stop.com, collapse ="\" \""),"\" already exist in argument \'data\'. \n",
             "These name(s) will be used internally and may conflict with existing data.\n")
    }
    index.outcome <- which(M.XY[,1] %in% object.vars[1])
    if(length(index.outcome)==0){
        stop("Argument \'repetition\' should specify how the outcome (Ychange) can be split into time-specific variable (Ybaseline, Yfollowup) \n",
             "Something like Ychange ~ Yfollowup - Ybaseline")
    }
    level.time <- rev(unlist(ls.X))
    name.outcome <- paste(vec.CS,collapse="")
    dataL <- stats::reshape(dataW, direction = "long",
                            varying = level.time,
                            v.names = name.outcome,
                            timevar = name.time,
                            idvar = var.cluster)
    index.X <- setdiff(object.vars[-1], M.XY[,1])
    if(length(repetition)==2){
        dataL[[name.time]] <- factor(dataL[[name.time]], labels = level.time)
        ## dataL[[name.time3]] <- as.numeric(dataL[[name.time]]==level.time[3])
        ## dataL[[name.time4]] <- as.numeric(dataL[[name.time]]==level.time[4])
        
    }
    ## ** fit LMM
    formula.txt <- paste0(name.outcome,"~",name.time)
    if(length(index.X)>0){
        formula.txt <- paste0(formula.txt,"+",name.time,"*(",paste(index.X, collapse = "+"),")")
    }
    repetition.txt <- paste0("~",name.time,"|",var.cluster)
    out <- lmmCC(object = stats::as.formula(formula.txt),
                 repetition = stats::as.formula(repetition.txt),
                 data = dataL,
                 lm.change = FALSE,
                 df = df, trace = trace, control = control)
    out$lm.change <- object
    
    ## ** export
    return(out)
}

## * .commonString
## find the common consecutive part between two strings
## copied from https://stackoverflow.com/questions/35381180/identify-a-common-pattern
## .commonString("aaachang2","aaabbb")
## .commonString("aaa235change2","aaachangebbb")
## .commonString("abcdef","xyz")
.commonString <- function(string1, string2){

    if(is.na(string1) || is.na(string2)){return(NA)}
    A <- strsplit(string1, "")[[1]]
    B <- strsplit(string2, "")[[1]]

    L <- matrix(0, length(A), length(B))
    ones <- which(outer(A, B, "=="), arr.ind = TRUE)
    ones <- ones[order(ones[, 1]), ,drop=FALSE]
    if(NROW(ones)==0){return(NA)}
    for(i in 1:NROW(ones)) {
        v <- ones[i, , drop = FALSE]
        L[v] <- ifelse(any(v == 1), 1, L[v - 1] + 1)
    }
    out <- paste0(A[(-max(L) + 1):0 + which(L == max(L), arr.ind = TRUE)[1]], collapse = "")
    return(out)
}
##----------------------------------------------------------------------
### lmmCC.R ends here
