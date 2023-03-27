### lmmCC.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 27 2023 (17:33) 
## Version: 
## Last-Updated: mar 27 2023 (18:39) 
##           By: Brice Ozenne
##     Update #: 47
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
##' 
##' @param formula [formula] Specify the model for the mean.
##' On the left hand side the outcome and on the right hand side the covariates affecting the mean value.
##' E.g. Y ~ Gender + Gene.
##' @param repetition [formula] Specify the structure of the data: the time/repetition variable and the grouping variable, e.g. ~ time|id.
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param lm.change [logical] Should a linear model on the change in outcome be estimated. Only possible with two repetitions.
##' Will match the mixed model if the later includes repetition-dependent effects for all covariates. 
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##' @param trace [interger, >0] Show the progress of the execution of the function.
##' @param control [list] Control values for the optimization method.
##' The element \code{optimizer} indicates which optimizer to use and additional argument will be pass to the optimizer.
##'
##' @return A \code{lmmCC} object, which inherits from \code{lmm}.

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
##' e.lm3 <- lm(change3 ~ X1 + X2, data = dW)
##' 
##' #### 2- fit Linear Mixed Model ####
##' dL <- reshape(dW[,c("id","X1","X2","Y1","Y2","Y3")],
##'               direction = "long", 
##'               varying = c("Y1","Y2","Y3"), sep = "", idvar = "id")
##' dL$time <- as.character(dL$time)
##' 
##' e.lmm2 <- lmmCC(Y ~ time + time*X1 + time*X2, repetition = ~time|id,
##'                 data = dL[dL$time %in% 1:2,])
##' e.lmm3 <- lmmCC(Y ~ time + time*X1 + time*X2, repetition = ~time|id,
##'                 data = dL[dL$time %in% c(1,3),])
##' e.lmm3.bis <- lmmCC(Y ~ time + time*X1 + time*X2, repetition = ~time|id,
##'                 data = dL[dL$time %in% c(1,3),], lm.change = TRUE)
##'
##' #### 3- compare ####
##' model.tables(e.lmm2)
##' summary(e.lm2)$coef
##' 
##' model.tables(e.lmm3)
##' summary(e.lm3)$coef
##' summary(e.lmm3.bis$lm.change)
##' head(e.lmm3.bis$lm.change$data)
##' head(e.lmm3.bis$data[order(e.lmm3.bis$data$id,e.lmm3.bis$data$time),])

## * lmmCC (code)
##' @export
lmmCC <- function(formula, repetition, data,
                  lm.change = FALSE,
                  df = NULL, trace = TRUE, control = NULL){

    ## ** normalize arguments
    ## formula
    if(!inherits(repetition,"formula")){
        stop("Argument \'repetition\' must be of class formula, something like: ~ time|cluster or strata ~ time|cluster. \n")
    }
    ## repetition
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
    
    ## ** complete case dataset
    Uvars <- unique(c(all.vars(formula),all.vars(repetition)))

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
    out <- lmm(formula = formula, repetition = repetition, data = data.NNA,
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
        Xvars <- all.vars(stats::delete.response(stats::terms(stats::as.formula(formula))))
        Xvars0 <- setdiff(Xvars,var.time)
        if(length(Xvars0)>0){
            test.cstX <- by(data.NNA[Xvars0],data.NNA[[var.cluster]],function(iData){NROW(unique(iData))})
            if(any(test.cstX!=1)){
                stop("Covariates should be constant within cluster to be able to fit a univariate linear model. \n",
                     "Consider setting argument \'lm.change\' to FALSE. \n")
            }
        }
        dataW.NNA <- reshape(data.NNA, direction = "wide", timevar = var.time,
                             idvar = var.cluster, v.names = var.outcome, sep = ".")
        dataW.NNA$change <- dataW.NNA[[paste(var.outcome,Utime[2],sep=".")]] - dataW.NNA[[paste(var.outcome,Utime[1],sep=".")]]

        if(length(Xvars0)>0){
            formula.lm <- paste0("change~",paste(Xvars0, collapse = "+"))
        }else{
            formula.lm <- "change~1"
        }
        out$lm.change <- lm(formula.lm, data = dataW.NNA)
        out$lm.change$data <- dataW.NNA
    }

    ## ** export
    class(out) <- append("lmmCC",class(out))
    return(out)
}

##----------------------------------------------------------------------
### lmmCC.R ends here
