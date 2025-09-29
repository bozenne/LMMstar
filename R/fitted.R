### fitted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  8 2021 (17:09) 
## Version: 
## Last-Updated: sep 29 2025 (13:23) 
##           By: Brice Ozenne
##     Update #: 438
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * fitted.lmm (documentation)
##' @title Fitted Outcome Value by a Linear Mixed Model.
##' @description Evaluate the expected outcome conditional on covariates
##' or the expected missing outcome value conditional on observed outcomes and covariates
##' based on a linear mixed model. 
##'
##' @param object a \code{lmm} object.
##' @param newdata [data.frame] the covariate values for each cluster.
##' @param type [character] by default fitted values are output (\code{NULL}).
##' Can also output the expected outcome (for missing outcomes) based on covariates and other outcome values from the same cluster (\code{"impute"}),
##' the change or expected change between baseline and each follow-up (\code{"change"}),
##' or the area under the curve of the outcome (\code{"auc"}, require a numeric repetition variable).
##' @param se [character] passed to \code{predict.lmm} to evaluate the standard error of the fitted value, expected outcome, change in expected outcome, or area under the curve.
##' @param df [logical] should a Student's t-distribution be used to model the distribution of the predicted mean. Otherwise a normal distribution is used.
##' @param format [character] should the prediction be output
##' in a matrix format with clusters in row and timepoints in columns (\code{"wide"}),
##' or in a data.frame/vector with as many rows as observations (\code{"long"})
##' @param keep.data [logical] should the dataset relative to which the predictions are evaluated be output along side the predicted values?
##' Only possible in the long format.
##' @param seed [integer, >0] random number generator (RNG) state used when starting imputation. If NULL no state is set.
##' @param simplify [logical] simplify the data format (vector instead of data.frame) and column names (no mention of the time variable) when possible.
##' @param ... additional argument passed the \code{\link{predict.lmm}}.
##'
##' @details Essentially a wrapper function of \code{\link{predict.lmm}} where the fitted value are stored in the outcome column of the dataset instead of in a separate column.
##' 
##' @return When \code{format="wide"}, a data.frame with as many rows as clusters.
##' When \code{format="long"}, a data.frame with as many rows as observations (\code{keep.data==TRUE})
##' or a vector of length the number of observations (\code{keep.data==TRUE}).
##' 
##' @keywords methods
##' 
##' @examples
##' #### single arm trial ####
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL <- gastricbypassL[order(gastricbypassL$id,gastricbypassL$visit),]
##' gastricbypassL$weight0 <- unlist(tapply(gastricbypassL$weight,gastricbypassL$id,
##' function(x){rep(x[1],length(x))}))
##' 
##' eUN.lmm <- lmm(glucagonAUC ~ visit + weight0, repetition = ~visit|id,
##'                data = gastricbypassL, df = FALSE)
##'
##' ## fitted mean (conditional on covariates only)
##' fitted(eUN.lmm)
##' fitted(eUN.lmm, newdata = data.frame(visit = "3", weight0 = 0))
##' fitted(eUN.lmm, newdata = data.frame(visit = "3", weight0 = 0),
##'        keep.data = TRUE)
##' 
##' ## fitted outcome value (conditional on covariates and covariates)
##' fitted(eUN.lmm, type = "outcome")
##' gastricbypassL.O <- fitted(eUN.lmm, type = "outcome", keep.data = TRUE)
##'
##' if(require(ggplot2)){
##' gg.outcome <- ggplot(gastricbypassL.O,
##'                      aes(x=time, y = glucagonAUC, color = impute, group = id))
##' gg.outcome <- gg.outcome + geom_point() + geom_line()## + facet_wrap(~id)
##' gg.outcome
##' }
##'
##' tapply(gastricbypassL.O$glucagonAUC, gastricbypassL.O$time, mean)
##' effects(eUN.lmm, variable = NULL)
##' 
##' ## fitted change value (conditional on covariates and covariates)
##' gastricbypassL.C <- fitted(eUN.lmm, type = "change", keep.data = TRUE)
##'
##' if(require(ggplot2)){
##' gg.change <- ggplot(gastricbypassL.C,
##'                     aes(x=time, y = glucagonAUC, color = impute, group = id))
##' gg.change <- gg.change + geom_point() + geom_line()
##' gg.change
##' }
##' 
##' tapply(gastricbypassL.C$glucagonAUC, gastricbypassL.O$time, mean)
##' effects(eUN.lmm, type = "change", variable = NULL)
##' 
##' ## fitted auc (conditional on covariates and covariates)
##' gastricbypassL.AUC <- fitted(eUN.lmm, type = "auc", keep.data = TRUE)
##'
##' if(require(ggplot2)){
##' gg.auc <- ggplot(gastricbypassL.AUC,
##'                     aes(x = "auc", y = glucagonAUC, color = impute))
##' gg.auc <- gg.auc + geom_point()
##' gg.auc
##' }
##'
##' mean(gastricbypassL.AUC$glucagonAUC)
##' effects(eUN.lmm, type = "auc", variable = NULL)
##' 
##' #### two arm trial ####
##' \dontrun{
##' if(require(nlmeU)){
##' data(armd.wide, package = "nlmeU")
##' armd.long <- reshape(armd.wide, direction  = "long",
##'       idvar = "subject", varying = paste0("visual",c(0,4,12,24,52)),
##'      times = c(0,4,12,24,52), timevar =  "week", v.names = "visual")
##' armd.long$week <- factor(armd.long$week, levels = c(0,4,12,24,52))
##' armd.longNNA <- armd.long[!is.na(armd.long$lesion),]
##' 
##' eUN2.lmm <- lmm(visual ~ treat.f*week + lesion,
##'                repetition = ~week|subject, structure = "UN",
##'                data = armd.longNNA)
##' 
##' ## fitted outcome value (conditional on covariates and covariates)
##' armd.O <- fitted(eUN2.lmm, type = "outcome", keep.data = TRUE)
##'
##' gg2.outcome <- ggplot(armd.O,
##'                      aes(x=week, y = visual, color = impute, group = subject))
##' gg2.outcome <- gg2.outcome + geom_point() + geom_line() + facet_wrap(~treat.f)
##' gg2.outcome
##' 
##' aggregate(visual ~ week + treat.f, FUN = mean, data = armd.O)
##' effects(eUN2.lmm, variable = "treat.f") ## mismatch due to adjustment on lesion
##' 
##' ## fitted change value (conditional on covariates and covariates)
##' armd.C <- fitted(eUN2.lmm, type = "change", keep.data = TRUE)
##'
##' gg.change <- ggplot(armd.C,
##'                     aes(x=week, y = visual, color = impute, group = subject))
##' gg.change <- gg.change + geom_point() + geom_line() + facet_wrap(~treat.f)
##' gg.change
##' 
##' coef(eUN2.lmm)
##' effects(eUN2.lmm, type = "change", variable = "treat.f")
##' effects(eUN2.lmm, effects = "difference", type = "change", variable = "treat.f")
##'
##' ## fitted auc (conditional on covariates and covariates)
##' armd.AUC <- fitted(eUN2.lmm, type = "auc", keep.data = TRUE)
##'
##' gg.auc <- ggplot(armd.AUC, aes(x = treat.f, y = visual, color = impute))
##' gg.auc <- gg.auc + geom_point()
##' gg.auc
##'
##' aggregate(visual ~ treat.f, data = armd.AUC, FUN = "mean")
##' effects(eUN2.lmm, type = "auc", variable = "treat.f") ## adjusted for lesion
##' effects(eUN2.lmm, effects = "difference", type = "auc", variable = "treat.f")
##' }}

## * fitted.lmm (code)
##' @export
fitted.lmm <- function(object, newdata = NULL, type = "mean", se = NULL, df = NULL,
                       keep.data = NULL, format = "long", seed = NULL, simplify = TRUE, ...){

    mycall <- match.call()

    ## ** extract from object
    outcome.var <- object$outcome$var
    cluster.var <- object$cluster$var
    time.var <- object$time$var
    object.data <- object$data ## with columns XXindexXX, XXclusterXX, ... and all observations (including NAs)
    object.data.original <- object$data.original 

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    
    ## *** type
    type <- match.arg(type, c("mean","outcome","impute","change","auc","auc-b"))
    type.prediction <- switch(type,
                              mean = "static",
                              outcome = "dynamic",
                              impute = "dynamic",
                              change = "change",
                              auc = "auc",
                              auc = "auc-b")

    ## *** newdata
    if(identical(newdata,"unique")){        
        if(format == "wide"){
            stop("Argument \'newdata\' equals \"unique\" is only valid when argument \'format\' equals \"long\". \n")
        }
        if(type != "mean"){
            stop("Argument \'newdata\' equals \"unique\" is only valid when argument \'type\' equals \"mean\". \n")
        }
        if(is.null(keep.data)){
            keep.data <- TRUE
        }
        if(is.null(se)){
            se <- FALSE
        }
        newdata <- stats::model.frame(object, type = c("unique","mean"), options = options)
    }else if(is.null(newdata)){
        newdata <- object$data.original
        index.na <- object$index.na
        newdata.index.cluster <- attr(object$design$index.cluster,"vectorwise")
        newdata.index.time <- attr(object$design$index.clusterTime,"vectorwise")                    
    }else if(format == "wide"){
        index.na <- NULL
        newdata.design <- stats::model.matrix(object, newdata = newdata, effects = "index", na.rm = FALSE, options = options)
        newdata.index.cluster <- attr(newdata.design$index.cluster, "vectorwise")
        newdata.index.time <- attr(newdata.design$index.clusterTime, "vectorwise")        
    }

    ## *** keep.data
    if(is.null(keep.data)){
        keep.data <- FALSE
    }
    index.NA <- which(is.na(newdata[[outcome.var]]))
    n.NA <- length(index.NA)

    if(n.NA==0 && type %in% c("outcome","impute")){
        message("No fit to return since there are no missing outcome values")
        if(format == "long"){
            if(keep.data){
                return(newdata)
            }else if(simplify){
                return(newdata[[outcome.var]])
            }else{
                return(newdata[outcome.var])
            }
        }else if(format == "wide"){
            M.pred <- cbind(newdata[[outcome.var]])
            colnames(M.pred) <- outcome.var
            out <- .reformat(M.pred, name = outcome.var, format = format, simplify = simplify,
                             keep.data = keep.data, data = newdata, index.na = index.na,
                             object.cluster = object$cluster, index.cluster = newdata.index.cluster,
                             object.time = object$time, index.time = newdata.index.time,                     
                             call = mycall)
            return(out)
        }
    }

    ## *** se and df
    if(type == "impute"){
        if(is.null(se)){
            se <- object$args$df
        }else if(all(se == FALSE)){
            stop("Argument \'se\' should not be FALSE when argument \'type\' equals to \"impute\". \n",
                 "The standard error of the fitted outcome is needed to perform imputation. \n")
        }
        if(is.null(df)){
            df <- FALSE
        }
    }else if(is.null(se)){
        se <- (!simplify || keep.data) && !(format == "wide")
        if(is.null(df)){
            df <- object$args$df & se[1]
        }
    }else if(is.null(df)){
        df <- object$args$df & se[1]
    }

    ## *** format
    if(format == "wide" && sum(se)>0 && type != "impute"){
        stop("Cannot export standard errors, confidence intervals, and p-values in the wide format. \n",
             "Consider setting the argument \'format\' to \"long\". \n")
    }
    if("imputed" %in% names(newdata)){
        stop("Argument \'newdata\' should not contain a column called \"imputed\". \n")
    }
    format[] <- match.arg(format, c("wide","long")) ## use 'format[] <-' instead of 'format <-' to keep the name that will be transferd to .reformat(

    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1 or FALSE/0. \n")
    }

    ## ** extract fitted values
    e.pred <- do.call(stats::predict,
                      args = c(list(object, newdata = newdata, se = se, type = type.prediction, df = df, keep.data = TRUE, format = "long", simplify = FALSE, options = options), dots))

    ## ** evaluate the different types
    if(type == "mean"){
        e.fit <- e.pred[setdiff(names(e.pred),outcome.var)]
        names(e.fit)[names(e.fit) == "estimate"] <- outcome.var
    }else if(type == "outcome"){
        e.fit <- cbind(e.pred[,setdiff(names(e.pred),"estimate")],
                       impute = FALSE)
        
        e.fit[index.NA,outcome.var] <- e.pred$estimate[index.NA]
        e.fit[index.NA,"impute"] <- TRUE
    }else if(type == "change"){
        if(sum(se)==0){
            e.fit <- cbind(e.pred[,setdiff(names(e.pred),"estimate")],
                           impute = attr(e.pred,"grad"))
        }else{
            e.fit <- cbind(e.pred[,setdiff(names(e.pred),"estimate")],
                           impute = !is.na(e.pred$se))
        }
        e.fit[[outcome.var]] <- e.pred$estimate
    }else if(type %in% c("auc","auc-b")){
        index.keep <- which(!duplicated(e.pred[[cluster.var]]))
        e.fit <- cbind(e.pred[index.keep,setdiff(names(e.pred),"estimate")],
                       impute = FALSE)

        e.fit[[outcome.var]] <- e.pred[index.keep,"estimate"]
        e.fit[e.fit[[cluster.var]] %in% e.pred[[cluster.var]][index.NA],"impute"] <- TRUE
    }else if(type == "impute"){
        e.fit <- cbind(e.pred[,setdiff(names(e.pred),c("estimate","se","df","lower","upper"))],
                       impute = FALSE)
        index.NA <- which(is.na(newdata[[outcome.var]]))
        n.NA <- length(index.NA)

        
        ## set seed for reproducibility
        if(n.NA>0 && !is.null(seed)){
            if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
                old <- .Random.seed # to save the current seed
                on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
            }else{
                on.exit(rm(.Random.seed, envir=.GlobalEnv))
            }
            set.seed(seed)        
        }

        e.fit[index.NA,outcome.var] <- stats::rnorm(n.NA, mean = e.pred[index.NA,"estimate"], sd = e.pred[index.NA,"se"])
        e.fit[index.NA,"impute"] <- TRUE
    }

    ## ** reformat
    if(format == "long"){
        if(keep.data){
            out <- e.fit
        }else{
            out <- e.fit[,intersect(names(e.fit),c(outcome.var,"se","df","lower","upper")),drop=FALSE]
        }
    }else if(format == "wide"){
        M.pred <- cbind(e.fit[[outcome.var]])
        colnames(M.pred) <- outcome.var

        out <- .reformat(M.pred, name = names(format), format = format, simplify = simplify,
                         keep.data = keep.data, data = newdata, index.na = index.na,
                         object.cluster = object$cluster, index.cluster = newdata.index.cluster,
                         object.time = object$time, index.time = newdata.index.time,                     
                         call = mycall)
    }

    
    if(simplify && (sum(se)==0 | type == "impute") && keep.data == FALSE && format == "long"){
        out <- out[[outcome.var]]
    }else{
        attr(out,"vcov") <- attr(e.pred,"vcov")
        if(simplify==FALSE){
            attr(out,"grad") <- attr(e.pred,"grad")
            attr(out,"args") <- c(attr(e.pred,"grad"),list(seed = seed))
            if(type=="impute"){
                attr(out,"args")$seed <- seed
            }
        }
    }

    ## ** export
    return(out)
    
}


## * fitted.mlmm (documentation)
##' @title Fitted Outcome Value by Multiple Linear Mixed Models.
##' @description Evaluate the expected mean conditional to covariates
##' or the expected missing outcome values conditional to observed outcome values and covariates.
##' based on group-specific linear mixed models.
##'
##' @param object a \code{mlmm} object.
##' @param newdata [data.frame] the covariate values for each cluster along with the \code{by} variable indicating which linear mixed model should be used.
##' @param simplify [logical] simplify the column names (no mention of the time variable) and remove the cluster, by, and time column when possible.
##' @param ... additional argument passed the \code{\link{predict.lmm}}.
##' 

## * fitted.mlmm (code)
##' @export
fitted.mlmm <- function(object, newdata = NULL, simplify = TRUE, ...){

    ## ** normalize user input

    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    
    ## *** newdata
    if(!is.null(newdata)){
        if(identical(newdata,"unique")){
            ls.newdata <- stats::setNames(lapply(names(object$model), function(iM){"unique"}), names(object$model))
        }else{
            by <- object$object$by
            if("by" %in% names(newdata)){
                stop("The by variable (here \"",by,"\") could not be found in argument \'newdata\'. \n")
            }else if(any(newdata[[by]] %in% names(object$model) == FALSE)){
                stop("Incorrect value for the by variable (here \"",by,"\"): \"",paste(setdiff(newdata[[by]], names(object$model)), collapse = "\", \""),"\". \n",
                     "Valid values: \"",paste(names(object$model), collapse = "\", \""),"\". \n")
            }
            ls.newdata <- split(newdata, newdata[[by]])
        }
    }else{
        ls.newdata <- stats::setNames(vector(mode = "list", length = length(object$model)), names(object$model))
    }

    ## ** extract
    ls.out <- lapply(names(ls.newdata), function(iBy){ ## iBy <- "A"
        do.call(stats::fitted, args = c(list(object$model[[iBy]], newdata = ls.newdata[[iBy]], options = options, simplify = simplify), dots))
    })

    ## ** reshape
    test.2D <- any(sapply(ls.out, inherits, "data.frame")) || any(sapply(ls.out, inherits, "matrix"))
    if(test.2D){
        out <- do.call(rbind,ls.out)
        if(simplify && is.data.frame(out)){
            object.manifest <- stats::variable.names(object)
            rm.manifest <- attributes(object.manifest)[setdiff(names(attributes(object.manifest)),c("by","cluster","time"))]
            out[unique(unlist(rm.manifest))] <- NULL
        }
    }else{
        out <- do.call(c,ls.out)
    }

    ## ** export
    return(out)

}

##----------------------------------------------------------------------
### fitted.R ends here
