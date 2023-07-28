### fitted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  8 2021 (17:09) 
## Version: 
## Last-Updated: jul 28 2023 (17:37) 
##           By: Brice Ozenne
##     Update #: 207
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * fitted.lmm (documentation)
##' @title Predicted Mean Value For Linear Mixed Model
##'
##' @param object a \code{lmm} object.
##' @param newdata [data.frame] the covariate values for each cluster.
##' @param impute [logical] Should the missing data in the outcome be imputed based on covariates and other outcome values from the same cluster.
##' @param se.impute [character] If \code{FALSE} the most likely value is imputed. Otherwise the imputed value is sampled from a normal distribution.
##' The value of the argument determine which standard deviation is used: all uncertainty about the predicted value (\code{"total"}),
##' only uncertainty related to the estimation of the model parameters (\code{"estimate"}), or only uncertainty related to the residual variance of the outcome (\code{"residual"}).
##' Passed to \code{predict.lmm}.
##' @param format [character] Should the prediction be output
##' in a matrix format with clusters in row and timepoints in columns (\code{"wide"}),
##' or in a data.frame/vector with as many rows as observations (\code{"long"})
##' @param keep.newdata [logical] Should the dataset relative to which the predictions are evaluated be output along side the predicted values?
##' Only possible in the long format.
##' @param simplify [logical] Simplify the data format (vector instead of data.frame) and column names (no mention of the time variable) when possible.
##' @param seed [integer, >0] Random number generator (RNG) state used when starting imputation.
##' If \code{NULL} no state is set. 
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return When \code{format="wide"}, a data.frame with as many rows as clusters.
##' When \code{format="long"} or \code{keep.newdata==TRUE}, a data.frame with as many rows as observations.
##' Otherwise: \itemize{
##' \item if \code{impute=FALSE} a vector of length the number of row of newdata containing the fitted values (i.e. based on the covariates only).
##' \item if \code{impute=TRUE} a vector of length the number of missing values in the outcome of newdata containing the cluster-specific conditional means
##' (i.e. based on the covariates and outcome measurements from the same cluster).
##' }
##' 
##' When \code{keep.newdata==TRUE}, a dataframe with an additional column containing the fitted values (i.e. based on the covariates only).
##' If \code{impute=TRUE}, the missing value in the outcome column are replaced by the cluster-specific conditional means
##' (i.e. based on the covariates and outcome measurements from the same cluster).
##' 
##' @keywords methods
##' 
##' @examples
##' #### simulate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' #### fit Linear Mixed Model ####
##' eCS.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id,
##'                structure = "CS", data = dL, df = FALSE)
##'
##' ## prediction
##' fitted(eCS.lmm)
##' fitted(eCS.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3))
##' fitted(eCS.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3), keep.newdata = TRUE)
##' 
##' #### fit Linear Mixed Model with missing data ####
##' dL2 <- dL
##' dL2[3,"Y"] <- NA
##' eCS2.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id,
##'                 structure = "CS", data = dL2, df = FALSE)
##' 
##' ## most likely value to impute
##' fitted(eCS2.lmm, impute = TRUE)
##' head(fitted(eCS2.lmm, impute = TRUE, keep.newdata = TRUE))
##' 
##' ## multiple imputation
##' dL2.imp1 <- data.frame(imp = "1",
##'     fitted(eCS2.lmm, impute = TRUE, se.impute = "total", keep.newdata = TRUE))
##' dL2.imp2 <- data.frame(imp = "2",
##'     fitted(eCS2.lmm, impute = TRUE, se.impute = "total", keep.newdata = TRUE))
##' head(dL2.imp1)
##' head(dL2.imp2)

## * fitted.lmm (code)
##' @export
fitted.lmm <- function(object, newdata = NULL, impute = FALSE, se.impute = FALSE, 
                       keep.newdata = FALSE, format = "long", simplify = TRUE, seed = NULL, ...){

    ## ** special case: no imputation
    if(impute == FALSE){
        return(stats::predict(object,
                              newdata = newdata,
                              se = FALSE,
                              type = "static",
                              df = FALSE,
                              keep.newdata = keep.newdata,
                              format = format,
                              simplify = simplify,
                              ...))
    }

    ## ** extract from object
    outcome.var <- object$outcome$var
    cluster.var <- object$cluster$var
    time.var <- object$time$var
    object.data <- object$data ## with columns XXindexXX, XXclusterXX, ... and all observations (including NAs)
    object.data.original <- object$data.original 

    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.logical(impute)){
        stop("Argument \'impute\' should be TRUE or FALSE. \n")
    }
    if("imputed" %in% names(newdata)){
        stop("Argument \'newdata\' should not contain a column called \"imputed\". \n")
    }
    type.prediction <- ifelse(impute,"dynamic", "static")

    format[] <- match.arg(format, c("wide","long")) ## use 'format[] <-' instead of 'format <-' to keep the name that will be transferd to .reformat(

    if(keep.newdata && format == "wide"){
        stop("Argument \'keep.newdata\' must be FALSE when using the wide format. \n")
    }    
    
    if(is.null(newdata)){
        newdata <- object$data.original
    }

    ## ** compute predictions
    e.pred <- stats::predict(object,
                             newdata = newdata,
                             se = se.impute,
                             type = type.prediction,
                             df = FALSE,
                             keep.newdata = TRUE,
                             format = "long",
                             simplify = FALSE)
    se.impute <- attr(e.pred,"args")$se
    
    ## ** update dynamic predictions
    index.NA <- which(is.na(newdata[[outcome.var]]))
    n.NA <- length(index.NA)

    ## set seed for reproducibility
    if(n.NA>0 && !is.null(se.impute) && !is.null(seed)){
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
        set.seed(seed)        
    }
    
    if(keep.newdata == FALSE){
        out <- newdata[[outcome.var]]
        if(n.NA>0){
            if(is.null(se.impute)){
                out[index.NA] <- e.pred[index.NA,"estimate"]
            }else{
                out[index.NA] <- stats::rnorm(n.NA, mean = e.pred[index.NA,"estimate"], sd = e.pred[index.NA,"se"])
            }
        }
        if(format == "long"){
            if(simplify){
                return(out)
            }else{
                out <- data.frame(out, impute = FALSE)
                names(out)[1] <- outcome.var
                out$impute[index.NA] <- TRUE
                return(out)
            }
        }
    }else if(keep.newdata){ ## must be the long format
        if(n.NA>0){
            if(is.null(se.impute)){
                newdata[index.NA,outcome.var] <- e.pred[index.NA,"estimate"]
            }else{
                newdata[index.NA,outcome.var] <- stats::rnorm(n.NA, mean = e.pred[index.NA,"estimate"], sd = e.pred[index.NA,"se"])
            }
            newdata$imputed <- FALSE
            newdata$imputed[index.NA] <- TRUE
        }else{
            newdata$imputed <- FALSE
        }
        return(newdata)
    }

    ## ** convert to wide format
    ## only case as otherwise the long format is dealt with before: keep.newdata == FALSE or keep.newdata == TRUE
    newdata.design <- model.matrix(object, data = newdata, effects = "index")
    newdata.index.cluster <- attr(newdata.design$index.cluster, "vectorwise")
    newdata.index.time <- attr(newdata.design$index.clusterTime, "vectorwise")
    Mout <- cbind(out)
    colnames(Mout) <- outcome.var

    out <- .reformat(Mout, name = names(format), format = format, simplify = simplify,
                     keep.data = keep.newdata, data = newdata, index.na = NULL,
                     object.cluster = object$cluster, index.cluster = newdata.index.cluster,
                     object.time = object$time, index.time = newdata.index.time,
                     call = match.call())


    ## ** export
    return(out)
    
}


## * fitted.mlmm (code)
##' @export
fitted.mlmm <- function(object, ...){
    stop("No \'fitted\' method for mlmm objects, consider using lmm instead of mlmm. \n")
}
##----------------------------------------------------------------------
### fitted.R ends here
