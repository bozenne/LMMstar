### fitted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  8 2021 (17:09) 
## Version: 
## Last-Updated: May 29 2022 (13:45) 
##           By: Brice Ozenne
##     Update #: 84
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
##' @param format [character] Should the predicted mean be output relative as a vector (\code{"long"}), or as a matrix with in row the clusters and in columns the outcomes (\code{"wide"}).
##' @param keep.newdata [logical] Should the argument \code{newdata} be output along side the predicted values? The output will then be a \code{data.frame}.
##' @param impute [logical] Should the missing data in the outcome be imputed based on covariates and other outcome values from the same cluster.
##' @param se.impute [character] If \code{FALSE} the most likely value is imputed. Otherwise the imputed value is sampled from a normal distribution.
##' The value of the argument determine which standard deviation is used: all uncertainty about the predicted value (\code{"total"}),
##' only uncertainty related to the estimation of the model parameters (\code{"estimate"}), or only uncertainty related to the residual variance of the outcome (\code{"residual"}).
##' Passed to \code{predict.lmm}.
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
fitted.lmm <- function(object, newdata = NULL, format = "long",
                       keep.newdata = FALSE, impute = FALSE, se.impute = FALSE, ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.logical(impute)){
        stop("Argument \'impute\' should be TRUE or FALSE. \n")
    }
    format <- match.arg(format, c("wide","long"))
    outcome.var <- object$outcome$var
    cluster.var <- object$cluster$var
    time.var <- object$time$var

    test.original.data <- is.null(newdata)
    if(test.original.data){
        newdata <- object$data.original
    }
        

    ## ** prepare for wide format transformation
    if(format == "wide"){
        keep.newdata <- TRUE
        if((cluster.var %in% names(newdata) == FALSE) || (time.var %in% names(newdata) == FALSE)){
            data2prepare <- as.data.frame(newdata)
            data2prepare$XXindexXX <- NULL
            data2prepare$XXclusterXX <- NULL
            data2prepare$XXcluster.indexXX <- NULL
            data2prepare$XXtimeXX <- NULL
            data2prepare$XXtime.indexXX <- NULL
            data2prepare$XXstrataXX <- NULL
            data2prepare$XXstrata.indexXX <- NULL
            newdata <- .prepareData(data = data2prepare, var.cluster = attr(cluster.var,"original"), var.time = attr(time.var,"original"), var.strata = NA, missing.repetition = NULL)
        }
    }
    
    ## ** impute missing values
    if(impute){
        e.impute <- stats::predict(object, newdata = newdata, type = "dynamic", se = se.impute, keep.newdata = TRUE)
        index.NA <- which(!is.na(e.impute$estimate))
        if(length(index.NA) > 0 && "se" %in% names(e.impute) == FALSE){
            newdata[index.NA,outcome.var] <- e.impute[index.NA,"estimate"]
        }else{
            newdata[index.NA,outcome.var] <- stats::rnorm(length(index.NA), mean = e.impute[index.NA,"estimate"], sd = e.impute[index.NA,"se"])
        }

        if(keep.newdata==FALSE){
            out <- newdata[index.NA,outcome.var]
        }else{
            if("imputed" %in% names(newdata)){
                stop("Argument \'newdata\' should not contain a column called \"imputed\". \n")
            }
            newdata$imputed <- FALSE
            newdata$imputed[index.NA] <- TRUE
            out <- newdata
        }
        value.var <- outcome.var
    }

    ## ** compute predictions
    ## design matrix
    if(!impute){
        if(test.original.data){
            if(length(object$index.na)>0){
                out <- rep(NA, NROW(newdata))
                out[-object$index.na] <- object$fitted
            }else{
                out <- as.vector(object$fitted)
            }
        }else{
            ## extract coefficients
            beta <- stats::coef(object, effects = "mean")
            name.beta <- names(beta)

            ## generate design matrix
            X.beta <- stats::model.matrix(object, data = newdata, effects = "mean")

            ## compute predictions
            out <- as.vector(X.beta %*% beta)
        }
        if(keep.newdata){
            if("estimate" %in% names(newdata)){
                stop("Argument \'newdata\' should not contain a column called \"estimate\". \n")
            }
            out <- cbind(newdata, estimate = out)
        }
        value.var <- "estimate"
    }

    ## ** convert to wide format
    if(format=="wide"){
        if(value.var!="estimate"){
            out[[value.var]] <- out[[value.var]]*c(NA,1)[1+out$imputed]
        }
        out <- reshape2::dcast(data = out,
                               formula = stats::as.formula(paste0(object$cluster$var,"~",object$time$var)), value.var = value.var)
    }

    ## ** export
    return(out)
    
}



##----------------------------------------------------------------------
### fitted.R ends here
