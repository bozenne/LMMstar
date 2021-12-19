### fitted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  8 2021 (17:09) 
## Version: 
## Last-Updated: Dec 15 2021 (14:05) 
##           By: Brice Ozenne
##     Update #: 56
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
##' @param keep.newdata [logical] Should the argument \code{newdata} be output along side the predicted values?
##' @param impute [logical] Should the missing data in the outcome be imputed based on covariates and other outcome values from the same cluster.
##' @param se.impute [character] If \code{FALSE} the most likely value is imputed. Otherwise the imputed value is sampled from a normal distribution.
##' The value of the argument determine which standard deviation is used: all uncertainty about the predicted value (\code{"total"}),
##' only uncertainty related to the estimation of the model parameters (\code{"estimate"}), or only uncertainty related to the residual variance of the outcome (\code{"residual"}).
##' Passed to \code{predict.lmm}.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return When \code{keep.newdata==FALSE}: \itemize{
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
fitted.lmm <- function(object, newdata = NULL, keep.newdata = FALSE, impute = FALSE, se.impute = FALSE, ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.logical(impute)){
        stop("Argument \'impute\' should be TRUE or FALSE. \n")
    }

    ## ** impute missing values
    if(impute){
        if(is.null(newdata)){
            newdata <- object$data.original
        }
        e.impute <- stats::predict(object, newdata = newdata, type = "dynamic", se = se.impute, keep.newdata = TRUE)
        index.NA <- which(!is.na(e.impute$estimate))
        if(length(index.NA) > 0 && "se" %in% names(e.impute) == FALSE){
            newdata[index.NA,object$outcome$var] <- e.impute[index.NA,"estimate"]
        }else{
            newdata[index.NA,object$outcome$var] <- stats::rnorm(length(index.NA), mean = e.impute[index.NA,"estimate"], sd = e.impute[index.NA,"se"])
        }

        if(keep.newdata==FALSE){
            return(newdata[index.NA,object$outcome$var])
        }else{
            if("imputed" %in% names(newdata)){
                stop("Argument \'newdata\' should not contain a column called \"imputed\". \n")
            }
            newdata$imputed <- FALSE
            newdata$imputed[index.NA] <- TRUE
        }
    }

    ## ** compute predictions
    ## design matrix
    if(is.null(newdata)){
        newdata <- object$data.original
            
        if(length(object$index.na)>0){
            prediction <- rep(NA, NROW(newdata))
            prediction[-object$index.na] <- object$fitted
        }else{
            prediction <- as.vector(object$fitted)
        }
    }else{
        ## extract coefficients
        beta <- coef(object, effects = "mean")
        name.beta <- names(beta)

        ## generate design matrix
        X.beta <- model.matrix(object, data = newdata, effects = "mean")

        ## compute predictions
        prediction <- as.vector(X.beta %*% beta)
    }


    ## ** export
    if(keep.newdata){
        return(cbind(newdata, estimate = prediction))
    }else{
        return(prediction)
    }
    
}



##----------------------------------------------------------------------
### fitted.R ends here
