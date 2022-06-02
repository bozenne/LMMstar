### fitted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  8 2021 (17:09) 
## Version: 
## Last-Updated: Jun  2 2022 (11:43) 
##           By: Brice Ozenne
##     Update #: 125
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

    ## ** extract from object
    outcome.var <- object$outcome$var
    cluster.var <- object$cluster$var
    time.var <- object$time$var
    object.fitted <- object$fitted
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
    if(impute && "imputed" %in% names(newdata)){
        stop("Argument \'newdata\' should not contain a column called \"imputed\". \n")
    }
    if(!impute){
        se.impute <- FALSE
    }
    type.prediction <- ifelse(impute,"dynamic", "static")

    format <- match.arg(format, c("wide","long"))
    
    test.original.data <- is.null(newdata)
    if(test.original.data){
        newdata <- object.data.original
    }

    if(format == "wide") {
        keep.newdata <- TRUE
    }

    ## ** compute predictions
    e.pred <- stats::predict(object,
                             newdata = newdata,
                             type = type.prediction,
                             se = se.impute,
                             keep.newdata = keep.newdata)
    
    ## ** store
    if(impute){ ## store dynamic predictions

            index.NA <- which(!is.na(e.pred$estimate))
            if(keep.newdata == FALSE){
                if("se" %in% names(e.pred) == FALSE){
                    out <- e.pred[index.NA,"estimate"]
                }else{
                    out <- stats::rnorm(length(index.NA),
                                        mean = e.pred[index.NA,"estimate"],
                                        sd = e.pred[index.NA,"se"])
                }
                
            }else if(length(index.NA) > 0){
                if("se" %in% names(e.pred) == FALSE){
                    newdata[index.NA,outcome.var] <- e.pred[index.NA,"estimate"]
                }else{
                    newdata[index.NA,outcome.var] <- stats::rnorm(length(index.NA),
                                                                  mean = e.pred[index.NA,"estimate"],
                                                                  sd = e.pred[index.NA,"se"])
                }
                newdata$imputed <- FALSE
                newdata$imputed[index.NA] <- TRUE
                out <- newdata
            }
            value.var <- c(outcome.var,"imputed")
    }else{ ## store static predictions

        if(keep.newdata == FALSE){
            out <- e.pred[["estimate"]]
        }else{
            out <- e.pred
        }
        value.var <- "estimate"

    }

    ## ** convert to wide format
    if(format=="wide"){
        out <- stats::reshape(data = out[,c(attr(time.var,"original"), attr(cluster.var,"original"), value.var),drop=FALSE], 
                              direction = "wide", timevar = attr(time.var,"original"), idvar = attr(cluster.var,"original"), v.names = value.var)
    }

    ## ** export
    return(out)
    
}



##----------------------------------------------------------------------
### fitted.R ends here
