### fitted.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jul  8 2021 (17:09) 
## Version: 
## Last-Updated: sep 23 2021 (20:53) 
##           By: Brice Ozenne
##     Update #: 12
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
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A vector of length the number of row of newdata
##'
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' ## prediction
##' fitted(eUN.lmm)
##' fitted(eUN.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3))
##' fitted(eUN.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3), keep.newdata = TRUE)

## * fitted.lmm (code)
##' @export
fitted.lmm <- function(object, newdata = NULL, keep.newdata = FALSE, ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    

    ## ** compute predictions

    ## design matrix
    if(is.null(newdata)){
        if(length(object$index.na)>0){
            prediction <- rep(NA, NROW(object$data.original))
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
        prediction <- X.beta %*% beta
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
