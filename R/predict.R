### predict.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: May 31 2021 (19:28) 
##           By: Brice Ozenne
##     Update #: 111
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
## * predict.lmm (documentation)
##' @title Predicted Mean Value For Linear Mixed Models
##'
##' @param object a \code{lmm} object.
##' @param newdata [data.frame] the covariate values for each cluster.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the predicted mean. Otherwise a normal distribution is used.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A data.frame with 5 columns:\itemize{
##' \item \code{estimate}: predicted mean.
##' \item \code{se}: uncertainty about the predicted mean.
##' \item \code{df}: degree of freedom
##' \item \code{lower}: lower bound of the confidence interval of the predicted mean
##' \item \code{upper}: upper bound of the confidence interval of the predicted mean
##' }
##'
##' @examples
##' ## simulate data in the long format
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' ## fit mixed model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' ## prediction
##' predict(eUN.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3))
##'
##' ## with Student's t-distribution
##' \dontrun{
##' predict(eUN.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3), df = TRUE)
##' }

## * predict.lmm (code)
##' @export
predict.lmm <- function(object, newdata, level = 0.95, df = !is.null(object$df), ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** compute predictions
    ff.mean <- stats::formula(object, effects = "mean")
    X.beta <- model.matrix(delete.response(terms(ff.mean)), newdata)
    n.obs <- NROW(X.beta)
    beta <- coef(object, effects = "mean")
    n.beta <- length(beta)
    name.beta <- names(beta)
    vcov.beta <- vcov(object, effects = "mean")
    
    prediction <- X.beta %*% beta
    prediction.vcov <- X.beta %*% vcov.beta %*% t(X.beta)
    prediction.se <- sqrt(diag(prediction.vcov))

    ## ** df
    if(df){
        vcov.param <- vcov(object, effects = "all", df = 2, transform.names = FALSE)
        dVcov <- attr(vcov.param,"dVcov")
        attr(vcov.param, "df") <- NULL
        attr(vcov.param, "dVcov") <- NULL
        prediction.df <- .dfX(X = X.beta, vcov.param = vcov.param, dVcov.param = dVcov)
    }else{
        prediction.df <- rep(Inf,n.obs)
    }


    ## ** export
    alpha <- 1-level
    out <- data.frame(estimate = prediction, se = prediction.se, df = prediction.df,
                      lower = prediction + stats::qt(alpha/2, df = prediction.df) *prediction.se,
                      upper = prediction + stats::qt(1-alpha/2, df = prediction.df)*prediction.se)
    return(out)
    
}

## * .dfX
.dfX <- function(X.beta, vcov.param, dVcov.param){
    
    if(!is.matrix(X.beta) && is.vector(X.beta)){
        X.beta <- rbind(X.beta)
    }
    n.obs <- NROW(X.beta)
    name.beta <- colnames(X.beta)
    n.beta <- length(name.beta)

    vcov.beta <- vcov.param[name.beta,name.beta,drop=FALSE]
    n.param <- NCOL(vcov.param)
    name.param <- colnames(vcov.param)
    prediction.vcov <- X.beta %*% vcov.beta %*% t(X.beta)
    
    ## ** variance covariance matrix of the variance of the mean parameters
    ## delta method:
    ## Cov[Sigma_\beta,Sigma_\beta'] = [dSigma_\beta/d\beta] [Sigma_\theta] [dSigma_\beta'/d\beta']
    n.beta2 <- n.beta^2
    Mname.beta2 <- expand.grid(name.beta,name.beta)
    name.beta2 <- interaction(Mname.beta2)
    
    vcov.vcovbeta <- matrix(0, nrow = n.beta2, ncol = n.beta2, dimnames = list(name.beta2,name.beta2))
    index.n0 <- which(apply(Mname.beta2,1,function(iP){any(abs(dVcov.param[iP[1],iP[2],])>1e-10)}))
    
    for(iP in index.n0){ ## iP <- index.n0[2]
        iParam <- name.beta2[iP]
        iProd <- as.double(dVcov.param[Mname.beta2[iP,1],Mname.beta2[iP,2],] %*% vcov.param)
        iIndex.n0 <- which(abs(iProd)>1e-10)

        for(iiP in index.n0[index.n0<=iP]){ ## iiP <- 2
            iiParam <-  name.beta2[iiP]
            vcov.vcovbeta[iParam,iiParam] <- sum(iProd[iIndex.n0] * dVcov.param[Mname.beta2[iiP,1],Mname.beta2[iiP,2],iIndex.n0])
            ## if(iParam != iiParam){
            ##     vcov.vcovbeta[iiParam,iParam] <- vcov.vcovbeta[iParam,iiParam]
            ## }
        }
    }
    vcov.vcovbeta[upper.tri(vcov.vcovbeta)]  <-  t(vcov.vcovbeta)[upper.tri(vcov.vcovbeta)]
    
    ## GS <- matrix(0, nrow = n.beta2, ncol = n.beta2, dimnames = list(name.beta2,name.beta2))
    ## for(iP in 1:n.beta2){ ## iP <- 1
    ##     for(iiP in 1:iP){ ## iiP <- 1
    ##         iParam <-  name.beta2[iP]
    ##         iiParam <-  name.beta2[iiP]
    ##         GS[iParam,iiParam] <- dVcov.param[Mname.beta2[iP,1],Mname.beta2[iP,2],] %*% vcov.param %*% dVcov.param[Mname.beta2[iiP,1],Mname.beta2[iiP,2],]
    ##         if(iParam != iiParam){
    ##             GS[iiParam,iParam] <- vcov.vcovbeta[iParam,iiParam]
    ##         }
    ##     }
    ## }

    ## ** gradient of the variance of the predictions relative to the variance of the mean parameters
    ## d(X\Sigma t(X))_ij is computed by X\delta_ij t(X)
    dXTX.dT <- matrix(NA, nrow = n.obs, ncol = n.beta2, dimnames = list(NULL,name.beta2))
        
    for(iP in 1:n.beta2){
        iDelta <- matrix(0, nrow = n.beta, ncol = n.beta, dimnames = list(name.beta,name.beta))
        iDelta[Mname.beta2[iP,1],Mname.beta2[iP,2]] <- 1
        dXTX.dT[,iP] <- diag(X.beta %*% iDelta %*% t(X.beta))
    }

    ## ** Satterthwaite approximation of the degrees of freedom
    ## NOTE: df(beta) is 2*diag(vcov.beta)^2 / diag(vcov.vcovbeta[Mname.beta2[,1]==Mname.beta2[,2],Mname.beta2[,1]==Mname.beta2[,2]])
    out <- sapply(1:n.obs, function(iObs){
        2 * prediction.vcov[iObs,iObs]^2 / (dXTX.dT[iObs,] %*% vcov.vcovbeta %*% dXTX.dT[iObs,])
    })


    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### predict.R ends here
