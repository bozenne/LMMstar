### predict.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: May 19 2021 (12:36) 
##           By: Brice Ozenne
##     Update #: 62
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
## * predict.lmm (code)
##' @export
predict.lmm <- function(object, newdata, level = 0.95, df = TRUE, ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** compute predictions
    X.beta <- model.matrix(object, effects = "mean")
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
        prediction.df <- .dfX(X = X.beta, vcov.param = vcov(object, effects = "all"), dVcov.param = object$dVcov)
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
    vcov.vcovbeta <- matrix(NA, nrow = n.beta2, ncol = n.beta2, dimnames = list(name.beta2,name.beta2))
    for(iP in 1:n.beta2){ ## iP <- 1
        for(iiP in 1:iP){ ## iiP <- 1
            iParam <-  name.beta2[iP]
            iiParam <-  name.beta2[iiP]
            vcov.vcovbeta[iParam,iiParam] <- dVcov.param[Mname.beta2[iP,1],Mname.beta2[iP,2],] %*% vcov.param %*% dVcov.param[Mname.beta2[iiP,1],Mname.beta2[iiP,2],]
            if(iParam != iiParam){
                vcov.vcovbeta[iiParam,iParam] <- vcov.vcovbeta[iParam,iiParam]
            }
        }
    }
    ## NOTE: df(beta) is 2*diag(vcov.beta)^2 / diag(vcov.vcovbeta[Mname.beta2[,1]==Mname.beta2[,2],Mname.beta2[,1]==Mname.beta2[,2]])
    
    ## ** gradient of the variance of the predictions relative to the variance of the mean parameters
    ## d(X\Sigma t(X))_ij is computed by X\delta_ij t(X)
    dXTX.dT <- matrix(NA, nrow = n.obs, ncol = n.beta2, dimnames = list(NULL,name.beta2))
        
    for(iP in 1:n.beta2){
        iDelta <- matrix(0, nrow = n.beta, ncol = n.beta, dimnames = list(name.beta,name.beta))
        iDelta[Mname.beta2[iP,1],Mname.beta2[iP,2]] <- 1
        dXTX.dT[,iP] <- diag(X.beta %*% iDelta %*% t(X.beta))
    }

    ## ** Satterthwaite approximation of the degrees of freedom
    out <- sapply(1:n.obs, function(iObs){
        2 * prediction.vcov[iObs,iObs]^2 / (dXTX.dT[iObs,] %*% vcov.vcovbeta %*% dXTX.dT[iObs,])
    })


    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### predict.R ends here
