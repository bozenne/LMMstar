### predict.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: Jul  9 2021 (10:05) 
##           By: Brice Ozenne
##     Update #: 152
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
## * predict.lmm (documentation)
##' @title Predicted Mean Value With Uncertainty For Linear Mixed Model
##'
##' @param object a \code{lmm} object.
##' @param newdata [data.frame] the covariate values for each cluster.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the predicted mean. Otherwise a normal distribution is used.
##' @param keep.newdata [logical] Should the argument \code{newdata} be output along side the predicted values?
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
##' ## fit Linear Mixed Model
##' eUN.lmm <- lmm(Y ~ X1 + X2 + X5, repetition = ~visit|id, structure = "UN", data = dL, df = FALSE)
##'
##' ## prediction
##' predict(eUN.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3))
##' predict(eUN.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3), keep.newdata = TRUE)
##'
##' ## with Student's t-distribution
##' \dontrun{
##' predict(eUN.lmm, newdata = data.frame(X1 = 1, X2 = 2, X5 = 3), df = TRUE)
##' }

## * predict.lmm (code)
##' @export
predict.lmm <- function(object, newdata, level = 0.95, df = !is.null(object$df), keep.newdata = FALSE, ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    

    ## ** compute predictions
    ## extract coefficients
    beta <- coef(object, effects = "mean")
    n.beta <- length(beta)
    name.beta <- names(beta)
    vcov.beta <- vcov(object, effects = "mean")

    ## extract design matrix
    X.beta <- .predict_model.matrix(object, newdata = newdata, name.beta = name.beta)
    n.obs <- NROW(X.beta)
    
    prediction <- X.beta %*% beta
    prediction.vcov <- X.beta %*% vcov.beta %*% t(X.beta)
    prediction.se <- sqrt(diag(prediction.vcov))

    ## ** df
    if(df){
        vcov.param <- vcov(object, effects = "all", df = 2, transform.names = FALSE)
        dVcov <- attr(vcov.param,"dVcov")
        attr(vcov.param, "df") <- NULL
        attr(vcov.param, "dVcov") <- NULL
        prediction.df <- .dfX(X.beta = X.beta, vcov.param = vcov.param, dVcov.param = dVcov)
    }else{
        prediction.df <- rep(Inf,n.obs)
    }


    ## ** export
    alpha <- 1-level
    out <- data.frame(estimate = prediction, se = prediction.se, df = prediction.df,
                      lower = prediction + stats::qt(alpha/2, df = prediction.df) *prediction.se,
                      upper = prediction + stats::qt(1-alpha/2, df = prediction.df)*prediction.se)
    if(keep.newdata){
        return(cbind(newdata,out))
    }else{
        return(out)
    }
    
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
    
    Mpair_dVcov.param <- do.call(rbind,lapply(1:n.beta2, function(iIndex){dVcov.param[Mname.beta2[iIndex,1],Mname.beta2[iIndex,2],]})) ## iIndex <- 240
    vcov.vcovbeta <- Mpair_dVcov.param %*% vcov.param %*% t(Mpair_dVcov.param)
    rownames(vcov.vcovbeta)  <-  name.beta2
    colnames(vcov.vcovbeta)  <-  name.beta2
    
    ## GS <- matrix(0, nrow = n.beta2, ncol = n.beta2, dimnames = list(name.beta2,name.beta2))
    ## for(iP in 1:n.beta2){ ## iP <- 15
    ##     for(iiP in 1:iP){ ## iiP <- 29
    ##         iParam <-  name.beta2[iP]
    ##         iiParam <-  name.beta2[iiP]
    ##         GS[iParam,iiParam] <- dVcov.param[Mname.beta2[iP,1],Mname.beta2[iP,2],] %*% vcov.param %*% dVcov.param[Mname.beta2[iiP,1],Mname.beta2[iiP,2],]
    ##         if(iParam != iiParam){
    ##             GS[iiParam,iParam] <- vcov.vcovbeta[iParam,iiParam]
    ##         }
    ##     }
    ## }
    ## GS[iP,iiP] - vcov.vcovbeta[iP,iiP]

    ## ** gradient of the variance of the predictions relative to the variance of the mean parameters
    ## d(X\Sigma t(X))_ij is computed by X\delta_ij t(X)
    dXTX.dT <- matrix(NA, nrow = n.obs, ncol = n.beta2, dimnames = list(NULL,name.beta2))

    ## matrix cookbook: A_ki B_lj = (A \delta_ij \trans{B})_kl
    ## so               A_ki A_kj = (A \delta_ij \trans{A})_kk
    ## so               t(A)_ik A_kj = (A \delta_ij \trans{A})_kk
    
    for(iP in 1:n.beta2){ ## iP <- 35
        ## iDelta <- matrix(0, nrow = n.beta, ncol = n.beta, dimnames = list(name.beta,name.beta))
        ## iDelta[Mname.beta2[iP,1],Mname.beta2[iP,2]] <- 1
        ## dXTX.dT[,iP] <- diag(X.beta %*% iDelta %*% t(X.beta))
        dXTX.dT[,iP] <- (X.beta[,Mname.beta2[iP,1]] * X.beta[,Mname.beta2[iP,2]])
    }

    ## ** Satterthwaite approximation of the degrees of freedom
    ## NOTE: df(beta) is 2*diag(vcov.beta)^2 / diag(vcov.vcovbeta[Mname.beta2[,1]==Mname.beta2[,2],Mname.beta2[,1]==Mname.beta2[,2]])
    out <- sapply(1:n.obs, function(iObs){
        2 * prediction.vcov[iObs,iObs]^2 / (dXTX.dT[iObs,] %*% vcov.vcovbeta %*% dXTX.dT[iObs,])
    })


    ## ** export
    return(out)
}

## * .predict_model.matrix
.predict_model.matrix <- function(object, newdata, name.beta){

    if(is.matrix(newdata)){
        if(all(name.beta == colnames(newdata))){
            X.beta <- newdata
        }else{
            stop("Argument \'newdata\' should contain a column relative to each coefficient when a matrix. \n")
        }
    }else if(is.data.frame(newdata)){
        ## use model.frame to be able to specify the na.action behavior
        ## use [,name.beta,] to make sure to only keep relevant columns (drop the columns drop when constraining baseline)
        if(object$strata$n==1){
            ff.mean <- stats::formula(object, effects = "mean")
            ff.meanRHS <- stats::delete.response(stats::terms(ff.mean))
            X.beta <- stats::model.matrix(ff.meanRHS, stats::model.frame(ff.meanRHS, newdata, na.action = stats::na.pass))[,name.beta,drop=FALSE]
        }else{
            ff.mean <- object$formula$mean.design
            ff.meanRHS <- stats::delete.response(stats::terms(ff.mean))
            X.beta <- matrix(NA, nrow = NROW(newdata), ncol = length(name.beta), dimnames = list(NULL, name.beta))
            for(iStrata in object$strata$levels){
                iIndex <- which(newdata[[object$strata$var]]==iStrata)
                iX.beta <- stats::model.matrix(ff.meanRHS, stats::model.frame(ff.meanRHS, newdata[iIndex,,drop=FALSE], na.action = stats::na.pass))
                colnames(iX.beta) <- paste0(colnames(iX.beta),":",iStrata)
                X.beta[iIndex,] <- 0
                X.beta[iIndex,colnames(iX.beta)] <- iX.beta
            }
        }
        rownames(X.beta) <- NULL
    }else{
        stop("Argument \'newdata\' should be a data.frame or a design matrix. \n")
    }
    return(X.beta)
}
##----------------------------------------------------------------------
### predict.R ends here
