### predict.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: okt  2 2021 (17:32) 
##           By: Brice Ozenne
##     Update #: 485
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
##' @param se [character] Type of uncertainty to be accounted for: estimation of the regression parameters (\code{"estimation"}), residual variance (\code{"residual"}), or both (\code{"total"}).
##' Can also be \code{NULL} to not compute standard error, p-values, and confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the predicted mean. Otherwise a normal distribution is used.
##' @param type [character] Should prediction be made conditional on the covariates only (\code{"static"}) or also on outcome values at other timepoints (\code{"dynamic"}).
##' @param keep.newdata [logical] Should the argument \code{newdata} be output along side the predicted values?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Static prediction are made using the linear predictor \eqn{X\beta} while dynamic prediction uses the conditional normal distribution of the missing outcome given the observed outcomes. So if outcome 1 is observed but not 2, prediction for outcome 2 is obtain by \eqn{X_2\beta + \sigma_{21}\sigma^{-1}_{22}(Y_1-X_1\beta)}. In that case, the uncertainty is computed as the sum of the conditional variance \eqn{\sigma_{22}-\sigma_{21}\sigma^{-1}_{22}\sigma_{12}} plus the uncertainty about the estimated conditional mean (obtained via delta method using numerical derivatives).
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
##' eUN.lmm <- lmm(Y ~ visit + X1 + X2 + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' ## prediction
##' newd <- data.frame(X1 = 1, X2 = 2, X5 = 3, visit = factor(1:3, levels = 1:3))
##' predict(eUN.lmm, newdata = newd)
##' predict(eUN.lmm, newdata = newd, keep.newdata = TRUE)
##'
##' ## dynamic prediction
##' newd.d1 <- cbind(newd, Y = c(NA,NA,NA))
##' predict(eUN.lmm, newdata = newd.d1, keep.newdata = TRUE, type = "dynamic")
##' newd.d2 <- cbind(newd, Y = c(6.61,NA,NA))
##' predict(eUN.lmm, newdata = newd.d2, keep.newdata = TRUE, type = "dynamic")
##' newd.d3 <- cbind(newd, Y = c(1,NA,NA))
##' predict(eUN.lmm, newdata = newd.d3, keep.newdata = TRUE, type = "dynamic")
##' newd.d4 <- cbind(newd, Y = c(1,1,NA))
##' predict(eUN.lmm, newdata = newd.d4, keep.newdata = TRUE, type = "dynamic")

## * predict.lmm (code)
##' @export
predict.lmm <- function(object, newdata, se = "estimation", df = !is.null(object$df), type = "static",
                        level = 0.95, keep.newdata = FALSE, ...){
    
    ## ** normalize user imput
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    options <- LMMstar.options()
    name.Y <- object$outcome$var
    name.time <- object$time$var
    if(object$strata$n>1){
        name.strata <- object$strata$var
    }else{
        name.strata <- NULL
    }
    U.time <- object$time$levels
    name.cluster <- object$cluster$var
    type.prediction <- match.arg(type, c("static","dynamic"))
    if(identical(se,FALSE)){
        se <- NULL
    }else if(!is.null(se)){
        if(identical(se,TRUE)){
            stop("Argument \'se\' should not be TRUE but one of \"estimation\", \"residual\", or \"total\".\n",
                 "Typically \"estimation\" will account for uncertainty about the fitted mean, \n",
                 "While \"total\" will reflect the prediction error. \n")
        }
        se <- match.arg(se, c("estimation","residual","total"))
    }
    if(is.null(newdata[[name.cluster]])){
        if(type=="static"){
            newdata[[name.cluster]] <- as.character(1:NROW(newdata))
        }else if(type == "dynamic"){
            if(any(duplicated(newdata[[name.time]]))){
                stop("Duplicated time values found in column ",name.time,".\n",
                     "Consider specifying the cluster variable in argument \'newdata\'. \n")
            }else{
                newdata[[name.cluster]] <- as.character("1")
            }
        }
    }else if(is.factor(newdata[[name.cluster]]) || is.numeric(newdata[[name.cluster]])){
        newdata[[name.cluster]] <- as.character(newdata[[name.cluster]])
    }
    if(!is.null(name.strata)){
        if(name.strata %in% names(newdata) == FALSE){
            stop("The strata column \"",name.strata,"\" in argument \'newdata\' is missing. \n")
        }
        if(is.factor(newdata[[name.strata]]) || is.numeric(newdata[[name.strata]])){
            newdata[[name.strata]] <- as.character(newdata[[name.strata]])
        }
        if(any(newdata[[name.strata]] %in% object$strata$level == FALSE)){
            stop("The strata column \"",name.strata,"\" in argument \'newdata\' contains unknown levels. \n")
        }
    } 
    
    if(type.prediction == "dynamic"){
        if(name.time %in% names(newdata) == FALSE){
            stop("The time column \"",name.time,"\" in argument \'newdata\' is missing and necessary when doing dynamic predictions. \n")
        }
        if(name.Y %in% names(newdata) == FALSE){
            stop("The outcome column \"",name.Y,"\" in argument \'newdata\' is missing and necessary when doing dynamic predictions. \n")
        }
        if(is.factor(newdata[[name.time]]) || is.numeric(newdata[[name.time]])){
            newdata[[name.time]] <- as.character(newdata[[name.time]])
        }
        if(any(newdata[[name.time]] %in% U.time == FALSE)){
            stop("The time column \"",name.time,"\" in argument \'newdata\' should match the existing times. \n",
                 "Existing times: \"",paste0(U.time, collapse = "\" \""),"\".\n")
        }
        test.duplicated <- tapply(newdata[[name.time]],newdata[[name.cluster]], function(iT){any(duplicated(iT))})
        if(any(test.duplicated)){
            stop("The time column \"",name.time,"\" in argument \'newdata\' should not have duplicated values within clusters. \n")
        }
        test.na <- tapply(newdata[[name.Y]],newdata[[name.cluster]], function(iY){any(is.na(iY))})
        if(any(test.na==FALSE)){
            stop("The outcome column \"",name.Y,"\" in argument \'newdata\' should contain at least one missing value for each cluster. \n",
                 "They are used to indicate the time at which prediction should be made. \n")
        }
        if(any("XXXindexXXX" %in% names(newdata))){
            stop("Argument \'newdata\' should not contain a column named \"XXXindexXXX\" as this name is used internally. \n")
        }
        newdata$XXXindexXXX <- 1:NROW(newdata)
    }else if(type.prediction == "static"){
        if(name.time %in% names(newdata) == FALSE){
            if(!is.null(se) && se %in% c("residual","total")){
                stop("The time column \"",name.time,"\" in missing from argument \'newdata\'. \n",
                     "It is necessary for uncertainty calculation involving the residual variance. \n")
            }else{
                newdata[[name.time]] <- U.time[1]
            }
        }
    }

    if(keep.newdata && any(c("estimate","se","df","lower","upper") %in% names(newdata))){
        stop("Argument \'newdata\' should not contain a column named \"estimate\", \"se\", \"lower\", \"upper\", or \"df\" as those names are used internally. \n")
    }

    if(!is.null(se)){
        factor.estimation <- se %in% c("estimation","total")
    }else{
        factor.estimation <- FALSE
    }
    if(!is.null(se)){
        factor.residual <- se %in% c("residual","total")
    }else{
        factor.residual <- FALSE
    }

    ## ** prepare id 
    seq.id <- unique(newdata[[name.cluster]])
    n.id <- length(seq.id)

    ## ** design matrix
    if(type.prediction == "dynamic" || factor.residual){
        newdesign <- model.matrix(object, data = newdata, effects = "all")
        if(any(is.na(newdesign$vcov$pattern.cluster))){
            stop("Could not identify covariance pattern for some clusters. \n")
        }
        X <- newdesign$mean
        pattern.cluster <- newdesign$vcov$X$pattern.cluster
        Upattern <- object$design$vcov$X$Upattern
    }else{
        X <- model.matrix(object, data = newdata, effects = "mean")
    }
    n.obs <- NROW(X)

    ## ** identify variance patterns
    if(type.prediction == "dynamic" || factor.residual){
        Omega <- object$Omega
        for(iO in 1:length(Omega)){
            dimnames(Omega[[iO]]) <- list(U.time[attr(Omega[[iO]],"time")],U.time[attr(Omega[[iO]],"time")])
        }
    }
    
    ## ** extract coefficients and variance
    beta <- coef(object, effects = "mean")
    vcov.beta <- vcov(object, effects = "mean")
        
    ## ** compute predictions
    if(type.prediction == "static"){
        ## compute predictions
        prediction <- X %*% beta
        if(keep.newdata){
            out <- cbind(newdata, estimate = prediction)
        }else{
            out <- data.frame(estimate = prediction)
        }
        ## compute uncertainty about the predictions
        ## prediction.se <- sqrt(diag(X %*% vcov.beta %*% t(X)))
        if(!is.null(se)){
            prediction.var <- rep(0,n.obs)
            if(factor.estimation){
                prediction.var <- prediction.var + rowSums((X %*% vcov.beta) *X)
            }
            if(factor.residual){
                ## find variance corresponding to each observation
                Omega.diag <- data.frame(value = unlist(lapply(Omega,diag)),
                                         pattern = unlist(lapply(1:length(Omega),function(iO){rep(names(Omega)[[iO]],NCOL(Omega[[iO]]))})),
                                         time = U.time[unlist(lapply(Omega,attr,"time"))])
                data.Omega <- data.frame(pattern = pattern.cluster[newdata[[name.cluster]]],
                                         time = newdata[[name.time]])
                index.value <- match(paste(data.Omega$pattern,data.Omega$time,sep="|"), paste(Omega.diag$pattern,Omega.diag$time,sep="|"))
                prediction.var <- prediction.var + Omega.diag$value[index.value]
            }
            out$se <- sqrt(prediction.var)
            out$df <- Inf
        }
    }else if(type.prediction == "dynamic"){
        ## prepare
        vcov.all <- vcov(object, effects = "all")
        param.all <- coef(object, effects = "all")
        name.beta <- colnames(X)
        
        newdata.order <- newdata[order(newdata[[name.time]]),,drop=FALSE]
        prediction <- rep(NA, n.obs)
        prediction.var <- rep(NA, n.obs)

        for(iId in 1:n.id){ ## iId <- 1
            iNewdata <- newdata.order[newdata.order[[name.cluster]] == seq.id[iId],,drop=FALSE] ## subset
            iIndex.con <- which(!is.na(iNewdata[[name.Y]]))
            iLevels.con <- iNewdata[[name.time]][iIndex.con]
            iIndex.pred <- which(is.na(iNewdata[[name.Y]])) ## position in the individal specific dataset
            iPos.pred <- iNewdata$XXXindexXXX[iIndex.pred] ## position in the original dataset
            iLevels.pred <- iNewdata[[name.time]][iIndex.pred]

            iX.con <- X[iNewdata$XXXindexXXX[iIndex.con],,drop=FALSE]
            iX.pred <- X[iPos.pred,,drop=FALSE]
            iOmega.pred <- Omega[[pattern.cluster[seq.id[iId]]]]

            if(length(iLevels.con)==0){ ## static prediction
                prediction[iPos.pred] <- iX.pred %*% beta
                ## iPred.var <- diag(iX.pred %*% vcov.beta %*% t(iX.pred)) + diag(iOmega.pred)
                if(factor.estimation){
                    prediction.var[iPos.pred] <- prediction.var[iPos.pred] + rowSums((iX.pred %*% vcov.beta) * iX.pred)
                }
                if(factor.residual){
                    prediction.var[iPos.pred] <- prediction.var[iPos.pred] + diag(iOmega.pred)
                }
            }else{ ## dynamic prediction
                iOmegaM1.con <- solve(iOmega.pred[iLevels.con,iLevels.con,drop=FALSE])
                iOmega.predcon <- iOmega.pred[iLevels.pred,iLevels.con,drop=FALSE]
                iOmega.conpred <- iOmega.pred[iLevels.con,iLevels.pred,drop=FALSE]

                ## extract current outcome
                iY <- iNewdata[iIndex.con,name.Y]
                ## compute normalized residuals
                iResiduals <- iOmegaM1.con %*% (iY - iX.con %*% beta)
                ## deduce prediction
                prediction[iPos.pred] <- iX.pred %*% beta + iOmega.predcon %*% iResiduals

                if(factor.estimation || factor.residual){
                    prediction.var[iPos.pred] <- 0
                }
                if(factor.estimation){
                    calcPred <- function(x){ ## x <- param.all
                        OO <- .calc_Omega(object = object$design$vcov, param = x, keep.interim = TRUE)[[pattern.cluster[seq.id[iId]]]]
                        dimnames(OO) <- list(U.time[attr(OO,"time")],U.time[attr(OO,"time")])
                        rr <- solve(OO[iLevels.con,iLevels.con,drop=FALSE]) %*% (iY - iX.con %*% x[name.beta])
                        pp <- iX.pred %*% x[name.beta] + OO[iLevels.pred,iLevels.con,drop=FALSE] %*% rr
                        return(pp)
                    }
                    ## calcPred(param.all)
                    iGrad <- numDeriv::jacobian(x = param.all, func = calcPred, method = options$method.numDeriv)
                    ## colnames(iGrad) <- names(param.all)
                    prediction.var[iPos.pred] <- prediction.var[iPos.pred] + rowSums((iGrad %*% vcov.all)*iGrad)
                }
                if(factor.residual){
                    ## iOmega.pred[2,2]*(1-attr(iOmega.pred,"cor")[1,2]^2)
                    prediction.var[iPos.pred] <- prediction.var[iPos.pred] + diag(iOmega.pred[iLevels.pred,iLevels.pred,drop=FALSE] - iOmega.predcon %*% iOmegaM1.con %*% iOmega.conpred)
                }
            }

            
        }

        if(keep.newdata){
            out <- cbind(newdata[,setdiff(colnames(newdata),"XXXindexXXX"),drop=FALSE], estimate = prediction)
            if(!is.null(se)){
                out$se <- sqrt(prediction.var)
                out$df <- ifelse(!is.na(prediction),Inf,NA)
            }
        }else{
            out <- data.frame(estimate = stats::na.omit(prediction))
            if(!is.null(se)){
                out$se <- sqrt(stats::na.omit(prediction.var))
                out$df <- Inf
            }
        }
    }

    ## ** df
    if(df && !is.null(se) && type.prediction == "static"){
        vcov.param <- vcov(object, effects = "all", df = 2, transform.names = FALSE)
        dVcov <- attr(vcov.param,"dVcov")
        attr(vcov.param, "df") <- NULL
        attr(vcov.param, "dVcov") <- NULL
        out$df <- .dfX(X.beta = X, vcov.param = vcov.param, dVcov.param = dVcov)     
    }

    ## ** export
    alpha <- 1-level
    if(!is.null(se)){
        out$lower <- out$estimate + stats::qt(alpha/2, df = out$df) * out$se
        out$upper <- out$estimate + stats::qt(1-alpha/2, df = out$df) * out$se
    }
    rownames(out) <- NULL
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

##----------------------------------------------------------------------
### predict.R ends here
