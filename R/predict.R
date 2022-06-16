### predict.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:39) 
## Version: 
## Last-Updated: jun 13 2022 (17:27) 
##           By: Brice Ozenne
##     Update #: 674
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
##' @description Predicted mean value conditional on covariates or on covariates and other outcome values.
##'
##' @param object a \code{lmm} object.
##' @param newdata [data.frame] the covariate values for each cluster.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the predictions. Only relevant if differs from the fitted values.
##' @param level [numeric,0-1] the confidence level of the confidence intervals.
##' @param se [character] Type of uncertainty to be accounted for: estimation of the regression parameters (\code{"estimation"}), residual variance (\code{"residual"}), or both (\code{"total"}).
##' Can also be \code{NULL} to not compute standard error, p-values, and confidence intervals.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the predicted mean. Otherwise a normal distribution is used.
##' @param type [character] Should prediction be made conditional on the covariates only (\code{"static"}) or also on outcome values at other timepoints (\code{"dynamic"}).
##' Can also output the model term (\code{"terms"}, similarly to \code{stats::predict.lm}.
##' @param keep.newdata [logical] Should the argument \code{newdata} be output along side the predicted values?
##' @param se.fit For internal use. When not missing mimic the output of \code{predict.se}. Overwrite argument \code{se}.
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Static prediction are made using the linear predictor \eqn{X\beta} while dynamic prediction uses the conditional normal distribution of the missing outcome given the observed outcomes. So if outcome 1 is observed but not 2, prediction for outcome 2 is obtain by \eqn{X_2\beta + \sigma_{21}\sigma^{-1}_{22}(Y_1-X_1\beta)}. In that case, the uncertainty is computed as the sum of the conditional variance \eqn{\sigma_{22}-\sigma_{21}\sigma^{-1}_{22}\sigma_{12}} plus the uncertainty about the estimated conditional mean (obtained via delta method using numerical derivatives).
##'
##' The model terms are computing by centering the design matrix around the mean value of the covariates used to fit the model.
##' Then the centered design matrix is multiplied by the mean coefficients and columns assigned to the same variable (e.g. three level factor variable) are summed together.
##' 
##' @return A data.frame with 5 columns:\itemize{
##' \item \code{estimate}: predicted mean.
##' \item \code{se}: uncertainty about the predicted mean.
##' \item \code{df}: degree of freedom
##' \item \code{lower}: lower bound of the confidence interval of the predicted mean
##' \item \code{upper}: upper bound of the confidence interval of the predicted mean
##' }
##' except when the argument \code{se.fit} is specified (see \code{predict.lm} for the output format).
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
##' predict(eUN.lmm, newdata = newd, keep.newdata = TRUE, se = "total")
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
predict.lmm <- function(object, newdata, p = NULL, se = "estimation", df = !is.null(object$df), type = "static",
                        level = 0.95, keep.newdata = FALSE, se.fit, ...){
    
    ## ** extract from object
    options <- LMMstar.options()
    name.Y <- object$outcome$var
    U.time <- object$time$levels
    name.time <- attr(object$time$var,"original")
    name.cluster <- attr(object$cluster$var,"original")
    object.mean <- object$design$mean

    table.param <- object$design$param
    name.mu <- table.param$name[table.param$type=="mu"]
    name.theta <- table.param$name

    ## ** normalize user imput
    dots <- list(...)    
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    ## type of prediction
    type.prediction <- match.arg(type, c("static0","static","terms","dynamic")) ## static0: no intercept
    if((type.prediction == "dynamic") && ("rho" %in% table.param$type == FALSE)){
        type.prediction <- "static"
        message("Move to static prediction as there is no correlation parameter in the model. \n")
    }
    if(type.prediction=="static0"){
        type.prediction <- "static"
        keep.intercept <- FALSE
    }else{
        keep.intercept <- TRUE
    }
    ## standard error
    if(!missing(se.fit)){se <- "estimation"}
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

    ## impute cluster when missing (if static) and unambiguous, i.e. no repeated times (id dynamic)
    if(!is.na(name.cluster)){
        if(is.null(newdata[[name.cluster]])){
            if(type.prediction=="static"){
                newdata[[name.cluster]] <- as.character(1:NROW(newdata))
            }else if(type.prediction == "dynamic" && all(!is.na(name.time)) && all(name.time %in% names(newdata))){
                if(any(duplicated(newdata[,name.time, drop=FALSE]))){
                    stop("Duplicated time values found in column ",name.time,".\n",
                         "Consider specifying the cluster variable in argument \'newdata\'. \n")
                }else{
                    newdata[[name.cluster]] <- "1"
                }
            }
        }else if((is.factor(newdata[[name.cluster]]) || is.numeric(newdata[[name.cluster]]))){
            newdata[[name.cluster]] <- as.character(newdata[[name.cluster]])
        }
    }

    ## check data
    if(type.prediction == "dynamic"){

        if(name.Y %in% names(newdata) == FALSE){
            stop("The outcome column \"",name.Y,"\" in argument \'newdata\' is missing and necessary when doing dynamic predictions. \n")
        }
        if(all(is.na(name.cluster))){
            stop("The cluster variable should be explicit when doing dynamic predictions. \n",
                 "Consider re-fitting lmm specifying the argument repetition as ~time|cluster. \n")
        }
        if(name.cluster %in% names(newdata) == FALSE){
            stop("The cluster column \"",name.cluster,"\" in argument \'newdata\' is missing and necessary when doing dynamic predictions. \n")
        }
        if(all(is.na(name.time))){
            stop("The time variable should be explicit when doing dynamic predictions. \n",
                 "Consider re-fitting lmm specifying the argument repetition as ~time|cluster. \n")
        }
        if(any(name.time %in% names(newdata) == FALSE)){
            stop("The time column \"",paste(name.time, collapse=", "),"\" in argument \'newdata\' is missing and necessary when doing dynamic predictions. \n")
        }

        test.level <- sapply(name.time, function(iVar){all(newdata[[iVar]] %in% attr(U.time,"original")[,iVar])})
        if(any(test.level == FALSE)){
            stop("The time column \"",names(which(test.level==FALSE)[1]),"\" in argument \'newdata\' should match the existing times. \n",
                 "Existing times: \"",paste0(attr(U.time,"original")[,names(which(test.level==FALSE)[1])], collapse = "\" \""),"\".\n")
        }
        test.duplicated <- by(newdata[,name.time,drop=FALSE],newdata[[name.cluster]], function(iT){any(duplicated(iT))})
        if(any(unlist(test.duplicated))){
            stop("The time column \"",paste(name.time, collapse = "\", \""),"\" in argument \'newdata\' should not have duplicated values within clusters. \n")
        }

    }else if(type.prediction == "static"){
        
        if(all(!is.na(name.time)) && any(name.time %in% names(newdata) == FALSE) && !is.null(se) && se %in% c("residual","total")){
            stop("The time column \"",paste(name.time[name.time %in% names(newdata) == FALSE], collapse = "\", \""),"\" in missing from argument \'newdata\'. \n",
                 "It is necessary for uncertainty calculation involving the residual variance. \n")
        }

    }

    if(keep.newdata && any(c("estimate","se","df","lower","upper") %in% names(newdata))){
        stop("Argument \'newdata\' should not contain a column named \"estimate\", \"se\", \"lower\", \"upper\", or \"df\" as those names are used internally. \n")
    }

    ## p
    if(!is.null(p) && any(name.theta %in% names(p) == FALSE)){
        stop("Incorrect argument \'p\' - it should be a vector with names containing all parameters. \n",
             "Missing parameter(s): \"",paste(name.theta[name.theta %in% names(p) == FALSE], collapse ="\" \""),"\"")
    }


    ## ** parameters
    if(is.null(p)){
        mu <- coef(object, effects = "mean")
        vcov.mu <- vcov(object, effects = "mean")
        theta <- coef(object, effects = "all")
        vcov.theta <- vcov(object, effects = "all")
    }else{
        theta <- p
        vcov.theta <- vcov(object, p = p, effects = "all")
        mu <- p[name.mu]
        vcov.mu <- vcov(object, p = p, effects = "mean")
    }


    ## ** design matrix
    X <- stats::model.matrix(object, data = newdata, effects = "mean")
    if(!keep.intercept && "(Intercept)" %in% colnames(X)){
        X[,"(Intercept)"] <- 0
    }
    n.obs <- NROW(X)

    ## ** terms        
    if(type.prediction == "terms"){
        Xmean <- colMeans(object.mean)
        Xc <- sweep(X, FUN = "-", MARGIN = 2, STATS = Xmean)
        Xcmu <- sweep(Xc, FUN = "*", MARGIN = 2,  STATS = mu)

            index.n0 <- which(attr(object.mean,"assign")!=0)
            if(length(index.n0)==0){
                Xterm <- matrix(nrow = NROW(newdata), ncol = 0)
            }else{
                Xterm <- do.call(cbind,tapply(index.n0,attr(object.mean,"assign")[index.n0],function(iCol){
                    list(rowSums(Xcmu[,iCol,drop=FALSE]))
                }))
                colnames(Xterm) <- attr(object.mean,"variable")                
            }

        if(any(attr(object.mean,"assign")==0)){
            attr(Xterm, "constant") <- sum(Xmean*mu)
        }
        return(Xterm)
    }

    ## ** identify variance patterns
    if(type.prediction == "dynamic" || factor.residual){
        Omega <- stats::sigma(object, cluster = newdata, simplifies = FALSE)
        index.cluster <- attr(Omega, "index.cluster")
        index.clusterTime <- attr(Omega, "index.clusterTime")
        U.cluster <- names(index.cluster)
        n.cluster <- length(U.cluster)


        if(!is.na(name.cluster) && any(newdata[[name.cluster]] %in% U.cluster == FALSE)){
            stop("Something went wrong when extracting the residual variance-covariance matrices from newdata. \n")
        }

    }
            
    ## ** compute predictions
    if(type.prediction == "static"){
        ## compute predictions
        prediction <- X %*% mu
        if(keep.newdata){
            out <- cbind(newdata, estimate = prediction[,1])
        }else{
            out <- data.frame(estimate = prediction[,1], stringsAsFactors = FALSE)
        }
        ## compute uncertainty about the predictions
        ## prediction.se <- sqrt(diag(X %*% vcov.mu %*% t(X)))
        if(!is.null(se)){
            prediction.var <- rep(0,n.obs)
            if(factor.estimation){
                prediction.var <- prediction.var + rowSums((X %*% vcov.mu) * X)
            }
            if(factor.residual){
                for(iCluster in U.cluster){ ## iCluster <- U.cluster[1]
                    prediction.var[index.cluster[[iCluster]]] <- prediction.var[index.cluster[[iCluster]]] + diag(Omega[[iCluster]])
                }
            }
            out$se <- sqrt(prediction.var)
            out$df <- Inf
        }
    }else if(type.prediction == "dynamic"){

        ## prepare
        prediction <- rep(NA, n.obs)
        prediction.var <- rep(NA, n.obs)
        
        for(iC in U.cluster){ ## iId <- 1
            iNewdata <- newdata[index.cluster[[iC]],,drop=FALSE] ## subset
            iIndex.con <- which(!is.na(iNewdata[[name.Y]]))
            iPos.con <- index.cluster[[iC]][iIndex.con]
            iLevels.con <- U.time[index.clusterTime[[iC]][iIndex.con]]
            
            iIndex.pred <- which(is.na(iNewdata[[name.Y]])) ## position in the individal specific dataset
            iPos.pred <- index.cluster[[iC]][iIndex.pred]
            iLevels.pred <- U.time[index.clusterTime[[iC]][iIndex.pred]]

            iX.con <- X[iPos.con,,drop=FALSE]
            iX.pred <- X[iPos.pred,,drop=FALSE]
            iOmega.pred <- Omega[[iC]]
            
            if(length(iPos.pred)>0 && length(iPos.con)==0){ ## static prediction

                prediction[iPos.pred] <- iX.pred %*% mu
                ## iPred.var <- diag(iX.pred %*% vcov.mu %*% t(iX.pred)) + diag(iOmega.pred)
                if(factor.estimation || factor.residual){
                    prediction.var[iPos.pred] <- 0
                }
                if(factor.estimation){
                    prediction.var[iPos.pred] <- prediction.var[iPos.pred] + rowSums((iX.pred %*% vcov.mu) * iX.pred)
                }
                if(factor.residual){
                    prediction.var[iPos.pred] <- prediction.var[iPos.pred] + diag(iOmega.pred)
                }

            }else if(length(iPos.pred)>0){ ## dynamic prediction
                iOmegaM1.con <- solve(iOmega.pred[iLevels.con,iLevels.con,drop=FALSE])
                iOmega.predcon <- iOmega.pred[iLevels.pred,iLevels.con,drop=FALSE]
                iOmega.conpred <- iOmega.pred[iLevels.con,iLevels.pred,drop=FALSE]

                ## extract current outcome
                iY <- iNewdata[iIndex.con,name.Y]
                ## compute normalized residuals
                iResiduals <- iOmegaM1.con %*% (iY - iX.con %*% mu)
                ## deduce prediction
                prediction[iPos.pred] <- iX.pred %*% mu + iOmega.predcon %*% iResiduals

                if(factor.estimation || factor.residual){
                    prediction.var[iPos.pred] <- 0
                }
                if(factor.estimation){
                    calcPred <- function(x){ ## x <- theta 
                        OO <- stats::sigma(object, p = x, cluster = iNewdata, simplifies = TRUE)
                        rr <- solve(OO[iLevels.con,iLevels.con,drop=FALSE]) %*% (iY - iX.con %*% x[name.mu])
                        pp <- iX.pred %*% x[name.mu] + OO[iLevels.pred,iLevels.con,drop=FALSE] %*% rr
                        return(pp)
                    }
                    ## calcPred(theta)
                    iGrad <- numDeriv::jacobian(x = theta, func = calcPred, method = options$method.numDeriv)
                    ## colnames(iGrad) <- names(name.theta)
                    prediction.var[iPos.pred] <- prediction.var[iPos.pred] + rowSums((iGrad %*% vcov.theta)*iGrad)
                }
                if(factor.residual){
                    ## iOmega.pred[2,2]*(1-attr(iOmega.pred,"cor")[1,2]^2)
                    prediction.var[iPos.pred] <- prediction.var[iPos.pred] + diag(iOmega.pred[iLevels.pred,iLevels.pred,drop=FALSE] - iOmega.predcon %*% iOmegaM1.con %*% iOmega.conpred)
                }
            }

            
        }
        if(keep.newdata){
            out <- cbind(newdata, estimate = prediction)
            if(!is.null(se)){
                out$se <- sqrt(prediction.var)
                out$df <- ifelse(!is.na(prediction),Inf,NA)
            }
        }else{
            out <- data.frame(estimate = stats::na.omit(prediction),stringsAsFactors = FALSE)
            if(!is.null(se)){
                if(NROW(out)>0){
                    out$se <- sqrt(stats::na.omit(prediction.var))
                    out$df <- Inf
                }else{
                    out$se <- numeric(0)
                    out$df <- numeric(0)
                }
            }
        }
    }

    ## ** df
    if(df && !is.null(se) && type.prediction == "static" && se == "estimation"){
        vcov.param <- vcov(object, effects = "all", df = 2, transform.names = FALSE)
        dVcov <- attr(vcov.param,"dVcov")
        attr(vcov.param, "df") <- NULL
        attr(vcov.param, "dVcov") <- NULL
        out$df <- .dfX(X.beta = X, vcov.param = vcov.param, dVcov.param = dVcov)
    }

    ## ** export
    if(!missing(se.fit)){
        return(list(fit = out$estimate,
                    se.fit = out$se,
                    residual.scale = NA,
                    df = out$df))
    }
    alpha <- 1-level
    if(!is.null(se)){
        if(NROW(out)>0){
            out$lower <- out$estimate + stats::qt(alpha/2, df = out$df) * out$se
            out$upper <- out$estimate + stats::qt(1-alpha/2, df = out$df) * out$se
        }else{
            out$lower <- numeric(0)
            out$upper <- numeric(0)
        }
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
    denum <- rowSums((dXTX.dT %*% vcov.vcovbeta) * dXTX.dT)
    out <- 2*diag(prediction.vcov)^2 / denum
    out[denum==0] <- Inf
    ## out <- sapply(1:n.obs, function(iObs){
    ##     2 * prediction.vcov[iObs,iObs]^2 / (dXTX.dT[iObs,] %*% vcov.vcovbeta %*% dXTX.dT[iObs,])
    ## })


    ## ** export
    return(out)
}

##----------------------------------------------------------------------
### predict.R ends here
