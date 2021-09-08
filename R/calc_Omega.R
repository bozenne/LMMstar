### calc_Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 21 2021 (18:12) 
## Version: 
## Last-Updated: sep  7 2021 (16:52) 
##           By: Brice Ozenne
##     Update #: 364
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .calc_Omega
.calc_Omega <- function(object, param, keep.interim = FALSE){

    X.var <- object$var
    X.cor <- object$cor
    Upattern <- object$pattern
    indicator <- object$indicator
    index.time <- object$index.time

    ## ** loop over covariance patterns
    out <- lapply(Upattern, function(iP){ ## iP <- Upattern[1]
        ## timepoint in this pattern
        iTime <- index.time[[iP]]
        iNtime <- length(iTime)
        ## design
        iUX.var <- X.var[[iP]]
        iUX.cor <- X.cor[[iP]]

        ## diagonal (standard error)
        Omega.sd <- unname(exp((iUX.var %*% log(param[colnames(iUX.var)]))))
        ## extra diagonal (correlation)
        if(length(iUX.cor)>0){
            Omega.cor <- matrix(iUX.cor %*% param[colnames(iUX.cor)], nrow = iNtime, ncol = iNtime)
        }else{
            Omega.cor <- diag(0, nrow = iNtime, ncol = iNtime)
        }
        ## assemble
        Omega <- diag(as.double(Omega.sd)^2, nrow = iNtime, ncol = iNtime) + Omega.cor * tcrossprod(Omega.sd)

        if(keep.interim){
            attr(Omega,"time") <- iTime
            attr(Omega,"sd") <- Omega.sd
            attr(Omega,"cor") <- Omega.cor
        }
        return(Omega)
    })

    ## ** export
    return(stats::setNames(out,Upattern))
}

## * .calc_dOmega
.calc_dOmega <-  function(object, param, type, Omega, Jacobian){

    ## ** prepare
    name.param <- names(param)
    n.param <- length(param)

    name.sigma <- name.param[type=="sigma"]
    name.k <- name.param[type=="k"]
    name.rho <- name.param[type=="rho"]
    name.paramVar <- intersect(name.param, c(name.sigma,name.k,name.rho))

    X.var <- object$var
    Upattern <- object$pattern
    Upattern.param <- object$param
    indicator <- object$indicator
    index.time <- object$index.time

    ## ** loop over covariance patterns
    out <- lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]

        iTime <- index.time[[iPattern]]
        iNtime <- length(iTime)
        iOmega.sd <- attr(Omega[[iPattern]],"sd")
        iOmega.cor <- attr(Omega[[iPattern]],"cor")
        iOmega <- Omega[[iPattern]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL; attr(iOmega,"time") <- NULL;

        iParam.sigma <- intersect(name.sigma, Upattern.param[[iPattern]])
        n.iParam.sigma <- length(iParam.sigma)
        iParam.k <- intersect(name.k, Upattern.param[[iPattern]])
        n.iParam.k <- length(iParam.k)
        iParam.rho <- intersect(name.rho, Upattern.param[[iPattern]])
        n.iParam.rho <- length(iParam.rho)
        iParamVar <- intersect(name.param, c(iParam.sigma, iParam.k, iParam.rho))
        
        iScore <- stats::setNames(vector(mode = "list", length = length(iParamVar)), iParamVar)

        if(n.iParam.sigma>0){
            for(iSigma in iParam.sigma){ ## iSigma <- iParam.sigma[1]
                ## positions where the sigma-parameter appears (not used since it appears in the whole matrix)
                ## ind.sigma <- indicator[[iPattern]][[iSigma]]
                ## compute derivative
                iParam.dsigma <- param[c(iParam.sigma,iParam.k)]
                iParam.dsigma[iSigma] <- 1
                ## propagate
                idOmega.sigma <- exp(X.var[[iPattern]] %*% log(iParam.dsigma))
                iScore[[iSigma]] <- diag(2*as.double(idOmega.sigma)*as.double(iOmega.sd), nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.sigma %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.sigma))
                ## iScore[[iSigma]] - 2*iOmega/param[iParam.dsigma]
            }
        }
        if(n.iParam.k>0){
            for(iK in iParam.k){ ## iK <- iParam.k[1]
                ## compute derivative
                iParam.dk <- param[c(iParam.sigma,iParam.k)]
                iParam.dk[iK] <- 1
                ## propagate
                idOmega.k <- exp(X.var[[iPattern]] %*% log(iParam.dk)) * X.var[[iPattern]][,iK]
                iScore[[iK]] <- diag(2*as.double(idOmega.k)*as.double(iOmega.sd), nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.k %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k))
                ## iScore[[iK]] - iOmega/param[iK] * (indicator[[iPattern]][[iK]] + tcrossprod(X.var[[iPattern]][,iK]))
            }
        }
        if(n.iParam.rho>0){
            for(iRho in iParam.rho){ ## iRho <- iParam.rho[1]
                ## derivative
                iScore[[iRho]] <- indicator[[iPattern]][[iRho]] * tcrossprod(iOmega.sd)
                ## iScore[[iRho]] - iOmega/param[iRho] * ind.rho
            }
        }

        ## apply transformation
        if(!is.null(Jacobian)){
            ## [dOmega_[11]/d theta_1] ... [dOmega_[11]/d theta_p] %*% Jacobian
            ## [dOmega_[ij]/d theta_1] ... [dOmega_[ij]/d theta_p] %*% Jacobian
            ## [dOmega_[mm]/d theta_1] ... [dOmega_[mm]/d theta_p] %*% Jacobian
            if(any(abs(Jacobian[iParamVar,setdiff(name.paramVar,iParamVar),drop=FALSE])>1e-10)){
                stop("Something went wrong when computing the derivative of the residual variance covariance matrix. \n",
                     "Contact the package manager with a reproducible example generating this error message. \n")
            }
            M.iScore <- do.call(cbind,lapply(iScore,as.double)) %*% Jacobian[iParamVar,iParamVar,drop=FALSE]
            iScore <- stats::setNames(lapply(1:NCOL(M.iScore), function(iCol){matrix(M.iScore[,iCol], nrow = iNtime, ncol = iNtime, byrow = FALSE)}), iParamVar)
        }
        return(iScore)
    })

    ## ** export
    out <- stats::setNames(out,Upattern)
    return(out)
}


## * .calc_d2Omega
.calc_d2Omega <- function(object, param, type, Omega, dOmega, pair, Jacobian, dJacobian){

    ## ** prepare
    name.param <- names(param)
    n.param <- length(param)

    name.sigma <- names(param)[type=="sigma"]
    name.k <- names(param)[type=="k"]
    name.rho <- names(param)[type=="rho"]
    name.paramVar <- intersect(name.param, c(name.sigma,name.k,name.rho))
    
    X.var <- object$var
    Upattern <- object$pattern
    Upattern.param <- object$param
    indicator <- object$indicator
    index.time <- object$index.time
        
    if(!is.null(Jacobian)){
        JacobianM1 <- solve(Jacobian)
    }
    ## ** loop over covariance patterns
    out <- lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]

        iTime <- index.time[[iPattern]]
        iNtime <- length(iTime)
        iOmega.sd <- attr(Omega[[iPattern]],"sd")
        iOmega.cor <- attr(Omega[[iPattern]],"cor")
        iOmega <- Omega[[iPattern]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL; attr(iOmega,"time") <- NULL;
        idOmega <- dOmega[[iPattern]]

        iParam.sigma <- intersect(name.sigma, Upattern.param[[iPattern]])
        n.iParam.sigma <- length(iParam.sigma)
        iParam.k <- intersect(name.k, Upattern.param[[iPattern]])
        n.iParam.k <- length(iParam.k)
        iParam.rho <- intersect(name.rho, Upattern.param[[iPattern]])
        n.iParam.rho <- length(iParam.rho)

        iPair <- pair[[iPattern]]
        n.iPair <- NCOL(iPair)

        iHess <- lapply(1:n.iPair, function(iPair){matrix(0, nrow = iNtime, ncol = iNtime)})

        for(iiPair in 1:n.iPair){ ## iiPair <- 1

            ## name of parameters
            iCoef1 <- iPair[1,iiPair]
            iCoef2 <- iPair[2,iiPair]

            ## type of parameters
            iType1 <- type[iCoef1]
            iType2 <- type[iCoef2]

            ## position where the parameters appear
            iIndicator <- indicator[[iPattern]][[iCoef1]] * indicator[[iPattern]][[iCoef2]]

            if(iType1 == "sigma"){
                if(iType2 == "sigma"){
                    if(iCoef1==iCoef2){
                        iHess[[iiPair]] <- 2 * iIndicator * iOmega / param[iCoef1]^2
                    }else{
                        stop("Cannot compute the Hessian with interacting sigma coefficients. \n")
                    }
                }else if(iType2 == "k"){
                    ## compute derivative
                    iParam.dksigma <- param[c(iParam.sigma,iParam.k)]
                    iParam.dksigma[c(iCoef1,iCoef2)] <- 1
                    ## propagate
                    idOmega.k <- exp(X.var[[iPattern]] %*% log(iParam.dksigma)) * X.var[[iPattern]][,iCoef2]
                    term1 <- diag(4*as.double(idOmega.k)*as.double(iOmega.sd), nrow = iNtime, ncol = iNtime)
                    term2 <- 2 * iOmega.cor * ((idOmega.k) %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k*X.var[[iPattern]][,iCoef2]))
                    iHess[[iiPair]] <- (term1 + term2) * iIndicator
                    ## iHess[[iiPair]] - (2*tcrossprod(X.var[[iPattern]][,iCoef2]) + 2) * iOmega/ (param[iCoef1] * param[iCoef2]) * ind.ksigma
                
                }else if(iType2 == "rho"){
                    iHess[[iiPair]] <- 2 * iIndicator * iOmega / (param[iCoef1] * param[iCoef2])
                }
            }else if(iType1 == "k"){                    
                if(iType2 == "k"){
                    ## compute derivative
                    iParam.dk1 <- param[c(iParam.sigma,iParam.k)]
                    iParam.dk1[iCoef1] <- 1
                    iParam.dk2 <- param[c(iParam.sigma,iParam.k)]
                    iParam.dk2[iCoef2] <- 1
                    ## propagate
                    idOmega.k1 <- exp(X.var[[iPattern]] %*% log(iParam.dk1)) * X.var[[iPattern]][,iCoef1]
                    idOmega.k2 <- exp(X.var[[iPattern]] %*% log(iParam.dk2)) * X.var[[iPattern]][,iCoef2]
                    iHess[[iiPair]] <- diag(2*as.double(idOmega.k1*idOmega.k2), nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.k1 %*% t(idOmega.k2) + idOmega.k2 %*% t(idOmega.k1))

                }else if(iType2 == "rho"){
                    ## compute derivative
                    iParam.dk <- param[c(iParam.sigma,iParam.k)]
                    iParam.dk[iCoef1] <- 1
                    ## propagate
                    idOmega.k <- exp(X.var[[iPattern]] %*% log(iParam.dk)) * X.var[[iPattern]][,iCoef1]
                    iHess[[iiPair]] <- iIndicator * (idOmega.k %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k))
                }
            }
        }

        ## apply transformation
        if(!is.null(Jacobian) || !is.null(dJacobian)){

            iParamVar <- Upattern.param[[iPattern]]
            n.iParamVar <- length(iParamVar)
            
            if(any(abs(JacobianM1[iParamVar,setdiff(name.paramVar,iParamVar),drop=FALSE])>1e-10) || any(abs(dJacobian[iParamVar,setdiff(name.paramVar,iParamVar),,drop=FALSE])>1e-10)){
                stop("Something went wrong when computing the derivative of the residual variance covariance matrix. \n",
                     "Contact the package manager with a reproducible example generating this error message. \n")
            }
            M.iScore <- do.call(cbind,lapply(dOmega[[iPattern]],as.double)) %*% JacobianM1[iParamVar,iParamVar,drop=FALSE]
            iHess2 <- vector(mode = "list", length = n.iPair)

            for(iP in 1:n.iParamVar){ ## iP <- 2
                ## d/d theta_1 = 
                ##  [dOmega_[11]/d2 theta_1] ... [dOmega_[11]/d theta_1 d theta_p] %*% Jacobian + [dOmega_[11]/d theta_1] ... [dOmega_[11]/d theta_p] %*% dJacobian/d theta_1
                ##  [dOmega_[ij]/d2 theta_1] ... [dOmega_[ij]/d theta_1 d theta_p] %*% Jacobian + [dOmega_[ij]/d theta_1] ... [dOmega_[ij]/d theta_p] %*% dJacobian/d theta_1
                ##  [dOmega_[mm]/d2 theta_1] ... [dOmega_[mm]/d theta_1 d theta_p] %*% Jacobian + [dOmega_[mm]/d theta_1] ... [dOmega_[mm]/d theta_p] %*% dJacobian/d theta_1
                iCoef1 <- Upattern.param[[iPattern]][iP]
                ## find all pairs involving the current parameter plus another parameter
                ## e.g. if the current parameter is sigma we want (sigma,sigma) (sigma,k2) (sigma,k3) (sigma k4) (sigma rho)
                ## e.g. if the current parameter is k4 we want (sigma,k4) (k2,k4) (k3,k4) (k4 k4) (k4 rho)
                iIndex.pair  <- which(colSums(iPair == iCoef1)>0)
                ## name of the coefficient of the other pair
                iCoef.pairCoef <- apply(iPair[,iIndex.pair,drop=FALSE], 2, function(iCol){
                    iOut <- setdiff(iCol,iCoef1)
                    if(length(iOut)==0){return(iCoef1)}else{return(iOut)} 
                })
                ## reorder the pairs according to the order of the parameters
                iIndex.pair <- iIndex.pair[match(iCoef.pairCoef, iParamVar)]
                ## apply transformation
                M.iHess2 <- (do.call(cbind,lapply(iHess[iIndex.pair],as.double)) %*% Jacobian[iParamVar,iParamVar,drop=FALSE] + M.iScore %*% dJacobian[iParamVar,iParamVar,iCoef1]) * Jacobian[iCoef1,iCoef1]
                ## convert back to time format
                iHess2[iIndex.pair] <- lapply(1:NCOL(M.iScore), function(iCol){matrix(M.iHess2[,iCol], nrow = iNtime, ncol = iNtime, byrow = FALSE)})
            }
            iHess <- iHess2
        }

        return(iHess)
        
    })

    ## ** export
    out <- stats::setNames(out,Upattern)
    return(out)
} 




##----------------------------------------------------------------------
### calc_Omega.R ends here
