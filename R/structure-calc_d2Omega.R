### calc_d2Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:18) 
## Version: 
## Last-Updated: okt  1 2021 (17:06) 
##           By: Brice Ozenne
##     Update #: 57
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calc_d2Omega
##' @title Second Derivative of the Residual Variance-Covariance Matrix
##' @description Second derivative of the residual variance-covariance matrix for given parameter values.
##' @noRd
##'
##' @param structure [structure]
##' @param param [named numeric vector] values of the parameters.
##' @param Omega [list of matrices] Residual Variance-Covariance Matrix for each pattern.
##' @param dOmega [list of matrices] First derivative of the residual Variance-Covariance Matrix for each pattern.
##' @param Jacobian [matrix] Jacobian of the reparametrisation.
##' @param dJacobian [array] First derivative of the Jacobian of the reparametrisation.
##'
##' @keywords internal
##' 
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$gender <- c("M","F")[as.numeric(gastricbypassL$id) %% 2+1]
##' dd <- gastricbypassL[!duplicated(gastricbypassL[,c("time","gender")]),]
##' 
##' ## independence
##' Sid1 <- skeleton(IND(~1, var.time = "time"), data = dd)
##' Sid4 <- skeleton(IND(~1|id, var.time = "time"), data = dd)
##' Sdiag1 <- skeleton(IND(~visit), data = dd)
##' Sdiag4 <- skeleton(IND(~visit|id), data = dd)
##' Sdiag24 <- skeleton(IND(~visit+gender|id, var.time = "time"), data = dd)
##' param24 <- setNames(c(1,2,2,3,3,3,5,3),Sdiag24$param$name)
##'
##' .calc_d2Omega(Sid4, param = c(sigma = 2))
##' .calc_d2Omega(Sdiag1, param = c(sigma = 1, k.visit2 = 2, k.visit3 = 3, k.visit4 = 4))
##' .calc_d2Omega(Sdiag4, param = c(sigma = 1, k.visit2 = 2, k.visit3 = 3, k.visit4 = 4))
##' .calc_d2Omega(Sdiag24, param = param24)
##' 
##' ## compound symmetry
##' Scs4 <- skeleton(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' Scs24 <- skeleton(CS(gender~time|id), data = gastricbypassL)
##' 
##' .calc_d2Omega(Scs4, param = c(sigma = 1,rho=0.5))
##' .calc_d2Omega(Scs4, param = c(sigma = 2,rho=0.5))
##' .calc_d2Omega(Scs24, param = c("sigma:F" = 2, "sigma:M" = 1,
##'                             "rho:F"=0.5, "rho:M"=0.25))
##' 
##' ## unstructured
##' Sun4 <- skeleton(UN(~visit|id), data = gastricbypassL)
##' param4 <- setNames(c(1,1.1,1.2,1.3,0.5,0.45,0.55,0.7,0.1,0.2),Sun4$param$name)
##' Sun24 <- skeleton(UN(gender~visit|id), data = gastricbypassL)
##' param24 <- setNames(c(param4,param4*1.1),Sun24$param$name)
##' 
##' .calc_d2Omega(Sun4, param = param4)
##' .calc_d2Omega(Sun24, param = param24)
`.calc_d2Omega` <-
    function(object, param, Omega, dOmega, Jacobian, dJacobian) UseMethod(".calc_d2Omega")

## * calc_d2Omega.UN
.calc_d2Omega.UN <- function(object, param, Omega, dOmega, Jacobian = NULL, dJacobian = NULL){

    ## ** prepare
    type <- stats::setNames(object$param$type,object$param$name) 
    name.sigma <- object$param$name[type=="sigma"]
    name.k <- object$param$name[type=="k"]
    name.rho <- object$param$name[type=="rho"]
    name.paramVar <- c(name.sigma,name.k,name.rho)

    param <- param[name.paramVar]
    name.param <- names(param)
    if(missing(Omega)){
        Omega <- .calc_Omega(object, param = param, keep.interim = TRUE)
    }
    if(missing(dOmega)){
        dOmega <- .calc_dOmega(object, param = param, Omega = Omega, Jacobian = Jacobian)
    }
    
    Upattern <- object$X$Upattern
    n.Upattern <- NROW(Upattern)
    pattern.cluster <- object$X$pattern.cluster
    X.var <- object$X$var
    X.cor <- object$X$cor
    if(!is.null(Jacobian)){
        test.nooffdiag <- all(abs(c(Jacobian[lower.tri(Jacobian,diag=FALSE)],Jacobian[upper.tri(Jacobian,diag=FALSE)]))<1e-10)
        if(test.nooffdiag){
            JacobianM1 <- Jacobian
            diag(JacobianM1) <- 1/diag(Jacobian)
            ## range(JacobianM1 - solve(Jacobian))
        }else{
            JacobianM1 <- solve(Jacobian)
        }
        
    }
    ## ** loop over covariance patterns
    out <- lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

        iPattern.var <- Upattern[iPattern,"var"]
        iPattern.cor <- Upattern[iPattern,"cor"]
        iTime <- Upattern[iPattern,"time"][[1]]
        iNtime <- length(iTime)
        iName.param <- Upattern[iPattern,"param"][[1]]

        iOmega.sd <- attr(Omega[[iPattern]],"sd")
        iOmega.cor <- attr(Omega[[iPattern]],"cor")
        iOmega <- Omega[[iPattern]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL; attr(iOmega,"time") <- NULL;

        iParam.sigma <- intersect(name.sigma, iName.param)
        n.iParam.sigma <- length(iParam.sigma)
        iParam.k <- intersect(name.k, iName.param)
        n.iParam.k <- length(iParam.k)
        iParam.rho <- intersect(name.rho, iName.param)
        n.iParam.rho <- length(iParam.rho)
        iParamVar <- c(iParam.sigma, iParam.k, iParam.rho)

        iScore <- stats::setNames(vector(mode = "list", length = length(iParamVar)), iParamVar)
        iX.var <- X.var[[iPattern.var]][,c(iParam.sigma,iParam.k),drop=FALSE]
        iX.cor <- X.cor[[iPattern.cor]][,c(iParam.rho),drop=FALSE]
        iIndicator <- c(attr(X.var[[iPattern.var]],"indicator.param"),attr(X.cor[[iPattern.cor]],"indicator.param"))

        iPair <- object$pair.varcoef[[Upattern[iPattern,"name"]]]
        n.iPair <- NCOL(iPair)

        iHess <- lapply(1:n.iPair, function(iPair){matrix(0, nrow = iNtime, ncol = iNtime)})

        for(iiPair in 1:n.iPair){ ## iiPair <- 2
            ## name of parameters
            iCoef1 <- iPair[1,iiPair]
            iCoef2 <- iPair[2,iiPair]

            ## type of parameters
            iType1 <- type[iCoef1]
            iType2 <- type[iCoef2]

            ## indicator
            iIndicator12 <- intersect(iIndicator[[iCoef1]],iIndicator[[iCoef2]])
            
            if(iType1 == "sigma"){
                if(iType2 == "sigma"){
                    if(iCoef1==iCoef2){
                        iHess[[iiPair]][iIndicator12] <- 2 * iOmega[iIndicator12] / param[iCoef1]^2
                    }else{
                        stop("Cannot compute the Hessian with interacting sigma coefficients. \n")
                    }
                }else if(iType2 == "k"){
                    ## compute derivative
                    iParam.dksigma <- param[c(iParam.sigma,iParam.k)]
                    iParam.dksigma[c(iCoef1,iCoef2)] <- 1
                    ## propagate
                    idOmega.k <- exp(iX.var %*% log(iParam.dksigma)) * iX.var[,iCoef2]
                    term1 <- diag(4*as.double(idOmega.k)*as.double(iOmega.sd), nrow = iNtime, ncol = iNtime)
                    term2 <- 2 * iOmega.cor * ((idOmega.k) %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k*iX.var[,iCoef2]))
                    iHess[[iiPair]][iIndicator12] <- (term1 + term2)[iIndicator12]
                    ## iHess[[iiPair]] - (2*tcrossprod(X.var[[iPattern.var]][,iCoef2]) + 2) * iOmega/ (param[iCoef1] * param[iCoef2]) * ind.ksigma
                
                }else if(iType2 == "rho"){
                    iHess[[iiPair]][iIndicator12] <- 2 * iOmega[iIndicator12] / (param[iCoef1] * param[iCoef2])
                }
            }else if(iType1 == "k"){                    
                if(iType2 == "k"){
                    ## compute derivative
                    iParam.dk1 <- param[c(iParam.sigma,iParam.k)]
                    iParam.dk1[iCoef1] <- 1
                    iParam.dk2 <- param[c(iParam.sigma,iParam.k)]
                    iParam.dk2[iCoef2] <- 1
                    ## propagate
                    idOmega.k1 <- exp(iX.var %*% log(iParam.dk1)) * iX.var[,iCoef1]
                    idOmega.k2 <- exp(iX.var %*% log(iParam.dk2)) * iX.var[,iCoef2]
                    iHess[[iiPair]] <- diag(2*as.double(idOmega.k1*idOmega.k2), nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.k1 %*% t(idOmega.k2) + idOmega.k2 %*% t(idOmega.k1))
                }else if(iType2 == "rho"){
                    ## compute derivative
                    iParam.dk <- param[c(iParam.sigma,iParam.k)]
                    iParam.dk[iCoef1] <- 1
                    ## propagate
                    idOmega.k <- exp(iX.var %*% log(iParam.dk)) * iX.var[,iCoef1]
                    iHess[[iiPair]][iIndicator12] <- (idOmega.k %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k))[iIndicator12]
                }
            }
        }

        ## apply transformation
        if(!is.null(Jacobian) || !is.null(dJacobian)){

            iParamVar <- Upattern[iPattern,"param"][[1]]
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
                iCoef1 <- iParamVar[iP]
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
    out <- stats::setNames(out,Upattern$name)
    return(out)
} 

## * calc_d2Omega.IND
.calc_d2Omega.IND <- .calc_d2Omega.UN

## * calc_d2Omega.CS
.calc_d2Omega.CS <- .calc_d2Omega.UN



##----------------------------------------------------------------------
### calc_d2Omega.R ends here
