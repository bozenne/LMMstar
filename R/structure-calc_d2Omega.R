### calc_d2Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:18) 
## Version: 
## Last-Updated: jul 26 2023 (14:27) 
##           By: Brice Ozenne
##     Update #: 217
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
##' @param transform.sigma,transform.k,transform.rho [character] Transformation used on the variance/correlation coefficients.
##' Only active if \code{"log"}, \code{"log"}, \code{"atanh"}: then the derivative is directly computed on the transformation scale instead of using the Jacobian.
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
    function(object, param, Omega, dOmega, Jacobian, dJacobian,
             transform.sigma, transform.k, transform.rho) UseMethod(".calc_d2Omega")

## * calc_d2Omega.ID
.calc_d2Omega.ID <- function(object, param, Omega, dOmega, Jacobian = NULL, dJacobian = NULL,
                              transform.sigma = NULL, transform.k = NULL, transform.rho = NULL){

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
        dOmega <- .calc_dOmega(object, param = param, Omega = Omega, Jacobian = Jacobian,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
    }

    Upattern <- object$Upattern
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern
    if(identical(transform.sigma,"log") && identical(transform.k,"log") && identical(transform.rho,"atanh")){
        Jacobian <- NULL
        dJacobian <- NULL
    }else{
        transform.sigma <- "none"
        transform.k <- "none"
        transform.rho <- "none"
    }
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
        iNtime <- Upattern[iPattern,"n.time"]
        iName.param <- Upattern[iPattern,"param"][[1]]
        if(is.null(iName.param)){return(NULL)}

        iOmega.sd <- attr(Omega[[iPattern]],"sd")
        iOmega.var <- tcrossprod(iOmega.sd)
        iOmega.cor <- attr(Omega[[iPattern]],"cor")
        iOmega <- Omega[[iPattern]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL; attr(iOmega,"time") <- NULL;

        iScore <- stats::setNames(vector(mode = "list", length = length(iName.param)), iName.param)

        iPair <- object$pair.vcov[[Upattern[iPattern,"name"]]]
        n.iPair <- NCOL(iPair)

        iHess <- lapply(1:n.iPair, function(iPair){matrix(0, nrow = iNtime, ncol = iNtime)})

        for(iiPair in 1:n.iPair){ ## iiPair <- 2

            ## name of parameters
            iCoef1 <- iPair[1,iiPair]
            iCoef2 <- iPair[2,iiPair]

            ## type of parameters
            iType1 <- type[iCoef1]
            iType2 <- type[iCoef2]

            ## indicators
            iMindicator.var1 <- attr(X.var[[iPattern.var]],"Mindicator.param")[[iCoef1]]
            iMindicator.var2 <- attr(X.var[[iPattern.var]],"Mindicator.param")[[iCoef2]]
            iIndicator.cor2 <- attr(X.cor[[iPattern.cor]],"indicator.param")[[iCoef2]]

            if(iType1 == "sigma"){
                if(iType2 == "sigma"){
                    if(iCoef1!=iCoef2){
                        stop("Cannot compute the Hessian with interacting sigma coefficients. \n")
                    }
                    if(transform.sigma == "log"){
                        iHess[[iiPair]] <- 4 * iOmega
                    }else{ ## no transformation  (other transformations are made through jacobian)
                        iHess[[iiPair]] <- 2 * iOmega / param[iCoef1]^2
                    }
                }else if(iType2 == "k"){
                    if(transform.k == "log"){
                        iHess[[iiPair]] <- iMindicator.var1 * iMindicator.var2 * iOmega
                    }else{ ## no transformation  (other transformations are made through jacobian)
                        iHess[[iiPair]] <- iMindicator.var1 * iMindicator.var2 * iOmega / (param[iCoef1]*param[iCoef2])
                    }
                }else if(iType2 == "rho"){
                    if(transform.rho == "atanh"){
                        iHess[[iiPair]][iIndicator.cor2] <- 2 * iOmega.var[iIndicator.cor2] * (1-param[iCoef2]^2)
                    }else{ ## no transformation (other transformations are made through jacobian)
                        iHess[[iiPair]][iIndicator.cor2] <- 2 * iOmega.var[iIndicator.cor2] / param[iCoef1]
                    }
                }
            }else if(iType1 == "k"){                    
                if(iType2 == "k"){
                    if(transform.k == "log"){
                        iHess[[iiPair]] <- iMindicator.var1 * iMindicator.var2 * iOmega
                    }else{ ## no transformation  (other transformations are made through jacobian)
                        if(iCoef1==iCoef2){
                            iHess[[iiPair]] <- iMindicator.var1*(iMindicator.var1-1) * iOmega / param[iCoef1]^2
                        }else{
                            iHess[[iiPair]] <- iMindicator.var1 * iMindicator.var2 * iOmega / (param[iCoef1]*param[iCoef2])
                        }
                    }
                }else if(iType2 == "rho"){
                    if(transform.k == "log"){
                        iHess[[iiPair]][iIndicator.cor2] <- iMindicator.var1[iIndicator.cor2] * iOmega.var[iIndicator.cor2] * (1-param[iCoef2]^2)
                    }else{ ## no transformation  (other transformations are made through jacobian)
                        iHess[[iiPair]][iIndicator.cor2] <- iMindicator.var1[iIndicator.cor2] * iOmega.var[iIndicator.cor2] / param[iCoef1]
                    }
                }
            }else if(transform.rho == "atanh" && iType1 == "rho" && iType2 == "rho" && iCoef1 == iCoef2){
                iHess[[iiPair]][iIndicator.cor2] <- - 2 * iOmega.var[iIndicator.cor2] * param[iCoef2] * (1 - param[iCoef2]^2)
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
.calc_d2Omega.IND <- .calc_d2Omega.ID

## * calc_d2Omega.CS
.calc_d2Omega.CS <- .calc_d2Omega.ID

## * calc_d2Omega.RE
.calc_d2Omega.RE <- .calc_d2Omega.ID

## * calc_d2Omega.TOEPLITZ
.calc_d2Omega.TOEPLITZ <- .calc_d2Omega.ID

## * calc_d2Omega.UN
.calc_d2Omega.UN <- .calc_d2Omega.ID

## * calc_d2Omega.CUSTOM
.calc_d2Omega.CUSTOM <- function(object, param, Omega, dOmega, Jacobian = NULL, dJacobian = NULL,
                                 transform.sigma = NULL, transform.k = NULL, transform.rho = NULL){

    Upattern <- object$Upattern
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern
    FCT.sigma <- object$FCT.sigma
    FCT.rho <- object$FCT.rho
    dFCT.sigma <- object$dFCT.sigma
    dFCT.rho <- object$dFCT.rho
    d2FCT.sigma <- object$d2FCT.sigma
    d2FCT.rho <- object$d2FCT.rho
    name.sigma <- object$param[object$param$type=="sigma","name"]
    name.rho <- object$param[object$param$type=="rho","name"]
    pair.varcoef <- object$pair.vcov

    if(!is.null(FCT.sigma) && is.null(d2FCT.sigma) || !is.null(FCT.rho) && is.null(d2FCT.rho) ){

        ## second derivative
        ## unlist(.calc_dOmega.CUSTOM(object, param = param, Omega = Omega))
        vec.dOmega <- numDeriv::jacobian(func = function(x){
            unlist(.calc_dOmega.CUSTOM(object, param = x, Omega = Omega))
        }, x = param[c(name.sigma,name.rho)])

        ## indicator of pattern
        vec.pattern <- unlist(lapply(names(dOmega), function(iName){ ## iName <- names(dOmega)[1]
            iParam <- names(dOmega[[iName]])
            iNtime <- Upattern[Upattern$name==iName,"n.time"]
            iOut <- lapply(iParam, function(iP){matrix(iName, nrow = iNtime, ncol = iNtime)})
            return(iOut)
        }))

        ## indicator of param
        vec.param <- unlist(lapply(names(dOmega), function(iName){ ## iName <- names(dOmega)[1]
            iParam <- names(dOmega[[iName]])
            iNtime <- Upattern[Upattern$name==iName,"n.time"]
            iOut <- lapply(iParam, function(iP){matrix(iP, nrow = iNtime, ncol = iNtime)})
            return(iOut)
        }))

        ## matrix to list of matrices according to param and pattern
        vec.patternXparam <- paste0(vec.pattern,"_",vec.param)
        list.d2Omega <- by(data = data.frame(pattern = vec.pattern, param = vec.param, vec.dOmega), INDICES = vec.patternXparam, FUN = function(idOmega){ ## idOmega <- vec.dOmega[1:16,]
            iOut <- apply(idOmega[,-(1:2),drop=FALSE], MARGIN = 2, simplify = FALSE, function(iVec){
                iNtime <- sqrt(length(iVec))
                matrix(iVec, nrow = iNtime, ncol = iNtime)
            })
            names(iOut) <- c(name.sigma,name.rho)
            return(iOut)
        }, simplify = FALSE)[unique(vec.patternXparam)]

        ## normalize to expected output
        out <- stats::setNames(vector(mode = "list", length = n.Upattern), Upattern$name)
        for(iPattern in 1:n.Upattern){
            
            iIndex.pattern <- vec.pattern[!duplicated(vec.patternXparam)]==Upattern$name[iPattern]
            iParam.pattern <- vec.param[!duplicated(vec.patternXparam)][iIndex.pattern]
            iList.d2Omega <- list.d2Omega[iIndex.pattern]

            out[[iPattern]] <- apply(pair.varcoef[[Upattern$name[iPattern]]], MARGIN = 2, simplify = FALSE, function(iCol){ ## iCol <- pair.varcoef[[Upattern$name[iPattern]]][,1]
                iList.d2Omega[[which(iParam.pattern==iCol[1])]][[iCol[2]]]
            })

        }
    }else{

        out <- stats::setNames(vector(mode = "list", length = n.Upattern), Upattern$name)
        for(iPattern in 1:n.Upattern){ ## iPattern <- 1

            iPattern.var <- Upattern$var[iPattern]
            iNtime <- Upattern$n.time[iPattern]
            iX.var <- X.var[[iPattern.var]]
            iOmega.sd <- attr(Omega[[iPattern]], "sd")
            idOmega.sd <- dFCT.sigma(p = param[name.sigma], n.time = iNtime, X = iX.var)
            id2Omega.sd <- d2FCT.sigma(p = param[name.sigma], n.time = iNtime, X = iX.var)

            if(iNtime > 1 && !is.na(Upattern$cor[iPattern])){
                iPattern.cor <- Upattern$cor[iPattern]
                iX.cor <- X.cor[[iPattern.cor]]
                iOmega.cor <- attr(Omega[[iPattern]], "cor")
                idOmega.cor <- dFCT.rho(p = param[name.rho], n.time = iNtime, X = iX.cor)
                id2Omega.cor <- d2FCT.rho(p = param[name.rho], n.time = iNtime, X = iX.cor)
            }

            out[[iPattern]] <- apply(pair.varcoef[[Upattern$name[iPattern]]], MARGIN = 2, simplify = FALSE, function(iCol){ ## iCol <- pair.varcoef[[Upattern$name[iPattern]]][,1]
                iDeriv <- matrix(0, iNtime, iNtime)
                if(iCol[1] %in% name.sigma && iCol[2] %in% name.sigma){
                    if(iCol[1]==iCol[2]){
                        iDeriv1 <- idOmega.sd[[iCol[1]]]
                        iDeriv2 <- id2Omega.sd[[iCol[1]]]
                    
                        ## diagonal sigma terms: f(a)^2 --> 2f(a)f'(a) --> 2[f'(a)f'(a)+f(a)f''(a)]
                        iDeriv <- iDeriv + 2*diag(iDeriv1^2 + iOmega.sd*iDeriv2, nrow = iNtime, ncol = iNtime)
                        ## off-diagonal sigma terms: f(a)f(b) --> f''(a)f(b)
                        if(iNtime > 1 && !is.null(X.cor)){
                            iDeriv <- iDeriv + iOmega.cor * (iDeriv2 %*% t(iOmega.sd) + iOmega.sd %*% t(iDeriv2))
                        }
                    }else if(iNtime > 1 && !is.null(X.cor)){
                        iDeriv1 <- idOmega.sd[[iCol[1]]]
                        iDeriv2 <- idOmega.sd[[iCol[2]]]
                    
                        ## off-diagonal sigma terms: f(a)f(b) --> f'(a)f'(b)
                        iDeriv <- iDeriv + iOmega.cor * (iDeriv1 %*% t(iDeriv2) + iDeriv2 %*% t(iDeriv1))
                    }
                }else if(iCol[1] %in% name.rho && iCol[2] %in% name.rho){
                    ## diagonal sigma terms: 0
                    ## off-diagonal sigma terms: f(a)f(b)f(c) --> f(a)f(b)f''(c)
                    if(iCol[1]==iCol[2]){
                        iDeriv <- iDeriv + id2Omega.cor[[iCol[1]]] * tcrossprod(iOmega.sd)
                    }
                }else{
                    iDeriv1 <- idOmega.sd[[iCol[iCol %in% name.sigma]]]
                    iDeriv2 <- idOmega.cor[[iCol[iCol %in% name.rho]]]
                    ## diagonal sigma terms: 0
                    ## off-diagonal sigma terms: f(a)f(b)f(c) --> f'(a)f(b)f'(c)
                    iDeriv <- iDeriv + iDeriv2 * (iDeriv1 %*% t(iOmega.sd) + iOmega.sd %*% t(iDeriv1))
                }

                return(iDeriv)
            })

        }
    }

    return(out)
}


##----------------------------------------------------------------------
### calc_d2Omega.R ends here
