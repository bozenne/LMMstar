### calc_d2Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:18) 
## Version: 
## Last-Updated: aug  4 2023 (17:05) 
##           By: Brice Ozenne
##     Update #: 307
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
    if(missing(Omega)){
        Omega <- .calc_Omega(object, param = param, keep.interim = TRUE)
    }
    if(missing(dOmega)){
        dOmega <- .calc_dOmega(object, param = param, Omega = Omega, Jacobian = Jacobian,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
    }

    object.param <- object$param[is.na(object$param$constraint),,drop=FALSE]
    type <- stats::setNames(object.param$type,object.param$name)
    name.sigma <- object.param[object.param$type=="sigma","name"]
    name.k <- object.param[object.param$type=="k","name"]
    name.rho <- object.param[object.param$type=="rho","name"]
    name.paramVar <- c(name.sigma,name.k,name.rho)

    Upattern <- object$Upattern
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern

    pair.vcov <- object$pair.vcov

    if(identical(transform.sigma,"log") && identical(transform.k,"log") && identical(transform.rho,"atanh")){
        Jacobian <- NULL
        dJacobian <- NULL
        transform <- TRUE
    }else{
        transform.sigma <- "none"
        transform.k <- "none"
        transform.rho <- "none"
        transform <- FALSE
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
        iName.param <- Upattern[iPattern,"param"][[1]]
        iN.time <- Upattern[iPattern,"n.time"]
        if(is.null(iName.param)){return(NULL)}

        iPair <- pair.vcov[[Upattern[iPattern,"name"]]]
        iPair.type <- attr(iPair,"typetype")
        iHess <- stats::setNames(lapply(1:NCOL(iPair), function(iP){matrix(0, nrow = iN.time, ncol = iN.time)}), colnames(iPair))
        
        iOmega <- Omega[[iPattern]]

        if("sigmak.sigmak" %in% names(iPair.type)){

            for(iiPair in iPair.type$sigmak.sigmak){ ## iiPair <- iPair.type$sigmak.sigmak[1]
                iDf.info <- attr(iPair,"index.pair")[[iiPair]]
                if(transform){
                    iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * iDf.info$value2 * iOmega[iDf.info$position]
                }else{ ## no transformation  (other transformations are made through jacobian)
                    if(iPair[1,iiPair]==iPair[2,iiPair]){
                        iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * (iDf.info$value1-1) * iOmega[iDf.info$position] / param[iPair[1,iiPair]]^2
                    }else if(iPair[1,iiPair]!=iPair[2,iiPair]){
                        iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * iDf.info$value2 * iOmega[iDf.info$position] / (param[iPair[1,iiPair]] * param[iPair[2,iiPair]])
                    }                    
                }
            }
            
        }

        if("sigmak.rho" %in% names(iPair.type)){

            for(iiPair in iPair.type$sigmak.rho){ ## iiPair <- iPair.type$sigmak.rho[1]
                iDf.info <- attr(iPair,"index.pair")[[iiPair]]
                if(transform){
                    iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * iDf.info$value2 * iOmega[iDf.info$position] * (1-param[iPair[2,iiPair]]^2) / param[iPair[2,iiPair]]
                }else{ ## no transformation (other transformations are made through jacobian)
                    iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * iDf.info$value2 * iOmega[iDf.info$position] / (param[iPair[1,iiPair]] * param[iPair[2,iiPair]])
                }
            }

        }

        if("rho.rho" %in% names(iPair.type)){
            for(iiPair in iPair.type$rho.rho){ ## iiPair <- iPair.type$rho.rho[1]
                iDf.info <- attr(iPair,"index.pair")[[iiPair]]
                if(transform){
                    if(iPair[1,iiPair]==iPair[2,iiPair]){
                        ## \rho = tanh(x), x = atanh(\rho), dx/d\rho = 1/(1-\rho^2) so d\rho/dx = 1-\rho^2
                        ## \dOmega/\dx = d\rho^a/dx = a \rho^{a-1} (1-\rho^2) = a (\rho^{a-1}-\rho^{a+1})
                        ## \d^2Omega/\dx^2 = a ((a-1)\rho^{a-2}-(a+1)\rho^{a})(1-\rho^2)=a(1-\rho^2)\rho^2((a-1)/\rho^2-(a+1))
                        iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * iOmega[iDf.info$position] * (1-param[iPair[2,iiPair]]^2) * ((iDf.info$value1-1)/param[iPair[1,iiPair]]^2 - (iDf.info$value1+1))
                    }else if(iPair[1,iiPair]!=iPair[2,iiPair]){
                        iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * iDf.info$value2 * iOmega[iDf.info$position] * (1-param[iPair[1,iiPair]]^2) * (1-param[iPair[2,iiPair]]^2) / (param[iPair[1,iiPair]] * param[iPair[2,iiPair]])
                    }
                    
                }else{ ## no transformation  (other transformations are made through jacobian)
                    if(iPair[1,iiPair]==iPair[2,iiPair]){
                        iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * (iDf.info$value1-1) * iOmega[iDf.info$position] / param[iPair[1,iiPair]]^2
                    }else if(iPair[1,iiPair]!=iPair[2,iiPair]){
                        iHess[[iiPair]][iDf.info$position] <- iDf.info$value1 * iDf.info$value2 * iOmega[iDf.info$position] / (param[iPair[1,iiPair]] * param[iPair[2,iiPair]])
                    }                    
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
        
        return(stats::setNames(iHess,colnames(iPair)))        
    })

    ## ** export
    out <- stats::setNames(out,Upattern$name)
    attr(out, "pair") <- pair.vcov
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
