### calc_d2Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:18) 
## Version: 
## Last-Updated: okt  8 2026 (14:54) 
##           By: Brice Ozenne
##     Update #: 490
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
##' @param param [named numeric vector] values of the parameters (transformed).
##' @param Omega [list of matrices] residual Variance-Covariance Matrix for each pattern.
##' @param transform.sigma,transform.k,transform.rho [character] transformation used on the variance/correlation coefficients.
##' Only active if \code{"log"}, \code{"log"}, \code{"atanh"}: then the derivative is directly computed on the transformation scale instead of using the Jacobian.
##' @param Upattern [data.frame] optional, used to only evaluate the second derivative of the residual variance-covariance with respect to a subset of patterns.
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
    function(object, param, Omega, transform.sigma, transform.k, transform.rho,
             Upattern) UseMethod(".calc_d2Omega")

## * calc_d2Omega.ID
.calc_d2Omega.ID <- function(object, param, Omega = NULL, transform.sigma, transform.k, transform.rho,
                             Upattern = NULL){

    ## ** prepare
    ## pattern
    if(is.null(Upattern)){
        Upattern <- object$Upattern
    }
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern
    
    ## param
    type <- stats::setNames(object$param$type, object$param$name)
    name.sigma <- object$param[type=="sigma" & is.na(object$param$constraint),"name"]
    name.k <- object$param[type=="k" & is.na(object$param$constraint),"name"]
    name.rho <- object$param[type=="rho" & is.na(object$param$constraint),"name"]
    
    ## Omega
    if(is.null(Omega)){
        Omega <- .calc_Omega(object, param = param,
                             transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, Upattern = Upattern, simplify = FALSE)
    }

    ## ** loop over covariance patterns
    out <- lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

        ## *** patterns
        iPattern.var <- Upattern[iPattern,"var"]
        iPattern.cor <- Upattern[iPattern,"cor"]
        iNtime <- Upattern[iPattern,"n.time"]

        ## *** Omega
        iOmega.sd <- attr(Omega[[iPattern]],"sd")
        iOmega.cor <- attr(Omega[[iPattern]],"cor")
        iOmega <- Omega[[iPattern]]; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL;

        ## *** relevant parameters
        iName.param <- Upattern[iPattern,"param"][[1]]
        ## special case with a single timepoint and therefore no correlation parameter
        ## so maybe no parameter at all if the variance is constrained
        if(is.null(iName.param)){
            return(NULL)
        }else{
            iPair <- object$pair.vcovvcov[which(object$pair.vcovvcov[[Upattern[iPattern,"name"]]]),c("name","param1","param2"),drop=FALSE]
            n.iPair <- NROW(iPair)
            iHess <- replicate(n = n.iPair, matrix(0, nrow = iNtime, ncol = iNtime), simplify = FALSE)
            names(iHess) <- iPair$name
        }
        
        ## *** loop over all pairs of parameters
        for(iiPair in 1:n.iPair){ ## iiPair <- 2

            ## name of parameters
            iCoef1 <- iPair[iiPair,"param1"]
            iCoef2 <- iPair[iiPair,"param2"]

            ## type of parameters
            iType1 <- type[iCoef1]
            iType2 <- type[iCoef2]

            ## first derivatives 
            if(iType1 == "sigma"){
                iDparam1 <- .dsigma_transform(value = param[iCoef1], Omega.sd = iOmega.sd, transform = transform.sigma, power = 1)
            }else if(iType1 == "k"){
                iDparam1 <- .dk_transform(value = param[iCoef1], indicator = (X.var[[iPattern.var]][,"k"] == iCoef1), Omega.sd = iOmega.sd, transform = transform.k, power = 1)
            } ## if iType1 == "rho" only requires second derivative as iCoef1 must be equal to iCoef2 otherwise derivative = 0

            if(iType2 == "sigma"){
                if(iCoef1 != iCoef2){
                    stop("Second Omega derivative cannot handle multiple sigma parameters in a single pattern. \n")
                }
                iDparam2 <- iDparam1
            }else if(iType2 == "k"){
                iDparam2 <- .dk_transform(value = param[iCoef2], indicator = (X.var[[iPattern.var]][,"k"] == iCoef2), Omega.sd = iOmega.sd, transform = transform.k, power = 1)
            }else if(iType2 == "rho" & transform.rho != "cov" & iType1 %in% c("sigma","k")){
                iDparam2 <- .drho_transform(value = param[iCoef2], indicator = (X.cor[[iPattern.cor]][,,"rho"] == iCoef2), transform = transform.rho, power = 1)
            }

            ## second derivative
            if(iType1 == "sigma" && iType2 == "sigma"){
                iDparam12 <- .dsigma_transform(value = param[iCoef1], Omega.sd = iOmega.sd, transform = transform.sigma, power = 2)
            }else if(iType1 == "sigma" && iType2 == "k"){
                iDparam12 <- .dsigmak_transform(value = param[c(iCoef1,iCoef2)], indicator = (X.var[[iPattern.var]][,"k"] == iCoef2), Omega.sd = iOmega.sd, transform = c(transform.sigma,transform.k), power = c(1,1))
            }else if(iType1 == "k" && iType2 == "k"){
                iDparam12 <- (iCoef1 == iCoef2) * .dk_transform(value = param[iCoef1], indicator = (X.var[[iPattern.var]][,"k"] == iCoef1), Omega.sd = iOmega.sd, transform = transform.sigma, power = 2)
            }else if(iType1 == "rho" && transform.rho == "atanh"){ ## for the first type to be rho it means that both types are rho
                iDparam12 <- (iCoef1 == iCoef2) * .drho_transform(value = param[iCoef1], indicator = (X.cor[[iPattern.cor]][,,"rho"] == iCoef1), transform = transform.rho, power = 2)
            }

            ## assemble
            if(iType1 %in% c("sigma","k") && iType2 %in% c("sigma","k")){
                
                ## d^2 u(x,y)u(x,y)^t / dxdy = d u(dx,y)u(x,y)^t + d u(x,y)u(dx,y)^t/dy
                ##                           = u(dx,dy)u(x,y)^t  + u(dx,y)u(x,dy)^t + u(x,dy)u(dx,y)^t + u(x,y)u(dx,dy)^t
                iHess[[iiPair]] <- (tcrossprod(iDparam12, iOmega.sd) + tcrossprod(iDparam1, iDparam2) + tcrossprod(iDparam2, iDparam1) + tcrossprod(iOmega.sd, iDparam12)) * iOmega.cor

            }else if(iType2 %in% c("rho")){

                if(transform.rho == "cov"){
                    ## iHess[[iiPair]] <- matrix(0, nrow = iNtime, ncol = iNtime) ## do nothing this is already the case.
                    ## first derivative with respect to cov was constant
                }else if(iType1 %in% c("sigma","k")){
                    iHess[[iiPair]] <- (tcrossprod(iOmega.sd, iDparam1) + tcrossprod(iDparam1, iOmega.sd))*iDparam2
                }else if(iType1 == "rho" && transform.rho == "atanh"){
                    iHess[[iiPair]] <- tcrossprod(iOmega.sd) * iDparam12
                }
            }
        }

        return(iHess)        
    })

    ## ** export
    return(stats::setNames(out,Upattern$name))
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

## * helper
## ** .dsigmak_transform
##' @param value [numeric vector] parameter value after transformation.
##' @param indicator [logical vector] TRUE where the k parameter is and FALSE otherwise.
##' For instance for [\sigma \sigma*k_1 \sigma*k_2 \sigma*k_3] it would be [FALSE TRUE FALSE FALSE] for k_1
##' @param Omega.sd [numeric vector] square root of the diagonal of the residual variance-covariance matrix.
##' @param transformation [character] transformation for the sigma parameter: \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}.
##' @param power [positive integer vector] order of each of the two derivative.
.dsigmak_transform <- function(value, indicator, Omega.sd, transform, power = c(1,1)){

    if(transform[2] == "none"){
        ## d^2 sigma sigma*k_1 sigma*k_2 / d sigma d k_2 = 0 1 0
        if(all(power==1)){
            out <- indicator
        }else{
            out <- rep(0, length(Omega.sd))
        }
    }else if(transform[2] == "log"){
        ## d^2 sigma exp(log(k)) / d log(k) d sigma = exp(log(k)) = k
        out <- indicator * value[2]
    }else if(transform[2] == "square"){
        ## d sigma sqrt(k^2) / d k^2 d sigma = 1 / (2 sqrt(k^2)) = 1 / (2 k)
        out <- indicator * prod(1/2 - seq(from = 0, to = power-1, by = 1)) / value[2]^power
    }else if(transform[2] == "logsquare"){
        ## d sigma exp(0.5 log(k^2)) / d log(k^2) d sigma = 0.5 exp(0.5 log(k^2)) = 0.5 k
        out <- indicator * value[2] / 2^power
    }else if(transform[2] %in% c("sd","logsd","var","logvar")){
        ## vector of values becomes [f(sigma) sd_2 sd_3 sd_4]
        ## so d^2 / d sd_2 d sigma ---> 0 because no elements contains sigma and sd_2
        out <- rep(0, length(Omega.sd))
    }

    if(transform[1] == "log"){
        ## d. / d log(x) = d. / d x * d x / d log(x) = d. / d x * 1 / (d log(x))/dx = x d . / d x
        out <- out * value[1]
    }else if(transform[1] == "square"){
        ## d. / d x^2 = d. / d x * d x / d x^2 = (1/2x) d . / d x
        out <- out * prod(1/2 - seq(from = 0, to = power-1, by = 1)) / value[1]^power
    }else if(transform[1] == "logsquare"){
        ## d. / d log(x^2) = d. / 2 d log(x) = x/2 d . / d x
        out <- out  * value[1] / 2^power
    }
    
    return(out)
}

##----------------------------------------------------------------------
### calc_d2Omega.R ends here
