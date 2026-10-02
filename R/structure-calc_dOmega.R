### calc_dOmega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:18) 
## Version: 
## Last-Updated: okt  2 2026 (16:38) 
##           By: Brice Ozenne
##     Update #: 413
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calc_dOmega
##' @title First Derivative of the Residual Variance-Covariance Matrix
##' @description First derivative of the residual variance-covariance matrix for given parameter values.
##' @noRd
##'
##' @param structure [structure]
##' @param param [named numeric vector] values of the parameters (transformed).
##' @param dOmega [list of matrices] first derivative of the residual Variance-Covariance Matrix for each pattern.
##' @param transform.sigma,transform.k,transform.rho [character] Transformation used on the variance/correlation coefficients.
##' @param Upattern [data.frame] Optional, used to only evaluate the derivative of the residual variance-covariance with respect to a subset of patterns.
##' 
##' @keywords internal
##' 
`.calc_dOmega` <-
    function(object, param, Omega, transform.sigma, transform.k, transform.rho, Upattern) UseMethod(".calc_dOmega")

## * calc_dOmega.ID
.calc_dOmega.ID <- function(object, param, Omega = NULL, transform.sigma, transform.k, transform.rho,
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
    type <- object$param$type

    name.sigma <- object$param$name[type=="sigma"]
    name.k <- object$param$name[type=="k"]
    name.rho <- object$param$name[type=="rho"]

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

        ## *** relevant parameters
        iName.param <- Upattern[iPattern,"param"][[1]]
        ## special case with a single timepoint and therefore no correlation parameter
        ## so maybe no parameter at all if the variance is constrained
        if(is.null(iName.param)){
            return(NULL)
        }else{
            iScore <- stats::setNames(vector(mode = "list", length = length(iName.param)), iName.param)
        }
        iParam.sigma <- intersect(iName.param, name.sigma)
        iParam.k <- intersect(iName.param, name.k)
        iParam.rho <- intersect(iName.param, name.rho)
        
        ## *** sigma & k parameter
        if(length(iParam.sigma)+length(iParam.k)>0){
            iX.var <- X.var[[Upattern[iPattern,"var"]]]
            
            iScore[c(iParam.sigma,iParam.k)] <- lapply(c(iParam.sigma,iParam.k), function(iParam){ ## iParam <- iParam.k[1]
                if(iParam %in% iParam.sigma){
                    iDparam.sd <- .dsigma_transform(value = param[iParam], Omega.sd = iOmega.sd, transform = transform.sigma, power = 1)
                }else{                
                    iDparam.sd <- .dk_transform(value = param[iParam], indicator = (iX.var[,"k"]==iParam), Omega.sd = iOmega.sd, transform = transform.k, power = 1)
                }
                if(transform.rho %in% c("none","atanh")){ ## all terms
                    ## (uu^t)' = u'u^t + uu'^t
                    iDomega.Omega <- (tcrossprod(iOmega.sd, iDparam.sd) + tcrossprod(iDparam.sd, iOmega.sd))*iOmega.cor
                }else if(transform.rho %in% c("cov")){ ## only diagonal terms
                    iDomega.Omega <- 2*diag(iOmega.sd * iDparam.sd)*iOmega.cor
                }
                return(iDomega.Omega)
            })
        }

        ## *** rho parameter
        ## transformations: none, atanh, cov
        if(length(iParam.rho)>0){
            iX.cor <- X.cor[[Upattern[iPattern,"cor"]]]
            
            iScore[iParam.rho] <- lapply(iParam.rho, function(iParam){ ## iParam <- iParam.rho[1]
                if(transform.rho == "cov"){
                    iDomega.Omega <- (iX.cor[,,"rho"]==iParam)
                }else{
                    iDomega.cor <- .drho_transform(value = param[iParam], indicator = (iX.cor[,,"rho"]==iParam), transform = transform.rho, power = 1)
                    iDomega.Omega <- tcrossprod(iOmega.sd, iOmega.sd)*iDomega.cor
                }                
                return(iDomega.Omega)
            })
        }

        return(iScore)
    })

    ## ** export
    return(stats::setNames(out,Upattern$name))
}

## * calc_dOmega.IND
.calc_dOmega.IND <- .calc_dOmega.ID

## * calc_dOmega.CS
.calc_dOmega.CS <- .calc_dOmega.ID

## * calc_dOmega.RE
.calc_dOmega.RE <- .calc_dOmega.ID

## * calc_dOmega.TOEPLITZ
.calc_dOmega.TOEPLITZ <- .calc_dOmega.ID

## * calc_dOmega.UN
.calc_dOmega.UN <- .calc_dOmega.ID

## * calc_dOmega.CUSTOM
.calc_dOmega.CUSTOM <- function(object, param, Omega, Jacobian = NULL,
                                transform.sigma = NULL, transform.k = NULL, transform.rho = NULL){

    ## ** prepare
    Upattern <- object$Upattern
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern
    FCT.sigma <- object$FCT.sigma
    FCT.rho <- object$FCT.rho
    dFCT.sigma <- object$dFCT.sigma
    dFCT.rho <- object$dFCT.rho
    name.sigma <- object$param[object$param$type=="sigma","name"]
    name.rho <- object$param[object$param$type=="rho","name"]

    ## ** Score
    if(!is.null(FCT.sigma) && is.null(dFCT.sigma) || !is.null(FCT.rho) && is.null(dFCT.rho) ){

        ## unlist(.calc_Omega.CUSTOM(object, param = param, simplify = TRUE))
        vec.dOmega <- numDeriv::jacobian(func = function(x){
            unlist(.calc_Omega.CUSTOM(object, param = x, simplify = TRUE))
        }, x = param[c(name.sigma,name.rho)])

        vec.pattern <- unlist(lapply(names(Omega), function(iName){
            iNtime <- Upattern[Upattern$name==iName,"n.time"]
            iOut <- matrix(iName, nrow = iNtime, ncol = iNtime)
        }))
        
        out <- by(data = vec.dOmega, INDICES = vec.pattern, FUN = function(idOmega){ ## idOmega <- vec.dOmega[1:16,]
            iOut <- apply(idOmega, MARGIN = 2, simplify = FALSE, function(iVec){
                iNtime <- sqrt(length(iVec))
                matrix(iVec, nrow = iNtime, ncol = iNtime)
            })
            names(iOut) <- c(name.sigma,name.rho)
            return(iOut)
        }, simplify = FALSE)
        class(out) <- "list"
        attr(out,"call") <- NULL

    }else{

        out <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

            ## derivative of sd with respect to the variance parameters
            iPattern.var <- Upattern$var[iPattern]
            iNtime <- Upattern$n.time[iPattern]
            iX.var <- X.var[[iPattern.var]]
            iOmega.sd <- attr(Omega[[iPattern]], "sd")
            idOmega.sd <- dFCT.sigma(p = param[name.sigma], n.time = iNtime, X = iX.var)

            ## derivative of rho with respect to the correlation parameters
            if(iNtime > 1 && !is.na(Upattern$cor[iPattern])){
                iPattern.cor <- Upattern$cor[iPattern]
                iX.cor <- X.cor[[iPattern.cor]]
                iOmega.cor <- attr(Omega[[iPattern]], "cor")
                idOmega.cor <- dFCT.rho(p = param[name.rho], n.time = iNtime, X = iX.cor)
            }

            ## derivative of Omega with respect to the variance and correlation parameters
            if(iNtime > 1 && !is.na(Upattern$cor[iPattern])){
                iOut <- c(
                    lapply(idOmega.sd, function(iDeriv){
                        iDeriv <- unname(iDeriv)
                        return(diag(2*iDeriv*iOmega.sd, nrow = iNtime, ncol = iNtime) + iOmega.cor * (iDeriv %*% t(iOmega.sd) + iOmega.sd %*% t(iDeriv)))
                    }),
                    lapply(idOmega.cor, function(iDeriv){
                        iDeriv <- unname(iDeriv)
                        return(iDeriv * tcrossprod(iOmega.sd))
                    })
                )
            }else{
                iOut <- lapply(idOmega.sd, function(iDeriv){
                    diag(2*as.double(iDeriv)*as.double(iOmega.sd), nrow = iNtime, ncol = iNtime)
                })
            }
        
            return(iOut)
        }), Upattern$name)
    }

    ## ** Jacobian: reparametrize derivative
    if(!is.null(Jacobian)){
        type <- object$param$type
        name.sigma <- object$param$name[type=="sigma"]
        name.k <- object$param$name[type=="k"]
        name.rho <- object$param$name[type=="rho"]
        name.paramVar <- c(name.sigma,name.k,name.rho)

        out <- lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

            iPattern.var <- Upattern[iPattern,"var"]
            iPattern.cor <- Upattern[iPattern,"cor"]
            iNtime <- Upattern[iPattern,"n.time"]
            iName.param <- Upattern[iPattern,"param"][[1]]

            ## [dOmega_[11]/d theta_1] ... [dOmega_[11]/d theta_p] %*% Jacobian
            ## [dOmega_[ij]/d theta_1] ... [dOmega_[ij]/d theta_p] %*% Jacobian
            ## [dOmega_[mm]/d theta_1] ... [dOmega_[mm]/d theta_p] %*% Jacobian
            if(any(abs(Jacobian[iName.param,setdiff(name.paramVar,iName.param),drop=FALSE])>1e-10)){
                stop("Something went wrong when computing the derivative of the residual variance covariance matrix. \n",
                     "Contact the package manager with a reproducible example generating this error message. \n")
            }
            M.iScore <- do.call(cbind,lapply(out[[iPattern]],as.double)) %*% Jacobian[iName.param,iName.param,drop=FALSE]
            iOut <- stats::setNames(lapply(1:NCOL(M.iScore), function(iCol){matrix(M.iScore[,iCol], nrow = iNtime, ncol = iNtime, byrow = FALSE)}), iName.param)
            return(iOut)
        })
        
        out <- stats::setNames(out,Upattern$name)
    }

    ## ** export
    return(out)
}

## * helper
## ** .dsigma_transform
##' @param value [numeric] parameter value after transformation.
##' @param Omega.sd [numeric vector] square root of the diagonal of the residual variance-covariance matrix.
##' @param transformation [character] transformation for the sigma parameter: \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}.
##' @param power [positive integer] order of the derivative.
.dsigma_transform <- function(value, Omega.sd, transform, power = 1){

    ## d sigma k / d sigma = k = sigma k / sigma
    if(transform == "none"){
        if(power == 1){
            out <- Omega.sd / value
        }else{
            out <- rep(0, length(Omega.sd))
        }
    }else if(transform == "log"){
        ## d exp(log(sigma))  k / d log(sigma) = exp(log(sigma)) k = sigma k
        out <- Omega.sd
    }else if(transform == "square"){
        ## d sqrt(sigma^2)  k / d sigma^2 =  k / (2 sqrt(sigma^2)) =  k sigma / (2 sigma^2)
        out <- Omega.sd * prod(1/2 - seq(from = 0, to = power-1, by = 1)) / (value^(power)) ## here value = sigma^2 due to transformation
    }else if(transform == "logsquare"){
        ## d exp(0.5 log(sigma^2)) k1 / d log(sigma^2) = 0.5 k1 exp(0.5 log(sigma^2)) = 0.5 k1 exp(log(sigma)) = 0.5 sigma k1
        out <- Omega.sd / 2^(power)
    }else{
        stop("Omega derivative cannot handle transformation ",transform," for sigma parameters. \n")
    }

    return(out)
}

## ** .dk_transform
##' @param value [numeric] parameter value after transformation.
##' @param indicator [logical vector] TRUE where the k parameter is and FALSE otherwise.
##' For instance for [\sigma \sigma*k_1 \sigma*k_2 \sigma*k_3] it would be [FALSE TRUE FALSE FALSE] for k_1
##' @param Omega.sd [numeric vector] square root of the diagonal of the residual variance-covariance matrix.
##' @param transformation [character] transformation for the k parameter: \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"}.
##' @param power [positive integer] order of the derivative.
.dk_transform <- function(value, indicator, Omega.sd, transform, power = 1){
    
    if(transform == "none"){
        ## d sigma k / d k = sigma = sigma k / k
        if(power == 1){
            out <- indicator * Omega.sd / value
        }else{
            out <- rep(0, length(Omega.sd))
        }
    }else if(transform == "log"){
        ## d sigma exp(log(k)) / d log(k) = sigma exp(log(k)) = sigma k
        out <- indicator * Omega.sd
    }else if(transform == "square"){
        ## d sigma sqrt(k^2) / d k^2 = sigma / (2 sqrt(k^2)) = sigma / (2 k) = sigma k / (2 k^2)
        ## d^2 sigma sqrt(k^2) / (d k^2)^2 =  d sigma / (2 sqrt(k^2))/ d k^2 = - sigma / (4 k^2 sqrt(k^2)) = - sigma k / (4 k^4)
        out <- indicator * Omega.sd * prod(1/2 - seq(from = 0, to = power-1, by = 1)) / value^(power) ## here value = k^2 due to transformation
    }else if(transform == "logsquare"){
        ## d sigma exp(0.5 log(k^2)) / d log(k^2) = 0.5 sigma exp(0.5 log(k^2)) = 0.5 sigma k
        out <- indicator * Omega.sd / 2^(power)
    }else if(transform == "sd"){
        ## d sigma k / d sigma k = 1
        if(power == 1){
            out <- indicator
        }else{
            out <- rep(0, length(Omega.sd))
        }
    }else if(transform == "logsd"){
        ## d exp(log(sigma k)) / d log(sigma k) = exp(log(sigma k)) = sigma k
        out <- indicator * Omega.sd
    }else if(transform == "var"){
        ## d sqrt((sigma k)^2) / d (sigma k)^2 = 1/(2 sqrt((sigma k)^2)) =  1/(2 sigma k)
        ## d^2 sqrt(sigma^2 k^2) / (d sigma^2 k^2)^2 =  d 1 / (2 sqrt(sigma^2 k^2))/ d sigma^2 k^2 = - 1 / (4 sigma^2 k^2 sqrt(sigma^2 k^2)) = - 1 / (4 sigma^3 k^3)
        out <- indicator * prod(1/2 - seq(from = 0, to = power-1, by = 1)) / Omega.sd^(2*power-1) ## here Omega.sd = sigma k as no transformation is applied to Omega.sd
    }else if(transform == "logvar"){
        ## d exp(0.5*log((sigma k)^2) / d log((sigma k)^2) = 0.5 exp(0.5*log((sigma k)^2)  =  0.5 sigma k
        out <- indicator * Omega.sd / 2^(power)
    }else{
        stop("Omega derivative cannot handle transformation ",transform," for k parameters. \n")
    }

    return(out)
}

## ** .drho_transform
##' @param value [numeric] parameter value after transformation.
##' @param indicator [logical matrix] TRUE where the rho parameter is in the residual variance-covariance matrix and FALSE otherwise.
##' @param Omega [numeric matrix] residual variance-covariance matrix.
##' @param transformation [character] transformation for the rho parameter: \code{"none"}, \code{"atanh"}, \code{"cov"}.
##' @param power [positive integer] order of the derivative.
.drho_transform <- function(value, indicator, Omega = NULL, transform, power = 1){

    if(transform == "none"){
        ## d rho / d rho = 1
        if(power == 1){
            if(is.null(Omega)){
                out <- indicator        
            }else{
                out <- indicator * Omega / value        
            }
        }else{
            out <- matrix(0, nrow = NROW(indicator), ncol = NCOL(indicator))
        }
    }else if(transform == "atanh"){
        ## d  tanh(atanh(rho)) / d atanh(rho) =  (1-tanh(atanh(rho))^2) =  (1-rho^2)
        ## d^2 tanh(atanh(rho)) / (d atanh(rho))^2 = d(1-tanh(atanh(rho))^2)/d atanh(rho)
        ##                                         = - 2 tanh(atanh(rho)) d tanh(atanh(rho))/d atanh(rho)
        ##                                         = - 2 tanh(atanh(rho)) (1 - tanh(atanh(rho))^2)
        if(power == 1){
            out <- indicator * (1 - tanh(value)^2)  ## here value = atanh(rho) due to transformation
        }else if(power == 2){
            out <- - indicator  * 2 * tanh(value) * (1 - tanh(value)^2)
        }else if(power == 3){
            out <- indicator *  (- 2 + 8 * tanh(value)^2 - 6 * tanh(value)^4)
        }else{
            stop("Derivative w.r.t. rho only implemented for order 1, 2, and 3 but not for order ",power,". \n")
        }
        if(!is.null(Omega)){
            out <- out * Omega / tanh(value)
        }
    }else if(transform == "cov"){
        if(power == 1){
            if(is.null(Omega)){
                warning("Derivative is taken w.r.t. the residual variance-covariance matrix, not the residual correlation matrix.\n")
            }
            out <- indicator
        }else{
            out <- matrix(0, nrow = NROW(indicator), ncol = NCOL(indicator))
        }
    }else{
        stop("Omega derivative cannot handle transformation ",transform," for rho parameters. \n")
    }

    return(out)
}

##----------------------------------------------------------------------
### calc_dOmega.R ends here
