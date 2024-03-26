### calc_dOmega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:18) 
## Version: 
## Last-Updated: Mar 24 2024 (14:58) 
##           By: Brice Ozenne
##     Update #: 222
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
##' @param param [named numeric vector] values of the parameters (untransformed).
##' @param Omega [list of matrices] Residual Variance-Covariance Matrix for each pattern.
##' @param Jacobian [matrix] Jacobian of the reparametrisation.
##' @param transform.sigma,transform.k,transform.rho [character] Transformation used on the variance/correlation coefficients.
##' Only active if \code{"log"}, \code{"log"}, \code{"atanh"}: then the derivative is directly computed on the transformation scale instead of using the Jacobian.
##' @keywords internal
##' 
##' @examples
##' data(gastricbypassL, package = "LMMstar")
##' gastricbypassL$gender <- c("M","F")[as.numeric(gastricbypassL$id) %% 2+1]
##' dd <- gastricbypassL[!duplicated(gastricbypassL[,c("time","gender")]),]
##' 
##' ## independence
##' Sid1 <- .skeleton(IND(~1, var.time = "time"), data = dd)
##' Sid4 <- .skeleton(IND(~1|id, var.time = "time"), data = dd)
##' Sdiag1 <- .skeleton(IND(~visit), data = dd)
##' Sdiag4 <- .skeleton(IND(~visit|id), data = dd)
##' Sdiag24 <- .skeleton(IND(~visit+gender|id, var.time = "time"), data = dd)
##' param24 <- setNames(c(1,2,2,3,3,3,5,3),Sdiag24$param$name)
##'
##' .calc_dOmega(Sid4, param = c(sigma = 2))
##' .calc_dOmega(Sdiag1, param = c(sigma = 1, k.visit2 = 2, k.visit3 = 3, k.visit4 = 4))
##' .calc_dOmega(Sdiag4, param = c(sigma = 1, k.visit2 = 2, k.visit3 = 3, k.visit4 = 4))
##' .calc_dOmega(Sdiag24, param = param24)
##' 
##' ## compound symmetry
##' Scs4 <- .skeleton(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' Scs24 <- .skeleton(CS(gender~time|id), data = gastricbypassL)
##' 
##' .calc_dOmega(Scs4, param = c(sigma = 1,rho=0.5))
##' .calc_dOmega(Scs4, param = c(sigma = 2,rho=0.5))
##' .calc_dOmega(Scs24, param = c("sigma:F" = 2, "sigma:M" = 1,
##'                             "rho:F"=0.5, "rho:M"=0.25))
##' 
##' ## unstructured
##' Sun4 <- .skeleton(UN(~visit|id), data = gastricbypassL)
##' param4 <- setNames(c(1,1.1,1.2,1.3,0.5,0.45,0.55,0.7,0.1,0.2),Sun4$param$name)
##' Sun24 <- .skeleton(UN(gender~visit|id), data = gastricbypassL)
##' param24 <- setNames(c(param4,param4*1.1),Sun24$param$name)
##' 
##' .calc_dOmega(Sun4, param = param4)
##' .calc_dOmega(Sun24, param = param24)
`.calc_dOmega` <-
    function(object, param, Omega, Jacobian, transform.sigma, transform.k, transform.rho) UseMethod(".calc_dOmega")

## * calc_dOmega.ID
.calc_dOmega.ID <- function(object, param, Omega, Jacobian = NULL,
                            transform.sigma = NULL, transform.k = NULL, transform.rho = NULL){

    ## ** prepare
    type <- object$param$type
    name.sigma <- object$param$name[type=="sigma"]
    name.k <- object$param$name[type=="k"]
    name.rho <- object$param$name[type=="rho"]
    name.paramVar <- c(name.sigma,name.k,name.rho)
    param <- param[name.paramVar]
    name.param <- names(param)
    if(missing(Omega)){
        Omega <- .calc_Omega(object, param = param, keep.interim = TRUE)
    }
    
    Upattern <- object$Upattern
    n.Upattern <- NROW(Upattern)
    pattern.cluster <- object$pattern
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern

    if(identical(transform.sigma,"log") && identical(transform.k,"log") && identical(transform.rho,"atanh")){
        Jacobian <- NULL
    }else{
        transform.sigma <- "none"
        transform.k <- "none"
        transform.rho <- "none"
    }
    log.param <- log(param[c(name.sigma,name.k)])

    ## ** loop over covariance patterns
    out <- lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

        iPattern.var <- Upattern[iPattern,"var"]
        iPattern.cor <- Upattern[iPattern,"cor"]
        iNtime <- Upattern[iPattern,"n.time"]
        iName.param <- Upattern[iPattern,"param"][[1]]
        if(is.null(iName.param)){return(NULL)}

        iOmega.sd <- attr(Omega[[iPattern]],"sd")
        iOmega.cor <- attr(Omega[[iPattern]],"cor")
        iOmega <- Omega[[iPattern]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL; attr(iOmega,"time") <- NULL;
        iParam.sigma <- intersect(iName.param, name.sigma)
        n.iParam.sigma <- length(iParam.sigma)
        iParam.k <- intersect(iName.param, name.k)
        n.iParam.k <- length(iParam.k)
        iParam.rho <- intersect(iName.param, name.rho)
        n.iParam.rho <- length(iParam.rho)

        iScore <- stats::setNames(vector(mode = "list", length = length(iName.param)), iName.param)
        iIndicator.cor <- attr(X.cor[[iPattern.cor]],"indicator.param")
        iMindicator.var <- attr(X.var[[iPattern.var]],"Mindicator.param")

        if(n.iParam.sigma>0){
            if(transform.sigma == "log"){
                iScore[[iParam.sigma]] <- 2 * iOmega
            }else{ ## no transformation  (other transformations are made through jacobian)
                iScore[[iParam.sigma]] <- 2 * iOmega / param[iParam.sigma]
            }
        }

        if(n.iParam.k>0){
            if(transform.k == "log"){
                iScore[names(iMindicator.var)] <- lapply(iMindicator.var, function(iM){iM * iOmega})
            }else{ ## no transformation  (other transformations are made through jacobian)
                iScore[names(iMindicator.var)] <- lapply(names(iMindicator.var), function(iParam){iMindicator.var[[iParam]] * iOmega / param[iParam]})
            }
        }

        if(n.iParam.rho>0){
            iOmega.var <- tcrossprod(iOmega.sd)
            
            for(iRho in iParam.rho){ ## iRho <- iParam.rho[1]
                iScore[[iRho]] <- diag(0, nrow = iNtime, ncol = iNtime)
                if(transform.rho == "atanh"){
                    iScore[[iRho]][iIndicator.cor[[iRho]]] <- iOmega.var[iIndicator.cor[[iRho]]] * (1-param[iRho]^2)
                }else{ ## no transformation (other transformations are made through jacobian)
                    iScore[[iRho]][iIndicator.cor[[iRho]]] <- iOmega.var[iIndicator.cor[[iRho]]]
                }
            }

        }
        ## apply transformation
        if(!is.null(Jacobian)){
            ## [dOmega_[11]/d theta_1] ... [dOmega_[11]/d theta_p] %*% Jacobian
            ## [dOmega_[ij]/d theta_1] ... [dOmega_[ij]/d theta_p] %*% Jacobian
            ## [dOmega_[mm]/d theta_1] ... [dOmega_[mm]/d theta_p] %*% Jacobian
            if(any(abs(Jacobian[iName.param,setdiff(name.paramVar,iName.param),drop=FALSE])>1e-10)){
                stop("Something went wrong when computing the derivative of the residual variance covariance matrix. \n",
                     "Contact the package manager with a reproducible example generating this error message. \n")
            }
            M.iScore <- do.call(cbind,lapply(iScore,as.double)) %*% Jacobian[iName.param,iName.param,drop=FALSE]
            iScore <- stats::setNames(lapply(1:NCOL(M.iScore), function(iCol){matrix(M.iScore[,iCol], nrow = iNtime, ncol = iNtime, byrow = FALSE)}), iName.param)
        }
        return(iScore)
    })

    ## ** export
    out <- stats::setNames(out,Upattern$name)
    return(out)
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

    if(!is.null(FCT.sigma) && is.null(dFCT.sigma) || !is.null(FCT.rho) && is.null(dFCT.rho) ){

        ## unlist(.calc_Omega.CUSTOM(object, param = param, keep.interim = FALSE))
        vec.dOmega <- numDeriv::jacobian(func = function(x){
            unlist(.calc_Omega.CUSTOM(object, param = x, keep.interim = FALSE))
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
    return(out)
}


##----------------------------------------------------------------------
### calc_dOmega.R ends here
