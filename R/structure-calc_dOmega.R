### calc_dOmega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:18) 
## Version: 
## Last-Updated: okt  1 2021 (17:07) 
##           By: Brice Ozenne
##     Update #: 53
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
##' @param param [named numeric vector] values of the parameters.
##' @param Omega [list of matrices] Residual Variance-Covariance Matrix for each pattern.
##' @param Jacobian [matrix] Jacobian of the reparametrisation.
##'
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
    function(object, param, Omega, Jacobian) UseMethod(".calc_dOmega")

## * calc_Omega.UN
.calc_dOmega.UN <- function(object, param, Omega, Jacobian = NULL){

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
    
    Upattern <- object$X$Upattern
    n.Upattern <- NROW(Upattern)
    pattern.cluster <- object$X$pattern.cluster
    X.var <- object$X$var
    X.cor <- object$X$cor

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
        iIndicator <- attr(X.cor[[iPattern.cor]],"indicator.param")

        if(n.iParam.sigma==1){
            ## compute derivative
            iParam.dsigma <- matrix(param[c(iParam.sigma,iParam.k)], nrow = n.iParam.sigma+n.iParam.k, ncol = n.iParam.sigma, byrow = FALSE,
                                    dimnames = list(c(iParam.sigma,iParam.k),iParam.sigma))
            iParam.dsigma[iParam.sigma,iParam.sigma] <- 1
            idOmega.sigma <- exp(iX.var %*% log(iParam.dsigma))
            ## propage
            iScore[[iParam.sigma]] <- unname(diag(2*as.double(idOmega.sigma)*as.double(iOmega.sd), nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.sigma %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.sigma)))
            ## iScore[[iParam.sigma]] - 2*iOmega/param[iParam.sigma]
        }else{
            stop("No sigma parameter in the structure. Something is wrong. \n")
        }
        
        if(n.iParam.k>0){
            ## compute derivative
            iParam.dk <- matrix(param[c(iParam.sigma,iParam.k)], nrow = n.iParam.sigma+n.iParam.k, ncol = n.iParam.k, byrow = FALSE,
                                dimnames = list(c(iParam.sigma,iParam.k),iParam.k))
            if(n.iParam.k==1){
                iParam.dk[iParam.k,iParam.k] <- 1
            }else{
                diag(iParam.dk[iParam.k,iParam.k]) <- 1
            }
            idOmega.k <- exp(iX.var %*% log(iParam.dk)) * iX.var[,iParam.k]
            for(iK in iParam.k){ ## iK <- iParam.k[1]
                iScore[[iK]] <- diag(2*as.double(idOmega.k[,iK])*as.double(iOmega.sd), nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.k[,iK] %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k[,iK]))
                ## iIndicator <- tcrossprod(X.var[[iPattern.var]][,iK], rep(1,NROW(X.var[[iPattern.var]][,iK]))) + t(tcrossprod(X.var[[iPattern.var]][,iK], rep(1,NROW(X.var[[iPattern.var]][,iK])))) > 0
                ## iScore[[iK]] - iOmega/param[iK] * (iIndicator + tcrossprod(X.var[[iPattern.var]][,iK]))
            }
        }
        if(n.iParam.rho>0){
            iOmega.var <- tcrossprod(iOmega.sd)
            
            for(iRho in iParam.rho){ ## iRho <- iParam.rho[1]
                iScore[[iRho]] <- diag(0, nrow = iNtime, ncol = iNtime)
                iScore[[iRho]][iIndicator[[iRho]]] <- iOmega.var[iIndicator[[iRho]]]

                ## iIndicator <- diag(0, nrow = iNtime, ncol = iNtime)
                ## iIndicator[attr(iX.cor,"index.vec2matrix")] <- iX.cor[,iRho]
                ## ## derivative
                ## iScore[[iRho]] <- iIndicator * tcrossprod(iOmega.sd)

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
    out <- stats::setNames(out,Upattern$name)
    return(out)
}

## * calc_Omega.IND
.calc_dOmega.IND <- .calc_dOmega.UN

## * calc_Omega.CS
.calc_dOmega.CS <- .calc_dOmega.UN



##----------------------------------------------------------------------
### calc_dOmega.R ends here
