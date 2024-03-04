### calc_Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 21 2021 (18:12) 
## Version: 
## Last-Updated: feb 28 2024 (13:55) 
##           By: Brice Ozenne
##     Update #: 613
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calc_Omega
##' @title Construct Residual Variance-Covariance Matrix
##' @description Construct residual variance-covariance matrix for given parameter values.
##' @noRd
##'
##' @param structure [structure]
##' @param param [named numeric vector] values of the parameters.
##' @param keep.interim [logical] should the correlation matrix and variance matrix be output
##' @param Upattern [data.frame] Optional, used to only evaluate the residual variance-covariance with respect to a subset of patterns.
##' Should contain the name of the pattern, the index of the variance pattern, the index of the correlation pattern.
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
##' .calc_Omega(Sid1, param = c(sigma = 1))
##' .calc_Omega(Sid4, param = c(sigma = 2))
##' .calc_Omega(Sdiag1, param = c(sigma = 1, k.visit2 = 2, k.visit3 = 3, k.visit4 = 4))
##' .calc_Omega(Sdiag4, param = c(sigma = 1, k.visit2 = 2, k.visit3 = 3, k.visit4 = 4))
##' .calc_Omega(Sdiag24, param = param24)
##' 
##' ## compound symmetry
##' Scs4 <- skeleton(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' Scs24 <- skeleton(CS(gender~time|id), data = gastricbypassL)
##' 
##' .calc_Omega(Scs4, param = c(sigma = 1,rho=0.5))
##' .calc_Omega(Scs4, param = c(sigma = 2,rho=0.5))
##' .calc_Omega(Scs24, param = c("sigma:F" = 2, "sigma:M" = 1,
##'                             "rho:F"=0.5, "rho:M"=0.25))
##' 
##' ## unstructured
##' Sun4 <- skeleton(UN(~visit|id), data = gastricbypassL)
##' param4 <- setNames(c(1,1.1,1.2,1.3,0.5,0.45,0.55,0.7,0.1,0.2),Sun4$param$name)
##' Sun24 <- skeleton(UN(gender~visit|id), data = gastricbypassL)
##' param24 <- setNames(c(param4,param4*1.1),Sun24$param$name)
##' 
##' .calc_Omega(Sun4, param = param4)
##' .calc_Omega(Sun24, param = param24, keep.interim = TRUE)
`.calc_Omega` <-
    function(object, param, keep.interim, Upattern) UseMethod(".calc_Omega")


## * calc_Omega.ID
.calc_Omega.ID <- function(object, param, keep.interim = FALSE, Upattern = NULL){

    if(is.null(Upattern)){
        Upattern <- object$Upattern
    }
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern

    Omega <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1
        iPattern.var <- Upattern[iPattern,"var"]
        iPattern.cor <- Upattern[iPattern,"cor"]
        iNtime <- Upattern[iPattern,"n.time"]

        if(length(X.var[[iPattern.var]])>0){
            Omega.sd <- unname(exp(X.var[[iPattern.var]] %*% log(param[colnames(X.var[[iPattern.var]])])))
        }else{
            Omega.sd <- rep(1, iNtime)
        }
        if(!is.null(X.cor) && !is.null(X.cor[[iPattern.cor]])){            
            Omega.cor <- attr(X.cor[[iPattern.cor]],"Omega.cor")
            iParam.cor <- attr(X.cor[[iPattern.cor]],"param")
            for(iiP in 1:length(iParam.cor)){ ## iiP <- 3
                Omega.cor[attr(X.cor[[iPattern.cor]],"indicator.param")[[iParam.cor[iiP]]]] <- param[iParam.cor[iiP]]
            }
            Omega <- Omega.cor * tcrossprod(Omega.sd)
        }else{
            Omega.cor <- NULL            
            Omega <- diag(as.double(Omega.sd)^2, nrow = iNtime, ncol = iNtime)
        }
        if(keep.interim){
            attr(Omega,"sd") <- Omega.sd
            attr(Omega,"cor") <- Omega.cor
            attr(Omega,"time") <- attr(X.var[[iPattern.var]], "index.time")
        }
        return(Omega)
    }), Upattern$name)

    ## print(Omega)
    return(Omega)
}

## * calc_Omega.IND
.calc_Omega.IND <- .calc_Omega.ID

## * calc_Omega.CS
.calc_Omega.CS <- .calc_Omega.ID

## * calc_Omega.RE
.calc_Omega.RE <- .calc_Omega.ID

## * calc_Omega.TOEPLITZ
.calc_Omega.TOEPLITZ <- .calc_Omega.ID

## * calc_Omega.UN
.calc_Omega.UN <- .calc_Omega.ID

## * calc_Omega.EXP
.calc_Omega.EXP <- function(object, param, keep.interim = FALSE, Upattern = NULL){

    if(is.null(Upattern)){
        Upattern <- object$Upattern
    }
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern
    regressor <- stats::setNames(object$param[object$param$type=="rho","code"],object$param[object$param$type=="rho","name"])
    
    Omega <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1
        iPattern.var <- Upattern[iPattern,"var"]
        iPattern.cor <- Upattern[iPattern,"cor"]
        iNtime <- Upattern[iPattern,"n.time"]

        if(length(X.var[[iPattern.var]])>0){
            Omega.sd <- unname(exp(X.var[[iPattern.var]] %*% log(param[colnames(X.var[[iPattern.var]])])))
        }else{
            Omega.sd <- rep(1, iNtime)
        }
        Omega.cor <- diag(0, nrow = iNtime, ncol = iNtime)
        
        if(!is.null(X.cor) && !is.null(X.cor[[iPattern.cor]])){
            iParam.cor <- attr(X.cor[[iPattern.cor]],"param")
            iTime.cor <- X.cor[[iPattern.cor]][,regressor[iParam.cor]]
            Omega.cor[attr(X.cor[[iPattern.cor]],"indicator.param")[[iParam.cor]]] <- exp(-param[iParam.cor] * iTime.cor)
        }
        Omega <- diag(as.double(Omega.sd)^2, nrow = iNtime, ncol = iNtime) + Omega.cor * tcrossprod(Omega.sd)
        
        if(keep.interim){
            attr(Omega,"sd") <- Omega.sd
            attr(Omega,"cor") <- Omega.cor
            attr(Omega,"time") <- attr(X.var[[iPattern.var]], "index.time")
        }
        return(Omega)
    }), Upattern$name)
    ## print(Omega)

    return(Omega)
    
    return(1)
}
## * calc_Omega.CUSTOM
.calc_Omega.CUSTOM <- function(object, param, keep.interim = FALSE, Upattern = NULL){

    if(is.null(Upattern)){
        Upattern <- object$Upattern
    }
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern
    FCT.sigma <- object$FCT.sigma
    FCT.rho <- object$FCT.rho
    name.sigma <- object$param[object$param$type=="sigma","name"]
    name.rho <- object$param[object$param$type=="rho","name"]

    Omega <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

        iPattern.var <- Upattern$var[iPattern]
        iNtime <- Upattern$n.time[iPattern]
        iX.var <- X.var[[iPattern.var]]
        iOmega.sd <- FCT.sigma(p = param[name.sigma], n.time = iNtime, X = iX.var)

        if(iNtime > 1 && !is.na(Upattern$cor[iPattern])){
            iPattern.cor <- Upattern$cor[iPattern]
            iX.cor <- X.cor[[iPattern.var]]
            iOmega.cor <- FCT.rho(p = param[name.rho], n.time = iNtime, X = iX.cor)
            diag(iOmega.cor) <- 0
            iOmega <- diag(as.double(iOmega.sd)^2, nrow = iNtime, ncol = iNtime) + iOmega.cor * tcrossprod(iOmega.sd)
        }else{
            iOmega.cor <- NULL
            iOmega <- diag(as.double(iOmega.sd)^2, nrow = iNtime, ncol = iNtime)
        }
        
        if(keep.interim){
            attr(iOmega,"sd") <- iOmega.sd
            attr(iOmega,"cor") <- iOmega.cor
        }
        return(iOmega)
    }), Upattern$name)

    return(Omega)
}


##----------------------------------------------------------------------
### calc_Omega.R ends here
