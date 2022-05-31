### calc_Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 21 2021 (18:12) 
## Version: 
## Last-Updated: May 30 2022 (23:35) 
##           By: Brice Ozenne
##     Update #: 487
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
    function(object, param, keep.interim) UseMethod(".calc_Omega")


## * calc_Omega.ID
.calc_Omega.ID <- function(object, param, keep.interim = FALSE){

    Upattern <- object$X$Upattern
    n.Upattern <- NROW(Upattern)
    pattern.cluster <- object$X$pattern.cluster
    X.var <- object$X$Xpattern.var
    X.cor <- object$X$Xpattern.cor

    Omega <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1
        iPattern.var <- Upattern[iPattern,"var"]
        iPattern.cor <- Upattern[iPattern,"cor"]
        iNtime <- Upattern[iPattern,"n.time"]
        Omega.sd <- unname(exp(X.var[[iPattern.var]] %*% log(param[colnames(X.var[[iPattern.var]])])))
        Omega.cor <- diag(0, nrow = iNtime, ncol = iNtime)
        if(!is.null(X.cor) && !is.null(X.cor[[iPattern.cor]])){
            iParam.cor <- attr(X.cor[[iPattern.cor]],"param")
            for(iiP in 1:length(iParam.cor)){
                iiParam <- iParam.cor[iiP]
                if(is.na(iiParam)){
                    iiParam <- which(is.na(names(attr(X.cor[[iPattern.cor]],"indicator.param")))) 
                    Omega.cor[attr(X.cor[[iPattern.cor]],"indicator.param")[[iiParam]]] <- NA
                }else{
                    Omega.cor[attr(X.cor[[iPattern.cor]],"indicator.param")[[iiParam]]] <- param[iiParam]
                }
            }
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
}

## * calc_Omega.IND
.calc_Omega.IND <- .calc_Omega.ID

## * calc_Omega.CS
.calc_Omega.CS <- .calc_Omega.ID

## * calc_Omega.UN
.calc_Omega.UN <- .calc_Omega.ID

## * calc_Omega.CS
.calc_Omega.CUSTOM <- function(object, param, keep.interim = FALSE){

    Upattern <- object$X$Upattern
    n.Upattern <- NROW(Upattern)
    pattern.cluster <- object$X$pattern.cluster
    X.var <- object$X$var
    X.cor <- object$X$cor
    FCT.sigma <- object$FCT.sigma
    FCT.rho <- object$FCT.rho
    name.sigma <- names(object$init.sigma)
    name.rho <- names(object$init.rho)

    Omega <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

        iPattern.var <- object$X$Upattern$var[iPattern]
        iNtime <- object$X$Upattern$n.time[iPattern]
        iX.var <- object$X$Xpattern.var[[iPattern.var]]
        iTime <- attr(iX.var, "index.time")
        iOmega.sd <- FCT.sigma(p = param[name.sigma], time = iTime, X = iX.var)
        
        if(iNtime > 1 && !is.null(X.cor)){
            iPattern.cor <- object$X$Upattern$cor[iPattern]
            iX.cor <- object$X$Xpattern[[iPattern.cor]]
            iOmega.cor <- FCT.rho(p = param[name.rho], time = iTime, X = iX.var)
            diag(iOmega.cor) <- 0
        }
        iOmega <- diag(as.double(iOmega.sd)^2, nrow = iNtime, ncol = iNtime) + iOmega.cor * tcrossprod(iOmega.sd)
        
        if(keep.interim){
            attr(iOmega,"sd") <- iOmega.sd
            attr(iOmega,"cor") <- iOmega.cor
            attr(iOmega,"time") <- iTime
        }
        return(iOmega)
    }), Upattern$name)

    return(Omega)
}


##----------------------------------------------------------------------
### calc_Omega.R ends here
