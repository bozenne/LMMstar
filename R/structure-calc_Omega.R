### calc_Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 21 2021 (18:12) 
## Version: 
## Last-Updated: okt  1 2021 (17:07) 
##           By: Brice Ozenne
##     Update #: 412
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


## * calc_Omega.UN
.calc_Omega.UN <- function(object, param, keep.interim = FALSE){
    Upattern <- object$X$Upattern
    n.Upattern <- NROW(Upattern)
    pattern.cluster <- object$X$pattern.cluster
    X.var <- object$X$var
    X.cor <- object$X$cor
    
    Omega <- stats::setNames(lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1
        iPattern.var <- Upattern[iPattern,"var"]
        iPattern.cor <- Upattern[iPattern,"cor"]
        iTime <- Upattern[iPattern,"time"][[1]]
        iNtime <- length(iTime)

        Omega.sd <- unname(exp(X.var[[iPattern.var]] %*% log(param[colnames(X.var[[iPattern.var]])])))
        Omega.cor <- diag(0, nrow = iNtime, ncol = iNtime)
        if(!is.null(X.cor) && !is.null(X.cor[[iPattern.cor]])){
            Omega.cor[attr(X.cor[[iPattern.cor]],"index.vec2matrix")] <- X.cor[[iPattern.cor]] %*% param[colnames(X.cor[[iPattern.cor]])]
        }
        Omega <- diag(as.double(Omega.sd)^2, nrow = iNtime, ncol = iNtime) + Omega.cor * tcrossprod(Omega.sd)
        
        if(keep.interim){
            attr(Omega,"time") <- iTime
            attr(Omega,"sd") <- Omega.sd
            attr(Omega,"cor") <- Omega.cor
        }
        return(Omega)
    }), Upattern$name)

    return(Omega)
}

## * calc_Omega.IND
.calc_Omega.IND <- .calc_Omega.UN

## * calc_Omega.CS
.calc_Omega.CS <- .calc_Omega.UN



##----------------------------------------------------------------------
### calc_Omega.R ends here
