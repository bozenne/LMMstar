### structure-calc_d3Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  4 2023 (17:08) 
## Version: 
## Last-Updated: aug  4 2023 (17:35) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calc_d3Omega
##' @title Third Derivative of the Residual Variance-Covariance Matrix
##' @description Third derivative of the residual variance-covariance matrix for given parameter values.
##' @noRd
##'
##' @param structure [structure]
##' @param param [named numeric vector] values of the parameters.
##' @param Omega [list of matrices] Residual Variance-Covariance Matrix for each pattern.
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
##' .calc_d3Omega(Sid4, param = c(sigma = 2))
##' .calc_d3Omega(Sdiag1, param = c(sigma = 1, k.visit2 = 2, k.visit3 = 3, k.visit4 = 4))
##' .calc_d3Omega(Sdiag4, param = c(sigma = 1, k.visit2 = 2, k.visit3 = 3, k.visit4 = 4))
##' .calc_d3Omega(Sdiag24, param = param24)
##' 
##' ## compound symmetry
##' Scs4 <- skeleton(CS(~1|id, var.time = "time"), data = gastricbypassL)
##' Scs24 <- skeleton(CS(gender~time|id), data = gastricbypassL)
##' 
##' .calc_d3Omega(Scs4, param = c(sigma = 1,rho=0.5))
##' .calc_d3Omega(Scs4, param = c(sigma = 2,rho=0.5))
##' .calc_d3Omega(Scs24, param = c("sigma:F" = 2, "sigma:M" = 1,
##'                             "rho:F"=0.5, "rho:M"=0.25))
##' 
##' ## unstructured
##' Sun4 <- skeleton(UN(~visit|id), data = gastricbypassL)
##' param4 <- setNames(c(1,1.1,1.2,1.3,0.5,0.45,0.55,0.7,0.1,0.2),Sun4$param$name)
##' Sun24 <- skeleton(UN(gender~visit|id), data = gastricbypassL)
##' param24 <- setNames(c(param4,param4*1.1),Sun24$param$name)
##' 
##' .calc_d3Omega(Sun4, param = param4)
##' .calc_d3Omega(Sun24, param = param24)
`.calc_d3Omega` <-
    function(object, param, Omega,
             transform.sigma, transform.k, transform.rho) UseMethod(".calc_d3Omega")

## * calc_d3Omega.ID
.calc_d3Omega.ID <- function(object, param, Omega,
                              transform.sigma = NULL, transform.k = NULL, transform.rho = NULL){

    ## ** prepare
    if(missing(Omega)){
        Omega <- .calc_Omega(object, param = param, keep.interim = TRUE)
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

    triplet.vcov <- object$triplet.vcov

    if(identical(transform.sigma,"log") && identical(transform.k,"log") && identical(transform.rho,"atanh")){
        transform <- TRUE
    }else if(identical(transform.sigma,"none") && identical(transform.k,"none") && identical(transform.rho,"none")){
        transform <- FALSE
    }else{
        stop("Method calc_d3Omega not available for the requested transformation. \n",
             "Only implemented for not transformation or log(sigma),log(k),atanh(rho). \n")
    }

    ## ** loop over covariance patterns
    out <- lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

        iPattern.var <- Upattern[iPattern,"var"]
        iPattern.cor <- Upattern[iPattern,"cor"]
        iName.param <- Upattern[iPattern,"param"][[1]]
        iN.time <- Upattern[iPattern,"n.time"]
        if(is.null(iName.param)){return(NULL)}
browser()
        iTriplet <- triplet.vcov[[Upattern[iPattern,"name"]]]
        iTriplet.type <- attr(iTriplet,"typetype")
        iDhess <- stats::setNames(lapply(1:NCOL(iTriplet), function(iP){matrix(0, nrow = iN.time, ncol = iN.time)}), colnames(iTriplet))
        
        iOmega <- Omega[[iPattern]]

        if("sigmak.sigmak.sigma.k" %in% names(iTriplet.type)){
        }

        if("sigmak.sigmak.rho" %in% names(iPair.type)){
        }

        if("sigmak.rho.rho" %in% names(iPair.type)){
        }

        if("rho.rho.rho" %in% names(iPair.type)){
        }

        return(stats::setNames(iHess,colnames(iPair)))        
    })

    ## ** export
    out <- stats::setNames(out,Upattern$name)
    attr(out, "pair") <- pair.vcov
    return(out)
} 

## * calc_d3Omega.IND
.calc_d3Omega.IND <- .calc_d3Omega.ID

## * calc_d3Omega.CS
.calc_d3Omega.CS <- .calc_d3Omega.ID

## * calc_d3Omega.RE
.calc_d3Omega.RE <- .calc_d3Omega.ID

## * calc_d3Omega.TOEPLITZ
.calc_d3Omega.TOEPLITZ <- .calc_d3Omega.ID

## * calc_d3Omega.UN
.calc_d3Omega.UN <- .calc_d3Omega.ID


##----------------------------------------------------------------------
### structure-calc_d3Omega.R ends here
