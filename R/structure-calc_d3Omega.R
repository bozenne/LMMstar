### calc_d3Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 16 2021 (13:18) 
## Version: 
## Last-Updated: okt  8 2026 (13:30) 
##           By: Brice Ozenne
##     Update #: 545
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
##' @param param [named numeric vector] values of the parameters (transformed).
##' @param Omega [list of matrices] residual Variance-Covariance Matrix for each pattern.
##' @param triplet [list of data.frame] first three columns contain triplets of variance-covaraince parameters.
##' Following columns indicate whether the triplet is present in each covariance pattern.
##' @param transform.sigma,transform.k,transform.rho [character] transformation used on the variance/correlation coefficients.
##' Only active if \code{"log"}, \code{"log"}, \code{"atanh"}: then the derivative is directly computed on the transformation scale instead of using the Jacobian.
##' @param Upattern [data.frame] optional, used to only evaluate the third derivative of the residual variance-covariance with respect to a subset of patterns.
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
    function(object, param, Omega, triplet, transform.sigma, transform.k, transform.rho,
             Upattern) UseMethod(".calc_d3Omega")

## * calc_d3Omega.ID
.calc_d3Omega.ID <- function(object, param, Omega, triplet, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL,
                             Upattern = NULL){

    ## ** prepare
    ## pattern
    Upattern <- object$Upattern
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern
    
    ## param
    type <- stats::setNames(object$param$type,object$param$name)
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
            iTriplet <- triplet[which(triplet[[Upattern[iPattern,"name"]]]),c("name","param1","param2","param3"),drop=FALSE]
            n.iTriplet <- NROW(iTriplet)
            iOut <- replicate(n = n.iTriplet, matrix(0, nrow = iNtime, ncol = iNtime), simplify = FALSE)
            names(iOut) <- iTriplet$name            
        }

        ## *** loop over all pairs of parameters
        for(iiTrip in 1:n.iTriplet){ ## iiTriplet <- 4

            ## name of parameters
            iCoef1 <- iTriplet[iiTrip,"param1"]
            iCoef2 <- iTriplet[iiTrip,"param2"]
            iCoef3 <- iTriplet[iiTrip,"param3"]

            ## type of parameters
            iType1 <- type[iCoef1]
            iType2 <- type[iCoef2]
            iType3 <- type[iCoef3]

            ## first derivatives 
            if(iType1 == "sigma"){
                iDparam1 <- .dsigma_transform(value = param[iCoef1], Omega.sd = iOmega.sd, transform = transform.sigma, power = 1)
            }else if(iType1 == "k"){
                iDparam1 <- .dk_transform(value = param[iCoef1], indicator = (X.var[[iPattern.var]][,"k"] == iCoef1), Omega.sd = iOmega.sd, transform = transform.k, power = 1)
            } ## if iType1 == "rho" only requires second derivative as iCoef1 must be equal to iCoef2 and iCoef3 otherwise derivative = 0

            if(iType2 == "sigma"){
                if(iCoef1 != iCoef2){
                    stop("Second Omega derivative cannot handle multiple sigma parameters in a single pattern. \n")
                }
                iDparam2 <- iDparam1
            }else if(iType2 == "k"){
                iDparam2 <- .dk_transform(value = param[iCoef2], indicator = (X.var[[iPattern.var]][,"k"] == iCoef2), Omega.sd = iOmega.sd, transform = transform.k, power = 1)
            } ## if iType2 == "rho" only requires second derivative as iCoef2 must be equal to iCoef1 and iCoef3 otherwise derivative = 0

            if(iType3 == "sigma"){
                if(iCoef1 != iCoef3){
                    stop("Second Omega derivative cannot handle multiple sigma parameters in a single pattern. \n")
                }
                iDparam3 <- iDparam1
            }else if(iType3 == "k"){
                iDparam3 <- .dk_transform(value = param[iCoef3], indicator = (X.var[[iPattern.var]][,"k"] == iCoef3), Omega.sd = iOmega.sd, transform = transform.k, power = 1)
            }else if(iType3 == "rho" & (transform.rho != "cov") & iType2 %in% c("sigma","k")){
                iDparam3 <- .drho_transform(value = param[iCoef3], indicator = (X.cor[[iPattern.cor]][,,"rho"] == iCoef3), transform = transform.rho, power = 1)
            }

            ## second derivative
            if(iType1 == "sigma" && iType2 == "sigma"){
                iDparam12 <- .dsigma_transform(value = param[iCoef1], Omega.sd = iOmega.sd, transform = transform.sigma, power = 2)
            }else if(iType1 == "sigma" && iType2 == "k"){
                iDparam12 <- .dsigmak_transform(value = param[c(iCoef1,iCoef2)], indicator = (X.var[[iPattern.var]][,"k"] == iCoef2), Omega.sd = iOmega.sd, transform = c(transform.sigma,transform.k), power = c(1,1))
            }else if(iType1 == "k" && iType2 == "k"){
                iDparam12 <- (iCoef1 == iCoef2) * .dk_transform(value = param[iCoef1], indicator = (X.var[[iPattern.var]][,"k"] == iCoef1), Omega.sd = iOmega.sd, transform = transform.sigma, power = 2)
            } ## no need for iDparam13 between sigma/k and rho terms
            
            if(iType1 == "sigma" && iType3 == "sigma"){
                iDparam13 <- .dsigma_transform(value = param[iCoef1], Omega.sd = iOmega.sd, transform = transform.sigma, power = 2)
            }else if(iType1 == "sigma" && iType3 == "k"){
                iDparam13 <- .dsigmak_transform(value = param[c(iCoef1,iCoef3)], indicator = (X.var[[iPattern.var]][,"k"] == iCoef3), Omega.sd = iOmega.sd, transform = c(transform.sigma,transform.k), power = c(1,1))
            }else if(iType1 == "k" && iType3 == "k"){
                iDparam13 <- (iCoef1 == iCoef3) * .dk_transform(value = param[iCoef1], indicator = (X.var[[iPattern.var]][,"k"] == iCoef1), Omega.sd = iOmega.sd, transform = transform.sigma, power = 2)
            } ## no need for iDparam13 between sigma/k and rho terms 

            if(iType2 == "sigma" && iType3 == "sigma"){
                iDparam23 <- .dsigma_transform(value = param[iCoef2], Omega.sd = iOmega.sd, transform = transform.sigma, power = 2)
            }else if(iType2 == "sigma" && iType3 == "k"){
                iDparam23 <- .dsigmak_transform(value = param[c(iCoef2,iCoef3)], indicator = (X.var[[iPattern.var]][,"k"] == iCoef3), Omega.sd = iOmega.sd, transform = c(transform.sigma,transform.k), power = c(1,1))
            }else if(iType2 == "k" && iType3 == "k"){
                iDparam23 <- (iCoef2 == iCoef3) * .dk_transform(value = param[iCoef2], indicator = (X.var[[iPattern.var]][,"k"] == iCoef2), Omega.sd = iOmega.sd, transform = transform.sigma, power = 2)
            }else if(iType2 == "rho" && (transform.rho != "cov") && (iType1 != "rho")){ ## for the first type to be rho it means that both types are rho
                iDparam23 <- (iCoef2 == iCoef3) * .drho_transform(value = param[iCoef2], indicator = (X.cor[[iPattern.cor]][,,"rho"] == iCoef2), transform = transform.rho, power = 2)
            } ## no need for iDparam13 between sigma/k and rho terms

            ## third derivative
            if(iType1 == "sigma" && iType2 == "sigma" && iType3 == "sigma"){
                iDparam123 <- .dsigma_transform(value = param[iCoef1], Omega.sd = iOmega.sd, transform = transform.sigma, power = 3)
            }else if(iType1 == "sigma" && iType2 == "sigma" && iType3 == "k"){
                iDparam123 <- .dsigmak_transform(value = param[c(iCoef1,iCoef3)], indicator = (X.var[[iPattern.var]][,"k"] == iCoef3),
                                                 Omega.sd = iOmega.sd, transform = c(transform.sigma,transform.k), power = c(2,1))
            }else if(iType1 == "sigma" && iType2 == "k" && iType3 == "k"){
                iDparam123 <- (iCoef2 == iCoef3) * .dsigmak_transform(value = param[c(iCoef1,iCoef3)], indicator = (X.var[[iPattern.var]][,"k"] == iCoef3),
                                                                      Omega.sd = iOmega.sd, transform = c(transform.sigma,transform.k), power = c(1,2))
            }else if(iType1 == "k" && iType2 == "k" && iType3 == "k"){
                iDparam123 <- (iCoef1 == iCoef2 && iCoef2 == iCoef3) * .dk_transform(value = param[iCoef1], indicator = (X.var[[iPattern.var]][,"k"] == iCoef1),
                                                                                     Omega.sd = iOmega.sd, transform = transform.sigma, power = 3)
            }else if(iType1 == "rho" && iType2 == "rho" && iType3 == "rho"){
                iDparam123 <- (iCoef1 == iCoef2 && iCoef2 == iCoef3) * .drho_transform(value = param[iCoef1], indicator = (X.cor[[iPattern.cor]][,,"rho"] == iCoef1),
                                                                                       transform = transform.rho, power = 3)
            }

            ## assemble
            if(iType1 %in% c("sigma","k") && iType2 %in% c("sigma","k") && iType3 %in% c("sigma","k")){

                term1 <- tcrossprod(iDparam123, iOmega.sd) + tcrossprod(iDparam12, iDparam3) + tcrossprod(iDparam13, iDparam2) + tcrossprod(iDparam1, iDparam23)
                term2 <- tcrossprod(iDparam23, iDparam1) + tcrossprod(iDparam2, iDparam13) + tcrossprod(iDparam3, iDparam12) + tcrossprod(iOmega.sd, iDparam123)
                iOut[[iiTrip]] <- (term1 + term2) * iOmega.cor

            }else if(iType3 %in% c("rho")){

                if(transform.rho == "cov"){
                    ## iOut[[iiTrip]] <- matrix(0, nrow = iNtime, ncol = iNtime) ## do nothing this is already the case.
                    ## first derivative with respect to cov was constant
                }else if(iType1 %in% c("sigma","k") && iType2 %in% c("sigma","k")){
                    iOut[[iiTrip]] <- (tcrossprod(iDparam12, iOmega.sd) + tcrossprod(iDparam1, iDparam2) + tcrossprod(iDparam2, iDparam1) + tcrossprod(iOmega.sd, iDparam12)) * iDparam3
                }else if(iType1 %in% c("sigma","k") && iType2 %in% c("rho")){
                    iOut[[iiTrip]] <- (tcrossprod(iOmega.sd, iDparam1) + tcrossprod(iDparam1, iOmega.sd)) * iDparam23
                }else if(iType1 %in% c("rho") && iType2 %in% c("rho") && transform.rho == "atanh"){
                    iOut[[iiTrip]] <- tcrossprod(iOmega.sd) * iDparam123
                }
            }
        }
        return(iOut)        
    })

    ## ** export
    return(stats::setNames(out,Upattern$name))
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


## * helper

##----------------------------------------------------------------------
### calc_d3Omega.R ends here
