### structure-calc_d3Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  4 2023 (17:08) 
## Version: 
## Last-Updated: aug  8 2023 (18:42) 
##           By: Brice Ozenne
##     Update #: 60
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

        iTriplet <- triplet.vcov[[Upattern[iPattern,"name"]]]
        iTriplet.type <- attr(iTriplet,"typetype")
        iDhess <- stats::setNames(lapply(1:NCOL(iTriplet), function(iP){matrix(0, nrow = iN.time, ncol = iN.time)}), colnames(iTriplet))
        iOmega <- Omega[[iPattern]]

        if(transform){

            if("sigmak.sigmak.sigmak" %in% names(iTriplet.type)){

                for(iiTriplet in iTriplet.type$sigmak.sigmak.sigmak){ ## iiTriplet <- iTriplet.type$sigmak.sigmak.sigmak[1]
                    iiPosition <- attr(iTriplet,"index.triplet")[[iiTriplet]]$position
                    iiValue1 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value1
                    iiValue2 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value2
                    iiValue3 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value3                
                    iDhess[[iiTriplet]][iiPosition] <- iiValue1 * iiValue2 * iiValue3 * iOmega[iiPosition]
                }
            
            }
                
            if("sigmak.sigmak.rho" %in% names(iTriplet.type)){

                for(iiTriplet in iTriplet.type$sigmak.sigmak.rho){ ## iiTriplet <- iTriplet.type$sigmak.sigmak.rho[1]
                    iiPosition <- attr(iTriplet,"index.triplet")[[iiTriplet]]$position
                    iiValue1 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value1
                    iiValue2 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value2
                    iiValue3 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value3
                    iiParam3 <- param[iTriplet[3,iiTriplet]]                
                    iDhess[[iiTriplet]][iiPosition] <- iiValue1 * iiValue2 * iiValue3 * iOmega[iiPosition] * (1-iiParam3^2) / iiParam3
                }

            }

            if("sigmak.rho.rho" %in% names(iTriplet.type)){

                for(iiTriplet in iTriplet.type$sigmak.rho.rho){ ## iiTriplet <- iTriplet.type$sigmak.rho.rho[1]
                    iiPosition <- attr(iTriplet,"index.triplet")[[iiTriplet]]$position
                    iiValue1 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value1
                    iiValue2 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value2
                    if(iTriplet[2,iiTriplet]==iTriplet[3,iiTriplet]){
                        iiParam2 <- param[iTriplet[2,iiTriplet]]                
                        iDhess[[iiTriplet]][iiPosition] <- iiValue1 * iiValue2 * (1-iiParam2^2) * iOmega[iiPosition] * ((iiValue2-1)/iiParam2^2 - (iiValue2+1))
                    }else if(iTriplet[2,iiTriplet]!=iTriplet[3,iiTriplet]){
                        iiValue3 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value3                    
                        iiParam23 <- param[iTriplet[2:3,iiTriplet]]                
                        iDhess[[iiTriplet]][iiPosition] <- iiValue1 * iiValue2 * iiValue3 * iOmega[iiPosition] * prod(1-iiParam23^2) / prod(iiParam23)
                    }
                }

            }

            if("rho.rho.rho" %in% names(iTriplet.type)){

                for(iiTriplet in iTriplet.type$rho.rho.rho){ ## iiTriplet <- iTriplet.type$rho.rho.rho[1]
                    iiPosition <- attr(iTriplet,"index.triplet")[[iiTriplet]]$position

                    if(iTriplet[1,iiTriplet]==iTriplet[2,iiTriplet] && iTriplet[1,iiTriplet]==iTriplet[3,iiTriplet]){
                        ## \rho = tanh(x), x = atanh(\rho), dx/d\rho = 1/(1-\rho^2) so d\rho/dx = 1-\rho^2
                        ## \dOmega/\dx = d\rho^a/dx = a \rho^{a-1} (1-\rho^2) = a (\rho^{a-1}-\rho^{a+1})
                        ## \d^2Omega/\dx^2 = a ((a-1)\rho^{a-2}-(a+1)\rho^{a})(1-\rho^2)=a(1-\rho^2)\rho^a((a-1)/\rho^2-(a+1))
                        ##                 = a (a-1) \rho^{a-2} - a (a-1) \rho^a - a (a+1) \rho^a + a (a+1) \rho^{a+2}
                        ##                 = a (a-1) \rho^{a-2} - 2 a^2 \rho^a + a (a+1) \rho^{a+2}
                        ## \d^3Omega/\dx^3 = (a (a-1) (a-2) \rho^{a-3} - 2 a^3  \rho^{a-1} + a (a+1) (a+2) \rho^{a+1} )(1-\rho^2)
                        ##                 = a \rho^{a} ((a-1) (a-2) \rho^{-3} - 2 a^2  \rho^{-1} + (a+1) (a+2) \rho )(1-\rho^2)
                        iiValue <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value1
                        iiParam <- param[iTriplet[1,iiTriplet]]
                        iDhess[[iiTriplet]][iiPosition] <- iiValue * iOmega[iiPosition] * ((iiValue-1)*(iiValue-2)/iiParam^3 - 2*iiValue^2/iiParam + (iiValue+1)*(iiValue+2)*iiParam)
                            
                    }else if(iTriplet[1,iiTriplet]==iTriplet[2,iiTriplet]){
                        iiValue1 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value1
                        iiValue3 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value3                    
                        iiParam1 <- param[iTriplet[1,iiTriplet]]                
                        iiParam3 <- param[iTriplet[3,iiTriplet]]                
                            
                        iDhess[[iiTriplet]][iiPosition] <- iiValue1 * iiValue3 * iOmega[iiPosition] * prod(1-c(iiParam1,iiParam3)^2) * ((iiValue1-1)/iiParam1^2 - (iiValue1+1)) / iiParam3
                    }else if(iTriplet[1,iiTriplet]==iTriplet[3,iiTriplet]){
                        iiValue1 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value1
                        iiValue2 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value2                    
                        iiParam1 <- param[iTriplet[1,iiTriplet]]                
                        iiParam2 <- param[iTriplet[2,iiTriplet]]                
                            
                        iDhess[[iiTriplet]][iiPosition] <- iiValue1 * iiValue2 * iOmega[iiPosition] * prod(1-c(iiParam1,iiParam2)^2) * ((iiValue1-1)/iiParam1^2 - (iiValue1+1)) / iiParam2
                    }else if(iTriplet[2,iiTriplet]==iTriplet[3,iiTriplet]){                            
                        iiValue1 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value1
                        iiValue2 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value2                   
                        iiParam1 <- param[iTriplet[1,iiTriplet]]                
                        iiParam2 <- param[iTriplet[2,iiTriplet]]                
                            
                        iDhess[[iiTriplet]][iiPosition] <- iiValue1 * iiValue2 * iOmega[iiPosition] * prod(1-c(iiParam1,iiParam2)^2) * ((iiValue2-1)/iiParam2^2 - (iiValue2+1)) / iiParam1
                    }else{ ## all different
                        iiValue1 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value1
                        iiValue2 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value3                    
                        iiValue3 <- attr(iTriplet,"index.triplet")[[iiTriplet]]$value3                    
                        iiParam <- param[iTriplet[,iiTriplet]]                
                    
                        iDhess[[iiTriplet]][iiPosition] <- iiValue1 * iiValue2 * iiValue3 * iOmega[iiPosition] * prod(1-iiParam^2) / prod(iiParam)
                    }
                }
                
            }

        }else{ ## no transformation

            for(iiTriplet in NCOL(iTriplet)){ ## iiTriplet <- 1
                iiPosition <- attr(iTriplet,"index.triplet")[[iiTriplet]]$position
                iiD2value <- attr(iTriplet,"index.triplet")[[iiTriplet]]$dvalue
                iiParam <- param[iTriplet[,iiTriplet]]    
                iDhess[[iiTriplet]][iiPosition] <- iiD2value * iOmega[iiPosition] / prod(iiParam)                
            }

        }

        return(stats::setNames(iDhess,colnames(iTriplet)))        
    })

    ## ** export
    out <- stats::setNames(out,Upattern$name)
    attr(out, "triplet") <- attr(triplet.vcov,"global")
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
