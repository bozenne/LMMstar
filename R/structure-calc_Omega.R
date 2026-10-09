### calc_Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 21 2021 (18:12) 
## Version: 
## Last-Updated: okt  7 2026 (11:39) 
##           By: Brice Ozenne
##     Update #: 688
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
##' @param param [named numeric vector] values of the parameters (transformed).
##' @param transform.sigma,transform.k,transform.rho [character] Transformation used on the variance/correlation coefficients.
##' @param Upattern [data.frame] Optional, used to only evaluate the residual variance-covariance with respect to a subset of patterns.
##' Should contain the name of the pattern, the index of the variance pattern, the index of the correlation pattern.
##' @param simplify [logical] should the correlation matrix and the vector of standard deviations be add to the output as attributes.
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
##' .calc_Omega(Sun24, param = param24, simplify = FALSE)
`.calc_Omega` <-
    function(object, param, transform.sigma, transform.k, transform.rho,
             Upattern, simplify) UseMethod(".calc_Omega")


## * calc_Omega.ID
.calc_Omega.ID <- function(object, param, transform.sigma, transform.k, transform.rho,
                           Upattern = NULL, simplify = TRUE){
    
    ## ** prepare
    ## pattern
    if(is.null(Upattern)){
        Upattern <- object$Upattern
    }
    n.Upattern <- NROW(Upattern)
    X.var <- object$var$Xpattern
    X.cor <- object$cor$Xpattern

    ## param
    param <- param[object$param$name] ## re-order and possibly remove mu parameters
    type <- object$param$type
    param1 <- c("one" = 1,param[type != "mu"])

    ## ** back-transform
    if(transform.sigma != "none"){
        name.sigma <- names(param)[type == "sigma"]
        if(transform.sigma  == "log"){
            param1[name.sigma] <- exp(param[name.sigma])
        }else if(transform.sigma  == "square"){
            param1[name.sigma] <- sqrt(param[name.sigma])
        }else if(transform.sigma  == "logsquare"){
            param1[name.sigma] <- exp(0.5*param[name.sigma])
        }        
    }   
    if(transform.k != "none"){
        name.k <- names(param)[type == "k"]
        if(transform.k  %in% c("log","logsd")){
            param1[name.k] <- exp(param[name.k])
        }else if(transform.k %in% c("square","var")){
            param1[name.k] <- sqrt(param[name.k])
        }else if(transform.k  %in% c("logsquare","logvar")){
            param1[name.k] <- exp(0.5*param[name.k])
        }        
    }
    if(transform.rho == c("atanh")){
        name.rho <- names(param)[type == "rho"]
        param1[name.rho] <- tanh(param[name.rho])
    }

    ## ** loop over covariance patterns
    out <- lapply(1:n.Upattern, function(iPattern){ ## iPattern <- 1

        ## *** patterns
        iX.var <- X.var[[Upattern[iPattern,"var"]]]
        iX.cor <- X.cor[[Upattern[iPattern,"cor"]]]
        iNtime <- Upattern[iPattern,"n.time"]

        ## *** variance
        ## convert vector of sigma: sigma sigma sigma sigma
        ##                of k    : 1     k2    k3    k4
        ## into values based of param1
        ## then take the product to get
        ##                       : sigma k2*sigma k3*sigma k4*sigma
        if(transform.k %in% c("none","log","square","logsquare")){
            Omega.sd <- apply(matrix(param1[iX.var], iNtime, NCOL(iX.var)), 1, prod)
        }else if(transform.k %in% c("sd","logsd","var","logvar")){
            Omega.sd <- param1[ifelse(iX.var[,"k"]=="one","sigma",iX.var[,"k"])]
        }

        ## *** correlation
        if(is.null(iX.cor) || iNtime == 1){
            Omega.cor <- diag(1, nrow = iNtime, ncol = iNtime)
        }else{
            ## convert matrix of rho: 1     rho12 rho13 rho14
            ##                        rho12     1 rho23 rho24
            ##                        rho13 rho23     1 rho34
            ##                        rho14 rho24 rho34     1
            ## into values and multiply them
            Omega.cor <- matrix(param1[iX.cor[,,"rho"]], nrow = iNtime, ncol = iNtime)
        }

        ## *** assemble
        if(transform.rho %in% c("none","atanh")){
            Omega <- tcrossprod(Omega.sd)*Omega.cor
        }else{
            Omega <- Omega.cor
            diag(Omega) <- Omega.sd^2
        }
        if(simplify == FALSE){
            attr(Omega,"sd") <- Omega.sd
            attr(Omega,"cor") <- Omega.cor
        }
        return(Omega)
    })

    ## print(Omega)
    return(stats::setNames(out,Upattern$name))
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
.calc_Omega.EXP <- function(object, param, Upattern = NULL, simplify = TRUE){

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
        
        if(simplify == FALSE){
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
.calc_Omega.CUSTOM <- function(object, param, Upattern = NULL, simplify = TRUE){

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
        
        if(simplify == FALSE){
            attr(iOmega,"sd") <- iOmega.sd
            attr(iOmega,"cor") <- iOmega.cor
        }
        return(iOmega)
    }), Upattern$name)

    return(Omega)
}


##----------------------------------------------------------------------
### calc_Omega.R ends here
