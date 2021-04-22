### calc_Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 21 2021 (18:12) 
## Version: 
## Last-Updated: Apr 22 2021 (10:29) 
##           By: Brice Ozenne
##     Update #: 67
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .calc_Omega
.calc_Omega <- function(object, sigma, k, rho, keep.interim = FALSE){

    ## ** loop over covariance patterns
    out <- lapply(1:attr(object,"nUpattern"), function(iP){
        ## timepoint in this pattern
        iTime <- attr(object,"Upattern.time")[[iP]]
        iNtime <- length(iTime)
        ## diagonal (standard error)
        Omega.sd <- unname(exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,k)))[,1]))
        ## extra diagonal (correlation)
        if(length(rho)>0){
            Omega.cor <- matrix(attr(object,"UX.cor")[[iP]] %*% rho, nrow = iNtime, ncol = iNtime)
        }else{
            Omega.cor <- diag(0, nrow = iNtime, ncol = iNtime)
        }
        ## assemble
        Omega <- diag(Omega.sd^2, nrow = iNtime, ncol = iNtime) + Omega.cor * tcrossprod(Omega.sd)
        
        if(keep.interim){
            attr(Omega,"time") <- iTime
            attr(Omega,"sd") <- Omega.sd
            attr(Omega,"cor") <- Omega.cor
        }
        return(Omega)
    })

    ## ** export
    return(setNames(out,attr(object,"Upattern")))
}

## * .calc_dOmega
.calc_dOmega <-  function(object, sigma, k, rho, Omega, transform){

    ## ** prepare
    name.sigma <- names(sigma)
    n.sigma <- length(sigma)
    name.k <- names(k)
    n.k <- length(k)
    name.rho <- names(rho)
    n.rho <- length(rho)
    allCoef <- c(name.sigma, name.k, name.rho)
    p <- length(allCoef)

    Upattern <- attr(object,"Upattern")
    Upattern.param <- attr(object,"Upattern.param")
    nUpattern <- attr(object,"nUpattern")
    Upattern.time <- attr(object,"Upattern.time")
    UX.var <- attr(object,"UX.var")
    indicator <- attr(object,"indicator")
    
    ## ** loop over covariance patterns
    out <- lapply(1:nUpattern, function(iP){ ## iP <- 1
        iScore <- setNames(vector(mode = "list", length = p), allCoef)
        iTime <- Upattern.time[[iP]]
        iNtime <- length(iTime)
        iOmega.sd <- attr(Omega[[iP]],"sd")
        iOmega.cor <- attr(Omega[[iP]],"cor")
        iOmega <- Omega[[iP]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL; attr(iOmega,"time") <- NULL;
     
        if(length(sigma)>0){
            for(iSigma in 1:n.sigma){ ## iSigma <- 1
                ## positions where the k-parameter appears
                ind.dsigma <- UX.var[[iP]][,name.sigma[iSigma]]
                ## compute derivative
                idsigma <- sigma
                idsigma[name.sigma[iSigma]] <- 1
                ## propagate
                idOmega.sigma <- exp((UX.var[[iP]] %*% log(c(idsigma,k)))[,1]) * ind.dsigma
                iScore[[name.sigma[iSigma]]] <- diag(2*idOmega.sigma*iOmega.sd, nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.sigma %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.sigma))
                ## iScore[[name.sigma[iSigma]]] - 2*iOmega/sigma
            }
        }
        if(length(k)>0){
            for(iK in 1:n.k){ ## iK <- 1
                ## position where the k-parameter appears
                ind.dk <- UX.var[[iP]][,name.k[iK]]
                ## compute derivative
                idk <- k
                idk[name.k[iK]] <- 1
                ## propagate
                idOmega.k <- exp((UX.var[[iP]] %*% log(c(sigma,idk)))[,1]) * ind.dk
                iScore[[name.k[iK]]] <- diag(2*idOmega.k*iOmega.sd, nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.k %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k))
                ## iScore[[name.k[iK]]] - (tcrossprod(ind.dk,rep(1,4))+tcrossprod(rep(1,4),ind.dk)) * iOmega/k[name.k[iK]]
            }
        }
        if(length(rho)>0){
            for(iRho in 1:n.rho){ ## iRho <- 1
                ## positions where the rho-parameter appears
                ind.drho <- indicator[[iP]][[name.rho[iRho]]]
                ## derivative
                iScore[[name.rho[iRho]]] <- ind.drho * tcrossprod(iOmega.sd)
                ## iOmega/iOmega.corn
            }
        }

        if(transform>0){
            M.transform <- .transformDeriv(transform, sigma = sigma, k = k, rho = rho, pattern.param = Upattern.param[[iP]])
            M.iScore <- do.call(cbind,lapply(iScore,as.double)) %*% M.transform
            iScore <- setNames(lapply(1:NCOL(M.iScore), function(iCol){matrix(M.iScore[,iCol], nrow = iNtime, ncol = iNtime, byrow = FALSE)}),allCoef)
        }
        return(iScore)
    })
    return(setNames(out,Upattern))
}


## * .calc_d2Omega
.calc_d2Omega <- function(object, sigma, k, rho, Omega, dOmega, pair, transform){

    ## ** prepare
    name.sigma <- names(sigma)
    n.sigma <- length(sigma)
    name.k <- names(k)
    n.k <- length(k)
    name.rho <- names(rho)
    n.rho <- length(rho)
    allCoef <- c(name.sigma, name.k, name.rho)
    p <- length(allCoef)
    n.pair  <- NCOL(pair)
        
    Upattern <- attr(object,"Upattern")
    Upattern.param <- attr(object,"Upattern.param")
    nUpattern <- attr(object,"nUpattern")
    Upattern.time <- attr(object,"Upattern.time")
    UX.var <- attr(object,"UX.var")
    indicator <- attr(object,"indicator")

    ## ** loop over covariance patterns
    out <- lapply(1:nUpattern, function(iP){
        iHess <- vector(mode = "list", length = n.pair)

        iTime <- Upattern.time[[iP]]
        iNtime <- length(iTime)
        iOmega.sd <- attr(Omega[[iP]],"sd")
        iOmega.cor <- attr(Omega[[iP]],"cor")
        iOmega <- Omega[[iP]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL; attr(iOmega,"time") <- NULL;
        idOmega <- dOmega[[iP]]

        for(iPair in 1:n.pair){ ## iPair <- 1
            ## name of parameters
            iCoef1 <- pair[1,iPair]
            iCoef2 <- pair[2,iPair]

            ## type of parameters
            iType1 <- c("sigma","k","rho")[which(c(iCoef1 %in% name.sigma,iCoef1 %in% name.k,iCoef1 %in% name.rho))]
            iType2 <- c("sigma","k","rho")[which(c(iCoef2 %in% name.sigma,iCoef2 %in% name.k,iCoef2 %in% name.rho))]

            ## positions where the parameters appears in the matrix
            ind1 <- indicator[[iP]][[iCoef1]]
            ind2 <- indicator[[iP]][[iCoef2]]
                    
            if(iType1 == "sigma"){
                if(iType2 == "sigma"){
                    if(iCoef1==iCoef2){
                        ## transformation
                        ## d2 sigma^2         -> 2 --> sigma^2 too much while 2 missing                                
                        ## d2 exp(2 logsigma) -> 4 exp(2 logsigma) = 4 sigma^2 --> 4 missing
                        ## d2 sigma2          -> 0
                        if(transform == 0){
                            iHess[[iPair]] <- 2 * ind1 * iOmega / sigma[iCoef1]^2
                        }else if(transform == 1){ 
                            iHess[[iPair]] <- 4 * ind1 * iOmega
                        }else if(transform == 2){ 
                            iHess[[iPair]] <- 0 * iOmega
                        }
                    }else if(all(abs(ind1 * ind2)<1e-6)){
                        iHess[[iPair]] <- 0*iOmega
                    }else{
                        stop("Cannot compute the Hessian with interacting sigma coefficients. \n")
                    }
                }else if(iType2 == "k"){
                    ## position where the parameters appear
                    ind.dksigma <- UX.var[[iP]][,iCoef1] * UX.var[[iP]][,iCoef2]
                    ## compute derivative
                    idsigma <- sigma; idsigma[iCoef1] <- 1
                    idk <- k; idk[iCoef2] <- 1
                    ## propagate
                    idOmega.k <- exp((UX.var[[iP]] %*% log(c(idsigma,idk)))[,1]) * ind.dksigma
                    iHess[[iPair]] <- diag(4*idOmega.k*iOmega.sd, nrow = iNtime, ncol = iNtime) + iOmega.cor * (iOmega.sd %*% t(idOmega.k) + iOmega.sd %*% t(idOmega.k))
                        
                        ## transformation
                        ## d sigma^2 d k^1or2 or d sigma^2 k   -> 4or2 sigma k^1or0         --> sigma k too much while 4/2 missing                                
                        ## d exp(2 logsigma) d exp(1or2 logk) -> 4or2 sigma^2 k^2or1     --> 4/2 missing
                        ## d (sigma2 k2) or sqrt(k2)         -> 1 or 1/(2 k)               --> sigma^2 k^2 too much while nothing or 1/(2k) missing
                        if(transform == 1){
                            iHess[[iPair]] <- iHess[[iPair]] * sigma[iCoef1] * k[iCoef2]
                        }else if(transform == 2){
                            browser()
                            iHess[[iPair]] <- iHess[[iPair]] / (4 * sigma[iCoef1] * k[iCoef2])
                        }
                        
                    }else if(iType2 == "rho"){
                        
                        ## transformation
                        ## d sigma^2 d rho                -> 2 sigma                                              --> sigma rho too much while 2 missing                                
                        ## d exp(2 logsigma) d atanh(rho) -> 2 exp(2 logsigma) (1-rho^2) = 2 sigma^2 * (1-rho^2) --> rho too much while 2*(1-rho^2) missing                                
                        ## d sigma2 d rho                -> 1                                                     --> sigma^2 too much
                        iHess[[iPair]] <- 2 * ind1 * ind2 * iOmega / (sigma[iCoef1] * rho[iCoef2])
                        if(transform == 1){ ## exp(2*logsigma) * \Omega
                            iHess[[iPair]] <- iHess[[iPair]] * sigma[iCoef1] * (1-rho[iCoef2]^2)
                        }else if(transform == 2){ ## sigma2 * \Omega
                            iHess[[iPair]] <- iHess[[iPair]] / (2 * sigma[iCoef1])
                        }
                    }
            }else if(iType1 == "k"){                    
                if(iType2 == "k"){
                    ## position where the parameter appear
                    ind.dk1 <- UX.var[[iP]][,iCoef1]
                    ind.dk2 <- UX.var[[iP]][,iCoef2]
                    ## compute derivative
                    idk1 <- k ; idk1[iCoef1] <- 1
                    idk2 <- k ; idk2[iCoef2] <- 1
                    ## propagate
                    idOmega.k1 <- exp((UX.var[[iP]] %*% log(c(sigma,idk1)))[,1]) * ind.dk1
                    idOmega.k2 <- exp((UX.var[[iP]] %*% log(c(sigma,idk2)))[,1]) * ind.dk2
                    iHess[[iPair]] <- diag(2*idOmega.k1*idOmega.k2, nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.k1 %*% t(idOmega.k2) + idOmega.k2 %*% t(idOmega.k1))
                        
                    ## transformation
                    ## d2 k^2  or d k d k'       = 2 k or 1
                    ## d2 exp(2 k)  or d exp(k) d exp(k') =  2 k^2 or k k'
                    ## d2 d (sigma2 k2) or d (sigma2 k2) d (sigma2' k2') = ????
                    if(transform == 1){ 
                        iHess[[iPair]] <- iHess[[iPair]] * (1+(iCoef1==iCoef2)) * k[iCoef1] * k[iCoef2]
                    }else if(transform == 2){
                        browser()
                    }
                    
                    }else if(iType2 == "rho"){
                        ## position where the parameters appear
                        ind.dk <- UX.var[[iP]][,iCoef1]
                        ## compute derivative
                        idk <- k; idk[iCoef1] <- 1
                        ## propagate
                        idOmega.k <- exp((UX.var[[iP]] %*% log(c(sigma,idk)))[,1]) * ind.dk
                        iHess[[iPair]] <- ind.2 * (idOmega.k %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k))

                        ## transformation
                        ## d k^2or1 d rho                -> 2or1 k^1or0                      
                        ## d exp(2or1 logk) d atanh(rho) -> 2or1 exp(2or1 logk) (1-rho^2) = 2or1 k^2or1 * (1-rho^2) --> missing k and 1-rho^2
                        ## d sigma2 k or sqrt(k)  d rho  -> ?????
                        if(transform == 1){ 
                            iHess[[iPair]] <- iHess[[iPair]] * k[iCoef1] * (1-rho[iCoef2]^2)
                        }else if(transform == 2){
                            browser()
                        }
                    }
            }else if(iType1 == "rho"){
                iHess[[iPair]] <- matrix(0, nrow = iNtime, ncol = iNtime, dimnames = list(U.time, U.time))
            }
        }
        return(iHess)
    })
    return(setNames(out,Upattern))
} 


## * .transformDeriv
.transformDeriv <- function(transform, sigma, k, rho, pattern.param){
    n.sigma <- length(sigma)
    n.k <- length(k)
    n.rho <- length(rho)

    sigma.pattern <- names(sigma)[names(sigma) %in% pattern.param]
    k.pattern <- names(k)[names(k) %in% pattern.param]
    rho.pattern <- names(rho)[names(rho) %in% pattern.param]
    n.sigma.pattern <- length(sigma.pattern)
    n.k.pattern <- length(k.pattern)
    n.rho.pattern <- length(rho.pattern)
    
    M.transform <- matrix(0, nrow = n.sigma+n.k+n.rho, ncol = n.sigma+n.k+n.rho,
                          dimnames = list(c(names(sigma),names(k),names(rho)),c(names(sigma),names(k),names(rho)))
                          )
    if(transform==1){
        ## \sigma = \exp(logsigma) so d\exp(logsigma) = \exp(logsigma) = \sigma
        M.transform[sigma.pattern,sigma.pattern] <- sigma[sigma.pattern]
        if(n.k.pattern>0){
            ## k = \exp(logk) so d\exp(logk) = \exp(logk) = k
            for(iK in k.pattern){
                M.transform[iK,iK] <- k[iK]
            }
        }
        if(n.rho.pattern>0){
            ## \rho = \tanh(atanhrho) so d\tanh(atanhrho) = (1-\tanh(atanhrho)^2) =  1-\rho^2
            for(iRho in rho.pattern){
                M.transform[iRho,iRho] <- 1-rho[iRho]^2
            }
        }
    }else if(transform==2){
        ## sqrt(\sigma1), sqrt(\sigma2/\sigma1), ...   so  1/(2 sqrt(\sigma1)) \sqrt(\sigma2)/(2 \sigma1^(3/2))       
        ##                                                      0              1 / (2*\sqrt(\sigma2)\sqrt(\sigma1))
        ##                            

        M.transform[sigma.pattern,sigma.pattern] <- 1/(2*sigma[sigma.pattern])
        if(n.k.pattern>0){
            M.transform[k.pattern,sigma.pattern] <- -k[k.pattern]/(2*sigma[sigma.pattern]^2)
            for(iK in k.pattern){
                M.transform[iK,iK] <- 1/(2*sigma[sigma.pattern]^2*k[iK])
            }
        }
    }else{
        stop("Unknown transformation. Should be 0, 1, or 2. \n")
    }
    return(M.transform)
}

## FCT_TRANS <- function(p){
##     c(sqrt(p[1]), sqrt(p[2]/p[1]))
## }
## x <- 2; y  <- 5; a <- sqrt(x); b <- sqrt(y/x)
## jacobian(FCT_TRANS, c(x,y))
## 1/(2*a); -b/(2*a^2); 1/(2*a^2*b)
## 1/(2*sqrt(x)); -sqrt(y)/(2*x^(3/2)); 1/(2*sqrt(x)*sqrt(y))


##----------------------------------------------------------------------
### calc_Omega.R ends here
