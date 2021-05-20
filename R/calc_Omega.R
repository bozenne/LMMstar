### calc_Omega.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 21 2021 (18:12) 
## Version: 
## Last-Updated: May 20 2021 (11:31) 
##           By: Brice Ozenne
##     Update #: 179
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
.calc_dOmega <-  function(object, sigma, k, rho, Omega, Jacobian){

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

        ## apply transformation
        if(!is.null(Jacobian)){
            ## [dOmega_[11]/d theta_1] ... [dOmega_[11]/d theta_p] %*% Jacobian
            ## [dOmega_[ij]/d theta_1] ... [dOmega_[ij]/d theta_p] %*% Jacobian
            ## [dOmega_[mm]/d theta_1] ... [dOmega_[mm]/d theta_p] %*% Jacobian
            M.iScore <- do.call(cbind,lapply(iScore,as.double)) %*% Jacobian
            iScore <- setNames(lapply(1:NCOL(M.iScore), function(iCol){matrix(M.iScore[,iCol], nrow = iNtime, ncol = iNtime, byrow = FALSE)}),allCoef)
        }
        return(iScore)
    })

    ## ** export
    out <- setNames(out,Upattern)
    return(out)
}


## * .calc_d2Omega
.calc_d2Omega <- function(object, sigma, k, rho, Omega, dOmega, pair, Jacobian, dJacobian){

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

    ## id of the coefficients associated to each column
    ## First coef: columns containing 1:p, second coef: columns containing p + (1:p), ...
    Mindex.pair <- do.call(cbind,lapply(1:NCOL(pair), function(iCol){
        iMatch <- match(pair[,iCol],allCoef)
        return(c(p*(iMatch[1]-1)+iMatch[2],p*(iMatch[2]-1)+iMatch[1]))
    }))
        
    Upattern <- attr(object,"Upattern")
    Upattern.param <- attr(object,"Upattern.param")
    nUpattern <- attr(object,"nUpattern")
    Upattern.time <- attr(object,"Upattern.time")
    UX.var <- attr(object,"UX.var")
    indicator <- attr(object,"indicator")

    if(!is.null(Jacobian)){
        JacobianM1 <- solve(Jacobian)
    }
    
    ## ** loop over covariance patterns
    out <- lapply(1:nUpattern, function(iP){

        iTime <- Upattern.time[[iP]]
        iNtime <- length(iTime)
        iOmega.sd <- attr(Omega[[iP]],"sd")
        iOmega.cor <- attr(Omega[[iP]],"cor")
        iOmega <- Omega[[iP]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL; attr(iOmega,"time") <- NULL;
        idOmega <- dOmega[[iP]]

        iHess <- lapply(1:n.pair, function(iPair){
            matrix(0, nrow = iNtime, ncol = iNtime, dimnames = list(iTime, iTime))
        })

        for(iPair in attr(pair,"index.hessian")){ ## iPair <- 1
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
                        iHess[[iPair]] <- 2 * ind1 * iOmega / sigma[iCoef1]^2
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
                    iHess[[iPair]] <- diag(4*idOmega.k*iOmega.sd, nrow = iNtime, ncol = iNtime) + 2 * iOmega.cor * (idOmega.k %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k))
                    ## 4 * sigma * k["k.Y2"]
                    ## 2 * sigma * c(rho["cor(Y2,Y1)"],k["k.Y3"]*rho["cor(Y3,Y2)"],k["k.Y4"]*rho["cor(Y4,Y2)"])
                }else if(iType2 == "rho"){
                    iHess[[iPair]] <- 2 * ind1 * ind2 * iOmega / (sigma[iCoef1] * rho[iCoef2])
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

                }else if(iType2 == "rho"){
                    ## position where the parameters appear
                    ind.dk <- UX.var[[iP]][,iCoef1]
                    ## compute derivative
                    idk <- k; idk[iCoef1] <- 1
                    ## propagate
                    idOmega.k <- exp((UX.var[[iP]] %*% log(c(sigma,idk)))[,1]) * ind.dk
                    iHess[[iPair]] <- ind2 * (idOmega.k %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k))

                }
            }
        }

        ## apply transformation
        if(!is.null(Jacobian) || !is.null(dJacobian)){
            M.iScore <- do.call(cbind,lapply(dOmega[[iP]],as.double)) %*% JacobianM1
            iScore <- setNames(lapply(1:NCOL(M.iScore), function(iCol){matrix(M.iScore[,iCol], nrow = iNtime, ncol = iNtime, byrow = FALSE)}),allCoef)

            iHess2 <- vector(mode = "list", length = n.pair)
                
            for(iC in 1:p){ ## iC <- 2
                ##  [dOmega_[11]/d2 theta_1] ... [dOmega_[11]/d theta_1 d theta_p] %*% Jacobian + [dOmega_[11]/d theta_1] ... [dOmega_[11]/d theta_p] %*% dJacobian/d theta_1
                ##  [dOmega_[ij]/d2 theta_1] ... [dOmega_[ij]/d theta_1 d theta_p] %*% Jacobian + [dOmega_[ij]/d theta_1] ... [dOmega_[ij]/d theta_p] %*% dJacobian/d theta_1
                ##  [dOmega_[mm]/d2 theta_1] ... [dOmega_[mm]/d theta_1 d theta_p] %*% Jacobian + [dOmega_[mm]/d theta_1] ... [dOmega_[mm]/d theta_p] %*% dJacobian/d theta_1
                iCoef <- allCoef[iC]
                iIndex.pairCoef <- which(colSums(apply(Mindex.pair, 2, `%in%`, p*(iC-1)+(1:p)))>0)
                iCoef.pairCoef <- apply(pair[,iIndex.pairCoef,drop=FALSE], 2, function(iCol){
                    iOut <- setdiff(iCol,iCoef)
                    if(length(iOut)==0){iCoef}else{iOut}
                })
                iIndex.pairCoef <- iIndex.pairCoef[match(iCoef.pairCoef, allCoef)]
                M.iHess2 <- (do.call(cbind,lapply(iHess[iIndex.pairCoef],as.double)) %*% Jacobian + M.iScore %*% dJacobian[,,iCoef]) * Jacobian[iC,iC]
                iHess2[iIndex.pairCoef] <- lapply(1:NCOL(M.iScore), function(iCol){matrix(M.iHess2[,iCol], nrow = iNtime, ncol = iNtime, byrow = FALSE)})
            }
            return(iHess2)
        }else{
            return(iHess)
        }
    })

    ## ** export
    out <- setNames(out,Upattern)
    ## test <- unlist(lapply(out[[1]],function(i){all(abs(i)<1e-5)}))==FALSE
    ## out.print <- out[[1]][test]
    ## print(setNames(out.print,interaction(pair[1,test],pair[2,test])))
    return(out)
} 




##----------------------------------------------------------------------
### calc_Omega.R ends here
