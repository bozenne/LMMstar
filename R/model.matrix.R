### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: mar 27 2021 (00:00) 
##           By: Brice Ozenne
##     Update #: 254
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * model.matrix.lmm (code)
##' @export
model.matrix.lmm <- function(object, data = NULL, type = "lmm-mean"){

    ## ** normalize user imput
    type <- match.arg(type, c("lmm","lmm-mean","lmm-variance","gls"))
    if(!is.null(data)){
        ff.allvars <- c(all.vars(x$formula$mean), all.vars(x$formula$var))
        if(any(ff.allvars %in% names(data) == FALSE)){
            stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
        }

        design <- .model.matrix.lmm(formula.mean = x$formula$mean.design,
                                    formula.var = x$formula$var.design,
                                    data = data,
                                    var.outcome = x$outcome$var,
                                    var.strata = x$strata$var, U.strata = x$strata$levels,
                                    var.time = x$time$var, U.time = x$time$levels,
                                    var.cluster = x$cluster$var,
                                    structure = x$structure
                                    )
    }
    
    if(type == "lmm"){
        if(is.null(data)){
            return(list(mean = object$design$X.mean,
                        variance = object$design$X.var))
        }else{
            
            return(list(mean = design$X.mean,
                        variance = design$X.var))
        }
    }else if(type == "lmm-mean"){
        if(is.null(data)){
            return(object$design$X.mean)
        }else{
            return(design$X.mean)
        }
    }else if(type == "lmm-variance"){
        if(is.null(data)){
            return(object$design$X.var)
        }else{
            return(design$X.var)
        }
    }else if(type=="gls"){
        if(object$strata$n==1){
            return(model.matrix(object$gls[[1]], data = data))
        }else{
            return(lapply(object$gls, model.matrix, data = data))
        }
    }
}

## * .model.matrix.lmm
.model.matrix.lmm <- function(formula.mean, formula.var, data, var.outcome,
                              var.strata, U.strata,
                              var.time, U.time,
                              var.cluster,
                              structure, transform){
    n.obs <- NROW(data)
    n.strata <- length(U.strata)
    n.time <- length(U.time)

    ## ** normalize data
    if("XXindexXX" %in% names(data) == FALSE){
        data$XXindexXX <- 1:NROW(data)
    }
    if(identical(var.strata, "XXstrata.indexXX") && var.strata %in% names(data) == FALSE){
        data$XXstrata.indexXX <- 1
    }else{
        data[[var.strata]] <- factor(data[[var.strata]], levels = U.strata)
    }
    data[[var.time]] <- factor(data[[var.time]], levels = U.time)
    
    ## ** mean
    if(n.strata==1){
        X.mean <- model.matrix(formula.mean, data)
    }else{
        ls.X.mean <- lapply(U.strata, function(iS){ ## iS <- U.strata[1]
            iX <- model.matrix(formula.mean, data[data[[var.strata]]==iS,])
            colnames(iX) <- paste0(colnames(iX),":",iS)
            attr(iX,"index") <- data[data[[var.strata]]==iS,"XXindexXX"]
            return(iX)
        })
        X.mean <- as.matrix(Matrix::bdiag(ls.X.mean))[order(unlist(lapply(ls.X.mean, attr, "index"))),]
        colnames(X.mean) <- unlist(lapply(ls.X.mean,colnames))
    }

    ## ** cluster
    U.cluster <- sort(unique(data[[var.cluster]]))
    n.cluster <- length(U.cluster)
    index.cluster <- match(data[[var.cluster]], U.cluster) ## ‘match’ returns a vector of the positions of (first) matches of its first argument in its second.
    
    ## ** variance
    ## *** parametrisation
    ## **** sigma
    if(n.strata==1){
        param.sigma <- "sigma"
    }else{
        param.sigma <- paste0("sigma",":",U.strata)
    }

    ## **** k
    if(structure == "CS"){
        if(n.strata==1){
            X.var <- model.matrix(~1, data)
            colnames(X.var) <- param.sigma            
        }else{
            X.var <- model.matrix(as.formula(paste0("~0+",var.strata)), data)
            colnames(X.var) <- param.sigma
        }
        param.k <- NULL
    }else if(structure == "UN"){
        if(n.strata==1){
            X.var <- model.matrix(formula.var, data)
            param.k <- colnames(X.var)[-1]
            for(iTime in 1:n.time){ ## iTime <- 2
                param.k <- gsub(pattern = paste0("^",var.time,U.time[iTime]), replacement = paste0(U.time[iTime]), x = param.k)
            }
            param.k <- paste("k",param.k,sep=".")
            colnames(X.var) <- c(param.sigma,param.k)
        }else{
            terms.var <- delete.response(terms(formula.var))
            formula2.var <- update(terms.var, paste0("~0+",var.strata,"+",var.strata,":.")) ## using ".:var.strata" does not work (it gives the same formula - does not invert . var.strata around the : sympbol)
            X.var <- model.matrix(formula2.var, data)
            index.tempo <- lapply(1:n.strata,function(iS){
                grep(paste0(paste0("^",var.strata,U.strata[iS],":")), colnames(X.var))
            })
            param.k <- colnames(X.var)[-(1:n.strata)]
            for(iStrata in 1:n.strata){ ## iStrata <- 1
                param.k <- sapply(strsplit(x = param.k, split = paste0("^",var.strata,U.strata[iStrata],":")), function(iVec){
                    if(length(iVec)==2){paste0(iVec[2],":",U.strata[iStrata])}else{iVec}
                })
            }
            for(iTime in 1:n.time){ ## iTime <- 2
                param.k <- gsub(pattern = paste0("^",var.time,U.time[iTime]), replacement = paste0(U.time[iTime]), x = param.k)
            }
            colnames(X.var) <- c(param.sigma, param.k)
            X.var[,sort(unlist(index.tempo))] <- X.var[,unlist(index.tempo)]
            param.k <- setdiff(colnames(X.var), param.sigma)
        }
    }
         
    attr(X.var,"level") <- as.numeric(droplevels(interaction(as.data.frame(X.var))))
    attr(X.var,"transform") <- transform ## should parameters be rescale (var: log transform, cor: atanh transform)

    ## **** rho
    if(n.time==1){
        param.rho <- NULL
    }else{
        Mtempo <- matrix(1:16,4,4)
        tMtempo <- t(Mtempo)
        if(structure == "CS"){
            param.rho <- "Rho"
            M.index <- cbind(index.lower = Mtempo[lower.tri(Mtempo)],
                             index.upper = tMtempo[lower.tri(Mtempo)])
        }else if(structure == "UN"){
            M.index <- cbind(which(lower.tri(diag(1, nrow = n.time, ncol = n.time)), arr.ind = TRUE),
                             index.lower = Mtempo[lower.tri(Mtempo)],
                             index.upper = tMtempo[lower.tri(Mtempo)])
            param.rho <- paste0("cor","(",M.index[,"row"],",",M.index[,"col"],")")
        }
        param.rho.save <- param.rho
        if(n.strata>1){
            param.rho <- unlist(lapply(U.strata, function(iS){paste(param.rho.save, iS, sep = ":")}))
        }
        
        M.indexAlltimes <- which(matrix(1,n.time, n.time)==1, arr.ind = TRUE)
        X.cor <- vector(mode = "list", length = n.strata)
        for(iStrata in 1:n.strata){ ## iStrata <- 1
            X.cor[[iStrata]] <- matrix(0, nrow = n.time^2, ncol = length(param.rho),
                                       dimnames = list(paste0("(",M.indexAlltimes[,"row"],",",M.indexAlltimes[,"col"],")"),param.rho))
            if(structure == "CS"){
                X.cor[[iStrata]][as.double(M.index[,c("index.lower","index.upper")]),param.rho[iStrata]] <- 1
            }else if(structure == "UN"){
                for(iParam in 1:length(param.rho.save)){ ## iParam <- 2
                    if(n.strata==1){
                        X.cor[[iStrata]][M.index[iParam,c("index.lower","index.upper")],param.rho[iParam]] <- 1
                    }else{
                        X.cor[[iStrata]][M.index[iParam,c("index.lower","index.upper")],paste(param.rho.save[iParam],U.strata[iStrata],sep=":")] <- 1
                    }
                }
            }
        }
    }

    ## *** design matrix
    attr(X.var,"pattern") <- tapply(attr(X.var,"level"), index.cluster, paste, collapse=".")
    attr(X.var,"Upattern") <- unname(sort(unique(attr(X.var,"pattern"))))
    attr(X.var,"nUpattern") <- length(attr(X.var,"Upattern"))
    indexPattern.tempo <- which(duplicated(attr(X.var,"pattern"))==FALSE) ## first pattern corresponding to each unique pattern
    indexCluster.tempo <- names(sort(attr(X.var,"pattern")[indexPattern.tempo])) ## corresponding cluster
    indexObs.tempo <- lapply(indexCluster.tempo, function(iCluster){which(index.cluster==iCluster)})
    attr(X.var,"UX.strata") <- setNames(sapply(indexObs.tempo, function(iObs){unique(data[iObs,var.strata])}), attr(X.var,"Upattern")) ## strata associated to each unique pattern
    attr(X.var,"UX.var") <- setNames(lapply(indexObs.tempo, function(iObs){ ## extract design matrix 
        iM <- X.var[iObs,,drop=FALSE]
        rownames(iM) <- data[[var.time]][iObs]
        return(iM)
    }), attr(X.var,"Upattern"))
    attr(X.var,"Upattern.time") <- lapply(attr(X.var,"UX.var"), rownames)
    if(!is.null(param.rho)){
        attr(X.var,"UX.cor") <- lapply(1:attr(X.var,"nUpattern"), function(iPattern){
            iIndex.time <- which(U.time %in% attr(X.var,"Upattern.time")[[iPattern]] == FALSE)
            M.test <- M.indexAlltimes
            if(length(iIndex.time)>0){
                M.test[M.indexAlltimes %in% iIndex.time] <- NA
            }
            iStrata <- which(U.strata == attr(X.var,"UX.strata")[iPattern])
            return(X.cor[[iStrata]][rowSums(is.na(M.test))==0,,drop=FALSE])
        })
    }
    index.vargroup <- setNames(match(attr(X.var,"pattern"),attr(X.var,"Upattern")),names(attr(X.var,"pattern")))
    
    ## ** prepare calculation of the variance matrices
    attr(X.var, "FUN.Omega") <- function(object, sigma, k, rho, keep.interim = FALSE){
        out <- lapply(1:attr(object,"nUpattern"), function(iP){

            iTime <- attr(object,"Upattern.time")[[iP]]
            iNtime <- length(iTime)
            Omega.var <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,k)))[,1])
            if(!is.null(rho)){
                Omega.cor <- matrix(attr(object,"UX.cor")[[iP]] %*% rho, nrow = iNtime, ncol = iNtime, dimnames = list(iTime,iTime))
            }else{
                Omega.cor <- diag(0, nrow = iNtime, ncol = iNtime)
            }
            Omega <- diag(Omega.var^2, nrow = iNtime, ncol = iNtime) + Omega.cor * tcrossprod(Omega.var)
            dimnames(Omega) <- list(iTime,iTime)
            if(keep.interim){
                attr(Omega,"var") <- Omega.var
                attr(Omega,"cor") <- Omega.cor
            }
            return(Omega)
        })
        return(setNames(out,attr(object,"Upattern")))
    }

    ## ** prepare calculation of the score
    ## indicator of the effect of each parameter on each element of the covariance
    attr(X.var, "indicator") <- setNames(lapply(1:attr(X.var,"nUpattern"), function(iP){ ## iP <- 1

        iIndicator <- setNames(vector(mode = "list", length = length(c(param.sigma,param.k,param.rho))), c(param.sigma,param.k,param.rho))
        iTime <- attr(X.var,"Upattern.time")[[iP]]
        iNtime <- length(iTime)
         
        if(!is.null(param.sigma)){
            for(iSigma in 1:length(param.sigma)){ ## iSigma <- 1
                ## positions where the sigma-parameter appears in the matrix
                dsigma <- rep(0,length(param.sigma))
                dsigma[iSigma] <- 1
                ind.var_dsigma <- attr(X.var,"UX.var")[[iP]] %*% c(dsigma,rep(0, length(param.k)))
                iIndicator[[param.sigma[iSigma]]] <- tcrossprod(ind.var_dsigma)
            }
        }
        if(!is.null(param.k)){
            for(iK in 1:length(param.k)){ ## iK <- 1
                ## positions where the k-parameter appears in the matrix
                dk <- rep(0,length(param.k))
                dk[iK] <- 1
                ind.var_dk <- attr(X.var,"UX.var")[[iP]] %*% c(rep(0,length(sigma)),dk)
                iIndicator[[param.k[iK]]] <- ((ind.var_dk %*% rep(1,iNtime)) + t(ind.var_dk %*% rep(1,iNtime))) > 0
            }
        }
        if(!is.null(param.rho)){
            for(iRho in 1:length(param.rho)){ ## iRho <- 1
                ## positions where the rho-parameter appears in the matrix
                iIndicator[[param.rho[iRho]]] <- matrix(attr(X.var,"UX.cor")[[iP]][,iRho] , nrow = iNtime, ncol = iNtime, dimnames = list(iTime,iTime))
            }
        }
        return(iIndicator)
    }), attr(X.var,"Upattern"))

    attr(X.var, "FUN.dOmega") <- function(object, sigma, k, rho, Omega){
        name.sigma <- names(sigma)
        n.sigma <- length(sigma)
        name.k <- names(k)
        n.k <- length(k)
        name.rho <- names(rho)
        n.rho <- length(rho)
        allCoef <- c(name.sigma, name.k, name.rho)
        p <- length(allCoef)
        transform <- attr(object,"transform")

        out <- lapply(1:attr(object,"nUpattern"), function(iP){
            iScore <- setNames(vector(mode = "list", length = p), allCoef)

            iTime <- attr(object,"Upattern.time")[[iP]]
            iNtime <- length(iTime)
            iOmega.var <- attr(Omega[[iP]],"var")
            iOmega.cor <- attr(Omega[[iP]],"cor")
            iOmega <- Omega[[iP]] ; attr(iOmega,"var") <- NULL; attr(iOmega,"cor") <- NULL;

            if(!is.null(sigma)){
                for(iSigma in 1:n.sigma){ ## iSigma <- 1
                    ## compute derivative
                    idsigma <- sigma
                    idsigma[iSigma] <- 1
                    ## propagate
                    idOmega.sigma <- exp((attr(object,"UX.var")[[iP]] %*% log(c(idsigma,k)))[,1])
                    iScore[[name.sigma[iSigma]]] <- diag(2*idOmega.sigma*iOmega.var, nrow = iNtime, ncol = iNtime) + iOmega.cor * 2 * tcrossprod(idOmega.sigma, iOmega.var)
                    ## log transform
                    if(transform){ 
                        iScore[[name.sigma[iSigma]]] <- iScore[[name.sigma[iSigma]]]/sigma[iSigma]
                    }
                    ## check
                    ## iScore[[name.sigma[iSigma]]] - attr(object,"indicator")[[iP]][[name.sigma[iSigma]]] * iOmega * 2 / sigma[iSigma]
                }
            }
            if(!is.null(k)){
                for(iK in 1:n.k){ ## iSigma <- 1
                    ## positions where the k-parameter appears in the matrix
                    ind.dk <- attr(object,"indicator")[[iP]][[name.k[iK]]]
                    ## compute derivative
                    idk <- k
                    idk[iK] <- 1
                    ## propagate
                    idOmega.k <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,idk)))[,1])
                    tcrossprod(iOmega.var)
                    ind.dk *tcrossprod(idOmega.k,iOmega.var) + ind.dk *tcrossprod(iOmega.var, idOmega.k)
                    
iScore[[name.sigma[iSigma]]] <- diag(2*idOmega.k*iOmega.var, nrow = iNtime, ncol = iNtime) + iOmega.cor * (tcrossprod(idOmega.k, ind.dk * iOmega.var) + tcrossprod(iOmega.var, idOmega.k))
                    browser()
                    ## dk^2 = 2k ; dk = k is obtained via (2k)/(2k)=1 and (2k)^2/(2k) = 2k
                    k2 <- k
                    k2[iK] <- 2*k[iK]
                    iOmega.var2 <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,k2)))[,1])
                    iScore[[name.k[iK]]] <- ind.dk * (diag(iOmega.var2^2, nrow = iNtime, ncol = iNtime) + iOmega.cor * tcrossprod(iOmega.var2)) / (2*k[iK])
                    if(transform){ ## log transform
                        iScore[[name.k[iK]]] <- iScore[[name.k[iK]]]/k[iK]
                    }
                }
            }
            if(!is.null(rho)){
                for(iRho in 1:n.rho){ ## iRho <- 1
                    ## positions where the rho-parameter appears in the matrix
                    ind.drho <- attr(object,"indicator")[[iP]][[name.rho[iRho]]]
                    ## derivative
                    iScore[[name.rho[iRho]]] <- ind.drho * tcrossprod(iOmega.var)
                    if(transform){ ## atanh transform
                        iScore[[name.rho[iRho]]] <- iScore[[name.rho[iRho]]]/(1-rho[iRho]^2)
                    }
                }
            }
            return(iScore)
        })
        return(setNames(out,attr(object,"Upattern")))
    }

    ## ** prepare calculation of the df
    attr(X.var, "FUN.d2Omega") <- function(object, sigma, k, rho, Omega){
        browser()
        name.sigma <- names(sigma)
        n.sigma <- length(sigma)
        name.k <- names(k)
        n.k <- length(k)
        name.rho <- names(rho)
        n.rho <- length(rho)
        allCoef <- c(name.sigma, name.k, name.rho)
        p <- length(allCoef)
        transform <- attr(object,"transform")
        
        out <- lapply(1:attr(object,"nUpattern"), function(iP){
            pair.varcoef <- .unorderedPairs(c(name.sigma,name.k, name.rho))
            npair.varcoef <- NCOL(pair.varcoef)
            iHess <- vector(mode = "list", length = npair.varcoef)

            iTime <- attr(object,"Upattern.time")[[iP]]
            iNtime <- length(iTime)
            iOmega.var <- attr(Omega[[iP]],"var")
            iOmega.cor <- attr(Omega[[iP]],"cor")
            iOmega <- Omega[[iP]] ; attr(iOmega,"var") <- NULL; attr(iOmega,"cor") <- NULL;

            for(iPair in 1:n.pair){ ## iPair <- 1
                ## name of parameters
                iCoef1 <- pair.varcoef[1,iPair]
                iCoef2 <- pair.varcoef[2,iPair]

                ## type of parameters
                iType1 <- c("sigma","k","rho")[which(c(iCoef1 %in% name.sigma,iCoef1 %in% name.k,iCoef1 %in% name.rho))]
                iType2 <- c("sigma","k","rho")[which(c(iCoef2 %in% name.sigma,iCoef2 %in% name.k,iCoef2 %in% name.rho))]

                ## positions where the parameters appears in the matrix
                ind1 <- attr(object,"indicator")[[iP]][[iCoef1]]
                ind2 <- attr(object,"indicator")[[iP]][[iCoef2]]
                    
                if(iType1 == "sigma"){
                    if(iType2 == "sigma"){
                        iHess[[iPair]] <- ind1 * ind2 * iOmega * 2 / (sigma[iCoef1] * sigma[iCoef2])
                    }else if(iType2 == "k"){
                        ## dk^2 = 2k ; dk = k is obtained via (2k)/(2k)=1 and (2k)^2/(2k) = 2k
                        k2 <- k
                        k2[iCoef2] <- 2*k[iCoef2]
                        iOmega.var2 <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,k2)))[,1])
                        iHess[[iPair]] <- ind1 * ind2 * (diag(iOmega.var2^2, nrow = iNtime, ncol = iNtime) + iOmega.cor * tcrossprod(iOmega.var2)) * 2 / (k[iCoef2] * sigma[iCoef1])
                    }else if(iType2 == "rho"){
                        iHess[[iPair]] <- ind1 * ind2 * tcrossprod(iOmega.var) * 2 / sigma[iCoef1]
                    }
                }else if(iType1 == "k"){
                    if(iType2 == "k"){
                        ## dk^2 = 2k ; dk = k is obtained via (2k)/(2k)=1 and (2k)^2/(2k) = 2k
                        iHess[[iPair]] <- ind1 * ind2 * iOmega / (k[iCoef1] * k[iCoef2])
                    }else if(iType2 == "rho"){
                        ## dk^2 = 2k ; dk = k is obtained via (2k)/(2k)=1 and (2k)^2/(2k) = 2k
                        k2 <- k
                        k2[iCoef2] <- 2*k[iCoef2]
                        iOmega.var2 <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,k2)))[,1])
                        iHess[[iPair]] <- ind1 * ind2 * tcrossprod(iOmega.var2) / k[iCoef1]
                    }
                }else if(iType1 == "rho"){
                    iHess[[iPair]] <- matrix(0, nrow = iNtime, ncol = iNtime, dimnames = list(U.time, U.time))
                }
            }
            return(iHess)
        })
        return(setNames(out,attr(object,"Upattern")))
    }

    ## ** export
    return(list(X.mean = X.mean,
                X.var = X.var,
                Y = data[[var.outcome]],
                index.cluster = index.cluster,
                index.vargroup = index.vargroup,
                index.strata = tapply(data[[var.strata]],data[[var.cluster]],unique),
                index.time = as.numeric(data[[var.time]]),
                cluster = list(n = n.cluster, levels = U.cluster, nobs = table(index.cluster)),
                param = list(mu = colnames(X.mean),sigma = param.sigma, k = param.k, rho = param.rho)))
}

##----------------------------------------------------------------------
### model.matrix.R ends here
