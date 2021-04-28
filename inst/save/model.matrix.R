### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: Apr 21 2021 (17:33) 
##           By: Brice Ozenne
##     Update #: 424
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

    if(transform %in% 0:2 == FALSE){
        stop("Argument \'transform\' must be 0, 1, or 2 in .model.matrix.lmm. \n")
    }

    ## ** normalize data
    if("XXindexXX" %in% names(data) == FALSE){
        data$XXindexXX <- 1:NROW(data)
    }
    if(identical(var.strata, "XXstrata.indexXX") && var.strata %in% names(data) == FALSE){
        data$XXstrata.indexXX <- 1
    }else{
        data[[var.strata]] <- factor(data[[var.strata]], levels = U.strata)
    }
    
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
    data[[var.time]] <- factor(data[[var.time]], levels = U.time)

    ## *** parametrisation
    ## **** sigma
    if(n.strata==1){
        param.sigma <- "sigma"
    }else{
        param.sigma <- paste0("sigma",":",U.strata)
    }

    ## **** k
    if(structure %in% c("CS","EXP")){
        if(n.strata==1){
            X.var <- model.matrix(~1, data)
            colnames(X.var) <- param.sigma            
        }else{
            X.var <- model.matrix(as.formula(paste0("~0+",var.strata)), data)
            colnames(X.var) <- param.sigma
        }
        param.k <- NULL
    }else if(structure == "UN"){
        data.relevel <- as.data.frame(data)
        if(n.strata==1){
            for(iVar in all.vars(formula.var)){ ## iVar <- "Gender"
                data.relevel[[iVar]] <- factor(data.relevel[[iVar]], levels = unique(data.relevel[[iVar]])) ## to match with gls which chooses the reference level according to the ordering
            }
            X.var <- model.matrix(formula.var, data.relevel)
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
    attr(X.var,"transform") <- transform ## should parameters be rescaled (0: sd, cor | 1 sd(log transform), cor(atanh transform), 2: var, cor)
    attr(X.var,"pattern") <- tapply(attr(X.var,"level"), index.cluster, paste, collapse=".") ## automatic re-ordering of the result
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
    attr(X.var,"Upattern.index.time") <- lapply(attr(X.var,"Upattern.time"), function(iVec){which(U.time %in% iVec)})
    index.vargroup <- setNames(match(attr(X.var,"pattern"),attr(X.var,"Upattern")),names(attr(X.var,"pattern")))

    ## **** rho
    if(n.time==1 || all(unlist(lapply(attr(X.var,"Upattern.time"),length))==1)){
        param.rho <- NULL
    }else{
        Mtempo <- matrix(1:16,4,4)
        tMtempo <- t(Mtempo)
        if(structure == "CS"){
            param.rho <- "Rho"
            M.index <- cbind(index.lower = Mtempo[lower.tri(Mtempo)],
                             index.upper = tMtempo[lower.tri(Mtempo)])
        }else if(structure == "EXP"){
            param.rho <- "range"
            browser()
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
                    browser()
                    if(n.strata==1){
                        X.cor[[iStrata]][M.index[iParam,c("index.lower","index.upper")],param.rho[iParam]] <- 1
                    }else{
                        X.cor[[iStrata]][M.index[iParam,c("index.lower","index.upper")],paste(param.rho.save[iParam],U.strata[iStrata],sep=":")] <- 1
                    }
                }
            }
        }

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

    ## ** prepare calculation of the variance matrices
    attr(X.var, "FUN.Omega") <- function(object, sigma, k, rho, keep.interim = FALSE){
        out <- lapply(1:attr(object,"nUpattern"), function(iP){

            iTime <- attr(object,"Upattern.time")[[iP]]
            iNtime <- length(iTime)
            Omega.sd <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,k)))[,1])
            if(length(rho)>0){
                Omega.cor <- matrix(attr(object,"UX.cor")[[iP]] %*% rho, nrow = iNtime, ncol = iNtime, dimnames = list(iTime,iTime))
            }else{
                Omega.cor <- diag(0, nrow = iNtime, ncol = iNtime)
            }
            Omega <- diag(Omega.sd^2, nrow = iNtime, ncol = iNtime) + Omega.cor * tcrossprod(Omega.sd)
            dimnames(Omega) <- list(iTime,iTime)
            if(keep.interim){
                attr(Omega,"sd") <- Omega.sd
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
         
        if(length(param.sigma)>0){
            for(iSigma in 1:length(param.sigma)){ ## iSigma <- 1
                ## positions where the sigma-parameter appears in the matrix
                dsigma <- rep(0,length(param.sigma))
                dsigma[iSigma] <- 1
                ind.var_dsigma <- attr(X.var,"UX.var")[[iP]] %*% c(dsigma,rep(0, length(param.k)))
                iIndicator[[param.sigma[iSigma]]] <- tcrossprod(ind.var_dsigma)
            }
        }
        if(length(param.k)>0){
            for(iK in 1:length(param.k)){ ## iK <- 1
                ## positions where the k-parameter appears in the matrix
                dk <- rep(0,length(param.k))
                dk[iK] <- 1
                ind.var_dk <- attr(X.var,"UX.var")[[iP]] %*% c(rep(0,length(sigma)),dk)
                iIndicator[[param.k[iK]]] <- ((ind.var_dk %*% rep(1,iNtime)) + t(ind.var_dk %*% rep(1,iNtime))) > 0
            }
        }
        if(length(param.rho)>0){
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
            iOmega.sd <- attr(Omega[[iP]],"sd")
            iOmega.cor <- attr(Omega[[iP]],"cor")
            iOmega <- Omega[[iP]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL;

            if(length(sigma)>0){
                for(iSigma in 1:n.sigma){ ## iSigma <- 1
                    ## positions where the k-parameter appears
                    ind.dsigma <- attr(object,"UX.var")[[iP]][,name.sigma[iSigma]]
                    ## compute derivative
                    idsigma <- sigma
                    idsigma[name.sigma[iSigma]] <- 1
                    idk <- k
                    if(transform==2 && length(k)>0){
                        iK <- which(colSums(object[object[,name.sigma[iSigma]]==1,names(k),drop=FALSE])>0)
                        idk[name.k[iK]] <- 1
                    }
                    ## propagate
                    idOmega.sigma <- exp((attr(object,"UX.var")[[iP]] %*% log(c(idsigma,idk)))[,1]) * ind.dsigma
                    iScore[[name.sigma[iSigma]]] <- diag(2*idOmega.sigma*iOmega.sd, nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.sigma %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.sigma))
                    ## iScore[[name.sigma[iSigma]]] - 2*iOmega/sigma

                    ## transformation
                    ## d sigma^2         -> 2 dsigma sigma
                    ## d exp(2 logsigma) -> 2 dlogsigma exp(2 logsigma) = 2 dlogsigma sigma^2 --> sigma missing
                    ## d sigma2          -> dsigma2 --> 2 sigma too much
                    if(transform == 1){
                        iScore[[name.sigma[iSigma]]] <- iScore[[name.sigma[iSigma]]] * sigma[name.sigma[iSigma]]
                    }else if(transform == 2){
                        ind.dsigma <- 1-do.call(prod,attr(object,"indicator")[[iP]][names(iK)])
                        iScore[[name.sigma[iSigma]]] <- ind.dsigma * iScore[[name.sigma[iSigma]]] / (2 * sigma[name.sigma[iSigma]])
                    }
                }
            }
            if(length(k)>0){
                for(iK in 1:n.k){ ## iK <- 1
                    ## position where the k-parameter appears
                    ind.dk <- attr(object,"UX.var")[[iP]][,name.k[iK]]
                    ## compute derivative
                    idk <- k
                    idk[name.k[iK]] <- 1
                    idsigma <- sigma
                    if(transform==2){
                        iSigma <- which(colSums(object[object[,name.k[iK]]==1,names(sigma),drop=FALSE])>0)
                        idsigma[name.sigma[iSigma]] <- 1
                    }
                    ## propagate
                    idOmega.k <- exp((attr(object,"UX.var")[[iP]] %*% log(c(idsigma,idk)))[,1]) * ind.dk
                    iScore[[name.k[iK]]] <- diag(2*idOmega.k*iOmega.sd, nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.k %*% t(iOmega.sd) + iOmega.sd %*% t(idOmega.k))
                    ## iScore[[name.k[iK]]] - (tcrossprod(ind.dk,rep(1,4))+tcrossprod(rep(1,4),ind.dk)) * iOmega/k[name.k[iK]]
                    
                    ## transformation
                    ## d k^2         -> 2 dk k
                    ## d exp(2 logk) -> 2 dlogk exp(2 logk) = 2 dlogk k^2 --> k missing
                    ## d k2sigma2 or d sqrt(k2sigma2)   -> 1 or 1/(2ksigma) --> 2 * k sigma2 too much
                    if(transform == 1){ 
                        iScore[[name.k[iK]]] <- iScore[[name.k[iK]]] * k[name.k[iK]]
                    }else if(transform == 2){ ## NOTE: the score relative to sigma has already been transformed to sigma2
                        iScore[[name.k[iK]]] <- iScore[[name.k[iK]]] / (2 * k[name.k[iK]] * sigma[name.sigma[iSigma]]) 
                    }
                }
            }
            if(length(rho)>0){
                for(iRho in 1:n.rho){ ## iRho <- 1
                    ## positions where the rho-parameter appears
                    ind.drho <- attr(object,"indicator")[[iP]][[name.rho[iRho]]]
                    ## derivative
                    iScore[[name.rho[iRho]]] <- ind.drho * tcrossprod(iOmega.sd)
                    ## iOmega/iOmega.corn

                    ## transformation
                    ## d rho            --> drho  --> rho too much                                
                    ## d tanh(atanhrho) --> datanhrho * (1-tanh(atanhrho)^2) = datanhrho(1-rho^2) --> rho too much while (1-rho^2) missing                                
                    if(transform == 1){ ## atanh transform
                        iScore[[name.rho[iRho]]] <- iScore[[name.rho[iRho]]] * (1-rho[name.rho[iRho]]^2)
                    }
                }
            }
            return(iScore)
       })
        return(setNames(out,attr(object,"Upattern")))
    }

    ## ** prepare calculation of the df
    attr(X.var, "FUN.d2Omega") <- function(object, sigma, k, rho, Omega, dOmega, pair){
        name.sigma <- names(sigma)
        n.sigma <- length(sigma)
        name.k <- names(k)
        n.k <- length(k)
        name.rho <- names(rho)
        n.rho <- length(rho)
        allCoef <- c(name.sigma, name.k, name.rho)
        p <- length(allCoef)
        transform <- attr(object,"transform")
        n.pair  <- NCOL(pair)
        
        out <- lapply(1:attr(object,"nUpattern"), function(iP){
            iHess <- vector(mode = "list", length = n.pair)

            iTime <- attr(object,"Upattern.time")[[iP]]
            iNtime <- length(iTime)
            iOmega.sd <- attr(Omega[[iP]],"sd")
            iOmega.cor <- attr(Omega[[iP]],"cor")
            iOmega <- Omega[[iP]] ; attr(iOmega,"sd") <- NULL; attr(iOmega,"cor") <- NULL;
            idOmega <- dOmega[[iP]]

            for(iPair in 1:n.pair){ ## iPair <- 1
                ## name of parameters
                iCoef1 <- pair[1,iPair]
                iCoef2 <- pair[2,iPair]

                ## type of parameters
                iType1 <- c("sigma","k","rho")[which(c(iCoef1 %in% name.sigma,iCoef1 %in% name.k,iCoef1 %in% name.rho))]
                iType2 <- c("sigma","k","rho")[which(c(iCoef2 %in% name.sigma,iCoef2 %in% name.k,iCoef2 %in% name.rho))]

                ## positions where the parameters appears in the matrix
                ind1 <- attr(object,"indicator")[[iP]][[iCoef1]]
                ind2 <- attr(object,"indicator")[[iP]][[iCoef2]]
                    
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
                        ind.dksigma <- attr(object,"UX.var")[[iP]][,iCoef1] * attr(object,"UX.var")[[iP]][,iCoef2]
                        ## compute derivative
                        idsigma <- sigma; idsigma[iCoef1] <- 1
                        idk <- k; idk[iCoef2] <- 1
                        ## propagate
                        idOmega.k <- exp((attr(object,"UX.var")[[iP]] %*% log(c(idsigma,idk)))[,1]) * ind.dksigma
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
                        ind.dk1 <- attr(object,"UX.var")[[iP]][,iCoef1]
                        ind.dk2 <- attr(object,"UX.var")[[iP]][,iCoef2]
                        ## compute derivative
                        idk1 <- k ; idk1[iCoef1] <- 1
                        idk2 <- k ; idk2[iCoef2] <- 1
                        ## propagate
                        idOmega.k1 <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,idk1)))[,1]) * ind.dk1
                        idOmega.k2 <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,idk2)))[,1]) * ind.dk2
                        iHess[[iPair]] <- diag(2*idOmega.k1*idOmega.k2, nrow = iNtime, ncol = iNtime) + iOmega.cor * (idOmega.k1 %*% t(idOmega.k2) + idOmega.k2 %*% t(idOmega.k1))
                        
                        ## transformation
                        ## d2 k^2  or d k d k'       = 2 k or 1
                        ## d2 exp(2 k)  or d exp(k) d exp(k') =  2 k^2 or k k'
                        ## d2 d (sigma2 k2) or d (sigma2 k2) d (sigma2' k2') = ????
                        if(transform == 1){ 
                            iHess[[iPair]] <- iHess[[iPair]] * (1+iCoef1==iCoef2) * k[iCoef1] * k[iCoef2]
                        }else if(transform == 2){
                            browser()
                        }
                    
                    }else if(iType2 == "rho"){
                        ## position where the parameters appear
                        ind.dk <- attr(object,"UX.var")[[iP]][,iCoef1]
                        ## compute derivative
                        idk <- k; idk[iCoef1] <- 1
                        ## propagate
                        idOmega.k <- exp((attr(object,"UX.var")[[iP]] %*% log(c(sigma,idk)))[,1]) * ind.dk
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
                param = list(mu = colnames(X.mean),sigma = param.sigma, k = param.k, rho = param.rho, pair.varcoef = .unorderedPairs(c(param.sigma,param.k,param.rho)))))
}

## * .unorderedPairs
## adapted from RecordLinkage package
.unorderedPairs <- function(x){
    n <- length(x)
    ls <- lapply(1:n, function(k){ rbind(x[k], x[k:n])})
    out <- array(unlist(ls), dim = c(2, n * (n + 1)/2))
    return(out)
}
##----------------------------------------------------------------------
### model.matrix.R ends here
