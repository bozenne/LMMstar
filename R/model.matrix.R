### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: Apr 24 2021 (22:47) 
##           By: Brice Ozenne
##     Update #: 447
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
model.matrix.lmm <- function(object, data = NULL, effects = "all", type.object = "lmm"){

    ## ** normalize user imput
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)

    ## ** update design matrix with new dataset
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
    }else{
        design <- object$design
    }
    
    ## ** update design matrix with new dataset
    if(type.object == "lmm"){
        if("mean" %in% effects && "variance" %in% effects){
            return(list(mean = design$X.mean,
                        variance = design$X.var))
        }else if("mean" %in% effects){
            return(design$X.mean)
        }else if("variance" %in% effects){
            return(design$X.var)
        }
    }else if(type.object == "gls"){
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
                              structure){
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
    attr(X.var,"Upattern.time") <- lapply(attr(X.var,"UX.var"), rownames) ## timepoints associated to each pattern
    attr(X.var,"Upattern.index.time") <- lapply(attr(X.var,"Upattern.time"), function(iVec){which(U.time %in% iVec)}) ## index of the timepoints associated to each pattern
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
                iIndicator[[param.sigma[iSigma]]] <- tcrossprod(ind.var_dsigma) > 0
            }
        }
        if(length(param.k)>0){
            for(iK in 1:length(param.k)){ ## iK <- 1
                ## positions where the k-parameter appears in the matrix
                dk <- rep(0,length(param.k))
                dk[iK] <- 1
                ind.var_dk <- attr(X.var,"UX.var")[[iP]] %*% c(rep(0,length(sigma)),dk)
                iIndicator[[param.k[iK]]] <- (ind.var_dk %*% rep(1,iNtime)) + t(ind.var_dk %*% rep(1,iNtime)) > 0
            }
        }
        if(length(param.rho)>0){
            for(iRho in 1:length(param.rho)){ ## iRho <- 1
                ## positions where the rho-parameter appears in the matrix
                iIndicator[[param.rho[iRho]]] <- matrix(attr(X.var,"UX.cor")[[iP]][,iRho] > 0, nrow = iNtime, ncol = iNtime, dimnames = list(iTime,iTime))
            }
        }
        return(iIndicator)
    }), attr(X.var,"Upattern"))

    ## which parameters are involved in the patterm    
    attr(X.var, "Upattern.param") <- setNames(lapply(attr(X.var,"indicator"), function(iP){ ## iP <- 1
        names(which(sapply(iP,any)))
    }), attr(X.var,"Upattern"))

    ## ** export
    return(list(X.mean = X.mean,
                X.var = X.var,
                Y = data[[var.outcome]],
                index.cluster = index.cluster,
                index.vargroup = index.vargroup,
                index.strata = tapply(data[[var.strata]],data[[var.cluster]],unique),
                index.time = as.numeric(data[[var.time]]),
                cluster = list(n = n.cluster, levels = U.cluster, nobs = table(index.cluster)),
                param = list(mu = colnames(X.mean), sigma = param.sigma, k = param.k, rho = param.rho,
                             pair.meanvarcoef = unname(t(expand.grid(colnames(X.mean),c(param.sigma,param.k,param.rho)))),
                             pair.varcoef = .unorderedPairs(c(param.sigma,param.k,param.rho)))))
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
