### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: May 27 2021 (15:46) 
##           By: Brice Ozenne
##     Update #: 653
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
model.matrix.lmm <- function(object, data = NULL, effects = "all", type.object = "lmm", ...){

    ## ** normalize user imput
    type.object <- match.arg(type.object, c("lmm","gls"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

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
    index.strata <- tapply(data[[var.strata]],data[[var.cluster]],unique)
    
    ## ** mean
    if(n.strata==1){
        X.mean <- model.matrix(formula.mean, data)
        strata.mu <- setNames(rep(1,NCOL(X.mean)), colnames(X.mean))
    }else{
        ls.X.mean <- lapply(U.strata, function(iS){ ## iS <- U.strata[1]
            iX <- model.matrix(formula.mean, data[data[[var.strata]]==iS,])
            colnames(iX) <- paste0(colnames(iX),":",iS)
            attr(iX,"index") <- data[data[[var.strata]]==iS,"XXindexXX"]
            return(iX)
        })

        X.mean <- as.matrix(Matrix::bdiag(ls.X.mean))[order(unlist(lapply(ls.X.mean, attr, "index"))),]
        colnames(X.mean) <- unlist(lapply(ls.X.mean,colnames))
        attr(X.mean, "assign") <- as.vector(do.call(cbind,lapply(ls.X.mean,attr,"assign")))

        strata.mu <- unlist(lapply(1:n.strata, function(iStrata){setNames(rep(iStrata, NCOL(ls.X.mean[[iStrata]])),colnames(ls.X.mean[[iStrata]]))}))
    }

    ## ** cluster
    U.cluster <- sort(unique(data[[var.cluster]]))
    n.cluster <- length(U.cluster)
    index.cluster <- match(data[[var.cluster]], U.cluster) ## ‘match’ returns a vector of the positions of (first) matches of its first argument in its second.

    ## ** variance
    data[[var.time]] <- factor(data[[var.time]], levels = U.time)
    index.time <- as.numeric(data[[var.time]])

    attr(index.cluster,"sorted") <- lapply(1:n.cluster, function(iId){
        iIndex <- which(index.cluster==iId)
        return(iIndex[order(index.time[iIndex])]) ## re-order observations according to the variance-covariance matrix
    })
    
    ## *** parametrisation
    ## **** sigma
    if(n.strata==1){
        param.sigma <- "sigma"
        strata.sigma <- setNames(1,param.sigma)
    }else{
        param.sigma <- paste0("sigma",":",U.strata)
        strata.sigma <- setNames(1:n.strata,param.sigma)
    }
    ## **** k
    if(structure %in% c("CS","EXP")){
        X.var <- model.matrix(formula.var, data)
        colnames(X.var) <- param.sigma
        attr(X.var,"assign") <- c(0,0)
        param.k <- NULL
        time.k <- NULL
        strata.k <- NULL
    }else if(structure == "UN"){
        data.relevel <- as.data.frame(data)
        if(n.strata>1){
            formula.var <- as.formula(paste0("~-1+",var.time,":",var.strata,"+",var.strata))
        }
        X.var <- model.matrix(formula.var, data.relevel)

        ## extract sigma columns
        index.sigma <- which(attr(X.var,"assign")<=(n.strata>1))
        if(n.strata==1){
            colnames(X.var)[index.sigma] <- param.sigma
        }else if(n.strata>1){
            ls.order.strata <- lapply(1:length(U.strata), function(iStrata){## iStrata <- 1
                c(which(colnames(X.var)==paste0(var.strata,U.strata[iStrata])),
                  grep(paste0("\\:",var.strata,U.strata[iStrata]),colnames(X.var)))
            })
            order.strata <-  unlist(ls.order.strata)
            colnames(X.var) <- gsub(paste0("^",var.strata),"",gsub(paste0("\\:",var.strata),"\\:",colnames(X.var)))
            colnames(X.var)[index.sigma] <- paste0("sigma:",colnames(X.var)[index.sigma])
        }

        ## extract k columns
        index.k <- which(attr(X.var,"assign")>(n.strata>1))
        time.k <- rep(as.character(NA), length(index.k))
        if(n.strata==1){
            for(iTime in 1:n.time){ ## iTime <- 2
                iPattern <- paste0("^",var.time,U.time[iTime],"$")
                time.k[grepl(pattern = iPattern, colnames(X.var)[index.k])] <- U.time[iTime]
            }
            colnames(X.var)[index.k] <- paste0("k.",gsub(paste0("^",var.time),"",colnames(X.var)[index.k]))
        }else if(n.strata>1){
            for(iTime in 1:n.time){ ## iTime <- 2
                iPattern <- paste0("^",var.time,U.time[iTime],"\\:")
                time.k[grepl(pattern = iPattern, colnames(X.var)[index.k])] <- U.time[iTime]
            }
            colnames(X.var)[index.k] <- paste0("k.",gsub(paste0("^",var.time),"",colnames(X.var)[index.k]))
        }
        param.k <- colnames(X.var)[index.k]

        if(n.strata>1){ ## re-order such that the coef for strata 1 appear first then those of strata 2
            strata.k <- unlist(lapply(1:n.strata, function(iStrata){setNames(rep(iStrata, length(ls.order.strata[[iStrata]])),colnames(X.var)[ls.order.strata[[iStrata]]])}))
            strata.k <- strata.k[names(strata.k) %in% param.sigma == FALSE]
            X.var <- X.var[,order.strata,drop=FALSE]
            param.k <- colnames(X.var)[colnames(X.var) %in% param.k]
        }else{
            strata.k <- setNames(rep(1,length(param.k)),param.k)
        }

    }

    ## uniaue variance levels
    level.variance <- as.numeric(droplevels(interaction(as.data.frame(X.var))))
    pattern <- tapply(level.variance, index.cluster, paste, collapse=".") ## automatic re-ordering of the result
    ## characterize each level
    Upattern <- unname(sort(unique(pattern)))
    n.pattern <- length(Upattern)
    
    indexPattern.tempo <- which(duplicated(pattern)==FALSE) ## first cluster corresponding to each unique pattern
    indexCluster.tempo <- names(sort(pattern[indexPattern.tempo])) ## corresponding cluster
    indexObs.tempo <- lapply(indexCluster.tempo, function(iCluster){which(index.cluster==iCluster)})

    strata.Upattern <- setNames(sapply(indexObs.tempo, function(iObs){unique(data[iObs,var.strata])}), Upattern) ## strata associated to each unique pattern
    X.Upattern <- setNames(lapply(indexObs.tempo, function(iObs){ ## extract design matrix
        iStrata <- which(U.strata %in% data[[var.strata]][iObs])
        iKeep.param <- c(names(strata.sigma)[strata.sigma==iStrata],names(strata.k)[strata.k==iStrata])
        iM <- X.var[iObs,iKeep.param,drop=FALSE]
        rownames(iM) <- data[[var.time]][iObs]
        return(iM)
    }), Upattern)
    time.Upattern <- lapply(X.Upattern, rownames) ## timepoints associated to each pattern
    indexTime.Upattern <- lapply(time.Upattern, function(iVec){which(U.time %in% iVec)}) ## index of the timepoints associated to each pattern

    ## **** rho
    if(n.time==1 || all(unlist(lapply(time.Upattern,length))==1)){
        param.rho <- NULL
        time.rho <- NULL
        strata.rho <- NULL
        X.cor <- NULL
    }else{
        if(structure == "CS"){
            param.rho <- "Rho"
            time.rho <- rbind(U.time[1],U.time[2])
        }else if(structure == "EXP"){
            param.rho <- "range"
            time.rho <- NULL
        }else if(structure == "UN"){
            Mtempo <- matrix(1:n.time^2, nrow = n.time, ncol = n.time,
                             dimnames = list(U.time,U.time))
            M.index <- cbind(which(lower.tri(diag(1, nrow = n.time, ncol = n.time)), arr.ind = TRUE),
                             index.lower = Mtempo[lower.tri(Mtempo)],
                             index.upper = t(Mtempo)[lower.tri(Mtempo)])
            param.rho <- paste0("cor","(",U.time[M.index[,"row"]],",",U.time[M.index[,"col"]],")")
            time.rho <- rbind(U.time[M.index[,"row"]],U.time[M.index[,"col"]])
            colnames(time.rho) <- param.rho
        }
        param.rho.save <- param.rho
        strata.rho <- setNames(rep(1,length(param.rho)),param.rho)

        if(n.strata>1){
            ls.param.rho <- lapply(U.strata, function(iS){paste(param.rho.save, iS, sep = ":")})
            param.rho <- unlist(ls.param.rho)
            strata.rho <- unlist(lapply(1:n.strata, function(iS){setNames(rep(iS,length(ls.param.rho[[iS]])),ls.param.rho[[iS]])}))
            time.rho <- do.call(cbind,lapply(U.strata, function(iS){time.rho}))
            colnames(time.rho) <- param.rho
        }
        
        M.indexAlltimes <- which(matrix(1,n.time, n.time)==1, arr.ind = TRUE)
        
        X.cor <- setNames(vector(mode = "list", length = n.pattern), Upattern)
        for(iPattern in Upattern){ ## iPattern <- Upattern[1]
            iStrata <- strata.Upattern[[iPattern]]
            iTime <- time.Upattern[[iPattern]]
            iParam.rho <- param.rho[U.strata[strata.rho]==iStrata]
            iN.time <- length(iTime)
            iM.indexAlltimes <-  M.indexAlltimes[(M.indexAlltimes[,"row"] %in% indexTime.Upattern[[iPattern]])*(M.indexAlltimes[,"col"] %in% indexTime.Upattern[[iPattern]])>0,,drop=FALSE]
            
            if(structure == "CS"){
                X.cor[[iPattern]] <- matrix(as.numeric(iM.indexAlltimes[,"row"] != iM.indexAlltimes[,"col"]), nrow = iN.time^2, ncol = length(iParam.rho),
                                            dimnames = list(paste0("(",iM.indexAlltimes[,"row"],",",iM.indexAlltimes[,"col"],")"), iParam.rho))
            }else if(structure == "UN"){
                X.cor[[iPattern]] <- matrix(0, nrow = iN.time^2, ncol = length(iParam.rho),
                                            dimnames = list(paste0("(",iM.indexAlltimes[,"row"],",",iM.indexAlltimes[,"col"],")"), iParam.rho))
                for(iParam in iParam.rho){ ## iParam <- iParam.rho[1]
                    X.cor[[iStrata]][paste0("(",paste(which(U.time %in% time.rho[,iParam]),collapse=","),")"),iParam] <- 1
                    X.cor[[iStrata]][paste0("(",paste(rev(which(U.time %in% time.rho[,iParam])),collapse=","),")"),iParam] <- 1
                }
            }
        }
    }

    ## ** prepare calculation of the score
    ## indicator of the effect of each parameter on each element of the covariance
    indicator.param <- setNames(lapply(Upattern, function(iP){ ## iP <- 1

        iStrata <- which(U.strata == strata.Upattern[iP])
        iParam.sigma <- names(strata.sigma)[strata.sigma==iStrata]
        iParam.k <- names(strata.k)[strata.k==iStrata]
        iParam.rho <- names(strata.rho)[strata.rho==iStrata]
        iParam.all <- c(iParam.sigma,iParam.k,iParam.rho)
        n.iParam.sigma <- length(iParam.sigma)
        n.iParam.k <- length(iParam.k)
        n.iParam.rho <- length(iParam.rho)
        n.iParam.all <- length(iParam.all)
        
        iIndicator <- setNames(vector(mode = "list", length = n.iParam.all), iParam.all)
        iTime <- time.Upattern[[iP]]
        iNtime <- length(iTime)

        if(n.iParam.sigma>0){
            for(iSigma in iParam.sigma){ ## iSigma <- 2
                ## positions where the sigma-parameter appears in the matrix
                ind.var_dsigma <- X.Upattern[[iP]][,iSigma,drop=FALSE]
                iIndicator[[iSigma]] <- tcrossprod(ind.var_dsigma) > 0
            }
        }
        if(n.iParam.k>0){
            for(iK in iParam.k){ ## iK <- iParam.k[1]
                ## positions where the k-parameter appears in the matrix
                ind.var_dk <- X.Upattern[[iP]][,iK,drop=FALSE]
                iIndicator[[iK]] <- (ind.var_dk %*% rep(1,iNtime)) + t(ind.var_dk %*% rep(1,iNtime)) > 0
            }
        }
        if(n.iParam.rho>0){
            for(iRho in iParam.rho){ ## iRho <- 1
                ## positions where the rho-parameter appears in the matrix
                iIndicator[[iRho]] <- matrix(X.cor[[iP]][,iRho] > 0, nrow = iNtime, ncol = iNtime, dimnames = list(iTime,iTime))
            }
        }
        return(iIndicator)
    }), Upattern)

    ## which parameters are involved in the patterm    
    param.Upattern <- setNames(lapply(indicator.param, function(iP){ ## iP <- 1
        names(which(sapply(iP,any)))
    }), Upattern)

    X.Upattern <- setNames(lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]
        iParam <- intersect(param.Upattern[[iPattern]], c(param.sigma, param.k))
        iXname <- colnames(X.Upattern[[iPattern]])
        iExtraName <- iXname[iXname %in% iParam == FALSE]
        if(identical(iParam, iXname)){
            return(X.Upattern[[iPattern]])
        }else if(all(iParam %in% iXname) && all(X.Upattern[[iPattern]][,iExtraName]==0) ){
            return(X.Upattern[[iPattern]][,iParam,drop=FALSE])
        }else{
            stop("Something went wrong when creating the design matrix for the variable parameters. \n",
                 "Contact the package manager with a reproducible example generating this error message. \n")
        }
    }), Upattern)
    attr(X.Upattern,"assign") <- attr(X.var,"assign")

    ## ** pairs
    pair.meanvarcoef <- setNames(lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]
        unname(t(expand.grid(names(strata.mu)[U.strata[strata.mu]==strata.Upattern[[iPattern]]], param.Upattern[[iPattern]])))
    }), Upattern)

    pair.varcoef <- setNames(lapply(Upattern, function(iPattern){## ## iPattern <- Upattern[1]
        .unorderedPairs(param.Upattern[[iPattern]])
    }),Upattern)

    if(length(param.rho)>0){
        attr(pair.varcoef,"test.hessian") <- setNames(lapply(pair.varcoef, function(iM){colSums(apply(iM, 2, `%in%`, param.rho))!=2}), Upattern)
    }else{
        attr(pair.varcoef,"test.hessian") <- setNames(lapply(pair.varcoef, function(iM){rep(TRUE,NCOL(iM))}), Upattern)
    }

    ## ** gather and export
    return(list(X.mean = X.mean,
                X.var = list(var = X.Upattern,
                             cor = X.cor,
                             pattern = Upattern,
                             strata =  strata.Upattern,
                             index.time = indexTime.Upattern,
                             cluster = setNames(match(pattern,Upattern),names(pattern)),
                             param =  param.Upattern,
                             indicator = indicator.param),
                Y = data[[var.outcome]],
                index.cluster = index.cluster,
                index.time = index.time,
                cluster = list(n = n.cluster, levels = U.cluster, nobs = table(index.cluster)),
                param = list(mu = colnames(X.mean), sigma = param.sigma, k = param.k, rho = param.rho,
                             strata.mu = strata.mu, strata.sigma = strata.sigma, strata.k = strata.k, strata.rho = strata.rho,
                             time.k = time.k, time.rho = time.rho,
                             pair.meanvarcoef = pair.meanvarcoef,
                             pair.varcoef = pair.varcoef)
                ))
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
