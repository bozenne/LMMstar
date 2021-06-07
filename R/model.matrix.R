### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: Jun  7 2021 (15:18) 
##           By: Brice Ozenne
##     Update #: 720
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
        ff.allvars <- c(all.vars(object$formula$mean), all.vars(object$formula$var))
        if(any(ff.allvars %in% names(data) == FALSE)){
            stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
        }

        design <- .model.matrix.lmm(formula.mean = object$formula$mean.design,
                                    formula.var = object$formula$var.design,
                                    data = data,
                                    var.outcome = object$outcome$var,
                                    var.strata = object$strata$var, U.strata = object$strata$levels,
                                    var.time = object$time$var, U.time = object$time$levels,
                                    var.cluster = object$cluster$var,
                                    structure = object$structure
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
        X.mean <- model.matrix_regularize(formula.mean, data)
        strata.mu <- stats::setNames(rep(1,NCOL(X.mean)), colnames(X.mean))
    }else{
        ls.X.mean <- lapply(U.strata, function(iS){ ## iS <- U.strata[1]
            iX <- model.matrix_regularize(formula.mean, data[data[[var.strata]]==iS,])
            colnames(iX) <- paste0(colnames(iX),":",iS)
            attr(iX,"index") <- data[data[[var.strata]]==iS,"XXindexXX"]
            return(iX)
        })

        X.mean <- as.matrix(Matrix::bdiag(ls.X.mean))[order(unlist(lapply(ls.X.mean, attr, "index"))),]
        colnames(X.mean) <- unlist(lapply(ls.X.mean,colnames))
        attr(X.mean, "assign") <- as.vector(do.call(cbind,lapply(ls.X.mean,attr,"assign")))

        strata.mu <- unlist(lapply(1:n.strata, function(iStrata){stats::setNames(rep(iStrata, NCOL(ls.X.mean[[iStrata]])),colnames(ls.X.mean[[iStrata]]))}))
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
        strata.sigma <- stats::setNames(1,param.sigma)
    }else{
        param.sigma <- paste0("sigma",":",U.strata)
        strata.sigma <- stats::setNames(1:n.strata,param.sigma)
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
            formula.var <- stats::as.formula(paste0("~-1+",var.time,":",var.strata,"+",var.strata))
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
            strata.k <- unlist(lapply(1:n.strata, function(iStrata){stats::setNames(rep(iStrata, length(ls.order.strata[[iStrata]])),colnames(X.var)[ls.order.strata[[iStrata]]])}))
            strata.k <- strata.k[names(strata.k) %in% param.sigma == FALSE]
            tempo <- attr(X.var,"assign")-1
            X.var <- X.var[,order.strata,drop=FALSE]
            attr(X.var,"assign") <- tempo
            param.k <- colnames(X.var)[colnames(X.var) %in% param.k]
        }else{
            strata.k <- stats::setNames(rep(1,length(param.k)),param.k)
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

    strata.Upattern <- stats::setNames(sapply(indexObs.tempo, function(iObs){unique(data[iObs,var.strata])}), Upattern) ## strata associated to each unique pattern
    X.Upattern <- stats::setNames(lapply(indexObs.tempo, function(iObs){ ## extract design matrix
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
        strata.rho <- stats::setNames(rep(1,length(param.rho)),param.rho)

        if(n.strata>1){
            ls.param.rho <- lapply(U.strata, function(iS){paste(param.rho.save, iS, sep = ":")})
            param.rho <- unlist(ls.param.rho)
            strata.rho <- unlist(lapply(1:n.strata, function(iS){stats::setNames(rep(iS,length(ls.param.rho[[iS]])),ls.param.rho[[iS]])}))
            time.rho <- do.call(cbind,lapply(U.strata, function(iS){time.rho}))
            colnames(time.rho) <- param.rho
        }
        
        M.indexAlltimes <- which(matrix(1,n.time, n.time)==1, arr.ind = TRUE)
        
        X.cor <- stats::setNames(vector(mode = "list", length = n.pattern), Upattern)
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
        attr(X.cor,"assign") <- rep(1,length(param.rho))
    }

    ## ** prepare calculation of the score
    ## indicator of the effect of each parameter on each element of the covariance
    indicator.param <- stats::setNames(lapply(Upattern, function(iP){ ## iP <- 1

        iStrata <- which(U.strata == strata.Upattern[iP])
        iParam.sigma <- names(strata.sigma)[strata.sigma==iStrata]
        iParam.k <- names(strata.k)[strata.k==iStrata]
        iParam.rho <- names(strata.rho)[strata.rho==iStrata]
        iParam.all <- c(iParam.sigma,iParam.k,iParam.rho)
        n.iParam.sigma <- length(iParam.sigma)
        n.iParam.k <- length(iParam.k)
        n.iParam.rho <- length(iParam.rho)
        n.iParam.all <- length(iParam.all)
        
        iIndicator <- stats::setNames(vector(mode = "list", length = n.iParam.all), iParam.all)
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
    param.Upattern <- stats::setNames(lapply(indicator.param, function(iP){ ## iP <- 1
        names(which(sapply(iP,any)))
    }), Upattern)

    X.Upattern <- stats::setNames(lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]
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
    pair.meanvarcoef <- stats::setNames(lapply(Upattern, function(iPattern){ ## iPattern <- Upattern[1]
        unname(t(expand.grid(names(strata.mu)[U.strata[strata.mu]==strata.Upattern[[iPattern]]], param.Upattern[[iPattern]])))
    }), Upattern)

    pair.varcoef <- stats::setNames(lapply(Upattern, function(iPattern){## ## iPattern <- Upattern[1]
        .unorderedPairs(param.Upattern[[iPattern]])
    }),Upattern)

    if(length(param.rho)>0){
        attr(pair.varcoef,"test.hessian") <- stats::setNames(lapply(pair.varcoef, function(iM){colSums(apply(iM, 2, `%in%`, param.rho))!=2}), Upattern)
    }else{
        attr(pair.varcoef,"test.hessian") <- stats::setNames(lapply(pair.varcoef, function(iM){rep(TRUE,NCOL(iM))}), Upattern)
    }

    ## ** gather and export
    return(list(X.mean = X.mean,
                X.var = list(var = X.Upattern,
                             cor = X.cor,
                             pattern = Upattern,
                             strata =  strata.Upattern,
                             index.time = indexTime.Upattern,
                             cluster = stats::setNames(match(pattern,Upattern),names(pattern)),
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
## form all combinations
.unorderedPairs <- function(x){
    n <- length(x)
    ls <- lapply(1:n, function(k){ rbind(x[k], x[k:n])})
    out <- array(unlist(ls), dim = c(2, n * (n + 1)/2))
    return(out)
}

## * model.matrix_regularize
## remove un-identifiable columns from the design matrix 
model.matrix_regularize <- function(formula, data){

    ## ** identify if there is an identifiability problem
    X <- stats::model.matrix(formula, data)
    X.qr <- qr(X)
    if(X.qr$rank==NCOL(X.qr$qr)){
        return(X)
    }else if(LMMstar.options()$drop.X==FALSE){
        stop("The design matrix does not have full rank according to the QR decomposition. \n")
    }
    
    ## ** prepare
    tt <- stats::delete.response(stats::terms(formula))
    tt.factors <- attr(tt,"factors")
    tt.term.labels <- attr(tt,"term.labels")
    tt.order <- attr(tt,"order")

    name.var <- all.vars(tt)
    test.factor <- sapply(name.var, function(iVar){is.character(data[[iVar]]) || is.factor(data[[iVar]])})
    var.factor <- name.var[test.factor]
    var.numeric <- name.var[!test.factor]
    name.factor <- name.var[sapply(name.var,function(iVar){is.factor(data[[iVar]])})]

    ## ** test 1: remove factor with empty level
    if(length(name.factor)>0){
        test.factor <- sapply(name.factor, function(iVar){all(levels(data[[iVar]]) %in% data[[iVar]])})
        if(any(test.factor==FALSE)){
            warning("Factor variable(s) with empty level: \"",paste(name.factor[test.factor==FALSE], collapse = "\" \""),"\"\n ",
                    "The empty level(s) will be remove internally. Consider applying droplevels to avoid this warning. \n")
            factor.drop <- name.factor[test.factor==FALSE]
            for(iFactor in factor.drop){
                data[[iFactor]] <- droplevels(data[[iFactor]])
            }
        }
    }

    ## ** create "naive" design matrix
    X <- .augmodel.matrix(tt,data)
    X.names <- colnames(X)
    X.assign <- tt.term.labels[attr(X,"assign")]
    X.Mlevel <- attr(X,"M.level")
    X.lslevel <- attr(X,"ls.level")
    X.reference <- attr(X,"reference")
    ## M.interactionName <- do.call(cbind,lapply(name.var, function(iVar){if(is.logical(X.Mlevel[,iVar])){return(c(NA,iVar)[X.Mlevel[,iVar]+1])}else{return(paste0(iVar,X.Mlevel[,iVar]))}}))
    ## vec.interactionName <- apply(M.interactionName, MARGIN = 1, FUN = function(iRow){paste(na.omit(iRow), collapse = ":")})

    ## ** test 2: form interaction and identify columns of the design matrix that are constant
    ls.rmX <- stats::setNames(lapply(which(tt.order>1), function(iInteraction){ ## iInteraction <- 4
        ## variables involved in the interactions
        iVar <- names(which(tt.factors[,iInteraction]>0)) 
        ## identify coefficient relative to this interaction 
        iNoVar <- setdiff(name.var,iVar)
        if(length(iNoVar)>0){ ## remove coefficients relative to other covariates
            iMlevel <- X.Mlevel[apply(X.Mlevel[,iNoVar,drop=FALSE],1,function(iParam){all(iParam==X.reference[iNoVar])}),,drop=FALSE]
        }else{
            iMlevel <- X.Mlevel
        }

        test.constant <- sapply(1:NROW(iMlevel), FUN = function(iParam){ ## iParam <- 1 
            iNumeric <- intersect(iVar, var.numeric)
            
            iValue <- do.call(cbind,lapply(var.factor, function(iFactor){data[[iFactor]]==iMlevel[iParam,iFactor]}))
            if(length(iNumeric)>0){
                iLs.ValueNum <- lapply(iNumeric, function(iN){if(iMlevel[iParam,iN]){data[[iN]]*iMlevel[iParam,iN]}else{rep(1,NROW(data))}})
                iValue <- cbind(iValue,do.call(cbind,iLs.ValueNum))
            }
            return(length(unique(apply(iValue, MARGIN = 1, prod)))<2)
        })
        return(rownames(iMlevel)[which(test.constant)])
    }), tt.term.labels[which(tt.order>1)])

    rmX <- unlist(ls.rmX)
    if(length(rmX)>0){
        warning("Constant values in the design matrix for interactions \"",paste(names(rmX), collapse = "\" \""),"\"\n ",
                "Coefficients \"",paste(rmX, collapse = "\" \""),"\" will be removed from the design matrix. \n",
                "Consider defining manually the interaction, e.g. via droplevels(interaction(.,.)) to avoid this warning. \n")
        X.old <- X
        test.keep <- colnames(X.old) %in% setdiff(colnames(X.old),rmX)
        X <- X.old[,test.keep]
        attr(X,"assign") <- attr(X.old,"assign")[test.keep] ## as.numeric(as.factor(attr(X.old,"assign")[test.keep])) - "(Intercept)" %in% colnames(X)
        attr(X,"contr.treatment") <- attr(X.old,"contr.treatment")
    }else{
        attr(X,"term.labels") <- NULL
        attr(X,"order") <- NULL
        attr(X,"ls.level") <- NULL
        attr(X,"M.level") <- NULL
        attr(X,"reference.level") <- NULL
    }
    
    ## ** export
    return(X)
}

## * .aumodel.matrix_match
## add more information about which variable contribute to which column in the design matrix
.augmodel.matrix <- function(formula, data){

    ## ** normalize formula
    formula.terms <- stats::delete.response(stats::terms(formula))
    var2term <- attr(formula.terms, "factors")
    
    ## ** create design matrix
    X <- stats::model.matrix(formula,data)
    p <- NCOL(X)
    
    ## ** identify factors and reference level
    all.variable <- all.vars(formula.terms)
    contrast.variable <- stats::setNames(lapply(all.variable, function(iVar){
        if(is.character(data[[iVar]])){
            out <- stats::contrasts(as.factor(data[[iVar]]))
            if(any(colSums(abs(out)>1e-12)>1)){
                stop("Cannot handle contrasts involving simultaneously several levels. \n")
            }
        }else if(is.factor(data[[iVar]])){
            out <- stats::contrasts(data[[iVar]])
            if(any(colSums(abs(out)>1e-12)>1)){
                stop("Cannot handle contrasts involving simultaneously several levels. \n")
            }
        }else{
            out <- NULL
        }
        return(out)
    }),all.variable)

    ls.reference <- stats::setNames(lapply(contrast.variable, function(iRef){
        if(!is.null(iRef)){
            return(rownames(iRef)[rowSums(abs(iRef)>0)==0])
        }else{
            return(NULL)
        }
    }),all.variable)
    reference <- unlist(ls.reference[!sapply(ls.reference,is.null)])
    reference2 <- data.frame(matrix(FALSE, nrow = 1, ncol = length(rownames(var2term)), dimnames = list(NULL, rownames(var2term))))
    reference2[,names(reference)] <- reference

    ## ** addtional informations
    X.names <- colnames(X)
    X.assign <- attr(X,"assign")
    X.order <- c(0,attr(formula.terms,"order"))[X.assign+1]
    X.term <- c("(Intercept)",attr(formula.terms,"term.labels"))[X.assign+1]
    X.level <- stats::setNames(vector(mode = "list", length = p), X.names)
    X.level2 <- stats::setNames(lapply(X.names, function(iName){reference2}), X.names)

    ## ** loop for each element of the design matrix and identify the right level
    for(iCol in 1:p){ ## iCol <- 5
        
        if(X.order[iCol]==0){
            ## reference level for intercept
            X.level[[iCol]] <- data.frame(reference)
            rownames(X.level[[iCol]]) <- NULL
        }else if(X.order[iCol]>0){
            ## variables involved
            iVar <- rownames(var2term)[var2term[,X.term[iCol]]==1]
            ## contrast for all factor/character variables involved
            iContrast <- contrast.variable[iVar]
            iContrast <- iContrast[!sapply(iContrast,is.null)]
            
            if(!is.null(iContrast) && length(iContrast)>0){
                ## re-create all possible names and identify the one matching the column name 
                iLs.factor <- stats::setNames(lapply(iVar,function(iName){rownames(iContrast[[iName]])}), iVar)
                iLs.name <- stats::setNames(lapply(iVar,function(iName){paste0(iName,iLs.factor[[iName]])}), iVar)
                iIndex  <- which(X.names[iCol] == interaction(iLs.name, sep=":"))
                ## deduce the factor variable
                iLs.index <- lapply(iVar,function(iName){if(is.null(iContrast[[iName]])){0}else{1:NROW(iContrast[[iName]])}})
                iIndex2 <- stats::setNames(as.numeric(unlist(strsplit(as.character(interaction(iLs.index, sep=":")[iIndex]), split = ":", fixed = TRUE))), iVar)
                X.level[[iCol]] <- as.data.frame(stats::setNames(lapply(iVar, function(iName){if(is.null(iLs.factor[[iName]])){as.numeric(NA)}else{iLs.factor[[iName]][iIndex2[[iName]]]}}), iVar))                
                X.level2[[iCol]][,names(X.level[[iCol]])] <- X.level[[iCol]]
                if(any(is.na(X.level[[iCol]]))){
                    X.level2[[iCol]][,is.na(X.level[[iCol]])] <- TRUE
                }
            }else{
                X.level[[iCol]] <- data.frame(matrix(as.numeric(NA), nrow = 1, ncol = length(iVar), dimnames = list(NULL, iVar)))
                X.level2[[iCol]][,iVar] <- TRUE
            }
        }
    }
    ## ** export
    attr(X,"term.labels") <- X.term
    attr(X,"order") <- X.order
    attr(X,"ls.level") <- X.level
    attr(X,"M.level") <- do.call(rbind,X.level2)
    attr(X,"reference.level") <- reference2
    return(X)
}


##----------------------------------------------------------------------
### model.matrix.R ends here
