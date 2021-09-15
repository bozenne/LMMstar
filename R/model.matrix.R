### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: sep 15 2021 (18:48) 
##           By: Brice Ozenne
##     Update #: 966
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
model.matrix.lmm <- function(object, data = NULL, effects = "mean", type.object = "lmm", ...){

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

    ## ** update design matrix with new datase
    if(type.object == "lmm"){
        if("mean" %in% effects && "variance" %in% effects){
            return(list(mean = design$mean,
                        variance = design$vcov))
        }else if("mean" %in% effects){
            return(design$mean)
        }else if("variance" %in% effects){
            return(design$vcov)
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
.model.matrix.lmm <- function(formula.mean, structure,
                              data, var.outcome,
                              var.strata, U.strata,
                              var.time, U.time,
                              var.cluster,
                              precompute.moments){
    n.obs <- NROW(data)
    n.strata <- length(U.strata)
    n.time <- length(U.time)
    
    ## ** normalize data
    ## index
    if("XXindexXX" %in% names(data) == FALSE){
        data$XXindexXX <- 1:NROW(data)
    }
    if(identical(var.strata, "XXstrata.indexXX") && var.strata %in% names(data) == FALSE){
        data$XXstrata.indexXX <- 1
    }else{
        data[[var.strata]] <- factor(data[[var.strata]], levels = U.strata)
    }
    index.strata <- tapply(data[[var.strata]],data[[var.cluster]],unique)

    ## time
    data[[var.time]] <- factor(data[[var.time]], levels = U.time)
    index.time <- as.numeric(data[[var.time]])

    ## cluster
    U.cluster <- sort(unique(data[[var.cluster]]))
    if(is.factor(U.cluster)){
        U.cluster <- as.character(U.cluster)
    }
    n.cluster <- length(U.cluster)
    index.cluster <- match(data[[var.cluster]], U.cluster) ## ‘match’ returns a vector of the positions of (first) matches of its first argument in its second.
    attr(index.cluster,"sorted") <- lapply(1:n.cluster, function(iId){
        iIndex <- which(index.cluster==iId)
        return(iIndex[order(index.time[iIndex])]) ## re-order observations according to the variance-covariance matrix
    })

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


    ## ** variance
    if(is.na(structure$name$strata)){
        structure$name$strata <- var.strata
    }else if(!identical(structure$name$strata,var.strata)){
        stop("The name of the strata variable in the residual variance-covariance structure does not match the actual strata variable. \n")
    }
    if(is.na(structure$name$time)){
        structure$name$time <- var.time
    }else if(!identical(structure$name$time,var.time)){
        stop("The name of the time variable in the residual variance-covariance structure does not match the actual time variable. \n")
    }
    if(is.na(structure$name$cluster)){
        structure$name$cluster <- var.cluster
    }else if(!identical(structure$name$cluster,var.cluster)){
        stop("The name of the cluster variable in the residual variance-covariance structure does not match the actual cluster variable. \n")
    }
    structure <- skeletonStructure(structure, data = data)
    
    ## ** prepare calculation of the score
    if(precompute.moments){
        precompute.XX <-  .precomputeXX(X = X.mean, pattern = structure$X$Upattern$name,
                                        pattern.time = structure$X$Upattern$time, pattern.cluster = structure$X$Upattern$cluster, index.cluster = attr(index.cluster,"sorted"))
        precompute.XY <-  .precomputeXR(X = precompute.XX$Xpattern, residuals = cbind(data[[var.outcome]]), pattern = structure$X$Upattern$name,
                                        pattern.time = structure$X$Upattern$time, pattern.cluster = structure$X$Upattern$cluster, index.cluster = attr(index.cluster,"sorted"))
    }else{
        precompute.XX <- NULL
        precompute.XY <- NULL
    }
    
    ## ** pairs
    pair.meanvarcoef <- stats::setNames(lapply(structure$X$Upattern$name, function(iPattern){ ## iPattern <- structure$X$Upattern$name[1]
        iParamMu <- names(strata.mu)[U.strata[strata.mu]==structure$X$Upattern$strata[[iPattern]]]
        iParamVar <- structure$X$Upattern$param[[iPattern]]
        iOut <- unname(t(expand.grid(iParamMu, iParamVar)))
        return(iOut)
    }), structure$X$Upattern$name)
    
    pair.varcoef <- stats::setNames(lapply(structure$X$Upattern$name, function(iPattern){## ## iPattern <- structure$X$Upattern$name[1]
        iParam <- structure$X$Upattern$param[[iPattern]]

        iOut <- .unorderedPairs(iParam)
        attr(iOut, "key") <- matrix(NA, nrow = length(iParam), ncol = length(iParam), dimnames = list(iParam,iParam))
        for(iCol in 1:NCOL(iOut)){
            attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
            attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol]] <- iCol
        }
        return(iOut)
    }),structure$X$Upattern$name)

    ## ** param
    skeleton.param <- list(mu = colnames(X.mean),
                           sigma = structure$param[structure$param$type=="sigma","name"],
                           k = structure$param[structure$param$type=="sigma","k"],
                           rho = structure$param[structure$param$type=="rho","name"],
                           strata.mu = strata.mu,
                           strata.sigma = structure$param[structure$param$type=="sigma","strata"],
                           strata.k = structure$param[structure$param$type=="k","strata"],
                           strata.rho = structure$param[structure$param$type=="rho","strata"],
                           time.sigma = structure$param[structure$param$type=="sigma","time"],
                           time.k = structure$param[structure$param$type=="k","time"],
                           time.rho = structure$param[structure$param$type=="rho","time"],
                           pair.meanvarcoef = pair.meanvarcoef,
                           pair.varcoef = pair.varcoef)
    name.param <- c(skeleton.param$mu,skeleton.param$sigma,skeleton.param$k,skeleton.param$rho)
    skeleton.param$type <- stats::setNames(c(rep("mu",length(skeleton.param$mu)),
                                             rep("sigma",length(skeleton.param$sigma)),
                                             rep("k",length(skeleton.param$k)),
                                             rep("rho",length(skeleton.param$rho))),
                                           name.param)
    skeleton.param$strata <- stats::setNames(c(skeleton.param$strata.mu,
                                               skeleton.param$strata.sigma,
                                               skeleton.param$strata.k,
                                               skeleton.param$strata.rho), name.param)


    ## ** gather and export
    out <- list(mean = X.mean,
                vcov = structure,
                Y = data[[var.outcome]],
                precompute.XX = precompute.XX,
                precompute.XY = precompute.XY,
                index.cluster = index.cluster,
                index.time = index.time,
                cluster = list(n = n.cluster, levels = U.cluster, nobs = table(index.cluster)),
                param = skeleton.param
                )
    return(out)
}

## * model.matrix_regularize
## remove un-identifiable columns from the design matrix 
model.matrix_regularize <- function(formula, data, augmodel = FALSE){

    ## ** identify if there is an identifiability problem
    X <- stats::model.matrix(formula, data)
    X.qr <- qr(X)
    if(X.qr$rank==NCOL(X.qr$qr)){
        if(augmodel){
            return(.augmodel.matrix(stats::delete.response(stats::terms(formula)),data))
        }else{
            return(X)
        }
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
    ## M.interactionName <- do.call(cbind,lapply(name.var, function(iVar){if(is.logical(X.Mlevel[,iVar])){return(c(NA,iVar)[X.Mlevel[,iVar]+1])}else{return(paste0(iVar,X.Mlevel[,iVar]))}}))
    ## vec.interactionName <- apply(M.interactionName, MARGIN = 1, FUN = function(iRow){paste(na.omit(iRow), collapse = ":")})
    
    ## ** test 2: remove column(s) corresponding to the same level
    index.keep <- which(!duplicated(attr(X,"M.level")))
    Xsave <- X
    X <- Xsave[,index.keep,drop=FALSE]
    attrX <- attributes(Xsave)[setdiff(names(attributes(Xsave)),names(attributes(X)))]
    attrX$assign <- attrX$assign[index.keep]
    attrX$term.labels <- attrX$term.labels[index.keep]
    attrX$order <- attrX$order[index.keep]
    attrX$ls.levels <- attrX$ls.levels[index.keep]
    attrX$M.level <- attrX$M.level[index.keep,,drop=FALSE]
    attributes(X) <- c(attributes(X),attrX)
    
    X.Mlevel <- attr(X,"M.level")
    X.reference <- attr(X,"reference")

    ## ** test 3: form interaction and identify columns of the design matrix that are constant
    ls.rmX <- stats::setNames(lapply(which(tt.order>1), function(iInteraction){ ## iInteraction <- 2 
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
        warning("Constant values in the design matrix in interactions \"",paste(names(ls.rmX), collapse = "\" \""),"\"\n ",
                "Coefficients \"",paste(rmX, collapse = "\" \""),"\" will be removed from the design matrix. \n",
                "Consider defining manually the interaction, e.g. via droplevels(interaction(.,.)) to avoid this warning. \n")
        X.old <- X
        test.keep <- colnames(X.old) %in% setdiff(colnames(X.old),rmX)
        X <- X.old[,test.keep]
        attr(X,"assign") <- attr(X.old,"assign")[test.keep] ## as.numeric(as.factor(attr(X.old,"assign")[test.keep])) - "(Intercept)" %in% colnames(X)
        attr(X,"contr.treatment") <- attr(X.old,"contr.treatment")
        if(augmodel){
            attr(X,"term.labels") <- attr(X.old,"term.labels")[test.keep]
            attr(X,"order") <- attr(X.old,"order")[test.keep]
            attr(X,"M.level") <- attr(X.old,"M.level")[test.keep,,drop=FALSE]
        }
    }else if(!augmodel){
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
    for(iCol in 1:p){ ## iCol <- 9

        if(X.order[iCol]==0){
            ## reference level for intercept
            X.level[[iCol]] <- data.frame(reference)
            rownames(X.level[[iCol]]) <- NULL
        }else if(X.order[iCol]>0){
            ## variables involved
            iVar <- rownames(var2term)[var2term[,X.term[iCol]]>=1]
            ## contrast for all factor/character variables involved
            iContrast <- contrast.variable[iVar]
            iContrast <- iContrast[!sapply(iContrast,is.null)]

            if(!is.null(iContrast) && length(iContrast)>0){
                ## re-create all possible names and identify the one matching the column name 
                iLs.factor <- stats::setNames(lapply(iVar,function(iName){rownames(iContrast[[iName]])}), iVar)
                iLs.grid <- expand.grid(iLs.factor)
                iLs.allnameInteraction <- interaction(stats::setNames(lapply(iVar,function(iName){paste0(iName,iLs.grid[[iName]])}), iVar), sep = ":")
                iIndex  <- which(X.names[iCol] == iLs.allnameInteraction)
                ## deduce the factor variable
                iValue.factor <- iLs.grid[iIndex,,drop=FALSE]
                X.level[[iCol]] <- as.data.frame(stats::setNames(lapply(iVar, function(iName){if(is.null(iValue.factor[[iName]])){as.numeric(NA)}else{iValue.factor[1,iName]}}), iVar))                
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

## * .unorderedPairs
## adapted from RecordLinkage package
## form all combinations
.unorderedPairs <- function(x, distinct = FALSE){
    n <- length(x)
    ls <- lapply(1:n, function(k){ rbind(x[k], x[k:n])})
    out <- array(unlist(ls), dim = c(2, n * (n + 1)/2))
    if(distinct){
        out <- out[,apply(out,2,function(iCol){all(!duplicated(iCol))}),drop=FALSE]
    }
    return(out)
}

## * .precompute
## ** .precomputeXX
.precomputeXX <- function(X, pattern, pattern.time, pattern.cluster, index.cluster){
    p <- NCOL(X)
    n.pattern <- length(pattern)
    n.time <- lapply(pattern.time,length)

    out <- list(pattern = stats::setNames(lapply(pattern, function(iPattern){array(0, dim = c(n.time[[iPattern]],n.time[[iPattern]],p*(p+1)/2))}), pattern),
                key = matrix(as.numeric(NA),nrow=p,ncol=p,dimnames=list(colnames(X),colnames(X))),
                Xpattern = stats::setNames(vector(mode = "list", length = n.pattern),pattern))

    ## key
    out$key[lower.tri(out$key,diag = TRUE)] <- 1:sum(lower.tri(out$key,diag = TRUE))
    out$key[upper.tri(out$key)] <- t(out$key)[upper.tri(out$key)]

    ## fill matrix
    for(iPattern in pattern){ ## iPattern <- pattern[1]
        iTime <- length(pattern.time[[iPattern]])

        if(iTime==1){
            out$Xpattern[[iPattern]] <- do.call(rbind,lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){X[iIndex,,drop=FALSE]}))
            iX.summary <- crossprod(out$Xpattern[[iPattern]])
            ## out$key[lower.tri(out$key,diag = TRUE)]
            out$pattern[[iPattern]][1,1,] <- iX.summary[lower.tri(iX.summary, diag = TRUE)]
        }else{
            out$Xpattern[[iPattern]] <- array(unlist(lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){X[iIndex,,drop=FALSE]})),
                                              dim = c(iTime,NCOL(X),length(index.cluster[pattern.cluster[[iPattern]]])))

            for(iCol1 in 1:p){ ## iCol1 <- 1
                for(iCol2 in 1:iCol1){ ## iCol2 <- 2
                    ## for(iId in pattern.cluster[[iPattern]]){
                    ##     out$pattern[[iPattern]][,,out$key[iCol1,iCol2]] <- out$pattern[[iPattern]][,,out$key[iCol1,iCol2]] + tcrossprod(X[index.cluster[[iId]],iCol1,drop=FALSE],X[index.cluster[[iId]],iCol2,drop=FALSE])
                    ## }
                    out$pattern[[iPattern]][,,out$key[iCol1,iCol2]] <- tcrossprod(out$Xpattern[[iPattern]][,iCol1,],out$Xpattern[[iPattern]][,iCol2,])
                }
            }
        }
    }
    return(out)
}

## ** .precomputeXR
.precomputeXR <- function(X, residuals, pattern, pattern.time, pattern.cluster, index.cluster){
    p <- NCOL(X[[1]])
    name.mucoef <- colnames(X)
    n.pattern <- length(pattern)
    n.time <- lapply(pattern.time,length)

    out <- stats::setNames(lapply(pattern, function(iPattern){array(0, dim = c(n.time[[iPattern]], n.time[[iPattern]], ncol = p), dimnames = list(NULL,NULL,name.mucoef))}), pattern)

    for(iPattern in pattern){ ## iPattern <- pattern[1]

        iTime <- length(pattern.time[[iPattern]])
        iResiduals <- do.call(cbind, lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){residuals[iIndex,,drop=FALSE]}))
        iX <- X[[iPattern]]

        if(iTime == 1){
            out[[iPattern]][1,1,] <- iResiduals %*% iX
        }else{
            for(iCol in 1:p){ ## iCol1 <- 1
                ## for(iId in 1:length(pattern.cluster[[iPattern]])){ ## iId <- 1
                ##     out[[iPattern]][,,iCol] <- out[[iPattern]][,,iCol] + tcrossprod(iX[[iId]][,iCol,drop=FALSE],iResiduals[[iId]])
                ## }
                out[[iPattern]][,,iCol] <- tcrossprod(iX[,iCol,], iResiduals)
            }
        }
    }

    return(out)
}

## ** .precomputeRR
.precomputeRR <- function(residuals, pattern.time, pattern, pattern.cluster, index.cluster){

    n.pattern <- length(pattern)
    n.time <- stats::setNames(lapply(pattern.time,length), pattern)
    out <- stats::setNames(lapply(pattern, function(iPattern){matrix(0, nrow = n.time[[iPattern]], ncol = n.time[[iPattern]])}), pattern)

    for(iPattern in pattern){ ## iPattern <- pattern[1]
        
        ## for(iId in pattern.cluster[[iPattern]]){
        ##     out[[iPattern]] <- out[[iPattern]] + tcrossprod(residuals[index.cluster[[iId]],,drop=FALSE])
        ## }
        out[[iPattern]] <- tcrossprod(do.call(cbind,lapply(index.cluster[pattern.cluster[[iPattern]]], function(iIndex){residuals[iIndex,,drop=FALSE]})))
    }
    return(out)
}

##----------------------------------------------------------------------
### model.matrix.R ends here
