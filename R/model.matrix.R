### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: sep 30 2021 (19:26) 
##           By: Brice Ozenne
##     Update #: 1247
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
        design <- list(mean = NULL,
                       vcov = NULL)

        ## prepare data
        if("variance" %in% effects){
            ff.allvars <- c(all.vars(object$formula$mean.design), all.vars(object$formula$var))
        }else{
            ff.allvars <- all.vars(object$formula$mean.design)
        }
        if(any(ff.allvars %in% names(data) == FALSE)){
            stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
        }
        ff.factor <- intersect(ff.allvars, names(object$xfactor))
        if(length(ff.factor)>0){
            for(iVar in ff.factor){ ## iVar <- ff.factor[2]
                if(any(data[[iVar]] %in% object$xfactor[[iVar]] == FALSE)){
                    Wf <- unique(data[[iVar]][data[[iVar]] %in% object$xfactor[[iVar]] == FALSE])
                    stop("Unknown factor(s) \"",paste0(Wf,collapse="\" \""),"\" for variable \"",iVar,"\".\n",
                         "Valid factors: \"",paste0(object$xfactor[[iVar]], collapse="\" \""),"\".\n")
                }
                data[[iVar]] <- factor(data[[iVar]], levels = object$xfactor[[iVar]])
            }
        }

        ## mean
        if("mean" %in% effects){
            design$mean <- .mean.matrix.lmm(formula = object$formula$mean.design, colnames = colnames(object$design$mean), data = data, 
                                            var.strata = object$strata$var, U.strata = object$strata$levels)
        }
    
        ## variance
        if("variance" %in% effects){
            indexData <- .extractIndexData(data = data, structure = object$design$vcov)
            design$vcov <- .vcov.matrix.lmm(structure = object$design$vcov, data = data,
                                           strata.var = indexData$strata.var, U.strata = indexData$U.strata,
                                           time.var = indexData$time.var, U.time = indexData$U.time,
                                           cluster.var = indexData$cluster.var, U.cluster = indexData$U.cluster, index.cluster = indexData$index.cluster,
                                           index.clusterTime = indexData$index.clusterTime, order.clusterTime = indexData$order.clusterTime)
            design$vcov$X$cor <- object$design$vcov$X$cor
            design$vcov$X$var <- object$design$vcov$X$var
        }
    }else{
        design <- object$design
    }

    ## ** update design matrix with new datase
    if(type.object == "lmm"){
        if("mean" %in% effects && "variance" %in% effects){
            return(list(mean = design$mean,
                        vcov = design$vcov))
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


## * .mean.matrix.lmm
.mean.matrix.lmm <- function(formula, colnames, data,
                             var.strata, U.strata){

    ## ** formula
    formula.model.matrix <- stats::as.formula(stats::delete.response(stats::terms(formula)))
    stratify.X <- length(setdiff(all.vars(formula),all.vars(formula.model.matrix)))>0

    ## ** design matrix
    if(stratify.X){
        ## *** normalize data
        ## index
        if("XXindexXX" %in% names(data) == FALSE){
            data$XXindexXX <- 1:NROW(data)
        }
        ## strata
        if(stratify.X && identical(var.strata, "XXstrata.indexXX") && var.strata %in% names(data) == FALSE){
            data$XXstrata.indexXX <- factor(1, levels = 1, labels = U.strata)
        }else{
            data[[var.strata]] <- factor(data[[var.strata]], levels = U.strata)
        }
        n.strata <- length(U.strata)

        ## *** generate design matrix for each strata
        ls.X.mean <- lapply(U.strata, function(iS){ ## iS <- U.strata[1]
            if(is.null(colnames)){
                iX <- .model.matrix_regularize(formula, data[data[[var.strata]]==iS,,drop=FALSE])
            }else{
                iX <- stats::model.matrix(formula, stats::model.frame(formula, data[data[[var.strata]]==iS,,drop=FALSE], na.action = stats::na.pass))
            }
            colnames(iX) <- paste0(colnames(iX),":",iS)
            attr(iX,"index") <- data[data[[var.strata]]==iS,"XXindexXX"]
            return(iX)
        })
        ## *** assemble
        X.mean <- as.matrix(Matrix::bdiag(ls.X.mean))[order(unlist(lapply(ls.X.mean, attr, "index"))),,drop=FALSE]
        colnames(X.mean) <- unlist(lapply(ls.X.mean,colnames))
        if(!is.null(colnames)){
            X.mean <- X.mean[,colnames,drop=FALSE]
        }
        attr(X.mean, "assign") <- as.vector(do.call(cbind,lapply(ls.X.mean,attr,"assign")))
        if(is.null(colnames)){
            strata.mu <- unlist(lapply(1:n.strata, function(iStrata){stats::setNames(rep(iStrata, NCOL(ls.X.mean[[iStrata]])),colnames(ls.X.mean[[iStrata]]))}))
        }else{
            attr.assign <- attr(X.mean, "assign")[match(colnames(X.mean),colnames)]
            X.mean <- X.mean[,colnames,drop=FALSE]
            attr(X.mean, "assign") <- attr.assign
        }
    }else{
        if(is.null(colnames)){
            X.mean <- .model.matrix_regularize(formula, data)
            strata.mu <- stats::setNames(rep(1,NCOL(X.mean)), colnames(X.mean))
        }else{
            X.mean <-  stats::model.matrix(formula, stats::model.frame(formula, data, na.action = stats::na.pass))[,colnames,drop=FALSE]
        }
    }

    ## ** export
    if(is.null(colnames)){
        attr(X.mean,"strata.mu") <- strata.mu
    }
    return(X.mean)

}

## * .vcov.matrix.lmm
.vcov.matrix.lmm <- function(structure, data, 
                            strata.var, U.strata,
                            time.var, U.time,
                            cluster.var, U.cluster, index.cluster, index.clusterTime, order.clusterTime){

    ## ** normalize data
    ## strata
    if(identical(strata.var, "XXstrata.indexXX") && strata.var %in% names(data) == FALSE){
        data$XXstrata.indexXX <- factor(1, levels = 1, labels = U.strata)
    }
    n.strata <- length(U.strata)
    
    ## time
    n.time <- length(U.time)

    ## formula
    formula.var <- structure$formula$var
    formula.cor <- structure$formula$cor

    ## ** check compatibility structre and data
    if(is.na(structure$name$strata)){
        structure$name$strata <- strata.var
    }else if(!identical(structure$name$strata,strata.var)){
        stop("The name of the strata variable in the residual variance-covariance structure does not match the actual strata variable. \n")
    }
    if(is.na(structure$name$time)){
        structure$name$time <- time.var
    }else if(!identical(structure$name$time,time.var)){
        stop("The name of the time variable in the residual variance-covariance structure does not match the actual time variable. \n")
    }
    if(is.na(structure$name$cluster)){
        structure$name$cluster <- cluster.var
    }else if(!identical(structure$name$cluster,cluster.var)){
        stop("The name of the cluster variable in the residual variance-covariance structure does not match the actual cluster variable. \n")
    }
    ## ** design matrix
    if(is.null(structure$param)){ ## structure
        out <- list(var = NULL, cor = NULL)
        out$var <- .colnameOrder(.model.matrix_regularize(formula.var, data = data, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)
        if(!is.null(formula.cor) && n.time>1 && any(sapply(order.clusterTime,length)>1)){  ## at least one individual with more than timepoint
            out$cor <- .colnameOrder(.model.matrix_regularize(formula.cor, data = data, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)
        }
    }else{ ## newdata
        out <- structure
        out$U.cluster <- U.cluster
        X.var <- model.matrix(formula.var, data = data)[,attr(structure$X$var,"original.colnames")]
        index.varPattern <- .buildVarPatterns(X.var = X.var, data = data, Upattern = structure$X$Upattern,
                                              U.cluster, index.cluster = index.cluster, index.clusterTime = index.clusterTime, order.clusterTime = order.clusterTime,
                                              U.strata = U.strata, strata.var = strata.var)

        if(!is.null(structure$X$cor)){
            X.cor <- model.matrix(formula.cor, data = data)[,attr(structure$X$cor,"original.colnames")]
            index.corPattern <- .buildCorPatterns(X.cor = X.cor, Upattern = structure$X$Upattern,
                                                  U.cluster, index.cluster = index.cluster, order.clusterTime = order.clusterTime)

            index.pattern <- match(paste(index.varPattern,index.corPattern,sep="."),paste(structure$X$Upattern$var,structure$X$Upattern$cor,sep="."))
            out$X$pattern.cluster <- stats::setNames(structure$X$Upattern$name[index.pattern],U.cluster)
        }else{
            out$X$pattern.cluster <- stats::setNames(structure$X$Upattern$name[match(index.varPattern,structure$X$Upattern$var)],U.cluster)
        }
        out$X$cluster.pattern <- stats::setNames(lapply(out$X$Upattern$name,function(iN){unname(which(iN==out$X$pattern.cluster))}), out$X$Upattern$name)
        out$X$Upattern$ncluster <- sapply(out$X$Upattern$name,function(iN){sum(iN==out$X$pattern.cluster)})
    }

    ## ** export
    return(out)
}

## * .model.matrix.lmm
.model.matrix.lmm <- function(formula.mean, structure,
                              data, var.outcome,
                              var.strata, U.strata,
                              var.time, U.time,
                              var.cluster,
                              precompute.moments){

    ## ** mean
    X.mean <- .mean.matrix.lmm(formula = formula.mean, colnames = NULL, data = data,
                               var.strata = var.strata, U.strata = U.strata)
    strata.mu <- attr(X.mean,"strata.mu")
    attr(X.mean,"strata.mu") <- NULL

    ## ** variance (update structure)
    structure <- .skeleton(structure = structure, data = data)

    ## ** cluster
    U.cluster <- sort(unique(data[[var.cluster]]))
    if(is.factor(U.cluster)){
        U.cluster <- as.character(U.cluster)
    }
    n.cluster <- length(U.cluster)
    index.cluster <- match(data[[var.cluster]], U.cluster) ## ‘match’ returns a vector of the positions of (first) matches of its first argument in its second.
    index.time <- as.numeric(data[[var.time]])
    attr(index.cluster,"sorted") <- lapply(1:n.cluster, function(iId){
        iIndex <- which(index.cluster==iId)
        return(iIndex[order(index.time[iIndex])]) ## re-order observations according to the variance-covariance matrix
    })

    ## ** prepare calculation of the score
    if(precompute.moments){
        precompute.XX <-  .precomputeXX(X = X.mean, pattern = structure$X$Upattern$name,
                                        pattern.time = structure$X$Upattern$time, pattern.cluster = structure$X$cluster.pattern, index.cluster = attr(index.cluster,"sorted"))
        precompute.XY <-  .precomputeXR(X = precompute.XX$Xpattern, residuals = cbind(data[[var.outcome]]), pattern = structure$X$Upattern$name,
                                        pattern.time = structure$X$Upattern$time, pattern.cluster = structure$X$cluster.pattern, index.cluster = attr(index.cluster,"sorted"))
    }else{
        precompute.XX <- NULL
        precompute.XY <- NULL
    }

    ## ** pairs
    Upattern.strata <- stats::setNames(structure$X$Upattern$strata,structure$X$Upattern$name)
    pair.meanvarcoef <- stats::setNames(lapply(structure$X$Upattern$name, function(iPattern){ ## iPattern <- structure$X$Upattern$name[1]
        iParamMu <- names(strata.mu)[U.strata[strata.mu]==Upattern.strata[iPattern]]
        iParamVar <- structure$X$Upattern$param[[iPattern]]
        iOut <- unname(t(expand.grid(iParamMu, iParamVar)))
        return(iOut)
    }), structure$X$Upattern$name)

    ## ** param
    skeleton.param <- list(mu = colnames(X.mean),
                           sigma = structure$param[structure$param$type=="sigma","name"],
                           k = structure$param[structure$param$type=="k","name"],
                           rho = structure$param[structure$param$type=="rho","name"],
                           strata.mu = strata.mu,
                           strata.sigma = structure$param[structure$param$type=="sigma","strata"],
                           strata.k = structure$param[structure$param$type=="k","strata"],
                           strata.rho = structure$param[structure$param$type=="rho","strata"],
                           time.sigma = unlist(structure$param[structure$param$type=="sigma","time"]),
                           time.k = unlist(structure$param[structure$param$type=="k","time"]),
                           time.rho = NULL,
                           pair.meanvarcoef = pair.meanvarcoef)
    if(length(skeleton.param$rho)>0){
        skeleton.param$time.rho <- base::t(do.call(rbind,structure$param[structure$param$type=="rho","time"]))
    }
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

## * helpers
## ** .model.matrix_regularize
## remove un-identifiable columns from the design matrix 
.model.matrix_regularize <- function(formula, data, augmodel = FALSE){

    ## ## ** test 0: remove variable(s) with single level in the formula
    test.1value <- sapply(all.vars(formula),function(iVar){length(unique(data[[iVar]]))})
    if(any(test.1value==1)){
        ## no warning because this is normal behavior when stratifying
        formula <- stats::update(formula, stats::as.formula(paste0("~.-",paste(names(test.1value)[test.1value==1],collapse="-"))))
    }
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

## ** .aumodel.matrix_match
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
                X.level[[iCol]] <- as.data.frame(stats::setNames(lapply(iVar, function(iName){if(is.null(iValue.factor[[iName]])){as.numeric(NA)}else{as.character(iValue.factor[1,iName])}}), iVar))                
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

    ## ** add attributes
    attr(X,"term.labels") <- X.term
    attr(X,"order") <- X.order
    attr(X,"ls.level") <- X.level
    attr(X,"M.level") <- do.call(rbind,X.level2)
    attr(X,"reference.level") <- reference2

    ## ** export
    return(X)
}

## ** .extractIndexData
.extractIndexData <- function(data, structure){
    
    ## *** find variable names
    time.var <- structure$name$time
    cluster.var <- structure$name$cluster
    strata.var <- structure$name$strata
    all.var <- stats::na.omit(c(time.var,cluster.var,strata.var,unlist(lapply(structure$formula,all.vars))))
    if(length(all.var)>0){
        for(iVar in all.var){
            if(is.factor(data[[iVar]]) && is.null(structure$X)){
                ## if the design matrix for the covariance pattern has been built
                ## then we may ask prediction for subsets so we don't want to drop levels
                data[[iVar]] <- droplevels(data[[iVar]])
            }else{
                data[[iVar]] <- as.factor(data[[iVar]])
            }
        }
    }

    ## *** find unique levels of each variable
    if(!is.na(cluster.var)){
        if(is.factor(data[[cluster.var]])){
            U.cluster <- levels(data[[cluster.var]])
        }else{
            U.cluster <- as.character(unique(data[[cluster.var]]))
        }
    }else{
        U.cluster <- paste0("id",1:NROW(data))
    }
    
    if(!is.na(time.var)){
        if(is.factor(data[[time.var]])){
            U.time <- levels(data[[time.var]])
        }else{
            U.time <- as.character(unique(data[[time.var]]))
        }
    }else if(!is.na(cluster.var)){
        U.time <- paste0("T",1:max(table(data[[cluster.var]])))
    }else{
        U.time <- "T1"
    }

    if(!is.na(strata.var)){
        if(is.factor(data[[strata.var]])){
            U.strata <- levels(data[[strata.var]])
        }else{
            U.strata <- as.character(unique(data[[strata.var]]))
        }
    }else{
        U.strata <- "S1"
    }

    ## *** find position of each cluster
    if(!is.na(cluster.var)){
        index.cluster <- tapply(1:NROW(data),data[[cluster.var]],function(iI){iI})
    }else{
        index.cluster <- stats::setNames(as.list(1:NROW(data)),U.cluster)
    }

    ## *** find time corresponding to each cluster
    if(is.na(cluster.var) && is.na(time.var)){
        index.clusterTime <- stats::setNames(as.list(rep(1,NROW(data))),U.cluster)
    }else if(!is.na(cluster.var) && !is.na(time.var)){
        index.clusterTime <- tapply(data[[time.var]],data[[cluster.var]],function(iT){as.numeric(factor(iT,levels = U.time))})
    }else if(!is.na(time.var)){
        index.clusterTime <- stats::setNames(as.list(as.numeric(factor(data[[time.var]], levels = U.time))),U.cluster)
    }else if(!is.na(cluster.var)){
        index.clusterTime <- tapply(data[[cluster.var]],data[[cluster.var]],function(iT){1:length(iT)})
    }
    order.clusterTime <- lapply(index.clusterTime, order)


    ## *** export
    out <- list(data = data,
                time.var = time.var,
                U.time = U.time,
                cluster.var = cluster.var,
                U.cluster = U.cluster,
                strata.var = strata.var,
                U.strata = U.strata,
                index.cluster = index.cluster,
                index.clusterTime = index.clusterTime,
                order.clusterTime = order.clusterTime
                )

    return(out)
}
## ** .buildVarPatterns
## convert design matrix into character string with associated strata and time
.buildVarPatterns <- function(X.var, data, Upattern,
                              U.cluster, index.cluster, index.clusterTime, order.clusterTime,
                              U.strata, strata.var){
    
    ## *** combine all covariate values at the observation levels and convert to positive integer
    level.var <- interaction(as.data.frame(X.var),drop=TRUE)
    if(!is.null(Upattern)){ ## restaure all levels
        level.var <- factor(level.var, levels = attr(Upattern,"levels.UpatternVar"))
    }
    Ulevel.var <- levels(level.var)
    numlevel.var <- as.numeric(level.var)

    ## *** structure the covariate pattern by cluster and add strata and time indexes
    patternVar.cluster <- stats::setNames(lapply(U.cluster, function(iC){ ## iC <- U.cluster[[1]]
        iIndex.obs <- index.cluster[[iC]][order.clusterTime[[iC]]]
        iLevel <- paste(numlevel.var[iIndex.obs], collapse=".")
        attr(iLevel,"index.level") <- numlevel.var[iIndex.obs]
        attr(iLevel,"index.alltime") <- index.clusterTime[[iC]][order.clusterTime[[iC]]]
        if(!is.na(strata.var)){
            attr(iLevel,"index.strata") <- which(unique(data[index.cluster[[iC]],strata.var])==U.strata)
        }else{
            attr(iLevel,"index.strata") <- 1
        }
        return(iLevel)
    }), U.cluster)
    attr(patternVar.cluster,"levels") <- Ulevel.var
    
    ## *** export
    if(is.null(Upattern)){
        return(patternVar.cluster)
    }else{
        return(.matchPatterns(cluster.pattern = patternVar.cluster,
                              Upattern = attr(Upattern,paste0("name.UpatternVar")), Upattern.levels = attr(Upattern,paste0("levels.UpatternVar")))
               )
    }
}

## ** .buildCorPatterns
## convert design matrix into character string with associated strata and time
.buildCorPatterns <- function(X.cor, Upattern,
                              U.cluster, index.cluster, order.clusterTime){

    ## *** deal with the case of no correlation
    if(is.null(X.cor)){
        if(is.null(Upattern)){
            patternCor.cluster <- stats::setNames(rep("1", length(U.cluster)),U.cluster)
            attr(patternCor.cluster,"levels") <- "1"
        }else{
            return(rep(1,length(U.cluster)))
        }
    }
    UpatternCor.cluster <- "1"
    patternCor.cluster <- stats::setNames(rep(1, length = length(U.cluster)), U.cluster)

    ## *** combine all covariate values at the observation levels and convert to positive integer
    level.cor <- interaction(as.data.frame(X.cor),drop=TRUE)
    if(!is.null(Upattern)){
        level.cor <- factor(level.cor, levels = attr(Upattern,"levels.UpatternCor"))
    }
    Ulevel.cor <- levels(level.cor)
    numlevel.cor <- as.numeric(level.cor)
    
    ## *** structure the covariate pattern by cluster
    patternCor.cluster <- stats::setNames(sapply(U.cluster, function(iC){
        paste(numlevel.cor[index.cluster[[iC]][order.clusterTime[[iC]]]],collapse=".")
    }), U.cluster)
    attr(patternCor.cluster,"levels") <- Ulevel.cor

    ## *** export
    if(is.null(Upattern)){
        return(patternCor.cluster)
    }else{
        return(.matchPatterns(cluster.pattern = patternCor.cluster,
                              Upattern = attr(Upattern,paste0("name.UpatternCor")), Upattern.levels = attr(Upattern,paste0("levels.UpatternCor")))
               )
    }
}

## ** .matchPatterns
.matchPatterns <- function(cluster.pattern, Upattern, Upattern.levels){

    vec.cluster.pattern <- unlist(cluster.pattern)
    
    if(all(vec.cluster.pattern %in% Upattern)){
        out <- match(vec.cluster.pattern, Upattern)
    }else{ ## some requested patterns may be a subset of the available ones. This tries to handle that.
        missingPattern <- setdiff(vec.cluster.pattern,Upattern)
        ls.UpatternLevels <- stats::setNames(lapply(strsplit(Upattern,split=".",fixed = TRUE),
                                             function(iPattern){
                                                 Upattern.levels[as.numeric(iPattern)]
                                             }),Upattern)
        substitutePattern <- sapply(missingPattern, function(iC){
            iLevel <- Upattern.levels[as.numeric(strsplit(iC,split=".",fixed = TRUE)[[1]])]
            names(ls.UpatternLevels)[sapply(ls.UpatternLevels, function(iPattern){all(iLevel %in% iPattern)})][1]
        })
        out <- as.numeric(factor(c(substitutePattern,stats::setNames(Upattern,Upattern))[vec.cluster.pattern], levels = Upattern))
    }
    return(stats::setNames(out, names(cluster.pattern)))
}


## ** .unorderedPairs
## adapted from RecordLinkage package
## form all combinations
## .unorderedPairs(1:5)
.unorderedPairs <- function(x, distinct = FALSE){
    n <- length(x)
    ls <- lapply(1:n, function(k){ rbind(x[k], x[k:n])})
    out <- do.call(cbind,ls)##array(unlist(ls), dim = c(2, n * (n + 1)/2))
    if(distinct){
        out <- out[,apply(out,2,function(iCol){all(!duplicated(iCol))}),drop=FALSE]
    }
    return(out)
}



##----------------------------------------------------------------------
### model.matrix.R ends here
