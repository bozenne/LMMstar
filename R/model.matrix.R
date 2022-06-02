### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: Jun  2 2022 (11:57) 
##           By: Brice Ozenne
##     Update #: 2259
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
model.matrix.lmm <- function(object, data = NULL, effects = "mean", simplifies = TRUE, ...){

    ## ** normalize user imput
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

        ## *** prepare output
        design <- list(mean = NULL,
                       vcov = NULL,
                       Y = NULL,
                       precompute.XX = NULL,
                       precompute.XY = NULL,
                       index.cluster = NULL,
                       index.clusterTime = NULL,
                       index.clusterStrata = NULL,
                       index.time = NULL,
                       time = object$design$time,
                       cluster = NULL,
                       param = object$design$param)
        cluster.var <- attr(object$cluster$var,"original")
        time.var <- attr(object$time$var,"original")
        strata.var <- attr(object$strata$var,"original")

        ## *** outcome
        if(object$outcome$var %in% names(data) && !simplifies){
            if(any(object$outcome$var %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing outcome variable \"",object$outcome$var,"\".\n")
            }
            design$Y <- data[[object$outcome$var]]
        }

        ## *** weights
        if(object$outcome$var %in% names(data) && !simplifies){
            if(!is.na(object$weight$var[1])){
                if(any(object$weight$var[1] %in% names(data) == FALSE)){
                    stop("Incorrect argument \'data\': missing weight variable \"",object$weight$var[1],"\".\n")
                }
                design$weight <- data.var[[object$weight$var[1]]]
            }
            if(!is.na(object$weight$var[2])){
                if(any(object$weight$var[2] %in% names(data) == FALSE)){
                    stop("Incorrect argument \'data\': missing weight variable \"",object$weight$var[2],"\".\n")
                }
                design$scale.Omega <- data.var[[object$weight$var[2]]]
            }
        }
        
        ## *** mean
        if("mean" %in% effects){
            data.mean <- data

            ## check dataset
            ff.allvars <- all.vars(object$formula$mean.design)
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) for the mean structure \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            ## convert to factor with the right levels
            ff.factor <- names(object$xfactor$mean)
            if(length(ff.factor)>0){
                for(iVar in ff.factor){ ## iVar <- ff.factor[1]
                    if(any(data[[iVar]] %in% object$xfactor$mean[[iVar]] == FALSE)){
                        Wf <- setdiff(unique(data[[iVar]]), iLevel)
                        stop("Unknown factor(s) \"",paste0(Wf,collapse="\" \""),"\" for variable \"",iVar,"\".\n",
                             "Valid factors: \"",paste0(object$xfactor$mean[[iVar]], collapse="\" \""),"\".\n")
                    }
                    data.mean[[iVar]] <- factor(data[[iVar]], levels = object$xfactor$mean[[iVar]])
                }
            }

            ## necessary to add the strata
            data.mean <- .prepareData(data.mean,
                                      var.cluster = NA,
                                      var.time = NA,
                                      var.strata = strata.var,
                                      missing.repetition = NULL,
                                      droplevels = TRUE)

            ## use stats::model.frame to handle spline
            design$mean  <- .mean.matrix.lmm(formula = object$formula$mean.design, colnames = colnames(object$design$mean),
                                             data = stats::model.frame(attr(object$design$mean,"terms"), data = data.mean , na.action = stats::na.pass), 
                                             stratify = (object$opt$name=="gls") && (object$strata$n>1), name.strata = object$strata$var, U.strata = object$strata$levels) ## only stratify mean if gls optimizer
        }

        ## *** variance-covariance
        if("variance" %in% effects){
            ## check dataset
            data.var <- data
            ff.factor <- unique(c(names(object$xfactor$var),names(object$xfactor$cor)))
            if(length(ff.factor)>0 && any(ff.factor %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) for the variance-covariance structure \"",paste(ff.factor[ff.factor %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            ## convert to factor with the right levels
            if(length(ff.factor)>0){
                for(iVar in ff.factor){ ## iVar <- ff.factor[1]
                    iLevel <- unique(c(object$xfactor$var[[iVar]],object$xfactor$cor[[iVar]]))
                    if(any(data[[iVar]] %in% iLevel == FALSE)){
                        Wf <- setdiff(unique(data[[iVar]]), iLevel)
                        stop("Unknown factor(s) \"",paste0(Wf,collapse="\" \""),"\" for variable \"",iVar,"\".\n",
                             "Valid factors: \"",paste0(iLevel, collapse="\" \""),"\".\n")
                    }
                    data.var[[iVar]] <- factor(data[[iVar]], levels = iLevel)
                }
            }

            if(cluster.var %in% names(data.var)){
                ## remove extra levels
                if(is.factor(data.var[[cluster.var]])){
                    data.var[[cluster.var]] <- droplevels(data.var[[cluster.var]])
                }
            }else if(all(time.var %in% names(data.var)) && all(!duplicated(data.var[,time.var,drop=FALSE]))){
                ## add cluster if no time repetition
                data.var[[cluster.var]] <- object$cluster$levels[1]
            }

            ## necessary to add the cluster index/time index/strata index
            data.var <- .prepareData(data.var,
                                     var.cluster = cluster.var,
                                     var.time = time.var,
                                     var.strata = strata.var,
                                     missing.repetition = NULL,
                                     droplevels = FALSE)

            ff.allvars <- unique(all.vars(object$formula$var,object$formula$cor))
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) for the covariance structure \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

                

            ## run part of .model.matrix.lmm
            outInit <- .extractIndexData(data = data.var, structure = object$design$vcov)
            U.cluster <- outInit$U.cluster
            U.time <- outInit$U.time
            U.strata <- outInit$U.strata
            design$index.cluster <- outInit$index.cluster
            design$index.clusterStrata <- outInit$index.clusterStrata
            design$index.clusterTime <- outInit$index.clusterTime
            design$cluster <- list(n = length(design$index.cluster), levels = levels(data.var$XXclusterXX), levels.original = NULL, nobs = sapply(design$index.cluster,length))

            ## copy structure 
            design$vcov <- object$design$vcov
            design$vcov$X <- list(var = NULL,
                                  cor = NULL,
                                  Upattern = NULL,
                                  Xpattern.var = NULL,
                                  Xpattern.cor = NULL,
                                  pattern.cluster = NULL,
                                  pair.varcoef = object$design$vcov$X$pair.varcoef)

            ## form design matrix    
            outDesign <- .vcov.matrix.lmm(structure = object$design$vcov, data = data.var, index.cluster = outInit$index.cluster)
            design$vcov$X$var <- outDesign$var
            design$vcov$X$cor <- outDesign$cor
            ## handle the case where structure is UN even though each cluster contain a single observation
            if(is.null(object$design$vcov$X$cor) && !is.null(outDesign$cor)){
                stop("Cannot build the design matrix for the covariance structure of clusters containing several obserations. \n",
                     "Whens fitting the model there were not replicates within cluster, so the correlation structure is unknown. \n")
            }

            ## find covariance patterns
            design$vcov <- .findUpatterns(design$vcov,
                                          index.clusterTime = outInit$index.clusterTime, U.time = U.time,
                                          index.cluster = outInit$index.cluster, U.cluster = U.cluster,
                                          index.clusterStrata = outInit$index.clusterStrata, U.strata = U.strata)            
        }
    }else{
        design <- object$design
    }

    ## ** export
    if(simplifies){
        if("mean" %in% effects && "variance" %in% effects){
            return(list(mean = design$mean,
                        vcov = design$vcov))
        }else if("mean" %in% effects){
            return(design$mean)
        }else if("variance" %in% effects){
            return(design$vcov)
        }
    }else{
        return(design)
    }
    
}


## * .mean.matrix.lmm
.mean.matrix.lmm <- function(formula, colnames, data,
                             stratify, name.strata, U.strata){


    ## ** design matrix
    if(stratify && length(U.strata)>1){
        n.strata <- length(U.strata)

        ## *** generate design matrix for each strata
        ls.X.mean <- lapply(U.strata, function(iS){ ## iS <- U.strata[1]
            if(is.null(colnames)){
                iX <- .model.matrix_regularize(formula, data[data$XXstrataXX==iS,,drop=FALSE])
            }else{
                iX <- stats::model.matrix(formula, data[data$XXstrataXX==iS,,drop=FALSE])
            }
            if(name.strata!="XXstrataXX"){
                colnames(iX) <- paste0(colnames(iX),":",name.strata,iS)
            }else{
                colnames(iX) <- paste0(colnames(iX),":",iS)
            }
            attr(iX,"index") <- data[data$XXstrataXX==iS,"XXindexXX"]
            return(iX)
        })

        ## *** assemble
        X.mean <- as.matrix(Matrix::bdiag(ls.X.mean))[order(unlist(lapply(ls.X.mean, attr, "index"))),,drop=FALSE]
        colnames(X.mean) <- unlist(lapply(ls.X.mean,colnames))
        if(!is.null(colnames)){
            X.mean <- X.mean[,colnames,drop=FALSE]
        }
        attr(X.mean, "assign") <- as.vector(do.call(cbind,lapply(ls.X.mean,attr,"assign")))
        attr(X.mean, "variable") <- attr(ls.X.mean[[1]],"variable")

        if(is.null(colnames)){
            strata.mu <- unlist(lapply(1:n.strata, function(iStrata){stats::setNames(rep(iStrata, NCOL(ls.X.mean[[iStrata]])),colnames(ls.X.mean[[iStrata]]))}))
        }else{
            attr.assign <- attr(X.mean, "assign")[match(colnames(X.mean),colnames)]
            attr.variable <- attr(X.mean, "variable")[match(colnames(X.mean),colnames)]
            X.mean <- X.mean[,colnames,drop=FALSE]
            attr(X.mean, "assign") <- attr.assign
            attr(X.mean, "variable") <- attr.assign
        }
    }else{
        if(is.null(colnames)){
            X.mean <- .model.matrix_regularize(formula, data)
            strata.mu <- stats::setNames(rep(U.strata,NCOL(X.mean)), colnames(X.mean))
        }else{
            X.mean <-  stats::model.matrix(formula, data)[,colnames,drop=FALSE]
        }
    }

    ## ** export
    if(is.null(colnames)){
        attr(X.mean,"strata.mu") <- strata.mu
    }
    return(X.mean)

}

## * .vcov.matrix.lmm
## output observation specific design matrix (but no covariance pattern)
.vcov.matrix.lmm <- function(structure, data, index.cluster){

    cluster.var <- structure$name$cluster
    time.var <- structure$name$time
    strata.var <- structure$name$strata
    U.cluster <- sort(unique(data[[cluster.var]]))
    U.time <- sort(unique(data[[time.var]]))
    U.strata <- sort(unique(data[[strata.var]]))
    
    ## ** normalize arguments
    ## strata
    if(is.na(strata.var)){
        strata.var <- "XXstrata.indexXX"
    }
    if(identical(strata.var, "XXstrata.indexXX") && strata.var %in% names(data) == FALSE){
        data$XXstrata.indexXX <- rep(1, NROW(data))
    }
    n.strata <- length(U.strata)
    
    ## time
    n.time <- length(U.time)

    ## formula
    formula.var <- structure$formula$var
    formula.cor <- structure$formula$cor

    ##  heterogeneous
    if(!is.null(structure$heterogeneous)){
        heterogeneous <- structure$heterogeneous
    }else{
        heterogeneous <- TRUE
    }

    ## data
    dataVar <- data
    if(length(all.vars(formula.var))>0 && structure$type %in% c("ID","IND","CS","UN")){
        for(iVar in all.vars(formula.var)){
            if(heterogeneous == FALSE){
                if(iVar == strata.var){
                    dataVar[[iVar]] <- as.factor(data[[iVar]])
                }else if(is.logical(data[[iVar]])){
                    dataVar[[iVar]] <- as.numeric(data[[iVar]]) + 1
                }else if(!is.numeric(data[[iVar]])){
                    dataVar[[iVar]] <- as.numeric(as.factor(data[[iVar]]))
                }else if(is.numeric(data[[iVar]])){
                    dataVar[[iVar]] <- data[[iVar]] - min(data[[iVar]]) + 1
                }
            }else if(heterogeneous){
                dataVar[[iVar]] <- as.factor(data[[iVar]])
            }
        }
    }
    dataCor <- data
    if(length(all.vars(formula.cor))>0 && structure$type %in% c("ID","IND","CS","UN")){
        for(iVar in all.vars(formula.cor)){
            if(heterogeneous == FALSE){
                if(iVar == strata.var){
                    dataCor[[iVar]] <- as.factor(data[[iVar]])
                }else if(is.logical(data[[iVar]])){
                    dataCor[[iVar]] <- as.numeric(data[[iVar]]) + 1
                }else if(!is.numeric(data[[iVar]])){
                    dataCor[[iVar]] <- as.numeric(as.factor(data[[iVar]]))
                }else if(is.numeric(data[[iVar]])){
                    dataCor[[iVar]] <- data[[iVar]] - min(data[[iVar]]) + 1
                }
            }else if(heterogeneous){
                dataCor[[iVar]] <- as.factor(data[[iVar]])
            }
        }
    }
    ## ** design matrix
    out <- list(var = NULL, cor = NULL, xfactor = list(var = NULL, cor = NULL))
    if(is.null(structure$param)){ ## structure
        out$var <- .colnameOrder(.model.matrix_regularize(formula.var, data = dataVar, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)
        out$xfactor$var <- stats::.getXlevels(stats::terms(formula.var),stats::model.frame(formula.var,dataVar))
        if(!is.null(formula.cor) && n.time>1 && any(sapply(index.cluster,length)>1)){  ## at least one individual with more than timepoint
            out$cor <- .colnameOrder(.model.matrix_regularize(formula.cor, data = dataCor, augmodel = TRUE), strata.var = strata.var, n.strata = n.strata)
            out$xfactor$cor <- stats::.getXlevels(stats::terms(formula.cor),stats::model.frame(formula.cor,dataCor)) 
        }
    }else{ ## newdata        
        out$var <- stats::model.matrix(formula.var, data = dataVar)[,attr(structure$X$var,"original.colnames"),drop=FALSE]
        colnames(out$var) <- colnames(structure$X$var)
        for(iAssign in setdiff(names(attributes(structure$X$var)), c("dim","dimnames"))){
            attr(out$var,iAssign) <- attr(structure$X$var,iAssign)
        }

        if(!is.null(formula.cor) && n.time>1 && any(sapply(index.cluster,length)>1)){  ## at least one individual with more than timepoint
            out$cor <- stats::model.matrix(formula.cor, data = dataCor)[,attr(structure$X$cor,"original.colnames"),drop=FALSE]
            colnames(out$cor) <- colnames(structure$X$cor)
            for(iAssign in setdiff(names(attributes(structure$X$cor)), c("dim","dimnames"))){
                attr(out$cor,iAssign) <- attr(structure$X$cor,iAssign)
            }
        }
    }

    ## ** export
    return(out)
}

## * .model.matrix.lmm
.model.matrix.lmm <- function(formula.mean, structure,
                              data, var.outcome, var.weights,
                              stratify.mean,
                              precompute.moments){

    ## ** indexes
    outInit <- .extractIndexData(data = data, structure = structure)
    U.cluster <- outInit$U.cluster
    U.time <- outInit$U.time
    U.strata <- outInit$U.strata
    index.cluster <- outInit$index.cluster
    index.clusterStrata <- outInit$index.clusterStrata
    index.clusterTime <- outInit$index.clusterTime

    ## ** mean
    ## use stats::model.frame to handle splines
    data.mf <- stats::model.frame(stats::update(formula.mean,~.+XXindexXX+XXtimeXX+XXclusterXX+XXstrataXX),data)
    X.mean <- .mean.matrix.lmm(formula = formula.mean, colnames = NULL, data = data.mf,
                               stratify = stratify.mean, name.strata = structure$name$strata, U.strata = U.strata)  ## only stratify mean if gls optimizer
    strata.mu <- attr(X.mean,"strata.mu")
    attr(X.mean,"strata.mu") <- NULL
    attr(X.mean,"terms") <- attr(data.mf,"terms")

    ## ** variance

    ## *** design matrix
    outDesign <- .vcov.matrix.lmm(structure = structure, data = data, index.cluster = outInit$index.cluster)

    structure$xfactor <- outDesign$xfactor
    structure$X <- list(var = outDesign$var,
                        cor = outDesign$cor)

    ## *** parametrization and patterns
    structure <- .skeleton(structure = structure, data = data, indexData = outInit)

    ## *** covariance pattern
    structure <- .findUpatterns(structure,
                                index.clusterTime = outInit$index.clusterTime, U.time = U.time,
                                index.cluster = outInit$index.cluster, U.cluster = U.cluster,
                                index.clusterStrata = outInit$index.clusterStrata, U.strata = U.strata)

    ## ** prepare calculation of the score
    if(precompute.moments){
        if(is.na(var.weights[1])){
            wX.mean <- X.mean
            wY <- cbind(data[[var.outcome]])
        }else{
            wX.mean <- sweep(X.mean, FUN = "*", MARGIN = 1, STATS = sqrt(data[[var.weights[1]]]))
            wY <- cbind(data[[var.outcome]]*sqrt(data[[var.weights[1]]]))
        }
        precompute.XX <-  .precomputeXX(X = wX.mean, pattern = structure$X$Upattern$name, 
                                        pattern.ntime = stats::setNames(structure$X$Upattern$n.time, structure$X$Upattern$name),
                                        pattern.cluster = structure$X$Upattern$index.cluster, index.cluster = index.cluster)
        precompute.XY <-  .precomputeXR(X = precompute.XX$Xpattern, residuals = wY, pattern = structure$X$Upattern$name,
                                        pattern.ntime = stats::setNames(structure$X$Upattern$n.time, structure$X$Upattern$name),
                                        pattern.cluster = structure$X$Upattern$index.cluster, index.cluster = index.cluster)
    }else{
        precompute.XX <- NULL
        precompute.XY <- NULL
    }

    ## ** pairs
    structure$X$pair.varcoef <- stats::setNames(lapply(structure$X$Upattern$name, function(iPattern){## iPattern <- structure$X$Upattern$name[1]
        iParamVar <- structure$X$Upattern$param[[iPattern]]

        iOut <- .unorderedPairs(iParamVar)
        attr(iOut, "key") <- matrix(NA, nrow = length(iParamVar), ncol = length(iParamVar), dimnames = list(iParamVar,iParamVar))
        for(iCol in 1:NCOL(iOut)){
            attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
            attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol]] <- iCol
        }
        return(iOut)
    }),structure$X$Upattern$name)

    pair.meanvarcoef <- stats::setNames(lapply(structure$X$Upattern$name, function(iPattern){ ## iPattern <- structure$X$Upattern$name[1]
        if(stratify.mean){
            iParamMu <- names(strata.mu)[strata.mu==structure$X$Upattern$index.strata[iPattern]]
        }else{
            iParamMu <- names(strata.mu)
        }
        iParamVar <- structure$X$Upattern$param[[iPattern]]
        iOut <- unname(t(expand.grid(iParamMu, iParamVar)))
        return(iOut)
    }), structure$X$Upattern$name)

    ## ** param
    ls.sub <- lapply(as.list(c("(Intercept)",attr(formula.mean,"term.labels"))[attr(X.mean,"assign")+1]),
                     function(iSub){
                         if(grepl(":",iSub,fixed = TRUE)){ ## wraper from a package, e.g. stats::poly()
                             unlist(strsplit(iSub,"::",fixed = TRUE))[-1]
                         }else if(grepl(":",iSub,fixed = TRUE)){ ## interaction
                             unlist(strsplit(iSub,":",fixed = TRUE))
                         }else{
                             return(iSub)
                         }
                     })
    mu.level <- mapply(sub = ls.sub,
                       name = colnames(X.mean),
                       FUN = function(sub,name){
                           for(iSub in sub){
                               name <- gsub(iSub,"",name, fixed = TRUE)
                           }
                           return(name)
                       })
    skeleton.param <- rbind(data.frame(name = colnames(X.mean),
                                       strata = NA,
                                       type = "mu",
                                       level = gsub("^:","",gsub(":$","",mu.level)),
                                       code = NA,
                                       code.x = NA,
                                       code.y = NA,
                                       sigma = NA,
                                       k.x = NA,
                                       k.y = NA),
                            structure$param)

    if(stratify.mean){
        skeleton.param$strata[skeleton.param$type=="mu"] <- strata.mu
        skeleton.param$strata <- as.list(skeleton.param$strata)
    }else{
        skeleton.param$strata[skeleton.param$type=="mu"] <- list(as.numeric(unique(stats::na.omit(skeleton.param$strata))))
    }
    rownames(skeleton.param) <- NULL
    attr(skeleton.param, "pair.meanvarcoef") <- pair.meanvarcoef

    order.type <- factor(skeleton.param$type, c("mu", "sigma", "k", "rho"))
    if(is.unsorted(order.type)){ ## in presence of strata, typically need to reorder by type
        skeleton.param <- skeleton.param[order(order.type),,drop=FALSE]
    }

    ## ** gather and export
    out <- list(mean = X.mean,
                vcov = structure,
                Y = data[[var.outcome]],
                precompute.XX = precompute.XX,
                precompute.XY = precompute.XY,
                index.cluster = index.cluster,
                index.clusterTime = index.clusterTime,
                index.clusterStrata = index.clusterStrata,
                time = list(n = max(unlist(index.clusterTime)), levels = levels(data$XXtimeXX), levels.original = NULL, nobs = table(unlist(index.clusterTime))),
                cluster = list(n = length(index.cluster), levels = levels(data$XXclusterXX), levels.original = NULL, nobs = sapply(index.cluster,length)),
                param = skeleton.param
                )
    if(!is.na(var.weights[1])){
        out$weights <- data[[var.weights[1]]]
        ## NOTE: only take first weight for each cluster as weights should be constant within cluster
        ## out$weights <- sapply(index.cluster, function(iIndex){data[iIndex[1],var.weights[1]]})
    }
    if(!is.na(var.weights[2])){
        out$scale.Omega <- data[[var.weights[2]]]
        ## NOTE: only take first weight for each cluster as weights should be constant within cluster
        ## out$scale.Omega <- sapply(index.cluster, function(iIndex){data[iIndex[1],var.weights[2]]})
    }
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

    if(X.qr$rank==NCOL(X)){
        if(augmodel){
            return(.augmodel.matrix(stats::delete.response(stats::terms(formula)),data))
        }else{
            attr(X,"variable") <- all.vars(formula)
            return(X)
        }
    }else{
        if(LMMstar.options()$drop.X==FALSE){
            stop("The design matrix does not have full rank according to the QR decomposition. \n")
        }
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
            message("Factor variable(s) with empty level: \"",paste(name.factor[test.factor==FALSE], collapse = "\" \""),"\"\n ",
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
    ## vec.interactionName <- apply(M.interactionName, MARGIN = 1, FUN = function(iRow){paste(stats::na.omit(iRow), collapse = ":")})

    ## ** test 2: remove column(s) corresponding to the same level
    index.keep <- which(!duplicated(attr(X,"M.level")))
    Xsave <- X
    X <- Xsave[,index.keep,drop=FALSE]
    attrX <- attributes(Xsave)[setdiff(names(attributes(Xsave)),names(attributes(X)))]
    attrX$assign <- attrX$assign[index.keep]
    attrX$variable <- attrX$variable
    attrX$term.labels <- attrX$term.labels[index.keep]
    attrX$order <- attrX$order[index.keep]
    attrX$ls.level <- attrX$ls.level[index.keep]
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
            iMlevel <- X.Mlevel[apply(X.Mlevel[,iNoVar,drop=FALSE],1,function(iParam){all(iParam==X.reference[iNoVar])}),iVar,drop=FALSE]
        }else{
            iMlevel <- X.Mlevel
        }

        ## create all possible values for the interaction and check if it is constant
        test.constant <- sapply(1:NROW(iMlevel), FUN = function(iParam){ ## iParam <- 1
            iValue <- do.call(cbind,lapply(var.factor, function(iFactor){data[[iFactor]]==iMlevel[iParam,iFactor]}))
            iNumeric <- intersect(iVar, var.numeric)
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
        X.old <- X
        test.keep <- colnames(X.old) %in% setdiff(colnames(X.old),rmX)
        X <- X.old[,test.keep,drop=FALSE]
        if(qr(X)$rank==X.qr$rank){
            message("Constant values in the design matrix in interactions \"",paste(names(ls.rmX), collapse = "\" \""),"\"\n ",
                    "Coefficients \"",paste(unique(rmX), collapse = "\" \""),"\" have been removed. \n")
        }else{
            warning("Constant values in the design matrix in interactions \"",paste(names(ls.rmX), collapse = "\" \""),"\"\n ",
                    "Coefficients \"",paste(rmX, collapse = "\" \""),"\" have been removed. \n")
        }
        attr(X,"assign") <- attr(X.old,"assign")[test.keep] ## as.numeric(as.factor(attr(X.old,"assign")[test.keep])) - "(Intercept)" %in% colnames(X)
        attr(X,"contrasts") <- attr(X.old,"contrasts")
        if(augmodel || X.qr$rank!=NCOL(X.qr$qr)){
            attr(X,"formula") <- attr(X.old,"formula")
            attr(X,"variable") <- attr(X.old,"variable")
            attr(X,"term.labels") <- attr(X.old,"term.labels")[test.keep]
            attr(X,"order") <- attr(X.old,"order")[test.keep]
            attr(X,"ls.level") <- attr(X.old,"ls.level")[test.keep]
            attr(X,"M.level") <- attr(X.old,"M.level")[test.keep,,drop=FALSE]
            attr(X,"reference.level") <- attr(X.old,"reference.level")
        }
    }

    ## ** identify if there is still an identifiability problem
    if(X.qr$rank==NCOL(X)){
        return(X)
    }
    if(any(attr(X,"order")>2)){
        stop("Cannot handle interaction involving more than two variables. \n")
    }

    ## ** test 4: remove columns one at a time, starting from the higher order term (i.e. the last columns)
    score.ref <- apply(attr(X,"M.level"),1,function(iRow){sum(iRow==attr(X,"reference"))})
    test.name <- colnames(X)[order(attr(X,"order")+score.ref/(length(attr(X,"reference"))+1), decreasing = TRUE)]
    
    rmX <- NULL
    for(iCol in test.name){
        iIndex <- colnames(X)!=iCol
        X.test <- X[,iIndex,drop=FALSE]
        if(X.qr$rank == qr(X.test)$rank){
            rmX <- c(rmX, iCol)
            keep.attr <- attributes(X)
            X <- X.test
            attr(X,"assign") <- keep.attr$assign[iIndex]
            attr(X,"contrast") <- keep.attr$contrasts
            if(augmodel){
                attr(X,"formula") <- keep.attr$formula
                attr(X,"variable") <- keep.attr$variable
                attr(X,"term.labels") <- keep.attr$term.labels[iIndex]
                attr(X,"order") <- keep.attr$order[iIndex]
                attr(X,"M.level") <- keep.attr$M.level[iIndex,,drop=FALSE]
                attr(X,"reference.level") <- keep.attr$reference.level
            }
        }
        if(NCOL(X)==X.qr$rank){break}
    }
    if(length(rmX)>0){
        message("Singular design matrix, coefficient",if(length(rmX)>1){"s"}," \"",paste(rmX, collapse = "\" \""),"\" has been removed. \n")
    } else if(!augmodel){
        attr(X,"formula") <- NULL
        attr(X,"term.labels") <- NULL
        attr(X,"order") <- NULL
        attr(X,"ls.level") <- NULL
        attr(X,"M.level") <- NULL
        attr(X,"reference.level") <- NULL
    }

    ## ** export
    return(X)
}

## ** .augmodel.matrix_match
## add more information about which variable contribute to which column in the design matrix
.augmodel.matrix <- function(formula, data){

    ## ** normalize formula
    formula.terms <- stats::delete.response(stats::terms(formula))
    var2term <- attr(formula.terms, "factors")

    ## if(NROW(var2term)>0){
    ##     test1fac <- sapply(data[,rownames(var2term),drop=FALSE], function(iVar){length(unique(iVar))})
    ##     if(any(test1fac==1)){
    ##         stop("Cannot create the design matrix due to a factor variable with a single level. \n",
    ##              "Variable: \"",paste(names(test1fac)[test1fac==1], collapse = "\" \""),"\". \n")
    ##     }
    ## }
    ## ** create design matrix
    X <- stats::model.matrix(formula,data)
    p <- NCOL(X)

    ## ** identify factors and reference level
    all.variable <- all.vars(formula.terms)
    contrast.variable <- stats::setNames(lapply(all.variable, function(iVar){
        if(is.character(data[[iVar]])){
            out <- stats::contrasts(as.factor(data[[iVar]]))
            if(any(colSums(abs(out)>1e-12)>1)){
                if(any(options()$contrast!=c("contr.treatment","contr.poly"))){
                    stop("Cannot handle contrasts involving simultaneously several levels. \n",
                         "Could be because options()$contrast has been modified to non-standard contrasts. \n")
                }else{
                    stop("Cannot handle contrasts involving simultaneously several levels. \n")
                }
            }
        }else if(is.factor(data[[iVar]])){
            out <- stats::contrasts(data[[iVar]])
            if(any(colSums(abs(out)>1e-12)>1)){
                if(any(options()$contrast!=c("contr.treatment","contr.poly"))){
                    stop("Cannot handle contrasts involving simultaneously several levels. \n",
                         "Could be because options()$contrast has been modified to non-standard contrasts. \n")
                }else{
                    stop("Cannot handle contrasts involving simultaneously several levels. \n")
                }
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
    reference2 <- data.frame(matrix(FALSE, nrow = 1, ncol = length(rownames(var2term)), dimnames = list(NULL, rownames(var2term))),stringsAsFactors = FALSE)
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
            X.level[[iCol]] <- data.frame(reference,stringsAsFactors = FALSE)
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
                iLs.grid <- expand.grid(iLs.factor[sapply(iLs.factor,length)>0])
                iLs.allnameInteraction <- interaction(stats::setNames(lapply(iVar,function(iName){paste0(iName,iLs.grid[[iName]])}), iVar), sep = ":")
                iIndex  <- which(X.names[iCol] == iLs.allnameInteraction)
                ## deduce the factor variable
                iValue.factor <- iLs.grid[iIndex,,drop=FALSE]
                X.level[[iCol]] <- as.data.frame(stats::setNames(lapply(iVar, function(iName){
                    if(iName %in% names(iValue.factor) == FALSE){
                        as.numeric(NA)
                    }else{
                        as.character(iValue.factor[1,iName])
                    }}), iVar))                
                X.level2[[iCol]][,names(X.level[[iCol]])] <- X.level[[iCol]]
                if(any(is.na(X.level[[iCol]]))){
                    X.level2[[iCol]][,is.na(X.level2[[iCol]])] <- TRUE
                }
            }else{
                X.level[[iCol]] <- data.frame(matrix(as.numeric(NA), nrow = 1, ncol = length(iVar), dimnames = list(NULL, iVar)),stringsAsFactors = FALSE)
                X.level2[[iCol]][,iVar] <- TRUE
            }
        }
    }
    ## ** add attributes
    attr(X,"formula") <- formula
    attr(X,"variable") <- colnames(reference2)
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

    time.var <- structure$name$time
    cluster.var <- structure$name$cluster
    strata.var <- structure$name$strata

    U.cluster <- sort(unique(data[[cluster.var]]))
    U.time <- sort(unique(data[[time.var]]))
    if(!is.na(strata.var)){
        U.strata <- sort(unique(data[[strata.var]]))
    }else{
        U.strata <- 1
    }

    ## *** find position of each cluster
    index.cluster <- tapply(1:NROW(data),data[[cluster.var]],function(iI){iI}, simplify = FALSE)
    
    ## *** find time corresponding to each cluster
    index.clusterTime <- tapply(data[[time.var]],data[[cluster.var]],function(iT){as.numeric(factor(iT,levels = U.time))})

    ## *** re-order according to time
    order.clusterTime <- lapply(index.clusterTime, order)
    index.cluster <- mapply(x = index.cluster, y = order.clusterTime, function(x,y){x[y]}, SIMPLIFY = FALSE)
    attr(index.cluster, "vectorwise") <- data[[cluster.var]]
    index.clusterTime <- mapply(x = index.clusterTime, y = order.clusterTime, function(x,y){x[y]}, SIMPLIFY = FALSE)
    attr(index.clusterTime, "vectorwise") <- as.numeric(factor(data[[time.var]], levels = U.time))

    ## *** find strata corresponding to each cluster
    if(!is.na(strata.var)){
        iMatch <- tapply(as.character(data[[strata.var]]),data[[cluster.var]],unique)
        index.clusterStrata <- stats::setNames(match(iMatch, U.strata), names(iMatch))
    }else{
        index.clusterStrata <- stats::setNames(rep(1, length(U.cluster)), U.cluster)
    }

    ## *** export
    out <- list(U.cluster = U.cluster,
                U.time = U.time,
                U.strata = U.strata,
                index.cluster = index.cluster,
                index.clusterTime = index.clusterTime,
                index.clusterStrata = index.clusterStrata
                )
    return(out)
}


##----------------------------------------------------------------------
### model.matrix.R ends here
