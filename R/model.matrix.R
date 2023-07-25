### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: jul 25 2023 (12:01) 
##           By: Brice Ozenne
##     Update #: 2872
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
model.matrix.lmm <- function(object, data = NULL, effects = "mean", simplify = TRUE, drop.X = NULL, ...){

    ## ** normalize user imput
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    if(is.null(drop.X)){
        drop.X <- LMMstar.options()$drop.X
    }
    if(!identical(effects,"index")){
        effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)
    }
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** update design matrix with new dataset
    if(!is.null(data)){

        data <- as.data.frame(data)

        ## *** prepare output
        design <- object$design
        design[setdiff(names(design),c("vcov","param","drop.X"))] <- list(NULL)
        design$vcov[c("var","cor","pattern","Upattern")] <- NULL
        var.outcome <- object$outcome$var
        var.cluster <- attr(object$cluster$var,"original")
        var.time <- attr(object$time$var,"original")
        var.strata <- attr(object$strata$var,"original")

        ## *** outcome
        if(!simplify){
            if(any(object$outcome$var %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing outcome variable \"",object$outcome$var,"\".\n")
            }
            design$Y <- data[[object$outcome$var]]
        }

        ## *** index
        structure.var <- stats::na.omit(c(attr(object$cluster$var,"original"),attr(object$time$var,"original"),attr(object$strata$var,"original")))
        if(("index" %in% effects) || ("variance" %in% effects) || ("correlation" %in% effects) || !simplify){

            ## check dataset
            if(any(structure.var %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) for the data structure \"",paste(structure.var[structure.var %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }
            ## generate cluster/time/strata variables
            data.Nindex <- .lmmNormalizeData(data = data[stats::na.omit(c(var.cluster,var.time,var.strata))],
                                             var.outcome = NA,
                                             var.cluster = var.cluster,
                                             var.time = var.time,
                                             var.strata = var.strata,
                                             droplevels = list(time = object$time$levels,
                                                               strata = object$strata$levels),
                                             initialize.cluster = object$design$vcov$ranef$crossed,
                                             initialize.time = setdiff(object$design$vcov$ranef$vars, object$design$vcov$ranef$var.cluster))$data
            
            ## extract indexes
            outInit <- .extractIndexData(data = data.Nindex, structure = design$vcov)
            if(identical(effects,"index")){
                return(outInit)
            }
            design$index.cluster <- outInit$index.cluster
            design$index.clusterStrata <- outInit$index.clusterStrata
            design$index.clusterTime <- outInit$index.clusterTime
            
        }
        
        ## *** weights
        if(!simplify){
            if(!is.na(object$weight$var[1])){
                if(any(object$weight$var[1] %in% names(data) == FALSE)){
                    stop("Incorrect argument \'data\': missing weight variable \"",object$weight$var[1],"\".\n")
                }
                design$weight <- data[[object$weight$var[1]]]
            }
        }
        if(!simplify){
            if(!is.na(object$weight$var[2])){
                if(any(object$weight$var[2] %in% names(data) == FALSE)){
                    stop("Incorrect argument \'data\': missing weight variable \"",object$weight$var[2],"\".\n")
                }
                design$scale.Omega <- data[[object$weight$var[2]]]
            }
        }
        
        ## *** mean
        if("mean" %in% effects){
            
            ## check dataset
            ff.allvars <- all.vars(object$formula$mean.design)
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) for the mean structure \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            ## convert to factor with the right levels
            data.mean <- .updateFactor(data, xfactor = object$xfactor$mean)

            ## use stats::model.frame to handle spline
            data.mf.mean <- stats::model.frame(object$formula$mean.design, data = data.mean, na.action = stats::na.pass)
            design$mean  <- stats::model.matrix(object$formula$mean.design, data.mf.mean)[,colnames(object$design$mean),drop=FALSE]
        }

        ## *** variance-covariance
        if("variance" %in% effects){
            ## check dataset
            ff.allvars <- c(all.vars(object$formula$var),all.vars(object$formula$cor))
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) for the variance-covariance structure \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            ## update levels
            design$vcov$var <- list(data = .updateFactor(data, xfactor = object$xfactor$var),
                                    lp2X = object$design$vcov$var$lp2X,
                                    lp2data = object$design$vcov$var$lp2data,
                                    X.attr = attributes(object$design$vcov$var$X))
            design$vcov$cor <- list(data = .updateFactor(data, xfactor = object$xfactor$cor),
                                    lp2X = object$design$vcov$cor$lp2X,
                                    lp2data = object$design$vcov$cor$lp2data,
                                    X.attr = attributes(object$design$vcov$cor$X))
            
            ## form design matrix
            outDesign <- .vcov.matrix.lmm(structure = design$vcov, data = data.Nindex, index.cluster = outInit$index.cluster, drop.X = drop.X)
            design$vcov$var <- outDesign$var
            design$vcov$cor <- outDesign$cor

            ## handle the case where structure is UN even though each cluster contain a single observation
            if(is.null(object$design$vcov$cor$X) && !is.null(outDesign$cor$X)){
                stop("Cannot build the design matrix for the covariance structure of clusters containing several obserations. \n",
                     "When fitting the model there were not replicates within cluster, so the correlation structure is unknown. \n")
            }

            ## find covariance patterns
            design$vcov <- .findUpatterns(design$vcov,
                                          index.clusterTime = outInit$index.clusterTime, U.time = outInit$U.time,
                                          index.cluster = outInit$index.cluster, U.cluster = outInit$U.cluster,
                                          index.clusterStrata = outInit$index.clusterStrata, U.strata = outInit$U.strata)            
        }
    }else{
        design <- object$design
    }

    ## ** export
    if(simplify){
        design2 <- list(mean = design$mean,
                        variance = design$vcov$var$X,
                        correlation = design$vcov$cor$X)
        attributes(design2$mean) <- attributes(design$mean)[c("dim","dimnames","assign")]
        attributes(design2$variance) <- attributes(design$variance)[c("dim","dimnames","assign")]
        attributes(design2$correlation) <- attributes(design$correlation)[c("dim","dimnames","assign")]
        if(length(effects)>1){
            return(design2[effects])
        }else{
            return(design2[[effects]])
        }
    }else{
        return(design)
    }
    
}

## * .vcov.matrix.lmm
## output observation specific design matrix (but no covariance pattern)
.vcov.matrix.lmm <- function(structure, data, index.cluster, drop.X){

    sep <- LMMstar.options()$sep["lp"]
    var.cluster <- structure$name$cluster
    var.time <- structure$name$time
    var.strata <- structure$name$strata
    U.cluster <- sort(unique(data[[var.cluster]]))
    U.time <- sort(unique(data[[var.time]]))
    structure.class <- structure$class
    structure.type <- structure$type
    formula <- structure$formula

    out <- list(var = list(), cor = list(),
                xfactor = list(var = NULL, cor = NULL))
    test.newdata <- !is.null(structure$param) ## if TRUE then model.matrix call on newdata otherwise call from lmm->.model.matrix on original data

    ## ** normalize dataset
    ## time
    n.time <- length(U.time)

    ## formula
    out$var$formula <- structure$formula$var
    out$cor$formula <- structure$formula$cor
    
    ## data: convert variable to numeric/factor according to the structure
    if(is.null(structure$var$data)){ ## from .model.matrix
        out$var$data <- data
        if(length(all.vars(out$var$formula))>0 && structure.class %in% c("ID","IND","CS","RE","UN","TOEPLITZ")){
            for(iVar in all.vars(out$var$formula)){
                out$var$data[[iVar]] <- as.factor(data[[iVar]])
            }
        }
    }else{ ## from model.matrix.lmm
        out$var$data <- structure$var$data
    }

    if(is.null(structure$cor$data)){ ## from .model.matrix
        out$cor$data <- data    
        if(length(all.vars(out$cor$formula))>0 && structure.class %in% c("CS","RE","UN","TOEPLITZ")){
            for(iVar in all.vars(out$cor$formula)){
                if((structure.class=="UN") || (structure.class=="CS" && structure.type %in% "heterogeneous") || (!is.na(var.strata) && iVar == var.strata)){
                    out$cor$data[[iVar]] <- as.factor(data[[iVar]])
                }else if(is.logical(data[[iVar]])){ ## "TOEPLITZ" or homogeneous "CS"
                    out$cor$data[[iVar]] <- as.numeric(data[[iVar]]) + 1
                }else if(!is.numeric(data[[iVar]])){ ## "TOEPLITZ" or homogeneous "CS"
                    out$cor$data[[iVar]] <- as.numeric(as.factor(data[[iVar]]))
                }else if(is.numeric(data[[iVar]])){ ## "TOEPLITZ" or homogeneous "CS"
                    out$cor$data[[iVar]] <- data[[iVar]] - min(data[[iVar]]) + 1
                }
            
            }
        }
    }else{ ## from model.matrix.lmm
        out$cor$data <- structure$cor$data
    }

    for(iMoment in c("var","cor")){ ## iMoment <- "var"

        iMoment.txt <- c(var = "variance", cor = "correlation")[iMoment]
        iData <- out[[iMoment]]$data
        iU.strata <- levels(iData[[var.strata]])
        iFormula <- out[[iMoment]]$formula

        ## ** design matrix
        if(test.newdata == FALSE){ ## structure
            if(!is.null(iFormula)){
                if(structure.class == "CUSTOM"){
                    out[[iMoment]]$X <- iData[,all.vars(iFormula),drop=FALSE]
                }else{
                    out[[iMoment]]$X <- .model.matrix_regularize(iFormula, data = iData, augmodel = TRUE, type = iMoment.txt, drop.X = drop.X)                    
                }
                attr(out[[iMoment]]$X, "original.colnames") <- colnames(out[[iMoment]]$X)
                out$xfactor[[iMoment]] <- stats::.getXlevels(stats::terms(iFormula),stats::model.frame(iFormula,iData))
            }

        }else{ ## newdata

            if(!is.null(iFormula)){
                iAttr.X <- structure[[iMoment]]$X.attr
                iKeep.attr <- setdiff(names(iAttr.X), c("dim","dimnames"))

                out[[iMoment]]$X <- stats::model.matrix(iFormula, data = iData)[,iAttr.X$original.colnames,drop=FALSE]
                colnames(out[[iMoment]]$X) <- iAttr.X$dimnames[[2]]
                for(iAssign in iKeep.attr){
                    attr(out[[iMoment]]$X,iAssign) <- iAttr.X[[iAssign]]
                }
            }
        }

        ## ** linear predictors & patterns
        if(!is.null(iFormula)){

            #### linear predictor for each observation
            iLp <- nlme::collapse(out[[iMoment]]$X, sep = sep, as.factor = !test.newdata)
            
            if(!test.newdata){
                #### position of the observations with distinct linear predictors
                iIndex.Ulp <- which(!duplicated(iLp))
                #### convert linear predictor to numeric, possible updating the order of the level to match "natural ordering"
                if(NCOL(out[[iMoment]]$X)>1){
                    if(is.na(var.strata)){
                        iNewlevel <- orderLtoR(out[[iMoment]]$X[iIndex.Ulp,,drop=FALSE], strata = NULL)
                    }else{                    
                        iNewlevel <- orderLtoR(out[[iMoment]]$X[iIndex.Ulp,,drop=FALSE], strata = factor(attr(out[[iMoment]]$X,"M.level")[[var.strata]], iU.strata))
                        ##  out[[iMoment]]$X[iIndex.Ulp[order(iNewlevel)],,drop=FALSE]
                    }
                    out[[iMoment]]$lp <- as.numeric(factor(iLp, levels = iLp[iIndex.Ulp][order(iNewlevel)]))
                }else{
                    iNewlevel <- levels(iLp)
                    out[[iMoment]]$lp <- as.numeric(iLp)
                }
                #### re-order position of the observations with distinct linear predictors to match increasing linear predictors
                if(length(iIndex.Ulp)>1){
                    iIndex.Ulp <- iIndex.Ulp[order(out[[iMoment]]$lp[iIndex.Ulp], decreasing = FALSE)]
                }
                #### re-express the linear predictor in term of design matrix or original data
                out[[iMoment]]$lp2X <- out[[iMoment]]$X[iIndex.Ulp,,drop=FALSE]
                rownames(out[[iMoment]]$lp2X) <- iLp[iIndex.Ulp]
                out[[iMoment]]$lp2data <- data[iIndex.Ulp,attr(out[[iMoment]]$X,"variable"),drop=FALSE]
                rownames(out[[iMoment]]$lp2data) <- iLp[iIndex.Ulp]
            }else{
                out[[iMoment]]$lp <- as.numeric(factor(iLp, levels = rownames(structure[[iMoment]]$lp2X)))
                out[[iMoment]]$lp2X <- structure[[iMoment]]$lp2X
                out[[iMoment]]$lp2data <- structure[[iMoment]]$lp2data
            }
            #### pattern for each cluster
            out[[iMoment]]$pattern <- as.numeric(as.factor(sapply(index.cluster, function(iC){paste(out[[iMoment]]$lp[iC], collapse = ".")})))
            ## Note: tapply(out[[iMoment]]$lp, attr(index.cluster,"vectorwise"), paste, collapse = ".")
            ## does not work when the time are not ordered in the initial dataset
            if(is.null(structure[[iMoment]]$pattern2lp)){
                #### position of the cluster with distinct patterns
                iIndex.pattern <- which(!duplicated(out[[iMoment]]$pattern))
                iIndex.pattern <- iIndex.pattern[order(out[[iMoment]]$pattern[iIndex.pattern], decreasing = FALSE)] ## re-order
                #### re-express the pattern in term of linear predictors
                out[[iMoment]]$pattern2lp <- unname(lapply(index.cluster[iIndex.pattern], function(iIndex){out[[iMoment]]$lp[iIndex]}))
            }else{
                out[[iMoment]]$pattern2lp <- structure[[iMoment]]$pattern2lp
            }
        }
    }

    ## ** export
    out$var$formula <- NULL
    out$var$data <- NULL
    out$cor$formula <- NULL
    out$cor$data <- NULL
    ## out$var
    ## out$cor
    return(out)
}

## * .model.matrix.lmm
.model.matrix.lmm <- function(formula.mean, structure,
                              data, var.outcome, var.weights,
                              drop.X, ## drop singular component of the design matrix
                              precompute.moments){

    ## ** indexes
    outInit <- .extractIndexData(data = data, structure = structure)
    U.cluster <- outInit$U.cluster
    U.time <- outInit$U.time
    U.strata <- outInit$U.strata
    n.strata <- length(U.strata)
    index.cluster <- outInit$index.cluster
    index.clusterStrata <- outInit$index.clusterStrata
    index.clusterTime <- outInit$index.clusterTime

    ## ** mean
    ## use stats::model.frame to handle splines
    data.mf <- stats::model.frame(stats::update(formula.mean,~.+XXindexXX+XXtimeXX+XXclusterXX+XXstrataXX),data)
    X.mean <- .model.matrix_regularize(formula.mean, data = data.mf, type = "mean", drop.X = drop.X)
    attr(X.mean,"terms") <- attr(data.mf,"terms")

    if(NCOL(X.mean)>0){
        ## gather all terms
        all.terms <- c("(Intercept)",attr(stats::terms(formula.mean),"term.labels"))[attr(X.mean,"assign")+1]
        ## identify level
        ls.sub <- sapply(all.terms, function(iSub){ ## iSub <- "X1:X2"
            if(grepl("::",iSub,fixed = TRUE)){ ## wraper from a package, e.g. stats::poly()
                return(unlist(strsplit(iSub,"::",fixed = TRUE))[-1])
            }else if(grepl(":",iSub,fixed = TRUE)){ ## interaction
                return(unlist(strsplit(iSub,":",fixed = TRUE)))
            }else{
                return(iSub)
            }
        }, simplify = FALSE)
        mu.level <- mapply(sub = ls.sub,
                           name = colnames(X.mean),
                           FUN = function(sub,name){
                               for(iSub in sub){
                                   name <- gsub(iSub,"",name, fixed = TRUE)
                               }
                               return(name)
                           })
        ## skeleton
        skeleton.mu <- data.frame(name = colnames(X.mean),
                                  index.strata = NA,
                                  type = "mu",
                                  constraint = NA,
                                  level = gsub("^:","",gsub(":$","",mu.level)),
                                  code = NA,
                                  lp.x = NA,
                                  lp.y = NA,
                                  sigma = NA,
                                  k.x = NA,
                                  k.y = NA)

        if(n.strata==1){
            skeleton.mu$index.strata <- n.strata
        }else if(structure$name$strata %in% attr(X.mean,"variable")){
            ## identify when stratifying the mean structure 
            X.mean2 <- .augmodel.matrix(stats::delete.response(stats::terms(formula.mean)),data.mf)
            M.factors <- attr(attr(X.mean2,"formula"),"factors")
            term.strata <- names(which(M.factors[structure$name$strata,]>0))
            col.strata <- intersect(colnames(X.mean2)[attr(X.mean2,"term.labels") %in% term.strata],
                                    colnames(X.mean))
            skeleton.mu$index.strata[match(col.strata,skeleton.mu$name)] <- match(attr(X.mean2,"M.level")[col.strata,structure$name$strata], U.strata)
        }
        attr(X.mean,"strata") <- stats::setNames(skeleton.mu$index.strata, skeleton.mu$name)
    }else{
        skeleton.mu <- NULL
    }

    ## ** variance
    ## *** design matrix
    outDesign <- .vcov.matrix.lmm(structure = structure, data = data, index.cluster = outInit$index.cluster, drop.X = drop.X)
    structure$xfactor <- outDesign$xfactor
    structure$var <- outDesign$var
    structure$cor <- outDesign$cor

    ## *** parametrization and patterns
    structure <- .skeleton(structure = structure, data = data, indexData = outInit)

    ## *** covariance pattern
    structure <- .findUpatterns(structure, 
                                index.clusterTime = outInit$index.clusterTime, U.time = U.time,
                                index.cluster = outInit$index.cluster, U.cluster = U.cluster,
                                index.clusterStrata = outInit$index.clusterStrata, U.strata = U.strata)
    ## ** prepare calculation of the score
    if(precompute.moments && NCOL(X.mean)>0){
        if(is.na(var.weights[1])){
            wX.mean <- X.mean
            wY <- cbind(data[[var.outcome]])
        }else{
            wX.mean <- sweep(X.mean, FUN = "*", MARGIN = 1, STATS = sqrt(data[[var.weights[1]]]))
            wY <- cbind(data[[var.outcome]]*sqrt(data[[var.weights[1]]]))
        }

        precompute.XX <-  .precomputeXX(X = wX.mean, pattern = structure$Upattern$name, 
                                        pattern.ntime = stats::setNames(structure$Upattern$n.time, structure$Upattern$name),
                                        pattern.cluster = attr(structure$pattern,"list"), index.cluster = index.cluster)

        precompute.XY <-  .precomputeXR(X = precompute.XX$Xpattern, residuals = wY, pattern = structure$Upattern$name,
                                        pattern.ntime = stats::setNames(structure$Upattern$n.time, structure$Upattern$name),
                                        pattern.cluster = attr(structure$pattern,"list"), index.cluster = index.cluster)
    }else{
        precompute.XX <- NULL
        precompute.XY <- NULL
    }

    ## ** find all pairs of coefficients
    structure$pair.vcov <- stats::setNames(lapply(structure$Upattern$name, function(iName){## iName <- structure$Upattern$name[1]
        iParamVcov <- structure$Upattern[structure$Upattern$name == iName,"param"][[1]]
        if(length(iParamVcov)==0){return(NULL)}

        iOut <- unorderedPairs(iParamVcov)
        attr(iOut, "key") <- matrix(NA, nrow = length(iParamVcov), ncol = length(iParamVcov), dimnames = list(iParamVcov,iParamVcov))
        for(iCol in 1:NCOL(iOut)){
            attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
            attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol]] <- iCol
        }
        return(iOut)
    }),structure$Upattern$name)

    structure$pair.meanvcov <- stats::setNames(lapply(structure$Upattern$name, function(iName){ ## iName <- structure$Upattern$name[1]
        iPattern <- structure$Upattern[structure$Upattern$name == iName,,drop=FALSE]
        if(length(iPattern$param[[1]])==0){return(NULL)}
        iParamMu <- c(skeleton.mu$name[which(skeleton.mu$index.strata == iPattern$index.strata)],skeleton.mu$name[is.na(skeleton.mu$index.strata)>0])
        iOut <- unname(t(expand.grid(iParamMu, iPattern$param[[1]])))
        return(iOut)
    }), structure$Upattern$name)

    ## ** param
    skeleton.param <- rbind(skeleton.mu,                            
                            structure$param[is.na(structure$param$constraint),,drop=FALSE] ## only keep free parameters
                            )
    rownames(skeleton.param) <- NULL

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
                param = skeleton.param,
                drop.X = drop.X
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
##' @description remove un-identifiable columns from the design matrix 
##' @noRd
.model.matrix_regularize <- function(formula, data, augmodel = FALSE, type, drop.X){

    ## ## ** test 0: remove variable(s) with single level in the formula
    test.1value <- sapply(all.vars(formula),function(iVar){length(unique(data[[iVar]]))})
    if(any(test.1value==1)){
        ## no warning because this is normal behavior when stratifying
        formula <- stats::update(formula, stats::as.formula(paste0("~.-",paste(names(test.1value)[test.1value==1],collapse="-"))))
    }

    ## ** identify if there is an identifiability problem
    X <- stats::model.matrix(formula, data)
    attr(X,"variable") <- all.vars(formula)
    X.qr <- qr(X)

    if(X.qr$rank==NCOL(X)){
        if(NCOL(X)==0){
            return(X)
        }else if(augmodel){
            return(.augmodel.matrix(stats::delete.response(stats::terms(formula)),data))
        }else{
            return(X)
        }
    }else{
        if(drop.X==FALSE){
            stop("The design matrix for the ",type," structure does not have full rank according to the QR decomposition. \n", sep = "")
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
                    "The empty level(s) will be remove internally for the ",type," structure.\n",
                    "Consider applying droplevels to avoid this warning. \n",
                    sep = "")
            factor.drop <- name.factor[test.factor==FALSE]
            for(iFactor in factor.drop){
                data[[iFactor]] <- droplevels(data[[iFactor]])
            }
        }
    }

    ## ** create "naive" design matrix
    X <- .augmodel.matrix(tt,data)

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

        txt <- paste0("Constant values in the design matrix for the ",type," structure.\n")
        if(length(unique(rmX))==1){
            txt <- paste0(txt, "Coefficient \"",paste(unique(rmX), collapse = "\" \""),"\" relative to interaction \"",paste(names(ls.rmX), collapse = "\" \""),"\" has been removed. \n")
        }else{
            txt <- paste0(txt, "Coefficients \"",paste(unique(rmX), collapse = "\" \""),"\" relative to interactions \"",paste(names(ls.rmX), collapse = "\" \""),"\" have been removed. \n")
        }
        if(qr(X)$rank==X.qr$rank){
            message(txt)
        }else{
            warning(txt)
        }
        attr(X,"assign") <- attr(X.old,"assign")[test.keep] ## as.numeric(as.factor(attr(X.old,"assign")[test.keep])) - "(Intercept)" %in% colnames(X)
        attr(X,"contrasts") <- attr(X.old,"contrasts")
        attr(X,"variable") <- attr(X.old,"variable")
        if(augmodel || X.qr$rank!=NCOL(X.qr$qr)){
            attr(X,"formula") <- attr(X.old,"formula")
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
            attr(X,"variable") <- keep.attr$variable
            if(augmodel){
                attr(X,"formula") <- keep.attr$formula
                attr(X,"term.labels") <- keep.attr$term.labels[iIndex]
                attr(X,"order") <- keep.attr$order[iIndex]
                attr(X,"M.level") <- keep.attr$M.level[iIndex,,drop=FALSE]
                attr(X,"reference.level") <- keep.attr$reference.level
            }
        }
        if(NCOL(X)==X.qr$rank){break}
    }
    if(length(rmX)>0){
        message("Design matrix for the ",type," structure is singular. \n",
                "Coefficient",if(length(rmX)>1){"s"}," \"",paste(rmX, collapse = "\" \""),"\" has been removed. \n")
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
##' @description add more information about which variable contribute to which column in the design matrix
##' @noRd
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
    ## order(data[,all.vars(formula)])
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
    for(iCol in 1:p){ ## iCol <- 1

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
                iLs.allnameInteraction <- nlme::collapse(stats::setNames(lapply(iVar,function(iName){paste0(iName,iLs.grid[[iName]])}), iVar), sep = ":", as.factor = TRUE)
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

    var.time <- structure$name$time
    var.cluster <- structure$name$cluster
    var.strata <- structure$name$strata

    ## *** find unique levels
    U.cluster <- sort(unique(data[[var.cluster]]))
    U.time <- sort(unique(data[[var.time]]))
    if(!is.na(var.strata)){
        U.strata <- sort(unique(data[[var.strata]]))
    }else{
        U.strata <- 1
    }

    ## *** find position of each cluster
    index.cluster <- unname(split(1:NROW(data),data[[var.cluster]]))

    ## *** find time corresponding to each cluster
    index.clusterTime <- unname(split(as.numeric(factor(data[[var.time]],levels = U.time)),data[[var.cluster]]))

    ## *** re-order according to time
    order.clusterTime <- lapply(index.clusterTime, order)
    index.cluster <- mapply(x = index.cluster, y = order.clusterTime, function(x,y){x[y]}, SIMPLIFY = FALSE)
    attr(index.cluster, "vectorwise") <- as.numeric(factor(data[[var.cluster]], levels = U.cluster))
    index.clusterTime <- mapply(x = index.clusterTime, y = order.clusterTime, function(x,y){x[y]}, SIMPLIFY = FALSE)
    attr(index.clusterTime, "vectorwise") <- as.numeric(factor(data[[var.time]], levels = U.time))

    ## *** find strata corresponding to each cluster
    if(!is.na(var.strata)){
        iMatch <- tapply(as.character(data[[var.strata]]),data[[var.cluster]],unique)
        index.clusterStrata <- match(iMatch, U.strata)
    }else{
        index.clusterStrata <- rep(1, length(U.cluster))
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

## ** .updatefactor
##' @description Update factor variables with previously defined levels
##' @noRd
.updateFactor <- function(data, xfactor){

    if(length(xfactor)==0){return(data)}
    
    ## convert to factor with the right levels
    ff.factor <- names(xfactor)
    if(length(ff.factor)>0){
        for(iVar in ff.factor){ ## iVar <- ff.factor[1]
            iLevel <- xfactor[[iVar]]
            if(any(data[[iVar]] %in% iLevel == FALSE)){
                Wf <- setdiff(unique(data[[iVar]]), iLevel)
                stop("Unknown factor(s) \"",paste0(Wf,collapse="\" \""),"\" for variable \"",iVar,"\".\n",
                     "Valid factors: \"",paste0(iLevel, collapse="\" \""),"\".\n")
            }
            data[[iVar]] <- factor(data[[iVar]], levels = iLevel)
        }
    }

    return(data)
}

##----------------------------------------------------------------------
### model.matrix.R ends here
