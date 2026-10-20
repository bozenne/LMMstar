### model.matrix.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (21:50) 
## Version: 
## Last-Updated: okt 20 2026 (16:37) 
##           By: Brice Ozenne
##     Update #: 4155
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.matrix.lmm (documentation)
##' @title Design Matrix for a Linear Mixed Model
##' @description Extract or construct design matrices for a linear mixed model.
##'
##' @param object an \code{lmm} object
##' @param newdata [data.frame] dataset relative to which the design matrix should be constructed.
##' @param effects [character] design matrix relative to the mean model (\code{"mean"}), variance model (\code{"variance"}), correlation model (\code{"correlation"}),
##' or all the previous (\code{"all"}).
##' Can also be \code{"index"} to only output the normalize data and the cluster, time, strata indexes.
##' @param simplify [logical] simplify the data format of the output (matrix instead of a list of matrix) when possible.
##' @param drop.X [logical] when the design matrix does not have full rank, should columns be dropped? 
##' @param na.rm [logical] Should row containing missing values for the variables used in the linear mixed model be removed?
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return When \code{simplify} is \code{FALSE}, a list with the followin elements: \itemize{
##' \item \code{mean}: design matrix for the mean model
##' \item \code{Y}: vector of outcome values
##' \item \code{vcov}: list of elements for the variance and correlation models.
##' \item \code{index.cluster}: list containing, for each cluster, the location of the corresponding observations in the processed dataset.
##' \item \code{index.clusterTime}: list containing, for each cluster, the repetition index corresponding observations.
##' \item \code{index.clusterStrata}: list containing, for each cluster, the strata index corresponding observations.
##' \item \code{param}: data.frame describing the modle parameters.
##' \item \code{drop.X}: logical value indicating whether columns in the design matrix should be dropped if it has not full rank.
##' \item \code{precompute.XX}, \code{precompute.XY}: moments of X and Y 
##' }
##'
##' When \code{simplify} is \code{TRUE}, this list will be simplified into a list with three elements: \itemize{
##' \item \code{mean}: design matrix for the mean model
##' \item \code{variance}: design matrix for the variance model
##' \item \code{correlation}: design matrix for the correlation model
##' }
##' or a single design matrix.
##' 
##' @keywords methods

## * model.matrix.lmm (code)
##' @export
model.matrix.lmm <- function(object, newdata = NULL, effects = "mean", simplify = TRUE, drop.X = NULL, na.rm = FALSE, ...){

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** effects
    if(!is.character(effects) || !is.vector(effects)){
        stop("Argument \'effects\' must be a character vector. \n")
    }
    valid.effects <- c("index","mean","fixed","variance","correlation","all")
    if(any(effects %in% valid.effects == FALSE)){
        stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
             "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
    }
    if(all("index" %in% effects)){
        if(length(effects)>1){
            stop("Argument \'effects\' must have length 1 when containing the element \"index\". \n")
        }
    }else if(all("all" %in% effects)){
        if(length(effects)>1){
            stop("Argument \'effects\' must have length 1 when containing the element \"all\". \n")
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else{    
        effects[effects== "fixed"] <- "mean"
    }

    ## *** drop.X
    if(is.null(drop.X)){
        drop.X <- options$drop.X
    }
    
    ## *** simplify
    if(!is.numeric(simplify) && !is.logical(simplify)){
        stop("Argument \'simplify\' must be numeric or logical. \n")
    }
    if(length(simplify)!=1){
        stop("Argument \'simplify\' must have length 1. \n")
    }
    if(simplify %in% c(0,0.5,1) == FALSE){
        stop("Argument \'simplify\' must be TRUE/1, 0.5, or FALSE/0. \n")
    }

    ## ** update design matrix with new dataset
    if(is.null(newdata) && (is.null(object$index.na) || na.rm)){

        design <- object$design

    }else{

        ## *** normalized data
        if(is.null(newdata)){
            newdata <- as.data.frame(object$data)
        }else{
            newdata <- as.data.frame(newdata)
        }

        ## *** detect missing values
        var.manifest <- stats::variable.names(object)
        var.cluster <- attr(object$cluster$var,"original")
        var.time <- attr(object$time$var,"original")
        var.manifest.newdata <- intersect(names(newdata),var.manifest)

        if(na.rm == TRUE){
            na.rm <- any(is.na(newdata[var.manifest.newdata]))
        }

        ## *** prepare output
        design <- object$design
        design[setdiff(names(design),c("vcov","param","drop.X"))] <- list(NULL)
        design$vcov[c("var","cor","pattern","Upattern")] <- NULL

        ## *** index and exclude missing values
        var.strata <- attr(object$strata$var,"original")
        var.outcome <- object$outcome$var
        structure.var <- stats::na.omit(c(var.cluster, var.time, var.strata))        
        keep.col_testNA <- intersect(names(newdata), var.manifest)

        if(any(c("index","variance","correlation") %in% effects) || (!simplify && all(structure.var %in% names(newdata))) || (na.rm && all(structure.var %in% names(newdata)))){

            ## check dataset
            if(any(structure.var %in% names(newdata) == FALSE)){
                stop("Incorrect argument \'newdata\': missing variable(s) for the data structure \"",paste(structure.var[structure.var %in% names(newdata) == FALSE], collapse = "\" \""),"\".\n")
            }
            if(any(!is.na(var.time)) && any(var.time %in% keep.col_testNA)){
                test <- lapply(intersect(keep.col_testNA, var.time), function(iTime){ ## iTime <- "time"
                    if(any(newdata[[iTime]] %in% attr(object$time$levels,"original")[[iTime]] == FALSE)){
                        stop("Incorrect argument \'newdata\': invalid time variable ",iTime,". \n",
                             "Proposed value: \"",paste(setdiff(newdata[[iTime]],attr(object$time$levels,"original")[[iTime]]), collapse ="\" \""),"\"\n",
                             "Valid values: \"",paste(attr(object$time$levels,"original")[[iTime]], collapse ="\" \""),"\"\n")
                    }
                })
                
            }

            ## generate cluster/time/strata variables
            newdata.norm <- .lmmNormalizeData(data = newdata[keep.col_testNA],
                                              var.outcome = NA,
                                              var.cluster = var.cluster,
                                              var.time = var.time,
                                              var.strata = var.strata,
                                              droplevels = list(time = object$time$levels,
                                                                strata = object$strata$levels),
                                              initialize.cluster = object$design$vcov$ranef$crossed,
                                              initialize.time = setdiff(object$design$vcov$ranef$vars, object$design$vcov$ranef$var.cluster),
                                              na.rm = na.rm)

            ## extract indexes
            outInit <- .extractIndexData(data = newdata.norm$data, structure = design$vcov)
            if(identical(effects,"index")){
                if(simplify){
                    return(outInit)
                }else{
                    return(c(outInit, newdata.norm))
                }
            }
            design$index.cluster <- outInit$index.cluster
            design$index.clusterStrata <- outInit$index.clusterStrata
            design$index.clusterTime <- outInit$index.clusterTime
            design$index.na <- newdata.norm$index.na
            newdata <- newdata.norm$data 
                        
        }else if(na.rm){
            ## same as before but does not require all cluster/time/strata variables as input

            if(all(is.na(var.cluster)) || all(var.cluster %in% names(newdata) == FALSE)){
                var.cluster.newdata <- NA
            }else{
                var.cluster.newdata <- intersect(var.cluster, names(newdata))
            }
            if(all(is.na(var.time)) || all(var.time %in% names(newdata) == FALSE)){
                var.time.newdata <- NA
            }else{
                var.time.newdata <- intersect(var.time, names(newdata))
            }
            if(all(is.na(var.strata)) || all(var.strata %in% names(newdata) == FALSE)){
                var.strata.newdata <- NA
            }else{
                var.strata.newdata <- intersect(var.strata, names(newdata))
            }

            data.Nindex <- .lmmNormalizeData(data = newdata[keep.col_testNA],
                                             var.outcome = NA,
                                             var.cluster = var.cluster.newdata,
                                             var.time = var.time.newdata,
                                             var.strata = var.strata.newdata,
                                             droplevels = TRUE,
                                             initialize.cluster = object$design$vcov$ranef$crossed,
                                             initialize.time = NULL,
                                             na.rm = TRUE)
            design$index.na <- data.Nindex$index.na
            newdata <- data.Nindex$data
        }

        ## *** outcome
        if(!simplify && var.outcome %in% names(newdata)){
            design$Y <- newdata[[var.outcome]]
            ## indexY.NNA <- which(!is.na(newdata[[var.outcome]]))
            ## design$Y[indexY.NNA] <- stats::model.response(stats::model.frame(object$formula$mean.outcome, newdata[indexY.NNA,var.outcome,drop=FALSE]))
        }

        ## *** weights
        if(!simplify){
            if(!is.na(object$weight$var["likelihood"]) && object$weight$var["likelihood"] %in% names(newdata)){
                design$weight.likelihood <- newdata[[object$weight$var["likelihood"]]]
            }
            if(!is.na(object$weight$var["Omega"]) && object$weight$var["Omega"] %in% names(newdata)){
                design$weight.Omega <- newdata[[object$weight$var["Omega"]]]
            }
        }
        
        ## *** mean
        if("mean" %in% effects){
            
            ## check dataset
            ff.allvars <- all.vars(object$formula$mean.design)
            if(any(ff.allvars %in% names(newdata) == FALSE)){
                stop("Incorrect argument \'newdata\': missing variable(s) for the mean structure \"",paste(ff.allvars[ff.allvars %in% names(newdata) == FALSE], collapse = "\" \""),"\".\n")
            }

            ## convert to factor with the right levels
            data.mean <- .updateFactor(newdata, xfactor = object$xfactor$mean)

            ## update formula with attributes from the design matrix
            
            data.mean$XXindexXX <- 1 ## add latent variables used when terms where defined
            data.mean$XXtimeXX <- 1
            data.mean$XXclusterXX <- 1
            data.mean$XXstrataXX <- 1

            ## use stats::model.frame, attr(object$design$mean,"terms") to handle spline (add attributes predvars with the position of the knots)
            data.mf.mean <- stats::model.frame(attr(object$design$mean,"terms"),
                                               data = data.mean, na.action = stats::na.pass)
            ## use object$formula$mean.design to handle interactions
            ## e.g. B:C + A:B:C does lead to the right names (otherwise can lead to names like B:C:A)
            design$mean  <- stats::model.matrix(object$formula$mean.design, data.mf.mean)[,colnames(object$design$mean),drop=FALSE]
        }

        ## *** variance-covariance
        if("variance" %in% effects){
            ## check dataset
            ff.allvars <- c(all.vars(object$formula$var),all.vars(object$formula$cor))
            if(any(ff.allvars %in% names(newdata) == FALSE)){
                stop("Incorrect argument \'newdata\': missing variable(s) for the variance-covariance structure \"",paste(ff.allvars[ff.allvars %in% names(newdata) == FALSE], collapse = "\" \""),"\".\n")
            }
            ## update levels
            design$vcov$var <- list(lp2X = object$design$vcov$var$lp2X,
                                    lp2data = object$design$vcov$var$lp2data,
                                    X.attr = attributes(object$design$vcov$var$X),
                                    xfactor = object$xfactor$var)
            design$vcov$cor <- list(lp2X = object$design$vcov$cor$lp2X,
                                    lp2data = object$design$vcov$cor$lp2data,
                                    X.attr = attributes(object$design$vcov$cor$X),
                                    xfactor = object$xfactor$cor)

            ## form design matrix
            outDesign <- .vcov.matrix.lmm(structure = design$vcov, data = newdata, index.cluster = outInit$index.cluster, drop.X = drop.X, sep = options$sep["lp"])
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
                                          index.clusterStrata = outInit$index.clusterStrata, U.strata = outInit$U.strata,
                                          sep = options$sep["pattern"])
        }
    }

    ## ** export
    if(simplify>0){
        design2 <- list(mean = design$mean,
                        variance = design$vcov$var$X,
                        correlation = design$vcov$cor$X)
        if(simplify >= 1){
            attributes(design2$mean) <- attributes(design$mean)[c("dim","dimnames","assign","contrasts")]
            attributes(design2$variance) <- attributes(design$vcov$var$X)[c("dim","dimnames","assign","contrasts")]
            attributes(design2$correlation) <- attributes(design$vcov$cor$X)[c("dim","dimnames","assign","contrasts")]
        }
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
.vcov.matrix.lmm <- function(structure, data, index.cluster, drop.X, sep){

    var.cluster <- structure$name$cluster
    var.time <- structure$name$time
    var.strata <- structure$name$strata
    U.cluster <- sort(unique(data[[var.cluster]]))
    
    structure.CUSTOM <- inherits(structure,"CUSTOM")
    structure.timeFormat <- attr(time(structure),"format")

    grid.timeFormat <- do.call(rbind,lapply(names(structure.timeFormat), function(iName){data.frame(moment = iName, format = structure.timeFormat[[iName]])}))
    rownames(grid.timeFormat) <- grid.timeFormat$moment
    grid.timeFormat$moment <- factor(grid.timeFormat$moment, levels = c("variance","correlation","correlation.cross"), labels = c("variance","correlation","correlation"))
    out <- list(variance = list(X = NULL, lp = NULL, lp2data = NULL, lp2X = NULL, xfactor = NULL, xnumeric = NULL))

    if(NROW(grid.timeFormat) == 1){ 
        ## do nothing 
    }else if(NROW(grid.timeFormat) == 2){
        out$correlation <-  list(X = NULL, lp = NULL, lp2data = NULL, lp2X = NULL, xfactor = NULL, xnumeric = NULL)
    }else if(NROW(grid.timeFormat)==3){
        
        if(grid.timeFormat["correlation","format"]==grid.timeFormat["correlation.cross","format"] & isTRUE(all.equal(structure$formula[["correlation"]],structure$formula[["correlation.cross"]]))){
            out$correlation <-  list(X = NULL, lp = NULL, lp2data = NULL, lp2X = NULL, xfactor = NULL, xnumeric = NULL)
            grid.timeFormat <- grid.timeFormat[c("variance","correlation"),] ## no need for specific design matrix for the cross blocks
        }else{
            out$correlation <-  list(X = NULL, lp = NULL, lp2data = NULL, lp2X = NULL, xfactor = NULL, xnumeric = NULL,
                                     X.cross = NULL, lp.cross = NULL, lp2data.cross = NULL, lp2X.cross = NULL, xfactor.cross = NULL, xnumeric.cross = NULL)
        }
    }else{
        stop("Something went wrong with the identification of the format of the time variable in the covariance structure. \n")
    }
    test.newdata <- !is.null(structure$param) ## if TRUE then model.matrix call on newdata otherwise call from lmm->.model.matrix on original data

    ## ** variance and correlation design lists
    for(iGrid in 1:NROW(grid.timeFormat)){ ## iGrid <- 1
        iMoment <- as.character(grid.timeFormat[iGrid,"moment"])
        iFormat <- grid.timeFormat[iGrid,"format"]
        iSuffix <- ifelse(iGrid==3,".cross","")
        iFormula <- structure$formula[[paste0(iMoment,iSuffix)]]
        iVars.strata <- stats::na.omit(structure$name$strata)
        iVars <- setdiff(all.vars(iFormula), iVars.strata)

        ## *** normalize dataset (e.g. set right levels to factor variables)
        if("xfactor" %in% names(structure[[iMoment]])){ ## from model.matrix.lmm
            iXfactor <- structure[[iMoment]][[paste0("xfactor",iSuffix)]]
            iXnumeric <- structure[[iMoment]][[paste0("xnumeric",iSuffix)]]
        }else if(iFormat == "original" || length(iVars)==0){
            iXfactor <- NULL
            iXnumeric <- NULL
        }else if(iFormat == "numeric"){ ## which variables to convert to numeric
            iXfactor <- NULL
            iXnumeric <- stats::na.omit(union(var.time,iVars))                        
        }else if(iFormat == "factor"){
            iXnumeric <- NULL
            iTestFac <- sapply(data[stats::na.omit(union(var.time,iVars))],is.factor)
            if(sum(iTestFac==FALSE)==0){
                iXfactor <- stats::.getXlevels(stats::terms(iFormula),data)
            }else{
                data.md <- stats::model.frame(iFormula,data)
                for(iiVar in names(iTestFac)[which(iTestFac==FALSE)]){
                    data.md[[iiVar]] <- as.factor(data.md[[iiVar]])
                }
                iXfactor <- stats::.getXlevels(stats::terms(iFormula),data.md)
            }
        }else if(iFormat == "factor0"){
            iXnumeric <- NULL
            if(is.na(var.time)){
                iXfactor <- NULL
            }else if(is.factor(data[[var.time]])){
                iXfactor <- stats::.getXlevels(stats::terms(iFormula),data)
            }else{
                data.md <- stats::setNames(data.frame(as.factor(data.md[[var.time]])), var.time)
                iXfactor <- stats::.getXlevels(stats::terms(iFormula), data.md)
            }            
        }

        if(length(iVars.strata)>0){
            iTestFac.strata <- sapply(data[iVars.strata],is.factor)
            iFormula.strata <- stats::reformulate(termlabels = iVars.strata)
            if(sum(iTestFac.strata==FALSE)==0){
                iXfactor <- c(iXfactor,stats::.getXlevels(stats::terms(iFormula.strata),data))
            }else{
                data.md.strata <- stats::model.frame(iFormula.strata,data)
                for(iiVar in names(iTestFac.strata)[which(iTestFac.strata==FALSE)]){
                    data.md.strata[[iiVar]] <- as.factor(data.md.strata[[iiVar]])
                }
                iXfactor <- c(iXfactor,stats::.getXlevels(stats::terms(iFormula.strata),data.md.strata))
            }
        }
        iData <- .updateFactor(data, xfactor = iXfactor, xnumeric = iXnumeric)

        ## *** design matrix
        if(structure.CUSTOM){
            iX <- iData[,iVars,drop=FALSE]                    
        }else if(!is.null(iFormula)){
            if(test.newdata == FALSE){ ## model fit (structure)
                iX <- .model.matrix_regularize(iFormula, data = iData, augmodel = TRUE, type = iMoment, drop.X = drop.X)                
                attr(iX, "original.colnames") <- colnames(iX)                
            }else{ ## newdata
                iAttr.X <- structure[[iMoment]]$X.attr
                iKeep.attr <- setdiff(names(iAttr.X), c("dim","dimnames"))
                iX <- stats::model.matrix(iFormula, data = iData)[,iAttr.X$original.colnames,drop=FALSE]
                colnames(iX) <- iAttr.X$dimnames[[2]]
                for(iAssign in iKeep.attr){
                attr(iX,iAssign) <- iAttr.X[[iAssign]]
                }
            }
        }
        ## *** linear predictors
        ## linear predictor for each observation
        iLp <- nlme::collapse(iX, sep = sep, as.factor = !test.newdata)
        
        if(!test.newdata){
            ## position of the observations with distinct linear predictors
            iIndex.Ulp <- which(!duplicated(iLp))
            ## convert linear predictor to numeric, possible updating the order of the level to match "natural ordering"
            if(NCOL(iX)>1){
                if(is.na(var.strata)){
                    iNewlevel <- orderLtoR(iX[iIndex.Ulp,,drop=FALSE], strata = NULL)
                }else if(structure.CUSTOM){
                    iNewOrder <- c(setdiff(colnames(iX),var.strata),var.strata)
                    iNewlevel <- orderLtoR(iX[iIndex.Ulp,iNewOrder,drop=FALSE], strata = NULL)
                }else{
                    iU.strata <- levels(iData[[var.strata]])
                    iNewlevel <- orderLtoR(iX[iIndex.Ulp,,drop=FALSE], strata = factor(attr(iX,"M.level")[[var.strata]], iU.strata))
                ##  iX[iIndex.Ulp[order(iNewlevel)],,drop=FALSE]
                }
                iLp.num <- as.numeric(factor(iLp, levels = iLp[iIndex.Ulp][order(iNewlevel)]))
            }else{
                iNewlevel <- levels(iLp)
                iLp.num <- as.numeric(iLp)
            }
            ## re-order position of the observations with distinct linear predictors to match increasing linear predictors
            if(length(iIndex.Ulp)>1){
                iIndex.Ulp <- iIndex.Ulp[order(iLp.num[iIndex.Ulp], decreasing = FALSE)]
            }
            ## re-express the linear predictor in term of design matrix or original data
            if(structure.CUSTOM){
                iLp2X <- NULL
                iLp2data <- iX[iIndex.Ulp,,drop=FALSE]
                rownames(iLp2data) <- iLp[iIndex.Ulp]
            }else{
                iLp2X <- iX[iIndex.Ulp,,drop=FALSE]
                rownames(iLp2X) <- iLp[iIndex.Ulp]
                iLp2data <- iData[iIndex.Ulp,attr(iX,"variable"),drop=FALSE]
                rownames(iLp2data) <- iLp[iIndex.Ulp]
            }
        }else{
            iLp.num <- as.numeric(factor(iLp, levels = rownames(structure[[iMoment]][[paste0("lp2X",iSuffix)]])))
            iLp2X <- structure[[iMoment]][[paste0("lp2X",iSuffix)]]
            iLp2data <- structure[[iMoment]][[paste0("lp2data",iSuffix)]]
        }

        ## *** pattern for each cluster
        iPattern <- as.numeric(as.factor(sapply(index.cluster, function(iC){paste(iLp.num[iC], collapse = ".")})))
        ## Note: tapply(out[[iMoment]]$lp, attr(index.cluster,"vectorwise"), paste, collapse = ".")
        ## does not work when the time are not ordered in the initial dataset
        if(!test.newdata){
            ## position of the cluster with distinct patterns
            iIndex.pattern <- which(!duplicated(iPattern))
            iIndex.pattern <- iIndex.pattern[order(iPattern[iIndex.pattern], decreasing = FALSE)] ## re-order
            ## re-express the pattern in term of linear predictors
            iPattern2lp <- unname(lapply(index.cluster[iIndex.pattern], function(iIndex){iLp.num[iIndex]}))
        }else{
            iPattern2lp <- structure[[iMoment]][[paste0("pattern2lp",iSuffix)]]
        }
        
        ## *** update
        out[[iMoment]][[paste0("xfactor",iSuffix)]] <- iXfactor
        out[[iMoment]][[paste0("xnumeric",iSuffix)]] <- iXnumeric
        out[[iMoment]][[paste0("X",iSuffix)]] <- iX 
        out[[iMoment]][[paste0("lp",iSuffix)]] <- iLp.num
        out[[iMoment]][[paste0("lp2X",iSuffix)]] <- iLp2X
        out[[iMoment]][[paste0("lp2data",iSuffix)]] <- iLp2data
        out[[iMoment]][[paste0("pattern",iSuffix)]] <- iPattern
        out[[iMoment]][[paste0("pattern2lp",iSuffix)]] <- iPattern2lp        
    }
    
    ## ** export
    return(out)
}

## * .model.matrix.lmm
.model.matrix.lmm <- function(formula.mean, structure,
                              data, var.outcome, var.weights,
                              drop.X, ## drop singular component of the design matrix
                              options, df){

    ## ** indexes
    outInit <- .extractIndexData(data = data, var.cluster = structure$name$cluster, var.time = structure$name$ordering, var.strata = structure$name$strata)
    U.cluster <- outInit$U.cluster
    U.time <- outInit$U.time
    U.strata <- outInit$U.strata
    n.strata <- length(U.strata)
    index.cluster <- outInit$index.cluster
    index.clusterStrata <- outInit$index.clusterStrata
    index.clusterTime <- outInit$index.clusterTime

    ## ** mean
    ## use stats::model.frame to handle splines
    data.mf <- stats::model.frame(formula.mean,data)
    outMean <- .model.matrix_regularize(formula.mean, data = data.mf, augmodel = TRUE, type = "mean", drop.X = drop.X)

    ## ** variance
    ## *** design matrix
    outDesign <- .vcov.matrix.lmm(structure = structure, data = data, index.cluster = index.cluster, drop.X = drop.X, sep = options$sep["lp"])
    structure$xfactor <- outDesign$xfactor
    structure$var <- outDesign$variance
    structure$cor <- outDesign$correlation

    ## *** parametrization and patterns
    structure <- .skeleton(structure = structure, data = data, indexData = outInit, options = options)

    ## *** covariance pattern
    structure <- .findUpatterns(structure, 
                                index.clusterTime = index.clusterTime, U.time = U.time,
                                index.cluster = index.cluster, U.cluster = U.cluster,
                                index.clusterStrata = index.clusterStrata, U.strata = U.strata,
                                sep = options$sep[c("rho.name","pattern")])
    n.pattern <- NROW(structure$Upattern)
    name.pattern <- structure$Upattern$name
    pattern.cluster <- structure$Upattern$index.cluster

    ## ** weights (cluster-specific)
    ## likelihood: use to rescale X and residuals when evaluating the moments in the log-likelihood, score, information, df
    ## Omega: use to rescale the residual variance covariance matrix

    ## *** prepare
    if(!is.na(var.weights["Omega"])){        
        test.unique <- tapply(data[[var.weights["Omega"]]], attr(index.cluster,"vectorwise"), FUN = function(iVec){sum(duplicated(iVec)==FALSE)})
        
        if(any(test.unique!=1)){

            rep.pattern <- stats::setNames(lapply(name.pattern, function(iPattern){ ## iPattern <- name.pattern[1]
                ## find unique repetition pattern
                iIndex.cluster <- index.cluster[pattern.cluster[[iPattern]]]
                iRep <- sapply(iIndex.cluster, function(iIndex){paste(data[[var.weights["Omega"]]][iIndex],collapse="_")})
                
                ## export
                factor(iRep, levels = unique(iRep))

            }), name.pattern)

            structure$pattern[unlist(pattern.cluster),"replicate"] <- unlist(lapply(rep.pattern, as.character))
            structure$pattern[unlist(pattern.cluster),"replicate.num"] <- unlist(lapply(rep.pattern, as.numeric))
            structure$pattern[unlist(pattern.cluster),"all"] <- paste(structure$pattern$vcov, structure$pattern$replicate.num,sep="_")

            multiplier <- stats::setNames(lapply(name.pattern, function(iPattern){ ## iPattern <- name.pattern[1]
                ## example cluster for each repetition pattern
                iIndex.cluster <- index.cluster[pattern.cluster[[iPattern]]]
                iRep <- rep.pattern[[iPattern]] 
                iIndex.cluster1 <- iIndex.cluster[which(duplicated(iRep)==FALSE)]
                
                ## number of repetitions per outcome for each repetition pattern
                iMweights <- do.call(rbind,lapply(iIndex.cluster1, function(iIndex){data[[var.weights["Omega"]]][iIndex]}))
                
                ## max number repetitions per pair of outcomes
                ## a) pmax approach
                iNtime <- structure$Upattern[structure$Upattern$name==iPattern,"n.time"]
                iOut <- array(NA, dim = c(iNtime,iNtime,NROW(iMweights)))
                for(iR in 1:NROW(iMweights)){
                    iOut[,,iR] <- pmax(matrix(iMweights[iR,], nrow = iNtime, ncol = iNtime, byrow = TRUE),
                                       matrix(iMweights[iR,], nrow = iNtime, ncol = iNtime, byrow = FALSE))
                }
                ## b) matrix product approach
                ## iMweights.sort <- do.call(rbind,apply(iMweights, MARGIN = 1, sort, simplify = FALSE)) ## same but sorted
                ## iTau <- lower.tri(matrix(0, nrow = iNtime, ncol = iNtime), diag = TRUE)
                ## iDelta <- apply(iMweights.sort - c(0,iMweights.sort[-iNtime]), MARGIN = 1, diag, simplify = FALSE)
                ## iP <- apply(cbind(iMweights,iMweights.sort), MARGIN = 1, FUN = function(iRow){
                ##     do.call(cbind,lapply(iRow[iNtime + (1:iNtime)], function(iW){(iW==iRow[(1:iNtime)])}))
                ## }, simplify = FALSE)
                ## lapply(1:length(levels(iRep)), function(iR){
                ##     (iP[[iR]] %*% iTau) %*% iDelta[[iR]] %*% t(iP[[iR]] %*% iTau)
                ## })
                return(iOut)
            }), name.pattern)

            structure$Upattern$n.replicate <- sapply(rep.pattern, function(iPattern){length(levels(iPattern))})            
        }
    }else{
        test.unique <- rep(1, length(index.cluster))
    }

    ## *** store
    index.cluster1 <- sapply(index.cluster,"[",1)        

    if(is.na(var.weights["likelihood"]) && (is.na(var.weights["Omega"]) || any(test.unique!=1))){

        weights <- list(XYR = NULL,
                        likelihood = rep(1, length(index.cluster)),
                        Omega = rep(1, length(index.cluster)))        

    }else if(!is.na(var.weights["likelihood"]) && (is.na(var.weights["Omega"]) || any(test.unique!=1))){

        weights <- list(XYR = data[[var.weights["likelihood"]]],
                        likelihood = data[[var.weights["likelihood"]]][index.cluster1],
                        Omega = rep(1, length(index.cluster)))

    }else if(is.na(var.weights["likelihood"]) && !is.na(var.weights["Omega"]) && all(test.unique==1)){ 

        weights <- list(XYR = data[[var.weights["Omega"]]],
                        likelihood = rep(1, length(index.cluster)),
                        Omega = data[[var.weights["Omega"]]][index.cluster1])
        
    }else if(!is.na(var.weights["likelihood"]) && !is.na(var.weights["Omega"]) && all(test.unique==1)){ 

        weights <- list(XYR = data[[var.weights["likelihood"]]] * data[[var.weights["Omega"]]],
                        likelihood = data[[var.weights["likelihood"]]][index.cluster1],
                        Omega = data[[var.weights["Omega"]]][index.cluster1])
    }
    if(any(test.unique!=1)){
        weights$multiplier <- multiplier
    }

    ## ** param
    browser()
    skeleton.param <- rbind(skeleton.mu,                            
                            structure$param[,names(skeleton.mu),drop=FALSE]
                            )    
    rownames(skeleton.param) <- NULL
    structure$param <- NULL

    order.type <- factor(skeleton.param$type, c("mu", "sigma", "k", "rho"))
    if(is.unsorted(order.type)){ ## in presence of strata, typically need to reorder by type
        skeleton.param <- skeleton.param[order(order.type),,drop=FALSE]
    }

    ## ** gather and export
    structure$pattern2param <- structure$Upattern[c("name","var","cor","index.strata","n.time")]
    structure$pattern2param$index.cluster1 <- sapply(structure$Upattern$index.cluster,"[",1)
    out <- list(param = skeleton.param,
                index.cluster = attr(index.cluster,"vectorwise"),
                index.clusterPattern = structure$pattern,
                index.clusterTime = attr(index.clusterTime,"vectorwise"),
                index.clusterStrata = index.clusterStrata,
                mean = X.mean, drop.X = drop.X,
                structure = structure[c("call","name","formula","twin","class","pattern2param")],
                variance = structure$var[c("Xpattern","pattern2lp","lp2data","lp2X")],
                correlation = structure$cor[c("Xpattern","pattern2lp","lp2data","lp2X")],                                
                weights = weights
                )

    return(out)
}

## * helpers
## ** .model.matrix_regularize
##' @description remove un-identifiable columns from the design matrix 
##' @noRd
.model.matrix_regularize <- function(formula, data, augmodel = FALSE, type, drop.X){

    ## ** initial check: is there is an identifiability problem?
    X <- stats::model.matrix(formula, data)
    X.qr <- qr(X)
    if(X.qr$rank==NCOL(X)){
        if(NCOL(X)==0){
            return(X)
        }else if(augmodel){
            return(.augmodel.matrix(formula,data))
        }else{
            return(X)
        }
    }else if(drop.X==FALSE){
        type <- switch(type, "var" = "variance", "cor" = "correlation", type)
        stop("The design matrix for the ",type," structure does not have full rank according to the QR decomposition. \n", sep = "")
    }
    regularize <- list()
    
    ## ** remedy 1: remove empty levels in factor or character/factor variable with a single level
    term.rhs <- stats::delete.response(stats::terms(formula))
    name.var <- all.vars(term.rhs)
    ## character ---> factor
    var.character <- name.var[sapply(data[name.var], is.character)]
    if(length(var.character)>0){
        for(iVar in name.character){
            data[name.character] <- as.factor(data[name.character])
        }
    }
    var.factor <- name.var[sapply(data[name.var], is.factor)]
    
    if(length(var.factor)>0){

        ## *** drop factor
        test.factor <- sapply(var.factor, function(iVar){all(levels(data[[iVar]]) %in% data[[iVar]])})
        if(any(test.factor==FALSE)){
            type <- switch(type, "var" = "variance", "cor" = "correlation", type)
            message("Factor variable(s) with empty level: \"",paste(name.factor[test.factor==FALSE], collapse = "\" \""),"\"\n ",
                    "The empty level(s) will be remove internally for the ",type," structure.\n",
                    "Consider applying droplevels to avoid this warning. \n",
                    sep = "")
            regularize$factor.drop <- var.factor[test.factor==FALSE]
            for(iFactor in regularize$factor.drop){
                data[[iFactor]] <- droplevels(data[[iFactor]])
            }
        }

        ## *** remove variable with single level
        test.1value <- sapply(var.factor, function(iVar){length(levels(data[[iVar]]))})
        if(any(test.1value==1)){
            ## no warning because this is normal behavior when stratifying
            regularize$var.drop <- names(test.1value)[test.1value==1]
            formula <- stats::update(formula, stats::as.formula(paste0("~.-",paste(regularize$var.drop,collapse="-"))))
        }
    }

    ## ** second check: is there still is an identifiability problem?
    X <- stats::model.matrix(formula, data)
    X.qr <- qr(X)
    if(X.qr$rank==NCOL(X)){
        if(augmodel){
            X <- .augmodel.matrix(formula, data)            
        }
        attr(X,"regularize") <- regularize
        return(X)        
    }else{ ## get information about the full design matrix
        outAug <- .augmodel.matrix(term.rhs, data)
    }

    ## ** intermediate check: linear dependencies order 1
    if(any(outAug$param$order==1)){
        X1 <- X[,rownames(outAug$param)[outAug$param$order==1],drop=FALSE]
        qr.order1 <- qr(X1)
        if(qr.order1$rank < NCOL(qr.order1$qr)){
            colnames(X1) <- paste0("X",1:NCOL(X1))
            lmDiag <- stats::lm(stats::reformulate(response = "Y", termlabels = colnames(X1), intercept = TRUE), data = data.frame(Y = 1:NROW(X1), X1))            
            stop("The design matrix for the ",type," structure does not have full rank according to the QR decomposition. \n", sep = "",
                 "Some variables are linearly dependent: \n", paste(capture.output(print(alias(lmDiag)$Complete)), collapse = "\n"),"\n")
        }
    }

    var.order1 <- unlist(outAug$param[outAug$param$order==1,"variable"])
    if(any(name.var %in% var.order1 ==  FALSE)){
        X1bis <- stats::model.matrix(reformulate(termlabels = name.var), data = data)
        qr.order1bis <- qr(X1bis)
        if(qr.order1bis$rank < NCOL(qr.order1bis$qr)){
            colnames(X1bis) <- paste0("X",1:NCOL(X1bis))
            lmDiag <- stats::lm(stats::reformulate(response = "Y", termlabels = colnames(X1bis), intercept = TRUE), data = data.frame(Y = 1:NROW(X1bis), X1bis))
            warning("Some variables are linearly dependent: \n", paste(capture.output(print(alias(lmDiag)$Complete)), collapse = "\n"),"\n")
        }
    }

    ## ** remedy 2: linear dependencies order 2 (or more)
    ## *** remove empty levels in factor or character/factor variable with a single level
    var.numeric <- name.var[sapply(data[name.var], is.numeric)]

    ls.rmX <- stats::setNames(lapply(which(outAug$param$order>1), function(iRow){ ## iRow <- 6
        ## variables involved in the interactions
        iVar <- outAug$param[iRow,"variable"][[1]]
        ## identify coefficient relative to this interaction 
        iNoVar <- setdiff(name.var,iVar)
        if(length(iNoVar)>0){ ## covariate combination relative to reference level for the 'other' covariates
            iSubset <- rowSums(outAug$param[,iNoVar,drop=FALSE] != do.call(rbind,replicate(outAug$param[which(outAug$param$reference),iNoVar,drop=FALSE], n = NROW(outAug$param))))
        }else{
            iSubset <- rep(0,NROW(outAug$param))
        }
        iM.level <- outAug$param[which(iSubset==0),,drop=FALSE]

        if(all(iVar %in% var.numeric == FALSE)){ ## all factor
            iTest.constant <- table(interaction(data[iVar]))
            iOut <- rownames(iM.level)[iTest.constant==0]
        }else if(all(iVar %in% var.numeric)){ ## all numeric
            ## no check needed
        }else{
            ## all combination of factor/character variables
            ## use as.numeric, as.factor, as.character to be robust against different ordering of the levels and special characters
            iM.level$group <- interaction(lapply(iM.level[,setdiff(iVar,var.numeric),drop=FALSE], FUN = function(iCol){as.numeric(as.factor(as.character(iCol)))}))
            ## find the term with the largest order 
            iM.level$group.highestOrder <- unlist(tapply(iM.level$order, iM.level$group, function(iOrder){iOrder == max(iOrder)}, simplify = FALSE))[order(order(iM.level$group))]
            ## split the data into groups according to the combination of factor/character variables
            iGroup <- factor(interaction(lapply(data[setdiff(iVar,var.numeric)], FUN = function(iCol){as.numeric(as.factor(as.character(iCol)))}), drop = TRUE),
                             levels = levels(iM.level$group))
            ## are all the numeric variables constant within groups?
            iTest.constant <- do.call(cbind,by(data[intersect(iVar,var.numeric)], iGroup, function(iDF){apply(iDF, MARGIN = 2, FUN = sd)}, simplify = FALSE))
            iGroup.constant <- which(colSums(iTest.constant==0)==(length(intersect(iVar,var.numeric))))
            ##
            iOut <- rownames(iM.level)[(iM.level$group %in% iGroup.constant) & iM.level$group.highestOrder]
        }
        return(iOut)
    }), outAug$param$term[outAug$param$order>1])

    regularize$drop.X <- unlist(unique(ls.rmX))
    if(length(regularize$drop.X)>0){
        attrAssign <- attr(outAug$X,"assign")[which(colnames(X) != regularize$drop.X)]
        attrContrast <- attr(outAug$X,"contrasts")[intersect(names(attr(outAug$X,"contrasts")),unlist(outAug$param[which(colnames(X) != regularize$drop.X),"variable"]))]
        outAug$X <- outAug$X[,which(colnames(X) != regularize$drop.X),drop=FALSE]
        attr(outAug$X,"assign") <- attrAssign
        attr(outAug$X,"contrasts") <- attrContrast
    }

    ## ** export
    X.qr <- qr(outAug$X)
    if(X.qr$rank==NCOL(outAug$X)){
        return(c(outAug, regularize))
    }else if(drop.X==FALSE){
        type <- switch(type, "var" = "variance", "cor" = "correlation", type)
        stop("The design matrix for the ",type," structure does not have full rank according to the QR decomposition. \n", sep = "")
    }    
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
                txt.error <- paste0("The column \"",iVar,"\" of argument \'data\' involves contrasts with multiple levels. \n",
                                    "This can be seen when running stats::contrasts(data$",iVar,"). \n")
                if(inherits(data[[iVar]],"ordered")){
                    stop(txt.error,
                         "Could be because \"",iVar,"\" is an ordered factor instead of a (usual) factor. \n",
                         "Considering applying as.factor to the variable before running lmm. \n")
                }else if(any(options()$contrast!=c("contr.treatment","contr.poly"))){
                    stop(txt.error,
                         "Could be because options()$contrast has been modified to non-standard contrasts. \n")
                }else{
                    stop(txt.error)
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
                iLs.grid <- expand.grid(iLs.factor[lengths(iLs.factor)>0])
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

    ## ** export
    M.level <- cbind(term = X.term,
                     order = X.order,
                     reference = NA,
                     variable = NA,
                     constraint = NA,
                     do.call(rbind,X.level2))
    M.level$reference <- rowSums(M.level[names(reference2)] != do.call(rbind, replicate(reference2, n = NROW(M.level), simplify = FALSE)))==0
    M.level$variable <- lapply(X.level, function(iDf){setdiff(names(iDf),"reference")})

    ## reference2
    return(list(X = X,
                param = M.level))
}

## ** .extractIndexData
.extractIndexData <- function(data, var.cluster, var.time, var.strata){

    ## *** find unique levels
    U.cluster <- sort(unique(data[[var.cluster]]))
    U.time <- sort(unique(data[[var.time]]))
    if(!is.na(var.strata)){
        U.strata <- sort(unique(data[[var.strata]]))
    }else{
        U.strata <- 1
    }

    ## *** find position of each cluster and each time
    data$XXindexXX <- 1:NROW(data)
    data.order <- data[order(data[[var.time]]),c("XXindexXX",var.cluster,var.time)]
    index.cluster <- unname(tapply(data.order$XXindexXX,data.order[[var.cluster]],base::identity, simplify = FALSE))
    attr(index.cluster, "vectorwise") <- as.numeric(factor(data[[var.cluster]], levels = U.cluster))
    index.clusterTime <- unname(tapply(data.order[[var.time]],data.order[[var.cluster]],base::identity, simplify = FALSE))
    attr(index.clusterTime, "vectorwise") <- as.numeric(factor(data[[var.time]], levels = U.time))

    ## index.cluster <- unname(split(1:NROW(data),data[[var.cluster]]))
    ## index.clusterTime <- unname(split(as.numeric(factor(data[[var.time]],levels = U.time)),data[[var.cluster]]))
    ## order.clusterTime <- lapply(index.clusterTime, order)
    ## index.cluster <- mapply(x = index.cluster, y = order.clusterTime, function(x,y){x[y]}, SIMPLIFY = FALSE)
    ## index.clusterTime <- mapply(x = index.clusterTime, y = order.clusterTime, function(x,y){x[y]}, SIMPLIFY = FALSE)

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
##' @description Update factor variables with previously defined levels and convert other variables to numeric
##' @param data [data.frame]
##' @param xfactor [list] each element of the list contain the levels of the variable to be converted to factor.
##' @param xnumeric [character vector] name of the variable(s) to be converted to numeric format.
##' @noRd
##'
##' @examples
##' df <- data.frame(age = 1:10, ageGroup = 1:10, sexInd = as.factor(c("F","M")))
##' str(.updateFactor(df, xfactor = list(ageGroup=1:10), xnumeric = "sexInd"))
.updateFactor <- function(data, xfactor = NULL, xnumeric = NULL){

    ## *** special case where nothing needs to be done
    xfactor <- na.omit(xfactor)
    xnumeric <- na.omit(xnumeric)
    if(length(xfactor)==0 && length(xnumeric)==0){
        return(data)
    }

    ## *** convert to factor
    ff.factor <- names(xfactor)    
    if(length(ff.factor)>0){
        for(iVar in ff.factor){ ## iVar <- ff.factor[1]
            iLevel <- xfactor[[iVar]]
            if(any(data[[iVar]] %in% c(iLevel,NA) == FALSE)){
                Wf <- setdiff(unique(data[[iVar]]), c(iLevel,NA))
                stop("Unknown factor(s) \"",paste0(Wf,collapse="\" \""),"\" for variable \"",iVar,"\".\n",
                     "Valid factors: \"",paste0(iLevel, collapse="\" \""),"\".\n")
            }
            data[[iVar]] <- factor(data[[iVar]], levels = iLevel)
        }
    }
    ## *** convert to numeric
    if(length(xnumeric)>0){
        for(iVar in xnumeric){ ## iVar <- ff.factor[1]
            if(is.logical(data[[iVar]])){ 
                data[[iVar]] <- as.numeric(data[[iVar]]) + 1
            }else if(!is.numeric(data[[iVar]])){ 
                data[[iVar]] <- as.numeric(as.factor(data[[iVar]]))
            }##else if(is.numeric(data[[iVar]])){ 
            ##    data[[iVar]] <- data[[iVar]] - min(data[[iVar]]) + 1
            ##}
        }
    }

    ## *** export
    return(data)
}


##----------------------------------------------------------------------
### model.matrix.R ends here
