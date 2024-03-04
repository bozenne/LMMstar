### model.frame.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (14:57) 
## Version: 
## Last-Updated: mar  4 2024 (19:21) 
##           By: Brice Ozenne
##     Update #: 94
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.frame.lmm (documentation)
##' @title Extracting the Model Frame from a Linear Mixed Model
##' @description Variables needed to fit the Linear Mixed Model.
##' @param formula [lmm] linear mixed model object
##' @param data [data.frame]
##' @param type [character] By default returns the processed dataset used to fit the Linear Mixed Model (\code{NULL}).
##' Can be used to add missing rows relative to missing repetition (\code{"add.NA"})
##' or obtain a dataset with unique sets of covariates (\code{"unique"}) with respect to the mean structure.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @keywords methods
##' 
##' @examples
##' data("armd.wide", package = "nlmeU")
##' e.lmH <- lmm(visual52 ~ lesion, structure = IND(~treat.f), data = armd.wide)
##' model.frame(e.lmH)
##' model.frame(e.lmH, type = "unique")
##'
##' data("gastricbypassL", package = "LMMstar")
##' dfL.NNA <- na.omit(gastricbypassL)
##' e.lmm <- lmm(glucagonAUC ~ visit, repetition = ~visit|id, data = dfL.NNA, df = FALSE)
##' model.frame(e.lmm, type = "add.NA")

## * model.frame.lmm (code)
##' @export
model.frame.lmm <- function(formula, data = NULL, type = NULL, ...){

    ## ** extract from object
    manifest.var <- lava::manifest(formula)
    outcome.var <- formula$outcome$var
    cluster.var <- formula$cluster$var
    time.var <- formula$time$var
    n.time <- formula$time$n
    strata.var <- formula$strata$var

    ## ** normalize user imput
    dots <- list(...)    
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(!is.null(type)){
        type[1] <- match.arg(type[1], c("add.NA","unique"))
    }
            
    ## ** generate data.frame
    if(is.null(type) || type[1]=="add.NA"){

        ## *** reformat data
        if(is.null(data)){
            newdata <- formula$data
            design <- formula$design
        }else{
            if(any(setdiff(manifest.var,outcome.var) %in% names(data) == FALSE)){
                missing.col <- setdiff(manifest.var,outcome.var)[setdiff(manifest.var,outcome.var) %in% names(data)]
                stop("Incorrect argument \'data\' due to missing column(s): \"",paste(missing.col, collapse = "\", \""),"\". \n")
            }
            design <- model.matrix(formula, data = data, effects = "index")
            newdata <- .lmmNormalizeData(data, var.outcome = outcome.var,
                                         var.time = attr(formula$time$var,"original"),
                                         var.cluster = attr(formula$cluster$var,"original"),                         
                                         var.strata = attr(formula$strata$var,"original"),                         
                                         droplevels = TRUE,
                                         initialize.cluster = formula$design$vcov$ranef$crossed,
                                         initialize.time = setdiff(formula$design$vcov$ranef$vars, formula$cluster$var))$data
        }

        n.clusterTime <- lengths(design$index.clusterTime)
        incomplete.cluster <- which(n.clusterTime!=n.time)
        if(!is.null(type) && length(incomplete.cluster)>0){
            newdata <- newdata[intersect(names(newdata),c(manifest.var,"XXcluster.indexXX","XXtime.indexXX"))]

            ## *** classify covariates
            Xbaseline.var <- .baselineVar.lmm(formula)
            Xtime.var <- union(setdiff(attr(manifest.var,"mean"),c(Xbaseline.var)),attr(time.var,"original"))
            if(length(Xtime.var)>0){ ## find covariate that are function of time
                n.time <- formula$time$n
                index.Utime <- which(!duplicated(formula$data[order(formula$data$XXtime.indexXX),"XXtime.indexXX"]))
                XtimeFCT.var <- stats::setNames(lapply(Xtime.var, function(iVar){
                    if(length(unique(formula$data[[iVar]]))>n.time){
                        return(numeric(0))
                    }
                    iTable <- table(formula$data$XXtime.indexXX,formula$data[[iVar]])
                    if(all(rowSums(iTable!=0)<=1)){
                        return(formula$data[index.Utime,iVar])
                    }else{
                        return(numeric(0))
                    }
                }), Xtime.var)
                XtimeFCT.var <- XtimeFCT.var[length(XtimeFCT.var)>0]
                Xtime.var <- setdiff(Xtime.var, names(XtimeFCT.var))
            }else{
                XtimeFCT.var <- NULL
            }                
            
            ## *** augment dataset
            Xbaseline.var2 <- unique(stats::na.omit(c(Xbaseline.var,attr(cluster.var,"original"),attr(strata.var,"original"),"XXcluster.indexXX")))

            ls.newrow <- lapply(incomplete.cluster, function(iC){ ## iC <- 5
                iDF <- newdata[newdata$XXcluster.indexXX == iC,,drop=FALSE]
                iN.missingTime <- n.time - NROW(iDF)
                    
                iOut <- stats::setNames(vector(mode = "list", length = NCOL(iDF)), colnames(iDF))
                iOut[[outcome.var]] <- rep(as.numeric(NA), iN.missingTime)
                iOut$XXtime.indexXX <- setdiff(1:n.time,iDF$XXtime.indexXX)
                
                iOut[Xbaseline.var2] <-  lapply(Xbaseline.var2, function(iVar){
                    rep(iDF[1,iVar], iN.missingTime)
                })
                iOut[Xtime.var] <-  lapply(Xtime.var, function(iVar){
                    rep(NA, iN.missingTime)
                })
                iOut[names(XtimeFCT.var)] <-  lapply(XtimeFCT.var, function(iVec){
                    iVec[iOut$XXtime.indexXX]
                })
                return(as.data.frame(iOut))
            })
            newdata <- rbind(newdata[manifest.var],
                             do.call(rbind,ls.newrow)[manifest.var])
        }
    }else if(type[1]=="unique"){

        if(length(type)==1){
            effects <- "mean"
        }else{
            effects <- type[-1]
            if(length(effects) == 1 && effects == "all"){
                effects <- c("mean","variance","correlation")            
            }else{
                effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)
            }
        }   

        if(is.null(data)){
            design <- formula$design
        }else{
            design <- model.matrix(formula, data = data, effects = "index")
        }
        cluster.var <- formula$cluster$var
        time.var <- formula$time$var

        
        ## *** reorder data by cluster and time
        data.index.cluster <- attr(design$index.cluster,"vectorwise")
        data.index.time <- attr(design$index.clusterTime,"vectorwise")                    
        vec.reorder <- order(data.index.cluster, data.index.time)
        if(is.null(data)){
            data <- formula$data[vec.reorder,,drop=FALSE]
        }

        newdata <- NULL
        if("variance" %in% effects || "correlation" %in% effects){
            keep.cluster <- sapply(attr(design$vcov$pattern,"list"),"[[",1)
            newdata.row <- which(data$XXcluster.indexXX %in% keep.cluster)
            vcov.var <- lava::manifest(formula, effects = setdiff(effects, "mean"))
            newdata.col <- stats::na.omit(unique(c(attr(cluster.var,"original"),unlist(vcov.var))))
        }else{
            newdata.row <- NULL
            newdata.col <- NULL
        }

        if("mean" %in% effects){
            mean.var <- lava::manifest(formula, effects = "mean")

            ## but putting first clusters (possibly) kept for the covariance covariates
            vec.reorder2 <- c(newdata.row,setdiff(vec.reorder,newdata.row))
            ## extract unique mean covariate sets
            newdata.row <- sort(union(newdata.row, which(!duplicated(data[,names(data) %in% mean.var,drop=FALSE]))))
            newdata.col <- union(newdata.col, mean.var)
        }

        ## *** subset
        newdata <- data[newdata.row,names(data) %in% newdata.col,drop=FALSE]
        attr(newdata,"index.cluster") <- data[newdata.row,"XXcluster.indexXX"]
        attr(newdata,"index.time") <- data[newdata.row,"XXtime.indexXX"]
        attr(newdata,"index.strata") <- data[newdata.row,"XXstrata.indexXX"]
    }
        

    ## ** export
    rownames(newdata) <- NULL
    return(newdata)

}

## * model.frame.lmmCC
##' @export
model.frame.lmmCC <- model.frame.lmm

##----------------------------------------------------------------------
### model.frame.R ends here
