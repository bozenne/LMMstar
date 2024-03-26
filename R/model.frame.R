### model.frame.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (14:57) 
## Version: 
## Last-Updated: Mar 26 2024 (09:46) 
##           By: Brice Ozenne
##     Update #: 124
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
##' @param newdata [data.frame] dataset relative to which the model frame should be constructed.
##' @param type [character] By default returns the processed dataset used to fit the Linear Mixed Model (\code{NULL}).
##' Can be used to add rows relative to missing repetitions (\code{"add.NA"})
##' or obtain a dataset with unique sets of covariates (\code{"unique"}) with respect to the mean structure.
##' @param add.index [logical] Should columns indexing the row number from the original dataset, time variable, cluster variable, strata variable
##' be added to the output?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details Column \code{"XXindexXX"} refers to the row of the original dataset (i.e. passed to argument \code{data} when calling \code{\link{lmm}}).
##' When adding rows relative to missing repetitions, since there is no row in the original dataset, a negative sign is used.
##' 
##' @keywords methods
##' 
##' @examples
##' data("armd.wide", package = "nlmeU")
##' e.lmH <- lmm(visual52 ~ lesion, structure = IND(~treat.f), data = armd.wide)
##' model.frame(e.lmH)
##' model.frame(e.lmH, add.index = TRUE)
##' model.frame(e.lmH, type = "unique")
##'
##' data("gastricbypassL", package = "LMMstar")
##' dfL.NNA <- na.omit(gastricbypassL)
##' e.lmm <- lmm(glucagonAUC ~ time, repetition = ~visit|id, data = dfL.NNA, df = FALSE)
##' model.frame(e.lmm, type = "unique")
##' model.frame(e.lmm, type = c("unique","correlation"))
##' model.frame(e.lmm, type = "add.NA", add.index = TRUE)

## * model.frame.lmm (code)
##' @export
model.frame.lmm <- function(formula, newdata = NULL, type = NULL, add.index = FALSE, ...){

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

        ## *** reformat newdata
        if(is.null(newdata)){
            newdata.norm <- formula$data
            design <- formula$design
        }else{
            if(any(setdiff(manifest.var,outcome.var) %in% names(newdata) == FALSE)){
                missing.col <- setdiff(manifest.var,outcome.var)[setdiff(manifest.var,outcome.var) %in% names(newdata)]
                stop("Incorrect argument \'newdata\' due to missing column(s): \"",paste(missing.col, collapse = "\", \""),"\". \n")
            }
            design <- model.matrix(formula, newdata = newdata, effects = "index")
            newdata.norm <- .lmmNormalizeData(newdata, var.outcome = outcome.var,
                                         var.time = attr(formula$time$var,"original"),
                                         var.cluster = attr(formula$cluster$var,"original"),                         
                                         var.strata = attr(formula$strata$var,"original"),                         
                                         droplevels = TRUE,
                                         initialize.cluster = formula$design$vcov$ranef$crossed,
                                         initialize.time = setdiff(formula$design$vcov$ranef$vars, formula$cluster$var))$data
        }

        ## *** add missing row (NA outcome)
        n.clusterTime <- lengths(design$index.clusterTime)
        incomplete.cluster <- which(n.clusterTime!=n.time)
        if(!is.null(type) && length(incomplete.cluster)>0){
            
            ## classify covariates
            Xbaseline.var <- .baselineVar.lmm(formula)
            Xtime.var <- union(setdiff(attr(manifest.var,"mean"),Xbaseline.var),attr(time.var,"original"))
            if(length(Xtime.var)>0){ ## find covariate that are function of time
                n.time <- formula$time$n
                object.data <- formula$data[order(formula$data$XXtime.indexXX),,drop=FALSE]
                index.Utime <- which(!duplicated(object.data$XXtime.indexXX))
                XtimeFCT.var <- stats::setNames(lapply(Xtime.var, function(iVar){ ## iVar <- "Time"
                    if(length(unique(formula$data[[iVar]]))>n.time){
                        return(numeric(0))
                    }
                    iTable <- table(formula$data$XXtime.indexXX,formula$data[[iVar]])
                    if(all(rowSums(iTable!=0)<=1)){
                        return(object.data[index.Utime,iVar])
                    }else{
                        return(numeric(0))
                    }
                }), Xtime.var)
                XtimeFCT.var <- XtimeFCT.var[length(XtimeFCT.var)>0]
                Xtime.var <- setdiff(Xtime.var, names(XtimeFCT.var))
            }else{
                XtimeFCT.var <- NULL
            }                
            ## augment dataset
            Xbaseline.var2 <- unique(stats::na.omit(c(Xbaseline.var,attr(cluster.var,"original"),attr(strata.var,"original"),
                                                      "XXclusterXX","XXstrataXX","XXcluster.indexXX","XXstrata.indexXX")))

            ls.newrow <- lapply(incomplete.cluster, function(iC){ ## iC <- 5
                iDF <- newdata.norm[newdata.norm$XXcluster.indexXX == iC,,drop=FALSE]
                iN.missingTime <- n.time - NROW(iDF)
                    
                iOut <- stats::setNames(vector(mode = "list", length = NCOL(iDF)), colnames(iDF))
                iOut[[outcome.var]] <- rep(as.numeric(NA), iN.missingTime)
                iOut$XXindexXX <- as.numeric(NA)
                iOut$XXtime.indexXX <- setdiff(1:n.time,iDF$XXtime.indexXX)
                iOut$XXtimeXX <- formula$time$levels[iOut$XXtime.indexXX]
                
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
            out <- rbind(newdata.norm, do.call(rbind,ls.newrow))
            out[is.na(out$XXindexXX),"XXindexXX"] <- -(1:sum(is.na(out$XXindexXX)))
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

        if(is.null(newdata)){
            design <- formula$design
        }else{
            design <- model.matrix(formula, data = newdata, effects = "index")
        }
        cluster.var <- formula$cluster$var
        time.var <- formula$time$var

        
        ## *** reorder data by cluster and time
        newdata.index.cluster <- attr(design$index.cluster,"vectorwise")
        newdata.index.time <- attr(design$index.clusterTime,"vectorwise")                    
        vec.reorder <- order(newdata.index.cluster, newdata.index.time)
        if(is.null(newdata)){
            newdata <- formula$newdata[vec.reorder,,drop=FALSE]
        }

        if("variance" %in% effects || "correlation" %in% effects){
            keep.cluster <- sapply(attr(design$vcov$pattern,"list"),"[[",1)
            newdata.row <- which(newdata$XXcluster.indexXX %in% keep.cluster)
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
            newdata.row <- sort(union(newdata.row, which(!duplicated(newdata[,names(newdata) %in% mean.var,drop=FALSE]))))
            newdata.col <- union(newdata.col, mean.var)
        }
        ## *** subset
        out <- newdata[newdata.row,,drop=FALSE]
   
    }

    ## ** export
    if(add.index==FALSE){
        out <- out[manifest.var]
    }
    return(out)

}

## * model.frame.lmmCC
##' @export
model.frame.lmmCC <- model.frame.lmm

##----------------------------------------------------------------------
### model.frame.R ends here
