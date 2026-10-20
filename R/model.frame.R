### model.frame.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun  7 2021 (14:57) 
## Version: 
## Last-Updated: okt 19 2026 (11:33) 
##           By: Brice Ozenne
##     Update #: 278
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
##' @description Contruct a data frame containing all the variables involved in a Linear Mixed Model.
##' @param formula a \code{lmm} object.
##' @param newdata [data.frame] dataset relative to which the model frame should be constructed.
##' @param type [character] dataset to be output: \itemize{
##' \item \code{"processed"}: dataset used to fit the Linear Mixed Model
##' \item \code{"original"}: original dataset provided by the user
##' \item \code{"unique"}: unique sets of covariates with respect to the mean structure.
##' \item \code{"extended"}: \code{"processed"} dataset where rows relative to missing repetitions have been added.
##' }.
##' @param add.index [logical] should columns indexing the row number from the original dataset, time variable, cluster variable, strata variable
##' be added to the output?
##' @param na.rm [logical] should rows containing missing values for the variables used in the linear mixed model be removed?
##' Not relevant when argument type is \code{"unique"}.
##' @param ... not used. For compatibility with the generic method.
##'
##' @details Column \code{"XXindexXX"} refers to the row of the original dataset (i.e. passed to argument \code{data} when calling \code{\link{lmm.formula}}).
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
model.frame.lmm <- function(formula, newdata = NULL, type = "processed", add.index = type %in% c("processed","extended"), na.rm = TRUE, ...){


    mycall <- match.call()

    ## ** extract from object
    table.variables <- stats::model.tables(formula, effects = "variables")
    var.outcome <- table.variables["outcome","var"]
    var.cluster <- table.variables["cluster","var"]
    var.time <- table.variables["time","var"]
    var.strata <- table.variables["strata","var"]
    n.time <- table.variables["time","n"]
    browser()
    index.cluster <- index(formula, variable = "cluster", type = "factor")

    ## ** normalize user imput
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

    ## *** newdata
    if(!is.null(newdata)){
        if(any(setdiff(var.manifest,var.outcome) %in% names(newdata) == FALSE)){
            missing.col <- setdiff(var.manifest,var.outcome)[setdiff(var.manifest,var.outcome) %in% names(newdata)]
            stop("Incorrect argument \'newdata\' due to missing column(s): \"",paste(missing.col, collapse = "\", \""),"\". \n")
        }
    }

    ## *** type
    if(length(type)!=1){
        stop("Argument \'type\' should have length 1. \n")
    }
    if(type %in% c("processed","original","unique","extended") == FALSE){
        stop("Argument \'type\' should be one of \"processed\", \"original\", \"unique\", \"extended\". \n",
             "Proposed value: \"",type,"\".")
    }

    ## *** na.rm
    if(type=="extended"){
        if("na.rm" %in% names(mycall) && mycall$na.rm){
            stop("Argument \'na.rm\' must be FALSE when arguement \'type\' equals \"add.NA\". \n")
        }
        na.rm <- FALSE        
    }

    ## ** retrieve original data.frame
    if(is.null(newdata)){
        test.newdata <- FALSE
        newdata <- formula$data
        index.na <- NULL
    }else{
        test.newdata <- TRUE
        index.na <- stats::na.omit(formula)
    }
    browser()
    ## ** add index
    if(add.index){
        if(test.newdata){
        }else{
            index.cluster <- formula$design$index.cluster

            index.nna <- setdiff(1:NROW(newdata), index.na)
            out$XXindexXX <- 1:NROW(newdata)
            browser()
            
            if(table.variables["strata","original"]){
                out$XXclusterXX <- newdata[[var.cluster]]
            }else{
                out$XXclusterXX <- 1
            }
            if(table.variables["time","original"]){
                out$XXtimeXX <- newdata[[var.time]]
            }else{
                out$XXtimeXX <- 1
            }            
            if(table.variables["strata","original"]){
                out$XXstrataXX <- newdata[[var.strata]]
            }else{
                out$XXstrataXX <- 1
            }
            out$XXcluster.indexXX <- 1
            out$XXtime.indexXX <- 1
            out$XXstrata.indexXX <-1             
        }
        
    }

    if(type == "original"){
        

        
        
    }else if(type %in% c("processed","extended")){

        ## *** reformat newdata
        if(is.null(newdata) && (is.null(formula$index.na) || na.rm)){
            out <- formula$data
            index.na <- formula$index.na
            index.data <- formula$design[c("index.cluster","index.clusterTime","index.clusterStrata")]
        }else{ ## new data or obtain model frame without removing the row with missing values
            
            if(is.null(newdata)){
                newdata <- formula$data
            }else if(any(setdiff(var.manifest,var.outcome) %in% names(newdata) == FALSE)){
                missing.col <- setdiff(var.manifest,var.outcome)[setdiff(var.manifest,var.outcome) %in% names(newdata)]
                stop("Incorrect argument \'newdata\' due to missing column(s): \"",paste(missing.col, collapse = "\", \""),"\". \n")
            }
            index.data <- stats::model.matrix(formula, newdata = newdata, effects = "index", na.rm = na.rm, simplify = FALSE)
            out <- index.data$data
            if(na.rm){
                index.na <- index.data$index.na
            }else{
                index.na <- NULL
            }
        }
        
    }else if(type=="unique"){

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
            design <- stats::model.matrix(formula, newdata = newdata, effects = "index")
        }
        var.cluster <- formula$cluster$var
        var.time <- formula$time$var

        
        ## *** reorder data by cluster and time
        newdata.index.cluster <- attr(design$index.cluster,"vectorwise")
        newdata.index.time <- attr(design$index.clusterTime,"vectorwise")                    
        vec.reorder <- order(newdata.index.cluster, newdata.index.time)
        if(is.null(newdata)){
            newdata <- formula$data[vec.reorder,,drop=FALSE]
        }

        ## *** keep unique mean and/or covariance pattern
        if("variance" %in% effects || "correlation" %in% effects){
            keep.cluster <- sapply(attr(design$vcov$pattern,"list"),"[[",1)
            newdata.row <- which(newdata$XXcluster.indexXX %in% keep.cluster)
            var.vcov <- stats::variable.names(formula, effects = setdiff(effects, "mean"))
            newdata.col <- stats::na.omit(unique(c(attr(var.cluster,"original"),unlist(var.vcov))))
        }else{
            newdata.row <- NULL
            newdata.col <- NULL
        }

        if("mean" %in% effects){
            var.mean <- stats::variable.names(formula, effects = "mean")

            ## but putting first clusters (possibly) kept for the covariance covariates
            vec.reorder2 <- c(newdata.row,setdiff(vec.reorder,newdata.row))
            ## extract unique mean covariate sets
            newdata.row <- sort(union(newdata.row, which(!duplicated(newdata[,names(newdata) %in% var.mean,drop=FALSE]))))
            newdata.col <- union(newdata.col, var.mean)
        }
        ## *** subset
        out <- newdata[newdata.row,,drop=FALSE]
        index.na <- NULL
        
    }

    ## ** add missing row (NA outcome)
    if(type == "extended"){

        n.clusterTime <- lengths(index.data$index.clusterTime)
        incomplete.cluster <- which(n.clusterTime!=n.time)

        if(length(incomplete.cluster)>0){

            if(all(is.na(attr(var.time,"original")))){
                stop("Cannot add rows relative to missing repetitions without time variable. \n",
                     "Consider specifying the argument \'repetiton\' when calling lmm() with something like ~time|id. \n")
            }

            ## *** classify covariates
            mean.type <- stats::variable.names(formula, effects = "mean.type")
            var.meanBaseline <- names(mean.type)[mean.type=="baseline"] ## baseline covariate
            var.meanTime <- names(mean.type)[mean.type=="time"] ## time covariate or simple transformation of time
            var.meanTimeVar <- names(mean.type)[mean.type=="timevar"] ## time varying covariate
            ls.tableTime.var <- attr(mean.type,"table") ## conversion from time variable to covariate 
            
            ## *** augment dataset
            var2.meanBaseline <- unique(stats::na.omit(c(var.meanBaseline, attr(var.cluster, "original"), attr(var.strata, "original"), "XXclusterXX", "XXstrataXX", "XXcluster.indexXX", "XXstrata.indexXX")))
            var2.meanTime <- unique(stats::na.omit(c(var.meanTime, attr(var.time, "original"))))
            if(any(var2.meanTime %in% names(ls.tableTime.var) == FALSE)){
                var.newtime <- var2.meanTime[var2.meanTime %in% names(ls.tableTime.var) == FALSE]
                ls.tableTime.var[var.newtime] <- lapply(var.newtime, function(iVar){ ## iVar <- var.newtime
                    attr(formula$time$level,"original")[,iVar]
                })
            }

            ls.newrow <- lapply(incomplete.cluster, function(iC){ ## iC <- 1
                iDF <- out[index.data$index.cluster[[iC]],,drop=FALSE]
                iN.missingTime <- n.time - NROW(iDF)

                iOut <- stats::setNames(vector(mode = "list", length = NCOL(iDF)), colnames(iDF))
                iOut[[var.outcome]] <- rep(as.numeric(NA), iN.missingTime)
                iOut$XXindexXX <- rep(as.numeric(NA), iN.missingTime)
                iOut$XXtime.indexXX <- setdiff(1:n.time,iDF$XXtime.indexXX)
                iOut$XXtimeXX <- formula$time$levels[iOut$XXtime.indexXX]
                
                iOut[var2.meanBaseline] <-  lapply(var2.meanBaseline, function(iVar){
                    rep(iDF[1,iVar], iN.missingTime)
                })
                iOut[var2.meanTime] <-  lapply(var2.meanTime, function(iVar){ ## iVar <- meanTime.var
                    ls.tableTime.var[[iVar]][iOut$XXtime.indexXX]
                })
                iOut[var.meanTimeVar] <-  lapply(var.meanTimeVar, function(iVar){
                    rep(NA, iN.missingTime)
                })
                return(as.data.frame(iOut))
            })
            out <- rbind(out, do.call(rbind,ls.newrow))
            out[is.na(out$XXindexXX),"XXindexXX"] <- -(1:sum(is.na(out$XXindexXX)))
        }
    }

    ## ** export
    if(length(index.na)>0){
        attr(out,"index.na") <- index.na
    }
    return(out)

}

## * model.frame.lmmCC
##' @export
model.frame.lmmCC <- model.frame.lmm

##----------------------------------------------------------------------
### model.frame.R ends here
