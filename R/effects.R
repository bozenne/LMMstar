### effects.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 29 2024 (09:47) 
## Version: 
## Last-Updated: jan 29 2024 (13:58) 
##           By: Brice Ozenne
##     Update #: 69
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## add conditional effect for interactions

effects.lmm <- function(object, type, variable, rhs = NULL, repetition = NULL, multivariate = FALSE, ...){

    call <- match.call()

    ## ** check arguments
    if(!is.character(type)){
        stop("Argument \'type\' must be a character value. \n")
    }
    if(length(type)!=1){
        stop("Argument \'type\' must have length 1. \n")
    }
    type <- tolower(type)
    valid.type <- c("mean","variance","correlation","ate","aco")
    if(type %in% valid.type == FALSE){
        stop("Argument \'type\' must be on of \"",paste(valid.type,collapse ="\", \""),"\". \n")
    }
    if(!is.character(variable)){
        stop("Argument \'variable\' must be a character value. \n")
    }
    if(length(variable)!=1){
        stop("Argument \'variable\' must have length 1. \n")
    }
    
    if(type %in% c("mean","ate","aco")){
        valid.variable <- attr(object$design$mean, "variable")
        rhs <- 0
    }else if(type=="variance"){
        valid.variable <- attr(object$design$vcov$var$X, "variable")
        transform.k <- list(...)$transform.k
        if(!is.null(transform.k) && transform.k %in% c("none","square")){
            rhs <- 1
        }else{
            rhs <- 0
        }
    }else if(type=="correlation"){
        valid.variable <- attr(object$design$vcov$cor$X, "variable")
        rhs <- 0
    }
    if(length(valid.variable) == 0){
        stop("Cannot evaluate the effect as no covariate was specified when fitting the model in the ",ifelse(type %in% c("ate","aco"),"mean",type)," structure. \n")
    }
    if(variable %in% valid.variable == FALSE){
        stop("Argument \'variable\' must be one of \"",paste(valid.variable,collapse ="\", \""),"\". \n")
    }

    ## ** normalized user input    
    if(type %in% c("mean","variance","correlation")){
        X <- switch(type,
                    "mean" = object$design$mean,
                    "variance" = object$design$vcov$var$X,
                    "correlation" = object$design$vcov$cor$X)

        assign <- which(attr(,"variable")==variable)
        labels <- colnames(object$design$vcov$var$X)[attr(object$design$vcov$var$X,"assign")==assign]
        effect <- paste0(labels,"=",rhs)
    }else if(type %in% c("aco","ate")){
        if(variable %in% names(object$xfactor$mean)){
            Ulevel.variable <- object$xfactor$mean[[variable]]
        }else{
            Ulevel.variable <- sort(unique(object$data[[variable]]))
        }
        if(length(Ulevel.variable)>2 && is.numeric(object$data[[variable]])){
            if(type == "aco"){
                stop("Cannot compute average treatment effects for numeric variables (unless they only have two levels). \n")
            }else if(type == "ate"){
                stop("Cannot compute average counterfactuals outcomes for numeric variables (unless they only have two levels). \n")
            }
        }

        if(type == "ate"){
            Upair.variable <- unorderedPairs(Ulevel.variable, distinct = TRUE)
        }
        ls.indexTime <- tapply(object$data$XXindexXX,object$data$XXtime.indexXX,identity)
        index.original <- model.matrix(object, data = object$data.original, effects = "index")

        if(type == "aco"){
            ls.C <- lapply(Ulevel.variable, function(iLevel){ ## iLevel <- 0
                iData <- object$data.original
                iData[[variable]][] <- iLevel
                iX <- model.matrix(object$formula$mean.design, iData)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint
                if(NROW(iX)==NROW(iData)){
                    iC <- do.call(rbind,by(iX,attr(index.original$index.clusterTime,"vectorwise"),colMeans))
                }else{ ## handle missing values in covariates, i.e. rows removed by model.matrix
                    iC <- do.call(rbind,by(iX,attr(index.original$index.clusterTime,"vectorwise")[as.numeric(rownames(iX))],colMeans))
                }
                rownames(iC) <- paste0(iLevel,"(",object$time$levels,")")
                return(iC)
            })
        }else if(type == "ate"){
            ls.C <- lapply(1:NCOL(Upair.variable), function(iPair){ ## iPair <- 1
                iData1 <- object$data.original
                iData2 <- object$data.original
                iData1[[variable]][] <- Upair.variable[1,iPair]
                iData2[[variable]][] <- Upair.variable[2,iPair]
                iX1 <- model.matrix(object$formula$mean.design, iData1)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint
                iX2 <- model.matrix(object$formula$mean.design, iData2)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint
                if(NROW(iX1)==NROW(iData1)){ ## handle missing values in covariates, i.e. rows removed by model.matrix
                    iC <- do.call(rbind,by(iX2-iX1,attr(index.original$index.clusterTime,"vectorwise"),colMeans))
                }else{
                    iC <- do.call(rbind,by(iX2-iX1,attr(index.original$index.clusterTime,"vectorwise")[as.numeric(rownames(iX))],colMeans))
                }
                rownames(iC) <- paste0(Upair.variable[2,iPair],"-",Upair.variable[1,iPair],"(",object$time$levels,")")
                return(iC)
            })
        }
        effect <- do.call(rbind,ls.C)
    }

    ## ** check arguments
    out <- anova(object, effect = effect, multivariate = multivariate, ...)
    out$args$effect <- type
    out$args$variable <- variable

    ## ** export
    attr(out,"class") <- append("effect_lmm",attr(out,"class"))
    return(out)
}
##----------------------------------------------------------------------
### effects.R ends here
