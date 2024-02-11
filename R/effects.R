### effects.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 29 2024 (09:47) 
## Version: 
## Last-Updated: Feb 11 2024 (23:26) 
##           By: Brice Ozenne
##     Update #: 165
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * effects.lmm (documentation)
##' @title Effects Derived For Linear Mixed Model
##' @description Estimate mean effects from linear mixed models relative to a variable.

##' @param object a \code{lmm} object. 
##' @param type [character] should all mean coefficients with respect to a variable be tested (\code{"mean"})
##' or the average counterfactual outcome (\code{"aoc"}) or average treatment effect (\code{"ate"}) be estimated?
##' @param variable [character] the variable relative to which the effect should be computed.
##' @param conditional [character] variable conditional to which the average conterfactual outcome or treatment effect should be computed.
##' @param rhs [numeric] the right hand side of the hypothesis.
##' @param repetition [character vector] repetition at which the effect should be assessed. By default it will be assessed at all repetitions.
##' @param multivariate [logical] should a multivariate Wald test be used to simultaneously test all null hypotheses.
##' @param prefix.time [character] When naming the estimates, text to be pasted before the value of the repetition variable.
##' Only relevant when \code{type = "aoc"} or \code{type = "ate"}.
##' @param prefix.var [logical] When naming the estimates, should the variable name be added or only the value?
##' @param sep.var [character] When naming the estimates, text to be pasted between the values to condition on.
##' Only relevant when \code{type = "aoc"} or \code{type = "ate"}.
##' @param ... Arguments passed to \code{\link{anova.lmm}}.
##'
##' @examples
##' #### simulate data in the long format ####
##' set.seed(10)
##' dL <- sampleRem(100, n.times = 3, format = "long")
##' 
##' #### fit Linear Mixed Model ####
##' eUN.lmm <- lmm(Y ~ visit + X1 + X2 + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' effects(eUN.lmm, type = "mean", variable = "X1")
##' effects(eUN.lmm, type = "aco", variable = "X1")
##' effects(eUN.lmm, type = "ate", variable = "X1")
##' effects(eUN.lmm, type = "ate", variable = "X1", repetition = "3")
##' 
##' #### fit Linear Mixed Model with interaction ####
##' dL$X1.factor <- as.factor(dL$X1)
##' dL$X2.factor <- as.factor(dL$X2)
##' eUN.lmmI <- lmm(Y ~ visit + X1.factor * X2.factor + X5,
##'                repetition = ~visit|id, structure = "UN", data = dL)
##'
##' effects(eUN.lmmI, type = "aco", variable = "X1.factor", conditional = c("X2.factor"), repetition = "3")
##' effects(eUN.lmmI, type = "ate", variable = "X1.factor", conditional = c("X2.factor"), repetition = "3")
##' effects(eUN.lmmI, type = "ate", variable = "X1.factor", conditional = list(X5=0:1), repetition = "3")
##' effects(eUN.lmmI, type = "ate", variable = "X1.factor", conditional = list(X5=0:5), repetition = "3")


## * effects.lmm (code)
##' @export
effects.lmm <- function(object, type, variable, conditional = NULL, rhs = NULL, repetition = NULL, multivariate = FALSE,
                        prefix.time = "t=", prefix.var = TRUE, sep.var = ",", ...){

    call <- match.call()

    ## ** check arguments
    if(!is.character(type)){
        stop("Argument \'type\' must be a character value. \n")
    }
    if(length(type)!=1){
        stop("Argument \'type\' must have length 1. \n")
    }
    type <- tolower(type)
    valid.type <- c("mean","ate","aco")
    if(type %in% valid.type == FALSE){
        stop("Argument \'type\' must be on of \"",paste(valid.type,collapse ="\", \""),"\". \n")
    }
    if(!is.character(variable)){
        stop("Argument \'variable\' must be a character value. \n")
    }
    valid.variable <- attr(object$design$mean, "variable")
    rhs <- 0

    if(length(valid.variable) == 0){
        stop("Cannot evaluate the effect as no covariate was specified when fitting the model in the ",ifelse(type %in% c("ate","aco"),"mean",type)," structure. \n")
    }
    if(length(variable)!=1){
        stop("Argument \'variable\' must have length 1. \n")
    }
    if(variable %in% valid.variable == FALSE){
        stop("Argument \'variable\' must be one of \"",paste(valid.variable,collapse ="\", \""),"\". \n")
    }
    ## ** normalized user input    
    if(type %in% c("mean")){
        if(!is.null(repetition)){
            warning("Argument \'repetition\' is ignored when argument \'type\' is \"mean\". \n")
        }

        X <- object$design$mean
        assign <- which(attr(X,"variable")==variable)
        labels <- colnames(X)[attr(X,"assign")==assign]
        effect <- paste0(labels,"=",rhs)
    }else if(type %in% c("aco","ate")){

        ## find unique levels of the exposure variable
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

        ## recover the time variable
        index.original <- model.matrix(object, data = object$data.original, effects = "index")
        time.original <- object$time$levels[attr(index.original$index.clusterTime,"vectorwise")]

        ## only consider requested times
        if(!is.null(repetition)){
            if(!all(repetition %in%  object$time$levels) &&  !all(repetition %in%  1:length(object$time$levels))){
                stop("Argument \'repetition\' incorrect. Should either be: \n",
                     "- any of \"",paste0(object$time$levels, collapse = "\", \""),"\" \n",
                     "- any integer from 1 to ",length(object$time$levels)," \n")
            }
            if(is.numeric(repetition)){
                repetition <- object$time$levels[repetition]
            }
            subset.repetition <- which(time.original %in% repetition)
            data.original <- cbind(object$data.original, XXtimeXX = time.original)[subset.repetition,,drop=FALSE]
        }else{
            subset.repetition <- NULL
            data.original <- cbind(object$data.original, XXtimeXX = time.original)
        }
        ## define conditioning set (normalize user input)
        if(is.character(conditional)){
            if(any(conditional %in% setdiff(names(object$xfactor$mean),variable) == FALSE)){
                stop("Incorrect value for argument \'conditional\'. \n",
                     "If a character should only include factor variables influencing the mean structure.\n",
                     "Valid values: \"",paste0(setdiff(names(object$xfactor$mean),variable), collapse = "\" \""),"\".\n")
            }
            grid.conditional <- expand.grid(object$xfactor$mean[conditional])
        }else if(is.list(conditional)){
            if(is.null(names(conditional))){
                stop("Missing names for argument \'conditional\'. \n",
                     "Valid values: \"",paste0(setdiff(valid.variable, variable), collapse = "\" \""),"\".\n")
            }
            if(any(names(conditional) %in% setdiff(valid.variable, variable) == FALSE)){
                stop("Incorrect names for argument \'conditional\'. \n",
                     "Valid values: \"",paste0(setdiff(valid.variable, variable), collapse = "\" \""),"\".\n")
            }
            grid.conditional <- expand.grid(conditional)
        }else if(!is.null(conditional)){
            stop("Argument \'conditional\' must either be NULL, or a character vector indicating factor variables, or a named list of covariate values. \n")
        }
        if(!is.null(conditional)){            
            key.conditional <- interaction(grid.conditional, sep = sep.var)
            key.original <- levels(key.conditional)[match(interaction(data.original[names(grid.conditional)], sep =  sep.var), key.conditional)]
            if(prefix.var){
                if(length(grid.conditional)==1){
                    data.original$XXstrataXX <- paste0(names(grid.conditional),"=",key.original)
                }else{
                    newlevel <- interaction(as.data.frame(lapply(names(grid.conditional), function(iName){paste0(iName,"=",grid.conditional[,iName])})), sep = sep.var)
                    data.original$XXstrataXX <- factor(key.original, levels = levels(key.conditional), labels = newlevel)
                }
            }else{
                data.original$XXstrataXX <- key.original
            }
            if(any(is.na(key.original))){
                data.original <- data.original[!is.na(key.original),,drop=FALSE]
            }
        }


        ## perform average
        if(type == "aco"){
            ls.C <- lapply(Ulevel.variable, function(iLevel){ ## iLevel <- 0
                iData <- data.original
                iData[[variable]][] <- iLevel
                iX <- model.matrix(object$formula$mean.design, iData)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint

                if(!is.null(conditional)){
                    iStrata <- paste0(iData$XXtimeXX, sep.var, iData$XXstrataXX)
                }else{
                    iStrata <- iData$XXtimeXX
                }

                if(NROW(iX)!=NROW(iData)){
                    ## handle missing values in covariates, i.e. rows removed by model.matrix
                    iStrata <- iStrata[as.numeric(rownames(iX))]
                }
                
                iC <- do.call(rbind,by(iX,iStrata,colMeans))
                rownames(iC) <- paste0(iLevel,"(",prefix.time,rownames(iC),")")
                if(prefix.var){
                    rownames(iC) <- paste0(variable,"=",rownames(iC))
                }
                return(iC)
            })
        }else if(type == "ate"){
            Upair.variable <- unorderedPairs(Ulevel.variable, distinct = TRUE)

            ls.C <- lapply(1:NCOL(Upair.variable), function(iPair){ ## iPair <- 1
                iData1 <- data.original
                iData2 <- data.original

                iData1[[variable]][] <- Upair.variable[1,iPair]
                iData2[[variable]][] <- Upair.variable[2,iPair]

                iX1 <- model.matrix(object$formula$mean.design, iData1)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint
                iX2 <- model.matrix(object$formula$mean.design, iData2)[,colnames(object$design$mean),drop=FALSE] ## remove uncessary columns in case of (baseline) constraint

                if(!is.null(conditional)){
                    iStrata <- paste0(Upair.variable[2,iPair],"-",Upair.variable[1,iPair],"(",prefix.time,iData1$XXtimeXX, sep.var, iData1$XXstrataXX,")")
                }else{
                    iStrata <- paste0(Upair.variable[2,iPair],"-",Upair.variable[1,iPair],"(",prefix.time,iData1$XXtimeXX,")")
                }

                if(NROW(iX1)!=NROW(iData1)){
                    ## handle missing values in covariates, i.e. rows removed by model.matrix
                    iStrata <- iStrata[as.numeric(rownames(iX1))]
                }

                iC <- do.call(rbind,by(iX2-iX1,iStrata,colMeans))
                if(prefix.var){
                    rownames(iC) <- paste0(variable,"=",rownames(iC))
                }
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

