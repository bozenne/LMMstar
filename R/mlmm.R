### mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 14 2022 (09:45) 
## Version: 
## Last-Updated: aug 25 2022 (15:35) 
##           By: Brice Ozenne
##     Update #: 154
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mlmm (documentation)
##' @title Fit Multiple Linear Mixed Model
##' @description Fit several linear mixed models, extract relevant coefficients, and combine them into a single table. 
##'
##' @param ... arguments passed to \code{\link{lmm}}.
##' @param data [data.frame] dataset (in the long format) containing the observations.
##' @param by [character] variable used to split the dataset. On each split a seperate linear mixed model is fit.
##' @param effects [character or numeric matrix] Linear combinations of coefficients relative to which Wald test should be computed.
##' Argument passed to \code{\link{anova.lmm}}.
##' Right hand side can be specified via an attribute \code{"rhs"}.
##' @param contrast.rbind [character or numeric matrix] Contrast to be be applied to compare the groups.
##' Argument passed to the argument \code{effects} of \code{\link{rbind.Wald_lmm}}.
##' Right hand side can be specified via an attribute \code{"rhs"}.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors.
##' Argument passed to \code{\link{anova.lmm}}.
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##' Argument passed to \code{\link{anova.lmm}}.
##' @param ci [logical] Should a confidence interval be output for each hypothesis?
##' Argument passed to \code{\link{anova.lmm}}.
##'
##' @details
##' \bold{Grouping variable} in argument repetition: when numeric, it will be converted into a factor variable, possibly adding a leading 0 to preserve the ordering.
##' This transformation may cause inconsistency when combining results between different \code{lmm} object. 
##' This is why the grouping variable should preferably be of type character or factor.
##' 
##' @examples
##'
##' #### univariate regression ####
##' if(require(lava)){
##' 
##' set.seed(10)
##' d1 <- cbind(sim(lvm(Y~0.5*X1), 25), group = "A")
##' d2 <- cbind(sim(lvm(Y~0.1*X1), 100), group = "B")
##' d3 <- cbind(sim(lvm(Y~0.01*X1), 1000), group = "C")
##' d1$id <- 1:NROW(d1)
##' d2$id <- 1:NROW(d2)
##' d3$id <- 1:NROW(d3)
##' 
##' d <- rbind(d1,d2,d3)
##' 
##' e.mlmm <- mlmm(Y~X1, data = d, by = "group", effects = "X1=0")
##' summary(e.mlmm, method = "single-step")
##' summary(e.mlmm, method = "bonferroni")
##' summary(e.mlmm, method = "single-step2")
##' ## summary(e.mlmm)
##' }
##' 
##' #### multivariate regression ####
##' set.seed(10)
##' dL <- sampleRem(250, n.times = 3, format = "long")
##'
##' e.mlmm <- mlmm(Y~X1+X2+X6, repetition = ~visit|id, data = dL,
##'                by = "X4", structure = "CS")
##' summary(e.mlmm, method = "none")
##' confint(e.mlmm, method = "none")
##' 
##' e.mlmmX1 <- mlmm(Y~X1+X2+X6, repetition = ~visit|id, data = dL,
##'                by = "X4", effects = "X1=0", structure = "CS")
##' summary(e.mlmmX1)
##' summary(e.mlmmX1, method = "single-step")
##' 

## * mlmm (code)
##' @export
mlmm <- function(..., data, by, contrast.rbind = NULL, effects = NULL, robust = NULL, df = TRUE, ci = TRUE,
                 transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE,
                 backtransform = NULL){

    ## ** normalizer user input
    options <- LMMstar.options()

    if(!inherits(data,"data.frame")){
        stop("Argument \'data\' must inherit from \"data.frame\". \n")
    }
    if(length(by)!=1){
        stop("Argument \'by\' must has length exactly 1. \n")
    }
    if(by %in% names(data) == FALSE){
        stop("Mismatch between argument \'by\' and \'data\'.\n",
             "Could not find column \"",by,"\" in data \n")
    }
    if(is.null(backtransform)){
        if(is.null(transform.sigma) && is.null(transform.k) && is.null(transform.rho)){
            backtransform <- options$backtransform.confint
        }else{
            backtransform <- FALSE
        }
    }else if(is.character(backtransform)){
        backtransform <-  eval(parse(text=backtransform))
    }
    ## used to decide on the null hypothesis of k parameters
    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = options$transform.sigma, x.transform.k = options$transform.k, x.transform.rho = options$transform.rho)
    
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    transform <- init$transform

    ## ** fit mixed models
    repetition <- list(...)$repetition
    if(!is.null(repetition)){
        res.split <- strsplit(deparse(repetition),"|", fixed = TRUE)[[1]]
        if(length(res.split)==2){
            var.cluster <- trimws(res.split[2], which = "both")
            if(var.cluster %in% names(data) && is.numeric(data[[var.cluster]])){
                if(all(data[[var.cluster]]>0) && all(data[[var.cluster]] %% 1 == 0)){
                    data[[var.cluster]] <- as.factor(sprintf(paste0("%0",ceiling(log10(max(data[[var.cluster]]))+0.1),"d"), data[[var.cluster]]))
                }else{
                    data[[var.cluster]] <- as.factor(data[[var.cluster]])
                }
            }
        }
    }

    ls.data <- base::split(data, data[[by]])
    ls.lmm <- lapply(ls.data, function(iData){
        lmm(..., data = iData, df = df)
    })
    name.lmm <- names(ls.lmm)
    n.lmm <- length(name.lmm)
    
    ## ** test linear combinations
    ls.name.allCoef <- lapply(ls.lmm, function(iLMM){names(coef(iLMM, effects = "all"))})
    ls.Cmat <- lapply(ls.name.allCoef, function(iName){matrix(0, nrow = length(iName), ncol = length(iName), dimnames = list(iName,iName))})

    ## *** identify contrast
    if(is.null(effects)){

        name.contrastCoef <- lapply(ls.lmm, function(iLMM){names(coef(iLMM, effects = options$effects))})

        for(iName in name.lmm){
            iNameCoef <- name.contrastCoef[[iName]]
            if(length(iNameCoef)==1){
                ls.Cmat[[iName]][iNameCoef,iNameCoef] <- 1
            }else{
                diag(ls.Cmat[[iName]][iNameCoef,iNameCoef]) <- 1
            }
            ls.Cmat[[iName]] <- ls.Cmat[[iName]][rowSums(abs(ls.Cmat[[iName]]))!=0,] ## remove lines with only 0
        }
        rhs <- NULL

    }else if(is.matrix(effects)){

        if(length(unique(sapply(ls.name.allCoef,length)))>1){
            stop("Cannot use matrix interface for argument \'effects\' when the number of model parameters varies over splits. \n")
        }
        if(any(sapply(ls.name.allCoef[-1],identical,ls.name.allCoef[[1]])==FALSE)){
            stop("Cannot use matrix interface for argument \'effects\' when the model parameters varies over splits. \n")
        }
        if(NCOL(effects) != length(ls.name.allCoef[[1]])){
            stop("Incorrect value for argument \'effects\'. \n",
                 "If a matrix, it should have as many columns as model parameters (",length(ls.name.allCoef[[1]]),")\n")
        }
        if(!is.null(colnames(effects))){
            if(!identical(sort(colnames(effects),sort(ls.name.allCoef[[1]])))){
                stop("Incorrect value for argument \'effects\'. \n",
                     "If a matrix, its column name should match the model parameters. \n",
                     "Model parameters: \"",paste(ls.name.allCoef[[1]], collapse="\" \""),"\".\n")
            }
            effects <- effects[,ls.name.allCoef[[1]],drop=FALSE]
        }else{
            colnames(effects) <- ls.name.allCoef[[1]]
        }
        ls.Cmat <- stats::setNames(lapply(1:n.lmm, function(iI){effects}), name.lmm)
        
        if(!is.null(attr(effects,"rhs"))){
            rhs <- stats::setNames(lapply(1:n.lmm, function(iI){rhs}), name.lmm)
        }else{
            rhs <- NULL
        }

    }else if(length(effects)==1 && is.character(effects) && grepl("=",effects,fixed=TRUE)){

        ls.Cmat <- stats::setNames(lapply(1:n.lmm, function(iI){effects}), name.lmm)
        rhs <- NULL

    }else if(all(is.character(effects))){
        valid.effects <- c("mean","fixed","variance","correlation","all")
        if(any(effects %in% valid.effects== FALSE)){
            stop("Incorrect value for argument \'effects\'. \n",
                 "When a character should either be a contrast (i.e. contrain the symbol = and have length one), \n",
                 "Or be among: \"",paste(valid.effects,collapse="\" \""),"\"\n")
        }

        name.contrastCoef <- lapply(ls.lmm, function(iLMM){names(coef(iLMM, effects = effects))})

        for(iName in name.lmm){
            iNameCoef <- name.contrastCoef[[iName]]
            if(length(iNameCoef)==1){
                ls.Cmat[[iName]][iNameCoef,iNameCoef] <- 1
            }else{
                diag(ls.Cmat[[iName]][iNameCoef,iNameCoef]) <- 1
            }
            ls.Cmat[[iName]] <- ls.Cmat[[iName]][rowSums(ls.Cmat[[iName]]!=0)>0,,drop=FALSE] 
        }

        
        rhs <- NULL

    }else{
        stop("Unknown value for argument \'effects\'. \n",
             "Can be a matrix, or a character encoding the contrast, or \"mean\", \"variance\", \"correlation\", \"all\".\n")
    }
    ## *** run 
    ls.anova <- stats::setNames(lapply(name.lmm, function(iName){ ## iName <- name.lmm[1]
        if(is.null(robust)){
            anova(ls.lmm[[iName]], effects = ls.Cmat[[iName]], rhs = rhs[[iName]], df = df, ci = ci,
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        }else{
            anova(ls.lmm[[iName]], effects = ls.Cmat[[iName]], rhs = rhs[[iName]], robust = robust, df = df, ci = ci,
                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        }
    }), name.lmm)

    ## ** joint inference
    name.model <- paste0(by,"=",name.lmm)

    if(!is.null(contrast.rbind)){
        rhs.by <- attr(contrast.rbind,"rhs")
        attr(contrast.rbind,"rhs") <- NULL
        sep <- ":"
    }else{
        rhs.by <- NULL
        sep <- ": "
    }
    out <- do.call("rbind.Wald_lmm",
                   args = c(list(model = ls.anova[[1]], effects = contrast.rbind, rhs = rhs.by, name = name.model, sep = sep), unname(ls.anova[-1]))
                   )
    out$model <- ls.lmm
        
    ## ** export
    attr(out,"call") <- match.call()
    class(out) <- append("mlmm", class(out))
    return(out)
}
##----------------------------------------------------------------------
### mlmm.R ends here
