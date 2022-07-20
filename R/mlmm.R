### mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 14 2022 (09:45) 
## Version: 
## Last-Updated: jul 20 2022 (15:13) 
##           By: Brice Ozenne
##     Update #: 104
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
##' @param effects [character] Linear combinations of coefficients relative to which Wald test should be computed.
##' @param robust [logical] Should robust standard errors (aka sandwich estimator) be output instead of the model-based standard errors. Argument passed to \code{anova.lmm}.
##' @param df [logical] Should the degree of freedom be computed using a Satterthwaite approximation?
##' @param ci [logical] Should a confidence interval be output for each hypothesis?
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
mlmm <- function(..., data, by, effects = NULL, robust = NULL, df = TRUE, ci = TRUE,
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

    if(is.null(effects) || all(effects %in% c("mean","fixed","variance","correlation","all"))){
        ls.anova <- lapply(ls.lmm, function(iLMM){ ## iLMM <- ls.lmm[[1]]
            iAllCoef <- names(coef(iLMM, effects = "all"))
            if(is.null(effects)){
                iAllCoef.effects <- names(coef(iLMM, effects = options$effects))
            }else{
                iAllCoef.effects <- names(coef(iLMM, effects = effects))
            }
            iC <- matrix(0, nrow = length(iAllCoef.effects), ncol = length(iAllCoef), dimnames = list(iAllCoef.effects,iAllCoef))
            if(length(iAllCoef.effects)==1){
                iC[iAllCoef.effects,iAllCoef.effects] <- 1
            }else{
                diag(iC[iAllCoef.effects,iAllCoef.effects]) <- 1
            }
            if(is.null(robust)){
                iOut <- anova(iLMM, effects = iC, df = df, ci = ci,
                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
            }else{
                iOut <- anova(iLMM, effects = iC, robust = robust, df = df, ci = ci,
                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
            }
            return(iOut)
        })
    }else{
        if(is.null(robust)){
            ls.anova <- lapply(ls.lmm, anova, effects = effects, df = df, ci = ci,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        }else{
            ls.anova <- lapply(ls.lmm, anova, effects = effects, robust = robust, df = df, ci = ci,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)
        }
    }

    ## ** joint inference
    name.model <- paste0(by,"=",unlist(lapply(ls.data, function(iData){iData[[by]][1]})))
    contrast <- attr(effects,"contrasts")
    rhs <- attr(effects,"rhs")
    out <- do.call("rbind.anova_lmm",
                   args = c(list(model = ls.anova[[1]], effects = contrast, rhs = rhs, name = name.model), unname(ls.anova[-1]))
                   )
    out$model <- ls.lmm
        
    ## ** export
    attr(out,"call") <- match.call()
    class(out) <- append("mlmm", class(out))
    return(out)
}
##----------------------------------------------------------------------
### mlmm.R ends here
