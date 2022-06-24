### mlmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 14 2022 (09:45) 
## Version: 
## Last-Updated: jun 24 2022 (18:03) 
##           By: Brice Ozenne
##     Update #: 61
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
##' @examples
##'
##' #### univariate regression ####
##' if(require(lava)){
##' library(LMMstar)
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
mlmm <- function(..., data, by, effects = NULL, robust = FALSE, df = TRUE, ci = TRUE){

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

    ## ** fit mixed models
    ls.data <- base::split(data, data[[by]])
    ls.lmm <- lapply(ls.data, function(iData){
        lmm(..., data = iData, df = df)
    })

    if(is.null(effects) || all(effects %in% c("mean","fixed","variance","correlation","all"))){
        ls.anova <- lapply(ls.lmm, function(iLMM){
            iAllCoef <- names(coef(iLMM, effects = "all"))
            iAllCoef.effects <- names(coef(iLMM, effects = options$effects))
            iC <- matrix(0, nrow = length(iAllCoef.effects), ncol = length(iAllCoef), dimnames = list(iAllCoef.effects,iAllCoef))
            if(length(iAllCoef.effects)==1){
                iC[iAllCoef.effects,iAllCoef.effects] <- 1
            }else{
                diag(iC[iAllCoef.effects,iAllCoef.effects]) <- 1
            }
            anova(iLMM, effects = iC, robust = robust, df = df, ci = ci)
        })
    }else{
        ls.anova <- lapply(ls.lmm, anova, effects = effects, robust = robust, df = df, ci = ci)
    }

    ## ** joint inference
    ## out <- do.call("rbind", args = c(list(model = ls.anova[[1]]), unname(ls.anova[-1])))
    name.model <- paste0(by,"=",unlist(lapply(ls.data, function(iData){iData[[by]][1]})))
    out <- do.call("rbind.anova_lmm",
                   args = c(list(model = ls.anova[[1]], name = name.model), unname(ls.anova[-1]))
                   )
    attr.callout <- list(df = attr(out,"df"),
                         test = attr(out,"test"),
                         robust = attr(out,"robust"))
    out$call <- match.call()
    attr(out$call,"df") <- attr.callout$df
    attr(out$call,"test") <- attr.callout$test
    attr(out$call,"robust") <- attr.callout$robust
    class(out) <- append("mlmm", class(out))
    return(out)
}
##----------------------------------------------------------------------
### mlmm.R ends here
