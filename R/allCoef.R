### tidy.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (14:58) 
## Version: 
## Last-Updated: okt 24 2020 (14:52) 
##           By: Brice Ozenne
##     Update #: 58
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * allCoef (documentation) 
#' @title Extract Model Coefficients With Confidence Intervals
#' 
#' @description Extract all model coefficients with confidence intervals.
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lmm} object.
#' @param conf.int [logical] Should the confidence interval be output
#' @param conf.level [numeric 0-1] Confidence level of the confidence intervals.
#' @param effects [character vector] Type of coefficient to be output.
#' Can be coefficients relative to the expectation of the outcome (\code{"mean"})
#' or to the variance-covariance structure of the residuals (\code{"variance"}).
#'
#' @examples
#' data(gastricbypassL)
#' library(nlme)
#' 
#' #### linear model ####
#' ## (wrong model as it does not account for repeated measurements)
#' e.lm <- lm(weight ~ time, data = gastricbypassL)
#' allCoef(e.lm)
#'
#' #### gls model ####
#' e.gls <- gls(weight ~ time,
#'              correlation = corSymm(form =~as.numeric(visit)|id),
#'              weights = varIdent(form =~1|visit),
#'              data = gastricbypassL)
#' allCoef(e.gls)
#' 
#' #### lme model ####
#' e.lme <- lme(weight ~ time,
#'              random = ~1|id,
#'              weights = varIdent(form =~1|visit),
#'              data = gastricbypassL)
#' allCoef(e.lme)
#'
#' #### lmm model ####
#' e.lmm <- lmm(weight ~ time, covariance = ~visit|id, data = gastricbypassL)
#' allCoef(e.lmm)

## * allCoef (code)
##' @export
`allCoef` <-
    function(object, conf.int, conf.level, effects) UseMethod("allCoef")

## * allCoef.lm (code)
##' @export
allCoef.lm <- function(object, conf.int = TRUE, conf.level = 0.95, effects = c("mean","variance")){
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    inter <- stats::confint(object, level = conf.level)

    out <- NULL
    if("mean" %in% effects){
        objectS <- summary(object, print = FALSE)$coef
        out <- rbind(out,data.frame(type = "mean",
                                    term = rownames(objectS),
                                    estimate = as.double(objectS[,"Estimate"]),
                                    std.error = as.double(objectS[,"Std. Error"]),
                                    statistic = as.double(objectS[,"t value"]),
                                    p.value = as.double(objectS[,"Pr(>|t|)"]),
                                    conf.low = as.double(inter[,1]),
                                    conf.high = as.double(inter[,2])
                                    ))
    }
    if("variance" %in% effects){
        out <- rbind(out,data.frame(type = "sigma",
                                    term = "sigma",
                                    estimate = as.double(stats::sigma(object)),
                                    std.error = as.numeric(NA),
                                    statistic = as.numeric(NA),
                                    p.value = as.numeric(NA),
                                    conf.low = as.numeric(NA),
                                    conf.high = as.numeric(NA)
                                    ))
    }
    return(out)
}

## * allCoef.gls (code)
##' @export
allCoef.gls <- function(object, conf.int = TRUE, conf.level = 0.95, effects = c("mean","variance")){
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    inter <- intervals(object, level = conf.level)

    out <- NULL
    if("mean" %in% effects){
        objectS <- summary(object, print = FALSE)$tTable
        out <- rbind(out,data.frame(type = "mean",
                                    term = rownames(inter$coef),
                                    estimate = as.double(inter$coef[,"est."]),
                                    std.error = as.double(objectS[,"Std.Error"]),
                                    statistic = as.double(objectS[,"t-value"]),
                                    p.value = as.double(objectS[,"p-value"]),
                                    conf.low = as.double(inter$coef[,"lower"]),
                                    conf.high = as.double(inter$coef[,"upper"])
                                    ))
        
    }
    if("variance" %in% effects){
        if(!is.null(object$modelStruct$corStruct)){
            out <- rbind(out,data.frame(type = "corStruct",
                                        term = rownames(inter$corStruct),
                                        estimate = as.double(inter$corStruct[,"est."]),
                                        std.error = NA,
                                        statistic = NA,
                                        p.value = NA,
                                        conf.low = as.double(inter$corStruct[,"lower"]),
                                        conf.high = as.double(inter$corStruct[,"upper"])
                                        ))
        }
        if(!is.null(object$modelStruct$varStruct)){
            out <- rbind(out,data.frame(type = "varStruct",
                                        term = rownames(inter$varStruct),
                                        estimate = as.double(inter$varStruct[,"est."]),
                                        std.error = NA,
                                        statistic = NA,
                                        p.value = NA,
                                        conf.low = as.double(inter$varStruct[,"lower"]),
                                        conf.high = as.double(inter$varStruct[,"upper"])
                                        ))
        }
        out <- rbind(out,data.frame(type = "sigma",
                                    term = "sigma",
                                    estimate = as.double(inter$sigma["est."]),
                                    std.error = NA,
                                    statistic = NA,
                                    p.value = NA,
                                    conf.low = as.double(inter$sigma["lower"]),
                                    conf.high = as.double(inter$sigma["upper"])
                                    ))
    }
    
    return(out)
}

## * summarize.lme (code)
##' @export
allCoef.lme <- function(object, conf.int = TRUE, conf.level = 0.95, effects = c("mean","variance")){
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    inter <- intervals(object, level = conf.level)

    out <- NULL
    if("mean" %in% effects){
        objectS <- summary(object, print = FALSE)$tTable
        out <- rbind(out,data.frame(type = "mean",
                                    term = rownames(inter$fixed),
                                    estimate = as.double(inter$fixed[,"est."]),
                                    std.error = as.double(objectS[,"Std.Error"]),
                                    statistic = as.double(objectS[,"t-value"]),
                                    p.value = as.double(objectS[,"p-value"]),
                                    conf.low = as.double(inter$fixed[,"lower"]),
                                    conf.high = as.double(inter$fixed[,"upper"])
                                    ))
        
    }
    if("variance" %in% effects){
        if(!is.null(object$modelStruct$corStruct)){
            out <- rbind(out,data.frame(type = "corStruct",
                                        term = rownames(inter$corStruct),
                                        estimate = as.double(inter$corStruct[,"est."]),
                                        std.error = NA,
                                        statistic = NA,
                                        p.value = NA,
                                        conf.low = as.double(inter$corStruct[,"lower"]),
                                        conf.high = as.double(inter$corStruct[,"upper"])
                                        ))
        }
        if(!is.null(object$modelStruct$varStruct)){
            out <- rbind(out,data.frame(type = "varStruct",
                                        term = rownames(inter$varStruct),
                                        estimate = as.double(inter$varStruct[,"est."]),
                                        std.error = NA,
                                        statistic = NA,
                                        p.value = NA,
                                        conf.low = as.double(inter$varStruct[,"lower"]),
                                        conf.high = as.double(inter$varStruct[,"upper"])
                                        ))
        }
        for(iTau in 1:length(inter$reStruct)){ ## iTau <- 1
            iNameTau <- names(inter$reStruct)[[iTau]]
            out <- rbind(out,data.frame(type = "random",
                                        term = paste0(iNameTau,"_",rownames(inter$reStruct[[iTau]])),
                                        estimate = as.double(inter$reStruct[[iTau]]["est."]),
                                        std.error = NA,
                                        statistic = NA,
                                        p.value = NA,
                                        conf.low = as.double(inter$reStruct[[iTau]]["lower"]),
                                        conf.high = as.double(inter$reStruct[[iTau]]["upper"])
                                        ))
        }
        out <- rbind(out,data.frame(type = "sigma",
                                    term = "sigma",
                                    estimate = as.double(inter$sigma["est."]),
                                    std.error = NA,
                                    statistic = NA,
                                    p.value = NA,
                                    conf.low = as.double(inter$sigma["lower"]),
                                    conf.high = as.double(inter$sigma["upper"])
                                    ))
    }
    
    return(out)
}

## * summarize.lmm (code)
##' @export
allCoef.lmm <- function(object, conf.int = TRUE, conf.level = 0.95, effects = c("mean")){
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    inter <- intervals(object, level = conf.level)

    out <- NULL
    if("mean" %in% effects){
        objectS <- summary(object, print = FALSE)[["mean"]]
        out <- rbind(out,data.frame(type = "mean",
                                    term = rownames(inter$coef),
                                    estimate = as.double(inter$coef[,"est."]),
                                    std.error = as.double(objectS[,"se"]),
                                    statistic = as.double(objectS[,"t-value"]),
                                    p.value = as.double(objectS[,"p-value"]),
                                    conf.low = as.double(inter$coef[,"lower"]),
                                    conf.high = as.double(inter$coef[,"upper"])
                                    ))
        
    }
    if("variance" %in% effects){
        if(!is.null(object$modelStruct$corStruct)){
            out <- rbind(out,data.frame(type = "corStruct",
                                        term = rownames(inter$corStruct),
                                        estimate = as.double(inter$corStruct[,"est."]),
                                        std.error = NA,
                                        statistic = NA,
                                        p.value = NA,
                                        conf.low = as.double(inter$corStruct[,"lower"]),
                                        conf.high = as.double(inter$corStruct[,"upper"])
                                        ))
        }
        if(!is.null(object$modelStruct$varStruct)){
            out <- rbind(out,data.frame(type = "varStruct",
                                        term = rownames(inter$varStruct),
                                        estimate = as.double(inter$varStruct[,"est."]),
                                        std.error = NA,
                                        statistic = NA,
                                        p.value = NA,
                                        conf.low = as.double(inter$varStruct[,"lower"]),
                                        conf.high = as.double(inter$varStruct[,"upper"])
                                        ))
        }
        out <- rbind(out,data.frame(type = "sigma",
                                    term = "sigma",
                                    estimate = as.double(inter$sigma["est."]),
                                    std.error = NA,
                                    statistic = NA,
                                    p.value = NA,
                                    conf.low = as.double(inter$sigma["lower"]),
                                    conf.high = as.double(inter$sigma["upper"])
                                    ))
    }
    
    return(out)
}

######################################################################
### tidy.lmm.R ends here
