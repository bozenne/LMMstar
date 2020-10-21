### tidy.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (14:58) 
## Version: 
## Last-Updated: okt 21 2020 (15:41) 
##           By: Brice Ozenne
##     Update #: 31
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * tidy.lmm (documentation) 
#' @title Extract Model Coefficients With Confidence Intervals
#' 
#' @description Extract model coefficients with confidence intervals.
#' @param x lmm object.
#' @param conf.int [logical] Should the confidence interval be output
#' @param conf.level [numeric 0-1] Confidence level of the confidence intervals.
#' @param effects [character vector] Type of coefficient to be output.
#' Can be coefficients relative to the expectation of the outcome (\code{"mean"})
#' or to the variance-covariance structure of the residuals (\code{"variance"}).
#' @param ... ignored. For compatibility with the generic method.
#'
#' @examples
#' library(broom)
#' data(gastricbypassL)
#' e.lmm <- lmm(weight ~ time, covariance = ~visit|id, data = gastricbypassL)
#' tidy(e.lmm)
#' tidy(e.lmm, effects = c("mean","variance"))

## * tidy.lmm (code) 
#' @export
tidy.lmm <- function(x, conf.int = TRUE, conf.level = 0.95, effects = c("mean"), ...){
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    inter <- intervals(x, level = conf.level)

    out <- NULL
    if("mean" %in% effects){
        xS <- summary(x, print = FALSE)[["mean"]]
        out <- rbind(out,data.frame(type = "mean",
                                    term = rownames(inter$coef),
                                    estimate = as.double(inter$coef[,"est."]),
                                    std.error = as.double(xS[,"se"]),
                                    statistic = as.double(xS[,"t-value"]),
                                    p.value = as.double(xS[,"p-value"]),
                                    conf.low = as.double(inter$coef[,"lower"]),
                                    conf.high = as.double(inter$coef[,"upper"])
                                    ))
        
    }
    if("variance" %in% effects){
        out <- rbind(out,data.frame(type = "corStruct",
                                    term = rownames(inter$corStruct),
                                    estimate = as.double(inter$corStruct[,"est."]),
                                    std.error = NA,
                                    statistic = NA,
                                    p.value = NA,
                                    conf.low = as.double(inter$corStruct[,"lower"]),
                                    conf.high = as.double(inter$corStruct[,"upper"])
                                    ))
        out <- rbind(out,data.frame(type = "varStruct",
                                  term = rownames(inter$varStruct),
                                  estimate = as.double(inter$varStruct[,"est."]),
                                  std.error = NA,
                                  statistic = NA,
                                  p.value = NA,
                                  conf.low = as.double(inter$varStruct[,"lower"]),
                                  conf.high = as.double(inter$varStruct[,"upper"])
                                  ))
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
