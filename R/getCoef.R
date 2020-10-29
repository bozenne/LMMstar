### tidy.lmm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (14:58) 
## Version: 
## Last-Updated: okt 29 2020 (15:28) 
##           By: Brice Ozenne
##     Update #: 159
##----------------------------------------------------------------------
## 
### Commentary: 
##
##
##
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getCoef (documentation) 
#' @title Extract Model Coefficients With Confidence Intervals
#' 
#' @description Extract all model coefficients with confidence intervals.
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lmm} object.
#' @param conf.int [logical] Should the confidence interval be output
#' @param conf.level [numeric 0-1] Confidence level of the confidence intervals.
#' @param effects [character vector] Type of coefficient to be output.
#' Can be coefficients relative to the expectation of the outcome (\code{"mean"})
#' or to the variance-covariance structure of the residuals (\code{"variance"}).
#' @param format [character] How the output should be shaped.
#' Can be \code{"default"}, \code{"publish"}, or \code{"SAS"}.
#' @param ... argument passed to the \code{publish} function (when \code{format="publish"}).
#'
#' 
#' @details Arugment \bold{format}: \cr
#' Setting the argument to \code{"default"} outputs a data.frame with columns type (mean or covariance),
#' term (name of the coefficient), estimate, std.error, statistic, p.value, conf.low, conf.high.
#'
#' Setting the argument to \code{"publish"} outputs a data.frame with columns Variable, Units Coefficients, CI, and p-value.
#' Call the function \code{publish} from the \code{publish} package.
#'                                     
#'
#' @examples
#' data(gastricbypassL, package = "repeated")
#' library(nlme)
#' 
#' #### linear model ####
#' ## (wrong model as it does not account for repeated measurements)
#' e.lm <- lm(weight ~ time, data = gastricbypassL)
#' 
#' getCoef(e.lm)
#' getCoef(e.lm, effects = "variance")
#' getCoef(e.lm, effects = "variance", format = "estimate")
#' if(require(Publish)){
#' getCoef(e.lm, format = "publish")
#' }
#' getCoef(e.lm, format = "SAS")
#'
#' #### gls model ####
#' e.gls <- gls(weight ~ time,
#'              correlation = corSymm(form =~as.numeric(visit)|id),
#'              weights = varIdent(form =~1|visit),
#'              data = gastricbypassL)
#' getCoef(e.gls)
#' getCoef(e.gls, effects = "variance")
#' getCoef(e.gls, effects = "variance", format = "estimate")
#' if(require(Publish)){
#' getCoef(e.gls, format = "publish")
#' }
#' getCoef(e.gls, format = "SAS")
#' 
#' #### lme model ####
#' e.lme <- lme(weight ~ time,
#'              random = ~1|id,
#'              weights = varIdent(form =~1|visit),
#'              data = gastricbypassL)
#' getCoef(e.lme)
#' getCoef(e.lme, effects = "variance")
#' getCoef(e.lme, effects = "variance", format = "estimate")
#' if(require(Publish)){
#' getCoef(e.lme, format = "publish")
#' }
#' getCoef(e.lme, format = "SAS")
#' 
#' #### lmm model ####
#' e.lmm <- lmm(weight ~ time, covariance = ~visit|id, data = gastricbypassL)
#' getCoef(e.lmm)
#' getCoef(e.lmm, effects = "variance")
#' getCoef(e.lmm, effects = "variance", format = "estimate")
#' if(require(Publish)){
#' getCoef(e.lmm, format = "publish")
#' }
#' getCoef(e.lmm, format = "SAS")

## * getCoef (code)
##' @export
`getCoef` <-
    function(object, conf.int, conf.level, effects, format, ...) UseMethod("getCoef")

## * getCoef.lm (code)
##' @export
getCoef.lm <- function(object, conf.int = TRUE, conf.level = 0.95, effects = c("mean"),
                       format = "default", ...){
    
    format <- match.arg(format, c("default","estimate","publish", "SAS"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(format=="publish"){
        return(Publish::publish(object, ...))
    }
    
    ## **  format=default
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
        ddf <- stats::df.residual(object)
        alpha <- 1 - conf.level

        out <- rbind(out,data.frame(type = "sigma",
                                    term = "sigma",
                                    estimate = as.double(stats::sigma(object)),
                                    std.error = as.numeric(NA), ## sqrt(2*stats::sigma(object)^4/ddf),
                                    statistic = as.numeric(NA),
                                    p.value = as.numeric(NA),
                                    conf.low = sqrt(ddf)*stats::sigma(object)/sqrt(stats::qchisq(1 - alpha/2, ddf)),
                                    conf.high = sqrt(ddf)*stats::sigma(object)/sqrt(stats::qchisq(alpha/2, ddf))
                                    ))
        if("mean" %in% effects == FALSE){
            out <- out[, c("type","term","estimate","conf.low","conf.high"),drop=FALSE]
        }
    }

    if(format == "SAS"){
        if("variance" %in% effects){
            stop("Argument \'effects\' must be \"mean\" when argument \'format\' is \"SAS\". \n")
        }
        X <- stats::model.matrix(object)
        return(.format2SAS(out, X = X, terms = object$terms))
    }else if(format == "estimate"){
        return(out[, c("type","term","estimate"),drop=FALSE])
    }else{
        return(out)
    }
}

## * getCoef.gls (code)
##' @export
getCoef.gls <- function(object, conf.int = TRUE, conf.level = 0.95, effects = c("mean"),
                        format = "default", ...){

    format <- match.arg(format, c("default","estimate","publish", "SAS"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(format=="publish"){
        return(Publish::publish(object, ...))
    }

    ## **  format=default
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
        if("mean" %in% effects == FALSE){
            out <- out[, c("type","term","estimate","conf.low","conf.high"),drop=FALSE]
        }
    }

    if(format == "SAS"){
        X <- stats::model.matrix(object, nlme::getData(object))
        return(.format2SAS(out, X = X, terms = stats::terms(stats::formula(object))))
    }else if(format == "estimate"){
        return(out[, c("type","term","estimate"),drop=FALSE])
    }else{
        return(out)
    }
}

## * summarize.lme (code)
##' @export
getCoef.lme <- function(object, conf.int = TRUE, conf.level = 0.95, effects = c("mean"),
                        format = "default", ...){

    format <- match.arg(format, c("default","estimate","publish", "SAS"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(format=="publish"){
        return(Publish::publish(object, ...))
    }

    ## **  format=default
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
        
        if("mean" %in% effects == FALSE){
            out <- out[, c("type","term","estimate","conf.low","conf.high"),drop=FALSE]
        }
    }
    
    if(format == "SAS"){
        X <- stats::model.matrix(object, nlme::getData(object))
        return(.format2SAS(out, X = X, terms = stats::terms(stats::formula(object))))
    }else if(format == "estimate"){
        return(out[, c("type","term","estimate"),drop=FALSE])
    }else{
        return(out)
    }
}

## * summarize.lmm (code)
##' @export
getCoef.lmm <- function(object, conf.int = TRUE, conf.level = 0.95, effects = c("mean"),
                        format = "default", ...){

    format <- match.arg(format, c("default","estimate","publish", "SAS"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(format=="publish"){
        object2 <- object
        class(object2) <- setdiff(class(object),"lmm")
        return(Publish::publish(object2, ...))
    }

    ## **  format=default
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
        if("mean" %in% effects == FALSE){
            out <- out[, c("type","term","estimate","conf.low","conf.high"),drop=FALSE]
        }

    }
    
    if(format == "SAS"){
        X <- stats::model.matrix(object, attr(object,"data"))
        return(.format2SAS(out, X = X, terms = stats::terms(stats::formula(object))))
    }else if(format == "estimate"){
        return(out[, c("type","term","estimate"),drop=FALSE])
    }else{
        return(out)
    }
    return(out)
}

## * .format2SAS
.format2SAS <- function(object, X, terms){

    ## ** check arguments
    if(is.null(assign)){
        stop("Argument \'X\' must contain an attribute \"assign\" \n")
    }

    ## ** initialize
    assign <- attr(X,"assign")
    terms.labels <- attr(terms,"term.labels")
    out <- data.frame(Effect = c("(Intercept)",terms.labels)[assign+1])
    
    ## ** add factors
    if(!is.null(attr(X,"contrasts"))){    
        ## type.X <-  attr(terms,"dataClasses")
        name.factor <- names(attr(X,"contrasts"))
        n.factor <- length(name.factor)

        Mfactor <- attr(terms,"factor")[name.factor,,drop=FALSE]
        terms.labels.factor <- names(which(colSums(Mfactor)>0))


        for(iFactor in 1:n.factor){## iFactor <- 1

            out[[name.factor[iFactor]]] <- ""
            iValue <- sapply(strsplit(colnames(X)[assign==iFactor],paste0("^",name.factor[iFactor])),"[[",2)
            
            for(iLabel in terms.labels){ ## iLabel <- "time"
                if(Mfactor[name.factor[iFactor],iLabel]>0){
                    out[out$Effect==iLabel,name.factor[iFactor]] <- iValue
                }
            }
        }
    }

    ## ** add estimates
    out$estimate <- object$estimate
    out$std.error <- object$std.error
    out$statistic <- object$statistic
    out$p.value <- object$p.value
    out$conf.low <- object$conf.low
    out$conf.high <- object$conf.high
    return(out)     
}

######################################################################
### tidy.lmm.R ends here
