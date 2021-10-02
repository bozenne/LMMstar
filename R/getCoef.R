### getCoef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (14:58) 
## Version: 
## Last-Updated: okt  1 2021 (16:51) 
##           By: Brice Ozenne
##     Update #: 229
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
#' @description Extract all model coefficients with confidence intervals.
#' 
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lmm} object.
#' @param conf.level [numeric 0-1] Confidence level of the confidence intervals.
#' @param effects [character vector] Type of coefficient to be output.
#' Can be coefficients relative to the expectation of the outcome (\code{"mean"} or \code{"fixed"})
#' or to the variance-covariance structure of the residuals (\code{"variance"}).
#' @param format [character] How the output should be shaped.
#' Can be \code{"default"}, \code{"estimate"}, \code{"publish"}, or \code{"SAS"}.
#' @param add.type [logical] Should the type of parameter be added.
#' @param ... argument passed to the \code{publish} function (when \code{format="publish"}).
#' 
#' @details Argument \bold{format}: \cr
#' Setting the argument to \code{"default"} outputs a data.frame with columns type (mean or covariance),
#' term (name of the coefficient), estimate, std.error, t.value, p.value, lower, upper.
#'
#' Setting the argument to \code{"publish"} outputs a data.frame with columns Variable, Units Coefficients, CI, and p-value.
#' Call the function \code{publish} from the \code{publish} package.
#' 
#' Setting the argument to \code{"estimate"} outputs a vector containing the estimated parameter values.
#'                                     
#' Argument \bold{add.type}: \cr
#' When \code{TRUE}, there can be 4 types of parameters in the output: \itemize{
#' \item \code{"mean"}: coefficients relative to the conditional mean of the outcome given the covariates. 
#' \item \code{"std.residual"}: (reference) residual standard deviation.
#' \item \code{"factor.std.residual"}: multiplicative factor to the residual standard deviation.
#' \item \code{"correlation"}: correlation coefficient between the residuals.
#' \item \code{"std.random"}: standard error of the random effects.
#' }
#'
#' @return A data.frame or a vector (see details section)
#' 
#' @examples
#' data(gastricbypassL, package = "LMMstar")
#' library(nlme)
#' 
#' #### linear model ####
#' ## (wrong model as it does not account for repeated measurements)
#' e.lm <- lm(weight ~ time, data = gastricbypassL)
#' 
#' getCoef(e.lm)
#' getCoef(e.lm, format = "estimate")
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

## * getCoef (code)
##' @export
`getCoef` <-
    function(object, conf.level, effects, format, add.type, ...) UseMethod("getCoef")

## * getCoef.lm (code)
##' @export
getCoef.lm <- function(object, conf.level = 0.95, effects = c("mean"),
                       format = "default", add.type = FALSE, ...){
    
    format <- match.arg(format, c("default","estimate","publish", "SAS"))
    if(identical(effects,"all")){
        effects <- c("mean","variance")
    }
    effects <- match.arg(effects, c("mean","fixed","variance"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"
    if(format=="publish"){
        return(Publish::publish(object, ...))
    }
    
    ## **  format=default
    inter <- stats::confint(object, level = conf.level)
    
    out <- NULL
    if("mean" %in% effects){
        objectS <- summary(object, print = FALSE)$coef

        iDF <- data.frame(type = "mean",
                          estimate = as.double(objectS[,"Estimate"]),
                          std.error = as.double(objectS[,"Std. Error"]),
                          t.value = as.double(objectS[,"t value"]),
                          p.value = as.double(objectS[,"Pr(>|t|)"]),
                          lower = as.double(inter[,1]),
                          upper = as.double(inter[,2])
                          )
        rownames(iDF) <- rownames(objectS)
        
        out <- rbind(out,iDF)
    }
    if("variance" %in% effects){
        ddf <- stats::df.residual(object)
        alpha <- 1 - conf.level

        iDF <- data.frame(type = "std.residual",
                          estimate = as.double(stats::sigma(object)),
                          std.error = as.numeric(NA), ## sqrt(2*stats::sigma(object)^4/ddf),
                          t.value = as.numeric(NA),
                          p.value = as.numeric(NA),
                          lower = sqrt(ddf)*stats::sigma(object)/sqrt(stats::qchisq(1 - alpha/2, ddf)),
                          upper = sqrt(ddf)*stats::sigma(object)/sqrt(stats::qchisq(alpha/2, ddf))
                          )
        rownames(iDF) <- "sigma"
        
        out <- rbind(out,iDF)
        if("mean" %in% effects == FALSE){
            out <- out[, c("type","estimate","lower","upper"),drop=FALSE]
        }
    }

    if(format == "SAS"){
        if("variance" %in% effects){
            stop("Argument \'effects\' must be \"mean\" when argument \'format\' is \"SAS\". \n")
        }
        X <- stats::model.matrix(object)
        return(.format2SAS(out, X = X, terms = object$terms))
    }else if(format == "estimate"){
        return(stats::setNames(as.double(out$estimate),rownames(out)))
    }else{
        if(add.type==FALSE){
            out <- out[,setdiff(names(out),"type"),drop=FALSE]
        }
        return(out)
    }
}

## * getCoef.gls (code)
##' @export
getCoef.gls <- function(object, conf.level = 0.95, effects = c("mean"),
                        format = "default", add.type = FALSE, ...){

    format <- match.arg(format, c("default","estimate","publish", "SAS"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(format=="publish"){
        return(Publish::publish(object, ...))
    }

    ## **  format=default
    inter <- .intervalsRobust(object, level = conf.level, effects = effects)

    out <- NULL
    if("mean" %in% effects){
        objectS <- summary(object, print = FALSE)$tTable
        
        iDF <- data.frame(type = "mean",
                          estimate = as.double(inter$coef[,"est."]),
                          std.error = as.double(objectS[,"Std.Error"]),
                          t.value = as.double(objectS[,"t-value"]),
                          p.value = as.double(objectS[,"p-value"]),
                          lower = as.double(inter$coef[,"lower"]),
                          upper = as.double(inter$coef[,"upper"])
                          )
        rownames(iDF) <- rownames(inter$coef)
        
        out <- rbind(out,iDF)
        
    }
    if("variance" %in% effects){
        if(!is.null(object$modelStruct$corStruct)){
            iDF <- data.frame(type = "correlation",
                              estimate = as.double(inter$corStruct[,"est."]),
                              std.error = NA,
                              t.value = NA,
                              p.value = NA,
                              lower = as.double(inter$corStruct[,"lower"]),
                              upper = as.double(inter$corStruct[,"upper"])
                              )
            rownames(iDF) <- rownames(inter$corStruct)
            
            out <- rbind(out,iDF)
        }
        if(!is.null(object$modelStruct$varStruct)){
            iDF <- data.frame(type = "factor.std.residual",
                              estimate = as.double(inter$varStruct[,"est."]),
                              std.error = NA,
                              t.value = NA,
                              p.value = NA,
                              lower = as.double(inter$varStruct[,"lower"]),
                              upper = as.double(inter$varStruct[,"upper"])
                              )
            rownames(iDF) <- rownames(inter$varStruct)
            
            out <- rbind(out,iDF)
        }
        iDF <- data.frame(type = "std.residual",
                          estimate = as.double(inter$sigma["est."]),
                          std.error = NA,
                          t.value = NA,
                          p.value = NA,
                          lower = as.double(inter$sigma["lower"]),
                          upper = as.double(inter$sigma["upper"])
                          )
        rownames(iDF) <- "sigma"
        out <- rbind(out,iDF)
        if("mean" %in% effects == FALSE){
            out <- out[, c("type","estimate","lower","upper"),drop=FALSE]
        }
    }

    if(format == "SAS"){
        X <- stats::model.matrix(object, nlme::getData(object))
        return(.format2SAS(out, X = X, terms = stats::terms(stats::formula(object))))
    }else if(format == "estimate"){
        return(stats::setNames(as.double(out$estimate),rownames(out)))
    }else{
        if(add.type==FALSE){
            out <- out[,setdiff(names(out),"type"),drop=FALSE]
        }
        return(out)
    }
}

## * getCoef.lme (code)
##' @export
getCoef.lme <- function(object, conf.level = 0.95, effects = c("mean"),
                        format = "default", add.type = FALSE, ...){

    format <- match.arg(format, c("default","estimate","publish", "SAS"))
    effects <- match.arg(effects, c("mean","variance"), several.ok = TRUE)
    if(format=="publish"){
        return(Publish::publish(object, ...))
    }

    ## **  format=default
    inter <- .intervalsRobust(object, level = conf.level, effects = effects)

    out <- NULL
    if("mean" %in% effects){
        objectS <- summary(object, print = FALSE)$tTable
        
        iDF <- data.frame(type = "mean",
                          estimate = as.double(inter$fixed[,"est."]),
                          std.error = as.double(objectS[,"Std.Error"]),
                          t.value = as.double(objectS[,"t-value"]),
                          p.value = as.double(objectS[,"p-value"]),
                          lower = as.double(inter$fixed[,"lower"]),
                          upper = as.double(inter$fixed[,"upper"])
                          )
        rownames(iDF) <- rownames(inter$fixed)
        
        out <- rbind(out,iDF)
        
    }
    if("variance" %in% effects){
        if(!is.null(object$modelStruct$corStruct)){
            iDF <- data.frame(type = "correlation",
                              estimate = as.double(inter$corStruct[,"est."]),
                              std.error = NA,
                              t.value = NA,
                              p.value = NA,
                              lower = as.double(inter$corStruct[,"lower"]),
                              upper = as.double(inter$corStruct[,"upper"])
                              )
            rownames(iDF) <- rownames(inter$corStruct)

            out <- rbind(out,iDF)
        }
        if(!is.null(object$modelStruct$varStruct)){
            iDF <- data.frame(type = "factor.std.residual",
                              estimate = as.double(inter$varStruct[,"est."]),
                              std.error = NA,
                              t.value = NA,
                              p.value = NA,
                              lower = as.double(inter$varStruct[,"lower"]),
                              upper = as.double(inter$varStruct[,"upper"])
                              )
            rownames(iDF) <- rownames(inter$varStruct)
            
            out <- rbind(out,iDF)
        }
        for(iTau in 1:length(inter$reStruct)){ ## iTau <- 1
            iNameTau <- names(inter$reStruct)[[iTau]]

            iDF <- data.frame(type = "std.random",
                              estimate = as.double(inter$reStruct[[iTau]]["est."]),
                              std.error = NA,
                              t.value = NA,
                              p.value = NA,
                              lower = as.double(inter$reStruct[[iTau]]["lower"]),
                              upper = as.double(inter$reStruct[[iTau]]["upper"])
                              )
            rownames(iDF) <- paste0(iNameTau,"_",rownames(inter$reStruct[[iTau]]))
            
            out <- rbind(out,iDF)
        }
        iDF <- data.frame(type = "std.residual",
                          estimate = as.double(inter$sigma["est."]),
                          std.error = NA,
                          t.value = NA,
                          p.value = NA,
                          lower = as.double(inter$sigma["lower"]),
                          upper = as.double(inter$sigma["upper"])
                          )
        rownames(iDF) <- "sigma"
        
        out <- rbind(out,iDF)
        
        if("mean" %in% effects == FALSE){
            out <- out[, c("type","estimate","lower","upper"),drop=FALSE]
        }
    }
    
    if(format == "SAS"){
        X <- stats::model.matrix(object, nlme::getData(object))
        return(.format2SAS(out, X = X, terms = stats::terms(stats::formula(object))))
    }else if(format == "estimate"){
        return(stats::setNames(as.double(out$estimate),rownames(out)))
    }else{
        if(add.type==FALSE){
            out <- out[,setdiff(names(out),"type"),drop=FALSE]
        }
        return(out)
    }
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
    out$t.value <- object$t.value
    out$p.value <- object$p.value
    out$lower <- object$lower
    out$upper <- object$upper
    return(out)     
}


## * .intervalsRobust
## intervals function but handle non-invertible vcov (i.e. only return point estimates)
.intervalsRobust <- function(object, level, effects){
    if(inherits(object,"gls")){
        inter <- intervals(object, level = level, which = "coef")
    }else if(inherits(object,"lme")){
        inter <- intervals(object, level = level, which = "fixed")
    }

    if("variance" %in% effects){
        tempo <- try(intervals(object, level = level, which = "var-cov"), silent = TRUE)

        if(inherits(tempo,"try-error")){
            if(inherits(object$modelStruct$corStruct,"corSymm") && length(object$modelStruct$reStruct)<=1){
            tempo <- list()

            if(!is.null(object$modelStruct$reStruct)){
                plen <- names(attr(object$modelStruct$reStruct, "plen"))
                ## undebug(nlme:::intervals.lme)
                tempo$reStruct <- lapply(object$modelStruct$reStruct, function(iStruct){
                    cbind(lower = NA, "est." = stats::coef(nlme::pdNatural(iStruct[[1]]), unconstrained = FALSE), "upper" = NA)
                })
                names(tempo$reStruct) <- names(object$modelStruct$reStruct)
            }

           
            mC <- attr(object$modelStruct$corStruct,"maxCov")
            M.index <- which(lower.tri(diag(mC)), arr.ind = TRUE)
            tempo$corStruct <- cbind(lower = NA, "est." = stats::coef(object$modelStruct$corStruct, unconstrained = FALSE), "upper" = NA)
            rownames(tempo$corStruct) <- paste0("cor(",M.index[,2],",",M.index[,1],")")

            if(!is.null(object$modelStruct$varStruct)){
                tempo$varStruct <- cbind(lower = NA, "est." = stats::coef(object$modelStruct$varStruct, unconstrained = FALSE), "upper" = NA)
            }
            
            tempo$sigma <- c("lower" = NA,
                             "est." = stats::sigma(object),
                             "upper" = NA)


            }else{
                stop(tempo)
            }
        }
        
        inter <- c(inter,tempo)
    }

    return(inter)
}


######################################################################
### getCoef.R ends here
