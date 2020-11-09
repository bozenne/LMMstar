### ufit.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 21 2020 (15:47) 
## Version: 
## Last-Updated: nov  9 2020 (11:43) 
##           By: Brice Ozenne
##     Update #: 103
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ufit (documentation)
#' @title Unique Fitted Values
#' 
#' @description Compute all the possible fitted values among the possible combinaisons of categorical covariables.
#' Continuous covariates are set to a single value.
#' 
#' @param object a \code{lmm} or \code{gls} object
#' @param value [named numeric vector] values used to fix the continuous covariates.
#' @param confint [logical] should confidence intervals for the fitted values be computed?
#' @param conf.quantile [numeric 0-1] the quantile used to compute the confidence intervals.

## * ufit (examples)
#' @examples
#' library(ggplot2)
#' library(data.table)
#' library(nlme)
#'
#' data(gastricbypassL, package = "repeated")
#' ## generate covariates
#' gastricbypassL <- as.data.table(gastricbypassL)
#' gastricbypassL[, baselineG := .SD$glucagon[1]<10000, by = "id"]
#' gastricbypassL[, age := round(rnorm(1, mean = 50, sd = 10)), by = "id"]
#' 
#' ## time evolution
#' e.lmm1 <- lmm(weight ~ time, covariance = ~visit|id, data = gastricbypassL)
#' e.ufit1 <- ufit(e.lmm1)
#' autoplot(e.ufit1)
#' 
#' e.gls2 <- gls(glucagon ~ time, correlation = corSymm(form=~as.numeric(visit)|id), weights = varIdent(form=~1|visit), data = gastricbypassL, na.action = na.exclude)
#' e.ufit2 <- ufit(e.gls2)
#' e.ufit2
#' 
#' e.gls3 <- gls(glucagon ~ time, correlation = corCompSymm(form=~1|id), data = gastricbypassL, na.action = na.exclude)
#' e.ufit3 <- ufit(e.gls3)
#' e.ufit3
#' 
#' ## time per group
#' e.lmm2 <- lmm(weight ~ time*baselineG, covariance = ~visit|id, data = gastricbypassL)
#' e.ufit2 <- ufit(e.lmm2)
#' autoplot(e.ufit2) + theme(axis.text.x = element_text(angle = 90))
#' autoplot(e.ufit2, x = "visit")
#' 
#' ## adding continuous covariates
#' e.lmm3 <- lmm(weight ~ time*baselineG + age, covariance = ~visit|id, data = gastricbypassL)
#' e.ufit45 <- ufit(e.lmm3, c("age" = 30))
#' e.ufit50 <- ufit(e.lmm3, c("age" = 60))
#' if(require(ggpubr)){
#' ggarrange(autoplot(e.ufit45, x = "visit") + coord_cartesian(ylim = c(80,150)),
#'           autoplot(e.ufit50, x = "visit") + coord_cartesian(ylim = c(80,150)),
#'           legend = "bottom")
#' }

## * ufit (code)
#' @export
ufit <- function(object, value = NULL, confint = TRUE, conf.quantile = stats::qnorm(0.975)){

    ## ** normalize arguments
    ff <- stats::formula(object)
    if(any(all.vars(ff) == "fit")){
        stop("Do not work when one of the variables in the model formula is called \"fit\". \n")
    }
    if(confint && any(all.vars(ff) %in% c("se","lower","upper"))){
        stop("Do not work when one of the variables in the model formula is called \"se\", \"lower\", or \"upper\". \n")
    }

    ## ** Update dataset according to argument values (i.e. neutralize continuous variables)
    newdata <- as.data.frame(getData(object))
    X <- stats::model.matrix(stats::terms(object), data = newdata)
    test.factor <- colSums((X == 0)+(X == 1))==NROW(X)

    name.contvar <- names(which(test.factor==FALSE))
    if(any(names(value) %in% name.contvar == FALSE)){
        stop("Incorrect variable in argument \'value\', must only contain values for the continuous covariates.\n",
             "Incorrect variable in argument \'value\': \"",paste0(setdiff(names(value),name.contvar), collapse = "\" \""),"\"\n")
    }
    if(!is.null(value)){
        if(is.null(names(value))){
            stop("Argument \'value\' should be named using the covariate names \n")
        }
        for(iContVar in names(value)){
            newdata[,iContVar] <- value[iContVar]
        }
        X <- stats::model.matrix(stats::terms(object), data = newdata)
    }
    test.continuous <- apply(X,2,function(x){length(unique(x))})
    
    if(any(test.continuous>2)){
        if(is.null(value)){
            stop("In presence of continuous covariates, the argument \'value\' must be specified.\n",
                 "It must be a vector containing one value for each continuous covariate, named with the covariate names. \n")
        }else{
            stop("In presence of continuous covariates, the argument \'value\' must contain one value for each continuous covariate. \n",
                 "Missing continuous variable in argument \'value\': \"",paste0(names(test.continuous[test.continuous>2]), collapse = "\" \""),"\"\n")
        }
    }
    index.unique <- which(!duplicated(X))
    if(length(object$na.action)!=0){
        index.unique <- (1:NROW(newdata))[-object$na.action][index.unique]
    }
    Unewdata <- newdata[index.unique,,drop=FALSE]

    ## ** compute predictions with confidence intervals
    if(confint){
        requireNamespace("AICcmodavg")
        tempo <- AICcmodavg::predictSE(object, newdata = Unewdata, se.fit = TRUE)
        Unewdata$fit <- tempo$fit
        Unewdata$se <- tempo$se.fit
        Unewdata$lower <- tempo$fit - conf.quantile * tempo$se.fit
        Unewdata$upper <- tempo$fit + conf.quantile * tempo$se.fit
    }else{
        Unewdata$fit <- stats::predict(object, newdata = Unewdata)
    }

    ## ** find visit variable
    if(!is.null(object$modelStruct$corStruct)){
        cor.var <- all.vars(stats::formula(object$modelStruct$corStruct))
        id.var <- utils::tail(all.vars(stats::formula(object$modelStruct$corStruct)),1)
        if(length(cor.var)==2){
            time.var <- cor.var[1]
        }else{
            Uid.var <- unique(newdata[[id.var]])

            if(!is.null(object$modelStruct$varStruct)){
                possible.var <- rhs.vars(stats::formula(object$modelStruct$varStruct))
                test.duplicated <- colSums(do.call(rbind,lapply(Uid.var, function(iId){apply(newdata[newdata[[id.var]]==iId,possible.var,drop=FALSE],2,duplicated)})))
                if(sum(test.duplicated==0)==1){
                    time.var <- as.character(possible.var[test.duplicated==0])
                }else{
                    time.var <- NULL
                }
            }else{
                time.var <- NULL
            }
        }
    }else{
        id.var <- NULL
        time.var <- NULL
    }
    
    ## ** export
    outcome <- as.character(lhs.vars(ff))
    Unewdata <- Unewdata[,setdiff(names(Unewdata), c(outcome, id.var))]
    attr(Unewdata,"outcome") <- outcome
    attr(Unewdata,"covariates") <- as.character(rhs.vars(ff))
    attr(Unewdata,"time.var") <- time.var
    attr(Unewdata,"confint") <- confint
    class(Unewdata) <- append("ufit",class(newdata))
    return(Unewdata)
}

## * autoplot.ufit
#' @title Display Unique Fitted Values
#' 
#' @description Display all the possible fitted values among the possible combinaisons of categorical covariables.
#' Continuous covariates are set to a single value.
#' 
#' @param object output of the \code{ufit} function.
#' @param x [character vector] time variable.
#' @param sep.strata [character] symbol used to separate the strata values in the caption
#' @param geom_confint [character] geometry used to display the confidence intervals.
#' Can be \code{"none"} for no display, \code{"errorbar"} to use intervals, or \code{"ribbon"} to use a shaded area.
#' @param alpha.ribbon [numeric 0-1] transparency parameter used to display the confidence intervals when using \code{geom_confint = "ribbon"}.
#' @param size.point [numeric >0] size argument used in \code{geom_point} to display point estimates.
#' @param size.line [numeric >0] size argument used in \code{geom_line} to connect point estimates belonging to the same strata.
#' @param size.ci [numeric >0] size argument used in \code{geom_errorbar} to display the confidence intervals.
#' @export
autoplot.ufit <- function(object, x = NULL,
                          sep.strata = ", ", geom_confint = "errorbar",
                          size.point = 3, size.line = 1.75, size.ci = 1, alpha.ribbon = 0.5){

    name.cov <- attr(object,"covariates")
    geom_confint <- match.arg(geom_confint, c("none","errorbar","ribbon"))

    ## ** find time variable
    if(is.null(x)){
        if(is.null(attr(object,"time.var"))){
            stop("Cannot guess what is the time variable. Consider specifying argument \'x\'. \n")
        }else{
            x <- attr(object,"time.var")
        }
        if(x %in% name.cov == FALSE){
            test.time <- sapply(name.cov, function(iCov){all(as.numeric(as.factor(object[[iCov]]))==as.numeric(as.factor(object[[x]])))})
            if(any(test.time)){
                x <- name.cov[test.time]
            }
        }
    }else{
        if(x %in% names(object) == FALSE){
            stop("Could not find the column corresponding to argument \'x\' in argument \'object\'. \n")
        }
        
        test.time <- sapply(name.cov, function(iCov){all(as.numeric(as.factor(object[[iCov]]))==as.numeric(as.factor(object[[x]])))})
        if(any(test.time)){
            name.cov <- name.cov[test.time==FALSE]
        }
    }
    
    ## ** convert x to factor
    level.x <- levels(object[[x]])
    object[[x]] <- as.numeric(object[[x]])
    
    ## ** find strata
    name.cov2 <- setdiff(name.cov,x)    
    if(length(name.cov2)==0){
        test.strata <- FALSE
    }else{
        test.strata <- TRUE
        if("strata" %in% names(object)){
            stop("Do not work when one of the variables in the dataset is called \"strata\". \n")
        }
        object$strata <- interaction(lapply(name.cov2, function(iCov){ ## iCov <- "age"
            return(paste0(iCov,"=",object[[iCov]]))
        }), sep = sep.strata)
    }
    
    ## ** display
    if(test.strata){
        gg <- ggplot2::ggplot(object, ggplot2::aes_string(x = x, y = "fit", group = "strata"))
        if(attr(object,"confint") && geom_confint != "none"){
            if(geom_confint == "errorbar"){
                gg <- gg + ggplot2::geom_errorbar(mapping = ggplot2::aes_string(ymin = "lower", ymax = "upper", color = "strata"), size = size.ci)
            }else if(geom_confint == "ribbon"){
                gg <- gg + ggplot2::geom_ribbon(mapping = ggplot2::aes_string(ymin = "lower", ymax = "upper", fill = "strata"), alpha = alpha.ribbon)
            }
        }
        gg <- gg + ggplot2::geom_point(ggplot2::aes_string(color = "strata"), size = size.point) + ggplot2::geom_line(ggplot2::aes_string(color = "strata"), size = size.line)
    }else{
        gg <- ggplot2::ggplot(object, ggplot2::aes_string(x = x, y = "fit"))
        if(attr(object,"confint") && geom_confint != "none"){
            if(geom_confint == "errorbar"){
                gg <- gg + ggplot2::geom_errorbar(mapping = ggplot2::aes_string(ymin = "lower", ymax = "upper"), size = size.ci)
            }else if(geom_confint == "ribbon"){
                gg <- gg + ggplot2::geom_ribbon(mapping = ggplot2::aes_string(ymin = "lower", ymax = "upper"), alpha = alpha.ribbon)
            }
        }
        gg <- gg + ggplot2::geom_point(size = size.point) + ggplot2::geom_line(size = size.line)
    }
    if(!is.null(level.x)){
        gg <- gg + ggplot2::scale_x_discrete(limits = level.x)
    }
    gg <- gg + ggplot2::ylab(attr(object,"outcome"))

    ## ** export
    return(gg)
}
