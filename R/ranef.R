### ranef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 26 2022 (11:18) 
## Version: 
## Last-Updated: mar 12 2024 (18:10) 
##           By: Brice Ozenne
##     Update #: 446
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ranef.lmm (documentation)
##' @title Estimate Random Effect From a Linear Mixed Model
##' @description Recover the random effects from the variance-covariance parameter of a linear mixed model.
##' @param object a \code{lmm} object.
##' @param effects [character] should the estimated random effects (\code{"mean"}) or the estimated variance/standard deviation of the random effects (\code{"variance"},\code{"std"}) be output?
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param ci [logical] should standard error and confidence intervals be evaluated using a delta method?
##' Will slow down the execution of the function.
##' @param format [character] should each type of random effect be output in a data.frame (\code{format="long"})
##' @param transform [logical] should confidence intervals for the variance estimates (resp. relative variance estimates) be evaluated using a log-transform (resp. atanh transformation)?
##' @param simplify [logical] when relevant will convert list with a single element to vectors and omit unessential output.
##' @param ... for internal use.
##'
##' @details Consider the following mixed model:
##' \deqn{Y = X\beta + \epsilon = X\beta + Z\eta + \xi}
##' where the variance of \eqn{\epsilon} is denoted \eqn{\Omega},
##' the variance of \eqn{\eta} is denoted \eqn{\Omega_{\eta}},
##' and the variance of \eqn{\xi} is \eqn{\sigma^2 I} with \eqn{I} is the identity matrix. \cr
##' The random effets are estimating according to:
##' \deqn{E[Y|\eta] = \Omega_{\eta} Z^{t} \Omega^{-1} (Y-X\beta)}
##' 
##' @keywords methods
##' 
##' @return A data.frame or a list depending on the argument \code{format}.
##' 
##' @examples
##' if(require(nlme)){
##' data(gastricbypassL, package = "LMMstar")
##' 
##' ## random intercept
##' e.RI <- lmm(weight ~ time + (1|id), data = gastricbypassL)
##' ranef(e.RI, effects = "mean")
##' ranef(e.RI, effects = "variance")
##'
##' }

## * ranef.lmm (code)
##' @export
ranef.lmm <- function(object, effects = "mean", ci = FALSE, transform = (effects %in% c("std","variance")),
                      p = NULL, format = "long", simplify = TRUE, ...){



    ## ** normalize user input
    mycall <- match.call()
    if(!inherits(object$design$vcov,"RE")){
        stop("Cannot estimate random effects linear mixed models defined by covariance structure (argument \'structure\'). \n",
             "Consider adding random effects in the argument \'formula\' instead. \n")
    }
    effects <- match.arg(effects, c("mean","std","variance"))
    format <- match.arg(format, c("wide","long"))

    param.name <- object$design$param$name
    if(!is.null(p)){
        if(any(duplicated(names(p)))){
            stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
        }
        if(any(param.name %in% names(p) == FALSE)){
            stop("Incorrect argument \'p\': missing parameter(s) \"",paste(param.name[param.name %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
        }
        p <- p[param.name]
    }else{
        p <- object$param
    }

    if(format == "wide" && "ci" %in% names(mycall) && ci == TRUE){
        message("Argument \'format\' ignored when argument \'ci\' is TRUE. \n")
        format <- "long"
    }
        
    if(transform && effects == "mean" && ci){
        stop("Argument \'transform\' should be FALSE when evaluating confidence intervals for the random effects. \n")
    }

    ## ** ci
    if(ci){
        ## by default do not compute degrees of freedom (not reliable)
        dots <- list(...)
        df <- dots$df
        if(is.null(df)){
            df <- FALSE
        }else{
            dots$df <- NULL
        }

        e.ranef <- nlme::ranef(object, effects = effects, ci = FALSE, p = p, format = format)
        e.delta <- lava::estimate(object, f = function(newp){
            iE <- nlme::ranef(object, effects = effects, ci = FALSE, p = newp, format = format)
            return(iE$estimate)
        }, df = df)
            
        if(transform){ ## recompute only CIs (backtransforming the se is not exact)
            eTrans.delta <- lava::estimate(object, f = function(newp){
                iE <- nlme::ranef(object, effects = effects, ci = FALSE, p = newp, format = format)
                iE[iE$type==effects,"estimate"] <- log(iE[iE$type==effects,"estimate"])
                iE[iE$type=="relative","estimate"] <- atanh(iE[iE$type=="relative","estimate"])
                return(iE$estimate)
            }, df = df)
            ## absolute
            e.delta$lower[e.ranef$type==effects] <- exp(eTrans.delta$lower[e.ranef$type==effects])
            e.delta$upper[e.ranef$type==effects] <- exp(eTrans.delta$upper[e.ranef$type==effects])
            ## relative
            e.delta$lower[e.ranef$type=="relative"] <- tanh(eTrans.delta$lower[e.ranef$type=="relative"])
            e.delta$upper[e.ranef$type=="relative"] <- tanh(eTrans.delta$upper[e.ranef$type=="relative"])
        }

        out <- cbind(e.ranef, e.delta[,c("se","df","lower","upper")])

        return(out)
    }else{
        dots <- list(...)
        if(length(dots)>0){
            stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
        }
    }

    ## param
    param.type <- stats::setNames(object$design$param$type,param.name)
    param.rho <- param.name[param.type=="rho"]
    param.strata <- unlist(object$design$param[param.type=="rho","index.strata"])

    ## cluster 
    var.cluster <- object$cluster$var
    n.cluster <- object$design$cluster$n
    index.cluster <- object$design$index.cluster
    Vindex.cluster <- attr(index.cluster, "vectorwise")
    
    ## strata
    var.strata <- object$strata$var
    index.clusterStrata <- object$design$index.clusterStrata
    U.strata <- object$strata$levels
    n.strata <- length(U.strata)

    ## design
    X.cor <- object$design$vcov$cor$X
    Xpattern.cor <- object$design$vcov$cor$Xpattern
    pattern.cluster <- object$design$vcov$pattern
    Upattern <- object$design$vcov$Upattern
    infoRanef <- object$design$vcov$ranef
    name.hierarchy <- unlist(infoRanef$hierarchy, use.names = FALSE)
    index.hierarchy <- unlist(lapply(1:length(infoRanef$hierarchy), function(iH){rep(iH,length(infoRanef$hierarchy[[iH]]))}))

    ## missing values
    index.na <- object$index.na

    ## all vars
    var.all <- unique(lava::manifest(object))
    
    ## ** converting correlation parameters into random effect variance
    cumtau <- coef(object, p = p, effects = "correlation", transform.rho = "cov", transform.names = FALSE)

    n.hierarchy <- length(infoRanef$param)
    index.hierarchy <- unlist(lapply(1:n.hierarchy, function(iH){rep(iH, length(infoRanef$param[[iH]]))}))
    name.RE <- unname(unlist(lapply(infoRanef$param, rownames)))
    n.RE <- length(name.RE)
    varRE <- matrix(NA, nrow = n.RE, ncol = n.strata,
                    dimnames =  list(name.RE, U.strata))
    for(iH in 1:n.hierarchy){ ## iH <- 1
        iHierarchy <- infoRanef$param[[iH]]
        varRE[rownames(iHierarchy),] <- apply(iHierarchy, MARGIN = 2, FUN = function(iName){
            cumtau[iName] - c(0,utils::head(cumtau[iName],-1))
        })
    }
    
    if(any(varRE<=0)){
        stop("Variance for the random effects is found to be negative - cannot estimate the random effects. \n")
    }
    if(effects %in% c("std","variance")){
        sigma2 <- coef(object, effects = "variance", transform.sigma = "square")

        if(format=="long"){
            out <- do.call(rbind,lapply(1:n.strata, function(iStrata){
                rbind(data.frame(variable = rownames(varRE),
                                 strata = U.strata[iStrata],
                                 type = effects,
                                 estimate = varRE[,iStrata]),
                      data.frame(variable = rownames(varRE),
                                 strata = U.strata[iStrata],
                                 type = "relative",
                                 estimate = varRE[,iStrata]/sigma2[iStrata])
                      )
            }))
            if(effects == "std"){
                out$estimate <- sqrt(out$estimate)
            }
        }else if(format=="wide"){
            out <- do.call(rbind,lapply(1:n.strata, function(iStrata){ ## iStrata <- 1
                iOut <- data.frame(variable = rownames(varRE),
                                   strata = U.strata[iStrata],
                                   variance = varRE[,iStrata],
                                   relative = varRE[,iStrata]/sigma2[iStrata])
                if(simplify == FALSE){
                    iOut <- rbind(data.frame(variable = "total", strata = U.strata[iStrata], variance = sigma2[iStrata], relative = 1),
                                  iOut,
                                  data.frame(variable = "residual", strata = U.strata[iStrata], variance = sigma2[iStrata]-sum(iOut$variance), relative = 1-sum(iOut$relative))
                                  )
                }
                return(iOut)
            }))
            if(effects == "std"){
                out$variance <- sqrt(out$variance)
                out$relative <- sqrt(out$relative)
            }
        }
        rownames(out) <- NULL
        if(simplify && n.strata==1){
            out$strata <- NULL
        }
        return(out)
    }

    ## ** extract normalized residuals
    ## head(stats::residuals(object, p = p, keep.data = TRUE, type = "response", format = "long"))
    df.epsilon <- stats::residuals(object, p = p, keep.data = TRUE, type = "normalized2", format = "long")
    if(object$strata$n==1){
        df.epsilon$XXstrata.indexXX <- U.strata
    }
    
    ## ** estimate random effects
    grid.ranef <- unlist(lapply(1:n.hierarchy, function(iH){ ## iH <- 1
        stats::setNames(lapply(1:length(infoRanef$hierarchy[[iH]]), function(iP){
            infoRanef$hierarchy[[iH]][1:iP]
        }), infoRanef$hierarchy[[iH]])
    }), recursive = FALSE)

    ls.out <- lapply(1:length(grid.ranef), function(iR){ ## iR <- 1
        do.call(rbind,by(df.epsilon, df.epsilon[grid.ranef[[iR]]], function(iDF){ ## iDF <- df.epsilon[df.epsilon$Subject == "F10",]
            iNames <- grid.ranef[[iR]]
            iTau <- varRE[names(grid.ranef)[iR],which(iDF[[var.strata]][1]==U.strata)]
            iOut <- data.frame(variable = NA,
                               strata = iDF[[var.strata]][1],
                               level = NA,
                               estimate = unname(sum(iDF$r.normalized)*iTau))
            iOut$variable <- list(iNames)
            iOut$level <- list(unlist(iDF[1,iNames,drop=FALSE]))
            return(iOut)
        }))
    })

    
    ## ** export
    out <- do.call(rbind,ls.out)
    rownames(out) <- NULL
    if(format == "wide"){
        out <- stats::setNames(lapply(1:n.hierarchy, function(iH){ ## iH <- 1
            iOut <- out[sapply(out$variable, utils::tail,1) %in% infoRanef$hierarchy[[iH]],,drop=FALSE]
            iOut$col <- sapply(iOut$level, function(iLevel){paste0(iLevel[-1], collapse = ":")})
            iOut$level <- sapply(iOut$level, utils::head, 1)
            iOutW <- stats::reshape(iOut, direction = "wide", 
                                    idvar = "level", timevar = "col", times = unique(iOut$col), drop = c("strata","variable"))
            colnames(iOutW)[2] <- "estimate"
            return(iOutW)
        }), sapply(infoRanef$hierarchy,"[",1))
        if(simplify && n.hierarchy == 1){
            names(out[[1]])[1] <- names(out)
            out <- out[[1]]            
        }
        
    }else if(format == "long"){
        if(simplify && n.strata == 1){
            out$strata <- NULL
        }
        if(simplify && all(lengths(infoRanef$hierarchy)==1)){
            out$variable <- unlist(out$variable)
            out$level <- unlist(out$level)
        }
    }
    return(out)
}

## * ranef.mlmm (code)
##' @export
ranef.mlmm <- function(object, ...){

    return(lapply(object$model, ranef, ...))
}

##----------------------------------------------------------------------
### ranef.R ends here
