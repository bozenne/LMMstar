### profile.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 16 2022 (15:19) 
## Version: 
## Last-Updated: jun 16 2022 (16:49) 
##           By: Brice Ozenne
##     Update #: 59
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * profile.lmm (documentation)
##' @title Display Contour of the log-Likelihood
##' @description Display the (restricted) log-likelihood when varying one parameter at a time.
##' @name profile
##'
##' @param fitted a \code{lmm} object.
##' @param effects [character vector] name of the parameters for which the (restricted) log-likelihood should be evaluated.
##' Alternatively can be the type of parameters, e.g. \code{"mean"}, \code{"variance"}, \code{"correlation"}, or \code{"all"}.
##' @param maxpts [integer, >0] the number of points use to discretize the likelihood.
##' @param conf.level [numeric, 0-1] the confidence level of the confidence intervals used to decide about the range of values for each parameter.
##' @param trace [logical] Show the progress of the execution of the function.
##' @param plot [logical] Should a graphical representation of the results be provided?
##' @param size [numeric vector of length 3] Size of the point for the MLE, the width of the line representing the likelihood, and the width of the corresponding quadratic approximation.
##' @param linetype [integer, >0] type of line used to represent the quadratic approximation of the likelihood.
##' @param shape [integer, >0] type of point used to represent the MLE.
##' @param scales,nrow,ncol argument passed to \code{ggplot2::facet_wrap}.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' 
##' @details WARNING: this differs from profile likelihood (which contrain one parameter and update the other accordingly).
##' In particular, it should not be used to obtain confidence intervals.
##'
##' @return A data.frame object containing for various parameter values, the log-likelihood. When argument \code{plot} is TRUE, also contain an attribute \code{"plot"} containing a \code{ggplot} object.
##' 
##' @examples
##' data(gastricbypassW)
##' e.lmm <- lmm(weight2 ~ weight1 + glucagonAUC1, data = gastricbypassW)
##' df.profile <- profile(e.lmm, effects = "all", maxpts = 10)
##' df.profile <- profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "log")

## * profile.lmm (code)
##' @export
profile.lmm <- function(fitted, effects = NULL, 
                        maxpts = 100, conf.level = 0.95, trace = FALSE,
                        plot = TRUE, size = c(3,2,1), linetype = 2, shape = 19, scales = "free", nrow = NULL, ncol = NULL,
                        transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE, ...){

    ## ** normalize user input
    call <- match.call()
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    p <- coef(fitted, effects = "all")
    name.p <- names(p)
    type.p <- setNames(fitted$design$param$type, names(p))

    p.trans <- confint(fitted, effects = "all", level = conf.level,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                       transform.names = transform.names)
    name.p.trans <- rownames(p.trans)

    if(is.null(effects)){
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation",names(p)), several.ok = TRUE)
    if(any(effects %in% names(p) == FALSE)){
        effects <- names(coef(fitted, effects = effects))
    }
    n.effects <- length(effects)

    ## ** profile likelihood
    if(trace>1){cat("Profile likelihood (",maxpts," points):\n")}
    ls.profile <- lapply(1:n.effects, function(iParam){ ## iParam <- 1

        iIndex <- which(name.p == effects[iParam])
        iName.trans <- name.p.trans[iIndex]
        iType <- unname(type.p[iIndex])
        iValue <- unname(p[iIndex])
        iValue.trans <- p.trans[iIndex,"estimate"]
        iLower.trans <- p.trans[iIndex,"lower"]
        iUpper.trans <- p.trans[iIndex,"upper"]

        if(trace>0){
            if(trace<=1){cat("*")}
            if(trace>1){cat(" - ",iName.trans," (between ",iLower.trans," and ",iUpper.trans,")\n",sep="")}
        }
        
        seqValue.trans <- c(iValue.trans,seq(iLower.trans, iUpper.trans, length.out = maxpts))
        seqValue <- .reparametrize(p = seqValue.trans, type = rep(iType,length(seqValue.trans)),
                                   transform.sigma = transform.sigma,
                                   transform.k = transform.k,
                                   transform.rho = transform.rho,
                                   transform.names = FALSE,
                                   inverse = TRUE, Jacobian = FALSE, dJacobian = FALSE)$p

        data.frame(param = iName.trans,
                   type = iType,
                   value = seqValue,
                   value.trans = seqValue.trans,
                   optimum = c(TRUE,rep(FALSE,maxpts)),
                   logLik = sapply(seqValue, function(iiValue){iP <- p; iP[iIndex] <- iiValue; logLik(fitted, p = iP)})
                   )
    })
    if(trace>0){cat("\n")}

    ## ** collect
    df.profile <- do.call(rbind,ls.profile)
df.profile$param <- factor(df.profile$param, levels = unique(df.profile$param))
    
    ## ** display
    if(plot){
        gg <- ggplot(df.profile[df.profile$optimum==FALSE,,drop=FALSE], aes(x=value.trans,y=logLik))
        gg <- gg + ylab(expression(paste("log-likelihood when varying a ",bold('single')," parameter")))
        gg <- gg + facet_wrap(~param, scales= scales, nrow = nrow, ncol = ncol)
        if(size[2]>0){
            gg <- gg + geom_line(size = size[2])
        }
        if(size[3]>0){
            gg <- gg + stat_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2), size = size[3], linetype = linetype, aes(color = "quadratic approximation"))
        }
        if(size[1]>0){
            gg <- gg  + geom_point(data = df.profile[df.profile$optimum==TRUE,,drop=FALSE], aes(color = "MLE"), size = size[1], shape = shape)
        }
        gg <- gg  + xlab("") + labs(color = "") + theme(legend.position = "bottom")
        if(plot>=1){
            print(gg)
        }
        attr(df.profile,"plot") <- gg
    }
    ## ** export
    return(df.profile)
}

##----------------------------------------------------------------------
### profile.R ends here
