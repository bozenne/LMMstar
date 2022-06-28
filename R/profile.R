### profile.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 16 2022 (15:19) 
## Version: 
## Last-Updated: jun 28 2022 (11:44) 
##           By: Brice Ozenne
##     Update #: 261
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
##' @description Display the (restricted) log-likelihood around Maximum Likelihood Estimate (MLE) under specific constrains.
##' @name profile
##'
##' @param fitted a \code{lmm} object.
##' @param effects [character vector] name of the parameters who will be constrained.
##' Alternatively can be the type of parameters, e.g. \code{"mean"}, \code{"variance"}, \code{"correlation"}, or \code{"all"}.
##' @param profile.likelihood [logical] should profile likelihood be performed? Otherwise varying one parameter at a time around the MLE while keeping the other constant).
##' @param maxpts [integer, >0] number of points use to discretize the likelihood, \code{maxpts} points smaller than the MLE and \code{maxpts} points higher than the MLE.
##' @param conf.level [numeric, 0-1] the confidence level of the confidence intervals used to decide about the range of values for each parameter.
##' @param trace [logical] Show the progress of the execution of the function.
##' @param plot [logical] Should a graphical representation of the results be provided?
##' @param ci [logical] Should a 95\% confidence intervals obtained from the Wald test (vertical lines) and Likelihood ratio test (horizontal line) be displayed?
##' @param size [numeric vector of length 4] Size of the point for the MLE,
##' width of the line representing the likelihood,
##' width of the corresponding quadratic approximation,
##' and width of the line representing the confidence intervals.
##' @param linetype [integer vector of length 2] type of line used to represent the quadratic approximation of the likelihood
##' and the confidence intervals.
##' @param shape [integer, >0] type of point used to represent the MLE.
##' @param scales,nrow,ncol argument passed to \code{ggplot2::facet_wrap}.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' 
##' @details Each parameter defined by the argument \code{effets} is treated separately:\itemize{
##' \item the confidence interval of a parameter is discretized with \code{maxpt} points,
##' \item this parameter is set to a discretization value.
##' \item the other parameters are either set to the (unconstrained) MLE (\code{profile.likelihood=FALSE})
##' or to constrained MLE  (\code{profile.likelihood=TRUE}). The latter case is much more computer intensive as it implies re-running the estimation procedure.
##' \item the (restricted) log-likelihood is evaluated.
##' }
##'
##' @return A data.frame object containing the log-likelihood for various parameter values.
##' When argument \code{plot} is TRUE, also contain an attribute \code{"plot"} containing a \code{ggplot} object.
##' 
##' @examples
##' data(gastricbypassW, package = "LMMstar")
##' e.lmm <- lmm(weight2 ~ weight1 + glucagonAUC1,
##'              data = gastricbypassW, control = list(optimizer = "FS"))
##'
##' ## profile logLiklihood
##' \dontrun{
##' profile(e.lmm, effects = "all", maxpts = 10, profile.likelihood = TRUE)
##' }
##' 
##' ## along a single parameter axis
##' profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "none")
##' profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "log")
##' 

## * profile.lmm (code)
##' @export
profile.lmm <- function(fitted, effects = NULL, profile.likelihood = FALSE,
                        maxpts = NULL, conf.level = 0.95, trace = FALSE,
                        plot = TRUE, ci = FALSE, size = c(3,2,1,1), linetype = c("dashed","dashed","dashed"), shape = 19, scales = "free", nrow = NULL, ncol = NULL,
                        transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    ## ** normalize user input
    call <- match.call()
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    
    p <- confint(fitted, effects = "all", level = conf.level,
                 transform.sigma = "none", transform.k = "none", transform.rho = "none")
    name.p <- rownames(p)
    type.p <- stats::setNames(fitted$design$param$type, name.p)

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = fitted$reparametrize$transform.sigma, x.transform.k = fitted$reparametrize$transform.k, x.transform.rho = fitted$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    p.trans <- confint(fitted, effects = "all", level = conf.level,
                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                       transform.names = transform.names)
    name.p.trans <- rownames(p.trans)

    if(is.null(effects)){
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation",name.p), several.ok = TRUE)
    if(any(effects %in% name.p == FALSE)){
        effects <- names(coef(fitted, effects = effects))
    }
    n.effects <- length(effects)

    if(fitted$opt$name!="FS" && profile.likelihood>0){
        stop("Argument \'profile.likelihood\' can only be TRUE when \"FS\" optimizer is used. \n",
             "Consider adding the argument control = list(optimizer = \"FS\") when fitting the mixed model with lmm. \n")
    }

    if(is.null(maxpts)){
        if(profile.likelihood>0){
            maxpts <- 15
        }else{
            maxpts <- 50
        }
        grid <- NULL
    }else if(length(maxpts)==1){
        grid <- NULL
    }else if(length(maxpts)>=1){
        grid <- maxpts
        maxpts <- length(maxpts)/2
    }

    if(ci && profile.likelihood==FALSE){
        stop("Can only display the confidence intervals when performing profile likelihood. \n")
    }

    ## ** profile likelihood
    if(trace>1){cat("Profile likelihood (",round(2*maxpts)," points):\n",sep="")}
    ls.profile <- lapply(1:n.effects, function(iParam){ ## iParam <- 4

        iIndex <- which(name.p == effects[iParam])
        iName.trans <- name.p.trans[iIndex]
        iType <- unname(type.p[iIndex])
        iValue <- unname(p[iIndex,"estimate"])
        iValue.trans <- p.trans[iIndex,"estimate"]
        iLower.trans <- p.trans[iIndex,"lower"]
        iUpper.trans <- p.trans[iIndex,"upper"]

        if(is.null(grid)){
            seqValue.trans <- c(seq(iLower.trans, iValue.trans, length.out = maxpts+1), seq(iValue.trans, iUpper.trans, length.out = maxpts+1)[-1])
            seqOptimum <- c(rep(FALSE,maxpts),TRUE,rep(FALSE,maxpts))
        }else{
            seqValue.trans <- sort(unique(c(grid, iValue.trans)))
            seqOptimum <- seqValue.trans == iValue.trans
        }
        if(trace>0){
            if(trace<=1){cat("*")}
            if(trace>1){cat(" - ",iName.trans," (between ",min(seqValue.trans)," and ",max(seqValue.trans),")\n",sep="")}
        }

        iMaxpts <- max(which(seqOptimum==TRUE)-1,length(seqOptimum)-which(seqOptimum==TRUE))
        
        iIndex.center <- which(seqOptimum==TRUE)
        seqValue <- .reparametrize(p = seqValue.trans, type = rep(iType,length(seqValue.trans)),
                                   transform.sigma = transform.sigma,
                                   transform.k = transform.k,
                                   transform.rho = transform.rho,
                                   transform.names = FALSE,
                                   inverse = TRUE, Jacobian = FALSE, dJacobian = FALSE)$p

        iOut <- data.frame(param = iName.trans,
                           type = iType,
                           value = seqValue,
                           value.trans = seqValue.trans,
                           optimum = seqOptimum,
                           logLik = NA,
                           cv = NA)
        iOut[iIndex.center,"logLik"] <- fitted$logLik
        iOut[iIndex.center,"cv"] <- TRUE

        keep.estimate <- NULL

        if(profile.likelihood>0){
            iInitInf <- stats::setNames(p[,"estimate"], name.p)
            iInitSup <- stats::setNames(p[,"estimate"], name.p)
            if(profile.likelihood>1){
                iOut[name.p] <- NA
                iOut[iIndex.center,name.p] <- p[,"estimate"]
            }

            for(iPts in 1:iMaxpts){ ## iPts <- 15 

                if(iIndex.center-iPts>0){
                    iResInf <- .constrain.lmm(fitted, effects = stats::setNames(seqValue[iIndex.center-iPts], effects[iParam]), init = iInitInf, trace = FALSE)
                    iOut[iIndex.center-iPts, c("logLik","cv")] <- c(logLik = iResInf$logLik, cv = iResInf$opt$cv)
                    iInitInf <- iResInf$param
                    if(profile.likelihood>1){
                        iOut[iIndex.center-iPts, name.p] <- iResInf$param
                    }
                }
                
                if(iIndex.center+iPts<=length(seqValue)){
                    iResSup <- .constrain.lmm(fitted, effects = stats::setNames(seqValue[iIndex.center+iPts], effects[iParam]), init = iInitSup, trace = FALSE)
                    iOut[iIndex.center+iPts, c("logLik","cv")] <- c(logLik = iResSup$logLik, cv = iResSup$opt$cv)
                    iInitSup <- iResSup$param
                    if(profile.likelihood>1){
                        iOut[iIndex.center+iPts, name.p] <- iResSup$param
                    }
                }
            }


        }else{
            iOut$logLik[-iIndex.center] <- sapply(seqValue[-iIndex.center], function(iiValue){
                iP <- stats::setNames(p[,"estimate"],name.p)
                iP[iIndex] <- iiValue
                return(logLik(fitted, p = iP))
            })
            iOut$cv <- TRUE
        }

        iOut$value.transC <- iOut$value.trans - iOut[iOut$optimum==TRUE,"value.trans"]
        iOut$likelihood.ratio <- exp(iOut$logLik - fitted$logLik)
        return(iOut)
    })
    if(trace>0){cat("\n")}

    ## ** collect
    df.profile <- do.call(rbind,ls.profile)
    df.profile$param <- factor(df.profile$param, levels = unique(df.profile$param))
    
    ## ** display
    if(plot){
        if(plot>=2){
            name.y <- "likelihood.ratio"
            reference <- 1
            legend.y <- "Likelihood relative to the maximum likelihood"
            fff <- likelihood.ratio ~ 0+I(value.trans-mean(value.trans)) + I((value.trans-mean(value.trans))^2)
        }else{
            name.y <- "logLik"
            reference <- fitted$logLik
            legend.y <- "Log-likelihood"
            fff <- logLik ~ 0+I(value.trans-mean(value.trans)) + I((value.trans-mean(value.trans))^2)
        }
        gg <- ggplot2::ggplot(df.profile, ggplot2::aes_string(x="value.trans",y=name.y))
        gg <- gg + ggplot2::ylab(legend.y)
        if(profile.likelihood>0){
            if(ci){
                gg <- gg + ggplot2::ggtitle(paste("Profile maximum likelihood estimation for parameter (95% CI):"))
            }else{
                gg <- gg + ggplot2::ggtitle(paste("Profile maximum likelihood estimation for parameter:"))
            }
        }else{
            gg <- gg + ggplot2::ggtitle(expression(paste("Varying a ",bold('single')," parameter:")))
        }
        gg <- gg + ggplot2::facet_wrap(~param, scales= scales, nrow = nrow, ncol = ncol)
        if(size[2]>0){
            gg <- gg + ggplot2::geom_line(size = size[2])
        }
        if(size[3]>0 && (plot<=1)){
            df.fit <- do.call(rbind,by(df.profile, df.profile$param, function(iDF){ ## iDF <- df.profile[df.profile$param=="sigma",]
                iDF$myset <- reference
                ## FOR CRAN test
                myset <- NULL
                iLM <- stats::lm(fff, data = iDF, offset = myset)
                iDF[[name.y]] <- stats::predict(iLM, newdata = iDF)
                return(iDF)
            }))
            df.fit$param <- factor(df.fit$param, levels = levels(df.profile$param))
            gg <- gg + ggplot2::geom_line(data = df.fit, size = size[3], linetype = linetype[1], ggplot2::aes(color = "quadratic approximation"))
        }
        if(size[1]>0){
            gg <- gg  + ggplot2::geom_point(data = df.profile[df.profile$optimum==TRUE,,drop=FALSE], ggplot2::aes(color = "MLE"), size = size[1], shape = shape)
        }
        if(ci){
            if(plot<=1){
                gg <- gg  + ggplot2::geom_abline(slope = 0, intercept = fitted$logLik - stats::qchisq(0.95, df = 1)/2, size = size[4], linetype = linetype[2])
            }else{
                gg <- gg  + ggplot2::geom_abline(slope = 0, intercept = exp(- stats::qchisq(0.95, df = 1)/2), size = size[4], linetype = linetype[2])
            }

            ## recompute ci at level 0.95
            p.trans2 <- confint(fitted, effects = "all", level = 0.95,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                               transform.names = transform.names)

            df.ci <- cbind(param = rownames(p.trans2), p.trans2)[which(name.p %in% effects),,drop=FALSE]
            df.ci$param <- factor(df.ci$param, levels = levels(df.profile$param))
            gg <- gg + ggplot2::geom_vline(data = df.ci, mapping = ggplot2::aes_string(xintercept = "lower"), size = size[4], linetype = linetype[2])
            gg <- gg + ggplot2::geom_vline(data = df.ci, mapping = ggplot2::aes_string(xintercept = "upper"), size = size[4], linetype = linetype[2])
        }
        gg <- gg  + ggplot2::xlab("") + ggplot2::labs(color = "") + ggplot2::theme(legend.position = "bottom")
        if(plot>=1){
            print(gg)
        }
        attr(df.profile,"plot") <- gg
    }
    ## ** export
    if(plot>=1){
        return(invisible(df.profile))
    }else{
        return(df.profile)
    }
}

##----------------------------------------------------------------------
### profile.R ends here
