### profile.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 16 2022 (15:19) 
## Version: 
## Last-Updated: aug  8 2024 (13:32) 
##           By: Brice Ozenne
##     Update #: 340
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * profile.lmm (documentation)
##' @title Log-Likelihood Contour For a Linear Mixed Model
##' @description Evaluate the (restricted) log-likelihood around Maximum Likelihood Estimate (MLE) of a linear mixed model.
##' It will vary one parameter and either keep the other parameters constant or set the other parameters at their MLE (profile likelihood).
##' Since a locally quadratic log-likelihood with an Hessian equivariant in law implies normally distributed estimates (Geyer 2013)
##' it can help trusting confidence intervals and p-values in small samples with a non-normally distributed outcome.
##' 
##' @param fitted a \code{lmm} object.
##' @param effects [character vector] name of the parameters who will be constrained.
##' Alternatively can be the type of parameters, e.g. \code{"mean"}, \code{"variance"}, \code{"correlation"}, or \code{"all"}.
##' @param profile.likelihood [logical] should the other parameters be set to the value maximizing the constrained likelihood or stay at their MLE value? 
##' @param maxpts [integer, >0] number of points use to discretize the likelihood, \code{maxpts} points smaller than the MLE and \code{maxpts} points higher than the MLE.
##' @param level [numeric, 0-1] the confidence level of the confidence intervals used to decide about the range of values for each parameter.
##' @param trace [logical] Show the progress of the execution of the function.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' 
##' @details Each parameter defined by the argument \code{effets} is treated separately:\itemize{
##' \item the confidence interval of a parameter is discretized with \code{maxpt} points,
##' \item this parameter is set to each discretization value.
##' \item the other parameters are either set to the (unconstrained) MLE (\code{profile.likelihood=FALSE})
##' or to constrained MLE (\code{profile.likelihood=TRUE}). The latter case is much more computer intensive as it implies re-running the estimation procedure.
##' \item the (restricted) log-likelihood is evaluated.
##' }
##'
##' @return A data.frame object containing the log-likelihood for various parameter values.
##'
##' @references Geyer, C. J. (2013). Asymptotics of maximum likelihood without the lln or clt or sample size going to infinity. In Advances in Modern Statistical Theory and Applications: A Festschrift in honor of Morris L. Eaton, pages 1–24. Institute of Mathematical Statistics.
##' @keywords htest
##' 
##' @examples
##' data(gastricbypassW, package = "LMMstar")
##' e.lmm <- lmm(weight2 ~ weight1 + glucagonAUC1,
##'              data = gastricbypassW, control = list(optimizer = "FS"))
##'
##' ## profile logLiklihood
##' \dontrun{
##' e.pro <- profile(e.lmm, effects = "all", maxpts = 10, profile.likelihood = TRUE)
##' head(e.pro)
##' plot(e.pro)
##' }
##' 
##' ## along a single parameter axis
##' e.sliceNone <- profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "none")
##' plot(e.sliceNone)
##' e.sliceLog <- profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "log")
##' plot(e.sliceLog)
##' 

## * profile.lmm (code)
##' @export
profile.lmm <- function(fitted, effects = NULL, profile.likelihood = FALSE,
                        maxpts = NULL, level = 0.95, trace = FALSE,                        
                        transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    table.param <- stats::model.tables(fitted, effects = "param")
    name.p <- table.param$name
    type.p <- stats::setNames(table.param$type, name.p)

    ## ** normalize user input
    ## *** dots
    dots <- list(...)
    if("options" %in% names(dots) && !is.null(dots$options)){
        options <- dots$options
    }else{
        options <- LMMstar.options()
    }
    dots$options <- NULL
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## *** transformations
    init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = fitted$reparametrize$transform.sigma, x.transform.k = fitted$reparametrize$transform.k, x.transform.rho = fitted$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    ## *** level
    p.trans <- stats::confint(fitted, effects = "all", level = level,
                              transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                              transform.names = transform.names, backtransform = FALSE, options = options)
    name.p.trans <- rownames(p.trans)
    rownames(p.trans) <- name.p

    ## *** effects
    if(is.null(effects)){
        if((is.null(transform.sigma) || identical(transform.sigma,"none")) && (is.null(transform.k) || identical(transform.k,"none")) && (is.null(transform.rho) || identical(transform.rho,"none"))){
            effects <- options$effects
        }else{
            effects <- c("mean","variance","correlation")
        }
    }else{
        if(!is.character(effects) || !is.vector(effects)){
            stop("Argument \'effects\' must be a character vector. \n")
        }
        valid.effects <- c("mean","fixed","variance","correlation","all",name.p)
        if(any(effects %in% valid.effects == FALSE)){
            stop("Incorrect value for argument \'effect\': \"",paste(setdiff(effects,valid.effects), collapse ="\", \""),"\". \n",
                 "Valid values: \"",paste(valid.effects, collapse ="\", \""),"\". \n")
        }
        if(all("all" %in% effects)){
            if(length(effects)>1){
                stop("Argument \'effects\' must have length 1 when containing the element \"all\". \n")
            }else{
                effects <- c("mean","variance","correlation")
            }
        }else{
            effects[effects == "fixed"] <- "mean"
        }
        if(any(effects %in% name.p == FALSE)){
            effects <- stats::model.tables(fitted, effects = c("param",effects))$name
        }
    
    }
    n.effects <- length(effects)

    ## *** profile.likelihood
    if(fitted$args$control$optimizer!="FS" && profile.likelihood>0){
        stop("Argument \'profile.likelihood\' can only be TRUE when \"FS\" optimizer is used. \n",
             "Consider adding the argument control = list(optimizer = \"FS\") when fitting the mixed model with lmm. \n")
    }

    ## *** maxpt
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

    ## ** profile likelihood
    if(trace>1){cat("Profile likelihood (",round(2*maxpts)," points):\n",sep="")}
    ls.profile <- lapply(1:n.effects, function(iParam){ ## iParam <- 4

        iIndex <- which(name.p == effects[iParam])
        iName <- name.p[iIndex]
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

            for(iPts in 1:iMaxpts){ ## iPts <- 13
                if(iIndex.center-iPts>0){
                    iResInf <- try(.constrain.lmm(fitted, effects = stats::setNames(seqValue[iIndex.center-iPts], effects[iParam]), init = iInitInf, trace = FALSE))
                    if(!inherits(iResInf,"try-error")){
                        iOut[iIndex.center-iPts, c("logLik","cv")] <- c(logLik = iResInf$logLik, cv = iResInf$opt$cv)
                        iInitInf <- iResInf$param
                        if(profile.likelihood>1){
                            iOut[iIndex.center-iPts, name.p] <- iResInf$param
                        }
                    }
                }
                
                if(iIndex.center+iPts<=length(seqValue)){
                    iResSup <- try(.constrain.lmm(fitted, effects = stats::setNames(seqValue[iIndex.center+iPts], effects[iParam]), init = iInitSup, trace = FALSE))
                    if(!inherits(iResSup,"try-error")){
                        iOut[iIndex.center+iPts, c("logLik","cv")] <- c(logLik = iResSup$logLik, cv = iResSup$opt$cv)                    
                        iInitSup <- iResSup$param
                        if(profile.likelihood>1){
                            iOut[iIndex.center+iPts, name.p] <- iResSup$param
                        }
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
    ## unique(df.profile$param)

    ## ** export
    attr(df.profile, "args") <- list(profile.likelihood = profile.likelihood,
                                     logLik = fitted$logLik,
                                     maxpts = maxpts,
                                     name.p = name.p,
                                     effects = effects,
                                     conf.level = level,
                                     ci = p.trans,
                                     transform.sigma = transform.sigma,
                                     transform.k = transform.k,
                                     transform.rho = transform.rho,
                                     transform.names = transform.names)
    rownames(attr(df.profile, "args")$ci) <- name.p.trans
    class(df.profile) <- append("profile_lmm", class(df.profile))
    return(df.profile)
}

##----------------------------------------------------------------------
### profile.R ends here
