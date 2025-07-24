### profile.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 16 2022 (15:19) 
## Version: 
## Last-Updated: jul 24 2025 (16:13) 
##           By: Brice Ozenne
##     Update #: 541
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
##' The values of a given parameter are varied over a pre-defined grid and the corresponding (contrained) likelihood w.r.t. each value is evaluated.
##' The other parameters are either kept constant or set to maximize the contrained likelihood (profile likelihood).
##' In the latter case, confidence intervals consistent with a likelihood ratio test (LRT) can be output.
##' @param fitted a \code{lmm} object.
##' @param effects [character vector] name of the parameters to be constrained.
##' Alternatively can be the type of parameters, e.g. \code{"mean"}, \code{"variance"}, \code{"correlation"}, or \code{"all"}.
##' @param profile.likelihood [FALSE,TRUE,"ci"] Should the unconstrained parameter(s) be kept at their (MLE) value (\code{FALSE}),
##' or set to the value maximizing the constrained likelihood (\code{TRUE}),
##' or should confidence intervals be computed for the parameters defined in argument \code{effects}. 
##' @param maxpts [integer, >0] number of points use to discretize the likelihood, \code{maxpts} points smaller than the MLE and \code{maxpts} points higher than the MLE.
##' @param level [numeric, 0-1] the confidence level of the confidence intervals.
##' Used to decide about the range of values for each parameter when argument \code{profile.likelihood} is \code{TRUE} or \code{FALSE}.
##' @param df [logical] Should a Student's t-distribution be used to model the distribution of the coefficients when evaluating the confidence intervals. Otherwise a normal distribution is used.
##' Ignored when \code{profile.likelihood} is \code{"ci"}.
##' @param trace [logical] Show the progress of the execution of the function.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' 
##' @details Each parameter defined by the argument \code{effets} is treated separately either to evaluate the constrained likelihood or compute a confidence interval (argument \code{profile.likelihood}).
##'
##' Confidence intervals are evaluating such that the lower and upper bound correspond to the same likelihood value,
##' aiming at having intervals where lack of coverage is equaly likely due to low or to high bounds.
##' This is performed using a root finding algorithm (\code{\link{uniroot}}).
##' 
##' The constrained likelihood is evaluated as follow:\itemize{
##' \item the confidence interval of a parameter is discretized with \code{maxpt} points. Increasing the confidence level will lead to a larger range of parameter values.
##' \item this parameter is set to each discretization value.
##' \item the other parameters are either set to the (unconstrained) MLE (\code{profile.likelihood=FALSE})
##' or to constrained MLE (\code{profile.likelihood=TRUE}). The latter case is much more computer intensive as it implies re-running the estimation procedure.
##' \item the (restricted) log-likelihood is evaluated.
##' }
##' Since a locally quadratic log-likelihood with an Hessian equivariant in law implies normally distributed estimates (Geyer 2013)
##' it can help trusting confidence intervals and p-values in small samples with a non-normally distributed outcome.
##'
##' @return [profile.likelihood = TRUE/FALSE] A data.frame object containing the log-likelihood for various parameter values. \cr
##' [profile.likelihood = "ci"] A data.frame object containing the REML or ML estimated parameter (\code{"estimate"}),
##' lower and upper bound of the profile-likelihood confidence interval (\code{"lower"}, \code{"upper"}),
##' and the discrepancy in log-likelihood between the bounds found by the root finding algorithm and the requested difference in log-likelihood reflecting the confidence level (\code{"error.lower"}, \code{"error.upper"}).
##'
##' @references Geyer, C. J. (2013). Asymptotics of maximum likelihood without the lln or clt or sample size going to infinity. In Advances in Modern Statistical Theory and Applications: A Festschrift in honor of Morris L. Eaton, pages 1–24. Institute of Mathematical Statistics.
##' 
##' @keywords htest
##' 
##' @examples
##' ### Linear regression ####
##' data(gastricbypassW, package = "LMMstar")
##' e.lmm <- lmm(weight2 ~ weight1 + glucagonAUC1, data = gastricbypassW)
##'
##' #### likelihood along a parameter axis (slice)
##' ## no transformation
##' e.sliceNone <- profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "none")
##' plot(e.sliceNone)
##' ## transformation
##' e.sliceLog <- profile(e.lmm, effects = "all", maxpts = 10, transform.sigma = "log")
##' plot(e.sliceLog)
##'
##' #### profile likelihood (local maxima of the likelihood - crest line)
##' \dontrun{
##' e.pro <- profile(e.lmm, effects = "all", profile.likelihood = TRUE)
##' plot(e.pro)
##' }
##' 
##' #### confidence interval based on profile likelihood
##' \dontrun{
##' e.PLCI <- profile(e.lmm, effects = c("weight1","sigma"), profile.likelihood = "ci")
##' e.PLCI
##' }
##'
##' ### Random intercept model ####
##' ## Data shown in Sahai and Ageel (2000) page 122.
##' ## The analysis of variance: fixed, random, and mixed models. Springer
##' df <- rbind(data.frame(Y = c(7.2, 7.7, 8, 8.1, 8.3, 8.4, 8.4, 8.5, 8.6, 8.7,
##'                              9.1, 9.1, 9.1, 9.8, 10.1, 10.3), type = "Hb SS"),
##'            data.frame(Y = c(8.1, 9.2, 10, 10.4, 10.6, 10.9, 11.1, 11.9, 12, 12.1),
##'                       type = "Hb S/beta"),
##'            data.frame(Y = c(10.7, 11.3, 11.5, 11.6, 11.7, 11.8, 12, 12.1, 12.3,
##'                             12.6, 12.6, 13.3, 13.3, 13.8, 13.9), type = "HB SC")
##' )
##' df$type <- factor(df$type, levels = unique(df$type))
##'
##' ## retrive first column of table I in the original publication
##' ## (doi:  10.1136/bmj.282.6260.283)
##' round(tapply(df$Y,df$type,mean),1)
##' round(tapply(df$Y,df$type,sd),1)
##'
##' e.RI <- lmm(Y ~ (1|type), data = df)
##'
##' #### delta method
##' confint(e.RI, effects = "correlation", df = FALSE) ## 0.76 [0.111; 0.955]
##'  
##' ranef(e.RI, effects = "variance", se = TRUE) ## 3.17 [0.429; 23.423]
##' ranef(e.RI, effects = "variance", se = TRUE, transform=FALSE) ## 3.17 [-3.167; 9.510]
##' ## same as confint(e.RI, effects = "correlation", transform.rho = "cov", df = FALSE)
##'
##' #### profile likelihood
##' \dontrun{
##' plot(profile(e.RI, effects = "correlation", profile.likelihood = TRUE, df = FALSE))
##' plot(profile(e.RI, effects = "correlation", profile.likelihood = TRUE, df = FALSE,
##'              transform.rho = "none", maxpts = seq(0.4,0.99,by=0.01)))
##'
##' profile(e.RI, effects = "correlation", profile.likelihood = "ci")
##' ## 0.760 [0.378; 0.983]
##' profile(e.RI, effects = "correlation", profile.likelihood = "ci", transform.rho = "cov")
##' }

## * profile.lmm (code)
##' @export
profile.lmm <- function(fitted, effects = NULL, profile.likelihood = FALSE,
                        maxpts = NULL, level = 0.95, df = NULL, trace = FALSE,                        
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

    ## *** profile.likelihood
    if(is.character(profile.likelihood) & tolower(profile.likelihood) %in% c("ci","interval","confidence interval","confidence-interval")){
        profile.likelihood <- 2
    }else if(profile.likelihood %in% 0:2 == FALSE){
        stop("Argument \'profile.likelihood\' should be FALSE, TRUE, or \"ci\". \n")
    }
    
    if(fitted$args$control$optimizer!="FS" && profile.likelihood>0){
        stop("Argument \'profile.likelihood\' can only be TRUE or \"ci\" when \"FS\" optimizer is used. \n",
             "Consider adding the argument control = list(optimizer = \"FS\") when fitting the mixed model with lmm. \n")
    }
    if(profile.likelihood == 2){
        delta.likelihood <- stats::qchisq(level,df=1)/2
        target.likelihood <- logLik(fitted) - stats::qchisq(level,df=1)/2

        if(is.null(transform.rho)){
            transform.rho <- "none"
        }
        if(is.null(transform.k) && (!is.null(transform.rho) && transform.rho == "none")){
            transform.k <- "none"
        }
        if(is.null(transform.sigma) && (!is.null(transform.k) && transform.k == "none") && (!is.null(transform.rho) && transform.rho == "none")){
            transform.sigma <- "none"
        } 
    }

    ## *** transformations
    init <- .init_transform(p = NULL, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = fitted$reparametrize$transform.sigma, x.transform.k = fitted$reparametrize$transform.k, x.transform.rho = fitted$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

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
    
    }
    if(any(effects %in% name.p == FALSE)){
        effects <- stats::model.tables(fitted, effects = c("param",effects))$name
    }
    n.effects <- length(effects)

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

    ## *** df
    if(is.null(df)){
        df <- !is.null(fitted$df)
    }

    ## ** prepare
    if(profile.likelihood == 2){
        p.trans <- stats::model.tables(fitted, effects = "all", level = level, 
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names,
                                       options = options)
        name.p.trans <- rownames(p.trans)        
        rownames(p.trans) <- name.p
    }else{
        p <- cbind(estimate = stats::coef(fitted, effects = "all", 
                                          transform.sigma = "none", transform.k = "none", transform.rho = "none", options = options))
        p.trans <- stats::confint(fitted, effects = "all", level = level, df = df,
                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                  transform.names = transform.names, backtransform = FALSE, options = options)
        name.p.trans <- rownames(p.trans)
        rownames(p.trans) <- name.p
    }

    ## ** profile likelihood
    if(trace>1){
        if(profile.likelihood==2){
            cat("Profile likelihood confidence intervals:\n",sep="")
        }else{
            cat("Profile likelihood (",round(2*maxpts)," points):\n",sep="")
        }
    }

    ls.profile <- lapply(1:n.effects, function(iParam){ ## iParam <- 4

        iIndex <- which(name.p == effects[iParam])
        iName <- name.p[iIndex]
        iType <- unname(type.p[iIndex])
        iName.trans <- name.p.trans[iIndex]
        iValue.trans <- p.trans[iIndex,"estimate"]
        iLower.trans <- p.trans[iIndex,"lower"]
        iUpper.trans <- p.trans[iIndex,"upper"]
        if(profile.likelihood!=2){
            iValue <- unname(p[iIndex,"estimate"])
        }

        ## *** evaluate CI
        if(profile.likelihood==2){

            ## Chebyshev's inequality inequality P(|X-muX|>k\sigma) < 1/k^2
            iBound <- .unirootBound(cbind(p.trans[iIndex,], type = iType), level = level,
                                    transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
            iFCT.dlogLik <- function(iiValue){ ## iiValue <- 3.2
                iLMM <- .constrain.lmm(fitted, effects = stats::setNames(iiValue, iName.trans), 
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                       trace = FALSE)
                if(iLMM$opt$cv<=0){
                    iDiff <- +Inf
                }else{
                    iDiff <- target.likelihood - logLik(iLMM)
                }
                return(iDiff)
            }
            ## optim
            resUp <- try(stats::uniroot(f = iFCT.dlogLik, lower = iValue.trans, upper = iBound[2], f.lower = logLik(fitted), extendInt = "upX"))
            resDown <- try(stats::uniroot(f = iFCT.dlogLik, lower = iBound[1], upper = iValue.trans, f.upper = logLik(fitted), extendInt = "downX"))
            
            ## export
            if(inherits(resUp,"try-error")){
                resUp <- list(root = NA, f.root = NA)
            }
            if(inherits(resDown,"try-error")){
                resDown <- list(root = NA, f.root = NA)
            }
            iDF <- data.frame(estimate = iValue.trans, lower = resDown$root, upper = resUp$root, error.lower = resDown$f.root, error.upper = resUp$f.root)
            rownames(iDF) <- iName
            return(iDF)
        }

        ## *** evaluate likelihood
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

        

        iOut <- data.frame(param = iName.trans,
                           type = iType,
                           value = NA,
                           value.trans = seqValue.trans,
                           optimum = seqOptimum,
                           logLik = NA,
                           cv = NA)
        iOut[iIndex.center,"logLik"] <- fitted$logLik
        iOut[iIndex.center,"cv"] <- TRUE
        iOut[iIndex.center,"value"] <- p[iName,"estimate"]

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
                    iResInf <- try(.constrain.lmm(fitted, effects = stats::setNames(seqValue.trans[iIndex.center-iPts], iName.trans), init = iInitInf,
                                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, trace = FALSE))

                    if(!inherits(iResInf,"try-error")){
                        iOut[iIndex.center-iPts, c("logLik","cv","value")] <- c(logLik = iResInf$logLik, cv = iResInf$opt$cv, value = iResInf$param[iName])
                        iInitInf <- iResInf$param
                        if(profile.likelihood>1){
                            iOut[iIndex.center-iPts, name.p] <- iResInf$param
                        }
                    }
                }

                if(iIndex.center+iPts<=length(seqValue.trans)){
                    iResSup <- try(.constrain.lmm(fitted, effects = stats::setNames(seqValue.trans[iIndex.center+iPts], iName.trans), init = iInitSup,
                                                  transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, trace = FALSE))
                    if(!inherits(iResSup,"try-error")){
                        iOut[iIndex.center+iPts, c("logLik","cv","value")] <- c(logLik = iResSup$logLik, cv = iResSup$opt$cv, value = iResSup$param[iName])                    
                        iInitSup <- iResSup$param
                        if(profile.likelihood>1){
                            iOut[iIndex.center+iPts, name.p] <- iResSup$param
                        }
                    }
                }
            }


        }else{
            iArgs.reparam <- list(p = seqValue.trans,
                                  type = rep(rep(iType),length(seqValue.trans))
                                  )
            if((iType == "k" && transform.k %in% c("sd","logsd","var","logvar")) || (iType == "rho" && transform.rho == "cov")){
                iArgs.reparam$p <- c(iArgs.reparam$p,stats::setNames(p.trans[table.param[iIndex,"sigma"],"estimate"],table.param[iIndex,"sigma"]))
                iArgs.reparam$type <- c(iArgs.reparam$type,"sigma")
                iArgs.reparam$sigma <- c(rep(table.param[iIndex,"sigma"],length(seqValue.trans)),NA)
            }
            if(iType == "rho" && transform.rho == "cov"){
                if(is.na(table.param[iIndex,"k.x"]) && is.na(table.param[iIndex,"k.y"])){
                    iArgs.reparam$k.x <- c(rep(table.param[iIndex,"k.x"],length(seqValue.trans)),NA)
                    iArgs.reparam$k.y <- c(rep(table.param[iIndex,"k.y"],length(seqValue.trans)),NA)
                }else{
                    iArgs.reparam$p <- c(iArgs.reparam$p, stats::setNames(p.trans[unlist(table.param[iIndex,c("k.x","k.y")]),"estimate"],table.param[iIndex,c("k.x","k.y")]))
                    iArgs.reparam$type <- c(iArgs.reparam$type,"k","k")
                    iArgs.reparam$sigma <- c(iArgs.reparam$sigma,rep(table.param[iIndex,"sigma"],2))
                    iArgs.reparam$k.x <- c(rep(table.param[iIndex,"k.x"],length(seqValue.trans)),NA,NA,NA)
                    iArgs.reparam$k.y <- c(rep(table.param[iIndex,"k.y"],length(seqValue.trans)),NA,NA,NA)
                }
            }
            
            seqValue <- .reparametrize(p = iArgs.reparam$p,
                                       type = iArgs.reparam$type,
                                       sigma = iArgs.reparam$sigma,
                                       k.x = iArgs.reparam$k.x,
                                       k.y = iArgs.reparam$k.y,                                       
                                       transform.sigma = transform.sigma,
                                       transform.k = transform.k,
                                       transform.rho = transform.rho,
                                       transform.names = FALSE,
                                       inverse = TRUE, Jacobian = FALSE, dJacobian = FALSE)$p[1:length(seqValue.trans)]

            iOut$value <- seqValue
            iOut$logLik[-iIndex.center] <- sapply(seqValue[-iIndex.center], function(iiValue){ ## iiValue <- -0.99984950
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
    ## unique(df.profile$param)

    ## ** export
    if(profile.likelihood != 2){
        df.profile$param <- factor(df.profile$param, levels = unique(df.profile$param))
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
    }
    return(df.profile)
}

## * .unirootBound
##' @description Determine conservative lower bounds and upper bounds for the confidence intervals to initialize the root finding algorithm
##'
##' @param object data.frame containing the estimate, its standard error, and its type (mu,sigma,k,rho).
##' @param level [numeric,0-1] confidence level.
##' @param tol [numeric] tolerance around the interval of definition of the parameter.
##' @noRd
##'
##' @details Based on Chebyshev's inequality inequality P(|X-muX|>k\sigma) < 1/k^2
.unirootBound <- function(object, level, tol = 1e-6,
                          transform.sigma, transform.k, transform.rho){

    bound <- c(object$estimate - object$se * 1/sqrt(1-level), object$estimate + object$se * 1/sqrt(1-level))
    if(object$type == "sigma" && transform.sigma %in% c("none","square")){
        bound[1] <- max(tol, bound[1])
    }else if(object$type == "k" && transform.k %in% c("none","square","sd","var")){
        bound[1] <- max(tol, bound[1])
    }else if(object$type == "rho" && transform.rho == "none"){
        bound <- c(max(-1+tol, bound[1]),min(1-tol, bound[2]))
    }

    return(bound)

}

##----------------------------------------------------------------------
### profile.R ends here
