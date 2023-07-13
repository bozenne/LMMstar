### ranef.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: May 26 2022 (11:18) 
## Version: 
## Last-Updated: jul 13 2023 (13:10) 
##           By: Brice Ozenne
##     Update #: 374
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
##' @param effects [character] should the estimated random effects (\code{"mean"}) or the estimated variance of the random effects (\code{"variance"}) be output?
##' @param p [numeric vector] value of the model coefficients to be used. Only relevant if differs from the fitted values.
##' @param ci [logical] should standard error and confidence intervals be evaluated using a delta method?
##' Will slow down the execution of the function.
##' @param format [character] should each type of random effect be output in a data.frame (\code{format="long"})
##' or using a list where the random effect are grouped by nesting factors (\code{format="wide"}).
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
##' data(gastricbypassL, package = "LMMstar")
##' 
##' ## random intercept
##' e.RI <- lmm(weight ~ time + (1|id), data = gastricbypassL)
##' ranef(e.RI)
##' 
##' ## nested random effects
##'
##' ## crossed random effects
##'
##'
##' #### random effect calculation in lme4 ####
##' if(require(lme4)){
##' data(Penicillin, package = "lme4")
##'
##' e.lmer <- lmer(diameter ~ (1|plate) + (1|sample), data = Penicillin)
##' y.lmer <- as.matrix(getME(e.lmer, "y") )
##' X.lmer <- as.matrix(getME(e.lmer, "X") )
##' Z.lmer <- as.matrix(getME(e.lmer, "Z") )
##' Lambda.lmer <- as.matrix(getME(e.lmer,"Lambda"))
##' sigma2.lmer <- getME(e.lmer,"sigma")^2
##' 
##' G.lmer <- sigma2.lmer*crossprod(t(Lambda.lmer))
##' Omega.lmer <- Z.lmer %*% G.lmer %*% t(Z.lmer) + sigma2.lmer*diag(NROW(Z.lmer))
##' epsilon.lmer <- Penicillin$diameter  - predict(e.lmer, re.form = ~0)
##'
##' ## solve system 17 in lme4 vignette
##' LHS <- c(t(Lambda.lmer) %*% t(Z.lmer) %*% y.lmer, t(X.lmer) %*% y.lmer)
##' B1 <- crossprod(Z.lmer %*% Lambda.lmer) + diag(1, NCOL(Z.lmer))
##' B2 <- t(X.lmer) %*% Z.lmer %*% Lambda.lmer
##' B3 <- crossprod(X.lmer) 
##' RHS <- rbind(cbind(B1,t(B2)), cbind(B2, B3))
##' hat <- solve(RHS) %*% cbind(LHS)
##' tail(hat,1) ## intercept
##' GS <- Lambda.lmer %*% head(hat,-1) ## random effects
##' GS - do.call(rbind,ranef(e.lmer))
##'
##' ## inuitive formula
##' test <- G.lmer %*% t(Z.lmer) %*% solve(Omega.lmer) %*% epsilon.lmer
##' (test / GS)
##'
##' ## in lmm
##' e.lmm <- lmm(diameter ~ (1|plate) + (1|sample), data = Penicillin, df = FALSE)
##' Tau <- ranef(e.lmm, effects = "variance")
##' OmegaM1epsilon.lmm <- residuals(e.lmm, type = "normalized2")
##' range(OmegaM1epsilon.lmm - solve(Omega.lmer) %*% epsilon.lmer)
##'
##' sum(OmegaM1epsilon.lmm[Penicillin$plate=="a"])*Tau[1]
##' sum(OmegaM1epsilon.lmm[Penicillin$sample=="F"])*Tau[2]
##' ranef(e.lmm)
##' ranef(e.lmer)
##' }

## * ranef.lmm (code)
##' @export
ranef.lmm <- function(object, effects = "mean", ci = FALSE, transform = (effects=="variance"),
                      p = NULL, format = "long", keep.strata = FALSE, ...){



    ## ** normalize user input
    mycall <- match.call()
    if(!inherits(object$design$vcov,"RE")){
        stop("Cannot estimate random effects linear mixed models defined by covariance structure (argument \'structure\'). \n",
             "Consider adding random effects in the argument \'formula\' instead. \n")
    }
    effects <- match.arg(effects, c("mean","variance"))
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

        e.ranef <- nlme::ranef(object, effects = effects, ci = FALSE, p = p, format = format, keep.strata = keep.strata)
        e.delta <- lava::estimate(object, f = function(newp){
            iE <- nlme::ranef(object, effects = effects, ci = FALSE, p = newp, format = format)
            return(iE$estimate)
        }, df = df)
            
        if(transform){ ## recompute only CIs (backtransforming the se is not exact)
            eTrans.delta <- lava::estimate(object, f = function(newp){
                iE <- nlme::ranef(object, effects = effects, ci = FALSE, p = newp, format = format)
                iE[iE$type=="variance","estimate"] <- log(iE[iE$type=="variance","estimate"])
                iE[iE$type=="relative","estimate"] <- atanh(iE[iE$type=="relative","estimate"])
                return(iE$estimate)
            }, df = df)
            ## absolute
            e.delta$lower[e.ranef$type=="variance"] <- exp(eTrans.delta$lower[e.ranef$type=="variance"])
            e.delta$upper[e.ranef$type=="variance"] <- exp(eTrans.delta$upper[e.ranef$type=="variance"])
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
    var.strata <- object$design$vcov$name$strata
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

    ## missing values
    index.na <- object$index.na

    ## all vars
    var.all <- unique(lava::manifest(object))
    
    ## ** converting correlation parameters into random effect variance
    cumtau <- coef(object, p = p, effects = "correlation", transform.rho = "cov", transform.names = FALSE)
    cumtau.strata <- tapply(cumtau,param.strata,identity, simplify = FALSE)

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
    if(effects == "variance"){
        sigma2 <- coef(object, effects = "variance", transform.sigma = "square")
        if(n.strata==1){
            if(format=="long"){
                out <- rbind(data.frame(variable = rownames(varRE),
                                        XXstrataXX = colnames(varRE),
                                        type = "variance",
                                        estimate = varRE[,1]),
                             data.frame(variable = rownames(varRE),
                                        XXstrataXX = colnames(varRE),
                                        type = "relative",
                                        estimate = varRE[,1]/sigma2)
                             )
            }else if(format=="wide"){
                out <- data.frame(variable = rownames(varRE),
                                  XXstrataXX = colnames(varRE),
                                  variance = varRE[,1],
                                  relative= varRE[,1]/sigma2)                
            }
            rownames(out) <- NULL
        }else{
            browser()
            out.wide <- data.frame(variable = rownames(varRE),
                                   varRE, check.names = FALSE)
            stats::reshape(out.wide, direction = "long", idvar = "variable", timevar = object$strata,
                           times = colnames(varRE))
        }
        if(keep.strata==FALSE){
            out$XXstrataXX <- NULL
        }else if(n.strata>1){
            browser()
        }
        return(out)
    }

    ## ** extract normalized residuals
    ## head(stats::residuals(object, p = p, keep.data = TRUE, type = "response", format = "long"))
    df.epsilon <- stats::residuals(object, p = p, keep.data = TRUE, type = "normalized2", format = "long")
    if(length(index.na)>0){
        df.epsilon <- df.epsilon[-index.na,,drop=FALSE]
    }
    
    ## ** estimate random effects
    out <- do.call(rbind,lapply(name.hierarchy, function(iName){ ## iName <- name.hierarchy[1]        
        iOut <- do.call(rbind,by(df.epsilon, df.epsilon[[iName]], function(iDF){ ## iDF <- df.epsilon[13:18,]
            iTau <- varRE[iName,iDF$XXstrataXX[1]]
            iOut <- cbind(variable = unname(iName),
                          XXstrataXX = iDF$XXstrataXX[1],
                          level = unname(iDF[1,iName,drop=FALSE]),
                          estimate = unname(sum(iDF$r.normalized)*iTau))
            return(iOut)
        }))
        return(iOut[match(unique(df.epsilon[[iName]]),iOut$level),,drop=FALSE])
    }))            
    
    ## ** export
    rownames(out) <- NULL
    if(format == "wide"){
        
        out <- stats::setNames(lapply(infoRanef$hierarchy, function(iH){
            iOut <- out[out$variable == iH[1],"estimate",drop=FALSE]
            rownames(iOut) <- out$level
            if(infoRanef$cross){
                browser()
            }
            return(iOut)
        }), sapply(infoRanef$hierarchy,"[",1))
        
    }else if(format == "long"){
        if(keep.strata == FALSE){
            out$XXstrataXX <- NULL
        }else if(n.strata>1){
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
