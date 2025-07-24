### moments.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (09:15) 
## Version: 
## Last-Updated: jul 24 2025 (13:31) 
##           By: Brice Ozenne
##     Update #: 727
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @param robust [0,1,2] 0: model-based s.e. are computed and df are relative to model-based s.e.
##'                       1: robust s.e. are computed but df are relative to model-based s.e.
##'                       2: robust s.e. are computed and df are relative to robust s.e.
##' @noRd

## * moments.lmm
moments.lmm <- function(x, effects = NULL, newdata = NULL, p = NULL,
                        logLik = TRUE, score = TRUE, information = TRUE, vcov = TRUE, df = TRUE,
                        indiv = FALSE, type.information = NULL, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

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
        valid.effects <- c("mean","fixed","variance","correlation","all")
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

    x.param <- stats::model.tables(x, effects = c("param",effects), transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)

    ## *** type.information
    if(is.null(type.information)){
        type.information <- x$args$type.information
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    ## *** transformation & p
    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho,
                            table.param = x$design$param)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    if(is.null(p)){
        theta <- x$param
    }else{
        theta <- init$p
    }
    
    ## ** extract or recompute information
    if(is.null(newdata) && is.null(p) && (indiv == FALSE) && test.notransform && x$args$type.information==type.information){
        keep.name <- stats::setNames(x.param$name, x.param$trans.name)    

        design <- x$design ## useful in case of NA
        out <- x$information[keep.name,keep.name,drop=FALSE]
        if(transform.names){
            dimnames(out) <- list(names(keep.name),names(keep.name))
        }
    }else{
         
        if(!is.null(newdata)){
            design <- stats::model.matrix(x, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- x$design
        }

    }

    ## ** evaluate moments 
    out <- .moments.lmm(value = theta, design = design, time = x$time, method.fit = x$args$method.fit, type.information = type.information,
                        transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                        logLik = logLik, score = score, information = information, vcov = vcov, df = df, indiv = indiv, effects = effects, robust = FALSE,
                        trace = FALSE, precompute.moments = !is.null(x$design$precompute.XX), method.numDeriv = options$method.numDeriv, transform.names = transform.names)

    ## ** export
    return(out)
}

## * .moments.lmm
.moments.lmm <- function(value, design, time, method.fit, type.information,
                         transform.sigma, transform.k, transform.rho,
                         logLik, score, information, vcov, df, indiv, effects, robust,
                         trace, precompute.moments, method.numDeriv, transform.names){

    param.value <- value[design$param$name]
    param.type <- design$param$type
    param.level <- design$param$level
    n.cluster <- length(design$index.cluster)
    df.analytic <- FALSE
    df.numeric <- df
    out <- list()

    ## ** 1- reparametrisation
    if(trace>=1){cat("- reparametrization \n")}
    name.allcoef <- names(param.value)
    index.var <- which(param.type %in% c("sigma","k","rho"))
    if(df.analytic){
        test.d2Omega <- TRUE
    }else if(vcov || information){
        test.d2Omega <- (method.fit == "REML" || type.information == "observed")
    }else{
        test.d2Omega  <- FALSE
    }

    out$reparametrize <- .reparametrize(p = param.value[index.var], type = param.type[index.var], level = param.level[index.var], 
                                        sigma = design$param$sigma[index.var], k.x = design$param$k.x[index.var], k.y = design$param$k.y[index.var],
                                        Jacobian = TRUE, dJacobian = 2*test.d2Omega, inverse = FALSE, ##  2 is necessary to export the right dJacobian
                                        transform.sigma = transform.sigma,
                                        transform.k = transform.k,
                                        transform.rho = transform.rho,
                                        transform.names = TRUE)

    newname.allcoef <- stats::setNames(name.allcoef, name.allcoef)
    if(out$reparametrize$transform==FALSE){
        out$reparametrize$newname <- NULL
        out$reparametrize$Jacobian <- NULL
        out$reparametrize$dJacobian <- NULL
    }else{
        newname.allcoef[names(out$reparametrize$p)] <- out$reparametrize$newname
    }
    if(score || information || vcov || df.analytic){
        type.effects <- c("mu","sigma","k","rho")[c("mean","variance","variance","correlation") %in% effects]
        attr(effects, "original.names") <- names(newname.allcoef[param.type %in% type.effects])
        attr(effects, "reparametrize.names") <- as.character(newname.allcoef[param.type %in% type.effects])        
    }

    ## ** 2- compute partial derivatives regarding the mean and the variance
    if(trace>=1){cat("- residuals \n")}
    out$fitted <- design$mean %*% param.value[colnames(design$mean)]
    out$residuals <- design$Y - out$fitted

    if(precompute.moments){

        wRR <- out$residuals
        if(attr(design$weights, "user-defined")){ 
            wRR <- sweep(wRR, FUN = "*", MARGIN = 1, STATS = sqrt(design$weights))
        } ## otherwise weights are set automatically to 1 but no need to update the residuals
        
        precompute <- list(weights = design$precompute.weights,
                           XX = design$precompute.XX,
                           RR = .precomputeRR(residuals = wRR, pattern = design$vcov$Upattern$name, 
                                              pattern.ntime = stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name),
                                              pattern.cluster = attr(design$vcov$pattern,"list"), index.cluster = design$index.cluster)                           
                           )

        if(score || information || vcov || df.analytic){
            wR <-  out$residuals
            if(attr(design$weights, "user-defined")){
                wR <- sweep(wR, FUN = "*", MARGIN = 1, STATS = design$weights)
            } ## otherwise weights are set automatically to 1 but no need to update the residuals
            precompute$XR  <-  .precomputeXR(X = design$mean, residuals = wR, pattern = design$vcov$Upattern$name,
                                             pattern.ntime = stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name),
                                             pattern.cluster = attr(design$vcov$pattern,"list"), index.cluster = design$index.cluster)
        }
        
    }else{
        precompute <- list()
    }

    if(trace>=1){cat("- Omega \n")}
    out$Omega <- .calc_Omega(object = design$vcov, param = param.value, simplify = FALSE)

    ## choleski decomposition
    Omega.chol <- lapply(out$Omega,function(iO){try(chol(iO),silent=TRUE)})
    if(any(sapply(Omega.chol,inherits,"try-error"))){
        index.error <- which(sapply(Omega.chol,inherits,"try-error"))
        attr(out,"error") <- c("Residuals variance-covariance matrix is not positive definite. Original error message:\n",
                               unique(unlist(Omega.chol[index.error])))
    }
    ## inverse
    out$OmegaM1 <- lapply(1:length(out$Omega),function(iP){ ## iP <- 1
        if(inherits(Omega.chol[[iP]],"try-error")){
            iOut <- try(solve(out$Omega[[iP]],silent=FALSE)) ## matrix may be negative definite, i.e., invertible but with negative eigenvalues
            attr(iOut,"vectorize") <- as.vector(iOut)
            iDet <- det(iOut)
            if(!is.na(iDet) & iDet>0){ ## handle negative determinant
                attr(iOut,"logdet") <- log(iDet)
            }else{
                attr(iOut,"logdet") <- NA
            }
        }else{
            iOut <- chol2inv(Omega.chol[[iP]])
            attr(iOut,"vectorize") <- as.vector(iOut)
            attr(iOut,"logdet") <- -2*sum(log(diag(Omega.chol[[iP]])))
        }
        return(iOut)
    })
    names(out$OmegaM1) <- names(out$Omega)

    ## log(sapply(out$OmegaM1,det))
    if(score || information || vcov || df.analytic){
        if(trace>=1){cat("- dOmega \n")}
        out$dOmega <- .calc_dOmega(object = design$vcov, param = param.value, Omega = out$Omega, 
                                   Jacobian = out$reparametrize$Jacobian,
                                   transform.sigma = transform.sigma,
                                   transform.k = transform.k,
                                   transform.rho = transform.rho)
    }

    if(test.d2Omega){
        if(trace>=1){cat("- d2Omega \n")}
        out$d2Omega <- .calc_d2Omega(object = design$vcov, param = param.value, Omega = out$Omega, dOmega = out$dOmega, 
                                     Jacobian = out$reparametrize$Jacobian, dJacobian = out$reparametrize$dJacobian,
                                     transform.sigma = transform.sigma,
                                     transform.k = transform.k,
                                     transform.rho = transform.rho)
    }

    ## ** 3- precompute
    ## *** require the full information whenever the information is not block diagonal
    ## all the matrix is need in order to get the inverse (vcov) 
    if((vcov && (method.fit=="REML"||type.information=="observed"))  || df.analytic){
        effects2 <- c("mean","variance","correlation")
        attr(effects2, "original.names") <- names(newname.allcoef)
        attr(effects2, "reparametrize.names") <- as.character(newname.allcoef)
    }else if(score || information || vcov || df.analytic){
        effects2 <- effects
    }

    ## *** matrix product between the residual variance-covariance matrix and its derivative
    precompute$Omega <- .precomputeOmega(precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega,
                                         effects = effects2, pair.vcov = design$vcov$pair.vcov,
                                         REML = method.fit=="REML", type.information = type.information,
                                         logLik = logLik, score = (score || (vcov && robust)), information = information, vcov = vcov, df = df.analytic)

    if(method.fit == "REML" && !is.null(precompute$XX)){
        precompute$REML <- .precomputeREML(precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, precompute = precompute, effects = effects2,
                                           logLik = logLik, score = (score || (vcov && robust)), information = information, vcov = vcov, df = df.analytic)
    }

    ## ** 4- compute likelihood derivatives

    ## *** log-likelihood
    if(logLik){
        if(trace>=1){cat("- log-likelihood \n")}
        out$logLik <- .logLik(X = design$mean, residuals = out$residuals, precision = out$OmegaM1, weights = design$weights,
                              pattern = design$vcov$pattern, index.cluster = design$index.cluster, 
                              indiv = indiv, REML = method.fit=="REML", precompute = precompute)
    }

    ## *** score
    if(score || (vcov && robust)){ 
        if(trace>=1){cat("- score \n")}

        Mscore <- .score(X = design$mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega, weights = design$weights, 
                         pattern = design$vcov$pattern, index.cluster = design$index.cluster, name.allcoef = name.allcoef,
                         indiv = indiv || (vcov && robust), REML = method.fit=="REML", effects = effects2, precompute = precompute)

        if(score){
            if(indiv){
                out$score <- Mscore[,attr(effects, "original.names"),drop=FALSE]
            }else{
                if(robust){
                    out$score <- colSums(Mscore[,attr(effects, "original.names"),drop=FALSE])
                }else{
                    out$score <- Mscore[attr(effects, "original.names")]
                }
            }
            if(transform.names && length(out$reparametrize$newname)>0){
                if(indiv){
                    colnames(out$score) <- newname.allcoef[colnames(out$score)]
                }else{
                    names(out$score) <- newname.allcoef[names(out$score)]
                }
            }
            attr(out$score,"message") <- attr(Mscore,"message")
        }
    }

    ## *** information
    if(information || vcov){
        if(trace>=1){cat("- information \n")}
        Minfo <- .information(X = design$mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, weights = design$weights, 
                              pattern = design$vcov$pattern, index.cluster = design$index.cluster, name.allcoef = name.allcoef,
                              pair.meanvcov = design$vcov$pair.meanvcov, pair.vcov = design$vcov$pair.vcov,
                              indiv = indiv && information, REML = (method.fit=="REML"), type.information = type.information, effects = effects2, 
                              precompute = precompute)

        if(information){
            if(indiv){
                out$information <- Minfo[,attr(effects, "original.names"),attr(effects, "original.names"),drop=FALSE]
            }else{
                out$information <- Minfo[attr(effects, "original.names"),attr(effects, "original.names"),drop=FALSE]
            }
            attr(out$information, "type.information") <- type.information
            attr(out$information,"message") <- attr(Minfo,"message")
            if(transform.names && length(out$reparametrize$newname)>0){
                if(indiv){
                    dimnames(out$information) <- list(NULL, attr(effects, "reparametrize.names"),attr(effects, "reparametrize.names"))
                }else{
                    dimnames(out$information) <- list(attr(effects, "reparametrize.names"),attr(effects, "reparametrize.names"))
                }
            }
        }
    }

    ## *** variance-covariance
    if(vcov){
        if(trace>=1){cat("- variance-covariance \n")}

        if(indiv && information){
            Minfo <- apply(Minfo, MARGIN = 2:3, sum)
        }

        if(is.invertible(Minfo, cov2cor = TRUE)){
            Mvcov <- solve(Minfo)
        }else{
            warning("Singular or nearly singular information matrix. \n")
            Mvcov <- try(solve(Minfo), silent = TRUE)
            if(inherits(Mvcov,"try-error")){
                Mvcov <- NA*Mvcov
            }
        }
        if(robust & !inherits(vcov,"try-error")){ 
            Mvcov <- Mvcov %*% crossprod(Mscore) %*% Mvcov
            attr(Mvcov,"message") <- attr(Mscore,"message")
        }

        if(vcov){
            out$vcov <- Mvcov[attr(effects, "original.names"),attr(effects, "original.names"),drop=FALSE]
            attr(out$vcov, "type.information") <- type.information
            attr(out$vcov, "robust") <- robust>0
            attr(out$vcov, "message") <- attr(Mvcov,"message")
            if(transform.names && length(out$reparametrize$newname)>0){
                dimnames(out$vcov) <- list(attr(effects, "reparametrize.names"),attr(effects, "reparametrize.names"))
            }
        }
    }

    if(df){
        if(trace>=1){cat("- degrees-of-freedom \n")}
        if(df.analytic){
            ## out$df2 <- .df_analytic(residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, Upattern.ncluster = Upattern.ncluster, vcov = out$vcov,
            ##                         pattern = design$vcov$pattern, index.clusterTime = design$index.time, index.cluster = design$index.cluster,
            ##                         name.varcoef = design$vcov$Upattern$param, name.allcoef = name.allcoef,
            ##                         pair.meanvarcoef = design$param$pair.meanvarcoef, pair.varcoef = design$vcov$pair.varcoef,
            ##                         indiv = indiv, REML = (method.fit=="REML"), type.information = type.information, name.effects = name.effects, robust = robust, diag = TRUE,
            ##                         precompute = precompute)
        }else if(df.numeric){
            ## require vcov for all parameters to compute df
            effects.all <- c("mean", "variance", "correlation")
            attr(effects.all, "original.names") <- names(newname.allcoef)
            attr(effects.all, "original.output") <- attr(effects, "original.names")
            attr(effects.all, "reparametrize.names") <- as.character(newname.allcoef)
            attr(effects.all, "reparametrize.output") <- attr(effects, "reparametrize.names")

            out$df <- .df_numDeriv(reparametrize = out$reparametrize,
                                   value = param.value, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   effects = effects.all, robust = (robust==2), ## if robust is 1 then robust s.e. are computed but df are relative to model-based s.e.
                                   precompute.moments = precompute.moments, method.numDeriv = method.numDeriv)
        }
        ## range(pmin(out$df2,10000)-pmin(out$df,10000))
        out$dVcov <- attr(out$df,"dVcov")
        attr(out$df,"dVcov") <- NULL
        if(transform.names && length(out$reparametrize$newname)>0){
            names(out$df) <- newname.allcoef[names(out$df)]
            dimnames(out$dVcov) <- list(newname.allcoef[dimnames(out$dVcov)[[1]]], newname.allcoef[dimnames(out$dVcov)[[2]]], newname.allcoef[dimnames(out$dVcov)[[3]]])
        }
    }

    ## ** 4- export
    return(out)
}

##----------------------------------------------------------------------
### moments.R ends here
