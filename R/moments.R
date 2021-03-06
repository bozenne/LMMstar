### moments.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (09:15) 
## Version: 
## Last-Updated: Jul  8 2021 (15:33) 
##           By: Brice Ozenne
##     Update #: 67
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

.moments.lmm <- function(param, design, time, method.fit, type.information,
                         transform.sigma, transform.k, transform.rho,
                         logLik, score, information, vcov, df, indiv, effects, robust,
                         trace, precompute.moments, method.numDeriv, transform.names){


    out <- list()

    ## ** 1- reparametrisation
    if(trace>=1){cat("- reparametrization \n")}
    name.allcoef <- names(param$value)
    
    index.var <- which(param$type %in% c("sigma","k","rho"))
    if(df){
        test.d2Omega <- TRUE
    }else if(vcov || information){
        test.d2Omega <- (method.fit == "REML" || type.information == "observed")
    }else{
        test.d2Omega  <- FALSE
    }
    out$reparametrize <- .reparametrize(p = param$value[index.var], type = param$type[index.var], strata = param$strata[index.var], time.levels = time$levels,
                                        time.k = design$param$time.k, time.rho = design$param$time.rho,
                                        Jacobian = TRUE, dJacobian = 2*test.d2Omega, inverse = FALSE, ##  2 is necessary to export the right dJacobian
                                        transform.sigma = transform.sigma,
                                        transform.k = transform.k,
                                        transform.rho = transform.rho,
                                        transform.names = TRUE)

    if(out$reparametrize$transform==FALSE){
        out$reparametrize$newname <- NULL
        out$reparametrize$Jacobian <- NULL
        out$reparametrize$dJacobian <- NULL
    }else{
        newname.allcoef <- stats::setNames(name.allcoef, name.allcoef)
        newname.allcoef[names(out$reparametrize$p)] <- out$reparametrize$newname
    }

    ## ** 2- compute partial derivatives regarding the mean and the variance
    if(trace>=1){cat("- residuals \n")}
    out$fitted <- design$X.mean %*% param$value[colnames(design$X.mean)]
    out$residuals <- design$Y - out$fitted
    if(precompute.moments){
        precompute <- list(XX = design$precompute.XX,
                           RR = .precomputeRR(residuals = out$residuals, pattern = design$X.var$pattern,
                                              pattern.time = design$X.var$index.time, pattern.cluster = attr(design$X.var$cluster, "index.byPattern"), index.cluster = attr(design$index.cluster,"sorted"))                           
                           )
        if(score || information || vcov || df){
            precompute$XR  <-  .precomputeXR(X = design$precompute.XX$Xpattern, residuals = out$residuals, pattern = design$X.var$pattern,
                                             pattern.time = design$X.var$index.time, pattern.cluster = attr(design$X.var$cluster, "index.byPattern"), index.cluster = attr(design$index.cluster,"sorted"))
        }

    }else{
        precompute <- NULL
    }

    if(trace>=1){cat("- Omega \n")}
    out$Omega <- .calc_Omega(object = design$X.var, param = param$value, keep.interim = TRUE)
    out$OmegaM1 <- lapply(out$Omega,solve)
    
    if(score || information || vcov || df){
        if(trace>=1){cat("- dOmega \n")}
        out$dOmega <- .calc_dOmega(object = design$X.var, param = param$value, type = param$type, Omega = out$Omega,
                                   Jacobian = out$reparametrize$Jacobian)
    }

    if(test.d2Omega){
        if(trace>=1){cat("- d2Omega \n")}
        out$d2Omega <- .calc_d2Omega(object = design$X.var, param = param$value, type = param$type,
                                     Omega = out$Omega, dOmega = out$dOmega, pair = design$param$pair.varcoef,
                                     Jacobian = out$reparametrize$Jacobian, dJacobian = out$reparametrize$dJacobian)
    }

    ## ** 3- compute likelihood derivatives
    if(logLik){
        if(trace>=1){cat("- log-likelihood \n")}
        out$logLik <- .logLik(X = design$X.mean, residuals = out$residuals, precision = out$OmegaM1,
                              index.variance = design$X.var$cluster, time.variance = design$index.time, index.cluster = design$index.cluster, 
                              indiv = indiv, REML = method.fit=="REML", precompute = precompute)
    }

    if(score){
        if(trace>=1){cat("- score \n")}
        out$score <- .score(X = design$X.mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega,
                            index.variance = design$X.var$cluster, time.variance = design$index.time, index.cluster = design$index.cluster,
                            name.varcoef = design$X.var$param, name.allcoef = name.allcoef,
                            indiv = indiv, REML = method.fit=="REML", effects = effects, precompute = precompute)

        if(transform.names && length(out$reparametrize$newname)>0){
            if(indiv){
                colnames(out$score) <- newname.allcoef[colnames(out$score)]
            }else{
                names(out$score) <- newname.allcoef[names(out$score)]
            }
        }
    }

    if(information || vcov){
        if(trace>=1){cat("- information \n")}
        out$information <- .information(X = design$X.mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, 
                                        index.variance = design$X.var$cluster, time.variance = design$index.time, index.cluster = design$index.cluster,
                                        name.varcoef = design$X.var$param, name.allcoef = name.allcoef,
                                        pair.meanvarcoef = design$param$pair.meanvarcoef, pair.varcoef = design$param$pair.varcoef,
                                        indiv = indiv, REML = (method.fit=="REML"), type.information = type.information, effects = effects, robust = robust,
                                        precompute = precompute)
        attr(out$information, "type.information") <- type.information
        attr(out$information, "robust") <- robust

        if(transform.names && length(out$reparametrize$newname)>0){
            if(indiv){
                dimnames(out$information) <- list(NULL,newname.allcoef[dimnames(out$information)[[2]]],newname.allcoef[dimnames(out$information)[[3]]])
            }else{
                dimnames(out$information) <- list(newname.allcoef[rownames(out$information)],newname.allcoef[colnames(out$information)])
            }
        }

    }

    if(vcov){
        if(trace>=1){cat("- variance-covariance \n")}
        out$vcov <- solve(out$information)
    }

    if(df){
        if(trace>=1){cat("- degrees of freedom \n")}
        out$df <- .df(param = param, reparametrize = out$reparametrize,
                      design, time, method.fit, type.information,
                      transform.sigma, transform.k, transform.rho, effects, robust = robust, diag = TRUE,
                      precompute.moments = precompute.moments, method.numDeriv = method.numDeriv)
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
