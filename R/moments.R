### moments.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (09:15) 
## Version: 
## Last-Updated: Jul  1 2022 (09:42) 
##           By: Brice Ozenne
##     Update #: 394
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * moment.lmm
.moments.lmm <- function(value, design, time, method.fit, type.information,
                         transform.sigma, transform.k, transform.rho,
                         logLik, score, information, vcov, df, indiv, effects, robust,
                         trace, precompute.moments, method.numDeriv, transform.names){

    param.value <- value[design$param$name]
    param.type <- design$param$type
    param.level <- design$param$level
    out <- list()

    ## ** 1- reparametrisation
    if(trace>=1){cat("- reparametrization \n")}
    name.allcoef <- names(param.value)
    index.var <- which(param.type %in% c("sigma","k","rho"))
    if(df){
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

    if(score || information || vcov || df){
        type.effects <- c("mu","sigma","k","rho")[c("mean","variance","variance","correlation") %in% effects]
        attr(effects, "original.names") <- names(newname.allcoef[param.type %in% type.effects])
        attr(effects, "reparametrize.names") <- as.character(newname.allcoef[param.type %in% type.effects])
    }
    
    ## ** 2- compute partial derivatives regarding the mean and the variance
    if(trace>=1){cat("- residuals \n")}
    out$fitted <- design$mean %*% param.value[colnames(design$mean)]
    out$residuals <- design$Y - out$fitted

    if(precompute.moments){
        if(is.null(design$weights)){
            wRR <- out$residuals
            wR <-  out$residuals
            Upattern.ncluster <- stats::setNames(design$vcov$X$Upattern$n.cluster,design$vcov$X$Upattern$name)
        }else{
            wR <- sweep(out$residuals, FUN = "*", MARGIN = 1, STATS = design$weights)
            wRR <- sweep(out$residuals, FUN = "*", MARGIN = 1, STATS = sqrt(design$weights))
            Upattern.ncluster <- tapply(tapply(design$weight,design$index.cluster,unique),design$vcov$X$pattern.cluster$pattern,sum)[design$vcov$X$Upattern$name]
        }
        weights.cluster <- NULL
        scale.cluster <- NULL

        precompute <- list(XX = design$precompute.XX,
                           RR = .precomputeRR(residuals = wRR, pattern = design$vcov$X$Upattern$name, 
                                              pattern.ntime = stats::setNames(design$vcov$X$Upattern$n.time, design$vcov$X$Upattern$name),
                                              pattern.cluster = design$vcov$X$Upattern$index.cluster, index.cluster = design$index.cluster)                           
                           )
        if(score || information || vcov || df){
            precompute$XR  <-  .precomputeXR(X = design$precompute.XX$Xpattern, residuals = wR, pattern = design$vcov$X$Upattern$name,
                                             pattern.ntime = stats::setNames(design$vcov$X$Upattern$n.time, design$vcov$X$Upattern$name),
                                             pattern.cluster = design$vcov$X$Upattern$index.cluster, index.cluster = design$index.cluster)
        }
        
        
    }else{
        precompute <- NULL
        Upattern.ncluster <- NULL
        if(is.null(design$weights)){
            weights.cluster <- stats::setNames(rep(1, design$cluster$n), design$cluster$levels)
        }else{
            weights.cluster <- stats::setNames(design$weights[sapply(design$index.cluster,"[[",1)], design$cluster$levels)
        }
        if(is.null(design$scale.Omega)){
            scale.cluster <- stats::setNames(rep(1, design$cluster$n), design$cluster$levels)
        }else{
            scale.cluster <- stats::setNames(design$scale.Omega[sapply(design$index.cluster,"[[",1)], design$cluster$levels)
        }
    }

    if(trace>=1){cat("- Omega \n")}
    out$Omega <- .calc_Omega(object = design$vcov, param = param.value, keep.interim = TRUE)

    ## choleski decomposition
    Omega.chol <- lapply(out$Omega,function(iO){try(chol(iO),silent=TRUE)})
    if(any(sapply(Omega.chol,inherits,"try-error"))){
        index.error <- which(sapply(Omega.chol,inherits,"try-error"))
        attr(out,"error") <- c("Residuals variance-covariance matrix is not positive definite. Original error message:\n",
                               unique(unlist(Omega.chol[index.error])))
    }
    ## inverse
    out$OmegaM1 <- lapply(Omega.chol,function(iChol){
        if(inherits(iChol,"try-error")){return(iChol)}else{return(chol2inv(iChol))}
    })
    ## determinant
    attr(out$OmegaM1,"logdet") <- sapply(Omega.chol, function(iChol){
        if(inherits(iChol,"try-error")){return(NA)}else{return(-2*sum(log(diag(iChol))))}
    })
    ## log(sapply(out$OmegaM1,det))
    if(score || information || vcov || df){
        if(trace>=1){cat("- dOmega \n")}
        out$dOmega <- .calc_dOmega(object = design$vcov, param = param.value, Omega = out$Omega,
                                   Jacobian = out$reparametrize$Jacobian,
                                   transform.sigma = transform.sigma,
                                   transform.k = transform.k,
                                   transform.rho = transform.rho)

        attr(out$dOmega, "ls.dOmega_OmegaM1") <- stats::setNames(lapply(design$vcov$X$Upattern$name, function(iPattern){
            lapply(out$dOmega[[iPattern]], function(iM){
                if(inherits(out$OmegaM1[[iPattern]],"try-error")){return(NA)}else{iM %*% out$OmegaM1[[iPattern]]}
            })
        }), design$vcov$X$Upattern$name)
        attr(out$dOmega, "ls.OmegaM1_dOmega_OmegaM1") <- stats::setNames(lapply(design$vcov$X$Upattern$name, function(iPattern){ ## iPattern <- "1:1"
            lapply(attr(out$dOmega, "ls.dOmega_OmegaM1")[[iPattern]], function(iM){
                if(inherits(out$OmegaM1[[iPattern]],"try-error")){return(NA)}else{out$OmegaM1[[iPattern]] %*% iM}
            })
        }), design$vcov$X$Upattern$name)
        attr(out$dOmega, "dOmega_OmegaM1") <- lapply(attr(out$dOmega, "ls.dOmega_OmegaM1"), function(iO){ do.call(cbind, lapply(iO,as.numeric)) })
        attr(out$dOmega, "OmegaM1_dOmega_OmegaM1") <- lapply(attr(out$dOmega, "ls.OmegaM1_dOmega_OmegaM1"), function(iO){ do.call(cbind, lapply(iO,as.numeric)) })
    }

    if(test.d2Omega){
        if(trace>=1){cat("- d2Omega \n")}
        out$d2Omega <- .calc_d2Omega(object = design$vcov, param = param.value, Omega = out$Omega, dOmega = out$dOmega, 
                                     Jacobian = out$reparametrize$Jacobian, dJacobian = out$reparametrize$dJacobian,
                                     transform.sigma = transform.sigma,
                                     transform.k = transform.k,
                                     transform.rho = transform.rho)
    }
    ## param.value
    ##     MM <- numDeriv::jacobian(func = function(iP){
    ##         as.vector(.calc_Omega(object = design$vcov, param = iP)[[1]])
    ##     }, x = param.value)
    ##     matrix(MM[,7],4,4)
    ##     matrix(MM[,8],4,4)

    ## MM <- numDeriv::jacobian(func = function(iP){
    ##     as.vector(.calc_dOmega(object = design$vcov, param = iP, 
    ##                            Jacobian = out$reparametrize$Jacobian)[[1]]$sigma)
    ## }, x = param.value)
    ## matrix(MM[,7],4,4)
    ## matrix(MM[,8],4,4)

    ## ** 3- compute likelihood derivatives
    if(logLik){
        if(trace>=1){cat("- log-likelihood \n")}
        out$logLik <- .logLik(X = design$mean, residuals = out$residuals, precision = out$OmegaM1,
                              Upattern.ncluster = Upattern.ncluster, weights = weights.cluster, scale.Omega = scale.cluster,
                              index.variance = design$vcov$X$pattern.cluster$pattern, time.variance = design$index.clusterTime, index.cluster = design$index.cluster, 
                              indiv = indiv, REML = method.fit=="REML", precompute = precompute)
    }

    if(score){ 
        if(trace>=1){cat("- score \n")}
        out$score <- .score(X = design$mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega,
                            Upattern.ncluster = Upattern.ncluster, weights = weights.cluster, scale.Omega = scale.cluster,
                            index.variance = design$vcov$X$pattern.cluster$pattern, time.variance = design$index.clusterTime, index.cluster = design$index.cluster,
                            name.allcoef = name.allcoef,
                            indiv = indiv, REML = method.fit=="REML", effects = effects, precompute = precompute)

        if(transform.names && length(out$reparametrize$newname)>0){
            if(indiv){
                colnames(out$score) <- newname.allcoef[colnames(out$score)]
            }else{
                names(out$score) <- newname.allcoef[names(out$score)]
            }
        }
    }

    if(information || vcov || df){## needed for finding the names of the coefficients and getting the variance-covariance matrix
        if(trace>=1){cat("- information \n")}
        if(vcov || df){ ## compute the full information otherwise the inverse (i.e. vcov) may not be the correct one
            effects2 <- c("mean","variance","correlation")
            attr(effects2, "original.names") <- names(newname.allcoef)
            attr(effects2, "reparametrize.names") <- as.character(newname.allcoef)
        }else{
            effects2 <- effects
        }
        Minfo <- .information(X = design$mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega,
                              Upattern.ncluster = Upattern.ncluster, weights = weights.cluster, scale.Omega = scale.cluster,
                              index.variance = design$vcov$X$pattern.cluster$pattern, time.variance = design$index.clusterTime, index.cluster = design$index.cluster,
                              name.allcoef = name.allcoef,
                              pair.meanvarcoef = attr(design$param, "pair.meanvarcoef"), pair.varcoef = design$vcov$X$pair.varcoef,
                              indiv = indiv, REML = (method.fit=="REML"), type.information = type.information, effects = effects2, robust = robust,
                              precompute = precompute)

        if(information){
            if(indiv){
                out$information <- Minfo[,attr(effects, "original.names"),attr(effects, "original.names"),drop=FALSE]
            }else{
                out$information <- Minfo[attr(effects, "original.names"),attr(effects, "original.names"),drop=FALSE]
            }
            attr(out$information, "type.information") <- type.information
            attr(out$information, "robust") <- robust
            if(transform.names && length(out$reparametrize$newname)>0){
                if(indiv){
                    dimnames(out$information) <- list(NULL, attr(effects, "reparametrize.names"),attr(effects, "reparametrize.names"))
                }else{
                    dimnames(out$information) <- list(attr(effects, "reparametrize.names"),attr(effects, "reparametrize.names"))
                }
            }
        }
    }
    if(vcov || df){
        if(trace>=1){cat("- variance-covariance \n")}
        if(robust && method.fit=="REML"){
            keep.cols <- intersect(names(which(rowSums(!is.na(Minfo))>0)),names(which(rowSums(!is.na(Minfo))>0)))
            Mvcov <- NA*Minfo
            Mvcov[keep.cols,keep.cols] <- solve(Minfo[keep.cols,keep.cols,drop=FALSE])
        }else{
            Mvcov <- solve(Minfo)
        }
        if(vcov){
            out$vcov <- Mvcov[attr(effects, "original.names"),attr(effects, "original.names"),drop=FALSE]
            if(transform.names && length(out$reparametrize$newname)>0){
                dimnames(out$vcov) <- list(attr(effects, "reparametrize.names"),attr(effects, "reparametrize.names"))
            }
        }
    }

    if(df){
        if(trace>=1){cat("- degrees of freedom \n")}
        ## system.time(
        out$df <- .df_numDeriv(value = param.value, reparametrize = out$reparametrize,
                               design = design, time = time, method.fit = method.fit, type.information = type.information,
                               transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, effects = effects2, 
                               robust = robust, diag = TRUE,
                               precompute.moments = precompute.moments, method.numDeriv = method.numDeriv)
        ## )
        ## system.time(
        ## out$df2 <- .df_analytic(residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, Upattern.ncluster = Upattern.ncluster, vcov = out$vcov,
        ##                         index.variance = design$vcov$X$pattern.cluster, time.variance = design$index.time, index.cluster = design$index.cluster,
        ##                         name.varcoef = design$vcov$X$Upattern$param, name.allcoef = name.allcoef,
        ##                         pair.meanvarcoef = design$param$pair.meanvarcoef, pair.varcoef = design$vcov$pair.varcoef,
        ##                         indiv = indiv, REML = (method.fit=="REML"), type.information = type.information, name.effects = name.effects, robust = robust, diag = TRUE,
        ##                         precompute = precompute)
        ## )
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
