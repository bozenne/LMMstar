### moments.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (09:15) 
## Version: 
## Last-Updated: aug  8 2023 (19:30) 
##           By: Brice Ozenne
##     Update #: 676
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

    ## ** 0- extract information
    param.value <- value[design$param$name]
    param.type <- design$param$type
    param.level <- design$param$level

    design.param <- design$param
    n.cluster <- length(design$index.cluster)
    pair.vcov <- design$vcov$pair.vcov
    pair.meanvcov <- design$vcov$pair.meanvcov
    out <- list()

    ## ** 1- reparametrisation
    if(trace>=1){cat("- reparametrization \n")}

    test.d2Omega <- (df || ((vcov || information) && (method.fit == "REML" || type.information == "observed")))
    
    name.allcoef <- names(param.value)
    index.var <- which(param.type %in% c("sigma","k","rho"))
    out$reparametrize <- .reparametrize(p = param.value[index.var], type = param.type[index.var], level = param.level[index.var], 
                                        sigma = design.param$sigma[index.var], k.x = design.param$k.x[index.var], k.y = design.param$k.y[index.var],
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
    test.freecoef <- is.na(design$param[match(design$param$name,names(newname.allcoef)),"constraint"]) 
    newname.freecoef <- newname.allcoef[test.freecoef]

    ## restrict effects to the ones meaningful for the model or those asked by the user
    if(score || information || vcov || df){
        param2effect <- sapply(unique(param.type), switch, "mu" = "mean", "sigma" = "variance", "k" = "variance", "rho" = "correlation")

        effects.all <-  names(param2effect)
        attr(effects.all, "original.names") <- names(newname.allcoef)
        attr(effects.all, "reparametrize.names") <- as.character(newname.allcoef)

        effects <- intersect(effects, param2effect)
        type.effects <- names(param2effect)[param2effect %in% effects]
        attr(effects, "original.names") <- names(newname.allcoef[param.type %in% type.effects])
        attr(effects, "reparametrize.names") <- as.character(newname.allcoef[param.type %in% type.effects])
    }

    ## ** 2- compute residual variance-covariance and its derivatives
    if(trace>=1){cat("- residuals \n")}
    out$fitted <- design$mean %*% param.value[colnames(design$mean)]
    out$residuals <- design$Y - out$fitted

    if(trace>=1){cat("- Omega \n")}
    out$Omega <- .calc_Omega(object = design$vcov, param = param.value, keep.interim = TRUE)

    if((score && (any(effects %in% c("variance","correlation")) || method.fit == "REML")) || information || vcov || df){
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
    
    if(df && (type.information == "observed" || method.fit == "REML")){
        if(trace>=1){cat("- d3Omega \n")}
        design$vcov$triplet.vcov <- .precomputeTripletParam(design$vcov)
        
        out$d3Omega <- .calc_d3Omega(object = design$vcov, param = param.value, Omega = out$Omega, 
                                     transform.sigma = transform.sigma,
                                     transform.k = transform.k,
                                     transform.rho = transform.rho)        
  }

    ## ** 3- Precompute useful quantities
    
    ## *** design matrix and residuals
    if(precompute.moments){
        precompute <- list(wXX = design$precompute.wXX)
        weights.cluster <- NULL ## weights are integrated in the residuals/number of clusters
        scale.cluster <- NULL ## does not handle scaled Omega
        
        if(is.null(design$weights)){
            wRR <- out$residuals
            wR <-  out$residuals
            Upattern.ncluster <- stats::setNames(design$vcov$Upattern$n.cluster,design$vcov$Upattern$name)
        }else{
            wR <- sweep(out$residuals, FUN = "*", MARGIN = 1, STATS = design$weights)
            wRR <- sweep(out$residuals, FUN = "*", MARGIN = 1, STATS = sqrt(design$weights))
            Upattern.ncluster <- tapply(tapply(design$weight,design$index.cluster,unique),design$vcov$pattern.cluster$pattern,sum)[design$vcov$Upattern$name]
        }

        precompute$wRR <- .precomputeRR(residuals = wRR, pattern = design$vcov$Upattern$name, 
                                        pattern.ntime = stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name),
                                        pattern.cluster = attr(design$vcov$pattern,"list"), index.cluster = design$index.cluster)                           
        if(!is.null(design$precompute.wXX) && (score || information || vcov || df)){
            precompute$wXR  <-  .precomputeXR(X = design$precompute.wXX$Xpattern, residuals = wR, pattern = design$vcov$Upattern$name,
                                              pattern.ntime = stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name),
                                              pattern.cluster = attr(design$vcov$pattern,"list"), index.cluster = design$index.cluster)
        }

        attr(precompute,"moments") <- TRUE
    }else{
        if(is.null(design$weights)){
            weights.cluster <- rep(1, n.cluster)
            Upattern.ncluster <- stats::setNames(design$vcov$Upattern$n.cluster,design$vcov$Upattern$name)
        }else{
            weights.cluster <- design$weights[sapply(design$index.cluster,"[[",1)]
            Upattern.ncluster <- tapply(tapply(design$weight,design$index.cluster,unique),design$vcov$pattern.cluster$pattern,sum)[design$vcov$Upattern$name]
        }
        if(is.null(design$scale.Omega)){
            scale.cluster <- rep(1, n.cluster)
        }else{
            scale.cluster <- design$scale.Omega[sapply(design$index.cluster,"[[",1)]
        }
        precompute <- list()
        attr(precompute,"moments") <- FALSE
    }

    ## *** Omega and its derivatives
    if(method.fit == "REML" && !precompute.moments){
        precompute$wX <- lapply(attr(design$vcov$pattern,"list"), function(iCs){
            lapply(iCs, function(iC){sqrt(weights.cluster[iC]) * design$mean[design$index.cluster[[iC]],,drop=FALSE]})
        })
    }
browser()
    resPrecomputeO <- .precomputeOmega(precompute = precompute,
                                       Omega = out$Omega, dOmega = out$dOmega, d2Omega = out$d2Omega, 
                                       method.fit = method.fit, type.information = type.information,
                                       score = score, information = information, vcov = vcov, df = df)
    precompute[names(resPrecomputeO)] <- resPrecomputeO

    
    ## ** 4- compute likelihood derivatives
    if(logLik){
        if(trace>=1){cat("- log-likelihood \n")}
        out$logLik <- .logLik(X = design$mean, residuals = out$residuals, precompute = precompute,
                              Upattern.ncluster = Upattern.ncluster, Upattern.ntime = design$vcov$Upattern$n.time, weights = weights.cluster, scale.Omega = scale.cluster,
                              pattern = design$vcov$pattern, index.cluster = design$index.cluster, 
                              indiv = indiv, REML = method.fit=="REML")
    }

    if(score){ 
        if(trace>=1){cat("- score \n")}
        
        out$score <- .score(X = design$mean, residuals = out$residuals, precompute = precompute,
                            Upattern.ncluster = Upattern.ncluster, weights = weights.cluster, scale.Omega = scale.cluster,
                            pattern = design$vcov$pattern, index.cluster = design$index.cluster, 
                            indiv = indiv, REML = method.fit=="REML", effects = effects)
        
        if(transform.names && length(out$reparametrize$newname)>0){
            if(indiv){
                colnames(out$score) <- newname.allcoef[colnames(out$score)]
            }else{
                names(out$score) <- newname.allcoef[names(out$score)]
            }
        }
    }

    if(information || vcov || df){
        if(trace>=1){cat("- information \n")}

        if(((vcov || robust) && type.information=="observed") || df){ 
            ## must compute the entire matrix to be able to inverse it 
            effects2 <- effects.all
        }else{
            ## can inverse the mean block / the variance block separately
            effects2 <- effects
        }
        Minfo <- .information(X = design$mean, residuals = out$residuals, precompute = precompute, pair.vcov = pair.vcov,
                              Upattern.ncluster = Upattern.ncluster, weights = weights.cluster, scale.Omega = scale.cluster,
                              pattern = design$vcov$pattern, index.cluster = design$index.cluster, 
                              indiv = indiv, REML = (method.fit=="REML"), type.information = type.information, effects = effects2, robust = robust)

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

    ## ** 5- Quantify uncertainty based on likelihood derivatives
    if(vcov || df){
        if(trace>=1){cat("- variance-covariance \n")}
        
        if(robust && method.fit=="REML"){
            keep.cols <- intersect(names(which(rowSums(!is.na(Minfo))>0)),names(which(rowSums(!is.na(Minfo))>0)))
            Mvcov <- NA*Minfo
            if(is.invertible(Minfo[keep.cols,keep.cols,drop=FALSE], cov2cor = TRUE)){
                Mvcov[keep.cols,keep.cols] <- solve(Minfo[keep.cols,keep.cols,drop=FALSE])
            }else{
                warning("Singular or nearly singular information matrix. \n")
                test <- try(solve(Minfo[keep.cols,keep.cols,drop=FALSE]), silent = TRUE)
                if(inherits(Mvcov,"try-error")){
                    df <- FALSE
                }else{
                    Mvcov[keep.cols,keep.cols] <- test
                }
            }
        }else{
            if(is.invertible(Minfo, cov2cor = TRUE)){
                Mvcov <- solve(Minfo)
            }else{
                warning("Singular or nearly singular information matrix. \n")
                Mvcov <- try(solve(Minfo), silent = TRUE)
                if(inherits(Mvcov,"try-error")){
                    Mvcov <- NA*Minfo
                    df <- FALSE
                }
            }
            
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
        ##                         pattern = design$vcov$pattern, index.clusterTime = design$index.time, index.cluster = design$index.cluster,
        ##                         name.varcoef = design$vcov$Upattern$param, name.allcoef = name.allcoef,
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

    ## ** 5- export
    return(out)
}

##----------------------------------------------------------------------
### moments.R ends here
