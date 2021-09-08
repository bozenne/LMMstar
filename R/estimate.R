### estimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 20 2021 (23:25) 
## Version: 
## Last-Updated: sep  8 2021 (13:09) 
##           By: Brice Ozenne
##     Update #: 117
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

.estimate <- function(design, time, method.fit, type.information,
                      transform.sigma, transform.k, transform.rho,
                      precompute.moments, init, n.iter, tol.score, tol.param, trace){

   
    if(!precompute.moments){
        stop("Only implemented when option \'precompute.moments\' is TRUE")
    }
    options <- LMMstar.options()
    if(is.null(n.iter)){
        n.iter <- options$param.optimizer["n.iter"]
    }
    if(is.null(tol.score)){
        tol.score <- options$param.optimizer["tol.score"]
    }
    if(is.null(tol.param)){
        tol.param <- options$param.optimizer["tol.param"]
    }
    if(is.null(trace)){
        trace <- FALSE
    }

    ## ** prepare
    param.mu <- design$param$mu
    param.sigma <- design$param$sigma
    param.k <- design$param$k
    param.rho <- design$param$rho
    param.Omega <- c(param.sigma,param.k,param.rho)
    param.name <- c(param.mu,param.Omega)
    n.param <- length(param.name)
    pattern <- design$X.var$pattern
    precompute.XY <- design$precompute.XY
    precompute.XX <- design$precompute.XX$pattern
    key.XX <- design$precompute.XX$key
    
    ## ** intialization
    if(is.null(init)){
        param.value <- setNames(rep(NA, n.param),param.name)

        start.OmegaM1 <- stats::setNames(lapply(pattern, function(iPattern){
            diag(1, nrow = length(design$X.var$index.time[[iPattern]]), ncol = length(design$X.var$index.time[[iPattern]]))
        }), pattern)
        param.value[param.mu] <- .estimateGLS(OmegaM1 = start.OmegaM1, pattern = pattern, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX)
        iResiduals.long <- design$Y - design$X.mean %*% param.value[param.mu]
        iResiduals.wide <- reshape2::dcast(data = data.frame(residuals = iResiduals.long, cluster = design$index.cluster, time = design$index.time),
                                           formula = cluster~time, value.var = "residuals")
        iM.rescor <- stats::cor(iResiduals.wide[,-1,drop=FALSE], use = "pairwise")
        diag(iM.rescor) <- NA
        iRescor <- mean(iM.rescor, na.rm = TRUE)
        iResvar <- mean(apply(iResiduals.wide[,-1,drop=FALSE], MARGIN = 2, FUN = stats::var, na.rm = TRUE))

        param.value[param.sigma] <- sqrt(iResvar)

        if(length(param.k)>0){
            param.value[param.k] <- 1
        }

        if(length(param.rho)>0){
            param.value[param.rho] <- iRescor
        }

    }else{
        param.value <- init[param.name]
    }
    ## microbenchmark(test =     .estimateGLS(OmegaM1 = start.OmegaM1, pattern = pattern, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX),
    ##                lm.fit = lm.fit(y = ncgsL$cholest[!is.na(ncgsL$cholest)], x = model.matrix(~time+highdose.time, data=ncgsL)[!is.na(ncgsL$cholest),]),
    ##                lm = lm(cholest~time+highdose.time, data=ncgsL)
    ##                )


    if(trace>1){
        cat("\nInitialization:\n")
        print(param.value)
    }
    
    ## ** loop
    cv <- FALSE
    param.valueM1 <- NULL
    if(trace>1){
        cat("\nLoop:\n")
    }

    for(iIter in 1:n.iter){ ## iIter <- 1
        param.valueM1 <- param.value
        
        ## *** moments
        outMoments <- .moments.lmm(value = param.value, design = design, time = time, method.fit = method.fit, type.information = "expected",
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   logLik = FALSE, score = TRUE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("variance","correlation"), robust = FALSE,
                                   trace = FALSE, precompute.moments = TRUE, transform.names = FALSE)

        if(all(abs(outMoments$score)<tol.score) && (iIter==1 || all(abs(param.valueM1 - param.value)<tol.param))){
            cv <- TRUE
            break
        }

        ## *** variance estimate
        param.newvalue.trans <- outMoments$reparametrize$p + as.double(outMoments$score %*% solve(outMoments$information))
        param.value[param.Omega] <- .reparametrize(param.newvalue.trans,
                                                   type = design$param$type[names(param.newvalue.trans)],
                                                   strata = design$param$strata[names(param.newvalue.trans)],
                                                   time.levels = time$levels, time.k = design$param$time.k, time.rho = design$param$time.rho,
                                                   Jacobian = FALSE, dJacobian = FALSE, inverse = TRUE,
                                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)$p
        

        ## *** mean estimate
        iOmega <- .calc_Omega(object = design$X.var, param = param.value, keep.interim = TRUE)
        param.value[param.mu] <- .estimateGLS(OmegaM1 = stats::setNames(lapply(iOmega, solve), names(iOmega)),
                                              pattern = pattern, precompute.XY = precompute.XY, precompute.XX = precompute.XX,
                                              key.XX = key.XX)

        
        if(trace==1){
            cat("*")
        }else if(trace>1){
            M.print <- rbind(estimate = param.value,
                             diff = param.value - param.valueM1,
                             score = c(rep(0, length(param.mu)),outMoments$score))
            rownames(M.print) <- paste0(rownames(M.print),".",iIter)
            print(M.print)
        }
        
    }

    if(trace==1){
        cat("\n")
    }else if(trace>1){
        if(cv){
            if(iIter==1){
                cat("Convergence after ",iIter," iteration: max score=",max(abs(outMoments$score)),"\n", sep = "")
            }else{
                cat("Convergence after ",iIter," iterations: max score=",max(abs(outMoments$score))," | max change in coefficient=",max(abs(param.valueM1 - param.value)),"\n", sep = "")
            }
        }else{
            if(iIter==1){
                cat("No convergence after ",iIter," iteration: max score=",max(abs(outMoments$score)),"\n")
            }else{
                cat("No convergence after ",iIter," iterations: max score=",max(abs(outMoments$score))," | max change in coefficient= ",max(abs(param.valueM1 - param.value)),"\n", sep = "")
            }
        }
    }

    ## ** export
    return(list(estimate = param.value,
                previous.estimate = param.valueM1,
                score = outMoments$score,
                n.iter = iIter,
                cv = cv))
}

.estimateGLS <- function(OmegaM1, pattern, precompute.XY, precompute.XX, key.XX){

    name.param <- colnames(key.XX)
    n.param <- length(name.param)
    max.key <- key.XX[n.param,n.param]
    numerator <- matrix(0, nrow = n.param, ncol = 1)
    denominator <- matrix(0, nrow = n.param, ncol = n.param)
    
    for(iPattern in pattern){ ## iPattern <- pattern[1]
        iVec.Omega <- as.double(OmegaM1[[iPattern]])
        iTime2 <- length(iVec.Omega)
        numerator  <- numerator + t(iVec.Omega %*%  matrix(precompute.XY[[iPattern]], nrow = iTime2, ncol = n.param))
        denominator  <- denominator + as.double(iVec.Omega %*%  matrix(precompute.XX[[iPattern]], nrow = iTime2, ncol = max.key))[key.XX]
    }

    out <- solve(denominator) %*% numerator
    return(stats::setNames(as.double(out), name.param))
}


##----------------------------------------------------------------------
### estimate.R ends here
