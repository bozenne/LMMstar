### estimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 20 2021 (23:25) 
## Version: 
## Last-Updated: Jun 21 2021 (22:00) 
##           By: Brice Ozenne
##     Update #: 51
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

.estimate <- function(param, design, time, method.fit, type.information,
                      transform.sigma, transform.k, transform.rho,
                      precompute.moments, init, n.iter, tol){

   
    if(!precompute.moments){
        stop("Only implemented when option \'precompute.moments\' is TRUE")
    }

    ## ** prepare
    param.name <- names(param$type)
    n.param <- length(param.name)
    param.mu <- param.name[param$type == "mu"]
    param.Omega <- param.name[param$type != "mu"]
    pattern <- design$X.var$pattern
    precompute.XY <- design$precompute.XY
    precompute.XX <- design$precompute.XX$pattern
    key.XX <- design$precompute.XX$key
    
    ## ** intialization
    iEstimate <- stats::setNames(rep(NA, n.param), param.name)
    if(is.null(init)){
        start.OmegaM1 <- stats::setNames(lapply(pattern, function(iPattern){
            diag(1, nrow = length(design$X.var$index.time[[iPattern]]), ncol = length(design$X.var$index.time[[iPattern]]))
        }), pattern)
        iEstimate[param.mu] <- .estimateGLS(OmegaM1 = start.OmegaM1, pattern = pattern, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX)
        iResiduals.long <- design$Y - design$X.mean %*% iEstimate[param.mu]
        iResiduals.wide <- reshape2::dcast(data = data.frame(residuals = iResiduals.long, cluster = design$index.cluster, time = design$index.time),
                                           formula = cluster~time, value.var = "residuals")
        iRescor <- mean(stats::cor(iResiduals.wide[,-1,drop=FALSE], use = "pairwise"))
        iResvar <- mean(apply(iResiduals.wide[,-1,drop=FALSE], MARGIN = 2, FUN = stats::var, na.rm = TRUE))
        
        if(transform.sigma %in% c("none")){
            iEstimate[param$type=="sigma"] <- sqrt(iResvar)
        }else if(transform.sigma %in% c("square")){
            iEstimate[param$type=="sigma"] <- iResvar
        }else if(transform.sigma %in% c("log")){
            iEstimate[param$type=="sigma"] <- log(sqrt(iResvar))
        }else if(transform.sigma %in% c("logsquare")){
            iEstimate[param$type=="sigma"] <- log(iResvar)
        }else{
            stop("No automatic initialization implemented for transform.sigma=",transform.sigma,".\n")
        }

        if(transform.k %in% c("none","sd")){
            iEstimate[param$type=="k"] <- 1
        }else if(transform.k %in% c("square","var")){
            iEstimate[param$type=="k"] <- 1
        }else if(transform.k %in% c("log","logsd")){
            iEstimate[param$type=="k"] <- 0
        }else if(transform.k %in% c("logsquare","logvar")){
            iEstimate[param$type=="k"] <- 0
        }else{
            stop("No automatic initialization implemented for transform.sigma=",transform.sigma,".\n")
        }

        if(transform.rho %in% c("none")){
            iEstimate[param$type=="rho"] <- iRescor
        }else if(transform.rho %in% c("atanh")){
            iEstimate[param$type=="rho"] <- atanh(iRescor)
        }else{
            stop("No automatic initialization implemented for transform.rho=",transform.rho,".\n")
        }

    }
    ## microbenchmark(test =     .estimateGLS(OmegaM1 = start.OmegaM1, pattern = pattern, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX),
    ##                lm.fit = lm.fit(y = ncgsL$cholest[!is.na(ncgsL$cholest)], x = model.matrix(~time+highdose.time, data=ncgsL)[!is.na(ncgsL$cholest),]),
    ##                lm = lm(cholest~time+highdose.time, data=ncgsL)
    ##                )
    ## ** loop
    iParam <- param
    cv <- FALSE
    print(iEstimate)
        
    for(iIter in 1:n.iter){ ## iIter <- 1
        
        ## *** moments
        iParam$value <- iEstimate
        outMoments <- .moments.lmm(param = iParam, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   logLik = FALSE, score = TRUE, information = TRUE, vcov = FALSE, df = FALSE, indiv = FALSE, effects = c("variance","correlation"), robust = FALSE,
                                   trace = FALSE, precompute.moments = TRUE, transform.names = FALSE)

        if(all(abs(outMoments$score)<tol)){
            cv <- TRUE
            next
        }
        
        ## *** variance estimate
        iEstimate[param.Omega] <- .estimateOmega(param = iEstimate[param.Omega], score = outMoments$score, information = outMoments$information)

        ## *** mean estimate
        iOmega <- .calc_Omega(object = design$X.var, param = iEstimate, keep.interim = TRUE)
        iEstimate[param.mu] <- .estimateGLS(OmegaM1 = stats::setNames(lapply(iOmega, solve), names(iOmega)), pattern = pattern, precompute.XY = precompute.XY, precompute.XX = precompute.XX, key.XX = key.XX)

        print(rbind(estimate = iEstimate,
                    score = c(rep(NA, length(param.mu)),outMoments$score)))
        
    }

    ## ** export
    return(list(estimate = iEstimate,
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

.estimateOmega <- function(param, score, information){
    return(param + score %*% solve(information))
}

##----------------------------------------------------------------------
### estimate.R ends here
