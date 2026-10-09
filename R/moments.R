### moments.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (09:15) 
## Version: 
## Last-Updated: okt  8 2026 (15:42) 
##           By: Brice Ozenne
##     Update #: 949
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @param df [TRUE,FALSE] should have an attribute method indicating how the degrees of freedom are to be computed.
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
                        trace = FALSE, method.numDeriv = options$method.numDeriv, transform.names = transform.names)

    ## ** export
    return(out)
}

## * .moments.lmm
.moments.lmm <- function(value, design, time, method.fit, type.information,
                         transform.sigma, transform.k, transform.rho,
                         logLik, score, information, vcov, df, indiv, effects, robust,
                         trace, method.numDeriv, transform.names){
    score <- information <- vcov <- df  <- TRUE
    
    out <- list()
    test.d2Omega <- df || ((vcov || information) & (method.fit == "REML" || type.information == "observed"))
    test.d3Omega <- df && (is.null(attr(df,"method")) || attr(df,"method")=="analytic")
    name.allcoef <- design$param$name
    nameparam.vcov <- design$param[design$param$type!="mu" & is.na(design$param$constraint),"name"]

    ## ** 1- compute partial derivatives regarding the mean and the variance
    if(trace>=1){cat("- residuals \n")}
    out$fitted <- design$mean %*% value[colnames(design$mean)]
    out$residuals <- design$Y - out$fitted

    wRR <- cbind(out$residuals)
    if(!is.null(design$weights.likelihood)){        
        wRR <- sweep(wRR, FUN = "*", MARGIN = 1, STATS = sqrt(design$weights.likelihood))
    }
    if(!is.null(design$weights.Omega)){        
        wRR <- sweep(wRR, FUN = "*", MARGIN = 1, STATS = sqrt(design$weights.Omega))
    }
    
    precompute <- list(weights = design$precompute.weights,
                       XX = design$precompute.XX,
                       RR = .precomputeRR(residuals = wRR, pattern = design$vcov$Upattern$name, 
                                          pattern.ntime = stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name),
                                          pattern.cluster = design$vcov$Upattern$index.cluster, index.cluster = design$index.cluster)                           
                       )

    if(score || information || vcov || df){

        wR <- out$residuals
        if(!is.null(design$weights.likelihood)){        
            wR <- sweep(out$residuals, FUN = "*", MARGIN = 1, STATS = design$weights.likelihood)
        }
        if(!is.null(design$weights.Omega)){        
            wR <- sweep(out$residuals, FUN = "*", MARGIN = 1, STATS = design$weights.Omega)
        }
        precompute$XR  <-  .precomputeXR(X = design$mean, residuals = wR, pattern = design$vcov$Upattern$name,
                                         pattern.ntime = stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name),
                                         pattern.cluster = design$vcov$Upattern$index.cluster, index.cluster = design$index.cluster)
    }
    
    if(trace>=1){cat("- Omega \n")}
    out$Omega <- .calc_Omega(object = design$vcov, param = value,
                             transform.sigma = transform.sigma,
                             transform.k = transform.k,
                             transform.rho = transform.rho,
                             simplify = FALSE)

    ## choleski decomposition
    Omega.chol <- lapply(out$Omega,function(iO){try(chol(iO),silent=TRUE)})
    if(any(sapply(Omega.chol,inherits,"try-error"))){
        index.error <- which(sapply(Omega.chol,inherits,"try-error"))
        attr(out,"error") <- c("Residuals variance-covariance matrix is not positive definite. Original error message:\n",
                               unique(unlist(Omega.chol[index.error])))
    }

    ## inverse
    out$OmegaM1 <- lapply(names(out$Omega),function(iPattern){ ## iP <- 1
        if(!is.null(design$weights.Omega)){
            iPattern.cluster <- design$vcov$Upattern[design$vcov$Upattern$name==iPattern,"index.cluster"][[1]]
            logdet_weights.Omega <- sapply(design$index.cluster[iPattern.cluster], function(iIndex){2*log(prod(design$weights.Omega[iIndex]))})
        }else{
            logdet_weights.Omega <- 0
        }
        if(inherits(Omega.chol[[iPattern]],"try-error")){
            iOut <- try(solve(out$Omega[[iPattern]],silent=FALSE)) ## matrix may be negative definite, i.e., invertible but with negative eigenvalues
            csiDet <- det(iOut)
            if(!is.na(iDet) & iDet>0){ ## handle negative determinant
                attr(iOut,"logdet") <- c(log(iDet),logdet_weights.Omega)
            }else{
                attr(iOut,"logdet") <- c(NA,logdet_weights.Omega)
            }
        }else{
            iOut <- chol2inv(Omega.chol[[iPattern]])
            attr(iOut,"logdet") <- c(-2*sum(log(diag(Omega.chol[[iPattern]]))),logdet_weights.Omega)
        }
        
        return(iOut)
    })
    names(out$OmegaM1) <- names(out$Omega)
        browser()

    ## log(sapply(out$OmegaM1,det))
    if(score || information || vcov || df){
        if(trace>=1){cat("- dOmega \n")}
        out$dOmega <- .calc_dOmega(object = design$vcov, param = value, Omega = out$Omega, 
                                   transform.sigma = transform.sigma,
                                   transform.k = transform.k,
                                   transform.rho = transform.rho)
    }

    if(test.d2Omega){
        if(trace>=1){cat("- d2Omega \n")}
        out$d2Omega <- .calc_d2Omega(object = design$vcov, param = value, Omega = out$Omega, 
                                     transform.sigma = transform.sigma,
                                     transform.k = transform.k,
                                     transform.rho = transform.rho)
    }

    if(test.d3Omega){
        if(trace>=1){cat("- d3Omega \n")}
        triplet <- .meanCovTriplet(structure = design$vcov, X = design$mean, index.cluster = design$index.cluster, index.clusterTime = design$index.clusterTime)

        out$d3Omega <- .calc_d3Omega(object = design$vcov, param = value, Omega = out$Omega, triplet = triplet$vcov3,
                                     transform.sigma = transform.sigma,
                                     transform.k = transform.k,
                                     transform.rho = transform.rho)
    }

    ## ** 2- precompute
    ## *** require the full information whenever the information is not block diagonal
    ## all the matrix is need in order to get the inverse (vcov)
    if((vcov && (method.fit=="REML"||type.information=="observed"))  || is.null(attr(df,"method")) || attr(df,"method")=="analytic"){
        effects2 <- c("mean","variance","correlation")
    }else{
        effects2 <- effects
    }

    ## *** matrix product between the residual variance-covariance matrix and its derivative
    precompute$Omega <- .precomputeOmega(precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, d3Omega = out$d3Omega,
                                         effects = effects2, pair.vcov = design$vcov$pair.vcovvcov, triplet.vcov = triplet$vcov3,
                                         REML = method.fit=="REML", type.information = type.information,
                                         logLik = logLik, score = (score || (vcov && robust)), information = information, vcov = vcov, df = df)

    if(method.fit == "REML" && !is.null(precompute$XX)){
        precompute$REML <- .precomputeREML(precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, d3Omega = out$d3Omega, 
                                           effects = effects2, param.vcov = nameparam.vcov, pair.vcov = design$vcov$pair.vcovvcov, triplet.vcov = triplet$vcov3, precompute = precompute, 
                                           logLik = logLik, score = (score || (vcov && robust)), information = information, vcov = vcov, df = df)
    }

    ## ** 3- compute likelihood derivatives

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
                out$score <- Mscore[,name.allcoef,drop=FALSE]
            }else{
                if(robust){
                    out$score <- colSums(Mscore[,name.allcoef,drop=FALSE])
                }else{
                    out$score <- Mscore[name.allcoef]
                }
            }
            attr(out$score,"message") <- attr(Mscore,"message")
        }
    }

    ## *** information
    if(information || vcov){
        if(trace>=1){cat("- information \n")}
        Minfo <- .information(X = design$mean, residuals = out$residuals, precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, weights = design$weights, 
                              pattern = design$vcov$pattern, index.cluster = design$index.cluster, name.allcoef = name.allcoef, pair.vcov = design$vcov$pair.vcovvcov,
                              indiv = indiv && information, REML = (method.fit=="REML"), type.information = type.information, effects = effects2, 
                              precompute = precompute)

        if(information){
            if(indiv){
                out$information <- Minfo[,name.allcoef,name.allcoef,drop=FALSE]
            }else{
                out$information <- Minfo[name.allcoef,name.allcoef,drop=FALSE]
            }
            attr(out$information, "type.information") <- type.information
            attr(out$information,"message") <- attr(Minfo,"message")            
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
            out$vcov <- Mvcov[name.allcoef,name.allcoef,drop=FALSE]
            attr(out$vcov, "type.information") <- type.information
            attr(out$vcov, "robust") <- robust>0
            attr(out$vcov, "message") <- attr(Mvcov,"message")            
        }
    }
browser()
    if(df){
        if(trace>=1){cat("- degrees-of-freedom \n")}
        if(is.null(attr(df,"method")) || attr(df,"method")=="analytic"){
            out$df <- .df_analytic(param = value, residuals = out$residuals, Omega = out$Omega, precision = out$OmegaM1, dOmega = out$dOmega, d2Omega = out$d2Omega, d3Omega = out$d3Omega,
                                   vcov = out$vcov, Upattern.ncluster = Upattern.ncluster, name.allcoef = name.allcoef,
                                   REML = (method.fit=="REML"), type.information = type.information, name.effects = name.effects, robust = robust, diag = TRUE,
                                   precompute = precompute, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        }else if(attr(df,"method")=="numeric"){
            out$df <- .df_numDeriv(reparametrize = out$reparametrize,
                                   value = param.value, design = design, time = time, method.fit = method.fit, type.information = type.information,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   effects = effects2, robust = (robust==2), ## if robust is 1 then robust s.e. are computed but df are relative to model-based s.e.
                                   method.numDeriv = method.numDeriv)
        }

        out$dVcov <- attr(out$df,"dVcov")
        attr(out$df,"dVcov") <- NULL        
    }

    ## ** 4- export
    return(out)
}

##----------------------------------------------------------------------
### moments.R ends here
