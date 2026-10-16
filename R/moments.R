### moments.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jun 18 2021 (09:15) 
## Version: 
## Last-Updated: okt 16 2026 (11:51) 
##           By: Brice Ozenne
##     Update #: 1140
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
    x.transform <- transform(x)
    init <- .init_transform(p = p, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x.transform$sigma, x.transform.k = x.transform$k, x.transform.rho = x.transform$rho,
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
                        updateGLS = FALSE, logLik = logLik, score = score, information = information, vcov = vcov, df = df, indiv = indiv, effects = effects, robust = FALSE,
                        trace = FALSE, method.numDeriv = options$method.numDeriv, transform.names = transform.names)

    ## ** export
    return(out)
}

## * .moments.lmm
.moments.lmm <- function(value, design, method.fit, type.information,
                         transform.sigma, transform.k, transform.rho,
                         updateGLS, logLik, score, information, vcov, df, indiv, effects, robust,
                         trace, method.numDeriv, transform.names){

    ## ** prepare
    out <- list()
    test.d2Omega <- df || ((vcov || information) & (method.fit == "REML" || type.information == "observed"))
    test.d3Omega <- df && (is.null(attr(df,"method")) || attr(df,"method")=="analytic")
    name.allcoef <- design$param$name
    attr(name.allcoef,"type") <- design$param$type
    nameparam.vcov <- design$param[design$param$type!="mu" & is.na(design$param$constraint),"name"]
    nameparam.mu <- design$param[design$param$type=="mu" & is.na(design$param$constraint),"name"]

    ## require the full information whenever the information is not block diagonal
    ## the whole matrix is needed in order to get the inverse (vcov)
    if((vcov && type.information=="observed")  || (df && attr(df,"method")=="analytic")){
        effects2 <- c("mean","variance","correlation")
    }else{
        effects2 <- effects
    }

    ## prepare for likelihood evaluation
    precompute <- list(weights = design$precompute.weights,
                       XX = design$precompute.XX)

    ## triplet of variance-covariance coefficients
    if(df && (is.null(attr(df,"method")) || attr(df,"method")=="analytic")){
        design$vcov$triplet.vcov <- .tripletVcov(structure = design$vcov, X = design$mean, index.cluster = design$index.cluster, index.clusterTime = design$index.clusterTime)
    }

    ## ** 1- compute partial derivatives regarding the variance
    if(trace>=1){cat("- Omega \n")}
    out$Omega <- .calc_Omega(object = design$vcov, param = value,
                             transform.sigma = transform.sigma,
                             transform.k = transform.k,
                             transform.rho = transform.rho,
                             simplify = FALSE)

    ## precompute matrix product between the residual variance-covariance matrix and its derivative
    precompute$Omega <- .precomputeOmega(object = design$vcov, Omega = out$Omega, precompute = precompute, multiplier = design$multiplier,
                                         param = value, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                         effects = effects2, REML = method.fit=="REML", type.information = type.information,
                                         logLik = logLik, score = (score || (vcov && robust)), information = information, vcov = vcov, df = df)
    out$OmegaM1 <- precompute$Omega$OmegaM1

    if(method.fit == "REML"){ ## precompute REML terms
        precompute$REML <- .precomputeREML(precompute = precompute, effects = effects2, param.vcov = nameparam.vcov, pair.vcov = design$vcov$pair.vcov, triplet.vcov = design$vcov$triplet.vcov, 
                                           logLik = logLik, score = score, information = information, vcov = vcov, df = df)
    }

    ## ** 2- GLS update of the mean parameters
    if(updateGLS){
        out$value <- value
        if(length(nameparam.mu)>0){
            if(trace>=1){cat("- GLS update \n")}
            value[nameparam.mu] <- .GLS(OmegaM1 = precompute$Omega$OmegaM1, precompute.XY = design$precompute.XY, precompute.XX = design$precompute.XX, precompute.weights = design$precompute.weights)[nameparam.mu]
        }
    }

    ## ** 3- compute partial derivatives regarding the mean
    if(trace>=1){cat("- residuals \n")}
    out$fitted <- design$mean %*% value[colnames(design$mean)]
    out$residuals <- design$Y - out$fitted

    if(is.null(design$weights$multiplier)){
        pattern.ntime <- stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name)
        pattern.cluster <- stats::setNames(design$vcov$Upattern$index.cluster, design$vcov$Upattern$name)
    }else{
        ## , design$precompute.weights$pattern)
        pattern.ntime <- stats::setNames(design$vcov$Upattern$n.time, design$vcov$Upattern$name)[design$precompute.weights$vcov]
        names(pattern.ntime) <- design$precompute.weights$pattern
        pattern.cluster <- tapply(1:NROW(design$vcov$pattern),factor(design$vcov$pattern$all, levels = design$precompute.weights$pattern), FUN = identity)
        names(pattern.cluster) <- design$precompute.weights$pattern
    }

    wRR <- cbind(out$residuals)
    if(!is.null(design$weights$XYR)){        
        wRR <- sweep(wRR, FUN = "*", MARGIN = 1, STATS = sqrt(design$weights$XYR))
    }
    precompute$RR <- .precomputeRR(residuals = wRR, pattern = design$precompute.weights$pattern, 
                                   pattern.ntime = pattern.ntime, pattern.cluster = pattern.cluster, index.cluster = design$index.cluster)                           

    if(score || information || vcov || df){

        wR <- out$residuals
        if(!is.null(design$weights$XYR)){        
            wR <- sweep(out$residuals, FUN = "*", MARGIN = 1, STATS = design$weights$XYR)
        }
        precompute$XR  <-  .precomputeXR(X = design$mean, residuals = wR, pattern = design$precompute.weights$pattern,
                                         pattern.ntime = pattern.ntime, pattern.cluster = pattern.cluster, index.cluster = design$index.cluster)
    }

    ## ** 4- compute likelihood derivatives

    ## *** log-likelihood
    if(logLik){
        if(trace>=1){cat("- log-likelihood \n")}
        out$logLik <- .logLik(X = design$mean, residuals = out$residuals, weights = design$weights, 
                              pattern = design$vcov$pattern, index.cluster = design$index.cluster, 
                              indiv = indiv, REML = method.fit=="REML", precompute = precompute)
    }

    ## *** score
    if(score || (vcov && robust)){ 
        if(trace>=1){cat("- score \n")}
        Mscore <- .score(X = design$mean, residuals = out$residuals, weights = design$weights, 
                         pattern = design$vcov$pattern, index.cluster = design$index.cluster, name.allcoef = name.allcoef,
                         indiv = indiv || (vcov && robust), REML = method.fit=="REML", effects = effects2, precompute = precompute)

        if(score){
            if(!indiv && (vcov && robust)){
                out$score <- colSums(Mscore)
            }else{
                out$score <- Mscore
            }
            attr(out$score,"message") <- attr(Mscore,"message")
        }
    }

    ## *** information
    if(information || vcov || (df & robust!=1)){
        if(trace>=1){cat("- information \n")}
        Minfo <- .information(X = design$mean, residuals = out$residuals, weights = design$weights, 
                              pattern = design$vcov$pattern, index.cluster = design$index.cluster, name.allcoef = name.allcoef, pair.vcov = design$vcov$pair.vcov,
                              indiv = indiv && information, REML = (method.fit=="REML"), type.information = type.information, effects = effects2, precompute = precompute)

        if(information){
            out$information <- Minfo
            attr(out$information, "type.information") <- type.information
            attr(out$information,"message") <- attr(Minfo,"message")            
        }
    }

    ## *** variance-covariance
    if(vcov || (df & robust!=1)){
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
            out$vcov <- Mvcov
            attr(out$vcov, "type.information") <- type.information
            attr(out$vcov, "robust") <- robust>0
            attr(out$vcov, "message") <- attr(Mvcov,"message")            
        }
    }else{
        Mvcov <- NULL
    }

    if(df){
        if(trace>=1){cat("- degrees-of-freedom \n")}
        if(is.null(attr(df,"method")) || attr(df,"method")=="analytic"){
            out$df <- .df_analytic(X = design$mean, residuals = out$residuals, weights = design$weights, vcov = out$vcov, robust = robust,
                                   pattern = design$vcov$pattern, index.cluster = design$index.cluster, name.allcoef = name.allcoef, pair.vcov = design$vcov$pair.vcov, triplet.vcov = design$vcov$triplet.vcov,
                                   REML = (method.fit=="REML"), type.information = type.information, effects = effects2, precompute = precompute, diag = TRUE)
        }else if(attr(df,"method")=="numeric"){
            out$df <- .df_numDeriv(param = value, design = design, method.fit = method.fit, type.information = type.information,
                                   transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                                   effects = effects2, vcov = Mvcov, robust = robust, method.numDeriv = method.numDeriv)
        }

        out$dVcov <- attr(out$df,"dVcov")
        attr(out$df,"dVcov") <- NULL        
    }

    ## ** 5- export
    return(out)
}

## * helper
## ** .optimGLS
## Implement GLS estimator i.e. \beta = (tX \OmegaM1 X)^{1} tX \OmegaM1 Y
.GLS <- function(OmegaM1, precompute.XY, precompute.XX, precompute.weights){

    param.mu <- colnames(precompute.XX$key)
    n.param <- length(param.mu)
    numerator <- matrix(0, nrow = n.param, ncol = 1)
    denominator <- matrix(0, nrow = n.param, ncol = n.param)

    for(iG in 1:NROW(precompute.weights)){ ## loop over covariance and replication patterns

        iPattern <- precompute.weights[iG,"pattern"]
        iVec.OmegaM1 <- as.double(OmegaM1[[iPattern]])
        numerator <- numerator + t(iVec.OmegaM1 %*%  precompute.XY[[iPattern]])
        denominator <- denominator + as.double(iVec.OmegaM1 %*%  precompute.XX$pattern[[iPattern]])[precompute.XX$key]

    }

    out <- solve(denominator, numerator)    
    return(stats::setNames(as.double(out), param.mu))
}

## ** .tripletVcov
.tripletVcov <- function(structure, X, index.cluster, index.clusterTime){
    
    n.pattern <- NROW(structure$Upattern)
    activeparam <- c(colnames(X),structure$param[is.na(structure$param$constraint),"name"])
    
    ## *** pairs of covariance parameters
    ## collect for each pattern the vector of parameters for a given cell in the residual variance-covariance matrix
    ls.Mparam.vcovvcov <- pairParamVcovVcov(structure)
    
    ## unique combinations of variance-covariance parameters across patterns
    Mparam.vcovvcov <- unique(do.call(rbind,ls.Mparam.vcovvcov))
    
    ## table with pairs for each pattern
    structure$pair.vcov <- structure$pair.vcov

    ## *** triplets of covariance parameters
    ## form all triplets in each pattern
    ls.triplet.vcovvcov <- apply(Mparam.vcovvcov, MARGIN = 1, FUN = function(iParam){
        unorderedTriplet(intersect(activeparam,iParam)) ## intersect will remove duplicates (e.g. sigma, one, sigma --> sigma, one) and keep a constant ordering
    }, simplify = FALSE)

    ## unique triplets across all patterns
    Utriplet.vcovvcov <- unique(t(do.call(cbind,ls.triplet.vcovvcov)))
    out <- data.frame(paste(Utriplet.vcovvcov[,1],Utriplet.vcovvcov[,2],Utriplet.vcovvcov[,3],sep = "_"),
                                    Utriplet.vcovvcov,
                                    matrix(NA, nrow = NROW(Utriplet.vcovvcov), ncol = n.pattern))
    names(out) <- c("name","param1","param2","param3",structure$Upattern$name)
    ## find which triplet appears in which pattern
    for(iPattern in structure$Upattern$name){ ## iP <- structure$Upattern$name[1]
        ## put parameter from the pattern first
        iParam.pattern <- structure$Upattern$param[[which(structure$Upattern$name == iPattern)]]
        iParam.active <- c(iParam.pattern,setdiff(activeparam,iParam.pattern))
        ## convert to numeric
        iUtriplet.vcov3.num <- matrix(as.numeric(factor(Utriplet.vcovvcov, levels = iParam.active)), nrow = NROW(Utriplet.vcovvcov), ncol = NCOL(Utriplet.vcovvcov))
        ## check all parameters belong to pattern
        out[,iPattern] <- rowSums(iUtriplet.vcov3.num <= length(iParam.pattern))==3
    }

    ## *** export
    return(out)
}
##----------------------------------------------------------------------
### moments.R ends here
