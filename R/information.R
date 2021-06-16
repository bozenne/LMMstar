### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: Jun 16 2021 (16:58) 
##           By: Brice Ozenne
##     Update #: 684
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * information.lmm (documentation)
##' @title Extract The Information From a Multivariate Gaussian Model
##' @description Extract or compute the (expected) second derivative of the log-likelihood of a multivariate gaussian model.
##' @name information
##' 
##' @param x a \code{lmm} object.
##' @param data [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the information be output? Otherwise output the sum of all clusters of the derivatives.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the information. Only relevant if differs from the fitted values.
##' @param effects [character] Should the information relative to all coefficients be output (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"}),
##' or only coefficients relative to the variance and correlation structure (\code{"variance"} or \code{"correlation"}).
##' @param type.information [character] Should the expected information be computed  (i.e. minus the expected second derivative) or the observed inforamtion (i.e. minus the second derivative).
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef} function.
##'
##' @return
##' When argument indiv is \code{FALSE}, a matrix with the value of the infroamtion relative to each pair of coefficient (in rows and columns) and each cluster (in rows).
##' When argument indiv is \code{TRUE}, a 3-dimensional array with the value of the information relative to each pair of coefficient (dimension 2 and 3) and each cluster (dimension 1).
##' 

## * information.lmm (code)
##' @rdname information
##' @export
information.lmm <- function(x, effects = "all", data = NULL, p = NULL, indiv = FALSE, type.information = NULL,
                            transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(is.null(type.information)){
        type.information <- attr(x$information,"type.information")
        test.detail <- FALSE
        robust <- FALSE
    }else{
        test.detail <- identical(attr(type.information,"detail"),TRUE)
        robust <- identical(attr(type.information,"robust"),TRUE)
        type.information <- match.arg(type.information, c("expected","observed"))
    }
    if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","variance","correlation"), several.ok = TRUE)

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    ## ** extract or recompute information
    if(is.null(data) && is.null(p) && (indiv == FALSE) && test.notransform && (robust==FALSE)){
        design <- x$design ## useful in case of NA
        out <- x$information
        if(transform.names){
            colnames(out)[match(names(x$reparametrize$p),colnames(out))] <- x$reparametrize$newname
            rownames(out)[match(names(x$reparametrize$p),rownames(out))] <- x$reparametrize$newname
        }
        if("mean" %in% effects == FALSE){
            out <- out[x$param$type!="mu",x$param$type!="mu",drop=FALSE]
        }else if("variance" %in% effects == FALSE && "correlation" %in% effects == FALSE){
            out <- out[x$param$type=="mu",x$param$type=="mu",drop=FALSE]
        }
        if(test.detail){
            param <- x$param
            param$value[names(x$reparametrize$p)] <- x$reparametrize$p
            attr(out,"detail") <- list(param = param, reparametrize = x$reparametrize, Y = x$design$Y, X.mean = x$design$X.mean, X.var = x$design$X.var,
                                       index.variance = x$design$X.var$cluster, time.variance = x$design$index.time, index.cluster = x$design$index.cluster, name.varcoef = x$design$X.var$param,
                                       pair.meanvarcoef = x$design$param$pair.meanvarcoef, pair.varcoef = x$design$param$pair.varcoef,
                                       REML = (x$method.fit=="REML"), precompute = list(XX = x$design$precompute.XX),
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        }
    }else{
        REML <- x$method.fit == "REML"
        test.precompute <- !is.null(x$design$precompute.XX) && !indiv
         
        if(!is.null(data)){
            ff.allvars <- c(all.vars(x$formula$mean), all.vars(x$formula$var))
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            design <- .model.matrix.lmm(formula.mean = x$formula$mean.design,
                                        formula.var = x$formula$var.design,
                                        data = data,
                                        var.outcome = x$outcome$var,
                                        var.strata = x$strata$var, U.strata = x$strata$levels,
                                        var.time = x$time$var, U.time = x$time$levels,
                                        var.cluster = x$cluster$var,
                                        structure = x$structure,
                                        precompute.moments = test.precompute)
        }else{
            design <- x$design
        }

        Y <- design$Y
        X <- design$X.mean
        index.vargroup <- design$X.var$cluster
        index.cluster <- design$index.cluster
        index.time <- design$index.time
        X.var <- design$X.var
        name.varcoef  <- design$X.var$param
        pair.varcoef  <- design$param$pair.varcoef
        pair.meanvarcoef  <- design$param$pair.meanvarcoef
        name.allcoef <- names(x$param$value)
        index.var <- x$param$type %in% c("sigma","k","rho")
        precompute <- design$precompute.XX

        if(!is.null(p) || (test.notransform == FALSE)){
            if(!is.null(p)){
                if(any(duplicated(names(p)))){
                    stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
                }
                if(any(names(x$param$type) %in% names(p) == FALSE)){
                    stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param$type)[names(x$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
                }
            }else{
                p <- x$param$value
            }

            reparametrize <- .reparametrize(p = p[index.var], type = x$param$type[index.var], strata = x$param$strata[index.var], time.levels = x$time$levels,
                                            time.k = x$design$param$time.k, time.rho = x$design$param$time.rho,
                                            Jacobian = TRUE, dJacobian = 2*(REML || type.information == "observed"), inverse = FALSE,
                                            transform.sigma = transform.sigma,
                                            transform.k = transform.k,
                                            transform.rho = transform.rho,
                                            transform.names = TRUE)

            if(reparametrize$transform==FALSE){
                reparametrize$newname <- NULL
                reparametrize$Jacobian <- NULL
                reparametrize$dJacobian <- NULL
            }

            Omega <- .calc_Omega(object = X.var, param = p, keep.interim = TRUE)
            precision <- lapply(Omega, solve)
            dOmega <- .calc_dOmega(object = X.var, param = p, type = x$param$type, Omega = Omega, Jacobian = reparametrize$Jacobian)
            if(REML || type.information == "observed"){
                d2Omega <- .calc_d2Omega(object = X.var, param = p, type = x$param$type,
                                         Omega = Omega, dOmega = dOmega, pair = pair.varcoef,
                                         Jacobian = reparametrize$Jacobian, dJacobian = reparametrize$dJacobian)
            }else{
                d2Omega <- NULL
            }
        }else{
            p <- x$param$value
            reparametrize <- x$reparametrize
            precision <- x$OmegaM1
            dOmega <- x$dOmega
            d2Omega <- x$d2Omega
        }

        out <- .information(X = X, residuals = Y - X %*% p[x$param$type=="mu"], precision = precision, dOmega = dOmega, d2Omega = d2Omega, robust = robust,
                            index.variance = index.vargroup, time.variance = index.time, index.cluster = index.cluster, name.varcoef = name.varcoef, name.allcoef = name.allcoef,
                            pair.meanvarcoef = pair.meanvarcoef, pair.varcoef = pair.varcoef, indiv = indiv, REML = REML, type.information = type.information,
                            effects = effects, precompute = precompute, X.var = X.var)

        if(test.detail){
            param <- x$param
            param$value <- p[names(param$type)]
            attr(out,"detail") <- list(param = param, reparametrize = reparametrize, Y = Y, X.mean = X, X.var = X.var,
                                       index.variance = index.vargroup, time.variance = index.time, index.cluster = index.cluster, name.varcoef = name.varcoef,
                                       pair.meanvarcoef = pair.meanvarcoef, pair.varcoef = pair.varcoef,
                                       REML = REML,
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        }

        if(transform.names && length(reparametrize$newname)>0){
            allnewname <- stats::setNames(name.allcoef,name.allcoef)
            allnewname[index.var] <- reparametrize$newname
            names(allnewname)[match(names(reparametrize$p), names(allnewname))] <- names(reparametrize$p)
            if(indiv){
                dimnames(out) <- list(NULL,allnewname[dimnames(out)[[2]]],allnewname[dimnames(out)[[3]]])
            }else{
                dimnames(out) <- list(allnewname[dimnames(out)[[1]]],allnewname[dimnames(out)[[2]]])
            }
        }
    }
    
    ## ** restaure NA
    if(length(x$index.na)>0 && indiv){ 
        iAdd <- .addNA(index.na = x$index.na, design = design, time = x$time)
        if(length(iAdd$missing.cluster)>0){
            out.save <- out
            out <- array(NA, dim = c(iAdd$n.allcluster, dim(out.save)[2:3]),
                         dimnames = c(list(NULL), dimnames(out.save)[2:3]))
            out[match(design$cluster$levels, iAdd$allcluster),,] <- out.save
        }
    }

    ## ** export
    return(out)
}

## * .information
## REML term
## d 0.5 tr[(X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 X)] = 0.5 tr[ (X \OmegaM1 d'\Omega \OmegaM1 X) (X \OmegaM1 X)^{-2} (X \OmegaM1 d\Omega \OmegaM1 X) ]
##                                                                 - 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d'\Omega \OmegaM1 d\Omega \OmegaM1 X) + (X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 d'\Omega \OmegaM1 X) ]
##                                                                 + 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d2\Omega \OmegaM1 X) ]
.information <- function(X, residuals, precision, dOmega, d2Omega, robust,
                         index.variance, time.variance, index.cluster, name.varcoef, name.allcoef,
                         pair.meanvarcoef, pair.varcoef, indiv, REML, type.information, effects, precompute, X.var){

    ## ** extract information
    test.loopIndiv <- indiv || is.null(precompute)
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    n.varcoef <- lapply(name.varcoef, length)
    n.allcoef <- length(name.allcoef)
    name.allvarcoef <- unique(unlist(name.varcoef))
    n.allvarcoef <- length(name.allvarcoef)
    U.pattern <- names(dOmega)
    n.pattern <- length(U.pattern)

    npair.meanvarcoef <- lapply(pair.meanvarcoef, NCOL)
    npair.varcoef <- lapply(pair.varcoef, NCOL)

    if(!is.null(precompute) && "RR" %in% names(precompute) == FALSE){
        precompute$RR <-  .precomputeRR(residuals = residuals, pattern = X.var$pattern,
                                        pattern.time = X.var$index.time, pattern.cluster = attr(X.var$cluster, "index.byPattern"), index.cluster = attr(index.cluster,"sorted"))
    }
    if(!is.null(precompute) && "XR" %in% names(precompute) == FALSE){
        precompute$XR <-  .precomputeRR(X = X, residuals = residuals, pattern = X.var$pattern,
                                        pattern.time = X.var$index.time, pattern.cluster = attr(X.var$cluster, "index.byPattern"), index.cluster = attr(index.cluster,"sorted"))
    }
    
    ## ** prepare output
    if("mean" %in% effects == FALSE){ ## compute information only for variance - correlation parameters
        if(test.loopIndiv && indiv){
            info <- array(0, dim = c(n.cluster, n.allvarcoef, n.allvarcoef),
                          dimnames = list(NULL, name.allvarcoef, name.allvarcoef))
        }else{
            info <- matrix(0, nrow = n.allvarcoef, ncol = n.allvarcoef,
                           dimnames = list(name.allvarcoef, name.allvarcoef)
                           )
        }    
        test.vcov <- TRUE
        test.mean <- FALSE
        if(REML && indiv){
            stop("Not possible to compute individual hessian for variance and/or correlation coefficients when using REML.\n")
        }
    }else if("variance" %in% effects == FALSE && "correlation" %in% effects == FALSE){ ## compute information only for mean parameters
        if(test.loopIndiv && indiv){
            info <- array(0, dim = c(n.cluster, n.mucoef, n.mucoef),
                          dimnames = list(NULL, name.mucoef, name.mucoef))
        }else{
            info <- matrix(0, nrow = n.mucoef, ncol = n.mucoef,
                           dimnames = list(name.mucoef, name.mucoef)
                           )
        }    
        test.vcov <- FALSE
        test.mean <- TRUE
    }else{ ## compute information only for all parameters
        if(test.loopIndiv && indiv){
            info <- array(0, dim = c(n.cluster, n.allcoef, n.allcoef),
                          dimnames = list(NULL, name.allcoef, name.allcoef))
        }else{
            info <- matrix(0, nrow = n.allcoef, ncol = n.allcoef,
                           dimnames = list(name.allcoef, name.allcoef)
                           )
        }    
        test.vcov <- TRUE
        test.mean <- TRUE
        if(REML && indiv){
            stop("Not possible to compute individual hessian for variance and/or correlation coefficients when using REML.\n")
        }
    }

    if(test.vcov && REML){
        REML.key <- matrix(NA, nrow = n.allvarcoef, ncol = n.allvarcoef, dimnames = list(name.allvarcoef, name.allvarcoef))
        maxkey <- sum(lower.tri(REML.key, diag = TRUE))
        REML.key[lower.tri(REML.key, diag = TRUE)] <- 1:maxkey
        REML.key[upper.tri(REML.key)] <- t(REML.key)[upper.tri(REML.key)]

        REML.denom <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
        ## REML.numerator1 <- stats::setNames(lapply(U.pattern, function(iPattern) { array(0, dim = c(n.mucoef, n.mucoef, n.varcoef[[iPattern]]), dimnames = list(name.mucoef, name.mucoef, name.varcoef[[iPattern]])) }), U.pattern)
        REML.numerator1 <- array(0, dim = c(n.mucoef, n.mucoef, n.allvarcoef), dimnames = list(name.mucoef, name.mucoef, name.allvarcoef))
        REML.numerator2 <- array(0, dim = c(n.mucoef, n.mucoef, maxkey), dimnames = list(name.mucoef, name.mucoef, NULL))
    }     
    
    ## ** compute information
    ## *** looping over individuals
    if(test.loopIndiv){
        
        ## ** precompute 
        if(test.vcov){
            OmegaM1_dOmega_OmegaM1 <- stats::setNames(vector(mode = "list", length = n.pattern), U.pattern)
            tr_OmegaM1_d2OmegaAndCo <- stats::setNames(lapply(1:n.pattern, function(iPattern){rep(NA, npair.varcoef[[iPattern]])}), U.pattern)
            if(REML || type.information == "observed"){
                OmegaM1_d2OmegaAndCo_OmegaM1 <- stats::setNames(lapply(1:n.pattern, function(iPattern){array(NA, dim = c(NCOL(precision[[iPattern]]),NCOL(precision[[iPattern]]), npair.varcoef[[iPattern]]))}), U.pattern)
            }

            for(iPattern in 1:n.pattern){ ## iPattern <- 3
                iOmega <- precision[[iPattern]]
                idOmega <- dOmega[[iPattern]]

                OmegaM1_dOmega_OmegaM1[[iPattern]] <- setNames(lapply(name.varcoef[[iPattern]], FUN = function(iVarcoef){iOmega %*% idOmega[[iVarcoef]] %*% iOmega}), name.varcoef[[iPattern]])

                ## loop over all pairs
                for(iPair in 1:npair.varcoef[[iPattern]]){ ## iPair <- 2
                    iCoef1 <- pair.varcoef[[iPattern]][1,iPair]
                    iCoef2 <- pair.varcoef[[iPattern]][2,iPair]

                    iTerm21 <- OmegaM1_dOmega_OmegaM1[[iPattern]][[iCoef2]] %*% idOmega[[iCoef1]]
                    
                    ## trace
                    if(type.information == "expected"){
                        tr_OmegaM1_d2OmegaAndCo[[iPattern]][iPair] <- tr(iTerm21)
                    }else if(type.information == "observed"){
                        tr_OmegaM1_d2OmegaAndCo[[iPattern]][iPair] <- - tr(iTerm21 - iOmega %*% d2Omega[[iPattern]][[iPair]])
                    
                    }

                    if(REML || type.information == "observed"){
                        iTerm12 <- OmegaM1_dOmega_OmegaM1[[iPattern]][[iCoef1]] %*% idOmega[[iCoef2]]
                        OmegaM1_d2OmegaAndCo_OmegaM1[[iPattern]][,,iPair] <- iOmega %*% d2Omega[[iPattern]][[iPair]] %*% iOmega - (iTerm12 + iTerm21) %*% iOmega
                    }
                    
                }
            }
        }
            
        ## loop
        for(iId in 1:n.cluster){ ## iId <- 7
            iPattern <- index.variance[iId]
            iIndex <- attr(index.cluster,"sorted")[[iId]]
            ## iIndex <- which(index.cluster==iId)
            ## iIndex <- iIndex[order(time.variance[iIndex])] ## re-order observations according to the variance-covariance matrix

            iX <- X[iIndex,,drop=FALSE]
            tiX <- t(iX)
            iOmega <- precision[[iPattern]]
            if(type.information == "observed"){
                iResidual <- residuals[iIndex,,drop=FALSE]
            }
        
            ## **** mean,mean
            iValue <- tiX %*% iOmega %*% iX
            if(test.mean){
                if(indiv){
                    info[iId,name.mucoef,name.mucoef] <- iValue
                }else{
                    info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + iValue
                }
            }
            if(REML && test.vcov){
                REML.denom <- REML.denom + iValue
                for(iVarcoef in 1:n.varcoef[[iPattern]]){ ## iVarcoef <- 1
                    REML.numerator1[,,iVarcoef] <- REML.numerator1[,,iVarcoef] + tiX %*% OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iX
                }
            }

            ## **** var,var
            if(test.vcov){

                for(iPair in 1:npair.varcoef[[iPattern]]){ ## iPair <- 1
                    iCoef1 <- pair.varcoef[[iPattern]][1,iPair]
                    iCoef2 <- pair.varcoef[[iPattern]][2,iPair]

                    iValue <- 0.5 * tr_OmegaM1_d2OmegaAndCo[[iPattern]][iPair]
                    ## 0.5 * tr(iOmega %*% idOmega$sigma %*% iOmega %*% idOmega$sigma)
                    if(type.information == "observed"){
                        iValue <- iValue - 0.5 * t(iResidual) %*% OmegaM1_d2OmegaAndCo_OmegaM1[[iPattern]][,,iPair] %*% iResidual
                    }
                    
                    if(indiv){
                        info[iId,iCoef1,iCoef2] <- iValue
                        if(iCoef1 != iCoef2){
                            info[iId,iCoef2,iCoef1] <- iValue
                        }
                    }else{
                        info[iCoef1,iCoef2] <- info[iCoef1,iCoef2] + iValue
                        if(iCoef1 != iCoef2){
                            info[iCoef2,iCoef1] <- info[iCoef2,iCoef1] + iValue
                        }
                    }

                    if(REML){
                        iKey <- REML.key[iCoef1,iCoef2]
                        REML.numerator2[,,iKey] <- REML.numerator2[,,iKey] + tiX %*% OmegaM1_d2OmegaAndCo_OmegaM1[[iPattern]][,,iPair] %*% iX
                    }
                }
            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){

                for(iPair in 1:npair.meanvarcoef[[iPattern]]){ ## iPair <- 1
                    iCoef1 <- pair.meanvarcoef[[iPattern]][1,iPair]
                    iCoef2 <- pair.meanvarcoef[[iPattern]][2,iPair]

                    iValue <- tiX[iCoef1,,drop=FALSE] %*% OmegaM1_dOmega_OmegaM1[[iPattern]][[iCoef2]] %*% iResidual 

                    if(indiv){
                        info[iId,iCoef1,iCoef2] <- iValue
                        info[iId,iCoef2,iCoef1] <- iValue
                    }else{
                        info[iCoef1,iCoef2] <- info[iCoef1,iCoef2] + iValue
                        info[iCoef2,iCoef1] <- info[iCoef2,iCoef1] + iValue
                    }
                }

            }
        }
    }

    ## *** looping over covariance patterns
    if(!test.loopIndiv){
        ## precompute
        ncluster.pattern <- sapply(attr(index.variance,"index.byPattern"),length)
        name.pattern <- names(ncluster.pattern)
    
        ## loop
        for (iPattern in name.pattern) { ## iPattern <- name.pattern[1]
            iOmega <- precision[[iPattern]]
            iTime <- NCOL(iOmega)
            iTime2 <- length(iOmega)
            iName.varcoef <- name.varcoef[[iPattern]]
            iN.varcoef <- length(iName.varcoef)

            iX <- matrix(unlist(precompute$XX$pattern[[iPattern]]), nrow = iTime2, ncol = dim(precompute$XX$pattern[[iPattern]])[3], byrow = FALSE)
                    
            ## **** mean,mean
            iValue <- (as.double(iOmega) %*% iX)[precompute$XX$key]
            if(test.mean){
                info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + iValue
            }

            ## **** var,var
            if(test.vcov){

                ## precompute
                iMat <- tblock(t(do.call(rbind, dOmega[[iPattern]]) %*% iOmega))
                dOmega_OmegaM1 <- matrix(iMat,
                                         nrow = iTime2, ncol = iN.varcoef, dimnames = list(NULL,iName.varcoef), byrow = FALSE)
                iOmegaM1_dOmega_OmegaM1 <- matrix(iOmega %*% iMat,
                                                  nrow = iTime2, ncol = iN.varcoef, dimnames = list(NULL,iName.varcoef), byrow = FALSE)
                iOmegaM1_d2Omega_OmegaM1 <- matrix(iOmega %*% tblock(t(do.call(rbind, d2Omega[[iPattern]]) %*% iOmega)),
                                                   nrow = iTime2, ncol = npair.varcoef[[iPattern]], byrow = FALSE)
                iOmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1 <- do.call(cbind,lapply(1:npair.varcoef[[iPattern]], function(iPair){ ## iPair <- 2
                    iCoef1 <- pair.varcoef[[iPattern]][1,iPair]
                    iCoef2 <- pair.varcoef[[iPattern]][2,iPair]
                    out <- matrix(dOmega_OmegaM1[,iCoef1], nrow = iTime, ncol = iTime) %*% matrix(iOmegaM1_dOmega_OmegaM1[,iCoef2], nrow = iTime, ncol = iTime)
                    return(as.double(out))
                }))

                if(REML){
                    iDouble2Mat <- as.vector(precompute$XX$key)
                    ## denominator
                    REML.denom <- REML.denom + (as.double(iOmega) %*% iX)[iDouble2Mat]
                    ## numerator 1
                    iX_OmegaM1_dOmega_OmegaM1_X <- t(iX) %*% iOmegaM1_dOmega_OmegaM1
                    for(iVarcoef in iName.varcoef){ ## iVarcoef <- iName.varcoef[1]
                        REML.numerator1[,,iVarcoef] <- REML.numerator1[,,iVarcoef] + iX_OmegaM1_dOmega_OmegaM1_X[iDouble2Mat,iVarcoef]
                    }
                    ## numerator 2
                    iX_OmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1_X <- t(iX) %*% iOmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1
                    for(iPair in 1:npair.varcoef[[iPattern]]){ ## iPair <- 1
                        REML.numerator2[,,iPair] <- REML.numerator2[,,iPair] + iX_OmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1_X[iDouble2Mat,iPair]
                    }
                }

                ## compute contribution
                id2Omega <- matrix(unlist(d2Omega[[iPattern]]), nrow = iTime2, ncol = npair.varcoef[[iPattern]])
                iTrace_d2Omega <- colSums(sweep(id2Omega, MARGIN = 1, FUN = "*", STATS = as.double(precision[[iPattern]])))
                iTrace_O_dO_O_dO <- sapply(1:npair.varcoef[[iPattern]], function(iPair){ ## iPair <- 2
                    return(crossprod(dOmega_OmegaM1[,pair.varcoef[[iPattern]][1,iPair]], dOmega_OmegaM1[,pair.varcoef[[iPattern]][2,iPair]]))
                })
                ## - 0.5 * tr(iOmega %*% dOmega[[iPattern]][[1]] %*% iOmega %*% dOmega[[iPattern]][[2]] - iOmega %*% d2Omega[[iPattern]][[iPair]])
                if(type.information == "expected"){
                    iValue <- 0.5 * iTrace_O_dO_O_dO
                }else if(type.information == "observed"){
                    iValue <- - 0.5 * iTrace_O_dO_O_dO + 0.5 * iTrace_d2Omega - 0.5 * as.double(as.double(precompute$RR[[iPattern]]) %*% iOmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1)
                }

                ## store
                info[iName.varcoef,iName.varcoef]  <- iValue[attr(pair.varcoef[[iPattern]],"key")]
            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){

                ## compute
                iValue <- t(iOmegaM1_dOmega_OmegaM1) %*% matrix(precompute$XR[[iPattern]], nrow = iTime2, ncol = n.mucoef, dimnames = list(NULL,name.mucoef))

                ## store
                info[iName.varcoef,name.mucoef] <- iValue
                info[name.mucoef,iName.varcoef] <- t(iValue)
            }
        }

    }

    ## ** export
    if(REML && test.vcov){
        REML.denomM1 <- solve(REML.denom)
        REML.numerator1.bis <- 0*REML.numerator2
        for(iKey in 1:maxkey){ ## iKey <- 1
            iIndex <- which(REML.key == iKey, arr.ind = TRUE)
            iCoef1 <- name.allvarcoef[iIndex[1,1]]
            iCoef2 <- name.allvarcoef[iIndex[1,2]]
            REML.numerator1.bis[,,iKey] <-  REML.numerator1[,,iCoef1] %*% REML.denomM1 %*% REML.numerator1[,,iCoef2]
        }
        REML.all <- as.double(REML.denomM1) %*% matrix(REML.numerator1.bis + REML.numerator2, nrow = prod(dim(REML.numerator2)[1:2]), ncol = dim(REML.numerator2)[3], byrow = FALSE)
        info[name.allvarcoef,name.allvarcoef] <- info[name.allvarcoef,name.allvarcoef] - 0.5 * REML.all[as.double(REML.key)]
        
    }

    if(robust){
        attr.info <- info
        attr.bread <- crossprod(.score(X = X, residuals = residuals, precision = precision, dOmega = dOmega, index.variance = index.variance, time.variance = time.variance, 
                                       index.cluster = index.cluster, name.varcoef = name.varcoef, name.allcoef = name.allcoef, indiv = TRUE, REML = REML, effects = effects, precompute = precompute) )
        info <- attr.info %*% solve(attr.bread) %*% attr.info
    }
    return(info)

}

##----------------------------------------------------------------------
### information.R ends here
