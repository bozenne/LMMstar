### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: okt  2 2021 (17:47) 
##           By: Brice Ozenne
##     Update #: 917
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * information.lmm (documentation)
##' @title Extract The Information From a Linear Mixed Model
##' @description Extract or compute the (expected) second derivative of the log-likelihood of a multivariate gaussian model.
##' @name information
##' 
##' @param x a \code{lmm} object.
##' @param data [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the information be output? Otherwise output the sum of all clusters of the derivatives.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the information. Only relevant if differs from the fitted values.
##' @param effects [character] Should the information relative to all coefficients be output (\code{"all"} or \code{"fixed"}),
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
information.lmm <- function(x, effects = NULL, data = NULL, p = NULL, indiv = FALSE, type.information = NULL,
                            transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    ## ** normalize user input
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(is.null(type.information)){
        type.information <- attr(x$information,"type.information")
        robust <- FALSE
    }else{
        type.information <- match.arg(type.information, c("expected","observed"))
        robust <- identical(attr(type.information,"robust"),TRUE)
    }
    if(is.null(effects)){
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform
    
    ## ** extract or recompute information
    if(is.null(data) && is.null(p) && (indiv == FALSE) && test.notransform && (robust==FALSE) && attr(x$information,"type.information")==type.information){
        keep.name <- stats::setNames(names(coef(x, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                                     names(coef(x, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))    

        design <- x$design ## useful in case of NA
        out <- x$information[keep.name,keep.name,drop=FALSE]
        if(transform.names){
            dimnames(out) <- list(names(keep.name),names(keep.name))
        }
    }else{
        test.precompute <- !is.null(x$design$precompute.XX) && !indiv
         
        if(!is.null(data)){
            ff.allvars <- c(all.vars(x$formula$mean), all.vars(x$formula$var))
            if(any(ff.allvars %in% names(data) == FALSE)){
                stop("Incorrect argument \'data\': missing variable(s) \"",paste(ff.allvars[ff.allvars %in% names(data) == FALSE], collapse = "\" \""),"\".\n")
            }

            design <- .model.matrix.lmm(formula.mean = x$formula$mean.design,
                                        structure = x$design$vcov,
                                        data = data,
                                        var.outcome = x$outcome$var,
                                        var.strata = x$strata$var, U.strata = x$strata$levels,
                                        var.time = x$time$var, U.time = x$time$levels,
                                        var.cluster = x$cluster$var,
                                        precompute.moments = test.precompute)
        }else{
            design <- x$design
        }

        if(!is.null(p)){
            if(any(duplicated(names(p)))){
                stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
            }
            if(any(names(x$param$type) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param$type)[names(x$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            p <- p[names(x$param$value)]
        }else{
            p <- x$param$value
        }
     
        out <- .moments.lmm(value = p, design = design, time = x$time, method.fit = x$method.fit, type.information = type.information,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                            logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects, robust = robust,
                            trace = FALSE, precompute.moments = test.precompute, transform.names = transform.names)$information
    }
    
    ## ** restaure NA
    if(length(x$index.na)>0 && indiv && is.null(data)){ 
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
.information <- function(X, residuals, precision, dOmega, d2Omega, Upattern.ncluster,
                         index.variance, time.variance, index.cluster, name.varcoef, name.allcoef,
                         pair.meanvarcoef, pair.varcoef, indiv, REML, type.information, effects, robust,
                         precompute){

    ## ** extract information
    test.loopIndiv <- indiv || is.null(precompute)
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    n.varcoef <- lapply(name.varcoef, length)
    n.allcoef <- length(name.allcoef)
    name.allvarcoef <- name.allcoef[name.allcoef %in% unique(unlist(name.varcoef))] ## make sure the ordering is correct
    n.allvarcoef <- length(name.allvarcoef)
    U.pattern <- names(dOmega)
    n.pattern <- length(U.pattern)
        
    npair.meanvarcoef <- lapply(pair.meanvarcoef, NCOL)
    npair.varcoef <- lapply(pair.varcoef, NCOL)

    ## ** prepare output
    name.effects <- attr(effects,"original.names")
    n.effects <- length(name.effects)
    if(test.loopIndiv && indiv){
        info <- array(0, dim = c(n.cluster, n.effects, n.effects),
                      dimnames = list(NULL, name.effects, name.effects))
    }else{
        info <- matrix(0, nrow = n.effects, ncol = n.effects,
                       dimnames = list(name.effects, name.effects)
                       )
    }    

    ## restrict to relevant parameters
    if(("variance" %in% effects == FALSE) && ("correlation" %in% effects == FALSE)){ ## compute hessian only for mean parameters
        test.vcov <- FALSE
        test.mean <- TRUE
    }else{
        if(REML && indiv){
            stop("Not possible to compute individual hessian for variance and/or correlation coefficients when using REML.\n")
        }
        if(("variance" %in% effects == FALSE) || ("correlation" %in% effects == FALSE)){ ## subset variance parameters
            name.varcoef <- stats::setNames(lapply(U.pattern,function(iPattern){intersect(name.effects,name.varcoef[[iPattern]])}),
                                     U.pattern)

            pair.meanvarcoef <- stats::setNames(lapply(U.pattern,function(iPattern){
                test.in <- (pair.meanvarcoef[[iPattern]][1,] %in% name.effects)+(pair.meanvarcoef[[iPattern]][2,] %in% name.effects)
                return(pair.meanvarcoef[[iPattern]][,test.in==2,drop=FALSE])
            }), U.pattern)

            pair.varcoef <- stats::setNames(lapply(U.pattern,function(iPattern){## iPattern <- 1
                test.in <- (pair.varcoef[[iPattern]][1,] %in% name.effects)+(pair.varcoef[[iPattern]][2,] %in% name.effects)
                iOut <- pair.varcoef[[iPattern]][,test.in==2,drop=FALSE]
                attr(iOut,"subset") <- which(test.in==2)
                attr(iOut,"key") <- matrix(NA, nrow = length(name.varcoef[[iPattern]]), ncol = length(name.varcoef[[iPattern]]), dimnames = list(name.varcoef[[iPattern]],name.varcoef[[iPattern]]))
                for(iCol in 1:NCOL(iOut)){
                    attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
                    attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol]] <- iCol
                }
                return(iOut)
            }), U.pattern)
            d2Omega <- stats::setNames(lapply(U.pattern,function(iPattern){
                
                return(d2Omega[[iPattern]][attr(pair.varcoef[[iPattern]],"subset")])
            }), U.pattern)

            n.varcoef <- lapply(name.varcoef, length)
            name.allvarcoef <- unique(unlist(name.varcoef))
            n.allvarcoef <- length(name.allvarcoef)
            npair.meanvarcoef <- lapply(pair.meanvarcoef, NCOL)
            npair.varcoef <- lapply(pair.varcoef, NCOL)
        }
        if("mean" %in% effects == FALSE){ ## compute hessian only for variance and/or correlation parameters
            if(REML && indiv){
                stop("Not possible to compute individual hessian for variance and/or correlation coefficients when using REML.\n")
            }

            test.vcov <- TRUE
            test.mean <- FALSE

        }else{ ## compute hessian all parameters
     
            test.vcov <- TRUE
            test.mean <- TRUE
        }
    }

    ## prepare REML term
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

            for(iPattern in 1:n.pattern){ ## iPattern <- 2
                iOmega <- precision[[iPattern]]
                idOmega <- dOmega[[iPattern]]

                OmegaM1_dOmega_OmegaM1[[iPattern]] <- stats::setNames(lapply(name.varcoef[[iPattern]], FUN = function(iVarcoef){iOmega %*% idOmega[[iVarcoef]] %*% iOmega}), name.varcoef[[iPattern]])

                ## loop over all pairs
                for(iPair in 1:npair.varcoef[[iPattern]]){ ## iPair <- 4
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
                for(iVarcoef in name.varcoef[[iPattern]]){ ## iVarcoef <- 1
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
    
        ## loop
        for (iPattern in U.pattern) { ## iPattern <- name.pattern[1]
            iOmega <- precision[[iPattern]]
            iTime <- NCOL(iOmega)
            iTime2 <- length(iOmega)
            iName.varcoef <- name.varcoef[[iPattern]]
            iN.varcoef <- length(iName.varcoef)

            iX <- matrix(unlist(precompute$XX$pattern[[iPattern]]), nrow = iTime2, ncol = dim(precompute$XX$pattern[[iPattern]])[3], byrow = FALSE)
                    
            ## **** mean,mean
            iValue <- (as.double(iOmega) %*% iX)[as.double(precompute$XX$key)]
            if(test.mean){
                info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + iValue
            }

            ## **** var,var
            if(test.vcov){

                ## precompute
                iMat <- tblock(t(do.call(rbind, dOmega[[iPattern]]) %*% iOmega))
                dOmega_OmegaM1 <- matrix(iMat,
                                         nrow = iTime2, ncol = iN.varcoef, dimnames = list(NULL,iName.varcoef), byrow = FALSE)
                tdOmega_OmegaM1 <- matrix(tblock(iMat),
                                          nrow = iTime2, ncol = iN.varcoef, dimnames = list(NULL,iName.varcoef), byrow = FALSE)
                iOmegaM1_dOmega_OmegaM1 <- matrix(iOmega %*% iMat,
                                                  nrow = iTime2, ncol = iN.varcoef, dimnames = list(NULL,iName.varcoef), byrow = FALSE)
                if(REML || type.information == "observed"){
                    iOmegaM1_d2Omega_OmegaM1 <- matrix(iOmega %*% tblock(t(do.call(rbind, d2Omega[[iPattern]]) %*% iOmega)),
                                                       nrow = iTime2, ncol = npair.varcoef[[iPattern]], byrow = FALSE)
                    iOmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1 <- do.call(cbind,lapply(1:npair.varcoef[[iPattern]], function(iPair){ ## iPair <- 4
                        iCoef1 <- pair.varcoef[[iPattern]][1,iPair]
                        iCoef2 <- pair.varcoef[[iPattern]][2,iPair]
                        out <- matrix(tdOmega_OmegaM1[,iCoef1], nrow = iTime, ncol = iTime) %*% matrix(iOmegaM1_dOmega_OmegaM1[,iCoef2], nrow = iTime, ncol = iTime) + matrix(tdOmega_OmegaM1[,iCoef2], nrow = iTime, ncol = iTime) %*% matrix(iOmegaM1_dOmega_OmegaM1[,iCoef1], nrow = iTime, ncol = iTime)
                        return(as.double(out))
                    }))
                    iOmegaM1_d2OmegaAndCo_OmegaM1 <- iOmegaM1_d2Omega_OmegaM1 - iOmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1
                    ## iOmega %*% d2Omega[[iPattern]][[iPair]] %*% iOmega - 2 * iOmega %*% dOmega[[iPattern]][[iCoef1]] %*% iOmega %*% dOmega[[iPattern]][[iCoef2]] %*% iOmega
                }
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
                    iX_OmegaM1_d2OmegaAndCo_OmegaM1_X <- t(iX) %*% iOmegaM1_d2OmegaAndCo_OmegaM1
                    for(iPair in 1:npair.varcoef[[iPattern]]){ ## iPair <- 1
                        iCoef1 <- pair.varcoef[[iPattern]][1,iPair]
                        iCoef2 <- pair.varcoef[[iPattern]][2,iPair]
                        REML.numerator2[,,REML.key[iCoef1,iCoef2]] <- REML.numerator2[,,REML.key[iCoef1,iCoef2]] + iX_OmegaM1_d2OmegaAndCo_OmegaM1_X[iDouble2Mat,iPair]
                    }
                    
                }

                ## compute contribution
                iTrace_O_dO_O_dO <- colSums(dOmega_OmegaM1[,pair.varcoef[[iPattern]][1,],drop=FALSE] * tdOmega_OmegaM1[,pair.varcoef[[iPattern]][2,],drop=FALSE])
                ## - 0.5 * tr(iOmega %*% dOmega[[iPattern]][[1]] %*% iOmega %*% dOmega[[iPattern]][[2]] - iOmega %*% d2Omega[[iPattern]][[iPair]])
                if(type.information == "expected"){
                    iValue <- 0.5 * Upattern.ncluster[iPattern] * iTrace_O_dO_O_dO
                }else if(type.information == "observed"){
                    id2Omega <- matrix(unlist(d2Omega[[iPattern]]), nrow = iTime2, ncol = npair.varcoef[[iPattern]])
                    iTrace_d2Omega <- colSums(sweep(id2Omega, MARGIN = 1, FUN = "*", STATS = as.double(precision[[iPattern]])))
                    iValue <- - 0.5 * Upattern.ncluster[iPattern] * (iTrace_O_dO_O_dO - iTrace_d2Omega) - 0.5 * as.double(as.double(precompute$RR[[iPattern]]) %*% iOmegaM1_d2OmegaAndCo_OmegaM1)
                    
                }
                
                ## store
                info[iName.varcoef,iName.varcoef]  <- info[iName.varcoef,iName.varcoef] + iValue[as.double(attr(pair.varcoef[[iPattern]],"key"))]
            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){

                ## compute
                iValue <- t(iOmegaM1_dOmega_OmegaM1) %*% matrix(precompute$XR[[iPattern]], nrow = iTime2, ncol = n.mucoef, dimnames = list(NULL,name.mucoef))

                ## store
                info[iName.varcoef,name.mucoef] <- info[iName.varcoef,name.mucoef] + iValue
                info[name.mucoef,iName.varcoef] <- info[name.mucoef,iName.varcoef] + t(iValue)
            }
        }

    }

    ## ** export
    if(REML && test.vcov){
        REML.denomM1 <- solve(REML.denom)
        REML.numerator2.bis <- array(NA, dim = dim(REML.numerator2), dimnames = dimnames(REML.numerator2))
        ls.REML.numerator1.denomM1 <- stats::setNames(lapply(1:dim(REML.numerator1)[3], FUN = function(iDim){REML.numerator1[,,iDim] %*% REML.denomM1}), dimnames(REML.numerator1)[[3]])
        ## ls.REML.numerator1.denomM1 <- apply(REML.numerator1, MARGIN = 3, FUN = function(iM){iM %*% REML.denomM1}, simplify = FALSE) ## only work on recent R versions
        for(iKey in 1:maxkey){ ## iKey <- 1
            iIndex <- which(REML.key == iKey, arr.ind = TRUE)
            REML.numerator2.bis[,,iKey] <-  ls.REML.numerator1.denomM1[[name.allvarcoef[iIndex[1,1]]]] %*% REML.numerator1[,,name.allvarcoef[iIndex[1,2]]]
        }

        REML.all <- as.double(REML.denomM1) %*% matrix(REML.numerator2.bis + REML.numerator2, nrow = prod(dim(REML.numerator2)[1:2]), ncol = dim(REML.numerator2)[3], byrow = FALSE)
        info[name.allvarcoef,name.allvarcoef] <- info[name.allvarcoef,name.allvarcoef] - 0.5 * REML.all[as.double(REML.key)]
        
    }

    if(robust){
        attr.info <- info
        attr.bread <- crossprod(.score(X = X, residuals = residuals, precision = precision, dOmega = dOmega, index.variance = index.variance, time.variance = time.variance, 
                                       index.cluster = index.cluster, name.varcoef = name.varcoef, name.allcoef = name.allcoef, indiv = TRUE, REML = REML, effects = effects,
                                       precompute = precompute) )
        info <- attr.info %*% solve(attr.bread) %*% attr.info
    }

    return(info)
}

##----------------------------------------------------------------------
### information.R ends here

