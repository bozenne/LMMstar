### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: May 20 2021 (12:17) 
##           By: Brice Ozenne
##     Update #: 376
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
##' @description Extract or compute the (expected) second derivative of the log-likelihood of a linear mixed model.
##' @name score
##' 
##' @param object a \code{lmm} object.
##' @param data [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the information be output? Otherwise output the sum of all clusters of the derivatives.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the information. Only relevant if differs from the fitted values.
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
information.lmm <- function(x, data = NULL, p = NULL, indiv = FALSE, type.information = NULL,
                            transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){
    options <- LMMstar.options()
    x.transform.sigma <- x$reparametrize$transform.sigma
    x.transform.k <- x$reparametrize$transform.k
    x.transform.rho <- x$reparametrize$transform.rho

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    if(is.null(type.information)){
        type.information <- options$type.information
        test.detail <- FALSE
    }else{
        test.detail <- identical(attr(type.information,"detail"),TRUE)
        type.information <- match.arg(type.information, c("expected","observed"))
    }

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, options = options,
                            x.transform.sigma = x.transform.sigma, x.transform.k = x.transform.k, x.transform.rho = x.transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    ## ** extract or recompute information
    if(is.null(data) && is.null(p) && (indiv == FALSE) && test.notransform){
        out <- x$information
    }else{
        REML <- x$method.fit == "REML"

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
                                        structure = x$structure)
            Y <- design$Y
            X <- design$X.mean
            index.vargroup <- design$index.vargroup
            index.cluster <- design$index.cluster
            index.time <- design$index.time
            X.var <- design$X.var
            pair.varcoef  <- design$param$pair.varcoef
            pair.meanvarcoef  <- design$param$pair.meanvarcoef
        }else{
            Y <- x$design$Y
            X <- x$design$X.mean
            index.vargroup <- x$design$index.vargroup
            index.cluster <- x$design$index.cluster
            index.time <- x$design$index.time
            X.var <- x$design$X.var
            pair.varcoef  <- x$design$param$pair.varcoef
            pair.meanvarcoef  <- x$design$param$pair.meanvarcoef
        }
        if(!is.null(p) || (test.notransform == FALSE)){
            if(!is.null(p)){
                if(any(duplicated(names(p)))){
                    stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
                }
                if(any(names(x$param$type) %in% names(p) == FALSE)){
                    stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param$type)[names(x$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
                }
            }else{
                p <- c(x$param$mu,x$param$sigma,x$param$k,x$param$rho)
            }

            beta <- p[names(x$param$mu)]
            name.pVar <- c(names(x$param$sigma),names(x$param$k),names(x$param$rho))

            reparametrize <- .reparametrize(p = p[name.pVar], type = x$param$type[name.pVar], strata = x$param$strata[name.pVar], time.levels = x$time$levels,
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

            Omega <- .calc_Omega(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$rho)], keep.interim = TRUE)
            precision <- lapply(Omega, solve)
            dOmega <- .calc_dOmega(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$rho)], Omega = Omega, Jacobian = reparametrize$Jacobian)
            if(REML || type.information == "observed"){
                d2Omega <- .calc_d2Omega(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$rho)],
                                         Omega = Omega, dOmega = dOmega, pair = pair.varcoef,
                                         Jacobian = reparametrize$Jacobian, dJacobian = reparametrize$dJacobian)
            }else{
                d2Omega <- NULL
            }
        }else{
            beta <- x$param$mu
            p <- c(beta,x$param$sigma,x$param$k,x$param$rho)
            name.pVar <- c(names(x$param$sigma),names(x$param$k),names(x$param$rho))

            reparametrize <- x$reparametrize

            precision <- x$OmegaM1
            dOmega <- x$dOmega
            d2Omega <- x$d2Omega
        }
        out <- .information(X = X, residuals = Y - X %*% beta, precision = precision, dOmega = dOmega, d2Omega = d2Omega,
                            index.variance = index.vargroup, time.variance = index.time, index.cluster = index.cluster, ## attr(X.var,"Upattern.index.time")
                            pair.meanvarcoef = pair.meanvarcoef, pair.varcoef = pair.varcoef, indiv = indiv, REML = REML, type.information = type.information)
        if(test.detail){
            attr(out,"detail") <- list(param = p, reparametrize = reparametrize, Y = Y, X.mean = X, X.var = X.var,
                                       index.variance = index.vargroup, time.variance = index.time, index.cluster = index.cluster,
                                       pair.meanvarcoef = pair.meanvarcoef, pair.varcoef = pair.varcoef,
                                       REML = REML,
                                       transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho)
        }

        if(transform.names && length(reparametrize$newname)>0){
            allnewname <- names(p)
            allnewname[match(name.pVar,colnames(out))] <- reparametrize$newname
            if(indiv){
                dimnames(out) <- list(NULL,allnewname,allnewname)
            }else{
                dimnames(out) <- list(allnewname,allnewname)
            }
        }
    }
    return(out)
}

## * .information
## REML term
## d 0.5 tr[(X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 X)] = 0.5 tr[ (X \OmegaM1 d'\Omega \OmegaM1 X) (X \OmegaM1 X)^{-2} (X \OmegaM1 d\Omega \OmegaM1 X) ]
##                                                                 - 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d'\Omega \OmegaM1 d\Omega \OmegaM1 X) + (X \OmegaM1 X)^{-1} (X \OmegaM1 d\Omega \OmegaM1 d'\Omega \OmegaM1 X) ]
##                                                                 + 0.5 tr[ (X \OmegaM1 X)^{-1} (X \OmegaM1 d2\Omega \OmegaM1 X) ]
.information <- function(X, residuals, precision, dOmega, d2Omega,
                         index.variance, time.variance, index.cluster,
                         pair.meanvarcoef, pair.varcoef, indiv, REML, type.information){

    if(indiv && REML){
        stop("Not possible to compute individual information when using REML.\n")
    }

    ## ** prepare
    ## coefficients
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- names(dOmega[[1]])
    n.varcoef <- length(name.varcoef)
    name.allcoef <- c(name.mucoef,name.varcoef)
    n.allcoef <- length(name.allcoef)
    var.pattern <- names(dOmega)

    npair.meanvarcoef <- NCOL(pair.meanvarcoef)
    npair.varcoef <- NCOL(pair.varcoef)

    ## store results
    if(indiv){
        info <- array(0, dim = c(n.cluster, n.allcoef, n.allcoef),
                     dimnames = list(NULL, name.allcoef, name.allcoef))
    }else{
        info <- matrix(0, nrow = n.allcoef, ncol = n.allcoef,
                      dimnames = list(name.allcoef, name.allcoef)
                      )
    }    

    ## precompute derivative for the variance
    dOmega.precomputed <- setNames(lapply(var.pattern, function(iPattern){
        iOut <- list(tr = vector(mode = "list", length = npair.varcoef))

        if(type.information == "observed"){
            iClusters <- which(var.pattern[index.variance]==iPattern)
            iIndex <- which(index.cluster %in% iClusters)
            iDF <- data.frame(id = index.cluster[iIndex],
                              res = residuals[iIndex,,drop=FALSE],
                              time =  time.variance[iIndex])
            iResidualsW <- reshape2::dcast(iDF, formula = id~time, value.var = "res")[,-1,drop=FALSE]

            iMisfit <- diag(1, NCOL(iResidualsW), NCOL(iResidualsW)) -  precision[[iPattern]] %*% crossprod(as.matrix(iResidualsW))/NROW(iResidualsW)
        }

        for(iPair in 1:npair.varcoef){
            term <- precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[1,iPair]]] %*% precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[2,iPair]]]
            if(type.information == "expected"){
                iOut$tr[[iPair]] <- 0.5 * tr(term)
            }else if(type.information == "observed"){
                ## iOut$tr[[iPair]] <- - 0.5 * tr(term) + tr(term %*% precision[[iPattern]] %*% crossprod(as.matrix(iResidualsW))/NROW(iResidualsW))
                iOut$tr[[iPair]] <- 0.5 * tr(term)
                if(!is.null(d2Omega[[iPattern]][[iPair]])){
                    iOut$tr[[iPair]] <- iOut$tr[[iPair]] - 0.5 * tr(precision[[iPattern]] %*% d2Omega[[iPattern]][[iPair]] %*% iMisfit)
                }
            }
            
        }

        if(REML || type.information == "observed"){
            iOut$cross <- setNames(lapply(1:n.varcoef, function(iCoef){            
                precision[[iPattern]] %*% dOmega[[iPattern]][[name.varcoef[iCoef]]] %*% precision[[iPattern]] 
            }), name.varcoef)
        }
        return(iOut)
    }), var.pattern)

    ## REML
    if(REML){
        REML.num <- lapply(1:npair.varcoef, function(i){
            list(numerator1a = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)),
                 numerator1b = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)),
                 numerator2 = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)),
                 numerator3 = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)))
        })
        REML.denom <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))

        REML.num.precomputed <- setNames(lapply(var.pattern, function(iPattern){
            iOut <- list(term2 = vector(mode = "list", length = npair.varcoef),
                         term3 = vector(mode = "list", length = npair.varcoef))                         

            for(iPair in 1:npair.varcoef){
                iOut$term2[[iPair]] <- precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[1,iPair]]] %*% precision[[iPattern]] %*% dOmega[[iPattern]][[pair.varcoef[2,iPair]]] %*% precision[[iPattern]]
                if(!is.null(d2Omega[[iPattern]][[iPair]])){
                    iOut$term3[[iPair]] <- precision[[iPattern]] %*% d2Omega[[iPattern]][[iPair]] %*% precision[[iPattern]]
                }
            }
            return(iOut)
        }), var.pattern)
    }

    ## ** compute information
    for(iId in 1:n.cluster){ ## iId <- 7
        iPattern <- index.variance[iId]
        iIndex <- attr(index.cluster,"sorted")[[iId]]
        ## iIndex <- which(index.cluster==iId)
        ## iIndex <- iIndex[order(time.variance[iIndex])] ## re-order observations according to the variance-covariance matrix

        iX <- X[iIndex,,drop=FALSE]
        tiX <- t(iX)
        iOmegaM1 <- precision[[iPattern]]

        ## *** mean,mean
        iValue <- tiX %*% iOmegaM1 %*% iX
        if(indiv){
            info[iId,name.mucoef,name.mucoef] <- iValue
        }else{
            ## info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + iValue
        }

        if(REML){
            REML.denom <- REML.denom + iValue
        }

        ## *** var,var
        for(iPair in 1:npair.varcoef){ ## iPair <- 1
            iCoef1 <- pair.varcoef[1,iPair]
            iCoef2 <- pair.varcoef[2,iPair]

            if(indiv){
                info[iId,iCoef1,iCoef2] <- dOmega.precomputed[[iPattern]]$tr[[iPair]]
                if(iCoef1 != iCoef2){
                    info[iId,iCoef2,iCoef1] <- dOmega.precomputed[[iPattern]]$tr[[iPair]]
                }
            }else{
                ## info[iCoef1,iCoef2] <- info[iCoef1,iCoef2] + dOmega.precomputed[[iPattern]]$tr[[iPair]]
                if(iCoef1 != iCoef2){
                    ## info[iCoef2,iCoef1] <- info[iCoef2,iCoef1] + dOmega.precomputed[[iPattern]]$tr[[iPair]]
                }
                if(REML){
                    REML.num[[iPair]]$numerator1a <- REML.num[[iPair]]$numerator1a + tiX %*% (dOmega.precomputed[[iPattern]]$cross[[iCoef1]]) %*% iX
                    REML.num[[iPair]]$numerator1b <- REML.num[[iPair]]$numerator1b + tiX %*% (dOmega.precomputed[[iPattern]]$cross[[iCoef2]]) %*% iX
                    REML.num[[iPair]]$numerator2 <- REML.num[[iPair]]$numerator2 + tiX %*% (REML.num.precomputed[[iPattern]]$term2[[iPair]]) %*% iX
                    if(!is.null(REML.num.precomputed[[iPattern]]$term3[[iPair]])){
                        REML.num[[iPair]]$numerator3 <- REML.num[[iPair]]$numerator3 + tiX %*% (REML.num.precomputed[[iPattern]]$term3[[iPair]]) %*% iX
                    }
                }
            }
        }

        ## *** mean,var
        if(type.information == "observed"){
            iResidual <- residuals[iIndex,,drop=FALSE]

            for(iPair in 1:npair.meanvarcoef){ ## iPair <- 1
                iCoef1 <- pair.meanvarcoef[1,iPair]
                iCoef2 <- pair.meanvarcoef[2,iPair]

                iValue <- tiX[iCoef1,,drop=FALSE] %*% dOmega.precomputed[[iPattern]]$cross[[iCoef2]] %*% iResidual 
                if(indiv){
                    info[iId,iCoef1,iCoef2] <- iValue
                    info[iId,iCoef2,iCoef1] <- iValue
                }else{
                    ## info[iCoef1,iCoef2] <- info[iCoef1,iCoef2] + iValue
                    ## info[iCoef2,iCoef1] <- info[iCoef2,iCoef1] + iValue
                }
            }

        }
    }

    ## ** export
    print(info)
    if(REML){
        REML.denom <- solve(REML.denom)
        for(iPair in 1:npair.varcoef){ ## iPair <- 1
            term1 <- tr(REML.denom %*% REML.num[[iPair]]$numerator1a %*% REML.denom %*% REML.num[[iPair]]$numerator1b)
            term2 <- tr(- 2 * REML.denom %*% REML.num[[iPair]]$numerator2)
            term3 <- tr(REML.denom %*% REML.num[[iPair]]$numerator3)

            iValue <- 0.5 * (term1 + term2 + term3) 
            info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] <- info[pair.varcoef[1,iPair],pair.varcoef[2,iPair]] + iValue
            if(pair.varcoef[1,iPair] != pair.varcoef[2,iPair]){
                info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] <- info[pair.varcoef[2,iPair],pair.varcoef[1,iPair]] + iValue
            }
        }
    }
    return(info)

}

##----------------------------------------------------------------------
### information.R ends here
