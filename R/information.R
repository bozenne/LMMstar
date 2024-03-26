### information.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2021 (22:13) 
## Version: 
## Last-Updated: Mar 26 2024 (09:39) 
##           By: Brice Ozenne
##     Update #: 1128
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
##' 
##' @param x a \code{lmm} object.
##' @param newdata [data.frame] dataset relative to which the information should be computed. Only relevant if differs from the dataset used to fit the model.
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
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef.lmm} function.
##'
##' @return
##' When argument indiv is \code{FALSE}, a matrix with the value of the infroamtion relative to each pair of coefficient (in rows and columns) and each cluster (in rows).
##' When argument indiv is \code{TRUE}, a 3-dimensional array with the value of the information relative to each pair of coefficient (dimension 2 and 3) and each cluster (dimension 1).
##'
##' @keywords methods

## * information.lmm (code)
##' @export
information.lmm <- function(x, effects = NULL, newdata = NULL, p = NULL, indiv = FALSE, type.information = NULL,
                            transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    ## ** normalize user input
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(is.null(type.information)){
        robust <- FALSE
        type.information <- x$args$type.information
    }else{
        robust <- identical(attr(type.information,"robust"),TRUE)
        type.information <- match.arg(type.information, c("expected","observed"))
    }
    if(is.null(effects)){
        effects <- options$effects
    }else if(identical(effects,"all")){
        effects <- c("mean","variance","correlation")
    }
    effects <- match.arg(effects, c("mean","fixed","variance","correlation"), several.ok = TRUE)
    effects[effects== "fixed"] <- "mean"

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
    if(is.null(newdata) && is.null(p) && (indiv == FALSE) && test.notransform && (robust==FALSE) && x$args$type.information==type.information){
        keep.name <- stats::setNames(names(coef(x, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                                     names(coef(x, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))    

        design <- x$design ## useful in case of NA
        out <- x$information[keep.name,keep.name,drop=FALSE]
        if(transform.names){
            dimnames(out) <- list(names(keep.name),names(keep.name))
        }
    }else{
        test.precompute <- !is.null(x$design$precompute.XX) && !indiv
         
        if(!is.null(newdata)){
            design <- stats::model.matrix(x, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- x$design
        }

        out <- .moments.lmm(value = theta, design = design, time = x$time, method.fit = x$args$method.fit, type.information = type.information,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                            logLik = FALSE, score = FALSE, information = TRUE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects, robust = robust,
                            trace = FALSE, precompute.moments = test.precompute, transform.names = transform.names)$information
    }

    ## ** restaure NAs and name
    if(indiv){

        if(!is.numeric(x$cluster$levels)){
            dimnames(out)[[1]] <- x$cluster$levels[match(1:dim(out)[[1]],x$cluster$index)]

        } 
        out <- restaureNA(out, index.na = x$index.na,
                          level = "cluster", cluster = x$cluster)
        
    }


    ## ** re-order values when converting to sd with strata (avoid sd0:0 sd0:1 sd1:0 sd1:1 sd2:0 sd2:1 ...)
    if("variance" %in% effects && transform.k %in% c("sd","var","logsd","logvar") && x$strata$n>1 && transform.names){
        out.name <- names(stats::coef(x, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = TRUE))
        if(indiv){
            out <- out[,out.name,out.name,drop=FALSE]
        }else{
            out <- out[out.name,out.name,drop=FALSE]
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
.information <- function(X, residuals, precision, dOmega, d2Omega,
                         Upattern.ncluster, weights, scale.Omega,
                         pattern, index.cluster, name.allcoef,
                         pair.meanvcov, pair.vcov, indiv, REML, type.information, effects, robust,
                         precompute){

    ## ** extract information
    test.loopIndiv <- indiv || is.null(precompute)
    n.obs <- length(index.cluster)
    n.cluster <- length(pattern)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- lapply(dOmega, names)
    n.varcoef <- lengths(name.varcoef)
    n.allcoef <- length(name.allcoef)
    name.allvarcoef <- name.allcoef[name.allcoef %in% unique(unlist(name.varcoef))] ## make sure the ordering is correct
    n.allvarcoef <- length(name.allvarcoef)
    U.pattern <- names(Upattern.ncluster)
    n.pattern <- length(U.pattern)
    npair.meanvcov <- lapply(pair.meanvcov, NCOL)
    npair.vcov <- lapply(pair.vcov, NCOL)

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
    if(any(is.na(attr(precision,"logdet")))){ ## non positive definite residual variance covariance
        return(info*NA)
    }

    ## restrict to relevant parameters
    if(("variance" %in% effects == FALSE) && ("correlation" %in% effects == FALSE)){ ## compute hessian only for mean parameters
        test.vcov <- FALSE
        test.mean <- n.mucoef>0
    }else{
        if(REML && indiv){
            stop("Not possible to compute individual hessian for variance and/or correlation coefficients when using REML.\n")
        }
        if(("variance" %in% effects == FALSE) || ("correlation" %in% effects == FALSE)){ ## subset variance parameters
            name.varcoef <- stats::setNames(lapply(U.pattern,function(iPattern){intersect(name.effects,name.varcoef[[iPattern]])}),
                                     U.pattern)

            pair.meanvcov <- stats::setNames(lapply(U.pattern,function(iPattern){
                test.in <- (pair.meanvcov[[iPattern]][1,] %in% name.effects)+(pair.meanvcov[[iPattern]][2,] %in% name.effects)
                return(pair.meanvcov[[iPattern]][,test.in==2,drop=FALSE])
            }), U.pattern)

            pair.vcov <- stats::setNames(lapply(U.pattern,function(iPattern){## iPattern <- 1
                test.in <- (pair.vcov[[iPattern]][1,] %in% name.effects)+(pair.vcov[[iPattern]][2,] %in% name.effects)
                iOut <- pair.vcov[[iPattern]][,test.in==2,drop=FALSE]
                attr(iOut,"subset") <- which(test.in==2)
                attr(iOut,"key") <- matrix(NA, nrow = length(name.varcoef[[iPattern]]), ncol = length(name.varcoef[[iPattern]]), dimnames = list(name.varcoef[[iPattern]],name.varcoef[[iPattern]]))
                for(iCol in 1:NCOL(iOut)){
                    attr(iOut, "key")[iOut[1,iCol],iOut[2,iCol]] <- iCol
                    attr(iOut, "key")[iOut[2,iCol],iOut[1,iCol]] <- iCol
                }
                return(iOut)
            }), U.pattern)
            d2Omega <- stats::setNames(lapply(U.pattern,function(iPattern){                
                return(d2Omega[[iPattern]][attr(pair.vcov[[iPattern]],"subset")])
            }), U.pattern)

            n.varcoef <- lengths(name.varcoef)
            name.allvarcoef <- unique(unlist(name.varcoef))
            n.allvarcoef <- length(name.allvarcoef)
            npair.meanvcov <- lapply(pair.meanvcov, NCOL)
            npair.vcov <- lapply(pair.vcov, NCOL)
        }
        if("mean" %in% effects == FALSE){ ## compute hessian only for variance and/or correlation parameters
            if(REML && indiv){
                stop("Not possible to compute individual hessian for variance and/or correlation coefficients when using REML.\n")
            }

            test.vcov <- any(n.varcoef>0)
            test.mean <- FALSE

        }else{ ## compute hessian all parameters
     
            test.vcov <- any(n.varcoef>0)
            test.mean <- n.mucoef>0
        }
    }

    ## prepare REML term
    if(test.vcov){

        ls.dOmega_OmegaM1 <- attr(dOmega,"ls.dOmega_OmegaM1")
        ls.OmegaM1_dOmega_OmegaM1 <- attr(dOmega,"ls.OmegaM1_dOmega_OmegaM1")
        dOmega_OmegaM1 <- attr(dOmega,"dOmega_OmegaM1")
        OmegaM1_dOmega_OmegaM1 <- attr(dOmega,"OmegaM1_dOmega_OmegaM1")

        if(REML){
            REML.key <- matrix(NA, nrow = n.allvarcoef, ncol = n.allvarcoef, dimnames = list(name.allvarcoef, name.allvarcoef))
            maxkey <- sum(lower.tri(REML.key, diag = TRUE))
            REML.key[lower.tri(REML.key, diag = TRUE)] <- 1:maxkey
            REML.key[upper.tri(REML.key)] <- t(REML.key)[upper.tri(REML.key)]

        REML.denom <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
        ## REML.numerator1 <- stats::setNames(lapply(U.pattern, function(iPattern) { array(0, dim = c(n.mucoef, n.mucoef, n.varcoef[[iPattern]]), dimnames = list(name.mucoef, name.mucoef, name.varcoef[[iPattern]])) }), U.pattern)
        REML.numerator1 <- array(0, dim = c(n.mucoef, n.mucoef, n.allvarcoef), dimnames = list(name.mucoef, name.mucoef, name.allvarcoef))
        REML.numerator2 <- array(0, dim = c(n.mucoef, n.mucoef, maxkey), dimnames = list(name.mucoef, name.mucoef, NULL))
        }
    }  

    ## ** compute information
    ## *** looping over individuals
    if(test.loopIndiv){
        
        ## ** precompute 
        if(test.vcov){
            tr_OmegaM1_d2OmegaAndCo <- stats::setNames(lapply(1:n.pattern, function(iPattern){rep(NA, npair.vcov[[iPattern]])}), U.pattern)
            if(REML || type.information == "observed"){
                OmegaM1_d2OmegaAndCo_OmegaM1 <- stats::setNames(lapply(1:n.pattern, function(iPattern){array(NA, dim = c(NCOL(precision[[iPattern]]),NCOL(precision[[iPattern]]), npair.vcov[[iPattern]]))}), U.pattern)
            }

            for(iPattern in 1:n.pattern){ ## iPattern <- 2
                iOmegaM1 <- precision[[iPattern]]
                idOmega <- dOmega[[iPattern]]

                ## loop over all pairs
                for(iPair in 1:npair.vcov[[iPattern]]){ ## iPair <- 4
                    iCoef1 <- pair.vcov[[iPattern]][1,iPair]
                    iCoef2 <- pair.vcov[[iPattern]][2,iPair]
                    iTerm21 <- ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iCoef2]] %*% idOmega[[iCoef1]]
                    
                    ## trace
                    if(type.information == "expected"){
                        tr_OmegaM1_d2OmegaAndCo[[iPattern]][iPair] <- tr(iTerm21)
                    }else if(type.information == "observed"){
                        tr_OmegaM1_d2OmegaAndCo[[iPattern]][iPair] <- - tr(iTerm21 - iOmegaM1 %*% d2Omega[[iPattern]][[iPair]])
                    }
                    if(REML || type.information == "observed"){
                        iTerm12 <- ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iCoef1]] %*% idOmega[[iCoef2]]
                        OmegaM1_d2OmegaAndCo_OmegaM1[[iPattern]][,,iPair] <- iOmegaM1 %*% d2Omega[[iPattern]][[iPair]] %*% iOmegaM1 - (iTerm12 + iTerm21) %*% iOmegaM1
                    }
                }
            }
        }
        
        ## loop
        for(iId in 1:n.cluster){ ## iId <- 7
            iPattern <- pattern[iId]
            iIndex <- index.cluster[[iId]]
            iWeight <- weights[iId]

            iX <- X[iIndex,,drop=FALSE]
            tiX <- t(iX)
            iOmegaM1 <- precision[[iPattern]] * scale.Omega[iId]
            if(type.information == "observed"){
                iResidual <- residuals[iIndex,,drop=FALSE]
            }
        
            ## **** mean,mean
            iValue <-  iWeight * (tiX %*% iOmegaM1 %*% iX)
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
                    REML.numerator1[,,iVarcoef] <- REML.numerator1[,,iVarcoef] + iWeight * (tiX %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iX) * scale.Omega[iId]
                }
            }

            ## **** var,var
            if(test.vcov){

                for(iPair in 1:npair.vcov[[iPattern]]){ ## iPair <- 1
                    iCoef1 <- pair.vcov[[iPattern]][1,iPair]
                    iCoef2 <- pair.vcov[[iPattern]][2,iPair]

                    iValue <- 0.5 * iWeight * tr_OmegaM1_d2OmegaAndCo[[iPattern]][iPair]
                    ## 0.5 * ntr(iOmegaM1 %*% idOmega$sigma %*% iOmegaM1 %*% idOmega$sigma)

                    if(type.information == "observed"){
                        iValue <- iValue - 0.5 * iWeight * (t(iResidual) %*% OmegaM1_d2OmegaAndCo_OmegaM1[[iPattern]][,,iPair] %*% iResidual) * scale.Omega[iId]
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
                        REML.numerator2[,,iKey] <- REML.numerator2[,,iKey] + iWeight * (tiX %*% OmegaM1_d2OmegaAndCo_OmegaM1[[iPattern]][,,iPair] %*% iX) * scale.Omega[iId]
                    }
                }
            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){

                for(iPair in 1:npair.meanvcov[[iPattern]]){ ## iPair <- 1
                    iCoef1 <- pair.meanvcov[[iPattern]][1,iPair]
                    iCoef2 <- pair.meanvcov[[iPattern]][2,iPair]

                    iValue <- iWeight * (tiX[iCoef1,,drop=FALSE] %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iCoef2]] %*% iResidual) * scale.Omega[iId]

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
            iOmegaM1 <- precision[[iPattern]]
            iTime <- NCOL(iOmegaM1)
            iTime2 <- length(iOmegaM1)
            iName.varcoef <- name.varcoef[[iPattern]]
            iN.varcoef <- length(iName.varcoef)
            iX <- precompute$XX$pattern[[iPattern]]

            ## **** mean,mean
            if(test.mean){
                info[name.mucoef,name.mucoef] <- info[name.mucoef,name.mucoef] + (as.double(iOmegaM1) %*% iX)[as.double(precompute$XX$key)]                
            }

            ## **** var,var
            if(test.vcov){
                
                ## precompute
                iLs.tdOmega_OmegaM1 <- lapply(ls.dOmega_OmegaM1[[iPattern]], FUN = base::t)
                tdOmega_OmegaM1 <- do.call(cbind,lapply(iLs.tdOmega_OmegaM1, as.numeric))
                iTrace_O_dO_O_dO <- colSums(dOmega_OmegaM1[[iPattern]][,pair.vcov[[iPattern]][1,],drop=FALSE] * tdOmega_OmegaM1[,pair.vcov[[iPattern]][2,],drop=FALSE])
                
                if(REML || type.information == "observed"){
                    iOmegaM1_d2Omega_OmegaM1 <- do.call(cbind,lapply(d2Omega[[iPattern]], function(iM){as.numeric(iOmegaM1 %*% iM %*% iOmegaM1)}))

                    iOmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1 <- do.call(cbind,lapply(1:npair.vcov[[iPattern]], function(iPair){ ## iPair <- 4
                        out <- iLs.tdOmega_OmegaM1[[pair.vcov[[iPattern]][1,iPair]]] %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[pair.vcov[[iPattern]][2,iPair]]]
                        return(as.double(out + t(out)))
                    }))
                    iOmegaM1_d2OmegaAndCo_OmegaM1 <- iOmegaM1_d2Omega_OmegaM1 - iOmegaM1_dOmega1_OmegaM1_dOmega2_OmegaM1
                }

                ## compute and store contribution
                if(type.information == "expected"){
                    iValue <- 0.5 * Upattern.ncluster[iPattern] * iTrace_O_dO_O_dO
                }else if(type.information == "observed"){
                    id2Omega <- matrix(unlist(d2Omega[[iPattern]]), nrow = iTime2, ncol = npair.vcov[[iPattern]])
                    iTrace_d2Omega <- colSums(sweep(id2Omega, MARGIN = 1, FUN = "*", STATS = as.double(precision[[iPattern]])))
                    iValue <- - 0.5 * Upattern.ncluster[iPattern] * (iTrace_O_dO_O_dO - iTrace_d2Omega) - 0.5 * as.double(as.double(precompute$RR[[iPattern]]) %*% iOmegaM1_d2OmegaAndCo_OmegaM1)
                }
                info[iName.varcoef,iName.varcoef]  <- info[iName.varcoef,iName.varcoef] + iValue[as.double(attr(pair.vcov[[iPattern]],"key"))]

                ## determinant
                if(REML){                    
                    iDouble2Mat <- as.vector(precompute$XX$key)
                    ## denominator
                    if(is.null(precompute$X.OmegaM1.X)){
                        REML.denom <- REML.denom + (as.double(iOmegaM1) %*% iX)[iDouble2Mat]
                    }else{
                        REML.denom <- REML.denom + precompute$X.OmegaM1.X[[iPattern]][iDouble2Mat]
                    }
                    ## numerator 1
                    iX_OmegaM1_dOmega_OmegaM1_X <- t(iX) %*% OmegaM1_dOmega_OmegaM1[[iPattern]]
                    for(iVarcoef in iName.varcoef){ ## iVarcoef <- iName.varcoef[1]
                        REML.numerator1[,,iVarcoef] <- REML.numerator1[,,iVarcoef] + iX_OmegaM1_dOmega_OmegaM1_X[iDouble2Mat,iVarcoef]
                    }
                    ## numerator 2
                    iX_OmegaM1_d2OmegaAndCo_OmegaM1_X <- t(iX) %*% iOmegaM1_d2OmegaAndCo_OmegaM1
                    for(iPair in 1:npair.vcov[[iPattern]]){ ## iPair <- 1
                        iCoef1 <- pair.vcov[[iPattern]][1,iPair]
                        iCoef2 <- pair.vcov[[iPattern]][2,iPair]
                        REML.numerator2[,,REML.key[iCoef1,iCoef2]] <- REML.numerator2[,,REML.key[iCoef1,iCoef2]] + iX_OmegaM1_d2OmegaAndCo_OmegaM1_X[iDouble2Mat,iPair]
                    }
                    
                }

            }

            ## **** mean,var
            if(type.information == "observed" && test.mean && test.vcov){

                ## compute
                iValue <- t(OmegaM1_dOmega_OmegaM1[[iPattern]]) %*% matrix(precompute$XR[[iPattern]], nrow = iTime2, ncol = n.mucoef, dimnames = list(NULL,name.mucoef))

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
        if(REML){
            effects2 <- "mean"
            attr(effects2,"original.names") <- attr(effects,"original.names")
            attr(effects2,"reparametrize.names") <- attr(effects,"reparametrize.names")
        }else{
            effects2 <- effects
        }
        if(is.null(weights)){
            weights <- rep(1, length(pattern))
        }
        if(is.null(scale.Omega)){
            scale.Omega <- rep(1, length(pattern))
        }
        attr.info <- info
        attr.bread <- crossprod(.score(X = X, residuals = residuals, precision = precision, dOmega = dOmega,
                                       weights = weights, scale.Omega = scale.Omega, pattern = pattern, 
                                       index.cluster = index.cluster, name.allcoef = name.allcoef, indiv = TRUE, REML = REML, effects = effects2,
                                       precompute = precompute) )

        
        
        if(any(c("mean","variance","correlation") %in% effects2 == FALSE)){
            keep.cols <- intersect(names(which(rowSums(abs(attr.bread))!=0)),names(which(rowSums(abs(attr.bread))!=0)))
            info <- NA*attr.info

            attr.infoM1 <- solve(attr.info)
            info[keep.cols,keep.cols] <- solve(attr.infoM1[keep.cols,keep.cols] %*% attr.bread[keep.cols,keep.cols,drop=FALSE] %*% attr.infoM1[keep.cols,keep.cols])
        }else{
            attr.infoM1 <- solve(attr.info)
            info <- solve(attr.infoM1 %*% attr.bread %*% attr.infoM1)
        }
    }
    return(info)
}

##----------------------------------------------------------------------
### information.R ends here

