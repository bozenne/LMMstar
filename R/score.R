### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:59) 
## Version: 
## Last-Updated: Jun 17 2021 (12:50) 
##           By: Brice Ozenne
##     Update #: 408
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * score.lmm (documentation)
##' @title Extract The Score From a Multivariate Gaussian Model
##' @description Extract or compute the first derivative of the log-likelihood of a multivariate gaussian model.
##' @name score
##' 
##' @param x a \code{lmm} object.
##' @param data [data.frame] dataset relative to which the score should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the score be output? Otherwise output the sum of all clusters of the derivatives.
##' @param effects [character] Should the score relative to all coefficients be output (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"}),
##' or only coefficients relative to the variance and correlation structure (\code{"variance"} or \code{"correlation"}).
##' @param p [numeric vector] value of the model coefficients at which to evaluate the score. Only relevant if differs from the fitted values.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef} function.
##'
##' @return
##' When argument indiv is \code{FALSE}, a vector with the value of the score relative to each coefficient.
##' When argument indiv is \code{TRUE}, a matrix with the value of the score relative to each coefficient (in columns) and each cluster (in rows).
##' 

## * score.lmm (code)
##' @rdname score
##' @export
score.lmm <- function(x, effects = "all", data = NULL, p = NULL, indiv = FALSE, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
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

    ## ** extract or recompute score
    if(is.null(data) && is.null(p) && (indiv == FALSE) && test.notransform){
        design <- x$design ## useful in case of NA
        out <- x$score
        if(transform.names && !is.null(x$reparametrize$newname)){
            names(out)[match(names(x$reparametrize$p),names(out))] <- x$reparametrize$newname
        }
        if("mean" %in% effects == FALSE){
            out <- out[x$param$type!="mu"]
        }else if("variance" %in% effects == FALSE && "correlation" %in% effects == FALSE){
            out <- out[x$param$type=="mu"]
        }
    }else{

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
            name.allcoef <- names(x$param$value)
        }
        Y <- design$Y
        X <- design$X.mean
        index.vargroup <- design$X.var$cluster
        index.cluster <- design$index.cluster
        index.time <- design$index.time
        X.var <- design$X.var
        name.varcoef <- design$X.var$param
        index.var <- x$param$type %in% c("sigma","k","rho")
        if(test.precompute){
            precompute <- list(XX = design$precompute.XX)
        }else{
            precompute <- NULL
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
                p <- x$param$value
            }
            beta <- p[x$param$type=="mu"]
            reparametrize <- .reparametrize(p = p[index.var], type = x$param$type[index.var], strata = x$param$strata[index.var], time.levels = x$time$levels,
                                            time.k = x$design$param$time.k, time.rho = x$design$param$time.rho,
                                            Jacobian = TRUE, dJacobian = FALSE, inverse = FALSE,
                                            transform.sigma = transform.sigma,
                                            transform.k = transform.k,
                                            transform.rho = transform.rho,
                                            transform.names = TRUE)
            if(reparametrize$transform==FALSE){
                reparametrize$newname <- NULL
                reparametrize$Jacobian <- NULL
                reparametrize$dJacobian <- NULL
            }
            newname.allcoef <- reparametrize$newname
            Omega <- .calc_Omega(object = X.var, param = p, keep.interim = TRUE)
            dOmega <- .calc_dOmega(object = X.var, param = p, type = x$param$type, Omega = Omega, Jacobian = reparametrize$Jacobian)
            precision <- lapply(Omega, solve)
        }else{
            newname.allcoef <- x$reparametrize$newname
            beta <- x$param$value[x$param$type=="mu"]
            precision <- x$OmegaM1
            dOmega <- x$dOmega
        }
        
        out <- .score(X = X, residuals = Y - X %*% beta, precision = precision, dOmega = dOmega,
                      index.variance = index.vargroup, time.variance = index.time, index.cluster = index.cluster, name.varcoef = name.varcoef, name.allcoef = name.allcoef,
                      indiv = indiv, REML = x$method.fit=="REML", effects = effects, precompute = precompute, X.var = X.var)

        if(transform.names && length(newname.allcoef)>0 && ("variance" %in% effects || "correlation" %in% effects)){
            if(indiv){
                index.var2 <- match(name.allcoef[index.var], colnames(out))
                colnames(out)[stats::na.omit(index.var2)] <- newname.allcoef[!is.na(index.var2)]
            }else{
                index.var2 <- match(name.allcoef[index.var], names(out))
                names(out)[stats::na.omit(index.var2)] <- newname.allcoef[!is.na(index.var2)]
            }
        }

    }

    ## ** restaure NA
    if(length(x$index.na)>0 && indiv){ 
        iAdd <- .addNA(index.na = x$index.na, design = design, time = x$time)
        if(length(iAdd$missing.cluster)>0){
            out.save <- out
            out <- matrix(NA, nrow = iAdd$n.allcluster, ncol = NCOL(out.save),
                          dimnames = list(NULL, colnames(out.save)))
            out[match(design$cluster$levels, iAdd$allcluster),] <- out.save
        }
    }

    ## ** export
    return(out)
}

## * .score
.score <- function(X, residuals, precision, dOmega,
                   index.variance, time.variance, index.cluster, name.varcoef, name.allcoef,
                   indiv, REML, effects, precompute, X.var){

    ## ** extract information
    test.loopIndiv <- indiv || is.null(precompute)
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    n.varcoef <- lapply(name.varcoef, length)
    name.allvarcoef <- unique(unlist(name.varcoef))
    n.allcoef <- length(name.allcoef)
    U.pattern <- names(dOmega)
    n.pattern <- length(U.pattern)

    if(!is.null(precompute) && "RR" %in% names(precompute) == FALSE){
        precompute$RR <-  .precomputeRR(residuals = residuals, pattern = X.var$pattern,
                                        pattern.time = X.var$index.time, pattern.cluster = attr(X.var$cluster, "index.byPattern"), index.cluster = attr(index.cluster,"sorted"))
    }
    if(!is.null(precompute) && "XR" %in% names(precompute) == FALSE){
        precompute$XR <-  .precomputeXR(X = X, residuals = residuals, pattern = X.var$pattern,
                                        pattern.time = X.var$index.time, pattern.cluster = attr(X.var$cluster, "index.byPattern"), index.cluster = attr(index.cluster,"sorted"))
    }

    ## ** prepare output
    if("mean" %in% effects == FALSE){ ## compute information only for variance - correlation parameters
        if(test.loopIndiv){
            Score <- matrix(0, nrow = n.cluster, ncol = length(name.allvarcoef),
                            dimnames = list(NULL, name.allvarcoef))
        }else{
            Score <- stats::setNames(rep(0, length(name.allvarcoef)), name.allvarcoef)
        }
        test.vcov <- TRUE
        test.mean <- FALSE
        if(REML && indiv){
            stop("Not possible to compute individual score for variance and/or correlation coefficients when using REML.\n")
        }
    }else if("variance" %in% effects == FALSE && "correlation" %in% effects == FALSE){ ## compute information only for mean parameters
        if(test.loopIndiv){
            Score <- matrix(0, nrow = n.cluster, ncol = n.mucoef,
                            dimnames = list(NULL, name.mucoef))
        }else{
            Score <- stats::setNames(rep(0, n.mucoef), name.mucoef)
        }
        test.vcov <- FALSE
        test.mean <- TRUE
    }else{ ## compute information only for all parameters
        if(test.loopIndiv){
            Score <- matrix(0, nrow = n.cluster, ncol = n.allcoef,
                            dimnames = list(NULL, name.allcoef))
        }else{
            Score <- stats::setNames(rep(0, n.allcoef), name.allcoef)
        }
        test.vcov <- TRUE
        test.mean <- TRUE
        if(REML && indiv){
            stop("Not possible to compute individual score for variance and/or correlation coefficients when using REML.\n")
        }
    }

    if(test.vcov && REML){
        REML.num <- array(0, dim = c(n.mucoef, n.mucoef, length(name.allvarcoef)), dimnames = list(name.mucoef,name.mucoef,name.allvarcoef))
        REML.denom <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
    }

    ## ** compute score
    ## *** looping over individuals
    if(test.loopIndiv){
        
        if(test.vcov){ ## precompute
            trOmegaM1_dOmega <- stats::setNames(vector(mode = "list", length = n.pattern), U.pattern)
            OmegaM1_dOmega_OmegaM1 <- stats::setNames(vector(mode = "list", length = n.pattern), U.pattern)

            for(iPattern in 1:n.pattern){ ## iPattern <- 4
                trOmegaM1_dOmega[[iPattern]]  <- stats::setNames(lapply(name.varcoef[[iPattern]], function(iVarcoef){tr(precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]])}), name.varcoef[[iPattern]])
                OmegaM1_dOmega_OmegaM1[[iPattern]] <- stats::setNames(lapply(name.varcoef[[iPattern]], function(iVarcoef){precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]] %*% precision[[iPattern]]}), name.varcoef[[iPattern]])
            }
        }
        
        ## loop
        for(iId in 1:n.cluster){ ## iId <- 7
            iPattern <- index.variance[iId]
            iIndex <- attr(index.cluster,"sorted")[[iId]]
            ## iIndex <- which(index.cluster==iId)
            ## iIndex <- iIndex[order(time.variance[iIndex])] ## re-order observations according to the variance-covariance matrix

            iResidual <- residuals[iIndex,,drop=FALSE]
            iX <- X[iIndex,,drop=FALSE]
            tiX <- t(iX)

            if(test.mean){
                Score[iId,name.mucoef] <- tiX %*% precision[[iPattern]] %*% iResidual
            }

            if(test.vcov){
                if(REML){
                    REML.denom <- REML.denom + tiX %*% precision[[iPattern]] %*% iX
                }

                for(iVarcoef in name.varcoef[[iPattern]]){ ## iVarcoef <- name.varcoef[1]
                    Score[iId,iVarcoef] <- -0.5 * trOmegaM1_dOmega[[iPattern]][[iVarcoef]] + 0.5 * t(iResidual) %*% OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iResidual

                    if(REML){
                        REML.num[,,iVarcoef] <- REML.num[,,iVarcoef] + tiX %*% OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iX
                    }
                }
            }
        }
        if(!indiv){
            Score <- colSums(Score)
        }
    }

    ## *** looping over covariance patterns
    if(!test.loopIndiv){
        ## precompute
        ncluster.pattern <- sapply(attr(index.variance,"index.byPattern"),length)
        name.pattern <- names(ncluster.pattern)

        ## loop
        for (iPattern in name.pattern) { ## iPattern <- name.pattern[1]
            iName.varcoef <- name.varcoef[[iPattern]]
            iOmega <- precision[[iPattern]]
            iTime2 <- length(iOmega)


            if(test.mean){
                ## X %*% iOmega^-1 %*% residual
                Score[name.mucoef] <- Score[name.mucoef] + as.double(iOmega) %*% matrix(unlist(precompute$XR[[iPattern]]), nrow = iTime2, ncol = n.mucoef, byrow = FALSE)
                ## Score[name.mucoef] <- Score[name.mucoef] + apply(precompute$XR[[iPattern]], MARGIN = 3, FUN = function(iM){sum(iM * iOmega)})
            }

            if(test.vcov){
                idOmega <- matrix(unlist(dOmega[[iPattern]]), nrow = iTime2, ncol = n.varcoef[[iPattern]], dimnames = list(NULL,name.varcoef[[iPattern]]))
                iTrace <- stats::setNames(colSums(sweep(idOmega, MARGIN = 1, FUN = "*", STATS = as.double(precision[[iPattern]]))), name.varcoef[[iPattern]])
                
                ## iOmega^-1 dOmega iOmega^-1
                iOmegaM1_dOmega_OmegaM1 <- matrix(iOmega %*% tblock(t(do.call(rbind, dOmega[[iPattern]]) %*% iOmega)),
                                                  nrow = iTime2, ncol = length(name.varcoef[[iPattern]]), dimnames = list(NULL,name.varcoef[[iPattern]]), byrow = FALSE)

                Score[iName.varcoef] <- Score[iName.varcoef] - 0.5 * ncluster.pattern[iPattern] * iTrace + 0.5 * as.double(precompute$RR[[iPattern]]) %*% iOmegaM1_dOmega_OmegaM1
                
                if(REML){
                    iX <- matrix(precompute$XX$pattern[[iPattern]], nrow = iTime2, ncol = dim(precompute$XX$pattern[[iPattern]])[3], byrow = FALSE)
                    iDouble2Mat <- as.vector(precompute$XX$key)
                    ## denominator
                    REML.denom <- REML.denom + (as.double(iOmega) %*% iX)[iDouble2Mat]
                    ## numerator
                    iX_OmegaM1_dOmega_OmegaM1_X <- t(iX) %*% iOmegaM1_dOmega_OmegaM1
                    for(iVarcoef in iName.varcoef){ ## iVarcoef <- iName.varcoef[1]
                        REML.num[,,iVarcoef] <- REML.num[,,iVarcoef] + iX_OmegaM1_dOmega_OmegaM1_X[iDouble2Mat,iVarcoef]
                    }
                }
            }
        }
    }

    ## ** export
    if(REML && test.vcov){
        REML.denomM1 <- solve(REML.denom)

        ## compute: 0.5 tr((X\OmegaM1X)^-1 (X\OmegaM1 d\Omega \OmegaM1 X)) in one go
        ## Score[name.allvarcoef] <-  Score[name.allvarcoef] + 0.5 * apply(REML.num, MARGIN = 3, function(x){tr(REML.denomM1 %*% x)})
        Score[name.allvarcoef] <-  Score[name.allvarcoef] + 0.5 * as.double(REML.denomM1) %*% matrix(REML.num, nrow = prod(dim(REML.num)[1:2]), ncol = dim(REML.num)[3], byrow = FALSE)
    }
    return(Score)
}


##----------------------------------------------------------------------
### score.R ends here
