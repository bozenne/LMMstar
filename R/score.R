### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:59) 
## Version: 
## Last-Updated: Mar 26 2024 (09:51) 
##           By: Brice Ozenne
##     Update #: 599
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * score.lmm (documentation)
##' @title Extract The Score From a Linear Mixed Model
##' @description Extract or compute the first derivative of the log-likelihood of a linear mixed model.
##' 
##' @param x a \code{lmm} object.
##' @param newdata [data.frame] dataset relative to which the score should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the score be output? Otherwise output the sum of all clusters of the derivatives.
##' @param effects [character] Should the score relative to all coefficients be output (\code{"all"}),
##' or only coefficients relative to the mean (\code{"mean"} or \code{"fixed"}),
##' or only coefficients relative to the variance and correlation structure (\code{"variance"} or \code{"correlation"}).
##' @param p [numeric vector] value of the model coefficients at which to evaluate the score. Only relevant if differs from the fitted values.
##' @param transform.sigma [character] Transformation used on the variance coefficient for the reference level. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"} - see details.
##' @param transform.k [character] Transformation used on the variance coefficients relative to the other levels. One of \code{"none"}, \code{"log"}, \code{"square"}, \code{"logsquare"}, \code{"sd"}, \code{"logsd"}, \code{"var"}, \code{"logvar"} - see details.
##' @param transform.rho [character] Transformation used on the correlation coefficients. One of \code{"none"}, \code{"atanh"}, \code{"cov"} - see details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details For details about the arguments \bold{transform.sigma}, \bold{transform.k}, \bold{transform.rho}, see the documentation of the \link[LMMstar]{coef.lmm} function.
##'
##' @return
##' When argument indiv is \code{FALSE}, a vector with the value of the score relative to each coefficient.
##' When argument indiv is \code{TRUE}, a matrix with the value of the score relative to each coefficient (in columns) and each cluster (in rows).
##'
##' @keywords methods

## * score.lmm (code)
##' @export
score.lmm <- function(x, effects = "mean", newdata = NULL, p = NULL, indiv = FALSE, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

    ## ** normalize user input
    dots <- list(...)
    options <- LMMstar.options()
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
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

    ## ** extract or recompute score
    if(is.null(newdata) && is.null(p) && (indiv == FALSE) && test.notransform){
        keep.name <- stats::setNames(names(coef(x, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                                     names(coef(x, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))    

        design <- x$design ## useful in case of NA
        out <- x$score[keep.name]
        if(transform.names){
            names(out) <- names(keep.name)
        }
    }else{
        test.precompute <- !is.null(x$design$precompute.XX) && !indiv

        if(!is.null(newdata)){
            design <- stats::model.matrix(x, newdata = newdata, effects = "all", simplify = FALSE)
        }else{
            design <- x$design
        }
        out <- .moments.lmm(value = theta, design = design, time = x$time, method.fit = x$args$method.fit,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                            logLik = FALSE, score = TRUE, information = FALSE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects,
                            trace = FALSE, precompute.moments = test.precompute, transform.names = transform.names)$score
    }

    ## ** name and restaure NAs
    if(indiv){

        if(!is.numeric(x$cluster$levels)){
            rownames(out) <- x$cluster$levels[match(1:NROW(out),x$cluster$index)]
        } 
        out <- restaureNA(out, index.na = x$index.na,
                          level = "cluster", cluster = x$cluster)        
        
    }

    ## ** export
    return(out)
}

## * .score
.score <- function(X, residuals, precision, dOmega,
                   Upattern.ncluster, weights, scale.Omega,
                   pattern, index.cluster, name.allcoef,
                   indiv, REML, effects,
                   precompute){


    ## ** extract information
    test.loopIndiv <- indiv || is.null(precompute)
    n.obs <- length(index.cluster)
    n.cluster <- length(pattern)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- lapply(dOmega,names)
    n.varcoef <- lengths(name.varcoef)
    name.allvarcoef <- unique(unlist(name.varcoef))
    U.pattern <- names(dOmega)
    n.pattern <- length(U.pattern)

    ## ** prepare output
    name.effects <- attr(effects,"original.names")
    n.effects <- length(name.effects)
    if(test.loopIndiv){
        Score <- matrix(0, nrow = n.cluster, ncol = n.effects,
                        dimnames = list(NULL, name.effects))
    }else{
        Score <- stats::setNames(rep(0, n.effects), name.effects)
    }
    if(any(is.na(attr(precision, "logdet")))){ ## non positive definite residual variance covariance
        return(Score*NA)
    }

    ## restrict to relevant parameters
    if(("variance" %in% effects == FALSE) && ("correlation" %in% effects == FALSE)){ ## compute score only for mean parameters
        test.vcov <- FALSE
        test.mean <- n.mucoef>0
    }else{
        if(REML && indiv){
            stop("Not possible to compute individual score for variance and/or correlation coefficients when using REML.\n")
        }
        if(("variance" %in% effects == FALSE) || ("correlation" %in% effects == FALSE)){
            name.varcoef <- stats::setNames(lapply(U.pattern,function(iPattern){intersect(name.effects,name.varcoef[[iPattern]])}),
                                     U.pattern)
            n.varcoef <- lengths(name.varcoef)
            name.allvarcoef <- unique(unlist(name.varcoef))
        }
        if("mean" %in% effects == FALSE){ ## compute score only for variance and/or correlation parameters
            if(REML && indiv){
                stop("Not possible to compute individual score for variance and/or correlation coefficients when using REML.\n")
            }

            test.vcov <- any(n.varcoef>0)
            test.mean <- FALSE

        }else{ ## compute score all parameters
     
            test.vcov <- any(n.varcoef>0)
            test.mean <- n.mucoef>0
        }
    }

    if(test.vcov){
        ls.OmegaM1_dOmega_OmegaM1 <- attr(dOmega,"ls.OmegaM1_dOmega_OmegaM1")
        OmegaM1_dOmega_OmegaM1 <- attr(dOmega,"OmegaM1_dOmega_OmegaM1")
                
        if(REML){
            REML.num <- array(0, dim = c(n.mucoef, n.mucoef, length(name.allvarcoef)), dimnames = list(name.mucoef,name.mucoef,name.allvarcoef))
            REML.denom <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))
        }
    }

    ## ** compute score
    ## *** looping over individuals
    if(test.loopIndiv){

        if(test.vcov){ ## precompute
            trOmegaM1_dOmega <- stats::setNames(vector(mode = "list", length = n.pattern), U.pattern)
            for(iPattern in 1:n.pattern){ ## iPattern <- 4
                trOmegaM1_dOmega[[iPattern]]  <- stats::setNames(lapply(name.varcoef[[iPattern]], function(iVarcoef){tr(precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]])}), name.varcoef[[iPattern]])
            }
        }
        
        ## loop
        for(iId in 1:n.cluster){ ## iId <- 7
            iPattern <- pattern[iId]
            iIndex <- index.cluster[[iId]]
            iWeight <- weights[iId]
            iOmegaM1 <- precision[[pattern[iId]]] * scale.Omega[iId]

            iResidual <- residuals[iIndex,,drop=FALSE]
            iX <- X[iIndex,,drop=FALSE]
            tiX <- t(iX)

            if(test.mean){
                Score[iId,name.mucoef] <- iWeight * (tiX %*% iOmegaM1 %*% iResidual)
            }

            if(test.vcov){
                if(REML){
                    REML.denom <- REML.denom + iWeight * (tiX %*% iOmegaM1 %*% iX)
                }

                for(iVarcoef in name.varcoef[[iPattern]]){ ## iVarcoef <- name.varcoef[1]
                    Score[iId,iVarcoef] <- -0.5 * iWeight * trOmegaM1_dOmega[[iPattern]][[iVarcoef]] + 0.5 * iWeight * (t(iResidual) %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iResidual) * scale.Omega[iId]

                    if(REML){
                        REML.num[,,iVarcoef] <- REML.num[,,iVarcoef] + iWeight * (tiX %*% ls.OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iX) * scale.Omega[iId]
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

        ## loop
        for (iPattern in U.pattern) { ## iPattern <- U.pattern[1]
            iName.varcoef <- name.varcoef[[iPattern]]
            iOmegaM1 <- precision[[iPattern]]
            iTime2 <- length(iOmegaM1)

            if(test.mean){
                ## X %*% iOmega^-1 %*% residual
                Score[name.mucoef] <- Score[name.mucoef] + as.double(iOmegaM1) %*% matrix(unlist(precompute$XR[[iPattern]]), nrow = iTime2, ncol = n.mucoef, byrow = FALSE)
                ## Score[name.mucoef] <- Score[name.mucoef] + apply(precompute$XR[[iPattern]], MARGIN = 3, FUN = function(iM){sum(iM * iOmegaM1)})
            }

            if(test.vcov){
                iTrace <- sapply(dOmega[[iPattern]], FUN = function(iO){sum(iO*precision[[iPattern]])})
                ## iOmega^-1 dOmega iOmega^-1
                Score[iName.varcoef] <- Score[iName.varcoef] - 0.5 * Upattern.ncluster[iPattern] * iTrace + 0.5 * as.double(precompute$RR[[iPattern]]) %*% OmegaM1_dOmega_OmegaM1[[iPattern]]
                
                if(REML){
                    iX <- precompute$XX$pattern[[iPattern]]
                    iDouble2Mat <- as.vector(precompute$XX$key)
                    ## denominator
                    if(is.null(precompute$X.OmegaM1.X)){
                        REML.denom <- REML.denom + (as.double(iOmegaM1) %*% iX)[iDouble2Mat]
                    }else{
                        REML.denom <- REML.denom + precompute$X.OmegaM1.X[[iPattern]][iDouble2Mat]
                    }
                    ## numerator
                    iX_OmegaM1_dOmega_OmegaM1_X <- t(iX) %*% OmegaM1_dOmega_OmegaM1[[iPattern]]
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
