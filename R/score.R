### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:59) 
## Version: 
## Last-Updated: aug  3 2023 (17:28) 
##           By: Brice Ozenne
##     Update #: 655
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
##' @param data [data.frame] dataset relative to which the score should be computed. Only relevant if differs from the dataset used to fit the model.
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
score.lmm <- function(x, effects = "mean", data = NULL, p = NULL, indiv = FALSE, transform.sigma = NULL, transform.k = NULL, transform.rho = NULL, transform.names = TRUE, ...){

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

    init <- .init_transform(transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, 
                            x.transform.sigma = x$reparametrize$transform.sigma, x.transform.k = x$reparametrize$transform.k, x.transform.rho = x$reparametrize$transform.rho)
    transform.sigma <- init$transform.sigma
    transform.k <- init$transform.k
    transform.rho <- init$transform.rho
    test.notransform <- init$test.notransform

    ## ** extract or recompute score
    if(is.null(data) && is.null(p) && (indiv == FALSE) && test.notransform){
            keep.name <- stats::setNames(names(coef(x, effects = effects, transform.sigma = "none", transform.k = "none", transform.rho = "none", transform.names = TRUE)),
                                                     names(coef(x, effects = effects, transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho, transform.names = transform.names)))    

            design <- x$design ## useful in case of NA
            out <- x$score[keep.name]
            if(transform.names){
                names(out) <- names(keep.name)
            }
    }else{

        test.precompute <- !is.null(x$design$precompute.XX) && !indiv

        if(!is.null(data)){
            design <- stats::model.matrix(x, data = data, effects = "all", simplify = FALSE)
        }else{
            design <- x$design
        }

        if(!is.null(p)){
            if(any(duplicated(names(p)))){
                stop("Incorrect argument \'p\': contain duplicated names \"",paste(unique(names(p)[duplicated(names(p))]), collapse = "\" \""),"\".\n")
            }
            if(any(names(x$param) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param)[names(x$param) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            p <- p[names(x$param)]
        }else{
            p <- x$param
        }
        out <- .moments.lmm(value = p, design = design, time = x$time, method.fit = x$args$method.fit,
                            transform.sigma = transform.sigma, transform.k = transform.k, transform.rho = transform.rho,
                            logLik = FALSE, score = TRUE, information = FALSE, vcov = FALSE, df = FALSE, indiv = indiv, effects = effects,
                            trace = FALSE, precompute.moments = test.precompute, transform.names = transform.names)$score
    }

    ## ** name and restaure NAs
    if(indiv){

        if(!is.numeric(x$cluster$levels)){
            rownames(out) <- x$cluster$levels[match(1:NROW(out),x$cluster$index)]
        } 
        out <- addNA(out, index.na = x$index.na,
                     level = "cluster", cluster = x$cluster)        
        
    }

    ## ** export
    return(out)
}

## * .score
.score <- function(X, residuals, precompute,
                   Upattern.ncluster, weights, scale.Omega,
                   pattern, index.cluster, 
                   indiv, REML, effects){


    ## ** extract information
    if(indiv && REML && c("variance","correlation") %in% effects){
        stop("Cannot compute cluster-specific score contributions relative to the variance-covariance parameters when using REML. \n")
    }
    test.loopIndiv <- indiv || !attr(precompute,"moments")
    
    name.effects <- attr(effects,"original.names")
    n.effects <- length(name.effects)
    if(!indiv && any(is.na(precompute$Omega.logdet))){ ## non positive definite residual variance covariance
        stats::setNames(rep(NA, n.effects), name.effects)
    }

    name.mucoef <- intersect(colnames(X), name.effects)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- lapply(precompute$dOmega_OmegaM1, function(iO){intersect(names(iO),name.effects)})
    name.allvarcoef <- unique(unlist(name.varcoef))
    n.allvarcoef <- length(name.allvarcoef)

    ## ** compute score
    ## *** looping over individuals
    if(test.loopIndiv){
        n.cluster <- length(pattern)
        Score <- matrix(0, nrow = n.cluster, ncol = n.effects,
                        dimnames = list(NULL, name.effects))
        ## loop
        for(iId in 1:n.cluster){ ## iId <- 7
            iPattern <- pattern[iId]
            iIndex <- index.cluster[[iId]]
            iOmegaM1 <- precompute$OmegaM1[[pattern[iId]]] * scale.Omega[iId]
            iResidual <- residuals[iIndex,,drop=FALSE]

            if(n.mucoef>0){
                Score[iId,name.mucoef] <- weights[iId] * (t(X[iIndex,name.mucoef,drop=FALSE]) %*% iOmegaM1 %*% iResidual)
            }

            if(n.allvarcoef>0){
                for(iVarcoef in name.varcoef[[iPattern]]){ ## iVarcoef <- name.varcoef[1]
                    term1 <- attr(precompute$dOmega_OmegaM1[[iPattern]],"tr")[iVarcoef]
                    term2 <- (t(iResidual) %*% precompute$OmegaM1_dOmega_OmegaM1[[iPattern]][[iVarcoef]] %*% iResidual) * scale.Omega[iId]
                    Score[iId,iVarcoef] <- 0.5 * weights[iId] * (- term1 + term2)
                }
            }
        }

        if(!indiv){
            Score <- colSums(Score)
        }
    }

    ## *** looping over covariance patterns
    if(!test.loopIndiv){
        U.pattern <- names(Upattern.ncluster)
        n.pattern <- length(U.pattern)
        Score <- stats::setNames(rep(0, n.effects), name.effects)

        ## loop
        for (iPattern in U.pattern) { ## iPattern <- U.pattern[1]
            
            iName.varcoef <- name.varcoef[[iPattern]]

            if(n.mucoef>0){
                ## X %*% iOmega^-1 %*% residual
                Score[name.mucoef] <- Score[name.mucoef] + (attr(precompute$wXR[[iPattern]],"vectorwise")[name.mucoef,,drop=FALSE] %*% attr(precompute$OmegaM1[[iPattern]],"vectorwise"))[,1]
                ## Score[name.mucoef] <- Score[name.mucoef] + apply(XR[[iPattern]], MARGIN = 3, FUN = function(iM){sum(iM * iOmegaM1)})
            }

            if(n.allvarcoef>0){
                term1 <- Upattern.ncluster[iPattern] * attr(precompute$dOmega_OmegaM1[[iPattern]],"tr")[iName.varcoef]
                term2 <- (attr(precompute$OmegaM1_dOmega_OmegaM1[[iPattern]],"vectorwise")[iName.varcoef,,drop=FALSE] %*% attr(precompute$wRR[[iPattern]],"vectorwise"))[,1]
                Score[iName.varcoef] <- Score[iName.varcoef] + 0.5 * (-term1 + term2)
            }
        }
    }

    ## *** add REML term
    if(REML && (n.allvarcoef>0)){
        ## compute: 0.5 tr((X\OmegaM1X)^-1 (X\OmegaM1 d\Omega \OmegaM1 X)) in one go
        ## Score[name.allvarcoef] <-  Score[name.allvarcoef] + 0.5 * sapply(wXOmegaM1dOmegaOmegaM1X, function(x){tr(wXOmegaM1X.M1 %*% x)})
        term1 <- attr(precompute$REML$wXOmegaM1dOmegaOmegaM1X,"vectorwise")[name.allvarcoef,,drop=FALSE] %*% attr(precompute$REML$wXOmegaM1X.M1,"vectorwise")
        Score[name.allvarcoef] <-  Score[name.allvarcoef] + 0.5 * term1[,1]
    }

    ## ** export
    return(Score)
}


##----------------------------------------------------------------------
### score.R ends here
