### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:59) 
## Version: 
## Last-Updated: Apr 22 2021 (18:11) 
##           By: Brice Ozenne
##     Update #: 171
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
##' @name score
##' 
##' @param object a \code{lmm} object.
##' @param data [data.frame] dataset relative to which the score should be computed. Only relevant if differs from the dataset used to fit the model.
##' @param indiv [logical] Should the contribution of each cluster to the score be output? Otherwise output the sum of all clusters of the derivatives.
##' @param p [numeric vector] value of the model coefficients at which to evaluate the score. Only relevant if differs from the fitted values.
##' @param transform [0,1,2] Transformation used on the variance coefficient. See details.
##' @param transform.names [logical] Should the name of the coefficients be updated to reflect the transformation that has been used?
##' @param ... Not used. For compatibility with the generic method.
##'
##' @details \bold{transform}: \cr
##' \itemize{
##' \item 0 means no transformation i.e. ouput stanrdard error, ratio of standard errors, and correlations.
##' \item 1 means log/atanh transformation i.e. ouput log(stanrdard error), log(ratio of standard errors), and atanh(correlations).
##' \item 2 ouput variance coefficients and correlations.
##' }
##'
##' @return
##' When argument indiv is \code{FALSE}, a vector with the value of the score relative to each coefficient.
##' When argument indiv is \code{FALSE}, a matrix with the value of the score relative to each coefficient (in columns) and each cluster (in rows).
##' 

## * score.lmm (code)
##' @rdname score
##' @export
score.lmm <- function(x, data = NULL, p = NULL, transform = NULL, indiv = FALSE, transform.names = TRUE, ...){
    options <- LMMstar.options()
    x.transform <- x$transform

    ## ** normalize user input
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    if(is.null(transform)){transform <- options$transform}

    ## ** extract or recompute score
    if(is.null(data) && is.null(p) && (indiv == FALSE) && (transform == x.transform)){
        out <- x$score
    }else{
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
                                        transform = transform
                                        )
            Y <- design$Y
            X <- design$X.mean
            index.vargroup <- design$index.vargroup
            index.cluster <- design$index.cluster
            index.time <- design$index.time
            X.var <- design$X.var
        }else{
            Y <- x$design$Y
            X <- x$design$X.mean
            index.vargroup <- x$design$index.vargroup
            index.cluster <- x$design$index.cluster
            index.time <- x$design$index.time
            X.var <- x$design$X.var
        }
        if(!is.null(p) || (transform != x.transform)){
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
            Omega <- .calc_Omega(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$rho)], keep.interim = TRUE)
            dOmega <- .calc_dOmega(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$rho)], Omega = Omega, transform = transform)
            precision <- lapply(Omega, solve)
        }else{
            beta <- x$param$mu
            precision <- x$OmegaM1
            dOmega <- x$dOmega
        }        
        out <- .score(X = X, residuals = Y - X %*% beta, precision = precision, dOmega = dOmega,
                      index.variance = index.vargroup, time.variance = index.time, index.cluster = index.cluster, ## attr(X.var,"Upattern.index.time")
                      indiv = indiv, REML = x$method.fit=="REML")

        if(transform>0 && transform.names){
            if(indiv){
                colnames(out) <- names(coef(e.lmm, transform = transform, effects = "all", transform.names = transform.names))
            }else{
                names(out) <- names(coef(e.lmm, transform = transform, effects = "all", transform.names = transform.names))
            }
        }

    }

    ## ** export
    return(out)
}

## * .score
.score <- function(X, residuals, precision, dOmega,
                   index.variance, time.variance, index.cluster, 
                   indiv, REML){

    if(indiv && REML){
        stop("Not possible to compute individual score when using REML.\n")
    }

    ## ** prepare
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    name.mucoef <- colnames(X)
    n.mucoef <- length(name.mucoef)
    name.varcoef <- names(dOmega[[1]])
    n.varcoef <- length(name.varcoef)
    name.allcoef <- c(name.mucoef,name.varcoef)
    n.allcoef <- length(name.allcoef)
    var.pattern <- names(dOmega)

    Score <- matrix(NA, nrow = n.cluster, ncol = n.allcoef,
                    dimnames = list(NULL, name.allcoef))

    ## precompute derivative for the variance
    dOmega.precomputed <- setNames(lapply(var.pattern, function(iPattern){
        iOut <- list(term1 = setNames(vector(mode = "list", length = n.varcoef), name.varcoef),
                     term2 = setNames(vector(mode = "list", length = n.varcoef), name.varcoef)
                     )
        for(iVarcoef in name.varcoef){
            iOut$term1[[iVarcoef]] <- tr(precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]])
            iOut$term2[[iVarcoef]] <- precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]] %*% precision[[iPattern]]
        }
        return(iOut)
    }), var.pattern)

    REML.num <- setNames(lapply(name.varcoef, function(i){
        matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef))
    }), name.varcoef)
    REML.denom <- matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef, name.mucoef))

    ## ** compute score
    for(iId in 1:n.cluster){ ## iId <- 7
        iPattern <- index.variance[iId]
        iIndex <- which(index.cluster==iId)
        iIndex <- iIndex[order(time.variance[iIndex])] ## re-order observations according to the variance-covariance matrix

        iResidual <- residuals[iIndex,,drop=FALSE]
        iX <- X[iIndex,,drop=FALSE]

        Score[iId,name.mucoef] <- t(iX) %*% precision[[iPattern]] %*% iResidual
        if(REML){
            REML.denom <- REML.denom + t(iX) %*% precision[[iPattern]] %*% iX
        }

        for(iVarcoef in name.varcoef){ ## iVarcoef <- name.varcoef[1]
            Score[iId,iVarcoef] <- -0.5 * dOmega.precomputed[[iPattern]]$term1[[iVarcoef]] + 0.5 * t(iResidual) %*% dOmega.precomputed[[iPattern]]$term2[[iVarcoef]] %*% iResidual
            ## Score[iId,iVarcoef] <- 0*Score[iId,iVarcoef]

            if(REML){
                REML.num[[iVarcoef]] <- REML.num[[iVarcoef]] + t(iX) %*% dOmega.precomputed[[iPattern]]$term2[[iVarcoef]] %*% iX
            }
        }
    }
    
    ## ** export
    if(indiv){
        return(Score)
    }else{
        Score <- colSums(Score)
        if(REML){
            REML.denom <- solve(REML.denom)
            Score[name.varcoef] <-  Score[name.varcoef] + 0.5 * sapply(REML.num, function(x){tr(REML.denom %*% x)})
            ## 0.5 tr((X\OmegaM1X)^-1 (X\OmegaM1 d\Omega \OmegaM1 X))
        }
        return(Score)
    }

}


##----------------------------------------------------------------------
### score.R ends here
