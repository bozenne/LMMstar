### score.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  5 2021 (12:59) 
## Version: 
## Last-Updated: mar 23 2021 (12:01) 
##           By: Brice Ozenne
##     Update #: 71
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * score.lmm (code)
##' @export
score.lmm <- function(x, data = NULL, p = NULL, transform = TRUE, indiv = FALSE, ...){

    if(is.null(data) && is.null(p) && indiv == FALSE && transform == TRUE){
        out <- x$score
    }else{
        if(!is.null(data) || transform == FALSE){
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
            index.variance <- design$index.vargroup
            index.cluster <- design$index.cluster
            X.var <- design$X.var
        }else{
            Y <- x$Y
            X <- x$design$X.mean
            index.variance <- x$design$index.vargroup
            index.cluster <- x$design$index.cluster
            X.var <- x$design$X.var
        }
        if(!is.null(p)){
            if(any(names(x$param$type) %in% names(p) == FALSE)){
                stop("Incorrect argument \'p\': missing parameter(s) \"",paste(names(x$param$type)[names(x$param$type) %in% names(p) == FALSE], collapse = "\" \""),"\".\n")
            }
            beta <- p[names(x$param$mu)]
            Omega <- attr(X.var,"FUN.Omega")(object = X.var, sigma = p[names(x$param$sigma)], k = p[names(x$param$k)], rho = p[names(x$param$cor)])
            precision <- lapply(Omega, solve)
        }else{
            beta <- x$param$mu
            precision <- x$OmegaM1
        }
        out <- .score(X = X, residuals = Y - X %*% beta, precision = precision,
                      index.variance = index.variance, index.cluster = index.cluster,
                      indiv = indiv, REML = object$method.fit=="REML")
    }
    
    return(out)
}

## * .score
.score <- function(X, residuals, precision, dOmega,
                   index.variance, index.cluster,
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
            iOut$term1[[iVarcoef]] <- -0.5 * tr(precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]])
            iOut$term2[[iVarcoef]] <- 0.5 * precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]] %*% precision[[iPattern]]
        }
        return(iOut)
    }), var.pattern)

    REML.det <- setNames(lapply(name.varcoef, function(i){
        list(numerator = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)),
             denominator = matrix(0, nrow = n.mucoef, ncol = n.mucoef, dimnames = list(name.mucoef,name.mucoef)))
    }), name.varcoef)
    
    ## ** compute score
    for(iId in 1:n.cluster){ ## iId <- 7
        iPattern <- index.variance[iId]
        iIndex <- which(index.cluster==iId)
        
        iResidual <- residuals[iIndex,,drop=FALSE]
        iX <- X[iIndex,,drop=FALSE]

        Score[iId,name.mucoef] <- t(iX) %*% precision[[iPattern]] %*% iResidual
        for(iVarcoef in name.varcoef){ ## iVarcoef <- name.varcoef[1]
            Score[iId,iVarcoef] <- dOmega.precomputed[[iPattern]]$term1[[iVarcoef]] + t(iResidual) %*% dOmega.precomputed[[iPattern]]$term2[[iVarcoef]] %*% iResidual

            ## Score[iId,iVarcoef] <- -0.5 * tr(precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]]) + 0.5 * t(iResidual) %*% precision[[iPattern]] %*% dOmega[[iPattern]][[iVarcoef]] %*% precision[[iPattern]] %*% iResidual
            if(REML){
                REML.det[[iVarcoef]]$numerator <- REML.det[[iVarcoef]]$numerator + t(iX) %*% (2*dOmega.precomputed[[iPattern]]$term2[[iVarcoef]]) %*% iX
                REML.det[[iVarcoef]]$denominator <- REML.det[[iVarcoef]]$denominator + t(iX) %*% precision[[iPattern]] %*% iX
            }
        }
    }
    
    ## ** export
    if(indiv){
        return(Score)
    }else{
        Score <- colSums(Score)
        if(REML){
            Score[name.varcoef] <-  Score[name.varcoef] + 0.5 * sapply(REML.det, function(x){tr(solve(x$denominator) %*% x$numerator)})
        }
        return(Score)
    }

}


##----------------------------------------------------------------------
### score.R ends here
